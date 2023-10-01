import csv
import pathlib


import click
import cyvcf2


from ... import util
from ...common import constants
from ...common import pcgr


@click.command(name='annotate')
@click.pass_context
@click.option('--tumor_name', required=True, type=str)
@click.option('--normal_name', required=True, type=str)

@click.option('--vcf_fp', required=True, type=click.Path(exists=True))

@click.option('--cancer_genes_fp', required=True, type=click.Path(exists=True))
@click.option('--annotations_dir', required=True, type=click.Path(exists=True))
@click.option('--pon_dir', required=True, type=click.Path(exists=True))

@click.option('--pcgr_data_dir', required=True, type=click.Path(exists=True))

@click.option('--pcgr_conda', required=False, type=str)
@click.option('--pcgrr_conda', required=False, type=str)

@click.option('--threads', required=False, default=4, type=int)
def entry(ctx, **kwargs):
    '''Annotate variants with information from several sources\f

    1. Set FILTER="PASS" for unfiltered variants
    2. Annotate with hotspot regions, high confidence regions, exclude regions, gnomAD, etc
    3. Annotate with panel of normal counts
    4. Prepare VCF for PCGR then annotate with clinical information using PCGR
    '''
    fn_prefix = get_filename_prefix(kwargs['vcf_fp'])

    # Set all FILTER="." to FILTER="PASS" as required by PURPLE
    filter_pass_fp = set_filter_pass(kwargs['vcf_fp'], fn_prefix)

    # Annotate with:
    #   - gnomAD r2.1 [INFO/gnomAD_AF]
    #   - Merged hotspots [INFO/HMF_HOTSPOT]
    #   - GIAB high confidence regions from HMF [INFO/HMF_GIAB_CONF]
    #   - ENCODE blocklist [INFO/ENCODE]
    #   - UCSC segmental duplications [INFO/SEGDUP]
    #   - GA4GH genome stratifications [INFO/TRICKY_*]
    #   - One other LCR file (competes with GA4GH LCR) [INFO/TRICKY_LCR]
    vcfanno_fp = general_annotations(filter_pass_fp, fn_prefix, kwargs['threads'], kwargs['annotations_dir'])

    # Annotate with UMCCR panel of normals [INFO/PON_COUNT]
    # NOTE(SW): done separately from above as the variant identity for the INDEL PON operates only
    # on position rather than position /and/ reference + allele
    pon_fp = panel_of_normal_annotations(vcfanno_fp, fn_prefix, kwargs['threads'], kwargs['pon_dir'])

    # Annotate with cancer-related and functional information from a range of sources using PCGR
    #   - Select variants to process - there is an upper limit for PCGR of around 500k
    #   - Set tumor and normal AF and DP in INFO for PCGR and remove all other annotations
    #   - Run PCGR on minimal VCF (pcgr_prep_fp)
    #   - Transfer selected PCGR annotations to unfiltered VCF (selected_fp)
    #       - PCGR ACMG TIER [INFO/PCGR_TIER]
    #       - VEP consequence [INFO/PCR_CSQ]
    #       - Known mutation hotspot [INFO/PCGR_MUTATION_HOTSPOT]
    #       - ClinVar clinical significant [INFO/PCGR_CLINVAR_CLNSIG]
    #       - Hits in COSMIC [INFO/PCGR_COSMIC_COUNT]
    #       - Hits in TCGA [INFO/PCGR_TCGA_PANCANCER_COUNT]
    #       - Hits in PCAWG [INFO/PCGR_ICGC_PCAWG_COUNT]
    selection_data = select_variants(pon_fp, fn_prefix, kwargs['cancer_genes_fp'])

    # NOTE(SW): here I set the correct filepath when filtering was not required
    if not (pcgr_prep_input_fp := selection_data.get('filtered')):
        pcgr_prep_input_fp = selection_data['selected']
    pcgr_prep_fp = pcgr.prepare_vcf_somatic(pcgr_prep_input_fp, fn_prefix, kwargs['tumor_name'], kwargs['normal_name'])
    pcgr_dir = pcgr.run_somatic(
        pcgr_prep_fp,
        kwargs['pcgr_data_dir'],
        threads=kwargs['threads'],
        pcgr_conda=kwargs['pcgr_conda'],
        pcgrr_conda=kwargs['pcgrr_conda'],
    )

    pcgr.transfer_annotations_somatic(
        selection_data['selected'],
        kwargs['tumor_name'],
        pcgr_dir,
        selection_data.get('filter_name'),
    )


def get_filename_prefix(vcf_fp):
    vcf_fn = pathlib.Path(vcf_fp).name
    return vcf_fn.replace('.vcf.gz', '')


def set_filter_pass(vcf_fp, fn_prefix):
    fp_out = f'{fn_prefix}.set_filter_pass.vcf.gz'

    in_fh = cyvcf2.VCF(vcf_fp)
    out_fh = cyvcf2.Writer(fp_out, in_fh, 'wz')

    for record in in_fh:
        if record.FILTER is None:
            record.FILTER = 'PASS'
        out_fh.write_record(record)

    return fp_out


def general_annotations(in_fp, fn_prefix, threads, annotations_dir):
    toml_fp = pathlib.Path(annotations_dir) / 'vcfanno_annotations.toml'
    out_fp = f'{fn_prefix}.general.vcf.gz'
    command = fr'''
        vcfanno -p {threads} -base-path $(pwd) {toml_fp} {in_fp} | bcftools view -o {out_fp} && \
            bcftools index -t {out_fp}
    '''
    util.execute_command(command)
    return out_fp


def panel_of_normal_annotations(in_fp, fn_prefix, threads, pon_dir):
    toml_snp_fp = pathlib.Path(pon_dir) / 'vcfanno_snps.toml'
    toml_indel_fp = pathlib.Path(pon_dir) / 'vcfanno_indels.toml'

    out_fp = f'{fn_prefix}.pon.vcf.gz'

    threads_quot, threads_rem = divmod(threads, 2)

    command = fr'''
        vcfanno -p {threads_quot+threads_rem} -base-path $(pwd) {toml_snp_fp} {in_fp} | \
            vcfanno -permissive-overlap -p {threads_quot} -base-path $(pwd) {toml_indel_fp} /dev/stdin | \
            bcftools view -o {out_fp} && \
            bcftools index -t {out_fp}
    '''
    util.execute_command(command)
    return out_fp


def select_variants(in_fp, fn_prefix, cancer_genes_fp):
    # Exclude variants until we hopefully move the needle below the threshold


    # TODO(SW): make clear logs of this data and process
    # TODO(SW): generate in tempdir and keep only the one that passes; probably avoid doing this
    # right now


    if util.count_vcf_records(in_fp) <= constants.MAX_SOMATIC_VARIANTS:
        return {'selected': in_fp}

    pass_fps = exclude_nonpass(in_fp, fn_prefix)
    if util.count_vcf_records(pass_fps.get('filtered')) <= constants.MAX_SOMATIC_VARIANTS:
        return pass_fps

    pop_fps = exclude_population_variants(in_fp, fn_prefix)
    if util.count_vcf_records(pop_fps.get('filtered')) <= constants.MAX_SOMATIC_VARIANTS:
        return pop_fps

    genes_fps = exclude_non_cancer_genes(in_fp, fn_prefix, cancer_genes_fp)
    if util.count_vcf_records(genes_fps.get('filtered')) <= constants.MAX_SOMATIC_VARIANTS:

        # TODO(SW): log warning

        pass

    return genes_fps


def generic_exclude(in_fp, fn_prefix, label, header_enum, process_fn, **kwargs):
    selected_fp = f'{fn_prefix}.{label}.vcf.gz'
    filtered_fp = f'{fn_prefix}.{label}.filtered.vcf.gz'

    in_fh = cyvcf2.VCF(in_fp)
    util.add_vcf_header_entry(in_fh, header_enum)

    selected_fh = cyvcf2.Writer(selected_fp, in_fh, 'wz')
    filtered_fh = cyvcf2.Writer(filtered_fp, in_fh, 'wz')

    for record in in_fh:
        process_fn(record, selected_fh=selected_fh, filtered_fh=filtered_fh, **kwargs)

    return {'selected': selected_fp, 'filtered': filtered_fp, 'filter_name': header_enum.value}


def exclude_nonpass(in_fp, fn_prefix):
    # Exclude variants that don't have FILTER=PASS and are not in a hotspot

    header_enum = constants.VcfFilter.MAX_VARIANTS_NON_PASS

    def process_fn(record, selected_fh, filtered_fh):
        # NOTES(SW): sanity check to ensure that if INFO/SAGE_HOTSPOT is present that FILTER=PASS
        if record.INFO.get(constants.VcfInfo.SAGE_HOTSPOT.value) is not None:
            assert not record.FILTER
        # Write to filtered_fp if passes criteria: in a hotspot or FILTER=PASS
        # Otherwise update FILTER appropriately then write to selected_fp
        if record.INFO.get(constants.VcfInfo.HMF_HOTSPOT.value) is not None or not record.FILTER:
            filtered_fh.write_record(record)
        else:
            assert 'PASS' not in record.FILTERS
            record.FILTER = ';'.join([*record.FILTERS, header_enum.value])
        selected_fh.write_record(record)

    return generic_exclude(in_fp, fn_prefix, 'pass_select', header_enum, process_fn)


def exclude_population_variants(in_fp, fn_prefix):
    # Set FILTER to MAX_VARIANTS_GNOMAD where INFO/gnomAD_AF is greater than 0.01 and variant is
    # not annotated as a hotspot

    header_enum = constants.VcfFilter.MAX_VARIANTS_GNOMAD

    def process_fn(record, selected_fh, filtered_fh):
        # NOTE(SW): for cases where gnomAD AF is not available we consider that the variant is not common
        gnomad_af = record.INFO.get(constants.VcfInfo.GNOMAD_AF.value, 0)
        is_population_common = float(gnomad_af) >= constants.MAX_SOMATIC_VARIANTS_GNOMAD_FILTER
        # Write to filtered_fp if passes criteria: in a hotspot or not a common population variant
        # Otherwise update FILTER appropriately then write to selected_fp
        if record.INFO.get(constants.VcfInfo.HMF_HOTSPOT.value) is not None or not is_population_common:
            filtered_fh.write_record(record)
        else:
            existing_filters = [e for e in record.FILTERS if e != 'PASS']
            record.FILTER = ';'.join([*existing_filters, header_enum.value])
        selected_fh.write_record(record)

    return generic_exclude(in_fp, fn_prefix, 'common_population', header_enum, process_fn)


def exclude_non_cancer_genes(in_fp, fn_prefix, bed_fp):
    # Set FILTER to MAX_VARIANTS_KEY_GENES where variant is not in a key cancer gene or annotated
    # as a hotspot

    header_enum = constants.VcfFilter.MAX_VARIANTS_NON_CANCER_GENES
    gene_variants = get_cancer_gene_variants(in_fp, bed_fp)

    def process_fn(record, selected_fh, filtered_fh, *, gene_variants):
        # Write to filtered_fp if passes criteria: in a hotspot or associated with a cancer gene
        # Otherwise update FILTER appropriately then write to selected_fp
        key = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        if record.INFO.get(constants.VcfInfo.HMF_HOTSPOT.value) is not None or key in gene_variants:
            filtered_fh.write_record(record)
        else:
            existing_filters = [e for e in record.FILTERS if e != 'PASS']
            record.FILTER = ';'.join([*existing_filters, header_enum.value])
        selected_fh.write_record(record)

    return generic_exclude(in_fp, fn_prefix, 'cancer_genes', header_enum, process_fn,
            gene_variants=gene_variants)


def get_cancer_gene_variants(in_fp, bed_fp):
    # NOTE(SW): current UMCCR cancer gene list (umccr_cancer_genes.hg38.ensembl107.sort.bed)
    # already contains 1,000 bp padding
    genes_variants_fp = pathlib.Path('genes_variants.vcf.gz')
    util.execute_command(f'bcftools view --regions-file {bed_fp} -o {genes_variants_fp} {in_fp}')

    gene_variants = set()
    for record in cyvcf2.VCF(genes_variants_fp):
        key = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        gene_variants.add(key)
    return gene_variants
