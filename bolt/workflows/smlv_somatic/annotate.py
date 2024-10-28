import logging
import os
import pathlib
import concurrent.futures
import click
import cyvcf2

from ... import util
from ...common import constants
from ...common import pcgr
from ...logging_config import setup_logging
logger = logging.getLogger(__name__)

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

@click.option('--output_dir', required=True, type=click.Path())

def entry(ctx, **kwargs):
    '''Annotate variants with information from several sources\f

    1. Set FILTER="PASS" for unfiltered variants
    2. Annotate with hotspot regions, high confidence regions, exclude regions, gnomAD, etc
    3. Annotate with panel of normal counts
    4. Prepare VCF for PCGR then annotate with clinical information using PCGR
    '''

    # Create output directory
    output_dir = pathlib.Path(kwargs['output_dir'])
    output_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

    setup_logging(output_dir)

    # Create a logger instance
    logger = logging.getLogger(__name__)

    # Set all FILTER="." to FILTER="PASS" as required by PURPLE
    filter_pass_fp = set_filter_pass(kwargs['vcf_fp'], kwargs['tumor_name'], output_dir)

    # Annotate with:
    #   - gnomAD [INFO/gnomAD_AF]
    #   - Hartwig hotspots [INFO/HMF_HOTSPOT]
    #   - ENCODE blocklist [INFO/ENCODE]
    #   - GIAB high confidence regions [INFO/GIAB_CONF]
    #   - Selected GA4GH/GIAB problem region stratifications [INFO/DIFFICULT_*]
    vcfanno_fp = general_annotations(
        filter_pass_fp,
        kwargs['tumor_name'],
        kwargs['threads'],
        kwargs['annotations_dir'],
        output_dir,
    )

    # Annotate with UMCCR panel of normals [INFO/PON_COUNT]
    # NOTE(SW): done separately from above as the variant identity for the INDEL PON operates only
    # on position rather than position /and/ reference + allele
    pon_fp = panel_of_normal_annotations(
        vcfanno_fp,
        kwargs['tumor_name'],
        kwargs['threads'],
        kwargs['pon_dir'],
        output_dir,
    )

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
    # Set selected data or full input

    # selection_data = select_variants(
        # pon_fp,
        # kwargs['tumor_name'],
        # kwargs['cancer_genes_fp'],
        # output_dir,
    # )
    pcgr_prep_fp = pcgr.prepare_vcf_somatic(
        pon_fp,
        kwargs['tumor_name'],
        kwargs['normal_name'],
        output_dir,
    )

    pcgr_output_dir = output_dir / 'pcgr'
    total_variants = util.count_vcf_records(pcgr_prep_fp)
    print(f"Total number of variants in the input VCF: {total_variants}")

    if total_variants > constants.MAX_SOMATIC_VARIANTS:
        vcf_chunks = util.split_vcf(
            pcgr_prep_fp,
            output_dir
        )
        pcgr_tsv_fp, pcgr_vcf_fp = run_somatic_chunck(
        vcf_chunks,
        kwargs['pcgr_data_dir'],
        output_dir,
        pcgr_output_dir,
        kwargs['threads'],
        kwargs['pcgr_conda'],
        kwargs['pcgrr_conda']
    )
    else:
        pcgr_tsv_fp, pcgr_vcf_fp = pcgr.run_somatic(
        pcgr_prep_fp,
        kwargs['pcgr_data_dir'],
        pcgr_output_dir,
        chunk_nbr=None,
        threads=kwargs['threads'],
        pcgr_conda=kwargs['pcgr_conda'],
        pcgrr_conda=kwargs['pcgrr_conda'],
    )

    # Transfer PCGR annotations to full set of variants
    pcgr.transfer_annotations_somatic(
        pon_fp,
        kwargs['tumor_name'],
        pcgr_vcf_fp,
        pcgr_tsv_fp,
        output_dir,
    )
    logger.info("Annotation process completed")


def run_somatic_chunck(vcf_chunks, pcgr_data_dir, output_dir, pcgr_output_dir, max_threads, pcgr_conda, pcgrr_conda):
    pcgr_tsv_files = []
    pcgr_vcf_files = []

    num_chunks = len(vcf_chunks)
    # Ensure we don't use more workers than available threads, and each worker has at least 2 threads
    max_workers = min(num_chunks, max_threads // 2)
    threads_quot, threads_rem = divmod(max_threads, num_chunks)
    threads_per_chunk = max(2, threads_quot)

    # Limit the number of workers to the smaller of num_chunks or max_threads // 2
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for chunk_number, vcf_file in enumerate(vcf_chunks, start=1):
            # Assign extra thread to the first 'threads_rem' chunks
            additional_thread = 1 if chunk_number <= threads_rem else 0
            total_threads = threads_per_chunk + additional_thread
            futures[executor.submit(pcgr.run_somatic, vcf_file, pcgr_data_dir, pcgr_output_dir, chunk_number, total_threads, pcgr_conda, pcgrr_conda)] = chunk_number

        for future in concurrent.futures.as_completed(futures):
            try:
                pcgr_tsv_fp, pcgr_vcf_fp = future.result()
                if pcgr_tsv_fp:
                    pcgr_tsv_files.append(pcgr_tsv_fp)
                if pcgr_vcf_fp:
                    pcgr_vcf_files.append(pcgr_vcf_fp)
            except Exception as e:
                print(f"Exception occurred: {e}")

    merged_vcf_fp, merged_tsv_fp = merging_pcgr_files(output_dir, pcgr_vcf_files, pcgr_tsv_files)
    return merged_tsv_fp, merged_vcf_fp

def merging_pcgr_files(output_dir, pcgr_vcf_files, pcgr_tsv_fp):
    # Step 3: Merge all chunk VCF files into a single file
    pcgr_dir = output_dir / 'pcgr/'
    pcgr_dir.mkdir(exist_ok=True)
    # Merge all TSV files into a single file in the pcgr directory    merged_tsv_fp = os.path.join(pcgr_dir, "nosampleset.pcgr_acmg.grch38.snvs_indels.tiers.tsv")
    merged_tsv_fp = os.path.join(pcgr_dir, "nosampleset.pcgr_acmg.grch38.snvs_indels.tiers.tsv")
    util.merge_tsv_files(pcgr_tsv_fp, merged_tsv_fp)
    # Step 5: Merge all VCF files into a single file in the pcgr directory
    merged_vcf_path = os.path.join(pcgr_dir, "nosampleset.pcgr_acmg.grch38")
    merged_vcf = util.merge_vcf_files(pcgr_vcf_files, merged_vcf_path)
    return merged_vcf, merged_tsv_fp

def set_filter_pass(input_fp, tumor_name, output_dir):
    output_fp = output_dir / f'{tumor_name}.set_filter_pass.vcf.gz'

    input_fh = cyvcf2.VCF(input_fp)
    output_fh = cyvcf2.Writer(output_fp, input_fh, 'wz')

    for record in input_fh:
        if record.FILTER is None:
            record.FILTER = 'PASS'
        output_fh.write_record(record)

    return output_fp


def general_annotations(input_fp, tumor_name, threads, annotations_dir, output_dir):
    toml_fp = pathlib.Path(annotations_dir) / 'vcfanno_annotations.toml'

    output_fp = output_dir / f'{tumor_name}.general_annotations.vcf.gz'

    command = fr'''
        vcfanno -p {threads} -base-path $(pwd) {toml_fp} {input_fp} | bcftools view -o {output_fp} && \
            bcftools index -t {output_fp}
    '''

    util.execute_command(command)
    return output_fp


def panel_of_normal_annotations(input_fp, tumor_name, threads, pon_dir, output_dir):
    toml_snp_fp = pathlib.Path(pon_dir) / 'vcfanno_snps.toml'
    toml_indel_fp = pathlib.Path(pon_dir) / 'vcfanno_indels.toml'

    output_fp = output_dir / f'{tumor_name}.pon.vcf.gz'

    threads_quot, threads_rem = divmod(threads, 2)

    command = fr'''
        vcfanno -p {threads_quot+threads_rem} -base-path $(pwd) {toml_snp_fp} {input_fp} | \
            vcfanno -permissive-overlap -p {threads_quot} -base-path $(pwd) {toml_indel_fp} /dev/stdin | \
            bcftools view -o {output_fp} && \
            bcftools index -t {output_fp}
    '''

    util.execute_command(command)
    return output_fp

def select_variants(input_fp, tumor_name, cancer_genes_fp, output_dir):
    # Exclude variants until we hopefully move the needle below the threshold


    # TODO(SW): make clear logs of this data and process
    # TODO(SW): generate in tempdir and keep only the one that passes; probably avoid doing this
    # right now


    if util.count_vcf_records(input_fp) <= constants.MAX_SOMATIC_VARIANTS:
        return {'selected': input_fp}

    pass_fps = exclude_nonpass(input_fp, tumor_name, output_dir)
    if util.count_vcf_records(pass_fps.get('filtered')) <= constants.MAX_SOMATIC_VARIANTS:
        return pass_fps

    pop_fps = exclude_population_variants(input_fp, tumor_name, output_dir)
    if util.count_vcf_records(pop_fps.get('filtered')) <= constants.MAX_SOMATIC_VARIANTS:
        return pop_fps

    genes_fps = exclude_non_cancer_genes(input_fp, tumor_name, cancer_genes_fp, output_dir)
    if util.count_vcf_records(genes_fps.get('filtered')) <= constants.MAX_SOMATIC_VARIANTS:

        # TODO(SW): log warning

        pass

    return genes_fps


def generic_exclude(input_fp, tumor_name, label, header_enum, process_fn, output_dir, **kwargs):
    selected_fp = output_dir / f'{tumor_name}.{label}.vcf.gz'
    filtered_fp = output_dir / f'{tumor_name}.{label}.filtered.vcf.gz'

    input_fh = cyvcf2.VCF(input_fp)
    util.add_vcf_header_entry(input_fh, header_enum)

    selected_fh = cyvcf2.Writer(selected_fp, input_fh, 'wz')
    filtered_fh = cyvcf2.Writer(filtered_fp, input_fh, 'wz')

    for record in input_fh:
        process_fn(record, selected_fh=selected_fh, filtered_fh=filtered_fh, **kwargs)

    return {'selected': selected_fp, 'filtered': filtered_fp, 'filter_name': header_enum.value}


def exclude_nonpass(input_fp, tumor_name, output_dir):
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

    return generic_exclude(input_fp, tumor_name, 'pass_select', header_enum, process_fn, output_dir)


def exclude_population_variants(input_fp, tumor_name, output_dir):
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

    return generic_exclude(input_fp, tumor_name, 'common_population', header_enum, process_fn, output_dir)


def exclude_non_cancer_genes(input_fp, tumor_name, bed_fp, output_dir):
    # Set FILTER to MAX_VARIANTS_KEY_GENES where variant is not in a key cancer gene or annotated
    # as a hotspot

    header_enum = constants.VcfFilter.MAX_VARIANTS_NON_CANCER_GENES
    gene_variants = get_cancer_gene_variants(input_fp, bed_fp, output_dir)

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

    return generic_exclude(input_fp, tumor_name, 'cancer_genes', header_enum, process_fn,
            output_dir, gene_variants=gene_variants)


def get_cancer_gene_variants(input_fp, bed_fp, output_dir):
    # Apply 1,000 bp padding to gene boundaries
    genes_variants_fp = output_dir / 'genes_variants.vcf.gz'
    util.execute_command(fr'''
        bcftools view \
            --regions-file <(awk 'BEGIN {{ OFS="\t" }} {{ $2-=1000; $3+=1000; print $0 }}' {bed_fp}) \
            --output {genes_variants_fp} \
            {input_fp}
    ''')

    gene_variants = set()
    for record in cyvcf2.VCF(genes_variants_fp):
        key = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        gene_variants.add(key)
    return gene_variants
