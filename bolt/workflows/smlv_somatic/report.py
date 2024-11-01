import collections
import csv
import json
import logging
import pathlib


import click
import cyvcf2
import yaml


from ... import util
from ...common import constants
from ...common import pcgr
from ...logging_config import setup_logging

logger = logging.getLogger(__name__)



@click.command(name='report')
@click.pass_context

@click.option('--tumor_name', required=True, type=str)
@click.option('--normal_name', required=True, type=str)

@click.option('--vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--vcf_filters_fp', required=True, type=click.Path(exists=True))
@click.option('--vcf_dragen_fp', required=True, type=click.Path(exists=True))

@click.option('--pcgr_conda', required=False, type=str)
@click.option('--pcgrr_conda', required=False, type=str)

@click.option('--pcgr_data_dir', required=False, type=str)
@click.option('--vep_dir', required=True, type=click.Path(exists=True))
@click.option('--purple_purity_fp', required=True, type=click.Path(exists=True))

@click.option('--cancer_genes_fp', required=True, type=click.Path(exists=True))
@click.option('--giab_regions_fp', required=True, type=click.Path(exists=True))
@click.option('--genome_fp', required=True, type=click.Path(exists=True))

@click.option('--threads', required=True, type=int, default=1)

@click.option('--output_dir', required=True, type=click.Path())

def entry(ctx, **kwargs):
    """Generate summary statistics and reports\f
    """

    # Create output directory
    output_dir = pathlib.Path(kwargs['output_dir'])
    output_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

    setup_logging(output_dir, "smlv_somatic_annotate.py")
    logger = logging.getLogger(__name__)


    # BCFtools stats
    bcftools_vcf_fp = bcftools_stats_prepare(kwargs['vcf_fp'], kwargs['tumor_name'], output_dir)
    run_bcftools_stats(bcftools_vcf_fp, kwargs['tumor_name'], output_dir)

    # Allele frequencies
    allele_frequencies(
        kwargs['vcf_fp'],
        kwargs['tumor_name'],
        kwargs['cancer_genes_fp'],
        kwargs['giab_regions_fp'],
        kwargs['genome_fp'],
        output_dir,
    )

    # Variant type counts
    # NOTE(SW): this is intended to preserve counts in the MultiQC report
    variant_counts_types_dragen = count_variant_types(kwargs['vcf_dragen_fp'])
    variant_counts_types_bolt= count_variant_types(kwargs['vcf_fp'])

    # NOTE(SW): using pass variants only for now

    variant_count_type_proportions = dict()
    for k in {*variant_counts_types_dragen['pass'], *variant_counts_types_bolt['pass']}:
        assert k not in variant_count_type_proportions
        if (count_dragen := variant_counts_types_dragen['pass'][k]) == 0:
            variant_count_type_proportions[k] = 0
        elif (count_bolt := variant_counts_types_bolt['pass'][k]) == 0:
            variant_count_type_proportions[k] = 0
        else:
            variant_count_type_proportions[k] = (count_dragen - count_bolt) / count_dragen * 100

    variant_count_data = {
        'snps': variant_counts_types_bolt['pass']['snps'],
        'indels': variant_counts_types_bolt['pass']['indels'],
        'others': variant_counts_types_bolt['pass']['others'],
        'filt_vars': variant_count_type_proportions['total'],
        'filt_snps': variant_count_type_proportions['snps'],
        'filt_indels': variant_count_type_proportions['indels'],
        'filt_others': variant_count_type_proportions['others'],
    }

    variant_counts_type_output_fn = f'{kwargs["tumor_name"]}.somatic.variant_counts_type.yaml'
    variant_counts_type_output_fp = output_dir / variant_counts_type_output_fn
    with variant_counts_type_output_fp.open('w') as fh:
        count_output = {
            'id': 'umccr',
            'data': { kwargs['tumor_name']: variant_count_data }
        }
        yaml.dump(count_output, fh, default_flow_style=False)

    # Variant process counts
    # NOTE(SW): this is intended to preserve counts in the Cancer Report
    variant_counts_process = count_variant_process(kwargs['vcf_filters_fp'])
    variant_counts_process_fn = f'{kwargs["tumor_name"]}.somatic.variant_counts_process.json'
    variant_counts_process_fp = output_dir / variant_counts_process_fn
    with variant_counts_process_fp.open('w') as fh:
        json.dump(variant_counts_process, fh, indent=4)
        fh.write('\n')

    # PCGR report
    if variant_counts_process['filter_pass'] <= constants.MAX_SOMATIC_VARIANTS:
        pcgr_input_vcf_fp = kwargs['vcf_fp']
    else:
        pcgr_input_vcf_fp = select_pcgr_variants(
            kwargs['vcf_fp'],
            kwargs['cancer_genes_fp'],
            kwargs['tumor_name'],
            output_dir,
        )

    purple_data = parse_purple_purity_file(kwargs['purple_purity_fp'])

    pcgr_prep_fp = pcgr.prepare_vcf_somatic(
        pcgr_input_vcf_fp,
        kwargs['tumor_name'],
        kwargs['normal_name'],
        output_dir,
    )

    pcgr_output_dir = output_dir / 'pcgr'
    pcgr.run_somatic(
        pcgr_prep_fp,
        kwargs['pcgr_data_dir'],
        pcgr_output_dir,
        threads=kwargs['threads'],
        pcgr_conda=kwargs['pcgr_conda'],
        pcgrr_conda=kwargs['pcgrr_conda'],
        purity=purple_data['purity'],
        ploidy=purple_data['ploidy'],
        sample_id=kwargs['tumor_name'],
    )


def bcftools_stats_prepare(input_fp, tumor_name, output_dir):
    input_fh = cyvcf2.VCF(input_fp)

    output_fp = output_dir / f'{tumor_name}.somatic.bcftools_stats.vcf.gz'
    output_fh = cyvcf2.Writer(output_fp, input_fh, 'wz')

    tumor_index = input_fh.samples.index(tumor_name)
    for record in input_fh:
        # NOTE(SW): SAGE and DRAGEN quality scores are not comparable; we only get stats of DRAGEN
        # FORMAT/SQ
        if (tumor_sq_value := record.format('SQ')) is not None:
            # Round SQ so that BCFtools stats uses integers on x-axis
            record.QUAL = round(tumor_sq_value[tumor_index,0])
        elif record.INFO.get('SAGE_NOVEL') is not None:
            record.QUAL = None
        else:
            assert False

        output_fh.write_record(record)

    return output_fp


def run_bcftools_stats(input_fp, tumor_name, output_dir):
    output_fp = output_dir / f'{tumor_name}.somatic.bcftools_stats.txt'
    command = fr'''
        bcftools stats {input_fp} | \
            sed '6 s#{input_fp}$#{tumor_name}#' > {output_fp}
    '''
    util.execute_command(command)


def allele_frequencies(input_fp, tumor_name, cancer_genes_fp, giab_regions_fp, genome_fp, output_dir):
    af_vcf_fp = output_dir / 'af_variants.vcf.gz'
    af_global_output_fp = output_dir / 'af_tumor.txt'
    af_keygenes_output_fp = output_dir / 'af_tumor_keygenes.txt'

    command = fr'''
        bcftools annotate -Ob -x 'INFO' {input_fp} | \
            bcftools view -Ob -s {tumor_name} -T <(gzip -cd {giab_regions_fp}) | \
            bcftools norm -Ob -m - -f {genome_fp} | \
            bcftools sort -o {af_vcf_fp}

        {{ echo af; bcftools query -f '[%AF]\n' {af_vcf_fp}; }} > {af_global_output_fp}

        {{
            echo -e 'chrom\tpos\tid\tref\talt\taf';
            bcftools view -f .,PASS {af_vcf_fp} | \
                bedtools intersect -a - -b {cancer_genes_fp} -header | \
                bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%AF]\n';
        }} > {af_keygenes_output_fp}
    '''

    util.execute_command(command)


def count_variant_types(vcf_fp):

    counts_templ = {'snps': 0, 'indels': 0, 'others': 0, 'total_incl_unfiltered': 0}

    counts = {
        'pass': counts_templ.copy(),
        'nonpass': counts_templ.copy(),
    }

    for record in cyvcf2.VCF(vcf_fp):

        variant_filter = str()
        if record.FILTER is None:
            variant_filter = 'pass'
        else:
            variant_filter = 'nonpass'

        if record.is_snp:
            counts[variant_filter]['snps'] += 1
        elif record.is_indel:
            counts[variant_filter]['indels'] += 1
        else:
            counts[variant_filter]['others'] += 1

    for variant_filter, data in counts.items():
        total = sum(v for v in data.values())
        counts[variant_filter]['total'] = total

    return counts


def count_variant_process(vcf_fp):
    # Set filter groups
    sage_add_info = {
        constants.VcfInfo.SAGE_RESCUE.value,
        constants.VcfInfo.SAGE_NOVEL.value,
    }

    bolt_annotation_filters = {
        constants.VcfFilter.MAX_VARIANTS_NON_PASS.value,
        constants.VcfFilter.MAX_VARIANTS_GNOMAD.value,
        constants.VcfFilter.MAX_VARIANTS_NON_CANCER_GENES.value,
    }

    # NOTE(SW): this assumes both all and /only/ bolt filters are defined in constants.VcfFilter
    bolt_filters = {f.value for f in constants.VcfFilter} - bolt_annotation_filters

    # Count table
    counts = {
        'dragen': 0,
        'sage': 0,
        'annotated': 0,
        'filter_pass': 0,
    }

    for record in cyvcf2.VCF(vcf_fp):

        # Restore existing filters to get accurate DRAGEN counts with rescued variants
        rescued_filters = record.INFO.get(constants.VcfInfo.RESCUED_FILTERS_EXISTING.value)
        if rescued_filters:
            assert not record.FILTER
            record.FILTER = rescued_filters.replace(',', ';')

        # Separate filters into useful groups
        record_filters_dragen = list()
        record_filters_bolt = list()
        record_filters_bolt_annotation = list()

        for filter_str in record.FILTERS:

            if filter_str == 'PASS':
                continue
            elif filter_str in bolt_annotation_filters:
                record_filters_bolt_annotation.append(filter_str)
            elif filter_str in bolt_filters:
                record_filters_bolt.append(filter_str)
            else:
                record_filters_dragen.append(filter_str)

        # Ensure no FILTERs for annotation 'max variants' have been missed
        assert not any(f.startswith('max_variants_') for f in record_filters_bolt)

        # Begin counting
        # Do not include novel SAGE calls in DRAGEN counts
        if record.INFO.get(constants.VcfInfo.SAGE_NOVEL.value) is None:
            counts['dragen'] += 1

        # All DRAGEN variants are passed to SAGE
        counts['sage'] += 1

        # Variants can be excluded for annotation when there are too many present
        if not record_filters_bolt_annotation:
            counts['annotated'] += 1

        # Filtering removes and rescues
        if not record.FILTER or rescued_filters:
            counts['filter_pass'] += 1

    return counts


def select_pcgr_variants(vcf_fp, cancer_genes_fp, tumor_name, output_dir):
    # Annotate variants in UMCCR somatic gene panel
    fp_annotated_out = output_dir / f'{tumor_name}.umccr_panel_variants_annotated.vcf.gz'
    util.execute_command(fr'''
        bcftools annotate \
            --annotations <(awk 'BEGIN {{ OFS="\t" }} {{ print $1, $2-2000, $3+2000, "1" }}' {cancer_genes_fp}) \
            --header-line '{util.get_vcf_header_line(constants.VcfInfo.PANEL)}' \
            --columns CHROM,FROM,TO,{constants.VcfInfo.PANEL.value} \
            --output {fp_annotated_out} \
            {vcf_fp}
    ''')

    # Set filter category for each variant
    variants_sorted = collections.defaultdict(list)
    for variant_count, variant in enumerate(cyvcf2.VCF(fp_annotated_out), 1):
        if any(variant.INFO.get(e) for e in constants.RETAIN_FIELDS_FILTERING):
            continue

        data = pcgr.get_variant_filter_data(variant)
        variant_filter = pcgr.determine_filter(data)
        assert variant_filter

        filter_category = (data['tier'], *variant_filter)
        variant_repr = pcgr.get_variant_repr(variant)
        variants_sorted[filter_category].append(variant_repr)


    # Determine the set of filter categories to come under the PCGR 500,000 variant threshold
    filter_sum = 0
    filter_categories = list()
    for key in pcgr.get_ordering():

        if (variant_count - filter_sum) <= constants.MAX_SOMATIC_VARIANTS:
            break

        filter_sum += len(variants_sorted.get(key, []))
        filter_categories.append(key)

    # Set FILTERS and write out records
    filter_variants = set()
    for key in filter_categories:
        filter_variants.update(variants_sorted[key])

    fh_in = cyvcf2.VCF(vcf_fp)
    util.add_vcf_header_entry(fh_in, constants.VcfFilter.PCGR_COUNT_LIMIT)

    # NOTE(SW): creating an additional VCF with all records for traceability
    fp_out = output_dir / f'{tumor_name}.pcgr_hypermutated.pass.vcf.gz'
    fp_set_out = output_dir / f'{tumor_name}.pcgr_hypermutated.filters_set.vcf.gz'

    fh_out = cyvcf2.Writer(fp_out, fh_in)
    fh_set_out = cyvcf2.Writer(fp_set_out, fh_in)

    for variant in fh_in:
        variant_repr = pcgr.get_variant_repr(variant)
        if variant_repr not in filter_variants:
            # Write only passing
            fh_out.write_record(variant)
        else:
            variant.FILTER = constants.VcfFilter.PCGR_COUNT_LIMIT.value
        # Write all variants including those with FILTER set
        fh_set_out.write_record(variant)

    fh_out.close()
    fh_set_out.close()

    return fp_out


def parse_purple_purity_file(fp):
    with open(fp, 'r') as fh:
        entries = list(csv.DictReader(fh, delimiter='\t'))
    assert len(entries) == 1
    [entry] = entries

    data = {
        'purity': entry['purity'],
        'ploidy': entry['ploidy'],
    }
    return data
