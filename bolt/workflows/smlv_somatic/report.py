import csv
import json
import pathlib


import click
import cyvcf2
import yaml


from ... import util
from ...common import constants
from ...common import pcgr


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

    script_name = pathlib.Path(__file__).stem
    setup_logging(output_dir, script_name)

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
    purple_data = parse_purple_purity_file(kwargs['purple_purity_fp'])

    pcgr_prep_fp = pcgr.prepare_vcf_somatic(
        kwargs['vcf_fp'],
        kwargs['tumor_name'],
        kwargs['normal_name'],
        output_dir,
    )

    pcgr.run_somatic(
        pcgr_prep_fp,
        kwargs['pcgr_data_dir'],
        output_dir,
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
