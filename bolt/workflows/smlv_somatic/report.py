import csv
import pathlib


import click
import cyvcf2
import yaml


from ... import util
from ...common import pcgr


@click.command(name='report')
@click.pass_context

@click.option('--tumor_name', required=True, type=str)
@click.option('--normal_name', required=True, type=str)

@click.option('--vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--vcf_unfiltered_fp', required=True, type=click.Path(exists=True))

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

    # Variant counts
    variant_counts_input = count_variant_types(kwargs['vcf_unfiltered_fp'])
    variant_counts_processed = count_variant_types(kwargs['vcf_fp'])

    # NOTE(SW): using pass variants only for now

    variant_count_proportions = dict()
    for k in {*variant_counts_input['pass'], *variant_counts_processed['pass']}:
        assert k not in variant_count_proportions
        if (count_input := variant_counts_input['pass'][k]) == 0:
            variant_count_proportions[k] = 0
        elif (count_processed := variant_counts_processed['pass'][k]) == 0:
            variant_count_proportions[k] = 0
        else:
            variant_count_proportions[k] = (count_input - count_processed) / count_input * 100

    variant_count_data = {
        'snps': variant_counts_processed['pass']['snps'],
        'indels': variant_counts_processed['pass']['indels'],
        'others': variant_counts_processed['pass']['others'],
        'filt_vars': variant_count_proportions['total'],
        'filt_snps': variant_count_proportions['snps'],
        'filt_indels': variant_count_proportions['indels'],
        'filt_others': variant_count_proportions['others'],
    }

    variant_counts_output_fp = output_dir / f'{kwargs["tumor_name"]}.somatic.variant_counts.yaml'
    with variant_counts_output_fp.open('w') as fh:
        count_output = {
            'id': 'umccr',
            'data': { kwargs['tumor_name']: variant_count_data }
        }
        yaml.dump(count_output, fh, default_flow_style=False)

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

    output_fp = f'{tumor_name}.somatic.bcftools_stats.vcf.gz'
    output_fh = cyvcf2.Writer(output_fp, input_fh, 'wz')

    tumor_index = input_fh.samples.index(tumor_name)
    for record in input_fh:
        # NOTE(SW): SAGE and DRAGEN quality scores are not comparable; we only get stats of DRAGEN
        # FORMAT/SQ
        if (tumor_sq_value := record.format('SQ')) is not None:
            record.QUAL = tumor_sq_value[tumor_index,0]
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
            sed '6 s/{input_fp}$/{tumor_name}/' > {output_fp}
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
