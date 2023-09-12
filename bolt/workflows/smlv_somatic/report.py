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

def entry(ctx, **kwargs):
    """Generate summary statistics and reports\f
    """

    # BCFtools stats
    bcftools_vcf_fp = bcftools_stats_prepare(kwargs['vcf_fp'], kwargs['tumor_name'])
    run_bcftools_stats(bcftools_vcf_fp, kwargs['tumor_name'])

    # Allele frequencies
    allele_frequencies(
        kwargs['tumor_name'],
        kwargs['vcf_fp'],
        kwargs['cancer_genes_fp'],
        kwargs['giab_regions_fp'],
        kwargs['genome_fp'],
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

    with open(f'{kwargs["tumor_name"]}.somatic.variant_counts.yaml', 'w') as fh:
        count_output = {
            'id': 'umccr',
            'data': { kwargs['tumor_name']: variant_count_data }
        }
        yaml.dump(count_output, fh, default_flow_style=False)

    # PCGR report
    purple_data = parse_purple_purity_file(kwargs['purple_purity_fp'])

    pcgr_prep_fp = pcgr.prepare_vcf_somatic(
        kwargs['vcf_fp'],
        f'{kwargs["tumor_name"]}.somatic',
        kwargs['tumor_name'],
        kwargs['normal_name'],
    )

    pcgr.run_somatic(
        pcgr_prep_fp,
        kwargs['pcgr_data_dir'],
        threads=kwargs['threads'],
        pcgr_conda=kwargs['pcgr_conda'],
        pcgrr_conda=kwargs['pcgrr_conda'],
        purity=purple_data['purity'],
        ploidy=purple_data['ploidy'],
        sample_id=kwargs['tumor_name'],
    )


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


def bcftools_stats_prepare(in_fp, tumor_name):
    in_fh = cyvcf2.VCF(in_fp)

    out_fp = f'{tumor_name}.somatic.bcftools_stats.vcf.gz'
    out_fh = cyvcf2.Writer(out_fp, in_fh, 'wz')

    tumor_index = in_fh.samples.index(tumor_name)
    for record in in_fh:
        # NOTE(SW): SAGE and DRAGEN quality scores are not comparable; we only get stats of DRAGEN
        # FORMAT/SQ
        if (tumor_sq_value := record.format('SQ')) is not None:
            record.QUAL = tumor_sq_value[tumor_index,0]
        elif record.INFO.get('SAGE_NOVEL') is not None:
            record.QUAL = None
        else:
            assert False

        out_fh.write_record(record)

    return out_fp


def run_bcftools_stats(vcf_fp, tumor_name):
    command = fr'''
        bcftools stats {vcf_fp} | \
            sed '6 s/{vcf_fp}$/{tumor_name}/' > {tumor_name}.somatic.bcftools_stats.txt
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


def allele_frequencies(tumor_name, vcf_fp, cancer_genes_fp, giab_regions_fp, genome_fp):
    command = fr'''
        bcftools annotate -Ob -x 'INFO' {vcf_fp} | \
            bcftools view -Ob -s {tumor_name} -T <(gzip -cd {giab_regions_fp}) | \
            bcftools norm -Ob -m - -f {genome_fp} | \
            bcftools sort -o af_variants.vcf.gz

        {{ echo af; bcftools query -f '[%AF]\n' af_variants.vcf.gz; }} > af_tumor.txt

        {{
            echo -e 'chrom\tpos\tid\tref\talt\taf';
            bcftools view -f .,PASS af_variants.vcf.gz | \
                bedtools intersect -a - -b {cancer_genes_fp} -header | \
                bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%AF]\n';
        }} > af_tumor_keygenes.txt
    '''

    util.execute_command(command)
