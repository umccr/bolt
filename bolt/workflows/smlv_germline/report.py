import pathlib
import yaml


import click


from ... import util
from ...common import pcgr


@click.command(name='report')
@click.pass_context
@click.option('--normal_name', required=True, type=str)

@click.option('--vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--vcf_unfiltered_fp', required=True, type=click.Path(exists=True))

@click.option('--pcgr_conda', required=False, type=str)
@click.option('--pcgrr_conda', required=False, type=str)

@click.option('--germline_panel_list_fp', required=True, type=click.Path(exists=True))
@click.option('--pcgr_data_dir', required=True, type=click.Path(exists=True))

@click.option('--threads', required=True, type=int, default=1)

def entry(ctx, **kwargs):
    """Generate summary statistics and reports\f
    """

    # BCFtools stats
    run_bcftool_stats(kwargs['vcf_unfiltered_fp'], kwargs['normal_name'])

    # Variant counts
    variant_counts_input = count_variants(kwargs['vcf_unfiltered_fp'])
    variant_counts_processed = count_variants(kwargs['vcf_fp'])

    variant_count_data = {
        'germline': variant_counts_input,
        'germline_predispose': variant_counts_processed,
    }

    with open(f'{kwargs["normal_name"]}.germline.variant_counts.yaml', 'w') as fh:
        count_output = {
            'id': 'umccr',
            'data': { kwargs['normal_name']: variant_count_data }
        }
        yaml.dump(count_output, fh, default_flow_style=False)


    # CPSR report
    vcf_norm_fp = split_multiallelic_records(kwargs['vcf_fp'])

    cpsr_prep_fp = pcgr.prepare_vcf_germline(
        vcf_norm_fp,
        f'{kwargs["normal_name"]}.germline',
        kwargs['normal_name'],
    )

    cpsr_dir = pcgr.run_germline(
        cpsr_prep_fp,
        kwargs['germline_panel_list_fp'],
        kwargs['pcgr_data_dir'],
        threads=kwargs['threads'],
        pcgr_conda=kwargs['pcgr_conda'],
        pcgrr_conda=kwargs['pcgrr_conda'],
        sample_id=kwargs['normal_name'],
    )

    pcgr.transfer_annotations_germline(
        vcf_norm_fp,
        kwargs['normal_name'],
        cpsr_dir,
    )


def run_bcftool_stats(vcf_fp, normal_name):
    command = fr'''
        bcftools stats -f PASS,. {vcf_fp} | \
            sed '6 s/{vcf_fp}$/{normal_name}/' > {normal_name}.germline.bcftools_stats.txt
    '''
    util.execute_command(command)


def count_variants(fp):
    process = util.execute_command(f'bcftools view -H -f PASS,. {fp} | wc -l')
    return process.stdout.strip()


def split_multiallelic_records(fp):
    fp_out = fp.replace('.vcf.gz', '.norm.vcf.gz')
    command = fr'''
        bcftools norm -m - -o {fp_out} {fp}
    '''
    util.execute_command(command)
    return fp_out
