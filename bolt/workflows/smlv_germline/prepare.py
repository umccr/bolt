import pathlib


import click
import cyvcf2


from ... import util
from ...common import constants


@click.command(name='prepare')
@click.pass_context

@click.option('--vcf_fp', required=True, type=click.Path(exists=True))

@click.option('--gene_panel_fp', required=True, type=click.Path(exists=True))

@click.option('--output_fp', required=True, type=click.Path())

def entry(ctx, **kwargs):
    '''Prepare germline variants for processing\f

    1. Select passing variants in the given gene panel
    '''

    pass_vcf = select_pass(kwargs['vcf_fp'])
    panel_vcf = select_gene_panel_variants(pass_vcf, kwargs['gene_panel_fp'], kwargs['output_fp'])


def select_pass(in_fp):
    in_fn = pathlib.Path(in_fp).name
    out_fp = f'{in_fn.replace(".vcf.gz", ".pass.vcf.gz")}'
    command = fr'''
        bcftools view -f PASS,. -o {out_fp} {in_fp} && \
            bcftools index -t {out_fp}
    '''
    util.execute_command(command)
    return out_fp


def select_gene_panel_variants(vcf_fp, panel_fp, output_fp):

    command = fr'''
        bcftools view --regions-file {panel_fp} -o {output_fp} {vcf_fp};
        bcftools index -t {output_fp};
    '''
    util.execute_command(command)

    return output_fp
