import pathlib


import click


from ... import util
from ...common import constants


@click.command(name='prepare')
@click.pass_context
@click.option('--tumor_name', required=True, type=str)
@click.option('--normal_name', required=True, type=str)

@click.option('--vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--somatic_vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--gene_panel_fp', required=True, type=click.Path(exists=True))
@click.option('--output_fp', required=True, type=click.Path())
def entry(ctx, **kwargs):
    '''Prepare germline variants for processing\f

    1. Add germline variants identified in somatic calls
    2. Select passing variants in the given gene panel
    '''
    pass_vcf = select_pass(kwargs['vcf_fp'])
    combined_vcf = add_germline_leakage(pass_vcf, kwargs['somatic_vcf_fp'], kwargs['tumor_name'], kwargs['normal_name'])
    panel_vcf = select_gene_panel_variants(combined_vcf, kwargs['gene_panel_fp'], kwargs['output_fp'])


def select_pass(in_fp):
    in_fn = pathlib.Path(in_fp).name
    out_fp = f'{in_fn.replace(".vcf.gz", ".pass.vcf.gz")}'
    command = fr'''
        bcftools view -f PASS,. -o {out_fp} {in_fp} && \
            bcftools index -t {out_fp}
    '''
    util.execute_command(command)
    return out_fp


def add_germline_leakage(germline_vcf, somatic_vcf, tumor_name, normal_name):

    # NOTE(SW): in order to traceback, must keep INFO/GNOMAD_COMMON and INFO/PON_COUNT; obvious also need to
    # include INFO/GERMLINE_LEAKAGE

    # NOTE(SW): we are purposefully including germline leakage variants that do not have a PASS,.
    # filter since this is what is done in Umccrise. I want to review/discuss the rationale driving
    # this aspect.

    info_default = [
        'INFO/DP',
        'INFO/END',
        'INFO/FractionInformativeReads',
        'INFO/MQ',
        'INFO/MQRankSum',
        'INFO/ReadPosRankSum',
        'INFO/hotspot',
    ]

    info_retain = [
        util.get_qualified_vcf_annotation(constants.VcfInfo.GERMLINE_LEAKAGE),
        util.get_qualified_vcf_annotation(constants.VcfInfo.PON_COUNT),
        'INFO/gnomAD_AF',
    ]

    format_default = [
        'FORMAT/AD',
        'FORMAT/AF',
        'FORMAT/DP',
        'FORMAT/F1R2',
        'FORMAT/F2R1',
        'FORMAT/GT',
        'FORMAT/MB',
        'FORMAT/PS',
        'FORMAT/SB',
        'FORMAT/SQ',
    ]

    out_fp = f'{germline_vcf.replace(".vcf.gz", ".include_leakage.vcf.gz")}'

    exclude_str = ',^'.join([*info_default, *info_retain, *format_default])

    command = fr'''
        bcftools view -s {tumor_name} -i 'INFO/GERMLINE_LEAKAGE=1' {somatic_vcf} | \
            sed '/^#CHROM/ s/{tumor_name}/{normal_name}/' | \
            bcftools annotate -x "{exclude_str}" -o germline_leakage.vcf.gz;
        bcftools index -t germline_leakage.vcf.gz;
    '''
    util.execute_command(command)

    command = fr'''
        bcftools concat \
            --allow-overlaps \
            {germline_vcf} \
            germline_leakage.vcf.gz | \
            bcftools sort --temp-dir ./bcftools.XXXXXX --output {out_fp};
        bcftools index -t {out_fp};
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
