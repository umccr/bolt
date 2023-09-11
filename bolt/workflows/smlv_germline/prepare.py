# NOTE(SW): somatic VCF is expected to be unfiltered to capture all germline leakage


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

    # NOTE(SW): the QUAL value for somatic but not germline variants is empty. DRAGEN uses INFO/SQ
    # in somatic variant records to indicate quality score, though this is scaled differently to
    # germline variant QUAL scores and so is not subsituted. Hence, germline leakage variants have
    # an empty QUAL field.

    # Retain all DRAGEN INFO fields
    info_default = [
        'INFO/DP',
        'INFO/FractionInformativeReads',
        'INFO/MQ',
    ]

    # These INFO fields are retained to enable traceback
    info_retain = [
        util.get_qualified_vcf_annotation(constants.VcfInfo.GERMLINE_LEAKAGE),
        util.get_qualified_vcf_annotation(constants.VcfInfo.PON_COUNT),
        util.get_qualified_vcf_annotation(constants.VcfInfo.GNOMAD_AF),
    ]

    # Keep all but FORMAT/SQ
    format_default = [
        'FORMAT/AD',
        'FORMAT/AF',
        'FORMAT/DP',
        'FORMAT/F1R2',
        'FORMAT/F2R1',
        'FORMAT/GT',
        'FORMAT/MB',
        'FORMAT/SB',
    ]

    exclude_str = ',^'.join([*info_default, *info_retain, *format_default])
    out_fp = f'{germline_vcf.replace(".vcf.gz", ".include_leakage.vcf.gz")}'

    # Process germline leakage
    # 1. select tumor call data
    # 2. rename as normal sample
    # 3. clear FILTERs
    # 4. remote unwanted INFO and FORMAT fields

    leakage_info_str = util.get_qualified_vcf_annotation(constants.VcfInfo.GERMLINE_LEAKAGE)
    command = fr'''
        bcftools view -s {tumor_name} -i '{leakage_info_str}=1' {somatic_vcf} | \
            sed '/^#CHROM/ s/{tumor_name}/{normal_name}/' | \
            awk 'BEGIN {{ OFS="\t" }}; /^[^#]/ {{ $7 = "PASS"; print; next }}; {{ print }}' | \
            bcftools annotate -x "{exclude_str}" -o germline_leakage.vcf.gz;
        bcftools index -t germline_leakage.vcf.gz;
    '''
    util.execute_command(command)

    # Combine prepared germline leakage with full normal calls
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
