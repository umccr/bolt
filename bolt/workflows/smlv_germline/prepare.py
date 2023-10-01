# NOTE(SW): somatic VCF is expected to be unfiltered to capture all germline leakage


import pathlib


import click
import cyvcf2


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

    # NOTE(SW): I am purposefully not selecting variants in the germline gene panel until all other
    # processing has been completed. Slow down is tolerable for the task, and negligible wrt to
    # overall workflow runtime.

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

    # Retain relevant DRAGEN INFO fields
    info_default = [
        'INFO/DP',
        'INFO/FractionInformativeReads',
        'INFO/MQ',
    ]

    # Retain specific INFO fields to enable traceback
    info_retain = [
        util.get_qualified_vcf_annotation(constants.VcfInfo.GERMLINE_LEAKAGE),
        util.get_qualified_vcf_annotation(constants.VcfInfo.PON_COUNT),
        util.get_qualified_vcf_annotation(constants.VcfInfo.GNOMAD_AF),
    ]

    # Keep all format fields except SQ (somatic quality score)
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
    fn_prefix = f'{germline_vcf.replace(".vcf.gz", "")}'

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

    # Identify calls that already exist in the germline call set
    gl_existing_vcf, gl_novel_vcf = get_existing_and_novel('germline_leakage.vcf.gz', germline_vcf)

    # Combine prepared novel germline leakage calls with the actual germline call set
    command = fr'''
        bcftools view -o {gl_novel_vcf}.gz {gl_novel_vcf};
        bcftools index -t {gl_novel_vcf}.gz;
        bcftools concat \
            --allow-overlaps \
            {germline_vcf} \
            {gl_novel_vcf}.gz | \
            bcftools sort --temp-dir ./bcftools.XXXXXX --output {fn_prefix}.leakage.vcf.gz;
    '''
    util.execute_command(command)

    # Annotate germline calls that also had a germline leakage call
    # Read in relevant germline leakage calls
    gl_existing_calls = set()
    for record in cyvcf2.VCF(gl_existing_vcf):
        key = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        assert key not in gl_existing_calls
        gl_existing_calls.add(key)

    # Iterate germline VCF and annotate
    fh_in = cyvcf2.VCF(f'{fn_prefix}.leakage.vcf.gz')
    util.add_vcf_header_entry(fh_in, constants.VcfInfo.GERMLINE_LEAKAGE_CALLED)
    fh_out = cyvcf2.Writer(f'{fn_prefix}.leakage.tagged.vcf.gz', fh_in, 'wz')

    for record in fh_in:
        key = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        if key in gl_existing_calls:
            record.INFO[constants.VcfInfo.GERMLINE_LEAKAGE_CALLED.value] = True
        fh_out.write_record(record)
    fh_out.close()

    # Index final output
    command = fr'''bcftools index -t {fn_prefix}.leakage.tagged.vcf.gz'''
    util.execute_command(command)

    return f'{fn_prefix}.leakage.tagged.vcf.gz'


def get_existing_and_novel(vcf_a, vcf_b):

    # NOTE(SW): overlap with smlv_somatic.rescue:get_sage_existing_and_novel

    dp = pathlib.Path('germline_leakage_isec/')
    util.execute_command(f'bcftools isec -p {dp} {vcf_a} {vcf_b}')

    # Shared records
    calls_existing = dp / '0003.vcf'

    # Records private to vcf_a
    calls_private = dp / '0000.vcf'

    return calls_existing, calls_private


def select_gene_panel_variants(vcf_fp, panel_fp, output_fp):

    command = fr'''
        bcftools view --regions-file {panel_fp} -o {output_fp} {vcf_fp};
        bcftools index -t {output_fp};
    '''
    util.execute_command(command)

    return output_fp
