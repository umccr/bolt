# NOTE(SW): this specifically takes the FILTERED SAGE VCF /and/ UNFILTER DRAGEN VCF, both with main
# assembly calls only i.e. chr1-22, chrXYM; expected to be done externally for now


import pathlib


import click
import cyvcf2


from ... import util
from ...common import constants


@click.command(name='rescue')
@click.pass_context
@click.option('--vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--sage_vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--hotspots_fp', required=True, type=click.Path(exists=True))
@click.option('--output_fp', required=True, type=click.Path())
def entry(ctx, **kwargs):
    '''Rescue variants using SAGE calls\f

    1. Identify existing and novel SAGE calls in hotspots
    2. Rescue/annotate passing variants otherwise update FILTER accordingly with existing SAGE calls
    3. Add new SAGE variants to input VCF
    '''
    # Identify variants that fall within hotspots

    # Select PASS SAGE variants in hotspots and then split into existing and novel calls
    sage_pass_vcf_fp = select_sage_pass_hotspot(kwargs['sage_vcf_fp'], kwargs['hotspots_fp'])
    sage_existing_vcf_fp, sage_novel_vcf_fp = get_sage_existing_and_novel(sage_pass_vcf_fp, kwargs['vcf_fp'])

    # Perform the following for variants re-called by SAGE:
    #  - SAGE FILTER=PASS, input FILTER=PASS:  set INFO/SAGE_HOTSPOT [annotate]
    #  - SAGE FILTER=PASS, input FILTER!=PASS: set INFO/SAGE_HOTSPOT, INFO/SAGE_RESCUE, FILTER=PASS [rescue]
    #  - SAGE FILTER!=PASS:                    append SAGE_lowconf to FILTER [exclude]
    # Additionally transfer SAGE FORMAT/AD, FORMAT/AF, FORMAT/DP, FORMAT/SB with the 'SAGE_' prefix for all
    # re-called variants regardless of FILTER
    vcf_anno_vcf_fp = annotate_existing_sage_calls(sage_existing_vcf_fp, kwargs['vcf_fp'])

    # Combine annotated, existing calls with SAGE novel calls
    combined_vcf_fp = combine_sage_novel(sage_novel_vcf_fp, vcf_anno_vcf_fp, kwargs['output_fp'])


def select_sage_pass_hotspot(in_fp, hotspots_fp):
    in_fn = pathlib.Path(in_fp).name
    out_fp = f'{in_fn.replace(".vcf.gz", ".hotspot_pass.vcf.gz")}'
    command = fr'''
        bcftools view -f PASS,. -R {hotspots_fp} -o {out_fp} {in_fp} && \
            bcftools index -t {out_fp}
    '''
    util.execute_command(command)
    return out_fp


def get_sage_existing_and_novel(sage_fp, vcf_fp):
    dp = pathlib.Path('sage_existing_and_novel_isec/')
    util.execute_command(f'bcftools isec -p {dp} {vcf_fp} {sage_fp}')
    sage_existing_fp = dp / '0003.vcf'  # records from sage_fp shared by vcf_fp
    sage_novel_fp = dp / '0001.vcf'  # records private to sage_fp
    return sage_existing_fp, sage_novel_fp


def annotate_existing_sage_calls(sage_fp, vcf_in_fp):
    # Read all SAGE calls into memory
    sage_calls = dict()
    for record in cyvcf2.VCF(sage_fp):
        key = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        assert key not in sage_calls
        sage_calls[key] = record

    # Get input file handle
    vcf_in_fh = cyvcf2.VCF(vcf_in_fp)

    # Add header entries so that they are included in the output file via templating done below
    util.add_vcf_header_entry(vcf_in_fh, constants.VcfFilter.SAGE_LOWCONF)

    util.add_vcf_header_entry(vcf_in_fh, constants.VcfInfo.SAGE_HOTSPOT)
    util.add_vcf_header_entry(vcf_in_fh, constants.VcfInfo.SAGE_RESCUE)

    # TODO(SW): check that defined header descriptions match those in the SAGE fp; collect as list
    # here and iterate to check and then add to vcf_in_fh header also in another loop

    util.add_vcf_header_entry(vcf_in_fh, constants.VcfFormat.SAGE_AD)
    util.add_vcf_header_entry(vcf_in_fh, constants.VcfFormat.SAGE_AF)
    util.add_vcf_header_entry(vcf_in_fh, constants.VcfFormat.SAGE_DP)
    util.add_vcf_header_entry(vcf_in_fh, constants.VcfFormat.SAGE_SB)

    # Open output file and use header from input file
    vcf_in_fn = pathlib.Path(vcf_in_fp).name
    vcf_out_fp = f'{vcf_in_fn.replace(".vcf.gz", ".anno.vcf.gz")}'
    vcf_out_fh = cyvcf2.Writer(vcf_out_fp, vcf_in_fh, 'wz')

    # Annotate SAGE recalls
    for record in vcf_in_fh:
        # Write unmodified records that were not re-called by SAGE
        key = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        if not (sage_record := sage_calls.get(key)):
            vcf_out_fh.write_record(record)
            continue

        # Perform the following for records re-called by SAGE:
        #  - SAGE FILTER=PASS, input FILTER=PASS:  set INFO/SAGE_HOTSPOT [annotate]
        #  - SAGE FILTER=PASS, input FILTER!=PASS: set INFO/SAGE_HOTSPOT, INFO/SAGE_RESCUE, FILTER=PASS [rescue]
        #  - SAGE FILTER!=PASS:                    append SAGE_lowconf to FILTER [exclude]
        if sage_record.FILTER is None:
            record.INFO[constants.VcfInfo.SAGE_HOTSPOT.value] = True
            if record.FILTER is not None:
                record.FILTER = 'PASS'
                record.INFO[constants.VcfInfo.SAGE_RESCUE.value] = True
        else:
            record.FILTER = ';'.join([*record.FILTERS, constants.VcfInfo.SAGE_LOWCONF.value])

        # NOTE(SW): previously in Umccrise, the FORMAT/AD and FORMAT/DP from SAGE calls were just used
        # to overwrite bcbio/DRAGEN equivalents but other data such as FORMAT/AF were not updated.
        # Here I instead retain these (and other) data by moving them into a new variables with the
        # form 'FORMAT/SAGE_<name>'.

        # Transfer some SAGE data to the output VCF
        record.set_format(constants.VcfFormat.SAGE_AD.value, sage_record.format('AD'))
        record.set_format(constants.VcfFormat.SAGE_AF.value, sage_record.format('AF'))
        record.set_format(constants.VcfFormat.SAGE_DP.value, sage_record.format('DP'))
        record.set_format(constants.VcfFormat.SAGE_SB.value, sage_record.format('SB'))

        vcf_out_fh.write_record(record)

    # Explicitly close to flush buffer then index output file
    vcf_out_fh.close()
    util.execute_command(f'bcftools index -t {vcf_out_fp}')

    return vcf_out_fp


def combine_sage_novel(sage_novel_fp, anno_fp, output_fp):

    # TODO(SW): check whether we need to reorder samples

    # We must rename some FORMAT fields to avoid namespace collision between the SAGE and DRAGEN
    # VCF; the input for the concat operation must also be bgzip compressed and have an index
    sage_novel_prep_fp = prepare_sage_novel(sage_novel_fp)

    # Add novel SAGE calls to annotated, existing calls
    command = fr'''
        bcftools concat -a {anno_fp} {sage_novel_prep_fp} | \
          bcftools sort -T $(pwd -P)/bcftools.XXXXXX -m 1G -o {output_fp}
    '''
    util.execute_command(command)


def prepare_sage_novel(in_fp):
    # Annotations to rename
    # NOTE(SW): here I only rename FORMAT/SB since SAGE measures this differently while others
    # should be sufficiently interchangeable
    annotations_rename = (
        #('FORMAT/AD', constants.VcfFormat.SAGE_AD.value),
        #('FORMAT/AF', constants.VcfFormat.SAGE_AF.value),
        #('FORMAT/DP', constants.VcfFormat.SAGE_DP.value),
        ('FORMAT/SB', constants.VcfFormat.SAGE_SB.value),
    )
    # Annotations to retain, all other SAGE annotations will be excluded
    annotations_retain_info = (
        util.get_qualified_vcf_annotation(constants.VcfInfo.SAGE_HOTSPOT),
        util.get_qualified_vcf_annotation(constants.VcfInfo.SAGE_NOVEL),
    )
    annotations_retain_format = (
        util.get_qualified_vcf_annotation(constants.VcfFormat.SAGE_SB),
        'FORMAT/AD',
        'FORMAT/AF',
        'FORMAT/DP',
    )

    # Write file rename annotations
    with open('rename_annotations.tsv', 'w') as fh:
        for ann, ann_new in annotations_rename:
            print(ann, ann_new, sep='\t', file=fh)

    # Get and write required header entries to file
    with open('header_entries.tsv', 'w') as fh:
        print(util.get_vcf_header_line(constants.VcfInfo.SAGE_HOTSPOT), file=fh)
        print(util.get_vcf_header_line(constants.VcfInfo.SAGE_NOVEL), file=fh)

    # NOTE(SW): the fill-tags BCFtools plugin would have been useful here instead of iterating with
    # awk but unforunately it doesn't seem to support anything other than integer and float INFO
    # types:
    #   * https://github.com/samtools/bcftools/blob/1.16/plugins/fill-tags.c#L423-L424

    sage_info_annotations = [constants.VcfInfo.SAGE_HOTSPOT.value, constants.VcfInfo.SAGE_NOVEL.value]
    sage_info_annotations_str = ';'.join(sage_info_annotations)

    # Build exclude spec: '^{comma delimited INFO},^{comma delimited FORMAT}'
    annotations_retain_info_str = '^' + ','.join(annotations_retain_info)
    annotations_retain_format_str = '^' + ','.join(annotations_retain_format)
    annotations_retain_str = f'{annotations_retain_info_str},{annotations_retain_format_str}'

    # NOTE(SW): the annotate commands must be split so that header lines for INFO/SAGE_HOTSPOT and
    # INFO/SAGE_NOVEL exist before applying the remove operation; doing in a single annotate
    # command fails for this reason

    # I prepare the novel SAGE calls as follows:
    #   1. rename existing SAGE annotations and add required header entries new SAGE annotations
    #   2. add new SAGE annotations to each call
    #   3. retrain only target annotations
    #   4. index output VCF
    out_fp = 'sage.novel.vcf.gz'
    command = fr'''
        bcftools annotate \
            --rename-annots rename_annotations.tsv \
            --header-lines header_entries.tsv \
            {in_fp} | \
            awk '
                BEGIN {{ OFS="\t" }}
                $1 ~ /^#/ {{ print }}
                $1 !~ /^#/ {{
                    if (length($8) > 1 ) {{
                        $8=$8";{sage_info_annotations_str}";
                    }} else {{
                        $8="{sage_info_annotations_str}";
                    }}
                    print
                }}
            ' | \
            bcftools annotate --remove '{annotations_retain_str}' | \
            bcftools view -o {out_fp} && \
            bcftools index -t {out_fp}
    '''
    util.execute_command(command)

    return out_fp
