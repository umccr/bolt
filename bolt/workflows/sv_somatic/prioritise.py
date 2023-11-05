import pathlib

import click
import cyvcf2

from ... import util
from ...external import prioritize_sv


@click.command(name='prioritise')
@click.pass_context

@click.option('--tumor_name', required=True, type=str)

@click.option('--sv_vcf', required=True, type=click.Path(exists=True))

@click.option('--refdata_known_fusion_pairs', required=True, type=click.Path(exists=True))
@click.option('--refdata_known_fusion_heads', required=True, type=click.Path(exists=True))
@click.option('--refdata_known_fusion_tails', required=True, type=click.Path(exists=True))
@click.option('--refdata_fusioncatcher_pairs', required=True, type=click.Path(exists=True))
@click.option('--refdata_key_genes', required=True, type=click.Path(exists=True))
@click.option('--refdata_key_tsgenes', required=True, type=click.Path(exists=True))

@click.option('--appris_fp', required=True, type=click.Path(exists=True))

@click.option('--output_dir', required=True, type=click.Path())

def entry(ctx, **kwargs):
    '''Prioritise SVs and CNVs\f
    '''

    # Create output directory
    output_dir = pathlib.Path(kwargs['output_dir'])
    output_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

    # Prioritise all variants
    prioritised_fp = output_dir / f'{kwargs["tumor_name"]}.prioritised.vcf'
    prioritize_sv.run(
        kwargs['sv_vcf'],
        kwargs['refdata_known_fusion_pairs'],
        kwargs['refdata_known_fusion_heads'],
        kwargs['refdata_known_fusion_tails'],
        kwargs['refdata_fusioncatcher_pairs'],
        kwargs['refdata_key_genes'],
        kwargs['refdata_key_tsgenes'],
        kwargs['appris_fp'],
        prioritised_fp,
    )

    sv_fp = split_records(prioritised_fp, kwargs['tumor_name'], output_dir, 'sv')
    cnv_fp = split_records(prioritised_fp, kwargs['tumor_name'], output_dir, 'cnv')

    create_sv_tsv(sv_fp, kwargs['tumor_name'], output_dir)
    create_cnv_tsv(cnv_fp, kwargs['tumor_name'], output_dir)


def create_sv_tsv(input_fp, tumor_name, output_dir):

    # NOTE(SW): the cancer report previously used Manta FORMAT/SR and FORMAT/PR as a diagnostic for
    # SV quality. However, GRIDSS handles read counts slightly different and has many measures of
    # various read support/non-support. For now I am using FORMAT/SR and FORMAT/RP.

    header = (
        'chrom',
        'start',
        'svtype',
        'SR_alt',
        'PR_alt',
        'SR_ref',
        'PR_ref',
        'QUAL',
        'tier',
        'annotation',
        'AF_PURPLE',
        'CN_PURPLE',
        'CN_change_PURPLE',
        'PURPLE_status',
        'ID',
        'MATEID',
        'ALT',
    )

    input_fh = cyvcf2.VCF(input_fp, samples=[tumor_name])
    output_fh = (output_dir / f'{tumor_name}.sv.prioritised.tsv').open('w')

    print(*header, sep='\t', file=output_fh)

    for record in input_fh:

        if record.FILTER and 'INFERRED' in record.FILTER:
            purple_status = 'INFERRED'
        elif record.INFO.get('RECOVERED') is not None:
            purple_status = 'RECOVERED'
        else:
            purple_status = ''

        data = (
            record.CHROM.replace('chr', ''),
            record.POS,
            record.INFO.get('EVENTTYPE', ''),
            parse_read_support_field(record, 'SR'),
            parse_read_support_field(record, 'RP'),
            parse_read_support_field(record, 'REF'),
            parse_read_support_field(record, 'REFPAIR'),
            record.QUAL,
            record.INFO.get('SV_TOP_TIER', 4),
            record.INFO['SIMPLE_ANN'],
            parse_info_field(record, 'PURPLE_AF'),
            parse_info_field(record, 'PURPLE_CN'),
            parse_info_field(record, 'PURPLE_CN_CHANGE'),
            purple_status,
            record.ID,
            parse_info_field(record, 'MATEID'),
            record.ALT[0],
        )

        print(*data, sep='\t', file=output_fh)


def create_cnv_tsv(input_fp, tumor_name, output_dir):

    input_fh = cyvcf2.VCF(input_fp, samples=[tumor_name])
    output_fh = (output_dir / f'{tumor_name}.cnv.prioritised.tsv').open('w')

    purple_fields = (
        'baf',
        'bafCount',
        'copyNumber',
        'depthWindowCount',
        'gcContent',
        'majorAlleleCopyNumber',
        'method',
        'minorAlleleCopyNumber',
        'segmentEndSupport',
        'segmentStartSupport',
    )

    header = (
        'chromosome',
        'start',
        'end',
        'svtype',
        *purple_fields,
        'sv_top_tier',
        'simple_ann',
    )

    print(*header, sep='\t', file=output_fh)

    for record in input_fh:

        purple_fields_data = list()
        for purple_field in purple_fields:
            purple_fields_data.append(record.INFO[f'PURPLE_{purple_field}'])

        data = (
            record.CHROM,
            record.POS,
            record.INFO['END'],

            record.ALT[0].strip('<>'),

            *purple_fields_data,

            record.INFO.get('SV_TOP_TIER', 4),
            record.INFO['SIMPLE_ANN'],

        )

        print(*data, sep='\t', file=output_fh)


def parse_info_field(record, field_name):
    data = record.INFO.get(field_name)
    if data is None:
        return str()
    elif isinstance(data, (list, tuple, set)):
        # Lazily assuming flat list/tuple/set
        return ','.join(str(e) for e in data)
    else:
        return data


def parse_read_support_field(record, field_name):
    if (data := record.format(field_name)):
        return data[0][0]
    else:
        return str()


def split_records(input_fp, tumor_name, output_dir, variant_type):
    if variant_type == 'sv':
        source = 'sv_gridss'
    elif variant_type == 'cnv':
        source = 'cnv_purple'
    else:
        assert False

    output_fp = output_dir / f'{tumor_name}.{variant_type}.prioritised.vcf.gz'
    command = fr'''
        bcftools view -i 'SOURCE="{source}"' {input_fp} | \
            bcftools annotate -x 'INFO/SOURCE' -o {output_fp}
    '''
    util.execute_command(command)

    return output_fp
