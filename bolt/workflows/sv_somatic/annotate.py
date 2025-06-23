import csv
import pathlib


import click
import cyvcf2
import pysam
import logging


from ... import util

logger = logging.getLogger(__name__)


@click.command(name='annotate')
@click.pass_context

@click.option('--tumor_name', required=True, type=str)

@click.option('--sv_fp', required=True, type=click.Path(exists=True))
@click.option('--cnv_fp', required=True, type=click.Path(exists=True))

@click.option('--reference_fasta_fp', required=True, type=click.Path(exists=True))
@click.option('--snpeff_database_dir', required=True, type=click.Path(exists=True))

@click.option('--output_dir', required=True, type=click.Path())

def entry(ctx, **kwargs):
    '''Annotate SVs and CNVs with functional information\f
    '''

    # Create output directory
    output_dir = pathlib.Path(kwargs['output_dir'])
    output_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

    # Create compiled VCF containing variants for annotation
    variants_fp = compile_variants(
        kwargs['sv_fp'],
        kwargs['cnv_fp'],
        kwargs['tumor_name'],
        kwargs['reference_fasta_fp'],
        output_dir,
    )

    # Annotate variants
    annotate_variants(variants_fp, kwargs['tumor_name'], kwargs['snpeff_database_dir'], output_dir)


def compile_variants(sv_fp, cnv_fp, tumor_name, ref_fp, output_dir):

    output_fp = output_dir / f'{tumor_name}.sv_cnv.vcf.gz'
    output_sorted_fp = output_dir / f'{tumor_name}.sv_cnv.sorted.vcf.gz'

    # Prepare VCF writer; use SV VCF as template
    sv_fh = cyvcf2.VCF(sv_fp)
    sv_fh.add_to_header('##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of variant">')
    sv_fh.add_to_header('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">')

    # TODO(SW): refactor when not under time pressure

    purple_fields = (
        ('Float', 'baf'),
        ('Integer', 'bafCount'),
        ('Float', 'copyNumber'),
        ('Integer', 'depthWindowCount'),
        ('Float', 'gcContent'),
        ('Float', 'majorAlleleCopyNumber'),
        ('String', 'method'),
        ('Float', 'minorAlleleCopyNumber'),
        ('String', 'segmentEndSupport'),
        ('String', 'segmentStartSupport'),
    )

    for field_type, purple_field in purple_fields:
        sv_fh.add_to_header(f'##INFO=<ID=PURPLE_{purple_field},Number=1,Type={field_type},Description="PURPLE {purple_field}">')

    output_fh = cyvcf2.Writer(output_fp, sv_fh, 'wz')

    # SVs
    for record in sv_fh:
        record.INFO['SOURCE'] = 'sv_esvee'
        output_fh.write_record(record)

    # CNVs
    # Collect information from TSV
    cnvs = list()
    with pathlib.Path(cnv_fp).open('r') as fh:
        for record in csv.DictReader(fh, delimiter='\t'):

            cn = float(record['copyNumber'])

            if cn >= 2.05:
                svtype = 'DUP'
            elif cn  <= 1.95:
                svtype = 'DEL'
            else:
                continue

            info_fields = dict()
            for field_type, purple_field in purple_fields:
                field_name = f'PURPLE_{purple_field}'
                assert field_name not in info_fields
                info_fields[field_name] = record[purple_field]

            info_fields['END'] = record['end']
            info_fields['SOURCE'] = 'cnv_purple'
            info_fields['SVTYPE'] = svtype

            record_data = {
                'chrom': record['chromosome'],
                'pos': record['start'],
                'alt': f'<{svtype}>',
                'info': info_fields,
            }

            cnvs.append(record_data)

    # Prepare as VCF record and write
    ref_fh = pysam.FastaFile(ref_fp)
    for i, record_data in enumerate(cnvs):
        ref_pos = int(record_data['pos'])
        ref = ref_fh.fetch(record_data['chrom'], ref_pos-1, ref_pos)
        record_variant_comps = (
            record_data['chrom'],
            record_data['pos'],
            f'cnv_purple_{i}',
            ref,
            record_data['alt'],
            '.',
        )
        info_str = ';'.join(f'{k}={v}' for k, v in record_data['info'].items())

        # NOTE(SW): being explicit here for format and sample columns
        format_column = '.'
        sample_a_column = '.'
        sample_b_column = '.'

        record_str_new = '\t'.join([
            *record_variant_comps,
            'PASS',
            info_str,
            format_column,
            sample_a_column,
            sample_b_column,
        ])

        record_new = output_fh.variant_from_string(record_str_new)
        output_fh.write_record(record_new)

    output_fh.close()

    util.execute_command(fr'''
        bcftools sort --temp-dir ./bcftools.XXXXXX --output {output_sorted_fp} {output_fp}
    ''')

    return output_sorted_fp


def annotate_variants(input_fp, tumor_name, snpeff_database, output_dir):

    output_fp = output_dir / f'{tumor_name}.annotated.vcf.gz'

    # NOTE(SW): snpEff installed via Conda requires an aboslute path for the database
    snpeff_database_abs = pathlib.Path(snpeff_database).resolve()

    snpeff_temp_dir = output_dir / 'snpeff_temp/'

    command = fr'''
        snpEff -Xms750m -Xmx4g -Djava.io.tmpdir={snpeff_temp_dir} \
            -dataDir {snpeff_database_abs} \
            -noStats \
            GRCh38.105 \
            {input_fp} | \
            sed -E -e 's/CHR([0-9]{{1,2}}|[XY])/chr\1/g' | \
            bcftools view -o {output_fp};

        bcftools index -t {output_fp};
    '''

    util.execute_command(command)
