import pathlib


import click


from ... import util


@click.command(name='cancer_report')
@click.pass_context
@click.option('--subject_name', required=True, type=str)
@click.option('--tumor_name', required=True, type=str)

@click.option('--af_global_fp', required=True, type=click.Path(exists=True))
@click.option('--af_keygenes_fp', required=True, type=click.Path(exists=True))

@click.option('--smlv_somatic_vcf_fp', required=True, type=click.Path(exists=True))

@click.option('--sv_somatic_tsv_fp', required=True, type=click.Path(exists=True))
@click.option('--sv_somatic_vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--cnv_somatic_tsv_fp', required=True, type=click.Path(exists=True))

@click.option('--purple_baf_plot_fp', required=True, type=click.Path(exists=True))

@click.option('--purple_dir', required=True, type=click.Path(exists=True))
@click.option('--virusbreakend_dir', required=True, type=click.Path(exists=True))

@click.option('--cancer_genes_fp', required=True, type=click.Path(exists=True))

@click.option('--output_dir', required=True, type=click.Path())
def entry(ctx, **kwargs):
    '''Generate UMCCR cancer report\f
    '''

    batch_name = f'{kwargs["subject_name"]}_{kwargs["tumor_name"]}'

    purple_plot_dir = pathlib.Path(kwargs['purple_dir']) / 'plot'

    output_dir = pathlib.Path(kwargs['output_dir'])
    output_image_dir = output_dir / 'img'
    output_table_dir = output_dir / 'cancer_report_tables'

    command = fr'''
        mkdir -p {output_image_dir}/;

        find -L \
            {purple_plot_dir} \
            {kwargs['purple_baf_plot_fp']} \
            -name '*png' -type f -exec cp {{}} {output_image_dir}/ \;;

        gpgr.R canrep \
            \
            --batch_name {batch_name} \
            --tumor_name {kwargs['tumor_name']} \
            \
            --af_global {kwargs['af_global_fp']} \
            --af_keygenes {kwargs['af_keygenes_fp']} \
            \
            --somatic_snv_vcf {kwargs['smlv_somatic_vcf_fp']} \
            \
            --somatic_sv_tsv {kwargs['sv_somatic_tsv_fp']} \
            --somatic_sv_vcf {kwargs['sv_somatic_vcf_fp']} \
            --purple_som_cnv_ann {kwargs['cnv_somatic_tsv_fp']} \
            \
            --purple_som_gene_cnv {kwargs['purple_dir']}/{kwargs['tumor_name']}.purple.cnv.gene.tsv \
            --purple_som_cnv {kwargs['purple_dir']}/{kwargs['tumor_name']}.purple.cnv.somatic.tsv \
            --purple_purity {kwargs['purple_dir']}/{kwargs['tumor_name']}.purple.purity.tsv \
            --purple_qc {kwargs['purple_dir']}/{kwargs['tumor_name']}.purple.qc \
            --purple_som_snv_vcf {kwargs['purple_dir']}/{kwargs['tumor_name']}.purple.somatic.vcf.gz \
            \
            --virusbreakend_tsv {kwargs['virusbreakend_dir']}/{batch_name}.virusbreakend.vcf.summary.tsv \
            --virusbreakend_vcf {kwargs['virusbreakend_dir']}/{batch_name}.virusbreakend.vcf \
            \
            --key_genes {kwargs['cancer_genes_fp']} \
            \
            --img_dir {output_image_dir}/ \
            --result_outdir {output_table_dir}/ \
            --out_file {output_dir}/report.html
    '''
    util.execute_command(command)
