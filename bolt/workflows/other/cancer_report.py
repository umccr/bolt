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
@click.option('--smlv_somatic_counts_process_fp', required=True, type=click.Path(exists=True))

@click.option('--sv_somatic_tsv_fp', required=True, type=click.Path(exists=True))
@click.option('--sv_somatic_vcf_fp', required=True, type=click.Path(exists=True))
@click.option('--cnv_somatic_tsv_fp', required=True, type=click.Path(exists=True))

@click.option('--purple_baf_plot_fp', required=True, type=click.Path(exists=True))

@click.option('--purple_dir', required=True, type=click.Path(exists=True))
@click.option('--virusbreakend_dir', required=True, type=click.Path(exists=True))

@click.option('--cancer_genes_fp', required=True, type=click.Path(exists=True))
@click.option('--oncokb_genes_fp', required=True, type=click.Path(exists=True))

@click.option('--output_dir', required=True, type=click.Path())

def entry(ctx, **kwargs):
    '''Generate UMCCR cancer report\f
    '''

    # Create output directory
    output_dir = pathlib.Path(kwargs['output_dir'])
    output_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

    # Normalise SAGE variants and remove duplicates that arise for MutationalPattern compatibility
    decomposed_snv_vcf = normalise_and_dedup_sage_variants(
        kwargs['smlv_somatic_vcf_fp'],
        kwargs['tumor_name'],
        output_dir,
    )

    # Prepare image directory as required by gpgr
    output_image_dir = prepare_gpgr_image_directory(
        kwargs['purple_dir'],
        kwargs['purple_baf_plot_fp'],
        output_dir,
    )

    # Set other required argument values
    batch_name = f'{kwargs["subject_name"]}_{kwargs["tumor_name"]}'
    output_table_dir = output_dir / 'cancer_report_tables'

    # Run gpgr canrep
    command = fr'''
        gpgr.R canrep \
            \
            --batch_name {batch_name} \
            --tumor_name {kwargs['tumor_name']} \
            \
            --af_global {kwargs['af_global_fp']} \
            --af_keygenes {kwargs['af_keygenes_fp']} \
            \
            --somatic_snv_vcf {decomposed_snv_vcf} \
            --somatic_snv_summary {kwargs['smlv_somatic_counts_process_fp']} \
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
            --virusbreakend_tsv {kwargs['virusbreakend_dir']}/{kwargs['tumor_name']}.virusbreakend.vcf.summary.tsv \
            --virusbreakend_vcf {kwargs['virusbreakend_dir']}/{kwargs['tumor_name']}.virusbreakend.vcf \
            \
            --key_genes {kwargs['cancer_genes_fp']} \
            --oncokb_genes {kwargs['oncokb_genes_fp']} \
            \
            --img_dir {output_image_dir}/ \
            --result_outdir {output_table_dir}/ \
            --out_file {output_dir}/{kwargs["tumor_name"]}.cancer_report.html
    '''
    util.execute_command(command)


def normalise_and_dedup_sage_variants(input_fp, tumor_name, output_dir):
    decomposed_snv_vcf = output_dir / f'{tumor_name}.snvs.normalised.vcf.gz'

    command = fr'''
        bcftools norm \
            --atomize \
            --remove-duplicates \
            --output {decomposed_snv_vcf} \
            {input_fp}
    '''
    util.execute_command(command)

    return decomposed_snv_vcf


def prepare_gpgr_image_directory(purple_dir, purple_baf_plot_fp, output_dir):
    purple_plot_dir = pathlib.Path(purple_dir) / 'plot'
    output_image_dir = output_dir / 'img'

    command = fr'''

        mkdir -p {output_image_dir}/;

        find -L \
            {purple_plot_dir} \
            {purple_baf_plot_fp} \
            -name '*png' -type f -exec cp {{}} {output_image_dir}/ \;;
    '''
    util.execute_command(command)

    return output_image_dir
