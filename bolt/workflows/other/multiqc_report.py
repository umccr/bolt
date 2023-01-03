import pathlib
import yaml


import click


from ... import util


@click.command(name='multiqc_report')
@click.pass_context
@click.option('--tumor_name', required=True, type=str)
@click.option('--normal_name', required=True, type=str)

@click.option('--input_dir', required=True, type=click.Path(exists=True))

@click.option('--output_dir', required=True, type=click.Path())
def entry(ctx, **kwargs):
    '''Generate MultiQC report\f
    '''

    # Create MultiQC config
    multiqc_conf = {
        'umccr': {
            'tumor_name': kwargs['tumor_name'],
            'normal_name': kwargs['normal_name'],
        },
    }

    multiqc_conf_fp = pathlib.Path('multiqc_config.yaml')
    with multiqc_conf_fp.open('w') as fh:
        yaml.dump(multiqc_conf, fh, default_flow_style=False)

    # Run MultiQC
    command = fr'''
        multiqc \
            --config {multiqc_conf_fp} \
            --outdir {kwargs["output_dir"]}/ \
            {kwargs["input_dir"]}/
    '''
    util.execute_command(command)
