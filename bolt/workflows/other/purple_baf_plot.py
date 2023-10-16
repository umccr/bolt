import pathlib


import click


from ... import util


@click.command(name='purple_baf_plot')
@click.pass_context

@click.option('--tumor_name', required=True, type=str)

@click.option('--purple_dir', required=True, type=click.Path(exists=True))
@click.option('--circos_conf_fp', required=True, type=click.Path(exists=True))
@click.option('--circos_gaps_fp', required=True, type=click.Path(exists=True))

@click.option('--output_dir', required=True, type=click.Path())

def entry(ctx, **kwargs):
    '''Render PURPLE β-allele frequency circos plot\f
    '''

    # Create output directory
    output_dir = pathlib.Path(kwargs['output_dir'])
    output_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

    # Get circos inputs
    purple_dir = pathlib.Path(kwargs['purple_dir'])
    circos_dir = purple_dir / 'circos'
    assert circos_dir.exists()

    circos_gaps_fp = pathlib.Path(kwargs['circos_gaps_fp'])

    # Generate plot
    command = fr'''
        sed 1>{output_dir}/circos_baf.conf \
          's/SAMPLE/'{kwargs["tumor_name"]}'/' \
          {kwargs["circos_conf_fp"]};

        for ftype in baf cnv map link; do
            src_fp={circos_dir}/{kwargs["tumor_name"]}.${{ftype}}.circos;
            cp ${{src_fp}} {output_dir}/;
        done;

        cp {circos_gaps_fp} {output_dir}/gaps.txt;

        circos \
            -nosvg \
            -conf {output_dir}/circos_baf.conf \
            -outputfile {kwargs["tumor_name"]}.circos_baf.png \
            -outputdir {output_dir};
    '''
    util.execute_command(command)
