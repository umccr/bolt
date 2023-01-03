import click


from .workflows import cli as workflow_cli


@click.group(cls=workflow_cli.ClickCustomGroupOrder)
@click.version_option()
def cli():
    pass


def entry():
    for cli_group in workflow_cli.CLI_GROUPS:
        cli.add_command(cli_group)
    cli()


if __name__ == '__main__':
    entry()
