import importlib
import pathlib


import click


from .. import util


# Desired sort oder of CLI groups; any group not present in this list will be placed equal last
CLI_ORDERING = [
    'smlv_somatic',
    'smlv_germline',
    'sv_somatic',
    'other',
]


def get_command_modules(package_dir_abs):
    root_path = util.get_project_root()
    package_dir_rel = package_dir_abs.relative_to(root_path)
    package_name = str(package_dir_rel).replace('/', '.')
    for command_fp in package_dir_abs.rglob('*.py'):

        if command_fp.name == '__init__.py':
            continue

        module_dir = command_fp.relative_to(package_dir_abs).parent
        module_name_base = str(module_dir).replace('/', '.')
        if module_name_base != '.':
            module_name_prefix = f'.{module_name_base}'
        else:
            module_name_prefix = ''
        module_name_suffix = command_fp.stem
        module_name = f'{module_name_prefix}.{module_name_suffix}'
        yield importlib.import_module(module_name, package_name)


def create_clis():
    clis = list()

    workflow_dir = pathlib.Path(__file__).parent
    workflow_fps = [fp for fp in workflow_dir.glob('*') if fp.name != 'cli.py']

    for fp in workflow_dir.glob('*'):
        if not fp.is_dir() or fp.name == '__pycache__':
            continue

        workflow_name = fp.name

        @click.group(name=workflow_name)
        def cli_group():
            pass

        for module in get_command_modules(fp):
            cli_group.add_command(module.entry)
        clis.append(cli_group)

    return clis


class ClickCustomGroupOrder(click.Group):

    def __init__(self, *args, name=None, commands=None, **kwargs):
        super().__init__(*args, **kwargs)
        if commands is None:
            commands = dict()
        self.commands = commands

    def list_commands(self, ctx):
        # Apply manual sorting to CLI groups, mostly so that 'other' appears last rather than first
        return {n: self.commands[n] for n in sorted(self.commands, key=self.command_sort)}

    def command_sort(self, name):
        # If we have defined sort order for a given CLI group, set position. Otherwise set to equal
        # last, presumably stable.
        return CLI_ORDERING.index(name) if name in CLI_ORDERING else len(CLI_ORDERING)


CLI_GROUPS = create_clis()
