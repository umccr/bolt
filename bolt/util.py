import pathlib
import subprocess
import sys
import textwrap
import cyvcf2


from .common import constants


# TODO(SW): create note that number this assumes location of `<root>/<package>/<file>`
def get_project_root():
    filepath = pathlib.Path(__file__).absolute()
    dirpath = filepath.parent
    project_root = dirpath.parent
    return project_root


def execute_command(command):
    command_prepared = command_prepare(command)

    print(command_prepared)

    process = subprocess.run(
        command_prepared,
        shell=True,
        executable='/bin/bash',
        capture_output=True,
        encoding='utf-8',
    )

    if process.returncode != 0:
        print(process)
        print(process.stderr)
        sys.exit(1)

    return process


def command_prepare(command):
    return f'set -o pipefail; {textwrap.dedent(command)}'


#def count_vcf_records(fp, exclude_args=None):
#    args = list()
#    if exclude_args:
#        args.append(f'-e \'{exclude_args}\'')
#
#    args_str = ' '.join(args)
#    command = f'bcftools view -H {args_str} {fp} | wc -l'
#
#    result = execute_command(command)
#    return int(result.stdout)


def count_vcf_records(fp):
    result = execute_command(f'bcftools view -H {fp} | wc -l')
    return int(result.stdout)


def add_vcf_header_entry(fh, anno_enum):
    header_entry = get_vcf_header_entry(anno_enum)
    if anno_enum in constants.VcfFilter:
        fh.add_filter_to_header(header_entry)
    elif anno_enum in constants.VcfInfo:
        fh.add_info_to_header(header_entry)
    elif anno_enum in constants.VcfFormat:
        fh.add_format_to_header(header_entry)
    else:
        assert False


def get_vcf_header_entry(anno_enum):
    header_entry = constants.VCF_HEADER_ENTRIES[anno_enum]
    return {'ID': anno_enum.value, **header_entry}


def get_vcf_header_line(anno_enum):
    header_entry = get_vcf_header_entry(anno_enum)

    if anno_enum in constants.VcfFilter:
        return f'##{anno_enum.namespace}=<ID={anno_enum.value},Description=\"{header_entry["Description"]}\">'
    elif anno_enum in constants.VcfInfo or anno_enum in constants.VcfFormat:
        return (
            f'##{anno_enum.namespace}=<'
            f'ID={anno_enum.value},'
            f'Number={header_entry["Number"]},'
            f'Type={header_entry["Type"]},'
            f'Description=\"{header_entry["Description"]}\">'
        )
    else:
        assert False

def get_qualified_vcf_annotation(anno_enum):
    assert anno_enum in constants.VcfInfo or anno_enum in constants.VcfFormat
    return f'{anno_enum.namespace}/{anno_enum.value}'


#def add_vcf_filter(record, filter_enum):
#    existing_filters = [e for e in record.FILTERS if e != 'PASS']
#    assert filter_enum.value not in existing_filters
#    return ';'.join([*existing_filters, filter_enum.value])

def check_annotation_headers(info_field_map, vcf_fp):
    # Ensure header descriptions from source INFO annotations match those defined here for the
    # output file; force manual inspection where they do not match
    vcf_fh = cyvcf2.VCF(vcf_fp)
    for header_dst, header_src in info_field_map.items():
        # Skip header lines that do not have an equivalent entry in the VCF
        try:
            header_src_entry = vcf_fh.get_header_type(header_src)
        except KeyError:
            continue

        header_dst_entry = get_vcf_header_entry(header_dst)
        # Remove leading and trailing quotes from source
        header_src_description_unquoted = header_src_entry['Description'].strip('"')
        try:
            assert  header_src_description_unquoted == header_dst_entry['Description']
        except AssertionError:
            print(f'Header description mismatch for {header_dst.value}')
            print(f'  src: {header_src_description_unquoted}')
            print(f'  dst: {header_dst_entry["Description"]}')
            sys.exit(1)