import gzip
import pathlib
import subprocess
import sys
import textwrap


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


def merge_tsv_files(tsv_files, merged_tsv_fp):
    """
    Merges all TSV files into a single TSV.
    """
    logger.info("Merging TSV files...")
    with gzip.open(merged_tsv_fp, 'wt') as merged_tsv:
        for i, tsv_file in enumerate(tsv_files):
            with gzip.open(tsv_file, 'rt') as infile:
                for line_number, line in enumerate(infile):
                    # Skip header except for the first file
                    if i > 0 and line_number == 0:
                        continue
                    merged_tsv.write(line)
    logger.info(f"Merged TSV written to: {merged_tsv_fp}")


def merge_vcf_files(vcf_files, merged_vcf_fp):
    """
    Merges multiple VCF files into a single sorted VCF file using bcftools.

    Parameters:
    - vcf_files: List of paths to VCF files to be merged.
    - merged_vcf_fp: Path to the output merged VCF file (without extension).

    Returns:
    - Path to the sorted merged VCF file.
    """
    logger.info("Merging VCF files...")
    merged_vcf_fp = pathlib.Path(merged_vcf_fp)
    merged_unsorted_vcf = merged_vcf_fp.with_suffix('.unsorted.vcf.gz')
    merged_vcf = merged_vcf_fp.with_suffix('.vcf.gz')

    # Prepare the bcftools merge command arguments
    command_args = [
        'bcftools merge',
        '-m all',
        '-Oz',
        f'-o {merged_unsorted_vcf}',
    ] + [str(vcf_file) for vcf_file in vcf_files]

    # Format the command for readability
    delimiter_padding = ' ' * 10
    delimiter = f' \\\n{delimiter_padding}'
    command_args_str = delimiter.join(command_args)

    command = f'''
    {command_args_str}
    '''

    # Run the bcftools merge command
    logger.info("Running bcftools merge...")
    execute_command(command)
    logger.info(f"Merged VCF written to: {merged_unsorted_vcf}")

    # Sort the merged VCF file
    sort_command_args = [
        'bcftools sort',
        '-Oz',
        f'-o {merged_vcf}',
        f'{merged_unsorted_vcf}'
    ]
    sort_command_args_str = delimiter.join(sort_command_args)
    sort_command = f'''
    {sort_command_args_str}
    '''

    logger.info("Sorting merged VCF file...")
    execute_command(sort_command)
    logger.info(f"Sorted merged VCF written to: {merged_vcf}")

    # Index the sorted merged VCF file
    index_command_args = [
        'bcftools index',
        '-t',
        f'{merged_vcf}'
    ]
    index_command_args_str = delimiter.join(index_command_args)
    index_command = f'''
    {index_command_args_str}
    '''

    logger.info("Indexing sorted merged VCF file...")
    execute_command(index_command)
    logger.info(f"Indexed merged VCF file: {merged_vcf}.tbi")

    # Optionally, remove the unsorted merged VCF file
    if merged_unsorted_vcf.exists():
        merged_unsorted_vcf.unlink()

    return merged_vcf