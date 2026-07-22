import gzip
import pathlib
import select
import subprocess
import sys
import textwrap
import logging
from types import SimpleNamespace

import cyvcf2

from .common import constants

# Set up logging
logger = logging.getLogger(__name__)

# TODO(SW): create note that number this assumes location of `<root>/<package>/<file>`
def get_project_root():
    filepath = pathlib.Path(__file__).absolute()
    dirpath = filepath.parent
    project_root = dirpath.parent
    return project_root


def execute_command(command, log_file_path=None):
    # set -e: exit on error, -u: exit on unset variable, -o pipefail: pipeline fails if any command fails
    prepared_command = f'set -euo pipefail; {textwrap.dedent(command)}'
    logger.info("Executing command: %s", command.strip())

    process = subprocess.Popen(
        prepared_command,
        shell=True,
        executable='/bin/bash',
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding='utf-8',
        bufsize=1,  # line buffered
    )

    stdout_lines = []
    stderr_lines = []
    stream_map = {
        process.stdout: (stdout_lines, logger.info),
        process.stderr: (stderr_lines, logger.warning),
    }
    open_streams = set(stream_map)
    log_file = log_file_path.open('a', encoding='utf-8') if log_file_path else None

    try:
        # select multiplexes stdout and stderr in a single thread, preserving arrival order
        # and preventing pipe buffer deadlock without threading races on log_file writes
        while open_streams:
            readable, _, _ = select.select(open_streams, [], [])
            for stream in readable:
                line = stream.readline()
                if line:
                    # Filter out bash libtinfo.so.6 warnings
                    if 'libtinfo.so.6: no version information available' not in line:
                        lines, log_fn = stream_map[stream]
                        log_fn(line.rstrip())
                        lines.append(line)
                        if log_file:
                            log_file.write(line)
                            log_file.flush()
                else:
                    open_streams.discard(stream)
    finally:
        process.wait()
        if log_file:
            log_file.close()

    if process.returncode != 0:
        logger.error("Command failed with return code %d: %s", process.returncode, command.strip())
        raise subprocess.CalledProcessError(
            process.returncode, command,
            output=''.join(stdout_lines),
            stderr=''.join(stderr_lines),
        )

    return SimpleNamespace(
        stdout=''.join(stdout_lines),
        stderr=''.join(stderr_lines),
        returncode=process.returncode,
        pid=process.pid,
        command=command,
    )

def count_vcf_records(fp):
    result = execute_command(f'bcftools view -H {fp} | wc -l')
    return int(result.stdout.strip())


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
    Merge gzipped TSV files into a single gzipped TSV.
    """

    with gzip.open(merged_tsv_fp, 'wt', encoding='utf-8') as merged_tsv:
        for i, tsv_file in enumerate(tsv_files):
            with gzip.open(tsv_file, 'rt', encoding='utf-8') as infile:
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
    merged_vcf_fp = pathlib.Path(merged_vcf_fp)
    merged_unsorted_vcf = merged_vcf_fp.parent / f'{merged_vcf_fp.name}.unsorted.vcf.gz'
    merged_vcf = merged_vcf_fp.parent / f'{merged_vcf_fp.name}.vcf.gz'

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
            assert header_src_description_unquoted == header_dst_entry['Description']
        except AssertionError:
            print(f'Header description mismatch for {header_dst.value}')
            print(f'  src: {header_src_description_unquoted}')
            print(f'  dst: {header_dst_entry["Description"]}')
            sys.exit(1)
