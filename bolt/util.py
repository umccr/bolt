import gzip
import pathlib
import subprocess
import textwrap
import logging
from types import SimpleNamespace


from .common import pcgr
from .common import constants


logger = logging.getLogger(__name__)


# TODO(SW): create note that number this assumes location of `<root>/<package>/<file>`
def get_project_root():
    filepath = pathlib.Path(__file__).absolute()
    dirpath = filepath.parent
    project_root = dirpath.parent
    return project_root

def execute_command(command, log_file_path=None):
    logger.info("Executing command: %s", command.strip())

    # Open the log file if provided
    log_file = log_file_path.open('a', encoding='utf-8') if log_file_path else None

    # Launch process with combined stdout and stderr streams, and line buffering enabled.
    process = subprocess.Popen(
        command,
        shell=True,
        executable='/bin/bash',
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding='utf-8',
        bufsize=1  # line buffered
    )

    output_lines = []
    # Iterate over each line as it becomes available
    with process.stdout:
        for line in iter(process.stdout.readline, ''):
            if line:
                logger.info(line.strip())
                output_lines.append(line)
                if log_file:
                    log_file.write(line)
                    log_file.flush()  # flush immediately for real-time logging
    process.wait()  # wait for the process to complete

    if log_file:
        log_file.close()

    result = SimpleNamespace(
        stdout=''.join(output_lines),
        returncode=process.returncode,
        pid=process.pid,
        command=command
    )

    return result

def command_prepare(command):
    return f'set -o pipefail; {textwrap.dedent(command)}'

def count_vcf_records(fp):
    result = subprocess.run(f'bcftools view -H {fp} | wc -l',
                            shell=True,
                            executable="/bin/bash",
                            capture_output=True,
                            text=True )
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

def split_vcf(input_vcf, output_dir):
    """
    Splits a VCF file into multiple chunks, each containing up to max_variants variants.
    Each chunk includes the VCF header.
    Ensures no overlapping positions between chunks.
    """
    output_dir = pathlib.Path(output_dir / "vcf_chunks")
    output_dir.mkdir(parents=True, exist_ok=True)

def merge_tsv_files(tsv_files, merged_tsv_fp):
    """
    Merges all TSV files into a single TSV.
    """
    with open(merged_tsv_fp, 'w') as merged_tsv:
        for i, tsv_file in enumerate(tsv_files):
            with open(tsv_file, 'r') as infile:
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