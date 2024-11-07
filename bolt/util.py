import pathlib
import subprocess
import textwrap
import logging


import cyvcf2


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
    """
    Executes a shell command, streaming output to console, logging, and capturing errors.

    Parameters:
    - command: Command to be executed as a formatted string.
    - log_file_path: Optional path to a log file to capture the command output.

    Returns:
    - (stdout, stderr): A tuple containing the complete stdout and stderr of the command.

    Raises:
    - subprocess.CalledProcessError if the command exits with a non-zero status.
    """
    logger.info(command.strip())

    # Open the log file if a path is specified
    log_file = None
    if log_file_path:
        log_file = log_file_path.open('a', encoding='utf-8')

    # Start the process
    with subprocess.Popen(command, shell=True, executable='/bin/bash',
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, encoding='utf-8') as process:

        # To collect the complete stdout and stderr
        full_stdout = []
        full_stderr = []

        # Stream stdout and stderr in real-time
        while True:
            stdout_line = process.stdout.readline()
            stderr_line = process.stderr.readline()

            # Handle stdout
            if stdout_line:
                full_stdout.append(stdout_line)
                if log_file:
                    log_file.write(stdout_line)
                    log_file.flush()  # Ensure immediate write
                logger.info(stdout_line.strip())

            # Handle stderr
            if stderr_line:
                full_stderr.append(stderr_line)
                if log_file:
                    log_file.write(stderr_line)
                    log_file.flush()
                logger.error(stderr_line.strip())

            # If the process has finished and there are no more lines, break
            if process.poll() is not None and not stdout_line and not stderr_line:
                break

        # Communicate to ensure we capture everything left in the buffers
        remaining_stdout, remaining_stderr = process.communicate()

        # Handle remaining stdout
        if remaining_stdout:
            full_stdout.append(remaining_stdout)
            if log_file:
                log_file.write(remaining_stdout)
                log_file.flush()
            logger.info(remaining_stdout.strip())

        # Handle remaining stderr
        if remaining_stderr:
            full_stderr.append(remaining_stderr)
            if log_file:
                log_file.write(remaining_stderr)
                log_file.flush()
            logger.error(remaining_stderr.strip())

    # Close the log file if it was opened
    if log_file:
        log_file.close()

    # Combine all the lines into final stdout and stderr strings
    stdout_str = ''.join(full_stdout)
    stderr_str = ''.join(full_stderr)
    
    return process


def command_prepare(command):
    return f'set -o pipefail; {textwrap.dedent(command)}'

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