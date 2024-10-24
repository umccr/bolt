import pathlib
import subprocess
import textwrap
import cyvcf2
import logging

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
    """
    Executes a shell command.

    Parameters:
    - command: Command to be executed as a formatted string.
    - log_file_path: Optional path to a log file to capture the command output.

    Raises:
    - SystemExit if the command fails.
    """
    logger.info('Executing command:')
    logger.info(command.strip())

    try:
        if log_file_path:
            with open(log_file_path, 'w') as log_file:
                process = subprocess.run(
                    command,
                    shell=True,
                    executable='/bin/bash',
                    check=True,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    text=True,
                    encoding='utf-8'
                )
        else:
            process = subprocess.run(
                command,
                shell=True,
                executable='/bin/bash',
                check=True,
                capture_output=True,
                text=True,
                encoding='utf-8'
            )
            if process.stdout:
                logger.info(process.stdout)
            if process.stderr:
                logger.error(process.stderr)
            return(process)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with error: {e}")
        raise SystemExit(e)

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

def split_vcf(input_vcf, output_dir):
    """
    Splits a VCF file into multiple chunks, each containing up to max_variants variants.
    Each chunk includes the VCF header.
    Ensures no overlapping positions between chunks.
    """
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    chunk_files = []
    chunk_number = 1
    variant_count = 0
    base_filename = pathlib.Path(input_vcf).stem
    chunk_filename = output_dir / f"{base_filename}_chunk{chunk_number}.vcf"
    chunk_files.append(chunk_filename)

    # Open the input VCF using cyvcf2
    vcf_in = cyvcf2.VCF(input_vcf)
    # Create a new VCF file for the first chunk
    vcf_out = cyvcf2.Writer(str(chunk_filename), vcf_in)

    last_position = None

    for record in vcf_in:
        current_position = record.POS
        # Check if we need to start a new chunk
        if variant_count >= 15000 and (last_position is None or current_position != last_position):
            # Close the current chunk file and start a new one
            vcf_out.close()
            chunk_number += 1
            chunk_filename = output_dir / f"{base_filename}_chunk{chunk_number}.vcf"
            chunk_files.append(chunk_filename)
            vcf_out = cyvcf2.Writer(str(chunk_filename), vcf_in)
            variant_count = 0

        # Write the record to the current chunk
        vcf_out.write_record(record)
        variant_count += 1
        last_position = current_position

    # Close the last chunk file
    vcf_out.close()
    vcf_in.close()

    logger.info(f"VCF file split into {len(chunk_files)} chunks.")

    return chunk_files

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
        f'bcftools merge',
        f'-m all',
        f'-Oz',
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
        f'bcftools sort',
        f'-Oz',
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
        f'bcftools index',
        f'-t',
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