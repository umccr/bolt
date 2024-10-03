import gzip
import pathlib
import subprocess
import sys
import os
import textwrap
import pysam


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

def execute_command(command):
    """
    Executes a shell command.

    Parameters:
    - command: Command to be executed as a formatted string.

    Raises:
    - SystemExit if the command fails.
    """
    print('Executing command:')
    print(command.strip())

    try:
        process = subprocess.run(
            command,
            shell=True,
            executable='/bin/bash',
            check=True,
            capture_output=True,
            text=True
        )
        if process.stdout:
            print(process.stdout)
        if process.stderr:
            print(process.stderr)
        return process
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        if e.stdout:
            print(f"Standard output:\n{e.stdout}")
        if e.stderr:
            print(f"Standard error:\n{e.stderr}")
        sys.exit(1)

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

#def add_vcf_filter(record, filter_enum):
#    existing_filters = [e for e in record.FILTERS if e != 'PASS']
#    assert filter_enum.value not in existing_filters
#    return ';'.join([*existing_filters, filter_enum.value])

def split_vcf(input_vcf, output_dir):
    """
    Splits a VCF file into multiple chunks, each containing up to max_variants variants.
    Each chunk includes the VCF header.
    """

    chunk_files = []
    chunk_number = 1
    variant_count = 0
    base_filename = input_vcf.stem
    chunk_filename = output_dir / f"{base_filename}_chunk{chunk_number}.vcf"
    chunk_files.append(chunk_filename)

    # Open the input VCF using pysam
    vcf_in = pysam.VariantFile(input_vcf, 'r')
    # Create a new VCF file for the first chunk
    vcf_out = pysam.VariantFile(chunk_filename, 'w', header=vcf_in.header)

    for record in vcf_in:
        if variant_count >= 300000: #constants.MAX_SOMATIC_VARIANTS
            # Close the current chunk file and start a new one
            vcf_out.close()
            chunk_number += 1
            chunk_filename = output_dir / f"{base_filename}_chunk{chunk_number}.vcf"
            chunk_files.append(chunk_filename)
            vcf_out = pysam.VariantFile(chunk_filename, 'w', header=vcf_in.header)
            variant_count = 0
        # Write the record to the current chunk
        vcf_out.write(record)
        variant_count += 1

    # Close the last chunk file
    vcf_out.close()
    vcf_in.close()

    print(f"VCF file split into {len(chunk_files)} chunks.")

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
    print(f"Merged TSV written to: {merged_tsv_fp}")


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
    print("Running bcftools merge...")
    execute_command(command)
    print(f"Merged VCF written to: {merged_unsorted_vcf}")

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

    print("Sorting merged VCF file...")
    execute_command(sort_command)
    print(f"Sorted merged VCF written to: {merged_vcf}")

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

    print("Indexing sorted merged VCF file...")
    execute_command(index_command)
    print(f"Indexed merged VCF file: {merged_vcf}.tbi")

    # Optionally, remove the unsorted merged VCF file
    if merged_unsorted_vcf.exists():
        merged_unsorted_vcf.unlink()

    return merged_vcf