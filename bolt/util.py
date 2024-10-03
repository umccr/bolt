import pathlib
import subprocess
import sys
import os
import textwrap
import pysam


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

def split_vcf(input_vcf, output_dir, max_variants=100000000):
    """
    Splits a VCF file into multiple chunks, each containing up to max_variants variants.
    Each chunk includes the VCF header.
    """
    # Ensure output_dir exists
    # Count total number of variants using util.count_vcf_records
    total_variants = count_vcf_records(input_vcf)
    print(f"Total number of variants in the input VCF: {total_variants}")


    chunk_files = []
    chunk_number = 1
    variant_count = 0
    base_filename = os.path.splitext(os.path.basename(input_vcf))[0]
    chunk_filename = os.path.join(output_dir, f"{base_filename}_chunk{chunk_number}.vcf")
    chunk_files.append(chunk_filename)

    # Open the input VCF using pysam
    vcf_in = pysam.VariantFile(input_vcf, 'r')
    # Create a new VCF file for the first chunk
    vcf_out = pysam.VariantFile(chunk_filename, 'w', header=vcf_in.header)

    for record in vcf_in:
        if variant_count >= max_variants:
            # Close the current chunk file and start a new one
            vcf_out.close()
            chunk_number += 1
            chunk_filename = os.path.join(output_dir, f"{base_filename}_chunk{chunk_number}.vcf")
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
    merged_unsorted_vcf = merged_vcf_fp + '.unsorted.vcf.gz'
    merged_vcf = merged_vcf_fp + '.vcf.gz'
    cmd = ['bcftools', 'merge', '-m', 'all', '-Oz', '-o', merged_unsorted_vcf] + vcf_files
    # Run the bcftools merge command
    try:
        print("Running bcftools merge...")
        subprocess.run(cmd, check=True)
        print(f"Merged VCF written to: {merged_unsorted_vcf}")
    except subprocess.CalledProcessError as e:
        print(f"Error merging VCF files with bcftools:\n{e}")
        raise

    # Sort the merged VCF file
    cmd_sort = ['bcftools', 'sort', '-Oz', '-o', merged_vcf, merged_unsorted_vcf]
    try:
        print(f"Sorting merged VCF file...")
        subprocess.run(cmd_sort, check=True)
        print(f"Sorted merged VCF written to: {merged_vcf}")
    except subprocess.CalledProcessError as e:
        print(f"Error sorting merged VCF file:\n{e}")
        raise
