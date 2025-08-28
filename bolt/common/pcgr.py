import collections
import csv
import functools
import gzip
import itertools
import pathlib
import re
import shutil
import tempfile

import cyvcf2


from .. import util
from ..common import constants


def prepare_vcf_somatic(input_fp, tumor_name, normal_name, output_dir):


    # TODO(SW): complete comments


    # PCGR requires tumor and normal FORMAT/DP and FORMAT/AF to be moved to the INFO field
    #   * <URL>

    # To generate the PCGR input, I exclude all unnecessary data from the VCF and set tumor
    # and normal FORMAT/AF and FORMAT/DP annotations as INFO annotations as required by PCGR.

    input_fh = cyvcf2.VCF(input_fp)

    tumor_index = input_fh.samples.index(tumor_name)
    normal_index = input_fh.samples.index(normal_name)
    assert tumor_name != normal_name
    assert tumor_index != normal_index

    output_fp = output_dir / f'{tumor_name}.pcgr_prep.vcf.gz'
    output_fh = cyvcf2.Writer.from_string(output_fp, get_minimal_header(input_fh), 'wz')

    for record in input_fh:
        # Collect tumor and normal FORMAT/AF and FORMAT/DP
        [tumor_dp] = record.format('DP')[tumor_index]
        [tumor_af] = record.format('AF')[tumor_index]
        [normal_dp] = record.format('DP')[normal_index]
        [normal_af] = record.format('AF')[normal_index]

        # Manually create INFO entries
        info_entries = (
            f'{constants.VcfInfo.TUMOR_AF.value}={tumor_af:.3}',
            f'{constants.VcfInfo.TUMOR_DP.value}={tumor_dp}',
            f'{constants.VcfInfo.NORMAL_AF.value}={normal_af:.3}',
            f'{constants.VcfInfo.NORMAL_DP.value}={normal_dp}',
        )
        info = ';'.join(info_entries)

        # Construct clean record containing new INFO data, also set FILTER=PASS and remove all
        # FORMAT and sample columns
        record_variant_comps = str(record).split('\t')[:6]  # CHROM, POS, ID, REF, ALT, QUAL
        record_str_new = '\t'.join([*record_variant_comps, 'PASS', info])
        record_new = output_fh.variant_from_string(record_str_new)
        output_fh.write_record(record_new)

    output_fh.close()

    # Index output
    command = fr'''bcftools index -t {output_fp}'''
    util.execute_command(command)

    return output_fp


def prepare_vcf_germline(input_fp, normal_name, output_dir):

    output_fp = output_dir / f'{normal_name}.cpsr.prep.vcf.gz'

    command = fr'''
        bcftools view -s {normal_name} {input_fp} | \
            bcftools annotate -x INFO,FILTER,FORMAT,^GT -o {output_fp};
            bcftools index -t {output_fp};
        '''

    util.execute_command(command)
    return output_fp


def get_minimal_header(input_fh):
    # Get a minimal VCF header for the PCGR input VCF
    # Filetype line
    filetype_line = '##fileformat=VCFv4.2'

    # Chromosome lines
    # NOTE(SW): the purpose of using an existing header is to obtain compatibile contig size for
    # each chromosome
    chrom_lines = list()
    chrom_prefix = [f'##contig=<ID={e},' for e in constants.CONTIGS_MAIN]
    header_lines = input_fh.raw_header.rstrip().split('\n')
    for line in header_lines:
        if any(line.startswith(c) for c in chrom_prefix):
            chrom_lines.append(line)

    # FORMAT lines
    format_lines = (
        util.get_vcf_header_line(constants.VcfInfo.TUMOR_AF),
        util.get_vcf_header_line(constants.VcfInfo.TUMOR_DP),
        util.get_vcf_header_line(constants.VcfInfo.NORMAL_AF),
        util.get_vcf_header_line(constants.VcfInfo.NORMAL_DP),
    )

    # Column header line
    column_line = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'

    # Construct and return
    return '\n'.join([filetype_line, *chrom_lines, *format_lines, column_line])


def run_somatic(input_fp, pcgr_refdata_dir, vep_dir, output_dir, chunk_nbr=None, threads=1, pcgr_conda=None, pcgrr_conda=None, purity=None, ploidy=None, sample_id=None):


    output_dir = output_dir / f"pcgr_{chunk_nbr}" if chunk_nbr is not None else output_dir

    if output_dir.exists():
        logger.warning(f"Output directory '{output_dir}' already exists and will be overwritten")
        shutil.rmtree(output_dir)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    if not sample_id:
        sample_id = 'nosampleset'

    command_args = [
        f'--sample_id {sample_id}',
        f'--input_vcf {input_fp}',
        f'--vep_dir {vep_dir}',
        f'--refdata_dir {pcgr_refdata_dir}',
        f'--tumor_dp_tag TUMOR_DP',
        f'--tumor_af_tag TUMOR_AF',
        f'--control_dp_tag NORMAL_DP',
        f'--control_af_tag NORMAL_AF',
        f'--genome_assembly grch38',
        f'--assay WGS',
        f'--estimate_signatures',
        f'--estimate_msi',
        f'--estimate_tmb',
        f'--vcfanno_n_proc {threads}',
        f'--vep_n_forks 4',
        f'--vep_pick_order biotype,rank,appris,tsl,ccds,canonical,length,mane_plus_clinical,mane_select',
    ]

    # NOTE(SW): VEP pick order is applied as a successive filter:
    #   * sort variants within the next category
    #   * select variants with lowest order value (i.e. best within category)
    #   * if only one selected variants return it, else proceed to next category
    #
    # Where pick categories have been exhausted and multiple variants remain tied, the first
    # variant in the sort array is returned.
    #
    # The 'biotype' category only sorts on basis of a variant impact being coding or non-coding.
    #
    # Cat sort: https://github.com/Ensembl/ensembl-vep/blob/105.0/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L718
    # Pick order: https://github.com/Ensembl/ensembl-vep/blob/105.0/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L760

    if pcgrr_conda:
        command_args.append(f'--pcgrr_conda {pcgrr_conda}')

    if purity:
        command_args.append(f'--tumor_purity {purity}')

    if ploidy:
        command_args.append(f'--tumor_ploidy {ploidy}')

    # NOTE(SW): placed here to always have output directory last
    command_args.append(f'--output_dir {output_dir}')

    delimiter_padding = ' ' * 10
    delimiter = f' \\\n{delimiter_padding}'

    command_args_str = delimiter.join(command_args)

    command = fr'''
        pcgr \
          {command_args_str}
    '''

    if pcgr_conda:
        command_conda = f'conda run -n {pcgr_conda} \\'
        command_formatting = '\n' + ' ' * 4
        command = command_formatting + command_conda + command

    # Log file path
    log_file_path = output_dir / "run_somatic.log"

    # Run the command and redirect output to the log file
    util.execute_command(command, log_file_path=log_file_path)

    pcgr_tsv_fp = pathlib.Path(output_dir) / f'{sample_id}.pcgr.grch38.snv_indel_ann.tsv.gz'
    pcgr_vcf_fp = pathlib.Path(output_dir) / f'{sample_id}.pcgr.grch38.pass.vcf.gz'

    # Check if both files exist
    if not pcgr_tsv_fp.exists():
        raise FileNotFoundError(f"Expected file {pcgr_tsv_fp} not found.")
    if not pcgr_vcf_fp.exists():
        raise FileNotFoundError(f"Expected file {pcgr_vcf_fp} not found.")

    return pcgr_tsv_fp, pcgr_vcf_fp


def run_germline(input_fp, panel_fp, pcgr_refdata_dir, vep_dir, output_dir, threads=1, pcgr_conda=None, pcgrr_conda=None, sample_id=None):

    if not sample_id:
        sample_id = 'nosampleset'


    cpsr_output_dir = output_dir / 'cpsr/'

    if cpsr_output_dir.exists():
        logger.warning(f"Output directory '{cpsr_output_dir}' already exists and will be overwritten")
        shutil.rmtree(cpsr_output_dir)

    # Create output directory
    cpsr_output_dir.mkdir(parents=True, exist_ok=True)

    command_args = [
        f'--sample_id {sample_id}',
        f'--input_vcf {input_fp}',
        f'--genome_assembly grch38',
        f'--custom_list {panel_fp}',
        f'--vep_dir {pcgr_refdata_dir}',
        f'--refdata_dir {pcgr_refdata_dir}',
        # NOTE(SW): probably useful to add versioning information here; weigh against maintainence
        # burden
        f'--custom_list_name umccr_germline_panel',
        f'--pop_gnomad global',
        f'--classify_all',
        f'--vcfanno_n_proc {threads}',
        f'--vep_pick_order biotype,rank,appris,tsl,ccds,canonical,length,mane_plus_clinical,mane_select',
    ]

    if pcgrr_conda:
        command_args.append(f'--pcgrr_conda {pcgrr_conda}')

    # NOTE(SW): placed here to always have output directory last
    command_args.append(f'--output_dir {cpsr_output_dir}')

    delimiter_padding = ' ' * 10
    delimiter = f' \\\n{delimiter_padding}'

    command_args_str = delimiter.join(command_args)


    command = fr'''
        cpsr \
          {command_args_str}
    '''
    if pcgr_conda:
        command_conda = f'conda run -n {pcgr_conda} \\'
        command_formatting = '\n' + ' ' * 4
        command = command_formatting + command_conda + command

    util.execute_command(command)

    return cpsr_output_dir


def transfer_annotations_somatic(input_fp, tumor_name, filter_name, pcgr_dir, output_dir):
    # Set destination INFO field names and source TSV fields
    info_field_map = {
        constants.VcfInfo.PCGR_MUTATION_HOTSPOT: 'MUTATION_HOTSPOT',
        constants.VcfInfo.PCGR_CLINVAR_CLASSIFICATION: 'CLINVAR_CLASSIFICATION',
        constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT: 'TCGA_PANCANCER_COUNT',
        constants.VcfInfo.PCGR_CSQ: 'CSQ',
    }

    pcgr_tsv_fp = pathlib.Path(output_dir) / 'nosampleset.pcgr_acmg.grch38.snvs_indels.tiers.tsv'
    pcgr_vcf_fp = pathlib.Path(output_dir) / 'nosampleset.pcgr_acmg.grch38.vcf.gz'

    # Enforce matching defined and source INFO annotations
    check_annotation_headers(info_field_map, pcgr_vcf_fp)

    # Gather PCGR annotation data for records
    pcgr_data = collect_pcgr_annotation_data(pcgr_tsv_fp, pcgr_vcf_fp, info_field_map)

    # Open filehandles, set required header entries
    input_fh = cyvcf2.VCF(input_fp)

    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_ACTIONABILITY_TIER)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_CSQ)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_MUTATION_HOTSPOT)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_CLINVAR_CLASSIFICATION)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT)

    output_fp = output_dir / f'{tumor_name}.annotations.vcf.gz'
    output_fh = cyvcf2.Writer(output_fp, input_fh, 'wz')

    # Transfer annotations and write to output
    for record in input_fh:
        # Do not process chrM since *snvs_indels.tiers.tsv does not include these annotations
        if record.CHROM == 'chrM':
            continue
        # Immediately print out variants that were not annotated
        if filter_name in record.FILTERS:
            output_fh.write_record(record)
            continue
        # Annotate and write
        record_ann = annotate_record(record, pcgr_data, allow_missing=True)
        output_fh.write_record(record_ann)


def transfer_annotations_germline(input_fp, normal_name, cpsr_dir, output_dir):
    # Set destination INFO field names and source TSV fields
    # Note: Only include fields that exist in CPSR v2.2.1 output
    info_field_map = {
        constants.VcfInfo.CPSR_CLINVAR_CLASSIFICATION: 'CLINVAR_CLASSIFICATION',
        constants.VcfInfo.CPSR_CSQ: 'CSQ',
    }

    cpsr_tsv_fp = pathlib.Path(cpsr_dir) / f'{normal_name}.cpsr.grch38.classification.tsv.gz'
    cpsr_vcf_fp = pathlib.Path(cpsr_dir) / f'{normal_name}.cpsr.grch38.vcf.gz'

    # Enforce matching defined and source INFO annotations
    check_annotation_headers(info_field_map, cpsr_vcf_fp)

    # Gather CPSR annotation data for records
    cpsr_data = collect_cpsr_annotation_data(cpsr_tsv_fp, cpsr_vcf_fp, info_field_map)

    # Open filehandles, set required header entries
    input_fh = cyvcf2.VCF(input_fp)

    util.add_vcf_header_entry(input_fh, constants.VcfInfo.CPSR_CLINVAR_CLASSIFICATION)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.CPSR_CSQ)

    output_fp = output_dir / f'{normal_name}.annotations.vcf.gz'
    output_fh = cyvcf2.Writer(output_fp, input_fh, 'wz')

    # Transfer annotations and write to output
    for record in input_fh:
        # Do not process chrM since *snvs_indels.tiers.tsv does not include these annotations
        if record.CHROM == 'chrM':
            continue
        # Annotate and write
        # NOTE(SW): allow missing CPSR annotations for input variants, CPSR seems to drop some
        record_ann = annotate_record(record, cpsr_data, allow_missing=True)
        output_fh.write_record(record_ann)


def check_annotation_headers(info_field_map, vcf_fp):
    # Ensure header descriptions from source INFO annotations match those defined here for the
    # output file; force manual inspection where they do not match
    vcf_fh = cyvcf2.VCF(vcf_fp)
    for header_dst, header_src in info_field_map.items():

        # Skip header lines that do not have an equivalent entry in the PCGR/CPSR VCF
        try:
            header_src_entry = vcf_fh.get_header_type(header_src)
        except KeyError:
            print(f"Missing in VCF: {header_src}")
            continue

        header_dst_entry = util.get_vcf_header_entry(header_dst)
        # Remove leading and trailing quotes from source
        header_src_description_unquoted = header_src_entry['Description'].strip('"')
        if header_src_description_unquoted != header_dst_entry['Description']:
            raise AssertionError(f"Mismatch for {header_src}:\nVCF: {header_src_description_unquoted}\nExpected: {header_dst_entry['Description']}")

def collect_pcgr_annotation_data(tsv_fp, vcf_fp, info_field_map):
    # Gather all annotations from TSV
    data_tsv = dict()
    with open(tsv_fp, 'r') as tsv_fh:
        for record in csv.DictReader(tsv_fh, delimiter='\t'):
            key, record_ann = get_annotation_entry_tsv(record, info_field_map)
            assert key not in data_tsv

            # Process PCGR_ACTIONABILITY_TIER
            # TIER 1, TIER 2, TIER 3, TIER 4, NONCODING
            record_ann[constants.VcfInfo.PCGR_ACTIONABILITY_TIER] = record['ACTIONABILITY_TIER'].replace(' ', '_')

            # Process PCGR_ACTIONABILITY_TIER
            # TIER 1, TIER 2, TIER 3, TIER 4, NONCODING
            record_ann[constants.VcfInfo.PCGR_ACTIONABILITY_TIER] = record['ACTIONABILITY_TIER']

            # Store annotation data
            data_tsv[key] = record_ann

    # Gather data from VCF
    data_vcf = get_annotations_vcf(vcf_fp, info_field_map)

    # Compile annotations, prefering TSV source
    return compile_annotation_data(data_tsv, data_vcf)


def collect_cpsr_annotation_data(tsv_fp, vcf_fp, info_field_map):
    # Gather annotations from TSV
    data_tsv = dict()
    gdot_re = re.compile('^(?P<chrom>[\dXYM]+):g\.(?P<pos>\d+)(?P<ref>[A-Z]+)>(?P<alt>[A-Z]+)$')
    with gzip.open(tsv_fp, 'rt') as tsv_fh:
        for record in csv.DictReader(tsv_fh, delimiter='\t'):
            # Decompose CPSR 'GENOMIC_CHANGE' field into CHROM, POS, REF, and ALT
            re_result = gdot_re.match(record['GENOMIC_CHANGE'])
            if not re_result:
                print(record['GENOMIC_CHANGE'])
                assert re_result
            record['CHROM'] = re_result.group('chrom')
            record['POS'] = re_result.group('pos')
            record['REF'] = re_result.group('ref')
            record['ALT'] = re_result.group('alt')

            key, record_ann = get_annotation_entry_tsv(record, info_field_map)
            assert key not in data_tsv
            data_tsv[key] = record_ann

    # Gather annotations from VCF
    data_vcf = get_annotations_vcf(vcf_fp, info_field_map)

    # Compile annotations, prefering TSV source
    return compile_annotation_data(data_tsv, data_vcf)

def parse_genomic_change(genomic_change):
    """
    Parse a genomic change string, e.g., "3:g.41224645T>C"
    Returns a tuple: (chrom, pos, ref, alt)
    """
    # Regular expression for the format "chrom:g.posRef>Alt"
    pattern = r'^(?P<chrom>\w+):g\.(?P<pos>\d+)(?P<ref>\w+)>(?P<alt>\w+)$'
    match = re.match(pattern, genomic_change)
    if not match:
        raise ValueError(f"Format not recognized: {genomic_change}")
    
    # Get values and format as needed
    chrom = f"chr{match.group('chrom')}"
    pos = int(match.group('pos'))
    ref = match.group('ref')
    alt = match.group('alt')
    return chrom, pos, ref, alt


def get_annotations_vcf(vcf_fp, info_field_map):
    data_vcf = dict()
    for record in cyvcf2.VCF(vcf_fp):
        # Set lookup key; PCGR strips leading 'chr' from contig names
        assert len(record.ALT) == 1
        [alt] = record.ALT
        key = (f'chr{record.CHROM}', record.POS, record.REF, alt)
        assert key not in data_vcf

        data_vcf[key] = dict()
        for info_dst, info_src in info_field_map.items():
            if (info_val := record.INFO.get(info_src)):
                data_vcf[key][info_dst] = info_val

    return data_vcf


def get_annotation_entry_tsv(record, info_field_map):
    # If GENOMIC_CHANGE is present, parse it for coordinates; otherwise, use separate fields.
    if "GENOMIC_CHANGE" in record and record["GENOMIC_CHANGE"]:
        chrom, pos, ref, alt = parse_genomic_change(record["GENOMIC_CHANGE"])

    if not chrom.startswith('chr'):
        chrom = f'chr{chrom}'

    key = (chrom, pos, ref, alt)

    record_ann = dict()
    for info_dst, info_src in info_field_map.items():
        if not (info_val := record.get(info_src)):
            continue
        elif info_val == 'NA':
            continue

        record_ann[info_dst] = info_val

    return key, record_ann


def compile_annotation_data(data_tsv, data_vcf):
    # Compile annotations, prefering TSV as source
    annotations = data_tsv
    for key, data_vcf_record in data_vcf.items():

        if key not in annotations:
            annotations[key] = dict()

        for info_name, info_val in data_vcf_record.items():

            if info_name in annotations[key]:
                continue

            annotations[key][info_name] = info_val
    return annotations


def annotate_record(record, annotations, *, allow_missing=False):
    # Get lookup key
    assert len(record.ALT) == 1
    [alt] = record.ALT
    key = (record.CHROM, record.POS, record.REF, alt)

    # Handle missing entries
    if key not in annotations:
        if allow_missing:
            return record
        else:
            assert key not in annotations

    # Transfer annotations
    for info_enum, v in annotations[key].items():
        record.INFO[info_enum.value] = v

    return record

def split_vcf(input_vcf, output_dir):
    """
    Splits a VCF file into multiple chunks, each containing up to max_variants variants.
    Each chunk includes the VCF header.
    Ensures no overlapping positions between chunks.
    """
    output_dir = pathlib.Path(output_dir / "vcf_chunks")
    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_files = []
    chunk_number = 1
    variant_count = 0
    base_filename = pathlib.Path(input_vcf).stem
    base_filename = input_vcf.stem
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
        if variant_count >= constants.MAX_SOMATIC_VARIANTS and (last_position is None or current_position != last_position):
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

def run_somatic_chunck(vcf_chunks, pcgr_data_dir, vep_dir, output_dir, pcgr_output_dir, max_threads, pcgr_conda, pcgrr_conda):
    pcgr_tsv_files = []
    pcgr_vcf_files = []
    num_chunks = len(vcf_chunks)
    # Ensure we don't use more workers than available threads, and each worker has at least 2 threads
    max_workers = min(num_chunks, max_threads // 2)
    threads_quot, threads_rem = divmod(max_threads, num_chunks)
    threads_per_chunk = max(2, threads_quot)
    # Limit the number of workers to the smaller of num_chunks or max_threads // 2
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for chunk_number, vcf_file in enumerate(vcf_chunks, start=1):
            # Assign extra thread to the first 'threads_rem' chunks
            additional_thread = 1 if chunk_number <= threads_rem else 0
            total_threads = threads_per_chunk + additional_thread
            futures[executor.submit(run_somatic, vcf_file, pcgr_data_dir, vep_dir, pcgr_output_dir, chunk_number, total_threads, pcgr_conda, pcgrr_conda)] = chunk_number
        for future in concurrent.futures.as_completed(futures):
            try:
                pcgr_tsv_fp, pcgr_vcf_fp = future.result()
                if pcgr_tsv_fp:
                    pcgr_tsv_files.append(pcgr_tsv_fp)
                if pcgr_vcf_fp:
                    pcgr_vcf_files.append(pcgr_vcf_fp)
            except Exception as e:
                print(f"Exception occurred: {e}")
    merged_vcf_fp, merged_tsv_fp = merging_pcgr_files(output_dir, pcgr_vcf_files, pcgr_tsv_files)
    return merged_tsv_fp, merged_vcf_fp

def merging_pcgr_files(output_dir, pcgr_vcf_files, pcgr_tsv_fp):
    pcgr_dir = pathlib.Path(output_dir) / 'pcgr'
    pcgr_dir.mkdir(exist_ok=True)

    # Merge all TSV files into a single file in the pcgr directory
    merged_tsv_fp = pcgr_dir / "nosampleset.pcgr_acmg.grch38.snvs_indels.tiers.tsv"
    util.merge_tsv_files(pcgr_tsv_fp, merged_tsv_fp)

    # Step 5: Merge all VCF files into a single file in the pcgr directory
    merged_vcf_path = pcgr_dir / "nosampleset.pcgr.grch38.pass.vcf.gz"
    merged_vcf = util.merge_vcf_files(pcgr_vcf_files, merged_vcf_path)

    return merged_vcf, merged_tsv_fp


def get_variant_filter_data(variant):
    attribute_names = (
        'tier',
        'difficult',
        'giab_conf',
        'intergenic',
        'intronic',
        'downstream',
        'upstream',
        'impacts_other',
    )

    data = {e: None for e in attribute_names}


    data['tier'] = variant.INFO['PCGR_TIER']


    info_keys = [k for k, v in variant.INFO]

    data['difficult'] = any(e.startswith('DIFFICULT') for e in info_keys)
    data['giab_conf'] = 'GIAB_CONF' in info_keys


    # NOTE(SW): GIAB_CONF always overrides DIFFICULT tags
    if data['giab_conf'] and data['difficult']:
        data['difficult']= False


    for impact in get_impacts(variant.INFO['PCGR_CSQ']):
        if impact == 'intergenic_variant':
            data['intergenic'] = True
        elif impact == 'intron_variant':
            data['intronic'] = True
        elif impact == 'downstream_gene_variant':
            data['downstream'] = True
        elif impact == 'upstream_gene_variant':
            data['upstream'] = True
        elif impact:
            data['impacts_other'] = True
        else:
            assert False

    return data


def get_impacts(csq_str_full):
    impacts = set()
    for csq_str in csq_str_full.split(','):
        csq_tokens = csq_str.split('|')
        impact_str = csq_tokens[1]
        impacts.update(impact_str.split('&'))
    return impacts


def determine_filter(data):

    for impact, region in get_ordering(tiers=False):

        # NOTE(SW): this is less efficient than nested loops since the outer block is reevaluated
        # within what would be the inner loop each cycle; taking this route for cleaner code

        impacts_higher = get_impacts_higher(impact)
        impact_filter = bool(data[impact]) and not any(bool(data[e]) for e in impacts_higher)

        region_filter = False
        if region == 'none':
            region_filter = not (data['difficult'] or data['giab_conf'])
        else:
            region_filter = data[region]

        if impact_filter and region_filter:
            return (impact, region)

    return False


def get_variant_repr(variant):
    return (variant.CHROM, variant.POS, variant.REF, tuple(variant.ALT))


@functools.cache
def get_ordering(tiers=True, impacts=True, regions=True):
    categories = [
        constants.PCGR_TIERS_FILTERING if tiers else None,
        constants.VEP_IMPACTS_FILTER if impacts else None,
        constants.GENOMIC_REGIONS_FILTERING if regions else None,
    ]

    # NOTE(SW): I'm not aware of any noncoding impacts for TIER_[1-4] other than TERT but keeping
    # in to be overly cautious
    ordering_iter = itertools.product(*(c for c in categories if c))

    return tuple(ordering_iter)


@functools.cache
def get_impacts_higher(impact):
    impact_index = constants.VEP_IMPACTS_FILTER.index(impact)
    if impact_index + 1 < len(constants.VEP_IMPACTS_FILTER):
        impacts_higher = constants.VEP_IMPACTS_FILTER[impact_index+1:len(constants.VEP_IMPACTS_FILTER)]
    else:
        impacts_higher = list()
    return impacts_higher
