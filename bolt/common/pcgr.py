import csv
import pathlib
import re
import shutil
import tempfile
import time

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


def run_somatic(input_fp, pcgr_refdata_dir, pcgr_output_dir, chunk_nbr=None, threads=1, pcgr_conda=None, pcgrr_conda=None, purity=None, ploidy=None, sample_id=None):

    # NOTE(SW): Nextflow FusionFS v2.2.8 does not support PCGR output to S3; instead write to a
    # temporary directory outside of the FusionFS mounted directory then manually copy across

    temp_dir = tempfile.TemporaryDirectory()
    pcgr_output_dir = pcgr_output_dir / f"{chunk_nbr}" if chunk_nbr is not None else pcgr_output_dir    # Check if the output directory already exists
    if pcgr_output_dir.exists():
        print(f"Warning: Output directory '{pcgr_output_dir}' already exists.")
        print("Waiting for 3 seconds before overwriting...")
        #time.sleep(3)  # Wait for 3 seconds

        # If the user chooses to continue, delete the existing directory
        shutil.rmtree(pcgr_output_dir)


    if not sample_id:
        sample_id = 'nosampleset'

    command_args = [
        f'--sample_id {sample_id}',
        f'--input_vcf {input_fp}',
        f'--tumor_dp_tag TUMOR_DP',
        f'--tumor_af_tag TUMOR_AF',
        f'--control_dp_tag NORMAL_DP',
        f'--control_af_tag NORMAL_AF',
        f'--pcgr_dir {pcgr_refdata_dir}',
        f'--genome_assembly grch38',
        f'--assay WGS',
        f'--estimate_signatures',
        f'--estimate_msi_status',
        f'--estimate_tmb',
        f'--show_noncoding',
        f'--vcfanno_n_proc {threads}',
        f'--vep_pick_order biotype,rank,appris,tsl,ccds,canonical,length,mane',
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
    command_args.append(f'--output_dir {temp_dir.name}')

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

    util.execute_command(command)

    shutil.copytree(temp_dir.name, pcgr_output_dir)

    return pcgr_output_dir


def run_germline(input_fp, panel_fp, pcgr_refdata_dir, output_dir, threads=1, pcgr_conda=None, pcgrr_conda=None, sample_id=None):

    if not sample_id:
        sample_id = 'nosampleset'

    # NOTE(SW): Nextflow FusionFS v2.2.8 does not support PCGR output to S3; instead write to a
    # temporary directory outside of the FusionFS mounted directory then manually copy across

    temp_dir = tempfile.TemporaryDirectory()
    cpsr_output_dir = output_dir / 'cpsr/'

    command_args = [
        f'--sample_id {sample_id}',
        f'--input_vcf {input_fp}',
        f'--genome_assembly grch38',
        f'--custom_list {panel_fp}',
        # NOTE(SW): probably useful to add versioning information here; weigh against maintainence
        # burden
        f'--custom_list_name umccr_germline_panel',
        f'--pop_gnomad global',
        f'--classify_all',
        f'--pcgr_dir {pcgr_refdata_dir}',
        f'--vcfanno_n_proc {threads}',
        f'--vep_pick_order biotype,rank,appris,tsl,ccds,canonical,length,mane',
    ]

    if pcgrr_conda:
        command_args.append(f'--pcgrr_conda {pcgrr_conda}')

    # NOTE(SW): placed here to always have output directory last
    command_args.append(f'--output_dir {temp_dir.name}')

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

    shutil.copytree(temp_dir.name, cpsr_output_dir)

    return cpsr_output_dir


def transfer_annotations_somatic(input_fp, tumor_name, filter_name, pcgr_dir, output_dir):
    # Set destination INFO field names and source TSV fields
    info_field_map = {
        constants.VcfInfo.PCGR_MUTATION_HOTSPOT: 'MUTATION_HOTSPOT',
        constants.VcfInfo.PCGR_CLINVAR_CLNSIG: 'CLINVAR_CLNSIG',
        constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT: 'TCGA_PANCANCER_COUNT',
        constants.VcfInfo.PCGR_CSQ: 'CSQ',
    }

    pcgr_tsv_fp = pathlib.Path(pcgr_dir) / 'nosampleset.pcgr_acmg.grch38.snvs_indels.tiers.tsv'
    pcgr_vcf_fp = pathlib.Path(pcgr_dir) / 'nosampleset.pcgr_acmg.grch38.vcf.gz'

    # Enforce matching defined and source INFO annotations
    check_annotation_headers(info_field_map, pcgr_vcf_fp)

    # Gather PCGR annotation data for records
    pcgr_data = collect_pcgr_annotation_data(pcgr_tsv_fp, pcgr_vcf_fp, info_field_map)

    # Open filehandles, set required header entries
    input_fh = cyvcf2.VCF(input_fp)

    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_TIER)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_CSQ)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_MUTATION_HOTSPOT)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_CLINVAR_CLNSIG)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_COSMIC_COUNT)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.PCGR_ICGC_PCAWG_COUNT)

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
        record_ann = annotate_record(record, pcgr_data)
        output_fh.write_record(record_ann)


def transfer_annotations_germline(input_fp, normal_name, cpsr_dir, output_dir):
    # Set destination INFO field names and source TSV fields
    info_field_map = {
        constants.VcfInfo.CPSR_FINAL_CLASSIFICATION: 'FINAL_CLASSIFICATION',
        constants.VcfInfo.CPSR_PATHOGENICITY_SCORE: 'CPSR_PATHOGENICITY_SCORE',
        constants.VcfInfo.CPSR_CLINVAR_CLASSIFICATION: 'CLINVAR_CLASSIFICATION',
        constants.VcfInfo.CPSR_CSQ: 'CSQ',
    }

    cpsr_tsv_fp = pathlib.Path(cpsr_dir) / f'{normal_name}.cpsr.grch38.snvs_indels.tiers.tsv'
    cpsr_vcf_fp = pathlib.Path(cpsr_dir) / f'{normal_name}.cpsr.grch38.vcf.gz'

    # Enforce matching defined and source INFO annotations
    check_annotation_headers(info_field_map, cpsr_vcf_fp)

    # Gather CPSR annotation data for records
    cpsr_data = collect_cpsr_annotation_data(cpsr_tsv_fp, cpsr_vcf_fp, info_field_map)

    # Open filehandles, set required header entries
    input_fh = cyvcf2.VCF(input_fp)

    util.add_vcf_header_entry(input_fh, constants.VcfInfo.CPSR_FINAL_CLASSIFICATION)
    util.add_vcf_header_entry(input_fh, constants.VcfInfo.CPSR_PATHOGENICITY_SCORE)
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
            continue

        header_dst_entry = util.get_vcf_header_entry(header_dst)
        # Remove leading and trailing quotes from source
        header_src_description_unquoted = header_src_entry['Description'].strip('"')
        assert  header_src_description_unquoted == header_dst_entry['Description']


def collect_pcgr_annotation_data(tsv_fp, vcf_fp, info_field_map):
    # Gather all annotations from TSV
    data_tsv = dict()
    with open(tsv_fp, 'r') as tsv_fh:
        for record in csv.DictReader(tsv_fh, delimiter='\t'):
            key, record_ann = get_annotation_entry_tsv(record, info_field_map)
            assert key not in data_tsv

            # Process PCGR_TIER
            # TIER_1, TIER_2, TIER_3, TIER_4, NONCODING
            record_ann[constants.VcfInfo.PCGR_TIER] = record['TIER'].replace(' ', '_')

            # Count COSMIC hits
            if record['COSMIC_MUTATION_ID'] == 'NA':
                cosmic_count = 0
            else:
                cosmic_count = len(record['COSMIC_MUTATION_ID'].split('&'))
            record_ann[constants.VcfInfo.PCGR_COSMIC_COUNT] = cosmic_count

            # Count ICGC-PCAWG hits by taking sum of affected donors where the annotation value has
            # the following format: project_code|tumor_type|affected_donors|tested_donors|frequency
            icgc_pcawg_count = 0
            if record['ICGC_PCAWG_OCCURRENCE'] != 'NA':
                for pcawg_hit_data in record['ICGC_PCAWG_OCCURRENCE'].split(','):
                    pcawrg_hit_data_fields = pcawg_hit_data.split('|')
                    affected_donors = int(pcawrg_hit_data_fields[2])
                    icgc_pcawg_count += affected_donors
                assert icgc_pcawg_count > 0
            record_ann[constants.VcfInfo.PCGR_ICGC_PCAWG_COUNT] = icgc_pcawg_count

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
    with open(tsv_fp, 'r') as tsv_fh:
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
    # Set lookup key; PCGR/CPSR strips leading 'chr' from contig names
    chrom = f'chr{record["CHROM"]}'
    pos = int(record['POS'])
    key = (chrom, pos, record['REF'], record['ALT'])

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
