import csv
import pathlib
import shutil
import tempfile


import cyvcf2


from .. import util
from ..common import constants


def prepare_vcf_somatic(in_fp, fn_prefix, tumor_name, normal_name):


    # TODO(SW): complete comments


    # PCGR requires tumor and normal FORMAT/DP and FORMAT/AF to be moved to the INFO field
    #   * <URL>

    # To generate the PCGR input, I exclude all unnecessary data from the VCF and set tumor
    # and normal FORMAT/AF and FORMAT/DP annotations as INFO annotations as required by PCGR.

    in_fh = cyvcf2.VCF(in_fp)

    tumor_index = in_fh.samples.index(tumor_name)
    normal_index = in_fh.samples.index(normal_name)
    assert tumor_name != normal_name
    assert tumor_index != normal_index

    out_fp = f'{fn_prefix}.pcgr_prep.vcf.gz'
    out_fh = cyvcf2.Writer.from_string(out_fp, get_minimal_header(in_fh), 'wz')

    for record in in_fh:
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
        record_new = out_fh.variant_from_string(record_str_new)
        out_fh.write_record(record_new)
    out_fh.close()

    # Index output
    command = fr'''bcftools index -t {out_fp}'''
    util.execute_command(command)

    return out_fp


def prepare_vcf_germline(in_fp, fn_prefix, normal_name):

    out_fp = f'{fn_prefix}.cpsr.prep.vcf.gz'
    command = fr'''
        bcftools view -s {normal_name} {in_fp} | \
            bcftools annotate -x INFO,FILTER,FORMAT,^GT -o {out_fp};
            bcftools index -t {out_fp};
        '''
    util.execute_command(command)

    return out_fp


def get_minimal_header(in_fh):
    # Get a minimal VCF header for the PCGR input VCF
    # Filetype line
    filetype_line = '##fileformat=VCFv4.2'

    # Chromosome lines
    # NOTE(SW): the purpose of using an existing header is to obtain compatibile contig size for
    # each chromosome
    chrom_lines = list()
    chrom_prefix = [f'##contig=<ID={e},' for e in constants.CONTIGS_MAIN]
    header_lines = in_fh.raw_header.rstrip().split('\n')
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


def run_somatic(in_fp, pcgr_dir, threads=1, pcgr_conda=None, pcgrr_conda=None, purity=None, ploidy=None, sample_id=None):

    # NOTE(SW): Nextflow FusionFS v0.6.5 does not support PCGR output to S3; instead write to a
    # temporary directory outside of the FusionFS mounted directory then manually copy across

    out_dir = 'pcgr_output/'
    tmp_dir = tempfile.TemporaryDirectory()

    if not sample_id:
        sample_id = 'nosampleset'

    command_args = [
        f'--sample_id {sample_id}',
        f'--input_vcf {in_fp}',
        f'--tumor_dp_tag TUMOR_DP',
        f'--tumor_af_tag TUMOR_AF',
        f'--control_dp_tag NORMAL_DP',
        f'--control_af_tag NORMAL_AF',
        f'--genome_assembly grch38',
        f'--assay WGS',
        f'--estimate_signatures',
        f'--estimate_msi_status',
        f'--estimate_tmb',
        f'--show_noncoding',
        f'--vcfanno_n_proc {threads}',
        f'--vep_pick_order biotype,rank,appris,tsl,ccds,canonical,length,mane',
    ]

    if pcgrr_conda:
        command_args.append(f'--pcgr_dir {pcgr_dir}')

    if purity:
        command_args.append(f'--tumor_purity {purity}')

    if ploidy:
        command_args.append(f'--tumor_ploidy {ploidy}')

    # NOTE(SW): placed here to always have output directory last
    command_args.append(f'--output_dir {tmp_dir.name}')

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

    shutil.copytree(tmp_dir.name, out_dir)

    return out_dir


def run_germline(in_fp, panel_fp, pcgr_dir, threads=1, pcgr_conda=None, pcgrr_conda=None, sample_id=None):

    if not sample_id:
        sample_id = 'nosampleset'

    # NOTE(SW): Nextflow FusionFS v0.6.5 does not support PCGR output to S3; instead write to a
    # temporary directory outside of the FusionFS mounted directory then manually copy across

    out_dir = 'cpsr_output/'
    tmp_dir = tempfile.TemporaryDirectory()

    command_args = [
        f'--sample_id {sample_id}',
        f'--input_vcf {in_fp}',
        f'--genome_assembly grch38',
        f'--custom_list {panel_fp}',
        # NOTE(SW): probably useful to add versioning information here; weigh against maintainence
        # burden
        f'--custom_list_name umccr_germline_panel',
        f'--pop_gnomad global',
        f'--classify_all',
        f'--pcgr_dir {pcgr_dir}',
        f'--vcfanno_n_proc {threads}',
        f'--vep_pick_order biotype,rank,appris,tsl,ccds,canonical,length,mane',
    ]

    if pcgrr_conda:
        command_args.append(f'--pcgr_dir {pcgr_dir}')

    # NOTE(SW): placed here to always have output directory last
    command_args.append(f'--output_dir {tmp_dir.name}')

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

    shutil.copytree(tmp_dir.name, out_dir)

    return out_dir


def transfer_annotations(src_fp, fn_prefix, pcgr_dir, filter_name):
    # Set destination INFO field names and source TSV fields
    info_field_map = {
        constants.VcfInfo.PCGR_MUTATION_HOTSPOT: 'MUTATION_HOTSPOT',
        constants.VcfInfo.PCGR_CLINVAR_CLNSIG: 'CLINVAR_CLNSIG',
        constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT: 'TCGA_PANCANCER_COUNT',
        constants.VcfInfo.PCGR_CSQ: 'CSQ',
    }

    pcgr_tsv_fp = pathlib.Path(pcgr_dir) / 'nosampleset.pcgr_acmg.grch38.snvs_indels.tiers.tsv'
    pcgr_vcf_fp = pathlib.Path(pcgr_dir) / 'nosampleset.pcgr_acmg.grch38.vcf.gz'

    # Ensure header descriptions from source INFO annotations match those defined here for the
    # output file; force manual inspection where they do not match
    pcgr_vcf_fh = cyvcf2.VCF(pcgr_vcf_fp)
    for header_dst, header_src in info_field_map.items():
        header_src_entry = pcgr_vcf_fh.get_header_type(header_src)
        header_dst_entry = util.get_vcf_header_entry(header_dst)
        # Remove leading and trailing quotes from source
        header_src_description_unquoted = header_src_entry['Description'].strip('"')
        assert  header_src_description_unquoted == header_dst_entry['Description']

    # Gather PCGR annotation data for records
    pcgr_data = collect_pcgr_annotation_data(pcgr_tsv_fp, pcgr_vcf_fp, info_field_map)

    # Open filehandles, set required header entries
    src_fh = cyvcf2.VCF(src_fp)

    util.add_vcf_header_entry(src_fh, constants.VcfInfo.PCGR_TIER)
    util.add_vcf_header_entry(src_fh, constants.VcfInfo.PCGR_CSQ)
    util.add_vcf_header_entry(src_fh, constants.VcfInfo.PCGR_MUTATION_HOTSPOT)
    util.add_vcf_header_entry(src_fh, constants.VcfInfo.PCGR_CLINVAR_CLNSIG)
    util.add_vcf_header_entry(src_fh, constants.VcfInfo.PCGR_COSMIC_COUNT)
    util.add_vcf_header_entry(src_fh, constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT)
    util.add_vcf_header_entry(src_fh, constants.VcfInfo.PCGR_ICGC_PCAWG_COUNT)

    out_fp = f'{fn_prefix}.annotations.vcf.gz'
    out_fh = cyvcf2.Writer(out_fp, src_fh, 'wz')

    # Transfer annotations and write to output
    for record in src_fh:


        # Do not process chrM since *snvs_indels.tiers.tsv does not include these annotations
        if record.CHROM == 'chrM':
            continue

        # Immediately print out variants that were not annotated
        if filter_name in record.FILTERS:
            out_fh.write_record(record)
            continue

        # Get lookup key
        assert len(record.ALT) == 1
        [alt] = record.ALT
        key = (record.CHROM, record.POS, record.REF, alt)
        # Annotate then write
        assert key in pcgr_data
        for info_enum, v in pcgr_data[key].items():
            record.INFO[info_enum.value] = v
        out_fh.write_record(record)


def collect_pcgr_annotation_data(tsv_fp, vcf_fp, info_field_map):
    # First get all CSQ from VCF (not available in the TSV)
    data_csq = dict()
    for record in cyvcf2.VCF(vcf_fp):
        # Set lookup key; PCGR strips leading 'chr' from contig names
        assert len(record.ALT) == 1
        [alt] = record.ALT
        key = (f'chr{record.CHROM}', record.POS, record.REF, alt)
        assert key not in data_csq
        data_csq[key] = record.INFO.get('CSQ')

    # Gather all annotations
    data = dict()
    with open(tsv_fp, 'r') as tsv_fh:
        for record in csv.DictReader(tsv_fh, delimiter='\t'):
            # Set lookup key; PCGR strips leading 'chr' from contig names
            chrom = f'chr{record["CHROM"]}'
            pos = int(record['POS'])
            key = (chrom, pos, record['REF'], record['ALT'])
            assert key not in data

            # Process INFO annotations that can be directly transferred for the current record
            data_rec = dict()
            for info_dst, info_src in info_field_map.items():
                if info_src == 'CSQ':
                    info_val = data_csq[key]
                else:
                    info_val = record[info_src]

                if info_val == 'NA':
                    continue
                data_rec[info_dst] = info_val

            # Process PCGR_TIER
            # TIER_1, TIER_2, TIER_3, TIER_4, NONCODING
            data_rec[constants.VcfInfo.PCGR_TIER] = record['TIER'].replace(' ', '_')

            # Count COSMIC hits
            if record['COSMIC_MUTATION_ID'] == 'NA':
                cosmic_count = 0
            else:
                cosmic_count = len(record['COSMIC_MUTATION_ID'].split('&'))
            data_rec[constants.VcfInfo.PCGR_COSMIC_COUNT] = cosmic_count

            # Count ICGC-PCAWG hits by taking sum of affected donors where the annotation value has
            # the following format: project_code|tumor_type|affected_donors|tested_donors|frequency
            icgc_pcawg_count = 0
            if record['ICGC_PCAWG_OCCURRENCE'] != 'NA':
                for pcawg_hit_data in record['ICGC_PCAWG_OCCURRENCE'].split(','):
                    pcawrg_hit_data_fields = pcawg_hit_data.split('|')
                    affected_donors = int(pcawrg_hit_data_fields[2])
                    icgc_pcawg_count += affected_donors
                assert icgc_pcawg_count > 0
            data_rec[constants.VcfInfo.PCGR_ICGC_PCAWG_COUNT] = icgc_pcawg_count

            # Store annotation data for the record
            data[key] = data_rec

    return data
