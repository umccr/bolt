import click
import pathlib


import cyvcf2


from ... import util
from ...common import constants
from ...logging_config import setup_logging


@click.command(name='filter')
@click.pass_context

@click.option('--tumor_name', required=True, type=str)

@click.option('--vcf_fp', required=True, type=click.Path(exists=True))

@click.option('--output_dir', required=True, type=click.Path())

def entry(ctx, **kwargs):
    '''Set and apply filters for variants\f
    '''

    # Create output directory
    output_dir = pathlib.Path(kwargs['output_dir'])
    output_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

    script_name = pathlib.Path(__file__).stem
    setup_logging("logs", script_name)

    # Open input VCF and set required header entries for output
    in_fh = cyvcf2.VCF(kwargs['vcf_fp'])
    header_filters = (
        constants.VcfFilter.MIN_AF,
        constants.VcfFilter.MIN_AD,
        constants.VcfFilter.MIN_AD_DIFFICULT,
        constants.VcfFilter.MIN_AD_NON_GIAB,
        constants.VcfFilter.PON,
        constants.VcfFilter.ENCODE,
        constants.VcfFilter.GNOMAD_COMMON,
        constants.VcfInfo.SAGE_HOTSPOT_RESCUE,
        constants.VcfInfo.PCGR_TIER_RESCUE,
        constants.VcfInfo.CLINICAL_POTENTIAL_RESCUE,
        constants.VcfInfo.RESCUED_FILTERS_EXISTING,
        constants.VcfInfo.RESCUED_FILTERS_PENDING,
    )
    for header_enum in header_filters:
        util.add_vcf_header_entry(in_fh, header_enum)

    # Apply FILTERs and annotate with other INFO data
    filters_fp = output_dir / f'{kwargs["tumor_name"]}.filters_set.vcf.gz'
    filters_fh = cyvcf2.Writer(filters_fp, in_fh, 'wz')

    tumor_index = in_fh.samples.index(kwargs['tumor_name'])
    for record in in_fh:
        set_filter_data(record, tumor_index)
        filters_fh.write_record(record)
    filters_fh.close()

    # Apply set FILTERs
    set_fp = output_dir / f'{kwargs["tumor_name"]}.pass.vcf.gz'
    command = fr'''
        bcftools view -f PASS,. -o {set_fp} {filters_fp}
        bcftools index -t {set_fp}
    '''
    util.execute_command(command)


def set_filter_data(record, tumor_index):
    # NOTE(SW): given the importance of the filtering and rescue logic I've decided to keep it all
    # inline under a single function to avoid complicating abstractions


    ########################
    ## Variant filtering  ##
    ########################
    filters = list()

    tumor_af = record.format('AF')[tumor_index,0]
    tumor_ad = record.format('AD')[tumor_index,1]

    ##
    # AF filter
    ##
    if tumor_af < constants.MIN_AF:
        filters.append(constants.VcfFilter.MIN_AF)

    ##
    # AD filter (general)
    ##
    if tumor_ad < constants.MIN_AD:
        filters.append(constants.VcfFilter.MIN_AD)

    ##
    # AD filter (degraded mappability)
    ##
    # If a variant falls within difficult to call regions (low sequence complexity, poor
    # mappability, segemental duplications, etc), increase required minimum allele depth in anticipation of
    # elevated alignment error rate
    difficult_region_tags_enums = (
        constants.VcfInfo.DIFFICULT_BAD_PROMOTER,
        constants.VcfInfo.DIFFICULT_GC15,
        constants.VcfInfo.DIFFICULT_GC70TO75,
        constants.VcfInfo.DIFFICULT_GC75TO80,
        constants.VcfInfo.DIFFICULT_GC80TO85,
        constants.VcfInfo.DIFFICULT_GC80,
        constants.VcfInfo.DIFFICULT_LOW_COMPLEXITY_DITR,
        constants.VcfInfo.DIFFICULT_LOW_COMPLEXITY_QUADTR,
        constants.VcfInfo.DIFFICULT_LOW_COMPLEXITY_TANDEMREPEATS,
        constants.VcfInfo.DIFFICULT_LOW_COMPLEXITY_TRITR,
        constants.VcfInfo.DIFFICULT_MAPPABILITY_NONUNIQUE,
        constants.VcfInfo.DIFFICULT_SEGDUP,
    )
    difficult_region_tags = {e.value for e in difficult_region_tags_enums}

    if tumor_ad < constants.MIN_AD_DIFFICULT_REGIONS:

        if any(record.INFO.get(e) is not None for e in difficult_region_tags):
            filters.append(constants.VcfFilter.MIN_AD_DIFFICULT)

        if record.INFO.get(constants.VcfInfo.GIAB_CONF.value) is None:
            filters.append(constants.VcfFilter.MIN_AD_NON_GIAB)

    # NOTE(SW): filter_somatic_vcf from umccr/vcf_stuff includes a mappability filter but the INFO
    # field used to evaluate does not exist in the input annotated VCF, so is not included

    ##
    # PON filter
    ##
    # NOTE(SW): 'max' is inclusive - keeps variants with 0 to n-1 PON hits; preserved from Umccrise
    pon_count = record.INFO.get(constants.VcfInfo.PON_COUNT.value, 0)
    if pon_count >= constants.PON_HIT_THRESHOLD:
        filters.append(constants.VcfFilter.PON)

    ##
    # ENCODE blocklist filter
    ##
    if record.INFO.get(constants.VcfInfo.ENCODE.value) is not None:
        filters.append(constants.VcfFilter.ENCODE)

    ##
    # Common population variant filter
    ##
    # NOTE(SW): rounding is essential here for accurate comparison; cyvcf2 floating-point error
    # means INFO/gnomAD_AF=0.01 can be represented as 0.009999999776482582
    gnomad_af = round(record.INFO.get(constants.VcfInfo.GNOMAD_AF.value, 0), 3)
    if gnomad_af >= constants.MAX_GNOMAD_AF:
        filters.append(constants.VcfFilter.GNOMAD_COMMON)


    ######################
    ##  Variant rescue  ##
    ######################
    # Attempt to rescue variants that would otherwise be filtered downstream

    # NOTE(SW): variants are only checked if they can be rescued (i.e. are set to be filtered)
    # later to improve readability below i.e. we may get rescue information for variants that are
    # already going to pass all filters

    # NOTE(SW): the logic to annotate rescued variants below requires that info_rescue is populated
    # with at least one entry for each successful rescue test
    info_rescue = list()

    ##
    # PCGR tier rescue
    ##
    pcgr_tier = record.INFO.get(constants.VcfInfo.PCGR_TIER.value)
    if pcgr_tier in constants.PCGR_TIERS_RESCUE:
        info_rescue.append(constants.VcfInfo.PCGR_TIER_RESCUE)

    ##
    # SAGE hotspot rescue
    ##
    # NOTE(SW): effectively reverts any FILTERs that may have been applied above
    if record.INFO.get(constants.VcfInfo.SAGE_HOTSPOT.value) is not None:
        info_rescue.append(constants.VcfInfo.SAGE_HOTSPOT_RESCUE)

    ##
    # Clinical potential rescue; hotspot, driver, otherwise known
    ##

    # TODO(SW): split these into driver, hotspot, and clinical significance; currently collapsed as
    # single CLINICAL_POTENTIAL_RESCUE flag

    # Get ClinVar clinical significance entries
    clinvar_clinsig = record.INFO.get(constants.VcfInfo.PCGR_CLINVAR_CLNSIG.value, '')
    clinvar_clinsigs = clinvar_clinsig.split(',')
    # Hit counts in relevant reference somatic mutation databases
    cosmic_count = record.INFO.get(constants.VcfInfo.PCGR_COSMIC_COUNT.value, 0)
    tcga_pancancer_count = record.INFO.get(constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT.value, 0)
    icgc_pcawg_count = record.INFO.get(constants.VcfInfo.PCGR_ICGC_PCAWG_COUNT.value, 0)
    if (
        record.INFO.get(constants.VcfInfo.HMF_HOTSPOT.value) is not None or
        record.INFO.get(constants.VcfInfo.PCGR_MUTATION_HOTSPOT.value) is not None or
        any(e in clinvar_clinsigs for e in constants.CLINVAR_CLINSIGS_RESCUE) or
        cosmic_count >= constants.MIN_COSMIC_COUNT_RESCUE or
        tcga_pancancer_count >= constants.MIN_TCGA_PANCANCER_COUNT_RESCUE or
        icgc_pcawg_count >= constants.MIN_ICGC_PCAWG_COUNT_RESCUE
    ):
        info_rescue.append(constants.VcfInfo.CLINICAL_POTENTIAL_RESCUE)

    ##
    # Apply valid rescues
    ##
    # Valid rescue occurs when a variant that would be excluded by either existing or new filters
    # meets at least one rescue criteria. A variant is rescued by setting FILTER=PASS and adding
    # rescue-specific INFO flags.
    rescue_candidate = filters or record.FILTER
    rescued_variant = bool(info_rescue)
    if rescue_candidate and rescued_variant:
        # Set rescue info
        for info_enum in info_rescue:
            assert record.INFO.get(info_enum.value) is None
            record.INFO[info_enum.value] = True
        # Add rescued filters
        if record.FILTER:
            record.INFO['RESCUED_FILTERS_EXISTING'] = record.FILTER.replace(';', ',')
        if filters:
            record.INFO['RESCUED_FILTERS_PENDING'] = ','.join(sorted(f.value for f in filters))
        # Clear filters
        record.FILTER = 'PASS'
        filters = list()


    #######################
    ## Apply new filters ##
    #######################
    if filters:
        filters_value = {e.value for e in filters}
        filters_existing = [e for e in record.FILTERS if e != 'PASS']
        assert all(e not in filters_existing for e in filters_value)
        record.FILTER = ';'.join([*filters_existing, *filters_value])
