import click
import pathlib


import cyvcf2


from ... import util
from ...common import constants


@click.command(name='filter')
@click.pass_context
@click.option('--vcf_fp', required=True, type=click.Path(exists=True))

@click.option('--tumor_name', required=True, type=str)

def entry(ctx, **kwargs):
    '''Set and apply filters for variants\f
    '''
    in_fh = cyvcf2.VCF(kwargs['vcf_fp'])

    header_filters = (
        constants.VcfFilter.MIN_AF,
        constants.VcfFilter.MIN_AD,
        constants.VcfFilter.MIN_AD_DIFFICULT,
        constants.VcfFilter.MIN_AD_NON_GIAB,
        constants.VcfFilter.PON,
        constants.VcfFilter.ENCODE,
        constants.VcfFilter.GNOMAD_COMMON,
        constants.VcfInfo.PCGR_TIER_RESCUE,
        constants.VcfInfo.CLINICAL_POTENTIAL_RESCUE,
        constants.VcfInfo.GERMLINE_LEAKAGE,
    )
    for header_enum in header_filters:
        util.add_vcf_header_entry(in_fh, header_enum)

    vcf_fn = pathlib.Path(kwargs['vcf_fp']).name
    vcf_prefix = vcf_fn.replace('.vcf.gz', '')

    filters_fp = f'{vcf_prefix}.filters_set.vcf.gz'
    filters_fh = cyvcf2.Writer(filters_fp, in_fh, 'wz')

    tumor_index = in_fh.samples.index(kwargs['tumor_name'])
    for record in in_fh:
        set_filter_data(record, tumor_index)
        filters_fh.write_record(record)
    filters_fh.close()

    # Apply set filters
    command = f'bcftools view -f PASS,. -o {vcf_prefix}.pass.vcf.gz {filters_fp}'
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
    # mappability, GIAB confidence, etc), increase required minimum allele depth in anticipation of
    # elevated alignment error rate
    difficult_region_tags_enums = (
        constants.VcfInfo.SEGDUP,
        constants.VcfInfo.TRICKY_LCR,
        constants.VcfInfo.TRICKY_GC15,
        constants.VcfInfo.TRICKY_GC70TO75,
        constants.VcfInfo.TRICKY_GC75TO80,
        constants.VcfInfo.TRICKY_GC80TO85,
        constants.VcfInfo.TRICKY_GC85,
        constants.VcfInfo.TRICKY_HENG_UM75_HS37D5,
        constants.VcfInfo.TRICKY_LOW_COMPLEXITY_51TO200BP,
        constants.VcfInfo.TRICKY_LOW_COMPLEXITY_GT200BP,
    )
    difficult_region_tags = {e.value for e in difficult_region_tags_enums}

    if tumor_ad < constants.MIN_AD_DIFFICULT_REGIONS:

        if any(record.INFO.get(e) for e in difficult_region_tags):
            filters.append(constants.VcfFilter.MIN_AD_DIFFICULT)

        if not record.INFO.get(constants.VcfInfo.GIAB_CONF.value):
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
    if record.INFO.get(constants.VcfInfo.ENCODE.value):
        filters.append(constants.VcfFilter.ENCODE)

    ##
    # Common population variant filter
    ##
    # NOTE(SW): gnomAD frequency is only checked in Umccrise (and still here) when no other FILTERs
    # are present so that it simplifies logic below for rescue - variants with only
    # FILTER/gnomAD_common are considered germline. I'd like to apply this filter to all variants
    # then check during germline leakage whether FILTER/gnomAD_common is the only FILTER set to
    # improve coherence of processing logic.
    if not record.FILTER and record.INFO.get(constants.VcfInfo.GNOMAD_AF.value, 0) >= constants.MAX_GNOMAD_AF:
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

    ###
    # SAGE hotspot rescue
    ###
    # NOTE(SW): effectively reverts any FILTERs that may have been applied above
    if record.INFO.get(constants.VcfInfo.SAGE_HOTSPOT.value):
        info_rescue.append(constants.VcfInfo.SAGE_HOTSPOT)

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
        record.INFO.get(constants.VcfInfo.HOTSPOT.value) or
        record.INFO.get(constants.VcfInfo.PCGR_MUTATION_HOTSPOT.value) or
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
        # Clear existing filters
        record.FILTER = 'PASS'
        filters = list()
        # Set rescue info
        for info_enum in info_rescue:
            assert not record.INFO.get(info_enum.value)
            record.INFO[info_enum.value] = True


    ######################
    ## Germline leakage ##
    ######################
    # Annotate variants thought to be germline leakage. Rescued variants are never treated as
    # germline leakage.

    # NOTE(SW): in Umccrise, where variants have set FILTER=PoN they are then immediately checked
    # if they are considered germline leakage but always fail the initial check for an empty
    # FILTER; germline leakage obtained via the gnomAD_AF filter is not affected by this issue

    # Germline leakage conditions
    # NOTE(SW): when implementing new conditions take care to include relevant filters in the below
    # conditional statement
    pon_gl = constants.VcfFilter.PON in filters and tumor_af >= constants.MIN_PON_GERMLINE_AF
    gnomad_gl = constants.VcfFilter.GNOMAD_COMMON in filters

    # Annotate variants where they:
    # 1. have not been rescued
    # 2. meet at least one germline leakage condition
    # 3. have no other filter set upstream
    if not rescued_variant and (pon_gl or gnomad_gl) and not record.FILTER:
        record.INFO[constants.VcfInfo.GERMLINE_LEAKAGE.value] = True


    #######################
    ## Apply new filters ##
    #######################
    if filters:
        filters_value = {e.value for e in filters}
        filters_existing = [e for e in record.FILTERS if e != 'PASS']
        assert all(e not in filters_existing for e in filters_value)
        record.FILTER = ';'.join([*filters_existing, *filters_value])
