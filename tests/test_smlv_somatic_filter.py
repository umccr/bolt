import tempfile
import unittest


import cyvcf2


import bolt.workflows.smlv_somatic.filter as smlv_somatic_filter
import bolt.common.constants as bolt_constants
import bolt.util as bolt_util


# TODO(SW): place this helper code to obtain cyvcf2 Variant classes somewhere else
HEADER_STR = (
    '##fileformat=VCFv4.2\n'
    '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="">\n'
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n'
    '##INFO=<ID=GIAB_CONF,Number=0,Type=Flag,Description="">\n'
    '##INFO=<ID=DIFFICULT_segdup,Number=0,Type=Flag,Description="">\n'
    '##INFO=<ID=PON_COUNT,Number=1,Type=Integer,Description="">\n'
    '##INFO=<ID=ENCODE,Number=0,Type=Flag,Description="">\n'
    '##INFO=<ID=HMF_HOTSPOT,Number=0,Type=Flag,Description="">\n'
    '##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="">\n'
    '##contig=<ID=chr1,length=248956422>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'
)


# Create VCF with required header and provided variant
def create_vcf(variant_str):
    fh = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    fh.write(HEADER_STR)
    fh.write(variant_str)
    fh.close()
    return fh.name

# Read in as cyvcf2 Variant record; first add all possible header entries
def get_record_from_str(variant_str):
    fp = create_vcf(variant_str)
    fh = cyvcf2.VCF(fp)
    for header_enum in bolt_constants.VCF_HEADER_ENTRIES:
        bolt_util.add_vcf_header_entry(fh, header_enum)
    return next(fh)


def get_record(
    chrom='chr1',
    pos='.',
    vid='.',
    ref='.',
    alt='.',
    qual='.',
    vfilter='.',
    info_data=None,
    format_data=None
):
    if info_data is None:
        info_data = dict()

    if format_data is None:
        format_data = {'GT': '0/0'}

    info_tokens = list()
    for k, v in info_data.items():
        info_tokens.append(f'{k}={v}' if v else k)
    info_str = ';'.join(info_tokens)

    format_str = ':'.join(format_data.keys())
    sample_str = ':'.join(format_data.values())

    variant_cmps = [chrom, pos, vid, ref, alt, qual, vfilter, info_str, format_str, sample_str]
    variant_str = '\t'.join(variant_cmps)

    return get_record_from_str(variant_str)


class TestSmlvSomaticFilter(unittest.TestCase):

    def setUp(self):
        self.records = {
            'pass_af10':        {'format_data': {'GT': '0/1', 'AD': '54,6',  'AF': '0.100'}},
            'pass_af20':        {'format_data': {'GT': '0/1', 'AD': '48,12', 'AF': '0.200'}},
            'filter_min_af9.9': {'format_data': {'GT': '0/1', 'AD': '91,10', 'AF': '0.099'}},
            'filter_min_ad3':   {'format_data': {'GT': '0/1', 'AD': '27,3',  'AF': '0.100'}},
            'filter_min_ad5':   {'format_data': {'GT': '0/1', 'AD': '45,5',  'AF': '0.100'}},
        }




    def test_min_af_filter(self):
        record = get_record(**self.records['filter_min_af9.9'])
        smlv_somatic_filter.set_filter_data(record, 0)
        assert record.FILTERS == [bolt_constants.VcfFilter.MIN_AF.value]


    def test_min_ad_filter(self):
        record = get_record(
            **self.records['filter_min_ad3'],
            info_data={'GIAB_CONF': ''},
        )
        smlv_somatic_filter.set_filter_data(record, 0)
        assert bolt_constants.VcfFilter.MIN_AD.value in record.FILTERS


    def test_min_ad_difficult_filter(self):
        record = get_record(
            **self.records['filter_min_ad5'],
            info_data={'DIFFICULT_segdup': '', 'GIAB_CONF': ''},
        )
        smlv_somatic_filter.set_filter_data(record, 0)

        assert record.FILTERS == [bolt_constants.VcfFilter.MIN_AD_DIFFICULT.value]


    def test_min_ad_difficult_filter_with_non_giab_conf(self):
        record = get_record(
            **self.records['filter_min_ad3'],
            info_data={'DIFFICULT_segdup': ''},
        )
        smlv_somatic_filter.set_filter_data(record, 0)

        filters_expected = [
            bolt_constants.VcfFilter.MIN_AD,
            bolt_constants.VcfFilter.MIN_AD_DIFFICULT,
            bolt_constants.VcfFilter.MIN_AD_NON_GIAB,
        ]
        assert all(fe.value in record.FILTERS for fe in filters_expected)


    def test_pon_filter(self):
        record = get_record(
            **self.records['pass_af10'],
            info_data={'PON_COUNT': 5},
        )
        smlv_somatic_filter.set_filter_data(record, 0)
        assert record.FILTERS == [bolt_constants.VcfFilter.PON.value]


    def test_encode_blocklist_filter(self):
        record = get_record(
            **self.records['pass_af10'],
            info_data={'ENCODE': ''},
        )
        smlv_somatic_filter.set_filter_data(record, 0)
        assert record.FILTERS == [bolt_constants.VcfFilter.ENCODE.value]


    def test_common_population_filter(self):
        record = get_record(
            **self.records['pass_af10'],
            info_data={'gnomAD_AF': 0.01},
        )
        smlv_somatic_filter.set_filter_data(record, 0)
        assert record.FILTERS == [bolt_constants.VcfFilter.GNOMAD_COMMON.value]




    def test_pcgr_tier_rescue(self):
        pcgr_tiers = [
            'TIER_1',
            'TIER_2',
        ]
        rescue_tag_str = bolt_constants.VcfInfo.PCGR_TIER_RESCUE.value

        for pcgr_tier in pcgr_tiers:
            record = get_record(
                **self.records['filter_min_af9.9'],
                info_data={'PCGR_TIER': pcgr_tier},
            )
            smlv_somatic_filter.set_filter_data(record, 0)
            assert not record.FILTER
            assert record.INFO.get(rescue_tag_str) is not None


    def test_sage_hotspot_rescue(self):
        record = get_record(
            **self.records['filter_min_af9.9'],
            info_data={'SAGE_HOTSPOT': ''},
        )
        smlv_somatic_filter.set_filter_data(record, 0)
        assert not record.FILTER
        assert record.INFO.get(bolt_constants.VcfInfo.SAGE_HOTSPOT_RESCUE.value) is not None


    def test_clinical_potential_rescue_general(self):
        info_data_sets = [
            {'HMF_HOTSPOT': ''},
            {'PCGR_MUTATION_HOTSPOT': ''},
            {'PCGR_COSMIC_COUNT': 11},
            {'PCGR_TCGA_PANCANCER_COUNT': 6},
            {'PCGR_ICGC_PCAWG_COUNT': 6},
        ]
        rescue_tag_str = bolt_constants.VcfInfo.CLINICAL_POTENTIAL_RESCUE.value

        for info_data_set in info_data_sets:
            record = get_record(
                **self.records['filter_min_af9.9'],
                info_data=info_data_set,
            )
            smlv_somatic_filter.set_filter_data(record, 0)
            assert not record.FILTER
            assert record.INFO.get(rescue_tag_str) is not None


    def test_clinical_potential_rescue_clinvar_clinsig(self):
        clinsigs = [
            'conflicting_interpretations_of_pathogenicity',
            'likely_pathogenic',
            'pathogenic',
            'uncertain_significance',
        ]
        rescue_tag_str = bolt_constants.VcfInfo.CLINICAL_POTENTIAL_RESCUE.value

        for clinsig in clinsigs:
            record = get_record(
                **self.records['filter_min_af9.9'],
                info_data={'PCGR_CLINVAR_CLNSIG': clinsig},
            )
            smlv_somatic_filter.set_filter_data(record, 0)
            assert not record.FILTER
            assert record.INFO.get(rescue_tag_str) is not None




    # NOTE(SW): ensuring that a passing variant with SAGE_HOTSPOT etc isn't tagged as being rescued
    def test_unmodified_sage_hotspot(self):
        record = get_record(
            **self.records['pass_af10'],
            info_data={'SAGE_HOTSPOT': ''},
        )
        smlv_somatic_filter.set_filter_data(record, 0)
        assert not record.FILTER
        assert record.INFO.get(bolt_constants.VcfInfo.SAGE_HOTSPOT_RESCUE.value) is None


    def test_stored_rescue_filters__pending_filters(self):
        filters_pending_expected = [
            bolt_constants.VcfFilter.MIN_AD,
            bolt_constants.VcfFilter.MIN_AD_NON_GIAB,
        ]
        filters_pending_expected_str = ','.join(sorted(f.value for f in filters_pending_expected))
        rescued_filters_str = bolt_constants.VcfInfo.RESCUED_FILTERS_PENDING.value

        for vfilter in ['.', 'PASS']:
            record = get_record(
                **self.records['filter_min_ad3'],
                vfilter=vfilter,
                info_data={'HOTSPOT': ''},
            )
            smlv_somatic_filter.set_filter_data(record, 0)

            assert not record.FILTER
            assert record.INFO.get(bolt_constants.VcfInfo.RESCUED_FILTERS_EXISTING.value) is None
            assert record.INFO.get(rescued_filters_str) == filters_pending_expected_str


    def test_stored_rescue_filters__existing_filter(self):
        rescued_filters_str = bolt_constants.VcfInfo.RESCUED_FILTERS_EXISTING.value

        record = get_record(
            **self.records['pass_af10'],
            vfilter='SEGDUP',
            info_data={'HOTSPOT': ''},
        )
        smlv_somatic_filter.set_filter_data(record, 0)

        assert not record.FILTER
        assert record.INFO.get(bolt_constants.VcfInfo.RESCUED_FILTERS_PENDING.value) is None
        assert record.INFO.get(rescued_filters_str) == 'SEGDUP'
