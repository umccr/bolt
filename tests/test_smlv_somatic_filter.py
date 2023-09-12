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
        self.record_min_af_filtered_data= {
            'format_data': {'GT': '0/1', 'AD': '100,10', 'AF': '0.010'},
        }


    def test_min_af_filter(self):
        record = get_record(**self.record_min_af_filtered_data)
        smlv_somatic_filter.set_filter_data(record, 0)

        assert record.FILTER is not None
        assert record.FILTERS == [bolt_constants.VcfFilter.MIN_AF.value]


    def test_clinvar_clinsig_rescue(self):
        clinsigs = [
            'conflicting_interpretations_of_pathogenicity',
            'likely_pathogenic',
            'pathogenic',
            'uncertain_significance',
        ]
        rescue_tag_str = bolt_constants.VcfInfo.CLINICAL_POTENTIAL_RESCUE.value

        for clinsig in clinsigs:
            record = get_record(
                **self.record_min_af_filtered_data,
                info_data={'PCGR_CLINVAR_CLNSIG': clinsig},
            )
            smlv_somatic_filter.set_filter_data(record, 0)

            assert record.FILTER is None
            assert record.INFO.get(rescue_tag_str) is not None
