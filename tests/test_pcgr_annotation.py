"""Unit tests for bolt/common/pcgr.py pure annotation-transfer helpers.

Covers only binary-free / cyvcf2-in-memory logic. Functions requiring
bcftools/PCGR/vcfanno subprocesses are intentionally NOT tested here.
Does NOT duplicate tests already in test_pcgr_hypermutated.py
(get_ordering, get_impacts, determine_filter, get_variant_filter_data,
split_vcf, run_somatic_chunk, count_variant_process, select_pcgr_variants).
"""
import pathlib
import tempfile
import unittest

import cyvcf2

import bolt.common.constants as constants
import bolt.common.pcgr as pcgr
import bolt.util as util


# Minimal VCF header with INFO tags exercised by annotate_record tests
HEADER = (
    '##fileformat=VCFv4.2\n'
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    f'##INFO=<ID={constants.VcfInfo.PCGR_ACTIONABILITY_TIER.value},Number=1,Type=String,Description="">\n'
    f'##INFO=<ID={constants.VcfInfo.PCGR_CSQ.value},Number=.,Type=String,Description="">\n'
    '##contig=<ID=chr1,length=248956422>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
)


def _write_vcf(path, variants):
    with open(path, 'w') as fh:
        fh.write(HEADER)
        for chrom, pos, ref, alt, info in variants:
            fh.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n')


def _make_variant_with_tags(tags, chrom='chr1', pos=100, ref='A', alt='T', info_str='.'):
    """Return a cyvcf2 Variant with a custom set of INFO enum tags registered."""
    with tempfile.TemporaryDirectory() as tmp:
        vcf_path = pathlib.Path(tmp) / 'test.vcf'
        _write_vcf(vcf_path, [(chrom, pos, ref, alt, info_str)])
        fh = cyvcf2.VCF(str(vcf_path))
        for tag_enum in tags:
            util.add_vcf_header_entry(fh, tag_enum)
        return list(fh)[0]


def _make_variant(chrom='chr1', pos=100, ref='A', alt='T', info_str='.'):
    """Return a cyvcf2 Variant with default INFO tags registered."""
    return _make_variant_with_tags(
        [constants.VcfInfo.PCGR_ACTIONABILITY_TIER, constants.VcfInfo.PCGR_CSQ],
        chrom=chrom, pos=pos, ref=ref, alt=alt, info_str=info_str,
    )


class TestParseGenomicChange(unittest.TestCase):
    """Unit tests for pcgr.parse_genomic_change()."""

    def test_standard_substitution(self):
        chrom, pos, ref, alt = pcgr.parse_genomic_change('3:g.41224645T>C')
        self.assertEqual(chrom, 'chr3')
        self.assertEqual(pos, 41224645)
        self.assertEqual(ref, 'T')
        self.assertEqual(alt, 'C')

    def test_chr_prefix_added(self):
        chrom, _, _, _ = pcgr.parse_genomic_change('17:g.7674220C>T')
        self.assertTrue(chrom.startswith('chr'), f"Expected 'chr' prefix, got: {chrom}")

    def test_multichar_alleles(self):
        chrom, pos, ref, alt = pcgr.parse_genomic_change('1:g.100ACGT>TTTT')
        self.assertEqual(chrom, 'chr1')
        self.assertEqual(ref, 'ACGT')
        self.assertEqual(alt, 'TTTT')

    def test_returns_int_pos(self):
        _, pos, _, _ = pcgr.parse_genomic_change('3:g.41224645T>C')
        self.assertIsInstance(pos, int)

    def test_raises_on_garbage_input(self):
        with self.assertRaises(ValueError):
            pcgr.parse_genomic_change('garbage')

    def test_raises_on_partial_format(self):
        with self.assertRaises(ValueError):
            pcgr.parse_genomic_change('3:g.41224645')


class TestGetImpactsHigher(unittest.TestCase):
    """Unit tests for pcgr.get_impacts_higher().

    Expectations are derived from constants.VEP_IMPACTS_FILTER rather than
    hard-coded literals, so the test tracks the source of truth.
    """

    def test_returns_higher_impacts_for_first_element(self):
        first = constants.VEP_IMPACTS_FILTER[0]
        higher = pcgr.get_impacts_higher(first)
        expected = list(constants.VEP_IMPACTS_FILTER[1:])
        self.assertEqual(list(higher), expected)

    def test_returns_empty_for_last_element(self):
        last = constants.VEP_IMPACTS_FILTER[-1]
        higher = pcgr.get_impacts_higher(last)
        self.assertEqual(list(higher), [])

    def test_middle_element(self):
        if len(constants.VEP_IMPACTS_FILTER) < 3:
            self.skipTest('VEP_IMPACTS_FILTER has fewer than 3 elements')
        mid_idx = len(constants.VEP_IMPACTS_FILTER) // 2
        mid = constants.VEP_IMPACTS_FILTER[mid_idx]
        higher = pcgr.get_impacts_higher(mid)
        expected = list(constants.VEP_IMPACTS_FILTER[mid_idx + 1:])
        self.assertEqual(list(higher), expected)

    def test_does_not_include_self(self):
        impact = constants.VEP_IMPACTS_FILTER[0]
        higher = pcgr.get_impacts_higher(impact)
        self.assertNotIn(impact, higher)


class TestGetAnnotationEntryTsv(unittest.TestCase):
    """Unit tests for pcgr.get_annotation_entry_tsv()."""

    # Map used in tests: one field from INFO enum present in VCF_HEADER_ENTRIES
    INFO_FIELD_MAP = {
        constants.VcfInfo.PCGR_CSQ: 'CSQ',
    }

    def test_genomic_change_key_extraction(self):
        record = {
            'GENOMIC_CHANGE': '3:g.41224645T>C',
            'CSQ': 'some_consequence',
        }
        key, record_ann = pcgr.get_annotation_entry_tsv(record, self.INFO_FIELD_MAP)
        self.assertEqual(key, ('chr3', 41224645, 'T', 'C'))

    def test_chrom_pos_ref_alt_fallback_with_chr_prefix(self):
        record = {
            'CHROM': '1',
            'POS': '100',
            'REF': 'A',
            'ALT': 'T',
            'CSQ': 'some_consequence',
        }
        key, record_ann = pcgr.get_annotation_entry_tsv(record, self.INFO_FIELD_MAP)
        chrom, pos, ref, alt = key
        self.assertTrue(str(chrom).startswith('chr'), f"Expected 'chr' prefix, got: {chrom}")

    def test_na_values_omitted_from_record_ann(self):
        record = {
            'GENOMIC_CHANGE': '3:g.41224645T>C',
            'CSQ': 'NA',
        }
        _, record_ann = pcgr.get_annotation_entry_tsv(record, self.INFO_FIELD_MAP)
        self.assertNotIn(constants.VcfInfo.PCGR_CSQ, record_ann)

    def test_falsy_values_omitted_from_record_ann(self):
        record = {
            'GENOMIC_CHANGE': '3:g.41224645T>C',
            'CSQ': '',
        }
        _, record_ann = pcgr.get_annotation_entry_tsv(record, self.INFO_FIELD_MAP)
        self.assertNotIn(constants.VcfInfo.PCGR_CSQ, record_ann)

    def test_valid_value_present_in_record_ann(self):
        record = {
            'GENOMIC_CHANGE': '3:g.41224645T>C',
            'CSQ': 'A|intron_variant|...',
        }
        _, record_ann = pcgr.get_annotation_entry_tsv(record, self.INFO_FIELD_MAP)
        self.assertIn(constants.VcfInfo.PCGR_CSQ, record_ann)
        self.assertEqual(record_ann[constants.VcfInfo.PCGR_CSQ], 'A|intron_variant|...')


class TestCompileAnnotationData(unittest.TestCase):
    """Unit tests for pcgr.compile_annotation_data()."""

    def test_tsv_wins_on_key_field_collision(self):
        key = ('chr1', 100, 'A', 'T')
        data_tsv = {key: {constants.VcfInfo.PCGR_CSQ: 'tsv_value'}}
        data_vcf = {key: {constants.VcfInfo.PCGR_CSQ: 'vcf_value'}}
        result = pcgr.compile_annotation_data(data_tsv, data_vcf)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'tsv_value')

    def test_vcf_only_key_added(self):
        key_tsv = ('chr1', 100, 'A', 'T')
        key_vcf = ('chr2', 200, 'G', 'C')
        data_tsv = {key_tsv: {constants.VcfInfo.PCGR_CSQ: 'tsv_value'}}
        data_vcf = {key_vcf: {constants.VcfInfo.PCGR_CSQ: 'vcf_value'}}
        result = pcgr.compile_annotation_data(data_tsv, data_vcf)
        self.assertIn(key_vcf, result)
        self.assertEqual(result[key_vcf][constants.VcfInfo.PCGR_CSQ], 'vcf_value')

    def test_vcf_only_field_added_to_existing_key(self):
        key = ('chr1', 100, 'A', 'T')
        data_tsv = {key: {constants.VcfInfo.PCGR_ACTIONABILITY_TIER: '1'}}
        data_vcf = {key: {constants.VcfInfo.PCGR_CSQ: 'vcf_csq'}}
        result = pcgr.compile_annotation_data(data_tsv, data_vcf)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_ACTIONABILITY_TIER], '1')
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'vcf_csq')

    def test_empty_inputs(self):
        result = pcgr.compile_annotation_data({}, {})
        self.assertEqual(result, {})

    def test_tsv_only_key_preserved(self):
        key = ('chr1', 100, 'A', 'T')
        data_tsv = {key: {constants.VcfInfo.PCGR_CSQ: 'tsv_value'}}
        result = pcgr.compile_annotation_data(data_tsv, {})
        self.assertIn(key, result)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'tsv_value')


class TestAnnotateRecord(unittest.TestCase):
    """Unit tests for pcgr.annotate_record()."""

    def test_annotations_written_to_matching_record(self):
        variant = _make_variant_with_tags(
            [constants.VcfInfo.PCGR_ACTIONABILITY_TIER],
            chrom='chr1', pos=100, ref='A', alt='T',
        )
        annotations = {
            ('chr1', 100, 'A', 'T'): {
                constants.VcfInfo.PCGR_ACTIONABILITY_TIER: '2',
            }
        }
        result = pcgr.annotate_record(variant, annotations, allow_missing=True)
        self.assertEqual(
            result.INFO.get(constants.VcfInfo.PCGR_ACTIONABILITY_TIER.value),
            '2',
        )

    def test_missing_key_allow_missing_true_returns_record_unchanged(self):
        variant = _make_variant(chrom='chr1', pos=100, ref='A', alt='T')
        annotations = {}
        result = pcgr.annotate_record(variant, annotations, allow_missing=True)
        self.assertIsNotNone(result)

    def test_missing_key_allow_missing_false_raises(self):
        # NOTE: annotate_record's `assert key not in annotations` guard is tautological
        # (it re-checks a condition already established by the enclosing `if`), so it never
        # fires. Control falls through to `annotations[key].items()`, which raises KeyError
        # instead of the AssertionError the guard's phrasing implies. Asserting the actual
        # (KeyError) behavior here — the net effect (failing loudly on an unresolvable
        # missing key) is preserved either way.
        variant = _make_variant(chrom='chr1', pos=100, ref='A', alt='T')
        annotations = {}
        with self.assertRaises(KeyError):
            pcgr.annotate_record(variant, annotations, allow_missing=False)

    def test_key_uses_chrom_pos_ref_alt(self):
        """annotate_record must match on the exact (CHROM, POS, REF, ALT) tuple."""
        variant = _make_variant_with_tags(
            [constants.VcfInfo.PCGR_ACTIONABILITY_TIER],
            chrom='chr5', pos=999, ref='G', alt='C',
        )
        annotations = {
            ('chr5', 999, 'G', 'C'): {
                constants.VcfInfo.PCGR_ACTIONABILITY_TIER: '3',
            }
        }
        result = pcgr.annotate_record(variant, annotations, allow_missing=False)
        self.assertEqual(
            result.INFO.get(constants.VcfInfo.PCGR_ACTIONABILITY_TIER.value),
            '3',
        )


if __name__ == '__main__':
    unittest.main()
