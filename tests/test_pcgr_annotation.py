"""Unit tests for bolt/common/pcgr.py pure annotation-transfer helpers.

Covers only binary-free / cyvcf2-in-memory logic. Functions requiring
bcftools/PCGR/vcfanno subprocesses are intentionally NOT tested here.
Does NOT duplicate tests already in test_pcgr.py / test_smlv_somatic_report.py
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

    def test_missing_key_allow_missing_false_raises_assertion(self):
        variant = _make_variant(chrom='chr1', pos=100, ref='A', alt='T')
        annotations = {}
        with self.assertRaises(AssertionError):
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


class TestGetAnnotationsVcf(unittest.TestCase):
    """Unit tests for pcgr.get_annotations_vcf() duplicate-key handling.

    PCGR strips the 'chr' prefix from its own output VCF, so fixtures here use
    bare contig names ('1', not 'chr1') to match real PCGR output.
    """

    PCGR_VCF_HEADER = (
        '##fileformat=VCFv4.2\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        f'##INFO=<ID={constants.VcfInfo.PCGR_CSQ.value},Number=.,Type=String,Description="">\n'
        '##contig=<ID=1,length=248956422>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )

    def _write_pcgr_vcf(self, path, rows):
        with open(path, 'w') as fh:
            fh.write(self.PCGR_VCF_HEADER)
            for chrom, pos, ref, alt, info in rows:
                fh.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n')

    def test_duplicate_vcf_key_keeps_first_no_crash(self):
        """A duplicate variant in PCGR's VCF output must be skipped, not raise AssertionError."""
        info_field_map = {constants.VcfInfo.PCGR_CSQ: 'PCGR_CSQ'}
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'pcgr.vcf'
            self._write_pcgr_vcf(vcf_fp, [
                ('1', 100, 'A', 'T', 'PCGR_CSQ=first'),
                ('1', 100, 'A', 'T', 'PCGR_CSQ=second'),  # duplicate transcript mapping
            ])
            result = pcgr.get_annotations_vcf(vcf_fp, info_field_map)

        key = ('chr1', 100, 'A', 'T')
        self.assertEqual(len(result), 1)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'first')

    def test_non_duplicate_vcf_keys_all_present(self):
        """Distinct variants are all retained unaffected by dedup handling."""
        info_field_map = {constants.VcfInfo.PCGR_CSQ: 'PCGR_CSQ'}
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'pcgr.vcf'
            self._write_pcgr_vcf(vcf_fp, [
                ('1', 100, 'A', 'T', 'PCGR_CSQ=v1'),
                ('1', 200, 'C', 'G', 'PCGR_CSQ=v2'),
            ])
            result = pcgr.get_annotations_vcf(vcf_fp, info_field_map)

        self.assertEqual(len(result), 2)
        self.assertEqual(result[('chr1', 100, 'A', 'T')][constants.VcfInfo.PCGR_CSQ], 'v1')
        self.assertEqual(result[('chr1', 200, 'C', 'G')][constants.VcfInfo.PCGR_CSQ], 'v2')


def _write_empty_pcgr_vcf(path):
    """Header-only sites-only VCF, so get_annotations_vcf() returns {}."""
    with open(path, 'w') as fh:
        fh.write(
            '##fileformat=VCFv4.2\n'
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            f'##INFO=<ID={constants.VcfInfo.PCGR_CSQ.value},Number=.,Type=String,Description="">\n'
            '##contig=<ID=1,length=248956422>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        )


def _write_tsv(path, header, rows):
    import csv as _csv
    with open(path, 'w', newline='') as fh:
        writer = _csv.DictWriter(fh, fieldnames=header, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


class TestCollectPcgrAnnotationData(unittest.TestCase):
    """Unit tests for pcgr.collect_pcgr_annotation_data() duplicate-key tier resolution."""

    HEADER = ['GENOMIC_CHANGE', 'CSQ', 'ACTIONABILITY_TIER']

    def _collect(self, tmp, rows):
        tsv_fp = pathlib.Path(tmp) / 'pcgr.tsv'
        vcf_fp = pathlib.Path(tmp) / 'pcgr.vcf'
        _write_tsv(tsv_fp, self.HEADER, rows)
        _write_empty_pcgr_vcf(vcf_fp)
        info_field_map = {constants.VcfInfo.PCGR_CSQ: 'CSQ'}
        return pcgr.collect_pcgr_annotation_data(tsv_fp, vcf_fp, info_field_map)

    def test_duplicate_key_keeps_more_actionable_tier_no_crash(self):
        """A second, more-actionable-tier row for the same variant replaces the first."""
        with tempfile.TemporaryDirectory() as tmp:
            result = self._collect(tmp, [
                {'GENOMIC_CHANGE': '1:g.100A>T', 'CSQ': 'transcript_a', 'ACTIONABILITY_TIER': 'TIER 3'},
                {'GENOMIC_CHANGE': '1:g.100A>T', 'CSQ': 'transcript_b', 'ACTIONABILITY_TIER': 'TIER 1'},
            ])
        key = ('chr1', 100, 'A', 'T')
        self.assertEqual(len(result), 1)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'transcript_b')
        self.assertEqual(result[key][constants.VcfInfo.PCGR_ACTIONABILITY_TIER], '1')

    def test_duplicate_key_skips_less_actionable_tier(self):
        """A second, less-actionable-tier row for the same variant is dropped."""
        with tempfile.TemporaryDirectory() as tmp:
            result = self._collect(tmp, [
                {'GENOMIC_CHANGE': '1:g.100A>T', 'CSQ': 'transcript_a', 'ACTIONABILITY_TIER': 'TIER 1'},
                {'GENOMIC_CHANGE': '1:g.100A>T', 'CSQ': 'transcript_b', 'ACTIONABILITY_TIER': 'TIER 3'},
            ])
        key = ('chr1', 100, 'A', 'T')
        self.assertEqual(len(result), 1)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'transcript_a')
        self.assertEqual(result[key][constants.VcfInfo.PCGR_ACTIONABILITY_TIER], '1')

    def test_non_duplicate_keys_all_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self._collect(tmp, [
                {'GENOMIC_CHANGE': '1:g.100A>T', 'CSQ': 'v1', 'ACTIONABILITY_TIER': 'TIER 1'},
                {'GENOMIC_CHANGE': '1:g.200C>G', 'CSQ': 'v2', 'ACTIONABILITY_TIER': 'TIER 2'},
            ])
        self.assertEqual(len(result), 2)
        self.assertIn(('chr1', 100, 'A', 'T'), result)
        self.assertIn(('chr1', 200, 'C', 'G'), result)


class TestCollectCpsrAnnotationData(unittest.TestCase):
    """Unit tests for pcgr.collect_cpsr_annotation_data() duplicate-key handling."""

    HEADER = ['GENOMIC_CHANGE', 'CSQ']

    def _collect(self, tmp, rows):
        tsv_fp = pathlib.Path(tmp) / 'cpsr.tsv.gz'
        vcf_fp = pathlib.Path(tmp) / 'cpsr.vcf'
        import gzip as _gzip
        import csv as _csv
        with _gzip.open(tsv_fp, 'wt', newline='') as fh:
            writer = _csv.DictWriter(fh, fieldnames=self.HEADER, delimiter='\t')
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        _write_empty_pcgr_vcf(vcf_fp)
        info_field_map = {constants.VcfInfo.PCGR_CSQ: 'CSQ'}
        return pcgr.collect_cpsr_annotation_data(tsv_fp, vcf_fp, info_field_map)

    def test_duplicate_key_keeps_first_no_crash(self):
        """A duplicate CPSR TSV row for the same variant must be skipped, not raise."""
        with tempfile.TemporaryDirectory() as tmp:
            result = self._collect(tmp, [
                {'GENOMIC_CHANGE': '1:g.100A>T', 'CSQ': 'first'},
                {'GENOMIC_CHANGE': '1:g.100A>T', 'CSQ': 'second'},
            ])
        key = ('chr1', 100, 'A', 'T')
        self.assertEqual(len(result), 1)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'first')

    def test_non_duplicate_keys_all_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self._collect(tmp, [
                {'GENOMIC_CHANGE': '1:g.100A>T', 'CSQ': 'v1'},
                {'GENOMIC_CHANGE': '1:g.200C>G', 'CSQ': 'v2'},
            ])
        self.assertEqual(len(result), 2)
        self.assertIn(('chr1', 100, 'A', 'T'), result)
        self.assertIn(('chr1', 200, 'C', 'G'), result)


if __name__ == '__main__':
    unittest.main()
