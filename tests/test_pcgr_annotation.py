"""Unit tests for bolt/common/pcgr.py pure annotation-transfer helpers.

Covers only pure / cyvcf2-in-memory logic. Functions requiring bcftools/PCGR/vcfanno
subprocesses (prepare_vcf_*, run_somatic*, transfer_annotations_*, get_variant_filter_data,
split_vcf, run_somatic_chunk, get_ordering, get_impacts, determine_filter,
select_pcgr_variants, count_variant_process) are NOT re-tested here — see
tests/test_pcgr_hypermutated.py for those.
"""
import pathlib
import tempfile
import unittest

import cyvcf2

import bolt.common.constants as constants
import bolt.common.pcgr as pcgr
import bolt.util as util


# Minimal VCF header — declares the PCGR_CSQ INFO tag via util.add_vcf_header_entry-compatible
# constants so annotate_record can write values onto records built from it.
HEADER = (
    '##fileformat=VCFv4.2\n'
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##contig=<ID=chr3,length=198295559>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
)


def _write_vcf(path, variants):
    with open(path, 'w') as fh:
        fh.write(HEADER)
        for chrom, pos, ref, alt in variants:
            fh.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n')


def _make_variant(chrom='chr3', pos=41224645, ref='T', alt='C', target_infos=()):
    """Return a (variant, writer_handle) built from a minimal VCF header.

    target_infos: iterable of VcfInfo enums to register on the header before parsing,
    so INFO writes in annotate_record succeed.
    """
    tmp_dir = tempfile.mkdtemp()
    vcf_path = pathlib.Path(tmp_dir) / 'test.vcf'
    _write_vcf(vcf_path, [(chrom, pos, ref, alt)])
    fh = cyvcf2.VCF(str(vcf_path))
    for info_enum in target_infos:
        util.add_vcf_header_entry(fh, info_enum)
    return list(fh)[0]


class TestParseGenomicChange(unittest.TestCase):
    """Unit tests for pcgr.parse_genomic_change()."""

    def test_valid_snv(self):
        result = pcgr.parse_genomic_change('3:g.41224645T>C')
        self.assertEqual(result, ('chr3', 41224645, 'T', 'C'))

    def test_valid_snv_different_chrom(self):
        result = pcgr.parse_genomic_change('X:g.100T>A')
        self.assertEqual(result, ('chrX', 100, 'T', 'A'))

    def test_unparseable_raises_valueerror(self):
        with self.assertRaises(ValueError):
            pcgr.parse_genomic_change('garbage')


class TestGetImpactsHigher(unittest.TestCase):
    """Unit tests for pcgr.get_impacts_higher()."""

    def test_first_impact_returns_remaining_tail(self):
        # get_impacts_higher slices the source tuple directly — result is a tuple, not a list.
        first = constants.VEP_IMPACTS_FILTER[0]
        expected = constants.VEP_IMPACTS_FILTER[1:]
        self.assertEqual(pcgr.get_impacts_higher(first), expected)

    def test_last_impact_returns_empty_list(self):
        last = constants.VEP_IMPACTS_FILTER[-1]
        self.assertEqual(pcgr.get_impacts_higher(last), [])

    def test_middle_impact_returns_correct_tail(self):
        idx = len(constants.VEP_IMPACTS_FILTER) // 2
        impact = constants.VEP_IMPACTS_FILTER[idx]
        expected = constants.VEP_IMPACTS_FILTER[idx + 1:]
        self.assertEqual(pcgr.get_impacts_higher(impact), expected)


class TestGetAnnotationEntryTsv(unittest.TestCase):
    """Unit tests for pcgr.get_annotation_entry_tsv()."""

    def setUp(self):
        self.info_field_map = {constants.VcfInfo.PCGR_CSQ: 'CSQ'}

    def test_genomic_change_field_parsed(self):
        record = {'GENOMIC_CHANGE': '3:g.41224645T>C', 'CSQ': 'some_csq_value'}
        key, record_ann = pcgr.get_annotation_entry_tsv(record, self.info_field_map)
        self.assertEqual(key, ('chr3', 41224645, 'T', 'C'))
        self.assertEqual(record_ann[constants.VcfInfo.PCGR_CSQ], 'some_csq_value')

    def test_chrom_pos_ref_alt_fallback_no_chr_prefix(self):
        record = {'CHROM': '3', 'POS': '41224645', 'REF': 'T', 'ALT': 'C', 'CSQ': 'v'}
        key, _ = pcgr.get_annotation_entry_tsv(record, self.info_field_map)
        self.assertEqual(key, ('chr3', 41224645, 'T', 'C'))

    def test_chrom_already_prefixed_not_double_prefixed(self):
        record = {'CHROM': 'chr3', 'POS': '41224645', 'REF': 'T', 'ALT': 'C', 'CSQ': 'v'}
        key, _ = pcgr.get_annotation_entry_tsv(record, self.info_field_map)
        self.assertEqual(key[0], 'chr3')

    def test_na_value_omitted(self):
        record = {'GENOMIC_CHANGE': '3:g.41224645T>C', 'CSQ': 'NA'}
        _, record_ann = pcgr.get_annotation_entry_tsv(record, self.info_field_map)
        self.assertNotIn(constants.VcfInfo.PCGR_CSQ, record_ann)

    def test_falsy_value_omitted(self):
        record = {'GENOMIC_CHANGE': '3:g.41224645T>C', 'CSQ': ''}
        _, record_ann = pcgr.get_annotation_entry_tsv(record, self.info_field_map)
        self.assertNotIn(constants.VcfInfo.PCGR_CSQ, record_ann)


class TestCompileAnnotationData(unittest.TestCase):
    """Unit tests for pcgr.compile_annotation_data()."""

    def test_tsv_wins_on_field_collision(self):
        key = ('chr3', 41224645, 'T', 'C')
        data_tsv = {key: {constants.VcfInfo.PCGR_CSQ: 'from_tsv'}}
        data_vcf = {key: {constants.VcfInfo.PCGR_CSQ: 'from_vcf'}}
        result = pcgr.compile_annotation_data(data_tsv, data_vcf)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'from_tsv')

    def test_vcf_only_key_added(self):
        key_tsv = ('chr3', 1, 'A', 'T')
        key_vcf = ('chr3', 2, 'A', 'T')
        data_tsv = {key_tsv: {constants.VcfInfo.PCGR_CSQ: 'tsv_val'}}
        data_vcf = {key_vcf: {constants.VcfInfo.PCGR_CSQ: 'vcf_val'}}
        result = pcgr.compile_annotation_data(data_tsv, data_vcf)
        self.assertIn(key_vcf, result)
        self.assertEqual(result[key_vcf][constants.VcfInfo.PCGR_CSQ], 'vcf_val')

    def test_vcf_only_field_added_to_existing_key(self):
        key = ('chr3', 41224645, 'T', 'C')
        data_tsv = {key: {constants.VcfInfo.PCGR_CSQ: 'tsv_csq'}}
        data_vcf = {key: {constants.VcfInfo.PCGR_MUTATION_HOTSPOT: 'vcf_hotspot'}}
        result = pcgr.compile_annotation_data(data_tsv, data_vcf)
        self.assertEqual(result[key][constants.VcfInfo.PCGR_CSQ], 'tsv_csq')
        self.assertEqual(result[key][constants.VcfInfo.PCGR_MUTATION_HOTSPOT], 'vcf_hotspot')


class TestAnnotateRecord(unittest.TestCase):
    """Unit tests for pcgr.annotate_record()."""

    def test_matching_key_writes_info(self):
        variant = _make_variant(
            chrom='chr3', pos=41224645, ref='T', alt='C',
            target_infos=[constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT],
        )
        key = ('chr3', 41224645, 'T', 'C')
        annotations = {key: {constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT: 7}}
        annotated = pcgr.annotate_record(variant, annotations, allow_missing=False)
        self.assertEqual(annotated.INFO.get(constants.VcfInfo.PCGR_TCGA_PANCANCER_COUNT.value), 7)

    def test_missing_key_allow_missing_returns_unchanged(self):
        variant = _make_variant(chrom='chr3', pos=41224645, ref='T', alt='C')
        annotations = {('chr3', 999, 'A', 'G'): {}}
        annotated = pcgr.annotate_record(variant, annotations, allow_missing=True)
        self.assertIs(annotated, variant)

    def test_missing_key_disallow_missing_raises(self):
        # NOTE: annotate_record's `assert key not in annotations` guard is tautological
        # (it re-checks a condition already established by the enclosing `if`), so it never
        # fires. Control falls through to `annotations[key].items()`, which raises KeyError
        # instead of the AssertionError one might expect from the guard's intent. Asserting
        # the actual (KeyError) behavior here rather than "fixing" the dead assert, since the
        # net effect — failing loudly on an unresolvable missing key — is preserved either way.
        variant = _make_variant(chrom='chr3', pos=41224645, ref='T', alt='C')
        annotations = {('chr3', 999, 'A', 'G'): {}}
        with self.assertRaises(KeyError):
            pcgr.annotate_record(variant, annotations, allow_missing=False)


if __name__ == '__main__':
    unittest.main()
