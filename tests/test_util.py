"""Unit tests for bolt/util.py — VCF header helpers and merge_tsv_files.

Covers only binary-free logic. Functions requiring bcftools (count_vcf_records,
merge_vcf_files, execute_command) are intentionally NOT tested here.
"""
import gzip
import pathlib
import tempfile
import unittest

import cyvcf2

import bolt.common.constants as constants
import bolt.util as util


# Minimal VCF header covering enum members exercised below
HEADER = (
    '##fileformat=VCFv4.2\n'
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    f'##FILTER=<ID={constants.VcfFilter.MIN_AF.value},Description="">\n'
    f'##INFO=<ID={constants.VcfInfo.SAGE_NOVEL.value},Number=0,Type=Flag,Description="">\n'
    '##contig=<ID=chr1,length=248956422>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
)


def _write_vcf(path, variants):
    with open(path, 'w') as fh:
        fh.write(HEADER)
        for pos, info in variants:
            fh.write(f'chr1\t{pos}\t.\tA\tT\t.\tPASS\t{info}\n')


def _make_vcf_handle(tmp_dir):
    """Return a cyvcf2.VCF handle opened on a minimal VCF written to tmp_dir."""
    vcf_path = pathlib.Path(tmp_dir) / 'test.vcf'
    _write_vcf(vcf_path, [(100, '.')])
    return cyvcf2.VCF(str(vcf_path))


class TestGetVcfHeaderEntry(unittest.TestCase):
    """Unit tests for util.get_vcf_header_entry()."""

    def test_id_matches_enum_value(self):
        entry = util.get_vcf_header_entry(constants.VcfInfo.SAGE_NOVEL)
        self.assertEqual(entry['ID'], constants.VcfInfo.SAGE_NOVEL.value)

    def test_merges_header_entries_fields(self):
        entry = util.get_vcf_header_entry(constants.VcfInfo.SAGE_NOVEL)
        expected = constants.VCF_HEADER_ENTRIES[constants.VcfInfo.SAGE_NOVEL]
        for key, value in expected.items():
            self.assertEqual(entry[key], value)

    def test_filter_enum_entry(self):
        entry = util.get_vcf_header_entry(constants.VcfFilter.MIN_AF)
        self.assertEqual(entry['ID'], constants.VcfFilter.MIN_AF.value)
        self.assertIn('Description', entry)


class TestGetVcfHeaderLine(unittest.TestCase):
    """Unit tests for util.get_vcf_header_line()."""

    def test_filter_line_format(self):
        line = util.get_vcf_header_line(constants.VcfFilter.MIN_AF)
        self.assertTrue(line.startswith('##FILTER=<'))
        self.assertIn(f'ID={constants.VcfFilter.MIN_AF.value}', line)

    def test_info_line_format(self):
        line = util.get_vcf_header_line(constants.VcfInfo.SAGE_NOVEL)
        self.assertTrue(line.startswith('##INFO=<'))
        self.assertIn(f'ID={constants.VcfInfo.SAGE_NOVEL.value}', line)
        self.assertIn('Number=', line)
        self.assertIn('Type=', line)

    def test_format_line_format(self):
        line = util.get_vcf_header_line(constants.VcfFormat.SAGE_AD)
        self.assertTrue(line.startswith('##FORMAT=<'))
        self.assertIn(f'ID={constants.VcfFormat.SAGE_AD.value}', line)


class TestGetQualifiedVcfAnnotation(unittest.TestCase):
    """Unit tests for util.get_qualified_vcf_annotation()."""

    def test_info_enum_qualified(self):
        result = util.get_qualified_vcf_annotation(constants.VcfInfo.SAGE_NOVEL)
        self.assertEqual(result, f'INFO/{constants.VcfInfo.SAGE_NOVEL.value}')

    def test_format_enum_qualified(self):
        result = util.get_qualified_vcf_annotation(constants.VcfFormat.SAGE_AD)
        self.assertEqual(result, f'FORMAT/{constants.VcfFormat.SAGE_AD.value}')

    def test_filter_enum_raises(self):
        with self.assertRaises(AssertionError):
            util.get_qualified_vcf_annotation(constants.VcfFilter.MIN_AF)


class TestAddVcfHeaderEntry(unittest.TestCase):
    """Unit tests for util.add_vcf_header_entry()."""

    def test_info_tag_added_to_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            fh = _make_vcf_handle(tmp)
            util.add_vcf_header_entry(fh, constants.VcfInfo.PCGR_ACTIONABILITY_TIER)
            header_type = fh.get_header_type(constants.VcfInfo.PCGR_ACTIONABILITY_TIER.value)
            self.assertEqual(header_type['ID'], constants.VcfInfo.PCGR_ACTIONABILITY_TIER.value)

    def test_filter_tag_added_to_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            fh = _make_vcf_handle(tmp)
            util.add_vcf_header_entry(fh, constants.VcfFilter.PON)
            # FILTER header lines live under BCF_HL_FLT (order=0); get_header_type
            # defaults to [INFO, FORMAT] (order=[1, 2]) so FILTER lookups need order=[0].
            header_type = fh.get_header_type(constants.VcfFilter.PON.value, order=[0])
            self.assertEqual(header_type['ID'], constants.VcfFilter.PON.value)

    def test_format_tag_added_to_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            fh = _make_vcf_handle(tmp)
            util.add_vcf_header_entry(fh, constants.VcfFormat.SAGE_DP)
            header_type = fh.get_header_type(constants.VcfFormat.SAGE_DP.value)
            self.assertEqual(header_type['ID'], constants.VcfFormat.SAGE_DP.value)


class TestMergeTsvFiles(unittest.TestCase):
    """Unit tests for util.merge_tsv_files()."""

    def _write_gz_tsv(self, path, rows):
        with gzip.open(path, 'wt', encoding='utf-8') as fh:
            for row in rows:
                fh.write(row + '\n')

    def test_single_header_row_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            tsv1 = tmp_path / 'a.tsv.gz'
            tsv2 = tmp_path / 'b.tsv.gz'
            self._write_gz_tsv(tsv1, ['CHROM\tPOS', 'chr1\t100'])
            self._write_gz_tsv(tsv2, ['CHROM\tPOS', 'chr1\t200'])

            merged_fp = tmp_path / 'merged.tsv.gz'
            util.merge_tsv_files([tsv1, tsv2], merged_fp)

            with gzip.open(merged_fp, 'rt', encoding='utf-8') as fh:
                lines = [line.rstrip('\n') for line in fh]

            self.assertEqual(lines.count('CHROM\tPOS'), 1)

    def test_all_data_rows_present_in_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            tsv1 = tmp_path / 'a.tsv.gz'
            tsv2 = tmp_path / 'b.tsv.gz'
            self._write_gz_tsv(tsv1, ['CHROM\tPOS', 'chr1\t100', 'chr1\t101'])
            self._write_gz_tsv(tsv2, ['CHROM\tPOS', 'chr1\t200'])

            merged_fp = tmp_path / 'merged.tsv.gz'
            util.merge_tsv_files([tsv1, tsv2], merged_fp)

            with gzip.open(merged_fp, 'rt', encoding='utf-8') as fh:
                lines = [line.rstrip('\n') for line in fh]

            self.assertEqual(
                lines,
                ['CHROM\tPOS', 'chr1\t100', 'chr1\t101', 'chr1\t200'],
            )

    def test_output_is_gzipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            tsv1 = tmp_path / 'a.tsv.gz'
            self._write_gz_tsv(tsv1, ['CHROM\tPOS', 'chr1\t100'])

            merged_fp = tmp_path / 'merged.tsv.gz'
            util.merge_tsv_files([tsv1], merged_fp)

            # gzip files start with the magic number 0x1f 0x8b
            with open(merged_fp, 'rb') as fh:
                magic = fh.read(2)
            self.assertEqual(magic, b'\x1f\x8b')


if __name__ == '__main__':
    unittest.main()
