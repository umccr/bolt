"""Unit tests for bolt/util.py — VCF header helpers and merge helpers.

Covers binary-free logic plus a bcftools-guarded integration test for
merge_vcf_files (see TestMergeVcfFiles). The remaining bcftools-dependent
functions (count_vcf_records, execute_command) are intentionally NOT tested
here.
"""
import gzip
import pathlib
import shutil
import subprocess
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


class TestCheckAnnotationHeaders(unittest.TestCase):
    """Unit tests for util.check_annotation_headers()."""

    def _write_vcf_with_sage_hotspot(self, path, description):
        with open(path, 'w') as fh:
            fh.write(
                '##fileformat=VCFv4.2\n'
                '##FILTER=<ID=PASS,Description="All filters passed">\n'
                f'##INFO=<ID=SAGE_HOTSPOT,Number=0,Type=Flag,Description="{description}">\n'
                '##contig=<ID=chr1,length=248956422>\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                'chr1\t100\t.\tA\tT\t.\tPASS\t.\n'
            )

    def _write_vcf_without_sage_fields(self, path):
        with open(path, 'w') as fh:
            fh.write(
                '##fileformat=VCFv4.2\n'
                '##FILTER=<ID=PASS,Description="All filters passed">\n'
                '##contig=<ID=chr1,length=248956422>\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                'chr1\t100\t.\tA\tT\t.\tPASS\t.\n'
            )

    def test_matching_description_returns_normally(self):
        expected = constants.VCF_HEADER_ENTRIES[constants.VcfInfo.SAGE_HOTSPOT]['Description']
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'matching.vcf'
            self._write_vcf_with_sage_hotspot(vcf_fp, expected)
            # No exception/SystemExit raised
            util.check_annotation_headers(
                {constants.VcfInfo.SAGE_HOTSPOT: 'SAGE_HOTSPOT'}, vcf_fp,
            )

    def test_mismatched_description_exits(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'mismatched.vcf'
            self._write_vcf_with_sage_hotspot(vcf_fp, 'a totally different description')
            with self.assertRaises(SystemExit):
                util.check_annotation_headers(
                    {constants.VcfInfo.SAGE_HOTSPOT: 'SAGE_HOTSPOT'}, vcf_fp,
                )

    def test_field_absent_from_target_vcf_is_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'no_sage.vcf'
            self._write_vcf_without_sage_fields(vcf_fp)
            # SAGE_HOTSPOT has no header entry in this VCF at all; must be
            # skipped rather than raising, so no exception/SystemExit here.
            util.check_annotation_headers(
                {constants.VcfInfo.SAGE_HOTSPOT: 'SAGE_HOTSPOT'}, vcf_fp,
            )


@unittest.skipUnless(shutil.which('bcftools'), 'bcftools not available')
class TestMergeVcfFiles(unittest.TestCase):
    """Integration tests for util.merge_vcf_files().

    merge_vcf_files reassembles PCGR hypermutated chunk outputs with
    `bcftools merge -m all`. These chunks are sites-only VCFs (no FORMAT or
    sample columns) because pcgr.prepare_vcf_somatic / get_minimal_header strip
    them. bcftools merge only fails with "Duplicate sample names" when inputs
    carry a same-named genotype column; on sites-only inputs it correctly
    produces the union. These tests lock that invariant in: any regression that
    reintroduces a sample column (which would break the merge) is caught here.
    """

    # Sites-only header (no FORMAT, no sample column) — mirrors get_minimal_header
    _HEADER = (
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1,length=248956422>\n'
        '##contig=<ID=chr2,length=242193529>\n'
        '##INFO=<ID=PCGR_TIER,Number=1,Type=String,Description="tier">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )

    def _write_chunk(self, path, records):
        """Write a sites-only VCF, then bgzip + tabix-index it via bcftools.

        `records` is an iterable of (chrom, pos, ref, alt) tuples. Positions
        within a chunk are sorted before writing so indexing succeeds; chunks
        may be mutually out of order to exercise the cross-chunk sort.
        """
        plain = pathlib.Path(f'{path}.plain.vcf')
        with open(plain, 'w') as fh:
            fh.write(self._HEADER)
            for chrom, pos, ref, alt in sorted(records, key=lambda r: (r[0], r[1])):
                fh.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tPCGR_TIER=1\n')
        subprocess.run(['bcftools', 'view', '-Oz', '-o', str(path), str(plain)], check=True)
        subprocess.run(['bcftools', 'index', '-t', str(path)], check=True)
        return records

    def _read_keys(self, vcf_fp):
        return [
            (record.CHROM, record.POS, record.REF, record.ALT[0])
            for record in cyvcf2.VCF(str(vcf_fp))
        ]

    def test_merge_is_lossless_and_sorted(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            # Disjoint positions, and chunks deliberately out of order relative
            # to each other so the merge must interleave/sort across chunks.
            chunk_a = self._write_chunk(
                tmp_path / 'chunk_a.vcf.gz',
                [('chr1', 300, 'G', 'A'), ('chr1', 100, 'A', 'T'), ('chr2', 50, 'C', 'G')],
            )
            chunk_b = self._write_chunk(
                tmp_path / 'chunk_b.vcf.gz',
                [('chr1', 200, 'C', 'G'), ('chr1', 400, 'T', 'C')],
            )

            merged_vcf = util.merge_vcf_files(
                [tmp_path / 'chunk_a.vcf.gz', tmp_path / 'chunk_b.vcf.gz'],
                tmp_path / 'merged.pass',
            )

            merged_keys = self._read_keys(merged_vcf)
            expected_keys = list(chunk_a) + list(chunk_b)

            # No loss, no duplication: exact multiset match
            self.assertEqual(len(merged_keys), len(expected_keys))
            self.assertCountEqual(merged_keys, expected_keys)
            # Position-sorted output (cross-chunk interleave)
            self.assertEqual(
                merged_keys,
                sorted(merged_keys, key=lambda k: (k[0], k[1])),
            )

    def test_merge_output_is_bgzipped_and_indexed(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            self._write_chunk(tmp_path / 'chunk_a.vcf.gz', [('chr1', 100, 'A', 'T')])
            self._write_chunk(tmp_path / 'chunk_b.vcf.gz', [('chr1', 200, 'C', 'G')])

            merged_vcf = util.merge_vcf_files(
                [tmp_path / 'chunk_a.vcf.gz', tmp_path / 'chunk_b.vcf.gz'],
                tmp_path / 'merged.pass',
            )

            # BGZF/gzip magic bytes
            with open(merged_vcf, 'rb') as fh:
                self.assertEqual(fh.read(2), b'\x1f\x8b')
            # merge_vcf_files tabix-indexes its output
            self.assertTrue(pathlib.Path(f'{merged_vcf}.tbi').exists())
            # Intermediate unsorted file is cleaned up
            self.assertFalse((tmp_path / 'merged.pass.unsorted.vcf.gz').exists())


if __name__ == '__main__':
    unittest.main()
