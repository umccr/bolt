"""Tests for count_variant_process counts."""

import pathlib
import tempfile
import unittest

import bolt.common.constants as constants
import bolt.workflows.smlv_somatic.report as report_mod


class TestCountVariantProcess(unittest.TestCase):
    """Verify count_variant_process counts.

    raw_pass counts DRAGEN raw PASS + SAGE novel, the caller-passed variants before bolt variant
    filters. Non-passing DRAGEN variants are excluded. The is_hypermutated and filtering_type
    fields are computed in the annotate step and merged in by entry; count_variant_process does
    not produce them.
    """

    COUNT_HEADER = (
        '##fileformat=VCFv4.2\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FILTER=<ID=DRAGEN_LOW_QUAL,Description="non-passing DRAGEN filter">\n'
        f'##FILTER=<ID={constants.VcfFilter.MIN_AF.value},Description="">\n'
        f'##FILTER=<ID={constants.VcfFilter.PON.value},Description="">\n'
        f'##FILTER=<ID={constants.VcfFilter.MAX_VARIANTS_NON_PASS.value},Description="">\n'
        f'##INFO=<ID={constants.VcfInfo.SAGE_NOVEL.value},Number=0,Type=Flag,Description="">\n'
        f'##INFO=<ID={constants.VcfInfo.RESCUED_FILTERS_EXISTING.value},Number=1,Type=String,Description="">\n'
        '##contig=<ID=chr1,length=248956422>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )

    def _write_count_vcf(self, path, rows):
        """Write a VCF with the given rows: list of (pos, filter, info) tuples."""
        with open(path, 'w') as fh:
            fh.write(self.COUNT_HEADER)
            for pos, filt, info in rows:
                fh.write(f'chr1\t{pos}\t.\tA\tT\t.\t{filt}\t{info}\n')

    def test_raw_pass_includes_variant_filtered_variants(self):
        """raw_pass counts PASS and bolt-variant-filtered variants; filter_pass counts PASS only."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            min_af = constants.VcfFilter.MIN_AF.value
            rows = [(i * 10, 'PASS', '.') for i in range(1, 4)]
            rows += [(i * 10 + 5, min_af, '.') for i in range(1, 3)]
            self._write_count_vcf(vcf_fp, rows)

            counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['raw_pass'], 5)
            self.assertEqual(counts['filter_pass'], 3)

    def test_raw_pass_excludes_nonpass_dragen(self):
        """Non-passing DRAGEN variants are excluded from raw_pass."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            rows = [(i * 10, 'PASS', '.') for i in range(1, 4)]
            rows += [(i * 10 + 5, 'DRAGEN_LOW_QUAL', '.') for i in range(1, 3)]
            self._write_count_vcf(vcf_fp, rows)

            counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['dragen'], 5)
            self.assertEqual(counts['raw_pass'], 3)

    def test_sage_novel_excluded_from_dragen_count(self):
        """SAGE_NOVEL variants are not DRAGEN variants, but are raw_pass."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            sage_novel = constants.VcfInfo.SAGE_NOVEL.value
            rows = [
                (10, 'PASS', '.'),
                (20, 'PASS', sage_novel),
                (30, 'PASS', '.'),
            ]
            self._write_count_vcf(vcf_fp, rows)

            counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['dragen'], 2)
            self.assertEqual(counts['sage'], 3)
            self.assertEqual(counts['raw_pass'], 3)

    def test_annotation_filter_excluded_from_annotated_count(self):
        """Variants with bolt annotation filters are excluded from the annotated count."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            annot_filter = constants.VcfFilter.MAX_VARIANTS_NON_PASS.value
            rows = [
                (10, 'PASS', '.'),
                (20, annot_filter, '.'),
            ]
            self._write_count_vcf(vcf_fp, rows)

            counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['annotated'], 1)
            self.assertEqual(counts['dragen'], 2)

    def test_bolt_variant_filtered_variant_is_raw_pass_but_not_annotated(self):
        """A bolt-variant-filtered variant (min_AF + max_variants_non_pass) is raw_pass, not annotated."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            min_af = constants.VcfFilter.MIN_AF.value
            annot_filter = constants.VcfFilter.MAX_VARIANTS_NON_PASS.value
            self._write_count_vcf(vcf_fp, [(10, f'{min_af};{annot_filter}', '.')])

            counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['raw_pass'], 1)
            self.assertEqual(counts['annotated'], 0)

    def test_rescue_variant_counted_in_raw_pass_and_filter_pass(self):
        """A rescued DRAGEN variant (RESCUED_FILTERS_EXISTING) is raw_pass and filter_pass."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            min_af = constants.VcfFilter.MIN_AF.value
            rescued_info = f'{constants.VcfInfo.RESCUED_FILTERS_EXISTING.value}={min_af}'
            rows = [
                (10, 'PASS', '.'),
                (20, 'PASS', rescued_info),
            ]
            self._write_count_vcf(vcf_fp, rows)

            counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['dragen'], 2)
            self.assertEqual(counts['raw_pass'], 2)
            self.assertEqual(counts['filter_pass'], 2)

    def test_raw_pass_keyed_before_filter_pass(self):
        """raw_pass is recorded before filter_pass so it appears first in the report."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            self._write_count_vcf(vcf_fp, [(10, 'PASS', '.')])

            counts = report_mod.count_variant_process(vcf_fp)

            keys = list(counts.keys())
            self.assertLess(keys.index('raw_pass'), keys.index('filter_pass'))


if __name__ == '__main__':
    unittest.main()
