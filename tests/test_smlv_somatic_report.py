"""Tests for bolt/workflows/smlv_somatic/report.py — hypermutated variant trimming."""
import pathlib
import shutil
import tempfile
import unittest
from unittest.mock import patch

from click.testing import CliRunner

import cyvcf2

import bolt.common.constants as constants
import bolt.common.pcgr as pcgr
import bolt.util as util
import bolt.workflows.smlv_somatic.report as report_mod

from tests.helpers import HEADER, _csq, _write_vcf


class TestSelectPcgrVariants(unittest.TestCase):
    """Integration tests for select_pcgr_variants() trimming logic."""

    def _run(self, variants, limit, tmp):
        """Run select_pcgr_variants with a small MAX_SOMATIC_VARIANTS limit."""
        vcf_fp = pathlib.Path(tmp) / 'input.vcf'
        _write_vcf(vcf_fp, variants)
        cancer_genes = pathlib.Path(tmp) / 'genes.bed'
        cancer_genes.write_text('chr1\t1\t9999999\n')

        # Mock bcftools annotate: copy input to the expected output path
        orig_execute = util.execute_command
        def fake_execute(cmd, **_):
            import re
            m = re.search(r'--output\s+(\S+)', cmd)
            if m and 'bcftools annotate' in cmd:
                shutil.copy(str(vcf_fp), m.group(1))
            else:
                orig_execute(cmd)

        with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', limit), \
             patch('bolt.util.execute_command', side_effect=fake_execute):
            out_fp = report_mod.select_pcgr_variants(
                vcf_fp, cancer_genes, 'TUMOR', pathlib.Path(tmp)
            )
        return sum(1 for _ in cyvcf2.VCF(str(out_fp)))

    def test_output_within_limit(self):
        """Output must never exceed MAX_SOMATIC_VARIANTS."""
        with tempfile.TemporaryDirectory() as tmp:
            # 15 variants: 2 hotspot + 5 TIER_1 + 5 TIER_3 + 3 NONCODING
            v = []
            for i in range(1, 3):    # hotspot
                v.append((i*10, f'HMF_HOTSPOT;PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'))
            for i in range(3, 8):    # TIER_1 intronic
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'))
            for i in range(8, 13):   # TIER_3 intronic
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=3;PCGR_CSQ={_csq("intron_variant")}'))
            for i in range(13, 16):  # NONCODING intergenic
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}'))
            count = self._run(v, limit=10, tmp=tmp)
            self.assertLessEqual(count, 10)

    def test_noncoding_dropped_before_tier1(self):
        """With limit = total - 3, the 3 NONCODING variants should be dropped (not TIER_1)."""
        with tempfile.TemporaryDirectory() as tmp:
            v = []
            for i in range(1, 6):   # 5 TIER_1 intronic
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'))
            for i in range(6, 9):   # 3 NONCODING intergenic
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}'))
            # limit=5: should drop the 3 NONCODING to get to 5
            count = self._run(v, limit=5, tmp=tmp)
            self.assertEqual(count, 5)

    def test_hotspots_never_dropped_by_tiered_filter(self):
        """Hotspot variants must survive tiered filtering.

        Note: HMF_HOTSPOT is not in RETAIN_FIELDS_FILTERING — these variants survive because
        they are TIER_1 (highest priority), not via the hotspot retention path.
        """
        with tempfile.TemporaryDirectory() as tmp:
            v = []
            for i in range(1, 3):   # 2 HMF_HOTSPOT TIER_1 variants (survive via tier priority)
                v.append((i*10, f'HMF_HOTSPOT;PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'))
            for i in range(3, 13):  # 10 NONCODING (should all be filtered)
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}'))
            # limit=2: only the 2 TIER_1 variants should remain
            count = self._run(v, limit=2, tmp=tmp)
            self.assertEqual(count, 2)

    def test_all_within_limit_nothing_filtered(self):
        """When total variants are below the limit, nothing is dropped."""
        with tempfile.TemporaryDirectory() as tmp:
            v = [(i*10, f'PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}')
                 for i in range(1, 6)]  # 5 NONCODING
            count = self._run(v, limit=10, tmp=tmp)
            self.assertEqual(count, 5)

    def test_retained_variants_bypass_tiered_filter(self):
        """Variants with PANEL or SAGE_HOTSPOT bypass tiered filtering and always survive."""
        with tempfile.TemporaryDirectory() as tmp:
            v = []
            for i in range(1, 3):   # 2 SAGE_HOTSPOT
                v.append((i*10, 'SAGE_HOTSPOT'))
            for i in range(3, 5):   # 2 PANEL-only
                v.append((i*10, 'PANEL'))
            for i in range(5, 15):  # 10 NONCODING that get dropped
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}'))
            count = self._run(v, limit=4, tmp=tmp)
            self.assertEqual(count, 4)

    def test_pcgr_mutation_hotspot_real_value_is_retained(self):
        """A real PCGR_MUTATION_HOTSPOT value (non-dot) must retain the variant."""
        with tempfile.TemporaryDirectory() as tmp:
            v = []
            for i in range(1, 3):   # 2 real hotspot variants — must survive
                v.append((i*10, f'PCGR_MUTATION_HOTSPOT=GRCH38_1_{i}_A_T;PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'))
            for i in range(3, 8):   # 5 NONCODING — dropped
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}'))
            count = self._run(v, limit=2, tmp=tmp)
            self.assertEqual(count, 2)

    def test_pcgr_mutation_hotspot_dot_not_treated_as_retained(self):
        """PCGR_MUTATION_HOTSPOT=. must not retain variants — '.' is a missing-value placeholder.

        cyvcf2 returns the string '.' (truthy) for String INFO fields written as '=.' by PCGR on
        every non-hotspot variant.  Without the fix, any(variant.INFO.get(e) ...) always returns
        True and ALL variants are treated as retained, so tiered filtering never drops anything and
        RuntimeError fires for any sample with >450k variants (sash #52 root cause).
        """
        with tempfile.TemporaryDirectory() as tmp:
            v = []
            for i in range(1, 3):   # 2 real SAGE_HOTSPOT — must be retained
                v.append((i*10, f'SAGE_HOTSPOT;PCGR_MUTATION_HOTSPOT=.;PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'))
            for i in range(3, 8):   # 5 NONCODING with PCGR_MUTATION_HOTSPOT=. — must be droppable
                v.append((i*10, f'PCGR_MUTATION_HOTSPOT=.;PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}'))
            # limit=2: only the 2 real SAGE_HOTSPOT variants survive; the 5 dot-placeholder ones are dropped
            count = self._run(v, limit=2, tmp=tmp)
            self.assertEqual(count, 2)

    def test_filters_set_vcf_marks_dropped_variants(self):
        """The traceability VCF marks filtered-out variants with PCGR_count_limit.

        The function drops entire categories, so we need two distinct categories:
        - 3 TIER_1 intronic (high priority — kept)
        - 2 NONCODING intergenic (lowest priority — dropped as a whole category)
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            vcf_fp = tmp_path / 'input.vcf'
            v = []
            for i in range(1, 4):   # 3 TIER_1 intronic — kept
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'))
            for i in range(4, 6):   # 2 NONCODING intergenic — dropped whole category
                v.append((i*10, f'PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}'))
            _write_vcf(vcf_fp, v)
            cancer_genes = tmp_path / 'genes.bed'
            cancer_genes.write_text('chr1\t1\t9999999\n')

            orig_execute = util.execute_command
            def fake_execute(cmd, **_):
                import re
                m = re.search(r'--output\s+(\S+)', cmd)
                if m and 'bcftools annotate' in cmd:
                    shutil.copy(str(vcf_fp), m.group(1))
                else:
                    orig_execute(cmd)

            # limit=3: the 2 NONCODING category is dropped, 3 TIER_1 survive
            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 3), \
                 patch('bolt.util.execute_command', side_effect=fake_execute):
                report_mod.select_pcgr_variants(vcf_fp, cancer_genes, 'TUMOR', tmp_path)

            filters_set_fp = tmp_path / 'TUMOR.pcgr_hypermutated.filters_set.vcf.gz'
            self.assertTrue(filters_set_fp.exists(), 'filters_set VCF not created')

            all_records = list(cyvcf2.VCF(str(filters_set_fp)))
            self.assertEqual(len(all_records), 5, 'filters_set VCF should contain all input variants')

            filter_tag = constants.VcfFilter.PCGR_COUNT_LIMIT.value
            dropped = [r for r in all_records if filter_tag in (r.FILTERS or [])]
            self.assertEqual(len(dropped), 2, 'Expected 2 NONCODING variants marked with PCGR_count_limit')


class TestSelectPcgrVariantsRaisesOnUnresolvableOverflow(unittest.TestCase):
    """select_pcgr_variants raises RuntimeError when all variants are retained (hotspots)."""

    def _fake_execute(self, vcf_fp):
        def _run(cmd, **_):
            import re
            m = re.search(r'--output\s+(\S+)', cmd)
            if m and 'bcftools annotate' in cmd:
                shutil.copy(str(vcf_fp), m.group(1))
            else:
                util.execute_command(cmd)
        return _run

    def test_raises_when_all_variants_are_hotspots(self):
        """All SAGE_HOTSPOT variants are RETAIN_FIELDS — tiered filtering cannot drop any; RuntimeError expected."""
        # RETAIN_FIELDS_FILTERING includes SAGE_HOTSPOT — use that flag, not HMF_HOTSPOT
        HOTSPOT_INFO = f'SAGE_HOTSPOT;PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            vcf_fp = tmp_path / 'input.vcf'
            variants = [(i * 10, HOTSPOT_INFO) for i in range(1, 6)]
            _write_vcf(vcf_fp, variants)
            cancer_genes = tmp_path / 'genes.bed'
            cancer_genes.write_text('chr1\t1\t9999999\n')

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 3), \
                 patch('bolt.util.execute_command', side_effect=self._fake_execute(vcf_fp)):
                with self.assertRaises(RuntimeError):
                    report_mod.select_pcgr_variants(vcf_fp, cancer_genes, 'TUMOR', tmp_path)


class TestCountVariantProcess(unittest.TestCase):
    """Verify count_variant_process counts and is_hypermutated flag (bolt #27).

    is_hypermutated must use the 'dragen' count (raw, pre-bolt-filter), not
    'filter_pass'. A sample with many DRAGEN variants that are mostly filtered
    away must still be flagged as hypermutated.
    """

    # Minimal header for count_variant_process: needs FILTER tags + SAGE_NOVEL INFO
    COUNT_HEADER = (
        '##fileformat=VCFv4.2\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        f'##FILTER=<ID={constants.VcfFilter.MIN_AF.value},Description="">\n'
        f'##FILTER=<ID={constants.VcfFilter.PON.value},Description="">\n'
        f'##FILTER=<ID={constants.VcfFilter.MAX_VARIANTS_NON_PASS.value},Description="">\n'
        f'##INFO=<ID={constants.VcfInfo.SAGE_NOVEL.value},Number=0,Type=Flag,Description="">\n'
        f'##INFO=<ID={constants.VcfInfo.RESCUED_FILTERS_EXISTING.value},Number=1,Type=String,Description="">\n'
        '##contig=<ID=chr1,length=248956422>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )

    def _write_count_vcf(self, path, rows):
        """rows: list of (pos, filter_str, info_str) tuples."""
        with open(path, 'w') as fh:
            fh.write(self.COUNT_HEADER)
            for pos, filt, info in rows:
                fh.write(f'chr1\t{pos}\t.\tA\tT\t.\t{filt}\t{info}\n')

    def test_is_hypermutated_uses_dragen_count(self):
        """is_hypermutated=True when dragen count > MAX_SOMATIC_VARIANTS even if filter_pass is below."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            min_af = constants.VcfFilter.MIN_AF.value
            # 3 DRAGEN PASS variants + 2 filtered by bolt (MIN_AF) — filter_pass=3, dragen=5
            rows = [(i * 10, 'PASS', '.') for i in range(1, 4)]
            rows += [(i * 10 + 5, min_af, '.') for i in range(1, 3)]
            self._write_count_vcf(vcf_fp, rows)

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 4):
                counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['dragen'], 5)
            self.assertEqual(counts['filter_pass'], 3)
            # dragen(5) > MAX(4) → hypermutated, even though filter_pass(3) ≤ MAX(4)
            self.assertTrue(counts['is_hypermutated'])

    def test_is_hypermutated_false_when_dragen_within_limit(self):
        """is_hypermutated=False when dragen count ≤ MAX_SOMATIC_VARIANTS."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            rows = [(i * 10, 'PASS', '.') for i in range(1, 4)]
            self._write_count_vcf(vcf_fp, rows)

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 10):
                counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['dragen'], 3)
            self.assertFalse(counts['is_hypermutated'])

    def test_sage_novel_excluded_from_dragen_count(self):
        """SAGE_NOVEL variants are not counted as DRAGEN variants."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            sage_novel_info = constants.VcfInfo.SAGE_NOVEL.value
            rows = [
                (10, 'PASS', '.'),              # dragen
                (20, 'PASS', sage_novel_info),  # sage novel — not dragen
                (30, 'PASS', '.'),              # dragen
            ]
            self._write_count_vcf(vcf_fp, rows)

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 100):
                counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['dragen'], 2)
            self.assertEqual(counts['sage'], 3)

    def test_annotation_filter_excluded_from_annotated_count(self):
        """Variants with bolt annotation filters are excluded from annotated count."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf_fp = pathlib.Path(tmp) / 'test.vcf'
            annot_filter = constants.VcfFilter.MAX_VARIANTS_NON_PASS.value
            rows = [
                (10, 'PASS', '.'),          # annotated
                (20, annot_filter, '.'),    # not annotated (bolt annotation filter)
            ]
            self._write_count_vcf(vcf_fp, rows)

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 100):
                counts = report_mod.count_variant_process(vcf_fp)

            self.assertEqual(counts['annotated'], 1)
            self.assertEqual(counts['dragen'], 2)


_PASS_COUNTS = {'pass': {'snps': 0, 'indels': 0, 'others': 0, 'total': 0}}


def _cli_args(dummy, output_dir):
    """Return CliRunner args list for report entry(); all file paths point to dummy."""
    d = str(dummy)
    return [
        '--tumor_name', 'TUMOR',
        '--normal_name', 'NORMAL',
        '--vcf_fp', d,
        '--vcf_filters_fp', d,
        '--vcf_dragen_fp', d,
        '--vep_dir', str(dummy.parent),
        '--purple_purity_fp', d,
        '--cancer_genes_fp', d,
        '--giab_regions_fp', d,
        '--genome_fp', d,
        '--threads', '1',
        '--output_dir', str(output_dir),
    ]


class TestEntrySkipsPcgrOnOverflow(unittest.TestCase):
    """entry() catches RuntimeError from select_pcgr_variants and skips PCGR entirely."""

    def test_run_somatic_not_called_when_cap_exceeded(self):
        """When select_pcgr_variants raises RuntimeError, run_somatic must not be called."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            dummy = tmp_path / 'dummy.vcf.gz'
            dummy.touch()

            with patch.object(report_mod, 'bcftools_stats_prepare', return_value=dummy), \
                 patch.object(report_mod, 'run_bcftools_stats'), \
                 patch.object(report_mod, 'allele_frequencies'), \
                 patch.object(report_mod, 'count_variant_types', return_value=_PASS_COUNTS), \
                 patch.object(report_mod, 'count_variant_process',
                              return_value={'filter_pass': constants.MAX_SOMATIC_VARIANTS + 1}), \
                 patch.object(report_mod, 'parse_purple_purity_file',
                              return_value={'purity': 0.8, 'ploidy': 2.0}), \
                 patch.object(report_mod, 'select_pcgr_variants',
                              side_effect=RuntimeError('595416 > 450000')), \
                 patch.object(pcgr, 'prepare_vcf_somatic') as mock_prep, \
                 patch.object(pcgr, 'run_somatic') as mock_run:
                result = CliRunner().invoke(report_mod.entry, _cli_args(dummy, tmp_path / 'out'))

            self.assertEqual(result.exit_code, 0, result.output)
            mock_prep.assert_not_called()
            mock_run.assert_not_called()

    def test_run_somatic_called_when_within_limit(self):
        """When PASS count is within limit, run_somatic must be called normally."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            dummy = tmp_path / 'dummy.vcf.gz'
            dummy.touch()
            fake_prep_output = tmp_path / 'prep.vcf.gz'
            fake_prep_output.touch()

            with patch.object(report_mod, 'bcftools_stats_prepare', return_value=dummy), \
                 patch.object(report_mod, 'run_bcftools_stats'), \
                 patch.object(report_mod, 'allele_frequencies'), \
                 patch.object(report_mod, 'count_variant_types', return_value=_PASS_COUNTS), \
                 patch.object(report_mod, 'count_variant_process',
                              return_value={'filter_pass': constants.MAX_SOMATIC_VARIANTS - 1}), \
                 patch.object(report_mod, 'parse_purple_purity_file',
                              return_value={'purity': 0.8, 'ploidy': 2.0}), \
                 patch.object(pcgr, 'prepare_vcf_somatic', return_value=fake_prep_output), \
                 patch.object(pcgr, 'run_somatic') as mock_run:
                result = CliRunner().invoke(report_mod.entry, _cli_args(dummy, tmp_path / 'out'))

            self.assertEqual(result.exit_code, 0, result.output)
            mock_run.assert_called_once()


if __name__ == '__main__':
    unittest.main()
