"""Tests for bolt/common/pcgr.py — tier ordering, filter categorisation, chunking."""
import pathlib
import tempfile
import unittest
from unittest.mock import patch

import bolt.common.constants as constants
import bolt.common.pcgr as pcgr

from tests.helpers import _csq, _count_vcf, _make_variant, _write_vcf


class TestTierOrdering(unittest.TestCase):
    """Verify the PCGR_TIERS_FILTERING fix: values and priority order."""

    def test_noncoding_filtered_before_tier1(self):
        """N (NONCODING) entries must all precede '1' (TIER_1) entries in get_ordering()."""
        ordering = pcgr.get_ordering()
        tiers = [key[0] for key in ordering]
        n_idx = [i for i, t in enumerate(tiers) if t == 'N']
        t1_idx = [i for i, t in enumerate(tiers) if t == '1']
        self.assertTrue(n_idx, 'No NONCODING (N) entries in get_ordering()')
        self.assertTrue(t1_idx, 'No TIER_1 (1) entries in get_ordering()')
        self.assertLess(max(n_idx), min(t1_idx),
                        'All NONCODING entries must precede all TIER_1 entries')

    def test_no_long_form_tier_values(self):
        """PCGR_TIERS_FILTERING must use short forms ('1'-'4', 'N'), not 'TIER_1' etc."""
        for v in constants.PCGR_TIERS_FILTERING:
            self.assertNotIn('TIER_', v,
                             f"Found long-form tier value '{v}' — must be short form")

    def test_priority_order(self):
        """Full ordering: N before 4 before 3 before 2 before 1."""
        expected = ('N', '4', '3', '2', '1')
        self.assertEqual(constants.PCGR_TIERS_FILTERING, expected)


class TestGetImpacts(unittest.TestCase):
    """Unit tests for pcgr.get_impacts() — CSQ string parsing."""

    def test_single_consequence(self):
        csq = _csq('intron_variant')
        self.assertEqual(pcgr.get_impacts(csq), {'intron_variant'})

    def test_multi_consequences_ampersand(self):
        """A single CSQ entry with two consequences joined by & returns both."""
        csq = _csq('intron_variant&upstream_gene_variant')
        self.assertEqual(pcgr.get_impacts(csq), {'intron_variant', 'upstream_gene_variant'})

    def test_multiple_csq_entries_union(self):
        """Comma-separated CSQ entries — returns the union of all consequences."""
        csq = f'{_csq("intron_variant")},{_csq("intergenic_variant")}'
        self.assertEqual(pcgr.get_impacts(csq), {'intron_variant', 'intergenic_variant'})


class TestDetermineFilter(unittest.TestCase):
    """Unit tests for pcgr.determine_filter() — filter category determination."""

    def _data(self, **overrides):
        base = {
            'tier': None,
            'difficult': False,
            'giab_conf': False,
            'intergenic': None,
            'intronic': None,
            'downstream': None,
            'upstream': None,
            'impacts_other': None,
        }
        base.update(overrides)
        return base

    def test_intergenic_difficult(self):
        data = self._data(intergenic=True, difficult=True)
        self.assertEqual(pcgr.determine_filter(data), ('intergenic', 'difficult'))

    def test_intergenic_no_region(self):
        data = self._data(intergenic=True, difficult=False, giab_conf=False)
        self.assertEqual(pcgr.determine_filter(data), ('intergenic', 'none'))

    def test_intergenic_giab_conf(self):
        data = self._data(intergenic=True, giab_conf=True)
        self.assertEqual(pcgr.determine_filter(data), ('intergenic', 'giab_conf'))

    def test_intronic_supersedes_intergenic(self):
        """When both intergenic and intronic are present, intronic wins (higher priority)."""
        data = self._data(intergenic=True, intronic=True, difficult=True)
        self.assertEqual(pcgr.determine_filter(data), ('intronic', 'difficult'))

    def test_impacts_other_highest_priority(self):
        """impacts_other is the last to be filtered — it wins over all other impacts."""
        data = self._data(
            intergenic=True, intronic=True, downstream=True,
            upstream=True, impacts_other=True, difficult=True,
        )
        self.assertEqual(pcgr.determine_filter(data), ('impacts_other', 'difficult'))

    def test_no_impact_returns_false(self):
        """A variant with no recognisable impact cannot be categorised."""
        data = self._data()  # all impacts None
        self.assertFalse(pcgr.determine_filter(data))

    def test_giab_conf_region(self):
        data = self._data(impacts_other=True, giab_conf=True)
        self.assertEqual(pcgr.determine_filter(data), ('impacts_other', 'giab_conf'))


class TestGetVariantFilterData(unittest.TestCase):
    """Unit tests for pcgr.get_variant_filter_data() — data extraction from VCF records."""

    def test_tier_extracted(self):
        info = f'PCGR_ACTIONABILITY_TIER=2;PCGR_CSQ={_csq("intron_variant")}'
        data = pcgr.get_variant_filter_data(_make_variant(info))
        self.assertEqual(data['tier'], '2')

    def test_intergenic_impact(self):
        info = f'PCGR_ACTIONABILITY_TIER=N;PCGR_CSQ={_csq("intergenic_variant")}'
        data = pcgr.get_variant_filter_data(_make_variant(info))
        self.assertTrue(data['intergenic'])
        self.assertFalse(data['intronic'])
        self.assertFalse(data['downstream'])
        self.assertFalse(data['upstream'])
        self.assertFalse(data['impacts_other'])

    def test_intronic_impact(self):
        info = f'PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'
        data = pcgr.get_variant_filter_data(_make_variant(info))
        self.assertTrue(data['intronic'])
        self.assertFalse(data['intergenic'])

    def test_giab_conf_overrides_difficult(self):
        """GIAB_CONF flag must clear the difficult flag even when DIFFICULT_* is also present."""
        info = f'GIAB_CONF;DIFFICULT_segdup;PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'
        data = pcgr.get_variant_filter_data(_make_variant(info))
        self.assertTrue(data['giab_conf'])
        self.assertFalse(data['difficult'])

    def test_difficult_without_giab(self):
        info = f'DIFFICULT_segdup;PCGR_ACTIONABILITY_TIER=1;PCGR_CSQ={_csq("intron_variant")}'
        data = pcgr.get_variant_filter_data(_make_variant(info))
        self.assertTrue(data['difficult'])
        self.assertFalse(data['giab_conf'])


class TestSplitVcf(unittest.TestCase):
    """Tests for pcgr.split_vcf() — chunking the annotation path for large VCFs.

    split_vcf() is the annotate-path strategy for hypermutated samples: it divides
    a VCF into ≤MAX_SOMATIC_VARIANTS chunks so each chunk can be run through PCGR
    independently. Tested 2026-05-13 with a synthetic 550k VCF: 550k → 450k + 100k.
    """

    def test_chunks_above_limit(self):
        """VCF exceeding the limit is split into correctly-sized chunks."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            vcf_fp = tmp_path / 'input.vcf'
            # 25 variants, limit=10 → expect 3 chunks (10, 10, 5)
            v = [(i * 10, f'PCGR_CSQ={_csq("intron_variant")}') for i in range(1, 26)]
            _write_vcf(vcf_fp, v)

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 10):
                chunks = pcgr.split_vcf(vcf_fp, tmp_path)

            self.assertEqual(len(chunks), 3)
            counts = [_count_vcf(c) for c in chunks]
            self.assertLessEqual(max(counts), 10)
            self.assertEqual(sum(counts), 25)

    def test_no_chunking_within_limit(self):
        """VCF within the limit produces a single chunk containing all variants."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            vcf_fp = tmp_path / 'input.vcf'
            v = [(i * 10, f'PCGR_CSQ={_csq("intron_variant")}') for i in range(1, 6)]
            _write_vcf(vcf_fp, v)

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 10):
                chunks = pcgr.split_vcf(vcf_fp, tmp_path)

            self.assertEqual(len(chunks), 1)
            self.assertEqual(_count_vcf(chunks[0]), 5)

    def test_chunks_are_gzipped(self):
        """Chunk files must be .vcf.gz — plain .vcf chunks violate CLAUDE.md and waste disk."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            vcf_fp = tmp_path / 'input.vcf'
            v = [(i * 10, f'PCGR_CSQ={_csq("intron_variant")}') for i in range(1, 26)]
            _write_vcf(vcf_fp, v)

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 10):
                chunks = pcgr.split_vcf(vcf_fp, tmp_path)

            for chunk in chunks:
                self.assertTrue(str(chunk).endswith('.vcf.gz'),
                                f'Expected .vcf.gz chunk, got: {chunk.name}')

    def test_chunks_are_tabix_indexed(self):
        """Each .vcf.gz chunk must have a .tbi index — PCGR v2.2.5 requires it."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            vcf_fp = tmp_path / 'input.vcf'
            v = [(i * 10, f'PCGR_CSQ={_csq("intron_variant")}') for i in range(1, 26)]
            _write_vcf(vcf_fp, v)

            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 10):
                chunks = pcgr.split_vcf(vcf_fp, tmp_path)

            for chunk in chunks:
                tbi = pathlib.Path(str(chunk) + '.tbi')
                self.assertTrue(tbi.exists(), f'Missing tabix index for {chunk.name}')


class TestRunSomaticChunkArgMapping(unittest.TestCase):
    """Regression test: run_somatic_chunk must forward args as keywords to run_somatic.

    Before the fix, run_somatic_chunk called run_somatic positionally (6 args),
    skipping pcgr_threads. This caused pcgr_conda ('pcgr') to land in the
    pcgr_threads slot → ValueError: invalid literal for int() with base 10: 'pcgr'.
    """

    def test_pcgr_conda_not_shifted_into_pcgr_threads(self):
        """pcgr_conda must reach run_somatic as pcgr_conda, not as pcgr_threads."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            vcf_fp = tmp_path / 'chunk.vcf'
            _write_vcf(vcf_fp, [(10, f'PCGR_CSQ={_csq("intron_variant")}')])

            captured = {}

            def fake_run_somatic(*args, **kwargs):
                captured['args'] = args
                captured['kwargs'] = kwargs
                return (None, None)

            with patch('bolt.common.pcgr.run_somatic', side_effect=fake_run_somatic), \
                 patch('bolt.common.pcgr.merging_pcgr_files',
                       return_value=(tmp_path / 'out.vcf', tmp_path / 'out.tsv')):
                pcgr.run_somatic_chunk(
                    [vcf_fp],
                    pcgr_data_dir=tmp_path / 'pcgr_data',
                    vep_dir=tmp_path / 'vep',
                    output_dir=tmp_path,
                    pcgr_output_dir=tmp_path / 'pcgr_output',
                    max_threads=4,
                    pcgr_conda='pcgr_env',
                    pcgrr_conda='pcgrr_env',
                )

            kw = captured['kwargs']
            self.assertEqual(kw.get('pcgr_conda'), 'pcgr_env',
                             'pcgr_conda was not forwarded — likely shifted into pcgr_threads')
            self.assertEqual(kw.get('pcgrr_conda'), 'pcgrr_env',
                             'pcgrr_conda was not forwarded correctly')
            self.assertEqual(kw.get('threads'), 4,
                             'threads (max_threads) was not forwarded correctly')
            self.assertEqual(kw.get('chunk_nbr'), 1,
                             'chunk_nbr was not forwarded correctly')

    def test_disable_estimates_passed_to_run_somatic(self):
        """run_somatic_chunk must pass disable_estimates=True to every run_somatic call.

        Chunked PCGR runs must not include --estimate_msi/--estimate_tmb per chunk —
        those flags produce per-chunk partial estimates that are meaningless after merging.
        Fixed in bolt 0.3.2 (umccr/sash#57 + disable_estimates wiring).
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            vcf_fp = tmp_path / 'chunk.vcf'
            _write_vcf(vcf_fp, [(10, f'PCGR_CSQ={_csq("intron_variant")}')])

            captured = {}

            def fake_run_somatic(*args, **kwargs):
                captured['kwargs'] = kwargs
                return (None, None)

            with patch('bolt.common.pcgr.run_somatic', side_effect=fake_run_somatic), \
                 patch('bolt.common.pcgr.merging_pcgr_files',
                       return_value=(tmp_path / 'out.tsv', tmp_path / 'out.vcf')):
                pcgr.run_somatic_chunk(
                    [vcf_fp],
                    pcgr_data_dir=tmp_path / 'pcgr_data',
                    vep_dir=tmp_path / 'vep',
                    output_dir=tmp_path,
                    pcgr_output_dir=tmp_path / 'pcgr_output',
                    max_threads=4,
                    pcgr_conda='pcgr_env',
                    pcgrr_conda='pcgrr_env',
                )

        self.assertTrue(
            captured['kwargs'].get('disable_estimates'),
            'run_somatic_chunk must forward disable_estimates=True — '
            'per-chunk MSI/TMB estimates are meaningless after merge',
        )


class TestRunSomaticCommandArgs(unittest.TestCase):
    """Verify run_somatic builds the correct PCGR command-line arguments."""

    def test_estimate_signatures_absent_from_command(self):
        """--estimate_signatures must not appear in the PCGR command (dropped in sash#57).

        --estimate_msi and --estimate_tmb must still be present for non-chunked runs.
        """
        captured = {}

        def fake_execute(cmd, **kwargs):
            captured['cmd'] = cmd
            # create the output files run_somatic expects to find after pcgr runs
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / 'nosampleset.pcgr.grch38.snv_indel_ann.tsv.gz').touch()
            (output_dir / 'nosampleset.pcgr.grch38.pass.vcf.gz').touch()

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / 'output'

            with patch('bolt.common.pcgr.util.execute_command', side_effect=fake_execute):
                pcgr.run_somatic(
                    input_fp=tmp_path / 'input.vcf.gz',
                    pcgr_refdata_dir=tmp_path / 'refdata',
                    vep_dir=tmp_path / 'vep',
                    output_dir=output_dir,
                )

        self.assertIn('cmd', captured, 'execute_command was not called')
        self.assertNotIn('--estimate_signatures', captured['cmd'],
                         '--estimate_signatures must be absent (dropped in sash#57)')
        self.assertIn('--estimate_msi', captured['cmd'],
                      '--estimate_msi must still be present for non-chunked runs')
        self.assertIn('--estimate_tmb', captured['cmd'],
                      '--estimate_tmb must still be present for non-chunked runs')

    def test_disable_estimates_suppresses_msi_tmb(self):
        """disable_estimates=True must suppress --estimate_msi and --estimate_tmb."""
        captured = {}

        def fake_execute(cmd, **kwargs):
            captured['cmd'] = cmd
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / 'nosampleset.pcgr.grch38.snv_indel_ann.tsv.gz').touch()
            (output_dir / 'nosampleset.pcgr.grch38.pass.vcf.gz').touch()

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            output_dir = tmp_path / 'output'

            with patch('bolt.common.pcgr.util.execute_command', side_effect=fake_execute):
                pcgr.run_somatic(
                    input_fp=tmp_path / 'input.vcf.gz',
                    pcgr_refdata_dir=tmp_path / 'refdata',
                    vep_dir=tmp_path / 'vep',
                    output_dir=output_dir,
                    disable_estimates=True,
                )

        self.assertNotIn('--estimate_msi', captured['cmd'],
                         '--estimate_msi must be absent when disable_estimates=True')
        self.assertNotIn('--estimate_tmb', captured['cmd'],
                         '--estimate_tmb must be absent when disable_estimates=True')


if __name__ == '__main__':
    unittest.main()
