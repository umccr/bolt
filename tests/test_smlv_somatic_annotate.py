"""Tests for select_variants filtering_type outcome."""

import pathlib
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import bolt.common.constants as constants
import bolt.workflows.smlv_somatic.annotate as annotate_mod


class TestSelectVariantsFilteringType(unittest.TestCase):
    """select_variants reports which max-variants ramp it applied.

    The four outcomes map to the return paths in select_variants:
    none, max_variants_non_pass, max_variants_gnomad, and max_variants_non_cancer_genes.
    If even the cancer-gene ramp leaves more than MAX_SOMATIC_VARIANTS variants, select_variants
    calls sys.exit(1) instead of returning.
    """

    HEADER = (
        '##fileformat=VCFv4.2\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FILTER=<ID=weak_evidence,Description="caller filter">\n'
        '##INFO=<ID=HMF_HOTSPOT,Number=0,Type=Flag,Description="">\n'
        '##INFO=<ID=SAGE_HOTSPOT,Number=0,Type=Flag,Description="">\n'
        '##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="">\n'
        '##contig=<ID=chr1,length=248956422>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )

    def _record(self, pos, filt='PASS', info='.'):
        return f'chr1\t{pos}\t.\tA\tT\t.\t{filt}\t{info}'

    def _write_vcf(self, fp, records):
        records = sorted(records, key=lambda r: int(r.split('\t')[1]))
        with open(fp, 'w') as fh:
            fh.write(self.HEADER)
            for rec in records:
                fh.write(rec + '\n')

    def _bgzip_index(self, vcf_fp, gz_fp):
        subprocess.run(['bcftools', 'view', '-Oz', '-o', str(gz_fp), str(vcf_fp)], check=True)
        subprocess.run(['bcftools', 'index', '-t', str(gz_fp)], check=True)

    def _run(self, records, max_variants, bed=None):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            output_dir = tmp / 'out'
            output_dir.mkdir()
            vcf_fp = tmp / 'input.vcf'
            gz_fp = tmp / 'input.vcf.gz'
            self._write_vcf(vcf_fp, records)
            self._bgzip_index(vcf_fp, gz_fp)

            bed_fp = None
            if bed is not None:
                bed_fp = tmp / 'genes.bed'
                bed_fp.write_text(bed)

            constants.MAX_SOMATIC_VARIANTS = max_variants
            return annotate_mod.select_variants(gz_fp, 'TUMOR', bed_fp, output_dir)

    def test_none_when_under_limit(self):
        records = [self._record(i * 10) for i in range(1, 5)]
        result = self._run(records, max_variants=10)
        self.assertEqual(result['filtering_type'], 'none')
        self.assertNotIn('filter_name', result)

    def test_non_pass_when_nonpass_ramp_resolves(self):
        records = [self._record(i * 10) for i in range(1, 5)]                      # 4 PASS
        records += [self._record(i * 10 + 5, 'weak_evidence') for i in range(1, 9)]  # 8 non-pass
        result = self._run(records, max_variants=10)
        self.assertEqual(result['filtering_type'], constants.VcfFilter.MAX_VARIANTS_NON_PASS.value)

    def test_gnomad_when_population_ramp_resolves(self):
        records = [self._record(i * 10, info='.') for i in range(1, 6)]             # 5 no-AF
        records += [self._record(60 + i * 10, info='gnomAD_AF=0.005') for i in range(2)]  # 2 low AF
        records += [self._record(90 + i * 10, info='gnomAD_AF=0.05') for i in range(3)]   # 3 common
        records += [self._record(130, 'weak_evidence', 'gnomAD_AF=0.05')]            # 1 common non-pass
        records += [self._record(140, info='HMF_HOTSPOT;gnomAD_AF=0.05')]            # 1 common hotspot
        result = self._run(records, max_variants=10)
        self.assertEqual(result['filtering_type'], constants.VcfFilter.MAX_VARIANTS_GNOMAD.value)

    def test_cancer_genes_when_gene_ramp_resolves(self):
        # chr1 10000-20000 gene region, padded to 9000-21000 by get_cancer_gene_variants
        bed = 'chr1\t10000\t20000\n'
        records = [self._record(15000 + i * 10) for i in range(3)]                  # 3 in-gene PASS
        records += [self._record(30000 + i * 10) for i in range(7)]                 # 7 outside PASS
        records += [self._record(40000 + i * 10, 'weak_evidence') for i in range(3)]  # 3 outside non-pass
        records += [self._record(50000, info='HMF_HOTSPOT')]                        # 1 outside hotspot
        result = self._run(records, max_variants=10, bed=bed)
        self.assertEqual(result['filtering_type'], constants.VcfFilter.MAX_VARIANTS_NON_CANCER_GENES.value)

    def test_unresolved_raises_when_gene_ramp_still_over_limit(self):
        records = [self._record(i * 10, info='HMF_HOTSPOT') for i in range(1, 13)]  # 12 hotspots
        with self.assertRaises(SystemExit):
            self._run(records, max_variants=10, bed='chr2\t1\t2\n')


class TestCountRawPass(unittest.TestCase):
    """count_raw_pass counts DRAGEN raw PASS + SAGE novel from the pre-filter VCF."""

    HEADER = (
        '##fileformat=VCFv4.2\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FILTER=<ID=weak_evidence,Description="caller filter">\n'
        '##INFO=<ID=SAGE_NOVEL,Number=0,Type=Flag,Description="">\n'
        '##contig=<ID=chr1,length=248956422>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )

    def _write_vcf(self, fp, records):
        with open(fp, 'w') as fh:
            fh.write(self.HEADER)
            for rec in records:
                fh.write(rec + '\n')

    def test_counts_pass_and_sage_novel_excludes_nonpass(self):
        records = [
            'chr1\t10\t.\tA\tT\t.\tPASS\t.',
            'chr1\t20\t.\tA\tT\t.\tweak_evidence\t.',
            'chr1\t30\t.\tA\tT\t.\tPASS\tSAGE_NOVEL',
        ]
        with tempfile.TemporaryDirectory() as tmp:
            fp = pathlib.Path(tmp) / 'input.vcf'
            self._write_vcf(fp, records)
            self.assertEqual(annotate_mod.count_raw_pass(fp), 2)

    def test_strict_greater_than_threshold(self):
        records = [f'chr1\t{i * 10}\t.\tA\tT\t.\tPASS\t.' for i in range(1, 5)]
        with tempfile.TemporaryDirectory() as tmp:
            fp = pathlib.Path(tmp) / 'input.vcf'
            self._write_vcf(fp, records)
            raw_pass = annotate_mod.count_raw_pass(fp)
            self.assertEqual(raw_pass, 4)
            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 4):
                self.assertFalse(raw_pass > constants.MAX_SOMATIC_VARIANTS)
            with patch('bolt.common.constants.MAX_SOMATIC_VARIANTS', 3):
                self.assertTrue(raw_pass > constants.MAX_SOMATIC_VARIANTS)


if __name__ == '__main__':
    unittest.main()
