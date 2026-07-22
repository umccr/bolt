"""Unit tests for bolt/workflows/smlv_somatic/rescue.py.

Covers only the SAGE VCF header-consistency check in
annotate_existing_sage_calls(), which is reachable without a bcftools binary
(the check runs before any subprocess call — a SystemExit from a mismatch
short-circuits execution before bcftools is ever invoked). The rest of
rescue.py orchestrates bcftools subprocesses end-to-end and is not
unit-tested here.
"""
import pathlib
import tempfile
import unittest

import bolt.common.constants as constants
import bolt.workflows.smlv_somatic.rescue as rescue


def _write_input_vcf_no_sage_fields(path):
    """A plain DRAGEN-style VCF with no SAGE_* header lines at all."""
    with open(path, 'w') as fh:
        fh.write(
            '##fileformat=VCFv4.2\n'
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
            '##contig=<ID=chr1,length=248956422>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n'
            'chr1\t100\t.\tA\tT\t.\tPASS\t.\tAD\t10,5\n'
        )


def _write_sage_vcf_mismatched_hotspot_description(path):
    """A SAGE-style VCF whose SAGE_HOTSPOT description does not match constants.py."""
    with open(path, 'w') as fh:
        fh.write(
            '##fileformat=VCFv4.2\n'
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            '##INFO=<ID=SAGE_HOTSPOT,Number=0,Type=Flag,Description="wrong description">\n'
            '##contig=<ID=chr1,length=248956422>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            'chr1\t100\t.\tA\tT\t.\tPASS\tSAGE_HOTSPOT\n'
        )


class TestAnnotateExistingSageCallsHeaderCheck(unittest.TestCase):
    """Regression test for rescue.py's SAGE VCF header-consistency check.

    annotate_existing_sage_calls(input_fp, tumor_name, sage_vcf_fp, output_dir)
    must validate header descriptions against sage_vcf_fp (the SAGE VCF), not
    input_fp (the DRAGEN VCF being annotated). Prior to the fix, the check was
    wired to input_fp, which has no SAGE_* header lines at all — so a mismatch
    in the real SAGE VCF's headers went undetected.
    """

    def test_header_mismatch_in_sage_vcf_triggers_exit(self):
        """A mismatched SAGE_HOTSPOT description in sage_vcf_fp must raise SystemExit.

        input_fp has no SAGE_* headers at all, so if the check were (incorrectly)
        run against input_fp instead, every field would be silently skipped and
        no SystemExit would be raised — proving the check now reads sage_vcf_fp.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            input_fp = tmp_path / 'input.vcf'
            sage_vcf_fp = tmp_path / 'sage.vcf'
            _write_input_vcf_no_sage_fields(input_fp)
            _write_sage_vcf_mismatched_hotspot_description(sage_vcf_fp)

            with self.assertRaises(SystemExit):
                rescue.annotate_existing_sage_calls(
                    input_fp, 'TUMOR', sage_vcf_fp, tmp_path,
                )


if __name__ == '__main__':
    unittest.main()
