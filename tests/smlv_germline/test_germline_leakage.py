import os
import pathlib
import sys
import unittest


import cyvcf2


import bolt.common.constants as bolt_constants
import bolt.workflows.smlv_germline.prepare as smlv_germline_prepare



def silence_stdouterr():
    sys.stdout = sys.stderr = open(os.devnull, 'w')


def restore_stdouterr():
    if sys.stderr.name == os.devnull:
        sys.stdout.close()
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__


class TestGermlineLeakage(unittest.TestCase):

    def setUp(self):
        self.assets_dir = pathlib.Path(__file__).parent / 'assets/'


    def test_add_germline_leakage(self):

        # NOTE(SW): avoiding decorator pattern for simplicity

        silence_stdouterr()
        final_vcf = smlv_germline_prepare.add_germline_leakage(
            str(self.assets_dir / 'gl.germline_calls.vcf.gz'),
            str(self.assets_dir / 'gl.somatic_calls.vcf.gz'),
            'sample_tumor',
            'sample_normal',
        )
        restore_stdouterr()

        records = dict()
        for record in cyvcf2.VCF(final_vcf):
            key = (record.CHROM, record.POS)
            assert key not in records
            records[key] = record

        expected_tags = {
            ('chr1', 10283239): bolt_constants.VcfInfo.GERMLINE_LEAKAGE_CALLED,
            ('chr7', 5996592): bolt_constants.VcfInfo.GERMLINE_LEAKAGE,
            ('chr7', 38352621): bolt_constants.VcfInfo.GERMLINE_LEAKAGE,
        }

        for key, tag_enum in expected_tags.items():
            assert key in records
            assert records[key].INFO.get(tag_enum.value) is not None
