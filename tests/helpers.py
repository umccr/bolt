"""Shared VCF-construction fixtures used across test modules."""
import pathlib
import tempfile

import cyvcf2


# Minimal CSQ: only tokens[1] (consequence) is read by get_impacts()
def _csq(consequence):
    return f'A|{consequence}|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.'


# Minimal VCF header with all INFO fields used by select_pcgr_variants
HEADER = (
    '##fileformat=VCFv4.2\n'
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##INFO=<ID=PCGR_ACTIONABILITY_TIER,Number=1,Type=String,Description="">\n'
    '##INFO=<ID=PCGR_CSQ,Number=.,Type=String,Description="">\n'
    '##INFO=<ID=HMF_HOTSPOT,Number=0,Type=Flag,Description="">\n'
    '##INFO=<ID=PCGR_MUTATION_HOTSPOT,Number=.,Type=String,Description="">\n'
    '##INFO=<ID=SAGE_HOTSPOT,Number=0,Type=Flag,Description="">\n'
    '##INFO=<ID=PANEL,Number=0,Type=Flag,Description="">\n'
    '##INFO=<ID=GIAB_CONF,Number=0,Type=Flag,Description="">\n'
    '##INFO=<ID=DIFFICULT_segdup,Number=0,Type=Flag,Description="">\n'
    '##contig=<ID=chr1,length=248956422>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
)


def _write_vcf(path, variants):
    with open(path, 'w') as fh:
        fh.write(HEADER)
        for pos, info in variants:
            fh.write(f'chr1\t{pos}\t.\tA\tT\t.\tPASS\t{info}\n')


def _count_vcf(fp):
    return sum(1 for _ in cyvcf2.VCF(str(fp)))


def _make_variant(info_str):
    """Return a cyvcf2 Variant built from info_str using the test VCF header."""
    with tempfile.TemporaryDirectory() as tmp:
        vcf_path = pathlib.Path(tmp) / 'test.vcf'
        _write_vcf(vcf_path, [(100, info_str)])
        return list(cyvcf2.VCF(str(vcf_path)))[0]
