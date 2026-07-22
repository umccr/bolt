# Testing

## Running tests

Canonical (matches CI):

```bash
python -m unittest discover tests/ --buffer
```

Also supported:

```bash
python -m pytest tests/ -v
```

Single test:

```bash
python -m pytest tests/test_smlv_somatic_filter.py::TestSmlvSomaticFilter::test_min_af_filter -v
```

Tests in `tests/` are pure Python + in-memory `cyvcf2` logic and must run
without any bioinformatics binary installed (no `bcftools`, `pcgr`, `cpsr`,
`vcfanno`, `snpEff`, `gpgr`, VEP). Functions that shell out to those tools are
either mocked/patched in tests or are not unit-tested (see below). The one
exception is `TestMergeVcfFiles`, an integration test that exercises the real
`merge_vcf_files` → `bcftools merge` path; it is guarded with
`@unittest.skipUnless(shutil.which('bcftools'), ...)`, so it runs in the conda
CI env and skips cleanly (never fails) where `bcftools` is absent.

## Test Coverage

| Module | Status | Functions covered | Test file |
|---|---|---|---|
| `bolt/util.py` | Partial | `get_vcf_header_entry`, `get_vcf_header_line`, `get_qualified_vcf_annotation`, `add_vcf_header_entry`, `merge_tsv_files`, `merge_vcf_files` (bcftools-guarded lossless/sorted integration test), `check_annotation_headers` | `tests/test_util.py` |
| `bolt/common/pcgr.py` | Partial | `get_ordering`, `get_impacts`, `determine_filter`, `get_variant_filter_data`, `split_vcf`, `run_somatic_chunk` (arg-mapping regression) | `tests/test_pcgr.py` |
| `bolt/common/pcgr.py` | Partial | `parse_genomic_change`, `get_impacts_higher`, `get_annotation_entry_tsv`, `compile_annotation_data`, `annotate_record`, `get_annotations_vcf` (duplicate-key regression), `collect_pcgr_annotation_data` (duplicate-key tier resolution), `collect_cpsr_annotation_data` (duplicate-key regression) | `tests/test_pcgr_annotation.py` |
| `bolt/workflows/smlv_somatic/filter.py` | Partial | `set_filter_data` | `tests/test_smlv_somatic_filter.py` |
| `bolt/workflows/smlv_somatic/report.py` | Partial | `select_pcgr_variants`, `count_variant_process`, `entry` overflow handling, `entry` `disable_estimates` branching | `tests/test_smlv_somatic_report.py` |
| `bolt/workflows/smlv_somatic/rescue.py` | Partial | `annotate_existing_sage_calls` (SAGE VCF header-consistency check only) | `tests/test_smlv_somatic_rescue.py` |

"Partial" means the module has meaningful test coverage for its pure/testable logic,
but not every function in the file is tested (see below for what is excluded and why).

## Untested / not unit-testable in CI

These require live bioinformatics binaries (bcftools, PCGR/CPSR, vcfanno, snpEff,
gpgr) or full end-to-end subprocess orchestration, so they are excluded from the
unit test suite:

| Function/module | Reason |
|---|---|
| `bolt/util.py: count_vcf_records` | Shells out to `bcftools view` |
| `bolt/util.py: execute_command` | Spawns real subprocesses via `/bin/bash` |
| `bolt/common/pcgr.py: prepare_vcf_somatic` / `prepare_vcf_germline` | Shells out to `bcftools index`/`bcftools view`/`bcftools annotate` |
| `bolt/common/pcgr.py: run_somatic` / `run_somatic_chunk` (execution path) / `run_germline` | Invoke `pcgr`/`cpsr` CLI directly |
| `bolt/common/pcgr.py: transfer_annotations_somatic` / `transfer_annotations_germline` | Depend on real PCGR/CPSR TSV+VCF output files |
| `bolt/common/pcgr.py: merging_pcgr_files` | Wraps `merge_vcf_files`/`merge_tsv_files` (bcftools-dependent) |
| `bolt/workflows/smlv_somatic/annotate.py` | Orchestrates vcfanno + PON + PCGR subprocess pipeline; no unit tests |
| `bolt/workflows/smlv_somatic/rescue.py` (all functions except the header check) | SAGE hotspot recall orchestrates `bcftools isec`/`concat`/`annotate` subprocesses end-to-end; no unit tests |
| `bolt/workflows/smlv_somatic/prepare.py` | bcftools-based VCF prep; no unit tests |
| `bolt/workflows/smlv_germline/prepare.py` | bcftools-based panel region selection; no unit tests |
| `bolt/workflows/smlv_germline/report.py` | bcftools stats + CPSR report generation; no unit tests |
| `bolt/workflows/sv_somatic/annotate.py` | snpEff subprocess annotation; no unit tests |
| `bolt/workflows/sv_somatic/prioritise.py` | Wraps `external/prioritize_sv.py`; no unit tests |
| `bolt/workflows/other/cancer_report.py` | Invokes `gpgr.R canrep` subprocess; no unit tests |
| `bolt/workflows/other/multiqc_report.py` | Invokes `multiqc` subprocess; no unit tests |
| `bolt/workflows/other/purple_baf_plot.py` | Invokes `circos` subprocess; no unit tests |
| `bolt/external/prioritize_sv.py` | Vendored third-party AstraZeneca SV annotation logic; no unit tests |

Contributions raising coverage for these are welcome, but will require mocking
subprocess calls (as done in `tests/test_smlv_somatic_report.py` for
`select_pcgr_variants`) or dedicated fixture VCFs/TSVs under `tests/fixtures/`.
