# bolt changelog

## 0.3.2

- [32](https://github.com/umccr/bolt/pull/32) - Fix `PCGR_MUTATION_HOTSPOT=.` (dot placeholder) treated as truthy in retention check — was preventing tiered filtering from running for any sample with >450k PASS variants ([sash#52](https://github.com/umccr/sash/issues/52))
- [32](https://github.com/umccr/bolt/pull/32) - Graceful PCGR skip when `select_pcgr_variants` cannot cap variants to `MAX_SOMATIC_VARIANTS` — logs warning and continues without cancer report; non-PCGR outputs still publish ([sash#52](https://github.com/umccr/sash/issues/52))
- Fix: disable `--estimate_msi`/`--estimate_tmb` in chunked PCGR annotation runs — estimates on partial VCFs are not meaningful
- Fix: pin `jlumbroso/free-disk-space` CI action to `v1.3.0` (was `@main`)
- Test: `TestSelectPcgrVariants` — 8 integration tests covering tiered trimming, hotspot retention, PANEL retention, NONCODING-first drop order, and the `PCGR_MUTATION_HOTSPOT=.` regression
- Test: `TestEntrySkipsPcgrOnOverflow` — 2 tests: entry() skips PCGR on `RuntimeError` from unresolvable overflow; entry() calls PCGR normally when within limit
- Test: `TestSelectPcgrVariantsRaisesOnUnresolvableOverflow` — asserts `RuntimeError` when retained variants alone exceed `MAX_SOMATIC_VARIANTS`
- Test: `TestRunSomaticCommandArgs` — 2 tests: `--estimate_signatures` absent from all `run_somatic` commands; `disable_estimates=True` suppresses `--estimate_msi`/`--estimate_tmb`
- Test: `TestRunSomaticChunkArgMapping.test_disable_estimates_passed_to_run_somatic` — `run_somatic_chunk` passes `disable_estimates=True` to every `run_somatic` call

## 0.3.1

- Fix `ModuleNotFoundError: No module named 'pkg_resources'` in `bolt:0.3.0-multiqc` — add `setuptools <81` to conda env
- Fix `merge_vcf_files` producing wrong output filename — `Path.with_suffix()` was stripping `.pass` component; use explicit path concatenation instead
- Fix VCF writers not closed in `transfer_annotations_somatic` and `transfer_annotations_germline` — BGZip output could be truncated
- Fix `split_vcf` writing uncompressed plain `.vcf` chunks — now uses `.vcf.gz` with `wz` mode
- Fix `PCGR_ACTIONABILITY_TIER` VCF header description — updated to match stored short-form values (`1`,`2`,`3`,`4`,`N`)
- Fix `split_vcf` chunks not tabix-indexed, causing PCGR to fail reading them
- Fix `build.yaml` and `Dockerfile.pcgr` build issues
- Remove unused `logging` import and fix `PCGR_MAX_SOMATIC_VARIANTS` header description in `constants.py`
- Bump `r-gpgr` to 2.3.1 in `Dockerfile.gpgr`
- Add CI smoke tests to catch Docker image startup failures before push
- [31](https://github.com/umccr/bolt/pull/31) - Drop `--estimate_signatures` from PCGR somatic invocation — signature analysis comes from gpgr/sigrap downstream; keep `--estimate_msi`/`--estimate_tmb` ([sash#57](https://github.com/umccr/sash/issues/57))
- Test: regression test for chunk file compression (`test_chunks_are_gzipped`)
- Test: regression test for chunk tabix indexing (`test_chunks_are_tabix_indexed`)
- Test: `TestTierOrdering` — 3 tests verifying `PCGR_TIERS_FILTERING` uses short forms (`N`,`4`,`3`,`2`,`1`) and NONCODING precedes TIER_1 in `get_ordering()`
- Test: `TestSplitVcf` — 4 tests: chunking above/below limit, `.vcf.gz` compression, `.tbi` indexing
- Test: `TestSelectPcgrVariants` — initial 6 integration tests for tiered trimming logic
- Test: `TestGetVariantFilterData`, `TestDetermineFilter`, `TestGetImpacts` — 14 unit tests covering variant attribute extraction and filter-category determination
- Test: `TestCountVariantProcess` — 4 tests: `is_hypermutated` flag, DRAGEN count, SAGE_NOVEL exclusion, annotation-filter exclusion
- Test: `TestRunSomaticChunkArgMapping` — asserts `pcgr_conda` is not shifted into `pcgr_threads` position on positional arg mapping

## 0.3.0

- [28](https://github.com/umccr/bolt/pull/28) - gpgr version bump to 2.2.12 for cancer report hypermutated flag fix

- [17](https://github.com/umccr/bolt/pull/17) - change dragen HRD file optional
- [14](https://github.com/umccr/bolt/pull/14) - gpgr version bump to 2.2.0
- [3](https://github.com/scwatts/bolt/pull/3) - Improve PCGR / CPSR argument handling
- [6](https://github.com/umccr/bolt/pull/6) - Change oncoanalyser v2.0.0 update, with switch sv caller from GRIPSS to eSVee
- [9](https://github.com/umccr/bolt/pull/9) Add hypermutation sample handling
