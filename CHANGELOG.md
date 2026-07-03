# bolt changelog

## 0.3.2

- [32](https://github.com/umccr/bolt/pull/32) - Graceful PCGR skip when `select_pcgr_variants` cannot cap variants to `MAX_SOMATIC_VARIANTS` — pipeline logs a warning and continues without the cancer report rather than hard-failing ([sash#52](https://github.com/umccr/sash/issues/52))
- Fix `PCGR_MUTATION_HOTSPOT=.` (dot placeholder) being treated as a retained hotspot variant during tiered filtering — was preventing tier-based filtering from running on any sample with >450k PASS variants
- Disable `--estimate_msi`/`--estimate_tmb` in chunked PCGR annotation runs (`split_vcf`) — estimates are not meaningful per-chunk
- Pin `jlumbroso/free-disk-space` CI action to `v1.3.0` (was `@main`)

## 0.3.1

- Fix `ModuleNotFoundError: No module named 'pkg_resources'` in `bolt:0.3.0-multiqc` — add `setuptools <81` to conda env
- Fix `merge_vcf_files` producing wrong output filename — `Path.with_suffix()` was stripping `.pass` component; use explicit path concatenation instead
- Fix VCF writers not closed in `transfer_annotations_somatic` and `transfer_annotations_germline` — BGZip output could be truncated
- Fix `split_vcf` writing uncompressed plain `.vcf` chunks — now uses `.vcf.gz` with `wz` mode
- Fix `PCGR_ACTIONABILITY_TIER` VCF header description — updated to match stored short-form values (`1`,`2`,`3`,`4`,`N`)
- Add regression test for chunk file compression (`test_chunks_are_gzipped`)
- Add CI smoke tests to catch Docker image startup failures before push
- Fix `build.yaml` and `Dockerfile.pcgr` build issues
- Remove unused `logging` import and fix `PCGR_MAX_SOMATIC_VARIANTS` header description in `constants.py`
- Fix `split_vcf` chunks not tabix-indexed, causing PCGR to fail reading them; add regression test
- Bump `r-gpgr` to 2.3.1 in `Dockerfile.gpgr`
- [31](https://github.com/umccr/bolt/pull/31) - Drop `--estimate_signatures` from PCGR somatic invocation — signature analysis comes from gpgr/sigrap downstream; keep `--estimate_msi`/`--estimate_tmb` ([sash#57](https://github.com/umccr/sash/issues/57))
- [32](https://github.com/umccr/bolt/pull/32) - Fix `PCGR_MUTATION_HOTSPOT=.` treated as truthy in retention check, preventing tiered filtering from running for all hypermutated samples ([sash#52](https://github.com/umccr/sash/issues/52))
- [32](https://github.com/umccr/bolt/pull/32) - Fix graceful PCGR skip when tiered filtering cannot bring PASS count below `MAX_SOMATIC_VARIANTS` — log warning and continue; non-PCGR outputs still publish ([sash#52](https://github.com/umccr/sash/issues/52))
- Fix: disable MSI/TMB estimates when running PCGR on annotation chunks — estimates on partial VCFs are meaningless
- Fix: pin `jlumbroso/free-disk-space` CI action to `v1.3.0` instead of `@main`

## 0.3.0

- [17](https://github.com/umccr/bolt/pull/17) - change dragen HRD file optional

- [14](https://github.com/umccr/bolt/pull/14) - gpgr version bump to 2.2.0

- [3](https://github.com/scwatts/bolt/pull/3) - Improve PCGR / CPSR argument handling

- [6](https://github.com/umccr/bolt/pull/6) - Change oncoanalyser v2.0.0 update, with switch sv caller from GRIPSS to eSVee

- [9](https://github.com/umccr/bolt/pull/9) Add hypermutation sample handling
