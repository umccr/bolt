# bolt changelog

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

## 0.3.0

- [17](https://github.com/umccr/bolt/pull/17) - change dragen HRD file optional

- [14](https://github.com/umccr/bolt/pull/14) - gpgr version bump to 2.2.0

- [3](https://github.com/scwatts/bolt/pull/3) - Improve PCGR / CPSR argument handling

- [6](https://github.com/umccr/bolt/pull/6) - Change oncoanalyser v2.0.0 update, with switch sv caller from GRIPSS to eSVee

- [9](https://github.com/umccr/bolt/pull/9) Add hypermutation sample handling
