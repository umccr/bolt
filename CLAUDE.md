# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install locally (requires conda env with binary deps for full functionality)
pip install -e .

# Run tests (canonical — matches CI)
python -m unittest discover tests/ --buffer

# Run tests with pytest (also works)
python -m pytest tests/ -v

# Run a single test
python -m pytest tests/test_smlv_somatic_filter.py::TestSmlvSomaticFilter::test_min_af_filter -v

# Run bolt CLI directly
python -m bolt smlv_somatic filter --help
```

CI triggers: push to any branch runs tests; push of a `v*.*.*` tag builds and pushes all 6 Docker images to `ghcr.io/umccr/bolt`.

## Architecture

`bolt` is a Click CLI toolkit for UMCCR WGS somatic/germline post-processing. Entry point: `bolt/__main__.py` → `bolt/workflows/cli.py` auto-discovers workflow groups from `bolt/workflows/*/`.

Each workflow directory (`smlv_somatic`, `smlv_germline`, `sv_somatic`, `other`) contains command modules. Each module exposes a `click` command named `entry`. Commands run as steps in the `sash` Nextflow pipeline, invoked via Docker images.

**Six Docker images** split heavy dependencies — each command maps to exactly one image (see README). CI builds all six on `v*.*.*` tag push.

## Key files

| File | Role |
|------|------|
| `bolt/common/constants.py` | All thresholds, VCF tag enums (`VcfFilter`, `VcfInfo`, `VcfFormat`), and header definitions. Authoritative — add new tags here first. |
| `bolt/common/pcgr.py` | PCGR/CPSR invocation, VCF splitting/merging for hypermutated chunking, annotation transfer |
| `bolt/util.py` | `execute_command`, `command_prepare`, `add_vcf_header_entry`, VCF counting/merging helpers |
| `bolt/workflows/smlv_somatic/annotate.py` | PCGR annotation with chunking for hypermutated samples |
| `bolt/workflows/smlv_somatic/report.py` | PCGR cancer report, hypermutated variant selection (`select_pcgr_variants`) |
| `bolt/workflows/smlv_somatic/filter.py` | Somatic variant filtering — representative of the VCF processing pattern |

## Coding patterns

**CLI:** Each command uses `entry(ctx, **kwargs)` with `@click.option` decorators. Call `setup_logging(output_dir, script_name)` early in `entry`.

**VCF writes:** Always use `cyvcf2.Writer(output_fp, in_fh, 'wz')` and close writers explicitly. Add new INFO/FILTER tags via `util.add_vcf_header_entry(fh, constants.VcfInfo.MY_TAG)` — never hardcode header strings.

**Shell commands:** Use `util.execute_command(cmd)` for single commands. For pipelines where `pipefail` matters, wrap with `util.command_prepare(cmd)` first.

## PCGR hypermutated handling

PCGR enforces a hard 500k variant limit in software (not an OOM — it raises an explicit error). `MAX_SOMATIC_VARIANTS = 450_000` (`constants.py`) is the safe threshold used throughout bolt.

- **Annotation path** (`annotate.py`): `split_vcf()` → chunks ≤450k → `run_somatic_chunk()` → merge. Handles any input size. ✅ Tested with synthetic 550k VCF — correctly splits into 450k + 100k chunks (2026-05-13).
- **Report path** (`report.py`): `select_pcgr_variants()` filters down to ≤450k via tiered priority (NONCODING first → TIER_1 last), then single `run_somatic()`. Requires `PCGR_ACTIONABILITY_TIER` (PCGR v2.2.1+) — if missing, variants have tier=None and fall into category `(None, impact, region)` which is absent from `get_ordering()`, so nothing gets dropped. If tiered filtering cannot reach ≤450k, raises `RuntimeError` — pipeline fails rather than hitting the PCGR hard software limit (not OOM).
  - `PCGR_TIERS_FILTERING` uses short forms `('N','4','3','2','1')` matching values written by `transfer_annotations_somatic()`.
- `RETAIN_FIELDS_FILTERING` = `('PANEL', *HOTSPOT_FIELDS_FILTERING)` marks variants never dropped by tiered filtering.

**Known failure mode (sash #52):** if a >500k sample has too few NONCODING/low-tier variants to drop, `select_pcgr_variants` raises `RuntimeError`. Fix not yet implemented — tracked in sash #52.

## Release chain

gpgr (conda publish) → `bolt/docker/Dockerfile.gpgr` (pin `r-gpgr ==X.Y.Z`) → bolt tag → 6 Docker images → sash module container tags. Runbook: `~/Documents/UMCCR/Runbooks/Release gpgr → bolt → sash.md`.

Dev image pattern: install gpgr branch via `install_github` in `Dockerfile.gpgr` (marked `# DEV: … revert before release`). Build with `--platform linux/amd64` on Apple Silicon — `linux-aarch64` conda channel lacks `bcftools==1.17`.
