# ADR-001: Cap somatic variants at 450k before PCGR (not 500k)

**Date:** 2026-07-22  
**Status:** Accepted  
**Context:** bolt #35, sash #52  
**Deciders:** Team (oral decision)

## Context

PCGR has a 500,000 variant threshold that triggers two problematic behaviours:

1. **Python side** (`pcgr/main.py`, `pcgr_vars.MAX_VARIANTS_FOR_REPORT = 500_000`):  
   When input variants exceed 500k, PCGR silently drops intergenic, intronic, upstream_gene, and downstream_gene variants. This makes the TMB calculation an underestimate and removes our control over which variants appear in the report.

2. **R side** (`pcgrr/R/main.R` ~line 954):  
   If variants are still ≥ 500k after the Python-side filtering, the HTML report is **silently not generated** — no error, just no output file.

Neither triggers a hard failure. Both are silent. This makes the 500k boundary dangerous: a sample could pass through bolt, enter PCGR, and produce either a misleading report (missing variants, wrong TMB) or no report at all — with no error in the logs.

## Decision

Set `MAX_SOMATIC_VARIANTS = 450_000` in bolt. This is the threshold used by `select_pcgr_variants()` for tiered filtering and by `count_variant_process()` for the `is_hypermutated` flag.

### Why 450k and not 500k?

1. **Bolt's filter is coarse-grained.** It drops entire variant categories at once (e.g. all NONCODING_INTERGENIC variants). Output is always ≤ `MAX_SOMATIC_VARIANTS`, but it can undershoot significantly.

2. **pcgrr's boundary check is `< 500000` (strict less-than).** If bolt outputs exactly 500,000 variants, pcgrr would skip HTML report generation. The 50k margin ensures we never land on this boundary.

3. **PCGR's Python filter uses `> 500000`.** So even 500,001 variants trigger silent intergenic/intronic removal. The margin gives a buffer against any minor variant-count inflation during PCGR's annotation pipeline.

### Why not lower (e.g. 400k)?

Lower thresholds would cause more samples to be flagged as hypermutated and have their variants filtered, potentially losing clinically relevant variants in TIER_1/2 categories for samples that PCGR could actually handle fine.

## Consequences

- Samples with > 450k PASS somatic variants get tiered filtering (lowest-priority categories dropped first until under 450k)
- If tiered filtering cannot bring the count below 450k (because retained variants — hotspots, panel genes — alone exceed it), PCGR is skipped entirely for that sample
- Real case: L2100242 had 595k PASS variants; retained variants alone exceeded 450k → PCGR skipped (sash #52)
- The PCGR HTML report, MAF output, and VCF2MAF are lost for skipped samples; cancer report (gpgr) is unaffected

## Alternatives Considered

| Option | Pros | Cons |
|--------|------|------|
| 500k (match PCGR exactly) | Fewer samples filtered | Risk of silent report loss at boundary; TMB underestimates |
| 450k (chosen) | Safe margin; bolt controls filtering priority | Slightly more aggressive filtering on borderline samples |
| 400k | Extra safety | Unnecessary — 50k margin already covers all known edge cases |
| No limit (let PCGR handle it) | Simplest | PCGR's filtering is indiscriminate (drops all intergenic regardless of tier); TMB affected; report may vanish |

## References

- PCGR source: `sigven/pcgr`, `pcgr/main.py` ~line 559, `pcgrr/R/main.R` ~line 954
- `pcgr_vars.MAX_VARIANTS_FOR_REPORT = 500_000`
- sash #52: hypermutated sample handling
- bolt #26: single-chunk PCGR merge guard
- Real failure: L2100242 (595k PASS variants, PANEL+hotspot retained variants alone > 450k)
