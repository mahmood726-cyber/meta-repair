# Code Review Findings: MetaRepair

**Reviewer**: Claude Opus 4.6 (1M context)
**Date**: 2026-04-03
**Files reviewed**: `pipeline.py` (655 lines), `index.html` (45 lines)

## P0 (Critical) -- 1 found

### P0-1: CSV output has no formula injection protection
- **File**: `pipeline.py`, lines 612-635
- **Issue**: CSV DictWriter output includes `review_id` and `pathologies` fields which could theoretically contain user-controllable strings (review IDs from filenames, pathology labels joined with `;`). While the current data is Cochrane IDs (e.g., "CD001234"), defensive protection should be added for portability.
- **Fix**: Add a `csv_safe()` helper that prepends `'` to cell values starting with `=`, `+`, `@`, `\t`, or `\r`.

## P1 (Important) -- 2 found

### P1-1: REML Fisher-scoring uses approximation, not true REML gradient
- **File**: `pipeline.py`, lines 199-214
- **Issue**: The REML update step uses `num = sum(wr^2 * (yi - tr)^2) - swr` and `den = sum(wr^2)`. This is a Paule-Mandel-like iteration, not the standard REML Fisher scoring (`d(loglik)/d(tau2)`). Results will be close but not identical to R's `metafor::rma(method="REML")`.
- **Impact**: Low -- the instability diagnostic uses this as one of 8 specs, so small deviations average out.

### P1-2: `safe_float` / zero handling
- **File**: `pipeline.py`, line 177
- **Issue**: `C = float(sw - np.sum(wi**2) / sw)` can be 0 when all studies have equal variance, leading to `tau2 = 0` by guard. This is correct behavior but the guard `if C > 0` means C=0 silently defaults to tau2=0, which is the correct DL behavior.
- **Status**: Acceptable.

## P2 (Minor) -- 2 found

### P2-1: EM algorithm uses fixed 20 iterations without convergence check
- **File**: `pipeline.py`, lines 56-73
- **Issue**: The GMM EM algorithm runs exactly 20 iterations. Adding convergence check (delta < epsilon) would be more robust.

### P2-2: PET-PEESE decision rule uses `abs(b0_pet) < 0.05`
- **File**: `pipeline.py`, line 321
- **Issue**: This uses an absolute threshold on the intercept magnitude, not a statistical test (t-test of b0). The standard PET-PEESE conditional uses statistical significance of the PET intercept.

## Summary
- P0: 1 | P1: 2 | P2: 2
- P0-1 FIXED: CSV injection protection added.
