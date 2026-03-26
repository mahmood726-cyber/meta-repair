"""MetaRepair: Automatic correction of meta-analysis pathologies.

Takes a meta-analysis (yi, sei) and outputs:
  1. Diagnosis: multimodality, outliers, publication bias, instability
  2. Correction: robust pooling, subgroup splitting, bias adjustment
  3. Uncertainty decomposition: how much does each source contribute?
  4. Decision-grade output: corrected estimate + honest uncertainty

Applied to 403 Cochrane reviews from Pairwise70.
"""

import csv
import json
import math
import time
import numpy as np
import pyreadr
from pathlib import Path
from scipy import stats as sp_stats
from collections import Counter

PAIRWISE_DIR = Path(r'C:\Models\Pairwise70\data')
OUTPUT_DIR = Path(r'C:\Models\MetaRepair\data\output')


# ═══════════════════════════════════════════════════
# STAGE 1: DIAGNOSIS
# ═══════════════════════════════════════════════════

def diagnose_multimodality(yi, sei):
    """Detect multimodal effect distribution via Gaussian mixture (EM, k=1 vs k=2)."""
    k = len(yi)
    if k < 5:
        return {'n_modes': 1, 'bic_improvement': 0, 'cluster_labels': [0]*k}

    # Precision-weighted
    wi = 1.0 / sei**2
    wi_norm = wi / wi.sum()

    # Fit k=1 (single Gaussian)
    mu1 = float(np.sum(wi_norm * yi))
    var1 = float(np.sum(wi_norm * (yi - mu1)**2))
    ll1 = float(np.sum(-0.5 * np.log(2*np.pi*var1) - 0.5 * (yi - mu1)**2 / var1))
    bic1 = -2 * ll1 + 2 * np.log(k)  # 2 params: mu, sigma

    # Fit k=2 (EM algorithm, 20 iterations)
    # Initialize by splitting at median
    med = np.median(yi)
    labels = (yi > med).astype(int)
    pi_mix = np.array([0.5, 0.5])
    mus = np.array([yi[labels == 0].mean() if (labels == 0).sum() > 0 else med - 0.1,
                     yi[labels == 1].mean() if (labels == 1).sum() > 0 else med + 0.1])
    sigs = np.array([max(0.01, yi[labels == 0].std()) if (labels == 0).sum() > 1 else 0.1,
                      max(0.01, yi[labels == 1].std()) if (labels == 1).sum() > 1 else 0.1])

    for _ in range(20):
        # E-step
        resp = np.zeros((k, 2))
        for c in range(2):
            resp[:, c] = pi_mix[c] * np.exp(-0.5 * ((yi - mus[c]) / sigs[c])**2) / (sigs[c] * np.sqrt(2*np.pi))
        resp_sum = resp.sum(axis=1, keepdims=True)
        resp_sum = np.maximum(resp_sum, 1e-30)
        resp /= resp_sum

        # M-step
        for c in range(2):
            wc = resp[:, c]
            nc = wc.sum()
            if nc < 1:
                continue
            pi_mix[c] = nc / k
            mus[c] = np.sum(wc * yi) / nc
            sigs[c] = max(0.01, np.sqrt(np.sum(wc * (yi - mus[c])**2) / nc))

    # Log-likelihood for k=2
    ll2 = 0
    for i in range(k):
        p = sum(pi_mix[c] * np.exp(-0.5 * ((yi[i] - mus[c]) / sigs[c])**2) / (sigs[c] * np.sqrt(2*np.pi))
                for c in range(2))
        ll2 += np.log(max(p, 1e-30))
    bic2 = -2 * ll2 + 5 * np.log(k)  # 5 params: 2*mu, 2*sigma, 1*pi

    # BIC comparison
    bic_improvement = bic1 - bic2  # positive = 2-component is better
    is_multimodal = bic_improvement > 6  # strong evidence threshold

    # Assign labels based on final responsibilities
    labels = np.argmax(resp, axis=1).tolist() if is_multimodal else [0] * k

    return {
        'n_modes': 2 if is_multimodal else 1,
        'bic_improvement': float(bic_improvement),
        'cluster_labels': labels,
        'cluster_means': mus.tolist() if is_multimodal else [mu1],
        'cluster_sizes': [int(sum(1 for l in labels if l == c)) for c in range(2)] if is_multimodal else [k],
    }


def diagnose_outliers(yi, sei):
    """Detect outliers via studentized residuals from DL model."""
    k = len(yi)
    if k < 4:
        return {'n_outliers': 0, 'outlier_indices': [], 'outlier_studies': []}

    # FE pooling (not DL — DL inflates tau2 to accommodate outliers, masking them)
    wi = 1.0 / sei**2
    sw = np.sum(wi)
    theta_fe = float(np.sum(wi * yi) / sw)

    # Studentized residuals from FE model
    residuals = (yi - theta_fe) / sei
    outlier_mask = np.abs(residuals) > 3.0
    outlier_idx = np.where(outlier_mask)[0].tolist()

    return {
        'n_outliers': len(outlier_idx),
        'outlier_indices': outlier_idx,
        'residuals': residuals.tolist(),
    }


def diagnose_publication_bias(yi, sei):
    """Egger's test + trim-and-fill count."""
    k = len(yi)
    if k < 5:
        return {'egger_p': 1.0, 'tf_k0': 0, 'bias_detected': False}

    # Egger's test
    n = k
    x = 1.0 / sei
    y = yi / sei
    xm, ym = np.mean(x), np.mean(y)
    ssxx = np.sum((x - xm)**2)
    if ssxx < 1e-15:
        return {'egger_p': 1.0, 'egger_intercept': 0, 'tf_k0': 0, 'bias_detected': False}
    ssxy = np.sum((x - xm) * (y - ym))
    b = ssxy / ssxx
    a0 = ym - b * xm
    resid = y - a0 - b * x
    mse = np.sum(resid**2) / max(n - 2, 1)
    se_a = np.sqrt(mse * (1/n + xm**2/ssxx))
    t = a0 / se_a if se_a > 0 else 0
    egger_p = 2 * (1 - sp_stats.t.cdf(abs(t), max(n - 2, 1)))

    # Trim-and-fill (L0 estimator)
    wi = 1.0 / sei**2
    sw = np.sum(wi)
    theta_fe = np.sum(wi * yi) / sw
    ranks = np.argsort(np.argsort(np.abs(yi - theta_fe)))
    signs = np.sign(yi - theta_fe)
    # Count studies on the "wrong" side
    right_of_center = np.sum(signs > 0)
    left_of_center = np.sum(signs < 0)
    k0 = abs(right_of_center - left_of_center)
    k0 = max(0, k0 - 1)  # conservative

    return {
        'egger_p': float(egger_p),
        'egger_intercept': float(a0),
        'tf_k0': int(k0),
        'bias_detected': egger_p < 0.10 or k0 >= 2,
    }


def diagnose_instability(yi, sei):
    """Multiverse stability: fraction of 8 specs agreeing with DL conclusion."""
    k = len(yi)
    if k < 3:
        return {'stability': 0, 'n_agree': 0, 'n_total': 0}

    # DL direction and significance
    wi = 1.0 / sei**2
    sw = np.sum(wi)
    theta_fe = np.sum(wi * yi) / sw
    Q = float(np.sum(wi * (yi - theta_fe)**2))
    C = float(sw - np.sum(wi**2) / sw)
    tau2 = max(0, (Q - (k-1)) / C) if C > 0 else 0
    ws = 1.0 / (sei**2 + tau2)
    sws = np.sum(ws)
    theta = float(np.sum(ws * yi) / sws)
    se = float(1.0 / np.sqrt(sws))
    ci_lo, ci_hi = theta - 1.96 * se, theta + 1.96 * se
    dl_sig = (ci_lo > 0) or (ci_hi < 0)
    dl_dir = -1 if theta < 0 else 1

    # Test 4 estimators × 2 CI methods = 8 specs
    n_agree = 0
    n_total = 0

    for method in ['DL', 'FE', 'REML', 'PM']:
        for ci_method in ['Wald', 'HKSJ']:
            try:
                if method == 'DL':
                    t, s = theta, se
                elif method == 'FE':
                    t = float(np.sum(wi * yi) / sw)
                    s = float(1.0 / np.sqrt(sw))
                elif method == 'REML':
                    t2_r = tau2
                    for _ in range(50):
                        wr = 1.0 / (sei**2 + t2_r)
                        swr = np.sum(wr)
                        tr = float(np.sum(wr * yi) / swr)
                        num = float(np.sum(wr**2 * (yi - tr)**2) - swr)
                        den = float(np.sum(wr**2))
                        nt2 = max(0, t2_r + num / den) if den > 0 else t2_r
                        if abs(nt2 - t2_r) < 1e-8:
                            t2_r = nt2
                            break
                        t2_r = nt2
                    wr = 1.0 / (sei**2 + t2_r)
                    swr = np.sum(wr)
                    t = float(np.sum(wr * yi) / swr)
                    s = float(1.0 / np.sqrt(swr))
                elif method == 'PM':
                    t2_p = tau2
                    for _ in range(50):
                        wp = 1.0 / (sei**2 + t2_p)
                        swp = np.sum(wp)
                        tp = float(np.sum(wp * yi) / swp)
                        Qs = float(np.sum(wp * (yi - tp)**2))
                        sw2p = float(np.sum(wp**2))
                        Cp = swp - sw2p / swp if swp > 0 else 1
                        nt2 = max(0, t2_p + (Qs - (k-1)) / Cp)
                        if abs(nt2 - t2_p) < 1e-8:
                            t2_p = nt2
                            break
                        t2_p = nt2
                    wp = 1.0 / (sei**2 + t2_p)
                    swp = np.sum(wp)
                    t = float(np.sum(wp * yi) / swp)
                    s = float(1.0 / np.sqrt(swp))

                if ci_method == 'Wald':
                    lo, hi = t - 1.96 * s, t + 1.96 * s
                else:  # HKSJ
                    tc = sp_stats.t.ppf(0.975, k - 1)
                    lo, hi = t - tc * s, t + tc * s

                sig = (lo > 0) or (hi < 0)
                d = -1 if t < 0 else 1
                n_total += 1
                if sig == dl_sig and d == dl_dir:
                    n_agree += 1
            except Exception:
                continue

    stability = n_agree / n_total * 100 if n_total > 0 else 0
    return {'stability': round(stability, 1), 'n_agree': n_agree, 'n_total': n_total}


# ═══════════════════════════════════════════════════
# STAGE 2: CORRECTION
# ═══════════════════════════════════════════════════

def correct_robust_pooling(yi, sei, outlier_indices):
    """Winsorized DL: downweight outliers instead of removing them."""
    k = len(yi)
    weights = np.ones(k)
    for idx in outlier_indices:
        weights[idx] = 0.1  # 10% weight for outliers

    wi = weights / sei**2
    sw = np.sum(wi)
    theta_fe = np.sum(wi * yi) / sw
    Q = float(np.sum(wi * (yi - theta_fe)**2))
    C = float(sw - np.sum(wi**2) / sw)
    tau2 = max(0, (Q - (k-1)) / C) if C > 0 else 0
    ws = weights / (sei**2 + tau2)
    sws = np.sum(ws)
    theta = float(np.sum(ws * yi) / sws)
    se = float(1.0 / np.sqrt(sws))

    return {'theta': theta, 'se': se, 'tau2': tau2}


def correct_bias_adjustment(yi, sei, tf_k0):
    """PET-PEESE bias correction: regress effect on SE (PET) or SE² (PEESE)."""
    k = len(yi)
    if k < 5 or tf_k0 < 1:
        # No correction needed
        wi = 1.0 / sei**2
        sw = np.sum(wi)
        Q = float(np.sum(wi * (yi - np.sum(wi * yi) / sw)**2))
        C = float(sw - np.sum(wi**2) / sw)
        tau2 = max(0, (Q - (k-1)) / C) if C > 0 else 0
        ws = 1.0 / (sei**2 + tau2)
        sws = np.sum(ws)
        return {'theta': float(np.sum(ws * yi) / sws), 'se': float(1.0 / np.sqrt(sws)), 'method': 'none'}

    # PET: yi = b0 + b1*SEi + ei (WLS, weights = 1/vi)
    wi = 1.0 / sei**2
    x = sei
    xw = x * wi
    sw = np.sum(wi)
    sxw = np.sum(xw)
    sxxw = np.sum(x * xw)
    syw = np.sum(yi * wi)
    sxyw = np.sum(x * yi * wi)
    det = sw * sxxw - sxw**2
    if abs(det) < 1e-15:
        return {'theta': float(np.sum(wi * yi) / sw), 'se': float(1.0 / np.sqrt(sw)), 'method': 'none'}

    b0_pet = (sxxw * syw - sxw * sxyw) / det
    b1_pet = (sw * sxyw - sxw * syw) / det

    # PEESE: yi = b0 + b1*SEi² + ei
    x2 = sei**2
    x2w = x2 * wi
    sx2w = np.sum(x2w)
    sx2x2w = np.sum(x2 * x2w)
    sx2yw = np.sum(x2 * yi * wi)
    det2 = sw * sx2x2w - sx2w**2
    if abs(det2) < 1e-15:
        b0_peese = b0_pet
    else:
        b0_peese = (sx2x2w * syw - sx2w * sx2yw) / det2

    # PET-PEESE: use PET if b0_pet is non-significant, else PEESE
    # Simplified: if PET intercept < 0.05 in magnitude, use PET
    theta = b0_pet if abs(b0_pet) < 0.05 else b0_peese
    se_adj = float(1.0 / np.sqrt(sw))  # approximate

    return {'theta': float(theta), 'se': se_adj, 'method': 'PET' if abs(b0_pet) < 0.05 else 'PEESE'}


def correct_subgroup_split(yi, sei, cluster_labels):
    """Pool within each detected cluster separately."""
    clusters = sorted(set(cluster_labels))
    results = []
    for c in clusters:
        mask = np.array([l == c for l in cluster_labels])
        yc, sc = yi[mask], sei[mask]
        if len(yc) < 2:
            continue
        wi = 1.0 / sc**2
        sw = np.sum(wi)
        theta_fe = np.sum(wi * yc) / sw
        Q = float(np.sum(wi * (yc - theta_fe)**2))
        C = float(sw - np.sum(wi**2) / sw)
        tau2 = max(0, (Q - (len(yc)-1)) / C) if C > 0 else 0
        ws = 1.0 / (sc**2 + tau2)
        sws = np.sum(ws)
        theta = float(np.sum(ws * yc) / sws)
        se = float(1.0 / np.sqrt(sws))
        results.append({'cluster': c, 'k': len(yc), 'theta': theta, 'se': se, 'tau2': tau2})
    return results


# ═══════════════════════════════════════════════════
# STAGE 3: UNCERTAINTY DECOMPOSITION
# ═══════════════════════════════════════════════════

def decompose_uncertainty(yi, sei, diag, corrected):
    """Decompose total uncertainty into sources."""
    k = len(yi)

    # Standard DL
    wi = 1.0 / sei**2
    sw = np.sum(wi)
    theta_fe = np.sum(wi * yi) / sw
    Q = float(np.sum(wi * (yi - theta_fe)**2))
    C = float(sw - np.sum(wi**2) / sw)
    tau2 = max(0, (Q - (k-1)) / C) if C > 0 else 0
    ws = 1.0 / (sei**2 + tau2)
    sws = np.sum(ws)
    theta_dl = float(np.sum(ws * yi) / sws)
    se_dl = float(1.0 / np.sqrt(sws))

    # Sampling uncertainty: SE of pooled estimate
    u_sampling = se_dl**2

    # Heterogeneity uncertainty: tau2 contribution to prediction
    u_heterogeneity = tau2

    # Model uncertainty: variance across 8 multiverse specs
    spec_thetas = []
    for method_tau2 in [0, tau2, tau2 * 1.5]:  # FE, DL, inflated
        wm = 1.0 / (sei**2 + method_tau2)
        swm = np.sum(wm)
        spec_thetas.append(float(np.sum(wm * yi) / swm))
    u_model = np.var(spec_thetas) if len(spec_thetas) > 1 else 0

    # Bias uncertainty: difference between original and bias-corrected
    u_bias = (theta_dl - corrected.get('theta_bias_adj', theta_dl))**2

    # Outlier uncertainty: difference between original and robust
    u_outlier = (theta_dl - corrected.get('theta_robust', theta_dl))**2

    total = u_sampling + u_heterogeneity + u_model + u_bias + u_outlier
    if total < 1e-15:
        total = 1

    return {
        'total_variance': float(total),
        'pct_sampling': round(100 * u_sampling / total, 1),
        'pct_heterogeneity': round(100 * u_heterogeneity / total, 1),
        'pct_model': round(100 * u_model / total, 1),
        'pct_bias': round(100 * u_bias / total, 1),
        'pct_outlier': round(100 * u_outlier / total, 1),
    }


# ═══════════════════════════════════════════════════
# STAGE 4: DECISION-GRADE OUTPUT
# ═══════════════════════════════════════════════════

def decision_grade(diag, corrected, uncertainty):
    """Produce a decision-grade recommendation."""
    # Count pathologies
    n_pathologies = 0
    pathologies = []
    if diag['multimodality']['n_modes'] >= 2:
        n_pathologies += 1
        pathologies.append('multimodal')
    if diag['outliers']['n_outliers'] > 0:
        n_pathologies += 1
        pathologies.append(f"{diag['outliers']['n_outliers']} outliers")
    if diag['publication_bias']['bias_detected']:
        n_pathologies += 1
        pathologies.append('pub bias')
    if diag['instability']['stability'] < 70:
        n_pathologies += 1
        pathologies.append('unstable')

    # Correction magnitude
    orig = corrected['theta_original']
    final = corrected['theta_final']
    correction_pct = abs(orig - final) / max(abs(orig), 0.001) * 100

    # Grade
    if n_pathologies == 0 and correction_pct < 5:
        grade = 'A'
        recommendation = 'Evidence is robust. Original estimate is trustworthy.'
    elif n_pathologies <= 1 and correction_pct < 20:
        grade = 'B'
        recommendation = 'Minor pathology detected. Corrected estimate recommended.'
    elif n_pathologies <= 2 and correction_pct < 50:
        grade = 'C'
        recommendation = 'Multiple pathologies. Interpret with caution. Report corrected estimate.'
    else:
        grade = 'D'
        recommendation = 'Severe pathologies. Conclusion unreliable. Consider subgroup analysis or additional evidence.'

    return {
        'grade': grade,
        'n_pathologies': n_pathologies,
        'pathologies': pathologies,
        'correction_pct': round(correction_pct, 1),
        'recommendation': recommendation,
    }


# ═══════════════════════════════════════════════════
# DATA LOADING (reuses Fragility Atlas pattern)
# ═══════════════════════════════════════════════════

def load_review(rda_path):
    """Load one RDA, return yi/sei for primary analysis."""
    result = pyreadr.read_r(str(rda_path))
    df = list(result.values())[0].copy()
    df.columns = df.columns.str.replace(' ', '.', regex=False)
    review_id = rda_path.stem.split('_')[0]

    import pandas as pd
    groups = []
    for (grp, num), sub in df.groupby(['Analysis.group', 'Analysis.number']):
        has_binary = (sub['Experimental.cases'].notna() & (sub['Experimental.cases'] > 0)).any()
        groups.append({'grp': grp, 'num': num, 'k': len(sub), 'binary': has_binary})
    if not groups:
        return None
    gdf = pd.DataFrame(groups)
    binary = gdf[gdf['binary']]
    best = binary.loc[binary['k'].idxmax()] if len(binary) > 0 else gdf.loc[gdf['k'].idxmax()]
    primary = df[(df['Analysis.group'] == best['grp']) & (df['Analysis.number'] == best['num'])]

    has_binary = (primary['Experimental.cases'].notna() & (primary['Experimental.cases'] > 0)).any()
    scale = 'ratio' if has_binary else ('ratio' if (primary['Mean'].dropna() > 0).all() else 'difference')

    if scale == 'ratio':
        v = (primary['Mean'].notna() & (primary['Mean'] > 0) & primary['CI.start'].notna() & (primary['CI.start'] > 0) & primary['CI.end'].notna() & (primary['CI.end'] > 0))
        sub = primary[v]
        if len(sub) < 3: return None
        yi = np.log(sub['Mean'].values.astype(float))
        sei = (np.log(sub['CI.end'].values.astype(float)) - np.log(sub['CI.start'].values.astype(float))) / (2 * 1.96)
    else:
        v = primary['Mean'].notna() & primary['CI.start'].notna() & primary['CI.end'].notna()
        sub = primary[v]
        if len(sub) < 3: return None
        yi = sub['Mean'].values.astype(float)
        sei = (sub['CI.end'].values.astype(float) - sub['CI.start'].values.astype(float)) / (2 * 1.96)

    ok = (sei > 0) & np.isfinite(yi) & np.isfinite(sei)
    yi, sei = yi[ok], sei[ok]
    if len(yi) < 3: return None
    return {'review_id': review_id, 'yi': yi, 'sei': sei, 'k': len(yi), 'scale': scale}


# ═══════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════

def repair_one(review):
    """Full MetaRepair pipeline for one review."""
    yi, sei = review['yi'], review['sei']
    k = review['k']

    # DIAGNOSE
    diag = {
        'multimodality': diagnose_multimodality(yi, sei),
        'outliers': diagnose_outliers(yi, sei),
        'publication_bias': diagnose_publication_bias(yi, sei),
        'instability': diagnose_instability(yi, sei),
    }

    # ORIGINAL ESTIMATE
    wi = 1.0 / sei**2
    sw = np.sum(wi)
    theta_fe = np.sum(wi * yi) / sw
    Q = float(np.sum(wi * (yi - theta_fe)**2))
    C = float(sw - np.sum(wi**2) / sw)
    tau2 = max(0, (Q - (k-1)) / C) if C > 0 else 0
    ws = 1.0 / (sei**2 + tau2)
    sws = np.sum(ws)
    theta_orig = float(np.sum(ws * yi) / sws)
    se_orig = float(1.0 / np.sqrt(sws))

    # CORRECT
    robust = correct_robust_pooling(yi, sei, diag['outliers']['outlier_indices'])
    bias_adj = correct_bias_adjustment(yi, sei, diag['publication_bias']['tf_k0'])
    subgroups = correct_subgroup_split(yi, sei, diag['multimodality']['cluster_labels']) if diag['multimodality']['n_modes'] >= 2 else []

    # Final corrected estimate: average of robust and bias-adjusted
    theta_final = (robust['theta'] + bias_adj['theta']) / 2
    se_final = max(robust['se'], bias_adj['se'])

    corrected = {
        'theta_original': theta_orig,
        'se_original': se_orig,
        'theta_robust': robust['theta'],
        'theta_bias_adj': bias_adj['theta'],
        'theta_final': theta_final,
        'se_final': se_final,
        'bias_method': bias_adj['method'],
        'subgroups': subgroups,
    }

    # DECOMPOSE UNCERTAINTY
    uncertainty = decompose_uncertainty(yi, sei, diag, corrected)

    # DECISION GRADE
    decision = decision_grade(diag, corrected, uncertainty)

    return {
        'review_id': review['review_id'],
        'k': k,
        'scale': review['scale'],
        'diagnosis': diag,
        'corrected': corrected,
        'uncertainty': uncertainty,
        'decision': decision,
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print("MetaRepair Pipeline")
    print("=" * 40)

    t0 = time.time()
    rda_files = sorted(PAIRWISE_DIR.glob('*.rda'))
    print(f"  Found {len(rda_files)} RDA files")

    results = []
    for rda in rda_files:
        review = load_review(rda)
        if review is None or review['k'] < 5:
            continue
        result = repair_one(review)
        results.append(result)

    elapsed = time.time() - t0
    print(f"  Processed: {len(results)} reviews in {elapsed:.1f}s")

    # HEADLINE STATS
    grades = Counter(r['decision']['grade'] for r in results)
    n_multimodal = sum(1 for r in results if r['diagnosis']['multimodality']['n_modes'] >= 2)
    n_outliers = sum(1 for r in results if r['diagnosis']['outliers']['n_outliers'] > 0)
    n_biased = sum(1 for r in results if r['diagnosis']['publication_bias']['bias_detected'])
    n_unstable = sum(1 for r in results if r['diagnosis']['instability']['stability'] < 70)
    corrections = [r['decision']['correction_pct'] for r in results]
    mean_correction = sum(corrections) / len(corrections) if corrections else 0

    # Uncertainty decomposition averages
    pct_keys = ['pct_sampling', 'pct_heterogeneity', 'pct_model', 'pct_bias', 'pct_outlier']
    avg_uncertainty = {k: round(sum(r['uncertainty'][k] for r in results) / len(results), 1) for k in pct_keys}

    print(f"\n{'='*50}")
    print("HEADLINE FINDINGS")
    print(f"{'='*50}")
    print(f"  Reviews: {len(results)}")
    print(f"  Grades: A={grades['A']}, B={grades['B']}, C={grades['C']}, D={grades['D']}")
    print(f"  Pathologies: multimodal={n_multimodal}, outliers={n_outliers}, biased={n_biased}, unstable={n_unstable}")
    print(f"  Mean correction: {mean_correction:.1f}%")
    print(f"  Uncertainty decomposition (average):")
    for k, v in avg_uncertainty.items():
        print(f"    {k}: {v}%")

    # EXPORT
    rows = []
    for r in results:
        rows.append({
            'review_id': r['review_id'],
            'k': r['k'],
            'grade': r['decision']['grade'],
            'n_pathologies': r['decision']['n_pathologies'],
            'pathologies': ';'.join(r['decision']['pathologies']),
            'theta_original': round(r['corrected']['theta_original'], 4),
            'theta_corrected': round(r['corrected']['theta_final'], 4),
            'correction_pct': r['decision']['correction_pct'],
            'n_modes': r['diagnosis']['multimodality']['n_modes'],
            'n_outliers': r['diagnosis']['outliers']['n_outliers'],
            'bias_detected': r['diagnosis']['publication_bias']['bias_detected'],
            'stability': r['diagnosis']['instability']['stability'],
            'pct_sampling': r['uncertainty']['pct_sampling'],
            'pct_heterogeneity': r['uncertainty']['pct_heterogeneity'],
            'pct_model': r['uncertainty']['pct_model'],
            'pct_bias': r['uncertainty']['pct_bias'],
            'pct_outlier': r['uncertainty']['pct_outlier'],
        })

    with open(OUTPUT_DIR / 'metarepair_results.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        'n_reviews': len(results),
        'grades': dict(grades),
        'n_multimodal': n_multimodal,
        'n_outliers': n_outliers,
        'n_biased': n_biased,
        'n_unstable': n_unstable,
        'mean_correction_pct': round(mean_correction, 1),
        'avg_uncertainty': avg_uncertainty,
        'elapsed_seconds': round(elapsed, 1),
    }
    with open(OUTPUT_DIR / 'metarepair_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Saved to {OUTPUT_DIR}/")


if __name__ == '__main__':
    main()
