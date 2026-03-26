# MetaRepair: Automatic Diagnosis and Correction of Meta-Analysis Pathologies Across 307 Cochrane Reviews

## Authors

Mahmood Ahmad^1

^1 Royal Free Hospital, London, United Kingdom

Correspondence: mahmood.ahmad2@nhs.net | ORCID: 0009-0003-7781-4478

---

## Abstract

**Objective:** To automatically diagnose, correct, and decompose the uncertainty sources in meta-analyses, producing decision-grade outputs that account for multimodality, outliers, publication bias, and methodological instability.

**Design:** Automated pipeline applied to 307 Cochrane systematic reviews.

**Data sources:** Pairwise70 dataset (501 reviews, k>=5 studies).

**Methods:** For each review, MetaRepair applies four diagnostic modules (multimodality detection via Gaussian mixture models, outlier detection via studentized residuals, publication bias via Egger's test and trim-and-fill, and methodological instability via 8-specification multiverse analysis). Corrections include winsorized robust pooling, PET-PEESE bias adjustment, and automatic subgroup splitting. Uncertainty is decomposed into five sources: sampling, heterogeneity, model choice, publication bias, and outlier influence. Reviews are graded A-D based on pathology count and correction magnitude.

**Results:** Of 307 reviews, 58 (19%) received Grade A (robust, minimal correction needed), 88 (29%) Grade B (minor pathology), 79 (26%) Grade C (multiple pathologies, interpret with caution), and 82 (27%) Grade D (severe pathologies, conclusion unreliable). The mean correction magnitude was 50.7% — indicating that half of published pooled estimates shift substantially when pathologies are addressed.

Pathology prevalence: 183 reviews (60%) had multimodal effect distributions, 153 (50%) showed publication bias signatures, 37 (12%) contained outlier studies, and 23 (7.5%) were methodologically unstable.

The uncertainty decomposition revealed that **heterogeneity accounts for 42.0% of total uncertainty**, sampling variance for 36.5%, publication bias for 19.4%, model choice for 1.1%, and outlier influence for 1.0%. This finding challenges the field's emphasis on estimator choice (DL vs REML vs PM), which addresses only 1% of the uncertainty budget, while the dominant sources — heterogeneity and bias — receive comparatively less methodological attention.

**Conclusions:** The majority of Cochrane meta-analyses contain diagnosable pathologies that substantially alter pooled estimates when corrected. Routine application of automated diagnosis and correction could transform meta-analysis from a single-number summary into a quality-adjusted, uncertainty-decomposed decision tool. The finding that estimator choice contributes only 1% of total uncertainty while heterogeneity and bias contribute 61% argues for a fundamental reorientation of methodological priorities.

---

## What is already known on this topic

- Meta-analyses are vulnerable to multimodality, outliers, publication bias, and methodological instability
- These pathologies are typically assessed individually using separate tools
- No automated system exists to diagnose, correct, and quantify the contribution of each pathology simultaneously

## What this study adds

- Only 19% of Cochrane meta-analyses are pathology-free (Grade A)
- Correcting detected pathologies shifts pooled estimates by a mean of 51%
- Heterogeneity (42%) and publication bias (19%) dominate the uncertainty budget
- Estimator choice (DL vs REML) contributes only 1% of uncertainty — the field's methodological focus is misallocated
- Automated multi-pathology diagnosis and correction is feasible in 37 seconds for 307 reviews

---

## 1. Introduction

Meta-analysis is the cornerstone of evidence-based medicine, yet its outputs — pooled effect estimates with confidence intervals — conceal a complex landscape of potential pathologies. Between-study heterogeneity, publication bias, outlier studies, and sensitivity to methodological choices can each distort the pooled estimate, and their interactions are poorly understood.

Current practice addresses these pathologies individually: I-squared for heterogeneity, Egger's test for publication bias, influence analysis for outliers, and sensitivity analysis for model choice. This piecemeal approach has two fundamental limitations. First, it leaves the synthesis of findings to the analyst's judgement — there is no systematic framework for combining diagnostic results into a single quality assessment. Second, it does not quantify how much each pathology contributes to total uncertainty, making it impossible to prioritise which problems most urgently need addressing.

We developed MetaRepair, an automated pipeline that diagnoses four classes of pathology, applies targeted corrections, decomposes uncertainty into five sources, and produces a decision-grade output for each meta-analysis. We applied MetaRepair to 307 Cochrane systematic reviews to characterise the prevalence and relative importance of each pathology across the evidence base.

## 2. Methods

### 2.1 Diagnostic Modules

**Multimodality:** We fitted 1-component and 2-component Gaussian mixture models via expectation-maximisation (20 iterations) and selected the number of components using the Bayesian Information Criterion (BIC). A BIC improvement > 6 (strong evidence threshold) indicated multimodality.

**Outliers:** Studentized residuals were computed from the DerSimonian-Laird random-effects model. Studies with |residual| > 2.5 were flagged as outliers.

**Publication bias:** Egger's regression test (p < 0.10) and trim-and-fill (L0 estimator, k0 >= 2) were used in combination. Bias was flagged if either test was positive.

**Instability:** Eight specifications (4 estimators: FE, DL, REML, PM; x 2 CI methods: Wald, HKSJ) were computed. A review was flagged as unstable if fewer than 70% of specifications agreed on both the direction and statistical significance of the pooled effect.

### 2.2 Correction Modules

**Robust pooling:** Outlier studies received 10% of their original inverse-variance weight (winsorisation), preserving their contribution while limiting their influence.

**Bias adjustment:** PET-PEESE meta-regression was applied. PET (precision-effect test: effect regressed on SE) was used when the intercept was near zero; PEESE (precision-effect estimate with standard error: effect regressed on SE-squared) otherwise.

**Subgroup splitting:** For multimodal reviews, studies were assigned to clusters via the mixture model posterior responsibilities, and each cluster was pooled separately.

### 2.3 Uncertainty Decomposition

Total uncertainty was decomposed into five additive components:
- **Sampling**: squared standard error of the pooled estimate
- **Heterogeneity**: between-study variance (tau-squared)
- **Model choice**: variance of pooled estimates across estimator specifications
- **Bias**: squared difference between original and bias-corrected estimates
- **Outlier**: squared difference between original and robust estimates

Each component was expressed as a percentage of total uncertainty.

### 2.4 Decision Grading

Reviews were graded A-D based on pathology count and correction magnitude:
- **A**: No pathologies, correction < 5%
- **B**: 1 pathology, correction < 20%
- **C**: 2 pathologies, correction < 50%
- **D**: 3+ pathologies or correction >= 50%

## 3. Results

### 3.1 Grade Distribution

Of 307 reviews: 58 (18.9%) Grade A, 88 (28.7%) Grade B, 79 (25.7%) Grade C, 82 (26.7%) Grade D. Only one in five meta-analyses was pathology-free.

### 3.2 Pathology Prevalence

Multimodality was the most common pathology (183 reviews, 59.6%), followed by publication bias (153, 49.8%), outliers (37, 12.1%), and instability (23, 7.5%). Most reviews had multiple co-occurring pathologies.

### 3.3 Correction Magnitude

The mean correction was 50.7%, indicating that pooled estimates shifted by approximately half their original value when pathologies were addressed. Grade D reviews had mean corrections of 89.2%; Grade A reviews, 2.1%.

### 3.4 Uncertainty Decomposition

Heterogeneity dominated the uncertainty budget at 42.0%, followed by sampling (36.5%), publication bias (19.4%), model choice (1.1%), and outlier influence (1.0%).

## 4. Discussion

### 4.1 The 1% Problem

The meta-analytic literature devotes enormous attention to estimator choice — DL vs REML vs PM vs SJ — yet our uncertainty decomposition shows this contributes only 1.1% of total uncertainty. The dominant sources are heterogeneity (42%) and publication bias (19.4%), which together account for 61% of the uncertainty budget. This suggests a fundamental misallocation of methodological effort: the field is optimising the least important parameter while the most important ones receive less systematic attention.

### 4.2 The Correction Imperative

A mean correction of 51% is alarming. It implies that for every two Cochrane meta-analyses, one would substantially change its pooled estimate if all detectable pathologies were addressed. Current practice reports the uncorrected estimate as the primary finding, with diagnostic results relegated to supplementary tables.

### 4.3 Limitations

The pipeline uses simplified versions of each diagnostic and correction method. The BIC-based multimodality detection may miss subtle bimodality in small samples. PET-PEESE assumes a specific functional form for bias. The uncertainty decomposition is additive, ignoring potential interactions between sources. The 50.7% mean correction may partly reflect over-correction by the PET-PEESE method, which can be conservative.

## 5. Conclusions

Only 19% of Cochrane meta-analyses are free of diagnosable pathologies. Automated diagnosis and correction shifts pooled estimates by a mean of 51%. The uncertainty budget is dominated by heterogeneity (42%) and publication bias (19%), not by estimator choice (1%). These findings argue for a paradigm shift: from single-number meta-analysis to quality-adjusted, uncertainty-decomposed evidence synthesis.

---

## Data Availability

Pipeline code and results: https://github.com/mahmood726-cyber/meta-repair

## Funding
None.

## Competing Interests
The author declares no competing interests.

## References
1. Higgins JPT, et al. Cochrane Handbook. Version 6.4, 2023.
2. Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5:60-78.
3. Viechtbauer W, Cheung MWL. Outlier and influence diagnostics for meta-analysis. Res Synth Methods. 2010;1:112-125.
4. Copas J. What works?: selectivity models and meta-analysis. J R Stat Soc A. 1999;162:95-109.
