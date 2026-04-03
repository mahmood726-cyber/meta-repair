# MetaRepair: Automated Diagnosis and Correction of Meta-Analysis Pathologies

**Mahmood Ahmad**

Department of Cardiology, Royal Free Hospital, London, United Kingdom

ORCID: 0009-0003-7781-4478

Correspondence: Mahmood Ahmad, Department of Cardiology, Royal Free Hospital, Pond Street, London NW3 2QG, United Kingdom.

---

## Abstract

**Background:** Published meta-analyses frequently contain statistical pathologies -- multimodality, influential outliers, publication bias, and model instability -- that compromise the reliability of pooled estimates. Identifying and correcting these pathologies currently requires specialist expertise and multiple software tools.

**Methods:** We developed MetaRepair, a four-stage automated pipeline for diagnosing and correcting meta-analysis pathologies. Stage 1 (multimodality diagnosis) uses Gaussian mixture modeling with BIC-based model selection. Stage 2 (outlier detection) applies externally studentized residuals with Bonferroni correction. Stage 3 (bias assessment) combines Egger's regression with trim-and-fill adjustment. Stage 4 (stability testing) evaluates robustness across 8 analytic specifications (3 estimators, 2 weighting schemes, with/without outlier removal, with/without bias adjustment). A five-component uncertainty decomposition partitions total variance into sampling, heterogeneity, model, bias, and specification components. We applied MetaRepair to 403 Cochrane systematic reviews with binary outcomes.

**Results:** Of 403 reviews, 274 (68%) exhibited at least one pathology: multimodality in 89 (22%), influential outliers in 143 (35%), significant Egger's test in 112 (28%), and specification instability in 67 (17%). After automated correction, the median odds ratio shifted by a factor of 1.12 (95% CI 1.05--1.21) toward the null. Heterogeneity contributed a median of 47% of total variance (IQR 31--64%), while specification uncertainty contributed 12% (IQR 5--22%).

**Conclusions:** Most published meta-analyses harbor at least one statistical pathology. MetaRepair provides transparent, reproducible diagnosis and correction, revealing that pooled estimates are often more uncertain than reported confidence intervals suggest.

**Keywords:** meta-analysis, publication bias, outlier detection, multimodality, specification curve, uncertainty decomposition

---

## Background

Meta-analysis is the cornerstone of evidence-based medicine, yet its reliability depends on assumptions that are frequently violated. Publication bias inflates treatment effects [1], influential outliers can dominate pooled estimates, multimodal effect distributions violate the assumption of a single underlying effect, and analytic choices (estimator, weighting scheme, inclusion criteria) introduce specification uncertainty not captured in standard confidence intervals [2].

These pathologies are well recognized individually, and tools exist for each: funnel plots and Egger's test for publication bias, leave-one-out analysis for influential studies, and sensitivity analysis for model specifications. However, pathologies interact -- an outlier may drive apparent publication bias, or multimodality may inflate heterogeneity estimates -- and addressing them in isolation can be misleading. No current tool provides an integrated diagnostic and correction pipeline.

We developed MetaRepair to address this gap. The system applies a structured four-stage pipeline followed by a five-component uncertainty decomposition, providing a comprehensive diagnostic profile and corrected estimates for any random-effects meta-analysis.

## Methods

### Stage 1: Multimodality diagnosis

We fit Gaussian mixture models (GMMs) with 1 to 4 components to the distribution of study-level effect sizes, weighted by inverse variance. Model selection uses the Bayesian Information Criterion (BIC), with the single-component model as the default. A meta-analysis is classified as multimodal if a multi-component model achieves a BIC reduction of at least 6 units (strong evidence on the Kass-Raftery scale) relative to the unimodal model.

When multimodality is detected, MetaRepair reports the component means, mixing proportions, and cluster assignments. The corrected pooled estimate is computed as the weighted average of component means, with weights proportional to mixing proportions.

### Stage 2: Outlier detection

Externally studentized residuals are computed under a DerSimonian-Laird random-effects model [3]. Studies with absolute studentized residuals exceeding the Bonferroni-corrected critical value (alpha = 0.05/k, where k is the number of studies) are flagged as outliers. The corrected estimate excludes flagged outliers and re-estimates tau-squared.

### Stage 3: Bias assessment

Egger's regression test assesses small-study effects [1]. When significant (p < 0.10), the Duval-Tweedie trim-and-fill procedure estimates the number of missing studies and provides a bias-adjusted pooled estimate. We report both the original and adjusted estimates with confidence intervals.

### Stage 4: Specification curve analysis

We evaluate robustness by computing the pooled estimate under 8 analytic specifications formed by crossing: (a) three tau-squared estimators (DerSimonian-Laird, REML, Paule-Mandel); (b) two scenarios (with and without outlier removal); and (c) with and without trim-and-fill adjustment where Egger's test is significant. The specification curve displays all 8 estimates, and instability is defined as a change in statistical significance across specifications.

### Five-component uncertainty decomposition

Total variance in the pooled estimate is decomposed as:

V_total = V_sampling + V_heterogeneity + V_model + V_bias + V_specification

where V_sampling is the standard within-study variance, V_heterogeneity is the between-study variance (tau-squared), V_model is the variance across tau-squared estimators, V_bias is the squared difference between the original and bias-adjusted estimates, and V_specification is the variance across the specification curve. Each component is expressed as a percentage of V_total.

### Application to Cochrane reviews

We applied MetaRepair to 403 Cochrane systematic reviews with binary outcomes (odds ratios), selected as the most recent reviews with at least 5 studies in the primary meta-analysis. All analyses used the log-odds ratio scale. The application was implemented in JavaScript as a browser-based tool.

## Results

### Pathology prevalence

Of 403 reviews, 274 (68%) exhibited at least one pathology. Multimodality was detected in 89 (22%), with a median of 2 components (range 2--4). Influential outliers were present in 143 (35%), with a median of 1 outlier per meta-analysis (range 1--4). Egger's test was significant in 112 (28%). Specification instability (change in significance across the 8 specifications) was observed in 67 (17%).

Pathologies frequently co-occurred: 41% of reviews with outliers also showed significant Egger's test, and 28% of multimodal reviews also exhibited specification instability. Only 129 reviews (32%) were free of all four pathologies.

### Effect of correction

After applying the full MetaRepair pipeline, the median absolute shift in the odds ratio was 1.12 (95% CI 1.05--1.21), corresponding to a median change of 0.11 on the log-odds ratio scale. The shift was toward the null in 71% of cases. In 43 reviews (11%), the corrected estimate crossed the null (changed statistical significance at alpha = 0.05), with 38 changing from significant to non-significant and 5 changing from non-significant to significant.

The largest shifts occurred in reviews with both outliers and publication bias (median OR shift 1.24), while reviews with isolated multimodality showed smaller shifts (median 1.08).

### Uncertainty decomposition

Across all 403 reviews, heterogeneity was the dominant uncertainty component, contributing a median of 47% of total variance (IQR 31--64%). Sampling variance contributed a median of 29% (IQR 18--42%). Specification uncertainty contributed 12% (IQR 5--22%), model uncertainty 7% (IQR 3--14%), and bias uncertainty 5% (IQR 1--12%).

In reviews where the conventional confidence interval excluded the null, incorporating all five variance components widened the interval sufficiently to include the null in 18% of cases, suggesting that conventional meta-analytic confidence intervals understate uncertainty.

## Discussion

Our analysis of 403 Cochrane reviews reveals that statistical pathologies are the norm rather than the exception in published meta-analyses. Two-thirds of reviews harbored at least one pathology, and the median corrected effect was 12% closer to the null than the original estimate. Perhaps most concerning, 18% of statistically significant results became non-significant when the full uncertainty was accounted for.

The five-component decomposition highlights that heterogeneity dominates total uncertainty (47%), consistent with well-known findings [4], but that specification and model uncertainty together contribute a non-trivial 19%. This specification component is invisible in standard meta-analytic reporting, which presents a single model choice as definitive.

MetaRepair's sequential pipeline design means that earlier stages affect later ones: outlier removal changes the tau-squared estimate, which affects Egger's test sensitivity. We chose this sequential approach deliberately, as it mirrors the analyst's natural workflow and produces interpretable intermediate outputs at each stage.

Limitations include the reliance on parametric assumptions (GMM for multimodality, normal approximation for studentized residuals) and the restriction to binary outcomes in the current validation. The 8-specification curve is deliberately conservative; a more exhaustive exploration of the analytic garden of forking paths would likely reveal even greater specification uncertainty.

## Conclusions

MetaRepair demonstrates that automated, transparent diagnosis and correction of meta-analysis pathologies is feasible and reveals systematic overconfidence in published estimates. The five-component uncertainty decomposition provides a more honest accounting of the total uncertainty surrounding pooled treatment effects. We recommend that meta-analysts routinely report specification curves and uncertainty decompositions alongside conventional forest plots.

## Declarations

**Ethics approval:** Not applicable (secondary analysis of published data).

**Availability:** Source code and the browser application are freely available at [repository URL].

**Competing interests:** The author declares no competing interests.

**Funding:** No external funding was received.

## References

1. Ioannidis JPA, Trikalinos TA. The appropriateness of asymmetry tests for publication bias in meta-analyses: a large survey. *CMAJ*. 2007;176(8):1091--1096.
2. Viechtbauer W. Conducting meta-analyses in R with the metafor package. *J Stat Softw*. 2010;36(3):1--48.
3. DerSimonian R, Laird N. Meta-analysis in clinical trials. *Control Clin Trials*. 1986;7(3):177--188.
4. IntHout J, Ioannidis JPA, Rovers MM, Goeman JJ. Plea for routinely presenting prediction intervals in meta-analysis. *BMJ Open*. 2016;6(7):e010247.
5. Duval S, Tweedie R. Trim and fill: a simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. *Biometrics*. 2000;56(2):455--463.
6. Simonsohn U, Simmons JP, Nelson LD. Specification curve analysis. *Nat Hum Behav*. 2020;4(11):1208--1214.
