Mahmood Ahmad
Tahir Heart Institute
author@example.com

MetaRepair: Automated Diagnosis and Correction of Meta-Analysis Pathologies

Can automated diagnosis and correction of meta-analytic pathologies produce uncertainty estimates decomposing variance into its constituent sources? We applied a four-stage pipeline to 403 Cochrane reviews: multimodality diagnosis via Gaussian mixture BIC comparison, outlier detection through studentized residuals, bias assessment using Egger regression and trim-and-fill, and stability testing across eight specifications combining four estimators with two CI methods. Corrections include Winsorized pooling for outliers, PET-PEESE for bias, and subgroup splitting for multimodal distributions, with five-component decomposition separating sampling, heterogeneity, model, bias, and outlier contributions. Across 403 reviews, the median OR shift after correction was 1.12 with 95% CI 1.05 to 1.21, and 68 percent exhibited at least one pathology. Uncertainty decomposition revealed heterogeneity contributed a median 47 percent of total variance, exceeding sampling uncertainty in most reviews. MetaRepair provides systematic diagnostic-then-correct workflow producing graded estimates with transparent uncertainty attribution. The limitation is that automated correction cannot replace judgment about whether pathologies reflect biological heterogeneity versus artifacts.

Outside Notes

Type: empirical
Primary estimand: Corrected pooled effect with uncertainty decomposition
App: MetaRepair Pipeline v1.0
Data: 403 Cochrane reviews from Pairwise70 dataset
Code: https://github.com/mahmood726-cyber/meta-repair
Version: 1.0
Validation: DRAFT

References

1. Guyatt GH, Oxman AD, Vist GE, et al. GRADE: an emerging consensus on rating quality of evidence and strength of recommendations. BMJ. 2008;336(7650):924-926.
2. Schunemann HJ, Higgins JPT, Vist GE, et al. Completing 'Summary of findings' tables and grading the certainty of the evidence. Cochrane Handbook Chapter 14. Cochrane; 2023.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
