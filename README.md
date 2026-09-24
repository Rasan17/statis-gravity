# Statis-Gravity — Clinical Biostatistics & Visualization Suite

**Statis-Gravity** is an interactive, zero-dependency clinical biostatistics and exploratory data analysis web application designed for biomedical researchers, clinicians, and neurosurgeons.

Conceived, supervised design and testing: **Dr G Narenthiran** MB ChB BSc(MedSci)(Hons) MRCS(Ed.) FEBNS FRCS(SN).  
Copyright, G Narenthiran FEBS FRCS(SN), `g_narenthiran@hotmail.com`.  
Dedicated to **Mrs Nirmaladevy Ganesalingam BSc (mother)**.

> ⚠️ **Notice**: AI was used to vibe code this WebApp. The App is still in testing phase. Not to use for clinical, research or decision making.

---

## 🚀 Key Clinical & Statistical Modules

1. **Descriptive Statistics & Normality**
   - Sample size $n$, Mean, Median, Mode, Range, Standard Deviation (SD), Variance, SEM.
   - Quartiles ($Q_1, Q_3$), Interquartile Range (IQR), and Tukey outlier detection fences.
   - 95% and 99% Student's $t$ confidence intervals for the mean.
   - Fisher-Pearson sample skewness, excess kurtosis, and Jarque-Bera normality assessment.
   - Interactive high-DPI frequency histogram with fitted normal density distribution overlay.

2. **Parametric & Non-Parametric Hypothesis Testing**
   - Independent Two-Sample **Welch's $t$-test** (Satterthwaite unequal variances).
   - Independent Two-Sample **Student's $t$-test** (Equal variances).
   - **Paired-Samples $t$-test** for dependent clinical pre/post measurements.
   - **Mann-Whitney $U$ test** (Wilcoxon rank-sum) for non-parametric rank comparisons.
   - Effect size: Cohen's $d$, rank-biserial correlation, and 95% CI of mean difference.
   - Comparative Box-and-Whisker plot with individual jittered patient observations.

3. **Multi-Group Comparisons (ANOVA)**
   - **One-Way Between-Subjects ANOVA**: $F$-statistic, exact $p$-value via incomplete beta / Snedecor $F$.
   - Effect sizes: Eta-squared ($\eta^2$) and Omega-squared ($\omega^2$).
   - Tukey HSD pairwise post-hoc multiple comparison tests.
   - Non-parametric **Kruskal-Wallis $H$ test** with Epsilon-squared ($\epsilon^2$).

4. **2×2 Contingency Tables & Epidemiological Risk Metrics**
   - Pearson's Chi-Square test ($\chi^2$) and Yates' continuity-corrected Chi-Square.
   - **Fisher's Exact Test** with exact two-tailed hypergeometric probabilities.
   - **Odds Ratio (OR)** with Woolf 95% confidence intervals.
   - **Relative Risk (RR)** with log-transformed 95% confidence intervals.
   - **Absolute Risk Reduction (ARR)** and **Number Needed to Treat (NNT)**.

5. **Correlation & Linear Regression**
   - **Pearson's product-moment correlation** ($r$, $R^2$, and $t$-test $p$-value).
   - **Spearman's rank correlation** ($\rho$ and significance).
   - Ordinary Least Squares (OLS) Linear Regression: $y = \beta_1 x + \beta_0$.
   - Standard errors of slope, intercept, and residuals ($s_{y.x}$).
   - Interactive scatter plot with regression trendline and **95% confidence bands**.

6. **Diagnostic Test Accuracy & ROC Curve Analysis**
   - Clinical 2×2 diagnostic matrix: Sensitivity, Specificity, PPV, NPV, Accuracy, Prevalence.
   - Positive Likelihood Ratio ($+LR$), Negative Likelihood Ratio ($-LR$), and **Youden's Index** $J$.
   - Wilson Score 95% confidence intervals for all clinical proportions.
   - Non-parametric **Receiver Operating Characteristic (ROC)** curve generation.
   - Area Under Curve (**AUC**) via trapezoidal integration with Hanley & McNeil standard error.
   - Automated optimal clinical decision threshold / cutoff determination.

7. **Sample Size & Statistical Power Calculation**
   - A-priori sample size calculation for two-group continuous mean differences ($d$, $\alpha$, $1-\beta$).
   - Total sample size ($N$) and required cohort allocation ($n_1, n_2$).
   - Grant and ethics committee-ready protocol justification summary.

---

## 🎨 Visualization & Publication Engine

- **High-DPI Retina Canvas Engine**: Custom, lightweight HTML5 canvas renderer supporting dynamic DPR scaling (`window.devicePixelRatio`).
- **One-Click Export**: Save high-resolution PNG charts directly to your desktop.
- **Publication Ready Reports**: Auto-generated statements formatted according to APA 7th edition and ICMJE medical journal standards, with 1-click clipboard copy.
- **Theme Support**: Seamless toggle between sleek clinical dark mode and publication light mode.

---

## 💻 How to Run Locally

Because Statis-Gravity is crafted with vanilla HTML5, CSS3, and ES6 JavaScript modules, it requires **no npm dependencies or build tools** to run!

```bash
# Simply open in Safari, Chrome, or Firefox:
open index.html

# Or serve via Python:
python3 -m http.server 8080
# then visit http://localhost:8080
```
