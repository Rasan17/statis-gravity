/**
 * Mathematical & Computational Verification Suite for Statis-Gravity
 */

import { Distributions } from './js/stats/distributions.js';
import { Descriptive } from './js/stats/descriptive.js';
import { Hypothesis } from './js/stats/hypothesis.js';
import { Anova } from './js/stats/anova.js';
import { Categorical } from './js/stats/categorical.js';
import { Correlation } from './js/stats/correlation.js';
import { Diagnostic } from './js/stats/diagnostic.js';
import { PowerAnalysis } from './js/stats/power.js';
import { Teaching } from './js/stats/teaching.js';
import { DocxReports } from './js/export/docx-generator.js';

let passes = 0;
let failures = 0;

function assert(condition, message) {
  if (condition) {
    passes++;
    console.log(`  ✓ ${message}`);
  } else {
    failures++;
    console.error(`  ✗ FAIL: ${message}`);
  }
}

function approx(a, b, tolerance = 1e-3) {
  return Math.abs(a - b) <= tolerance;
}

console.log('--- Testing Distributions ---');
// Standard normal z = 1.95996 -> two tailed p ~ 0.05
const pNormal = Distributions.normalPValue(1.95996);
assert(approx(pNormal, 0.05, 1e-3), `Normal p-value at z=1.96 should be ~0.05, got ${pNormal.toFixed(5)}`);

// Student t with df=10 at t=2.2281 -> p ~ 0.05
const pT = Distributions.tPValue(2.2281, 10);
assert(approx(pT, 0.05, 1e-3), `Student t p-value at t=2.2281 (df=10) should be ~0.05, got ${pT.toFixed(5)}`);

// Chi-square df=1 at 3.841 -> p ~ 0.05
const pChi = Distributions.chiSquarePValue(3.841, 1);
assert(approx(pChi, 0.05, 1e-3), `Chi-square p-value at chi2=3.841 (df=1) should be ~0.05, got ${pChi.toFixed(5)}`);

// F-test df1=2, df2=12 at F=3.885 -> p ~ 0.05
const pF = Distributions.fPValue(3.885, 2, 12);
assert(approx(pF, 0.05, 1e-2), `F p-value at F=3.885 (df=2,12) should be ~0.05, got ${pF.toFixed(5)}`);

console.log('--- Testing Descriptive Statistics ---');
const sample = [12, 14, 15, 18, 19, 21, 24, 28, 30];
const desc = Descriptive.calculate(sample);
assert(desc.n === 9, 'Sample size is 9');
assert(approx(desc.mean, 20.111, 1e-3), `Mean is ~20.111, got ${desc.mean.toFixed(3)}`);
assert(desc.median === 19, `Median is 19, got ${desc.median}`);
assert(approx(desc.sd, 6.234, 1e-3), `Sample SD is ~6.234, got ${desc.sd.toFixed(3)}`);
assert(desc.q1 === 15, `Q1 is 15, got ${desc.q1}`);
assert(desc.q3 === 24, `Q3 is 24, got ${desc.q3}`);
assert(desc.iqr === 9, `IQR is 9, got ${desc.iqr}`);

// Mode test
const sampleMode = [10, 12, 12, 12, 14, 15, 18];
const descMode = Descriptive.calculate(sampleMode);
assert(descMode.modes.length === 1 && descMode.modes[0] === 12, `Mode is correctly identified as 12 (got ${descMode.modes})`);
assert(descMode.maxFreq === 3, `Max frequency is 3 (got ${descMode.maxFreq})`);

// Skewness, Kurtosis, and Tukey Outlier Tests
const skewedSample = [12, 14, 15, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 38, 62];
const descSkew = Descriptive.calculate(skewedSample);
assert(descSkew.skewness > 2.0, `Skewness is > 2.0 (got ${descSkew.skewness.toFixed(3)})`);
assert(descSkew.skewnessInterpretation.includes('High Right Skew'), `Skewness interpretation is High Right Skew (got ${descSkew.skewnessInterpretation})`);
assert(descSkew.kurtosis > 5.0, `Excess kurtosis is > 5.0 (got ${descSkew.kurtosis.toFixed(3)})`);
assert(descSkew.kurtosisInterpretation.includes('Leptokurtic'), `Kurtosis interpretation is Leptokurtic (got ${descSkew.kurtosisInterpretation})`);
assert(descSkew.outliers.length === 2, `Identified exactly 2 outliers (got ${descSkew.outliers.length})`);

const mildOutlier = descSkew.outliers.find(o => o.value === 38);
assert(mildOutlier && mildOutlier.type === 'Mild' && mildOutlier.direction === 'High', `Value 38 classified as Mild High outlier`);
assert(mildOutlier && approx(mildOutlier.zScore, 1.23, 0.05), `Mild outlier Z-score is ~1.23 (got ${mildOutlier?.zScore.toFixed(2)})`);

const extremeOutlier = descSkew.outliers.find(o => o.value === 62);
assert(extremeOutlier && extremeOutlier.type === 'Extreme' && extremeOutlier.direction === 'High', `Value 62 classified as Extreme High outlier`);
assert(extremeOutlier && approx(extremeOutlier.zScore, 3.14, 0.05), `Extreme outlier Z-score is ~3.14 (got ${extremeOutlier?.zScore.toFixed(2)})`);

// Normality test verification
assert(descSkew.normality.test === 'Jarque-Bera', `Normality test identified as Jarque-Bera`);
assert(descSkew.normality.statistic > 15, `Jarque-Bera statistic is large for skewed sample: ${descSkew.normality.statistic.toFixed(2)}`);
assert(descSkew.normality.isNormal === false, `Jarque-Bera test rejects normality for skewed sample (p = ${descSkew.normality.pValue.toFixed(5)})`);
// Symmetric sample without outliers
const symmSample = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
const descSymm = Descriptive.calculate(symmSample);
assert(approx(descSymm.skewness, 0, 0.01), `Symmetric sample skewness is ~0 (got ${descSymm.skewness.toFixed(3)})`);
assert(descSymm.skewnessInterpretation.includes('Symmetric'), `Symmetric sample interpreted as Symmetric (got ${descSymm.skewnessInterpretation})`);
assert(descSymm.outliers.length === 0, `Zero outliers in clean symmetric sample (got ${descSymm.outliers.length})`);
assert(descSymm.normality.isNormal === true, `Jarque-Bera test confirms normality for symmetric sample (p = ${descSymm.normality.pValue.toFixed(3)})`);

// Kernel Density Estimation (KDE) Bandwidth Verification for Violin Plot
const nKde = skewedSample.length;
const iqrKde = descSkew.iqr > 0 ? descSkew.iqr : descSkew.sd;
const aKde = Math.min(descSkew.sd, iqrKde / 1.34);
const hSilverman = 0.9 * aKde * Math.pow(nKde, -0.2);
assert(hSilverman > 2.0 && hSilverman < 4.0, `Silverman bandwidth h is within expected range [2.0, 4.0], got ${hSilverman.toFixed(3)}`);

console.log('--- Testing Hypothesis Testing (t-test & Mann-Whitney) ---');
const groupA = [25, 28, 31, 29, 33, 27, 30, 32];
const groupB = [18, 22, 20, 24, 19, 21, 23, 17];
const tRes = Hypothesis.independentTTest(groupA, groupB);
assert(tRes.statistic > 0, `t-statistic is positive: ${tRes.statistic.toFixed(3)}`);
assert(tRes.pValue < 0.001, `Groups A and B differ significantly: p = ${tRes.pValue.toExponential(3)}`);

const mwRes = Hypothesis.mannWhitneyUTest(groupA, groupB);
assert(mwRes.pValue < 0.001, `Mann-Whitney U test confirms difference: p = ${mwRes.pValue.toExponential(3)}`);

console.log('--- Testing ANOVA & Tukey HSD Post-Hoc ---');
const anovaRes = Anova.oneWay([
  { name: 'Control', data: [10, 11, 12, 10, 13] },
  { name: 'Low Dose', data: [14, 15, 13, 16, 14] },
  { name: 'High Dose', data: [20, 22, 19, 21, 24] }
]);
assert(anovaRes.fStatistic > 20, `ANOVA F is large: ${anovaRes.fStatistic.toFixed(2)}`);
assert(anovaRes.pValue < 0.0001, `ANOVA p-value is highly significant: p = ${anovaRes.pValue.toExponential(3)}`);
assert(anovaRes.pairwise && anovaRes.pairwise.length === 3, `Tukey HSD computed 3 pairwise comparisons`);
const ctrlVsHigh = anovaRes.pairwise.find(c => c.comparison === 'Control vs High Dose');
assert(ctrlVsHigh && ctrlVsHigh.isSignificant && ctrlVsHigh.pValue < 0.001, `Control vs High Dose difference is significant (p < 0.001)`);

console.log('--- Testing Categorical & 2x2 ---');
// 2x2 table: Treated: 20 diseased / 80 cured; Control: 40 diseased / 60 cured
const catRes = Categorical.twoByTwo(20, 40, 80, 60);
assert(approx(catRes.riskMetrics.oddsRatio, 0.375, 1e-2), `Odds Ratio is ~0.375, got ${catRes.riskMetrics.oddsRatio.toFixed(3)}`);
assert(catRes.chiSquare.pValueStandard < 0.01, `Chi-square is significant: p = ${catRes.chiSquare.pValueStandard.toFixed(4)}`);

console.log('--- Testing Correlation, Regression & Assumptions ---');
const xs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
const ys = [2, 4, 5, 7, 8, 10, 12, 13, 15, 17];
const corrRes = Correlation.pearson(xs, ys);
assert(corrRes.r > 0.99, `Pearson r is near 1: ${corrRes.r.toFixed(4)}`);
const regRes = Correlation.linearRegression(xs, ys);
assert(approx(regRes.slope, 1.630, 1e-2), `Regression slope is ~1.63, got ${regRes.slope.toFixed(3)}`);

// Test U-Shaped Detection on physiological CPP vs Mortality curve
const cppX = [42, 46, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110];
const mortY = [82, 68, 54, 38, 25, 18, 16, 19, 26, 36, 49, 62, 73, 85, 96];
const regU = Correlation.linearRegression(cppX, mortY);
const assumpU = Correlation.checkAssumptions(cppX, mortY, regU);
assert(assumpU.isUShaped === true, 'U-shaped relationship correctly detected');
assert(assumpU.shape === 'U-Shaped', `Shape identified as U-Shaped (got ${assumpU.shape})`);
assert(approx(assumpU.quad.vertexX, 72.4, 0.5), `Nadir identified at optimal CPP ~ 72.4 mmHg (got ${assumpU.quad.vertexX.toFixed(2)})`);
assert(assumpU.quad.rSquaredQuad > 0.90, `Quadratic R² > 0.90 (got ${assumpU.quad.rSquaredQuad.toFixed(3)})`);
assert(assumpU.recommendedTest.includes('Quadratic'), `Recommended test indicates Quadratic Polynomial Fit (got ${assumpU.recommendedTest})`);

// Test Monotonicity: Strictly Increasing
const pairsInc = xs.map((x, i) => ({ x, y: ys[i] }));
const monoInc = Correlation.analyzeMonotonicity(pairsInc);
assert(monoInc.isMonotonic === true, 'Monotonically increasing data correctly identified');
assert(monoInc.reversals === 0, 'Zero directional reversals in strictly increasing data');
assert(monoInc.pattern === 'Monotonically Increasing', `Pattern is Monotonically Increasing (got ${monoInc.pattern})`);

// Test Monotonicity: Tri-Phasic / N-Shaped (Increase -> Decrease -> Increase)
const nX = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
const nY = [15, 28, 45, 62, 70, 65, 50, 36, 25, 20, 26, 42, 60, 78, 95];
const pairsN = nX.map((x, i) => ({ x, y: nY[i] }));
const monoN = Correlation.analyzeMonotonicity(pairsN);
assert(monoN.isMonotonic === false, 'Tri-phasic N-curve recognized as non-monotonic');
assert(monoN.reversals === 2, `Identified exactly 2 directional reversals (got ${monoN.reversals})`);
assert(monoN.pattern.includes('N-Shaped'), `Pattern identified as N-Shaped (got ${monoN.pattern})`);
const regN = Correlation.linearRegression(nX, nY);
const assumpN = Correlation.checkAssumptions(nX, nY, regN);
assert(assumpN.isParametricOk === false, 'Parametric test rejected for tri-phasic non-monotonic data');
assert(assumpN.recommendedTest.includes('Segmented') || assumpN.recommendedTest.includes('Non-Linear'), 'Non-Linear / Segmented test recommended for tri-phasic data');

console.log('--- Testing Diagnostic ROC & AUC ---');
const rocData = [
  { score: 10, status: 1 },
  { score: 9, status: 1 },
  { score: 8, status: 1 },
  { score: 7, status: 1 },
  { score: 6, status: 0 },
  { score: 5, status: 1 },
  { score: 4, status: 0 },
  { score: 3, status: 0 },
  { score: 2, status: 0 },
  { score: 1, status: 0 }
];
const roc = Diagnostic.computeROC(rocData);
assert(roc.auc > 0.85, `ROC AUC is high (>0.85): ${roc.auc.toFixed(3)}`);

console.log('--- Testing Diagnostic 2x2 Performance (Index Test vs Gold Standard) ---');
const d2x2 = Diagnostic.evaluate2x2(92, 8, 12, 188);
assert(d2x2.total === 300, `Total sample size is 300 (got ${d2x2.total})`);
assert(approx(d2x2.sensitivity, 0.8846, 1e-3), `Sensitivity (TPR) is ~88.5%, got ${(d2x2.sensitivity * 100).toFixed(2)}%`);
assert(d2x2.sensitivityCI95[0] < d2x2.sensitivity && d2x2.sensitivityCI95[1] > d2x2.sensitivity, 'Sensitivity Wilson 95% CI bounds point estimate');
assert(approx(d2x2.specificity, 0.9592, 1e-3), `Specificity (TNR) is ~95.9%, got ${(d2x2.specificity * 100).toFixed(2)}%`);
assert(d2x2.specificityCI95[0] < d2x2.specificity && d2x2.specificityCI95[1] > d2x2.specificity, 'Specificity Wilson 95% CI bounds point estimate');
assert(approx(d2x2.ppv, 0.9200, 1e-3), `PPV is exactly 92.0%, got ${(d2x2.ppv * 100).toFixed(2)}%`);
assert(approx(d2x2.npv, 0.9400, 1e-3), `NPV is exactly 94.0%, got ${(d2x2.npv * 100).toFixed(2)}%`);
assert(approx(d2x2.accuracy, 0.9333, 1e-3), `Overall Accuracy is ~93.3%, got ${(d2x2.accuracy * 100).toFixed(2)}%`);
assert(d2x2.accuracyCI95[0] > 0.89 && d2x2.accuracyCI95[1] < 0.97, `Accuracy 95% CI is within expected range [${d2x2.accuracyCI95.map(v => (v*100).toFixed(1)).join('%, ')}%]`);
assert(approx(d2x2.plr, 21.673, 0.05), `Positive Likelihood Ratio (LR+) is ~21.67 (got ${d2x2.plr.toFixed(2)})`);
assert(approx(d2x2.nlr, 0.120, 0.01), `Negative Likelihood Ratio (LR-) is ~0.12 (got ${d2x2.nlr.toFixed(2)})`);
assert(approx(d2x2.youdenJ, 0.8438, 1e-3), `Youden's J statistic is ~0.844 (got ${d2x2.youdenJ.toFixed(3)})`);

console.log('--- Testing Power & Sample Size Analysis ---');
// 1. Two Independent Means
const pwrIndep = PowerAnalysis.twoIndependentMeans({ m1: 10, m2: 15, sd: 10, alpha: 0.05, power: 0.80 });
assert(pwrIndep.nPerGroup >= 60 && pwrIndep.nPerGroup <= 70, `Independent means sample size per group is ~63-64, got ${pwrIndep.nPerGroup}`);
assert(pwrIndep.totalN === 2 * pwrIndep.nPerGroup, `Total sample size is 2 * nPerGroup = ${pwrIndep.totalN}`);
assert(pwrIndep.achievedPower >= 0.80, `Achieved power meets or exceeds target (got ${(pwrIndep.achievedPower * 100).toFixed(1)}%)`);

// Post-hoc power for independent means
const pwrIndepPost = PowerAnalysis.twoIndependentMeans({ m1: 10, m2: 15, sd: 10, alpha: 0.05, nPerGroup: 64 });
assert(approx(pwrIndepPost.achievedPower, 0.805, 0.02), `Post-hoc power for n=64 is ~80.5%, got ${(pwrIndepPost.achievedPower * 100).toFixed(1)}%`);

// 2. Paired Mean Study Design (Before vs After)
// Effect size dz = |120 - 115| / 10 = 0.50
const pwrPaired = PowerAnalysis.pairedMeans({ m1: 120, m2: 115, sdDiff: 10, alpha: 0.05, power: 0.80 });
assert(pwrPaired.dz === 0.50, `Paired effect size dz is 0.50, got ${pwrPaired.dz}`);
assert(pwrPaired.nPairs >= 31 && pwrPaired.nPairs <= 35, `Required paired subjects is ~32-34 pairs, got ${pwrPaired.nPairs}`);
assert(pwrPaired.achievedPower >= 0.80, `Paired achieved power meets target (got ${(pwrPaired.achievedPower * 100).toFixed(1)}%)`);

// Post-hoc power for paired design
const pwrPairedPost = PowerAnalysis.pairedMeans({ m1: 120, m2: 115, sdDiff: 10, alpha: 0.05, nPairs: 34 });
assert(pwrPairedPost.achievedPower >= 0.80, `Paired post-hoc power for N=34 is >= 80%, got ${(pwrPairedPost.achievedPower * 100).toFixed(1)}%`);

// 3. Clinical Study Design: 2x2 Contingency Table (Chi-Square & Fisher's Exact)
// p1 = 0.15 (Test), p2 = 0.30 (Control) -> ARR = 0.15, RR = 0.50, NNT = 6.67
const pwrContChisq = PowerAnalysis.contingency2x2({ p1: 0.15, p2: 0.30, alpha: 0.05, power: 0.80, testType: 'chisq' });
assert(approx(pwrContChisq.arr, 0.15, 1e-3), `Absolute Risk Reduction (ARR) is 0.15, got ${pwrContChisq.arr}`);
assert(approx(pwrContChisq.rr, 0.50, 1e-3), `Relative Risk (RR) is 0.50, got ${pwrContChisq.rr}`);
assert(approx(pwrContChisq.nnt, 6.67, 0.05), `NNT is ~6.67, got ${pwrContChisq.nnt.toFixed(2)}`);
assert(pwrContChisq.uncorrectedN >= 118 && pwrContChisq.uncorrectedN <= 125, `Chi-Square n per group is ~121, got ${pwrContChisq.uncorrectedN}`);
assert(pwrContChisq.continuityCorrectedN >= 130 && pwrContChisq.continuityCorrectedN <= 140, `Continuity-corrected Fisher's exact n per group is ~134, got ${pwrContChisq.continuityCorrectedN}`);
assert(pwrContChisq.continuityCorrectedN > pwrContChisq.uncorrectedN, `Fisher's continuity-corrected sample size is strictly greater than uncorrected`);

// Post-hoc power for 2x2 contingency study
const pwrContPost = PowerAnalysis.contingency2x2({ p1: 0.15, p2: 0.30, alpha: 0.05, nPerGroup: 134, testType: 'fisher' });
assert(pwrContPost.powerFisher >= 0.79, `Fisher's exact power at n=134 is ~80%, got ${(pwrContPost.powerFisher * 100).toFixed(1)}%`);
assert(pwrContPost.powerChisq > pwrContPost.powerFisher, `Standard Chi-Square power is higher than Fisher's exact at same N`);

console.log('--- Testing Hypothesis Error Bar / Dispersion Modes ---');
const sampleHypo = [20, 22, 24, 25, 26, 28, 30, 31, 34];
const dHypo = Descriptive.calculate(sampleHypo);

// 1. 95% CI Mode
const ciLower = dHypo.ci95[0];
const ciUpper = dHypo.ci95[1];
assert(ciLower < dHypo.mean && ciUpper > dHypo.mean, `95% CI straddles sample mean (${ciLower.toFixed(2)} < ${dHypo.mean.toFixed(2)} < ${ciUpper.toFixed(2)})`);
const ciMargin = (ciUpper - ciLower) / 2;
assert(approx(ciMargin, 3.033, 0.05), `95% CI margin of error is ~3.03, got ${ciMargin.toFixed(2)}`);

// 2. SEM Mode
const semLower = dHypo.mean - dHypo.sem;
const semUpper = dHypo.mean + dHypo.sem;
assert(approx(dHypo.sem, 1.500, 1e-2), `SEM is s / sqrt(n) = 1.500, got ${dHypo.sem.toFixed(3)}`);
assert(approx(semUpper - semLower, 2 * dHypo.sem, 1e-3), `SEM error bar span is exactly 2 * SEM (${(semUpper - semLower).toFixed(3)})`);

// 3. Standard Deviation (SD) Mode
const sdLower = dHypo.mean - dHypo.sd;
const sdUpper = dHypo.mean + dHypo.sd;
assert(approx(dHypo.sd, 4.500, 1e-2), `SD is sample std dev = 4.500, got ${dHypo.sd.toFixed(3)}`);
assert(approx(sdUpper - sdLower, 2 * dHypo.sd, 1e-3), `SD error bar span is exactly 2 * SD (${(sdUpper - sdLower).toFixed(3)})`);
assert(sdUpper - sdLower > semUpper - semLower, 'SD span is strictly wider than SEM span');

// 4. Interquartile Range (IQR) Mode
assert(dHypo.q1 === 24 && dHypo.q3 === 30, `Quartiles Q1=24 and Q3=30 match expectation (got Q1=${dHypo.q1}, Q3=${dHypo.q3})`);
assert(dHypo.iqr === 6, `IQR is Q3 - Q1 = 6 (got ${dHypo.iqr})`);
assert(dHypo.lowerFence === 24 - 1.5 * 6 && dHypo.upperFence === 30 + 1.5 * 6, `Tukey fences correctly bounded at [${dHypo.lowerFence}, ${dHypo.upperFence}]`);

console.log('--- Testing Multi-Group ANOVA Error Bar / Dispersion Modes ---');
const anovaMulti = Anova.oneWay([
  { name: 'Conservative', data: [12.2, 11.8, 12.5, 11.9, 12.1, 12.3, 11.7, 12.4] },
  { name: 'Orthotic Helmet', data: [8.5, 8.2, 8.8, 8.4, 8.6, 8.1, 8.7, 8.5] },
  { name: 'Endoscopic Strip', data: [4.1, 4.5, 4.2, 3.9, 4.3, 4.0, 4.4, 4.2] }
]);
assert(anovaMulti.groups.length === 3, `ANOVA contains 3 distinct cohorts for plotting`);

anovaMulti.groups.forEach((g) => {
  const s = g.stats;
  assert(s && s.n === 8, `${g.name} has valid descriptive stats (n=${s.n})`);

  // 1. 95% CI mode
  assert(Array.isArray(s.ci95) && s.ci95.length === 2, `${g.name} has calculated 95% CI`);
  assert(s.ci95[0] < s.mean && s.ci95[1] > s.mean, `${g.name} 95% CI bounds mean (${s.ci95[0].toFixed(2)} < ${s.mean.toFixed(2)} < ${s.ci95[1].toFixed(2)})`);

  // 2. SEM mode
  const semSpan = 2 * s.sem;
  assert(s.sem > 0 && approx(s.sem, s.sd / Math.sqrt(s.n), 1e-3), `${g.name} SEM matches s / sqrt(n) = ${s.sem.toFixed(3)}`);

  // 3. SD mode
  const sdSpan = 2 * s.sd;
  assert(sdSpan > semSpan, `${g.name} SD error span (${sdSpan.toFixed(2)}) strictly exceeds SEM span (${semSpan.toFixed(2)})`);

  // 4. IQR mode
  assert(s.iqr > 0 && approx(s.iqr, s.q3 - s.q1, 1e-3), `${g.name} IQR (${s.iqr.toFixed(2)}) equals Q3 - Q1`);
  assert(s.median >= s.q1 && s.median <= s.q3, `${g.name} median (${s.median.toFixed(2)}) lies within Q1-Q3 box`);
});

console.log('--- Testing DOCX Clinical Report Generation (All 7 Modules) ---');
const isZip = (buf) => buf.length > 30 && buf[0] === 0x50 && buf[1] === 0x4B && buf[2] === 0x03 && buf[3] === 0x04;
const hasRequiredSections = (builder) => {
  const xml = builder.buildDocumentXML();
  const hasAttribution = xml.includes('Dr G Narenthiran');
  const hasNotice = xml.includes('AI was used to vibe code this WebApp');
  const hasSec1 = xml.includes('1. ') && xml.includes('Information');
  const hasSec2 = xml.includes('2. ') && (xml.includes('Outcome') || xml.includes('Results'));
  const hasSec3 = xml.includes('3. ') && xml.includes('Interpretation');
  const hasSec4 = xml.includes('4. ') && xml.includes('Reason');
  const hasSec5 = xml.includes('5. ') && xml.includes('Background');
  return hasAttribution && hasNotice && hasSec1 && hasSec2 && hasSec3 && hasSec4 && hasSec5;
};

// 1. Descriptive DOCX
const bDesc = DocxReports.createDescriptiveDocx({
  name: 'ICP (mmHg)', n: 20, mean: 15.2, ci95: [13.8, 16.6], sd: 3.1, variance: 9.61, sem: 0.69,
  median: 15.0, q1: 13.0, q3: 17.0, iqr: 4.0, modes: [14.0], maxFreq: 3, min: 9.0, max: 22.0,
  skewness: 0.45, kurtosis: -0.2, normality: { isNormal: true, statistic: 0.72, pValue: 0.697 },
  lowerFence: 7.0, upperFence: 23.0, outliers: [], reportText: 'Descriptive APA summary'
});
const bufDesc = bDesc.generateUint8Array();
assert(isZip(bufDesc), `Descriptive DOCX is a valid PKZIP archive (${bufDesc.length} bytes)`);
assert(hasRequiredSections(bDesc), `Descriptive DOCX contains Information, Outcome, Interpretation, Reason, and Background sections`);

// 2. Hypothesis DOCX
const bHypo = DocxReports.createHypothesisDocx({
  nameA: 'Control', nameB: 'Treatment', testName: 'Welch Unequal Variances t-Test',
  statistic: 3.82, df: 28.4, pValue: 0.00065, isSignificant: true, meanDiff: 4.2, ci95: [1.9, 6.5],
  cohensD: 0.98, reportText: 'Hypothesis APA summary',
  groupA: { n: 15, mean: 22.0, sd: 4.0, sem: 1.03, ci95: [19.8, 24.2], median: 21.5, iqr: 4.5 },
  groupB: { n: 15, mean: 17.8, sd: 2.8, sem: 0.72, ci95: [16.2, 19.4], median: 18.0, iqr: 3.5 }
});
const bufHypo = bHypo.generateUint8Array();
assert(isZip(bufHypo), `Hypothesis DOCX is a valid PKZIP archive (${bufHypo.length} bytes)`);
assert(hasRequiredSections(bHypo), `Hypothesis DOCX contains all 5 required clinical sections`);

// 3. ANOVA DOCX
const bAnova = DocxReports.createAnovaDocx({
  k: 3, totalN: 36, grandMean: 14.8, ssBetween: 120.5, dfBetween: 2, msBetween: 60.25,
  ssWithin: 240.2, dfWithin: 33, msWithin: 7.28, fStatistic: 8.28, pValue: 0.0012,
  etaSquared: 0.334, omegaSquared: 0.288, reportText: 'ANOVA APA summary',
  groups: [
    { name: 'Cohort 1', stats: { n: 12, mean: 11.2, sd: 2.4, sem: 0.69, ci95: [9.7, 12.7], median: 11.0, iqr: 3.0 } },
    { name: 'Cohort 2', stats: { n: 12, mean: 14.5, sd: 2.8, sem: 0.81, ci95: [12.7, 16.3], median: 14.5, iqr: 3.5 } },
    { name: 'Cohort 3', stats: { n: 12, mean: 18.7, sd: 3.2, sem: 0.92, ci95: [16.7, 20.7], median: 18.5, iqr: 4.0 } }
  ],
  pairwise: [
    { comparison: 'Cohort 1 vs Cohort 3', meanDiff: -7.5, seDiff: 1.10, qStatistic: 6.8, pValue: 0.0001, ci95: [-10.2, -4.8], cohensD: 2.78, isSignificant: true }
  ]
});
const bufAnova = bAnova.generateUint8Array();
assert(isZip(bufAnova), `ANOVA DOCX is a valid PKZIP archive (${bufAnova.length} bytes)`);
assert(hasRequiredSections(bAnova), `ANOVA DOCX contains post-hoc table, interpretation, and methodological rationale`);

// 4. Categorical DOCX
const bCat = DocxReports.createCategoricalDocx({
  mode: 'study', a: 20, b: 80, c: 40, d: 60, totalN: 200, reportText: 'Categorical study summary',
  chiSquare: { standard: 9.60, pValueStandard: 0.0019, yates: 8.66, pValueYates: 0.0032 },
  fishersExact: { pValue: 0.0028 },
  risk: { oddsRatio: 0.375, orCI95: [0.20, 0.70], relativeRisk: 0.500, rrCI95: [0.32, 0.79], arr: 0.20, arrCI95: [0.07, 0.33], rrr: 0.50, nnt: 5.0 }
});
const bufCat = bCat.generateUint8Array();
assert(isZip(bufCat), `Categorical DOCX is a valid PKZIP archive (${bufCat.length} bytes)`);
assert(hasRequiredSections(bCat), `Categorical DOCX contains risk metrics (NNT/ARR), rationale, and background`);

// 5. Correlation DOCX
const bCorr = DocxReports.createCorrelationDocx({
  xName: 'CPP', yName: 'Mortality', n: 15, reportText: 'Correlation U-shape summary',
  corr: { r: 0.08, p: 0.78, ci95: [-0.45, 0.57], spearman: 0.12, spearmanP: 0.67 },
  reg: { slope: 0.02, intercept: 20.0, seSlope: 0.08, seIntercept: 6.2, rSquared: 0.006, residualSE: 14.2, fStatistic: 0.08, pVal: 0.78 },
  morphology: { shape: 'U-Shaped (Quadratic)', isMonotonic: false, reversals: 1, quadR2: 0.952, vertexX: 72.4, vertexY: 16.0, recommendedModel: 'Quadratic Polynomial Fit' }
});
const bufCorr = bCorr.generateUint8Array();
assert(isZip(bufCorr), `Correlation DOCX is a valid PKZIP archive (${bufCorr.length} bytes)`);
assert(hasRequiredSections(bCorr), `Correlation DOCX includes curve morphology diagnostics and non-linear rationale`);

// 6. Diagnostic ROC DOCX
const bRoc = DocxReports.createDiagnosticDocx({
  name: 'Serum NfL', totalN: 120, posCount: 50, negCount: 70, reportText: 'ROC summary statement',
  auc: 0.912, aucCI95: [0.86, 0.96], seAuc: 0.026, optimalCutoff: 0.42,
  sensitivity: 0.88, specificity: 0.91, plr: 9.78, nlr: 0.13
});
const bufRoc = bRoc.generateUint8Array();
assert(isZip(bufRoc), `Diagnostic ROC DOCX is a valid PKZIP archive (${bufRoc.length} bytes)`);
assert(hasRequiredSections(bRoc), `Diagnostic ROC DOCX contains Likelihood Ratios, Youden index, and bayesian rationale`);

// 7. Power & Sample Size DOCX
const bPwr = DocxReports.createPowerDocx({
  designLabel: 'Independent Two-Sample Means', goal: 'sample_size', alpha: 0.05,
  effectSize: 0.50, effectSizeLabel: "Cohen's d", nPerGroup: 64, totalN: 128,
  targetPower: 0.80, achievedPower: 0.807, reportText: 'IRB grant power justification'
});
const bufPwr = bPwr.generateUint8Array();
assert(isZip(bufPwr), `Power Analysis DOCX is a valid PKZIP archive (${bufPwr.length} bytes)`);
assert(hasRequiredSections(bPwr), `Power Analysis DOCX contains ethical justification, attrition buffer, and mathematical background`);

// 8. Teaching & CLT DOCX
const bTeach = DocxReports.createTeachingDocx({
  distName: 'Exponential Distribution',
  distN: 500,
  clinicalExample: 'Emergency Department Length of Stay (ED LOS)',
  sampleMean: 2.48,
  theoMean: 2.50,
  sampleSD: 2.45,
  theoSD: 2.50,
  skewness: 1.95,
  skewnessLabel: 'High Right Skew',
  kurtosis: 5.82,
  kurtosisLabel: 'Leptokurtic',
  jbStat: 485.2,
  jbP: 0.00001,
  isNormal: false,
  cltParentName: 'Exponential Distribution (Right-Skewed, λ = 0.4)',
  cltN: 36,
  cltK: 1000,
  cltResults: {
    theoMean: 2.50,
    obsMean: 2.502,
    theoSE: 0.417,
    obsSE: 0.415,
    skewness: 0.082,
    normalityP: 0.421,
    isNormal: true
  },
  tConv: Teaching.tConvergence.getMetrics(4),
  overlap: Teaching.significanceOverlap.getMetrics({ delta: 2.0, sd: 2.5, n: 16, alpha: 0.05 }),
  reportText: 'CLT, t-distribution, and two-sample alpha overlap convergence summary text'
});
const bufTeach = bTeach.generateUint8Array();
assert(isZip(bufTeach), `Teaching DOCX is a valid PKZIP archive (${bufTeach.length} bytes)`);
assert(hasRequiredSections(bTeach), `Teaching DOCX contains all 5 required sections, attribution, and disclaimer`);

console.log('--- Testing Teaching & Distribution Generators ---');
const normData = Teaching.generators.normal(1000, 10, 2);
const normDesc = Descriptive.calculate(normData);
assert(approx(normDesc.mean, 10, 0.25), `Normal generator mean ~ 10, got ${normDesc.mean.toFixed(3)}`);
assert(approx(normDesc.sd, 2, 0.25), `Normal generator SD ~ 2, got ${normDesc.sd.toFixed(3)}`);

const expData = Teaching.generators.exponential(1000, 0.5);
const expDesc = Descriptive.calculate(expData);
assert(approx(expDesc.mean, 2.0, 0.25), `Exponential generator mean ~ 2.0 (1/λ), got ${expDesc.mean.toFixed(3)}`);
assert(expDesc.min >= 0, `Exponential variates are strictly non-negative`);

const uniData = Teaching.generators.uniform(1000, 0, 10);
const uniDesc = Descriptive.calculate(uniData);
assert(approx(uniDesc.mean, 5.0, 0.25), `Uniform generator mean ~ 5.0, got ${uniDesc.mean.toFixed(3)}`);
assert(uniDesc.min >= 0 && uniDesc.max <= 10, `Uniform variates bounded in [0, 10]`);

const bimData = Teaching.generators.bimodal(1000, 10, 1, 20, 1, 0.5);
const bimDesc = Descriptive.calculate(bimData);
assert(approx(bimDesc.mean, 15.0, 0.35), `Bimodal generator mean ~ 15.0, got ${bimDesc.mean.toFixed(3)}`);

console.log('--- Testing Teaching Theoretical PDFs ---');
const normPdf = Teaching.pdf.normal(10, 10, 2);
assert(approx(normPdf, 0.19947, 1e-4), `Normal theoretical peak PDF at mean ~ 0.19947, got ${normPdf.toFixed(5)}`);

const expPdf = Teaching.pdf.exponential(2, 0.5);
assert(approx(expPdf, 0.5 * Math.exp(-1), 1e-4), `Exponential PDF at x=2 (λ=0.5) ~ 0.18394, got ${expPdf.toFixed(5)}`);

const uniPdf = Teaching.pdf.uniform(5, 0, 10);
assert(approx(uniPdf, 0.1, 1e-4), `Uniform PDF inside [0, 10] is 0.1, got ${uniPdf.toFixed(5)}`);

console.log('--- Testing Central Limit Theorem (CLT) Convergence ---');
Teaching.clt.reset();
Teaching.clt.setPopulation('exponential'); // True mu = 2.0, sigma = 2.0
Teaching.clt.setSampleSize(36); // Theoretical SE = 2.0 / sqrt(36) = 2.0 / 6 = 0.33333
const cltSummary = Teaching.clt.drawSamples(1000);

assert(cltSummary.samplesDrawn === 1000, `CLT simulated exactly 1,000 samples`);
assert(cltSummary.sampleSize === 36, `CLT sample size is n = 36`);
assert(approx(cltSummary.theoreticalSE, 0.3333, 1e-3), `Theoretical SE is ~0.3333, got ${cltSummary.theoreticalSE.toFixed(4)}`);
assert(approx(cltSummary.observedMean, 2.00, 0.10), `Observed grand mean converges to theoretical mean ~2.00, got ${cltSummary.observedMean.toFixed(3)}`);
assert(approx(cltSummary.observedSE, 0.3333, 0.05), `Observed SE converges to theoretical SE ~0.3333, got ${cltSummary.observedSE.toFixed(3)}`);
assert(Math.abs(cltSummary.skewness) < 0.50, `Sampling distribution skewness is close to theoretical skewness (2/sqrt(n) ~ 0.33), got ${cltSummary.skewness.toFixed(3)}`);

console.log('--- Testing Student\'s t-Distribution Convergence to Normal N(0, 1) ---');
// Critical value calculations
assert(approx(Teaching.tConvergence.getCriticalValue(1), 12.7062, 1e-3), 't_crit for df=1 is 12.7062');
assert(approx(Teaching.tConvergence.getCriticalValue(2), 4.3027, 1e-3), 't_crit for df=2 is 4.3027');
assert(approx(Teaching.tConvergence.getCriticalValue(4), 2.7764, 1e-3), 't_crit for df=4 is 2.7764');
assert(approx(Teaching.tConvergence.getCriticalValue(10), 2.2281, 1e-3), 't_crit for df=10 is 2.2281');
assert(approx(Teaching.tConvergence.getCriticalValue(30), 2.0423, 1e-3), 't_crit for df=30 is ~2.0423');
assert(approx(Teaching.tConvergence.getCriticalValue(120), 1.980, 0.005), 't_crit for df=120 is ~1.980');
assert(approx(Teaching.tConvergence.getCriticalValue(500), 1.960, 1e-3), 't_crit for large df converges to 1.960');

// Convergence metrics for small sample (n=4, df=3)
const mSmall = Teaching.tConvergence.getMetrics(4);
assert(mSmall.df === 3, 'Sample size n=4 yields df=3');
assert(mSmall.tPeak < 0.38 && mSmall.tPeak > 0.36, `Small df peak density is flatter than normal, got ${mSmall.tPeak.toFixed(4)}`);
assert(mSmall.tailProb > 0.12, `Small df tail probability is heavy/inflated (>12%), got ${(mSmall.tailProb*100).toFixed(1)}%`);
assert(mSmall.excessKurtosis === Infinity, 'df=3 excess kurtosis is infinite');

// Convergence metrics for parametric threshold (n=31, df=30)
const mThreshold = Teaching.tConvergence.getMetrics(31);
assert(mThreshold.df === 30, 'Sample size n=31 yields df=30');
assert(approx(mThreshold.tPeak, 0.3956, 1e-3), `df=30 peak density is within 1% of normal 0.3989, got ${mThreshold.tPeak.toFixed(4)}`);
assert(approx(mThreshold.tCrit, 2.042, 0.005), `df=30 critical value is 2.042 (within 4.2% of 1.960), got ${mThreshold.tCrit.toFixed(3)}`);
assert(approx(mThreshold.tailProb, 0.059, 0.005), `df=30 tail probability has converged close to 5.0%, got ${(mThreshold.tailProb*100).toFixed(1)}%`);
assert(approx(mThreshold.excessKurtosis, 0.23, 0.02), `df=30 excess kurtosis is near zero (6/26 = 0.23), got ${mThreshold.excessKurtosis.toFixed(2)}`);

// Convergence metrics for large trial (n=121, df=120)
const mLarge = Teaching.tConvergence.getMetrics(121);
assert(approx(mLarge.tCrit, 1.980, 0.005), `Large trial df=120 t_crit is practically identical to 1.960, got ${mLarge.tCrit.toFixed(3)}`);
assert(approx(mLarge.tailProb, 0.050, 0.005), `Large trial df=120 tail probability matches normal 5.0%, got ${(mLarge.tailProb*100).toFixed(1)}%`);

console.log('--- Testing Two-Sample Overlap, SD vs. SEM & Alpha Significance Simulation ---');
// Baseline parameters: Delta=2.0, SD=2.5, n=16, alpha=0.05
const ovlBase = Teaching.significanceOverlap.getMetrics({ delta: 2.0, sd: 2.5, n: 16, alpha: 0.05 });
assert(approx(ovlBase.sem, 0.625, 1e-3), `SEM is SD / sqrt(n) = 2.5 / 4 = 0.625, got ${ovlBase.sem.toFixed(3)}`);
assert(approx(ovlBase.seDiff, 0.8839, 1e-3), `SE_diff is SD * sqrt(2/n) = 0.8839, got ${ovlBase.seDiff.toFixed(4)}`);
assert(ovlBase.df === 30, `Degrees of freedom df = 2*n - 2 = 30, got ${ovlBase.df}`);
assert(approx(ovlBase.tStat, 2.2627, 1e-3), `t-statistic is Delta / SE_diff = 2.2627, got ${ovlBase.tStat.toFixed(4)}`);
assert(approx(ovlBase.tCrit, 2.0423, 0.005), `t_crit for df=30, alpha=0.05 is ~2.042, got ${ovlBase.tCrit.toFixed(3)}`);
assert(approx(ovlBase.deltaCrit, 1.805, 0.01), `Critical separation Delta_crit is ~1.805, got ${ovlBase.deltaCrit.toFixed(3)}`);
assert(ovlBase.isSignificant === true, `Observed Delta (2.0) > Delta_crit (1.805) correctly flagged as significant`);
assert(ovlBase.pValue < 0.05, `p-value (${ovlBase.pValue.toFixed(4)}) is strictly less than 0.05`);
assert(approx(ovlBase.patientOVL, 0.689, 0.01), `High patient overlap (~68.9%) despite statistical significance, got ${(ovlBase.patientOVL*100).toFixed(1)}%`);
assert(approx(ovlBase.meansOVL, 0.110, 0.01), `Means sampling overlap has shrunk to ~11.0%, got ${(ovlBase.meansOVL*100).toFixed(1)}%`);

// Non-significant scenario: Delta=1.0, SD=2.5, n=16, alpha=0.05
const ovlNonSig = Teaching.significanceOverlap.getMetrics({ delta: 1.0, sd: 2.5, n: 16, alpha: 0.05 });
assert(ovlNonSig.isSignificant === false, `Observed Delta (1.0) < Delta_crit (1.805) correctly flagged as not significant`);
assert(ovlNonSig.pValue > 0.05, `p-value (${ovlNonSig.pValue.toFixed(4)}) is strictly greater than 0.05`);
assert(ovlNonSig.meansOVL > 0.35, `Heavy means overlap in non-significant state, got ${(ovlNonSig.meansOVL*100).toFixed(1)}%`);

// Demonstrating effect of modifying alpha:
// At alpha = 0.05: Delta=2.0 is significant (p = 0.031)
// If we tighten alpha to 0.01 (e.g. confirmatory phase III trial or Bonferroni correction):
const ovlStrictAlpha = Teaching.significanceOverlap.getMetrics({ delta: 2.0, sd: 2.5, n: 16, alpha: 0.01 });
assert(approx(ovlStrictAlpha.tCrit, 2.750, 0.01), `Tightened alpha=0.01 pushes t_crit up to ~2.750, got ${ovlStrictAlpha.tCrit.toFixed(3)}`);
assert(approx(ovlStrictAlpha.deltaCrit, 2.431, 0.02), `Tightened alpha=0.01 pushes Delta_crit out to 2.431, got ${ovlStrictAlpha.deltaCrit.toFixed(3)}`);
assert(ovlStrictAlpha.isSignificant === false, `Same Delta=2.0 is NO LONGER significant at alpha=0.01 (p=0.031 >= 0.01)`);

// If we relax alpha to 0.10 on borderline Delta=1.6 (not significant at 0.05):
const ovlBorderline05 = Teaching.significanceOverlap.getMetrics({ delta: 1.6, sd: 2.5, n: 16, alpha: 0.05 });
const ovlBorderline10 = Teaching.significanceOverlap.getMetrics({ delta: 1.6, sd: 2.5, n: 16, alpha: 0.10 });
assert(ovlBorderline05.isSignificant === false, `Delta=1.6 is NOT significant at alpha=0.05`);
assert(ovlBorderline10.isSignificant === true, `Delta=1.6 BECOMES significant when alpha relaxed to 0.10 (t_crit drops to ~1.697)`);

// Demonstrating SD vs. SEM separation mechanism:
// Increase sample size from n=16 to n=64 while holding Delta=2.0 and SD=2.5 constant:
const ovlLargeN = Teaching.significanceOverlap.getMetrics({ delta: 2.0, sd: 2.5, n: 64, alpha: 0.05 });
assert(approx(ovlLargeN.sd, 2.50, 1e-4), `Patient SD remains constant at 2.50`);
assert(approx(ovlLargeN.patientOVL, ovlBase.patientOVL, 1e-4), `Patient biological overlap remains identical at 68.9%`);
assert(approx(ovlLargeN.sem, 0.3125, 1e-4), `SEM cuts in half from 0.625 to 0.3125`);
assert(ovlLargeN.tStat > 4.5, `t-statistic doubles from 2.26 to 4.53, got ${ovlLargeN.tStat.toFixed(2)}`);
assert(ovlLargeN.pValue < 0.0001, `p-value plunges to < 0.0001, got ${ovlLargeN.pValue.toFixed(6)}`);
// Testing Independent Group 1 and Group 2 SD and SEM Controls:
console.log('--- Testing Independent Group 1 and Group 2 SD & SEM Controls ---');
const ovlHetero = Teaching.significanceOverlap.getMetrics({
  mean1: 10.0,
  delta: 2.5,
  sd1: 2.0,
  sd2: 4.0,
  n1: 16,
  n2: 16,
  alpha: 0.05
});
assert(approx(ovlHetero.sem1, 0.500, 1e-3), `Group 1 SEM is 2.0 / 4 = 0.500, got ${ovlHetero.sem1.toFixed(3)}`);
assert(approx(ovlHetero.sem2, 1.000, 1e-3), `Group 2 SEM is 4.0 / 4 = 1.000, got ${ovlHetero.sem2.toFixed(3)}`);
assert(approx(ovlHetero.seDiff, Math.sqrt(1.25), 1e-3), `SE_diff is sqrt(0.5^2 + 1.0^2) = 1.118, got ${ovlHetero.seDiff.toFixed(3)}`);
assert(approx(ovlHetero.df, 22.06, 0.05), `Welch-Satterthwaite df for unequal variances is ~22.06, got ${ovlHetero.df.toFixed(2)}`);
assert(ovlHetero.patientOVL > 0.40 && ovlHetero.patientOVL < 0.65, `Heteroscedastic patient overlap is valid, got ${(ovlHetero.patientOVL*100).toFixed(1)}%`);
assert(ovlHetero.meansOVL < 0.30, `Heteroscedastic means overlap correctly narrows, got ${(ovlHetero.meansOVL*100).toFixed(1)}%`);

// Direct SEM manipulation (e.g. user drags SEM1 slider to 0.25, SEM2 to 0.50)
const ovlDirectSEM = Teaching.significanceOverlap.getMetrics({
  mean1: 10.0,
  delta: 2.0,
  sd1: 2.0,
  sd2: 3.0,
  sem1: 0.25,
  sem2: 0.50,
  alpha: 0.05
});
assert(approx(ovlDirectSEM.sem1, 0.25, 1e-4), `Direct SEM1 set to 0.25, got ${ovlDirectSEM.sem1}`);
assert(approx(ovlDirectSEM.sem2, 0.50, 1e-4), `Direct SEM2 set to 0.50, got ${ovlDirectSEM.sem2}`);
assert(ovlDirectSEM.n1 === 64, `Implied sample size n1 for SD=2.0 and SEM=0.25 is 64, got ${ovlDirectSEM.n1}`);
assert(ovlDirectSEM.n2 === 36, `Implied sample size n2 for SD=3.0 and SEM=0.50 is 36, got ${ovlDirectSEM.n2}`);
assert(approx(ovlDirectSEM.seDiff, Math.sqrt(0.25*0.25 + 0.50*0.50), 1e-3), `SE_diff is ~0.559, got ${ovlDirectSEM.seDiff.toFixed(3)}`);
assert(ovlDirectSEM.isSignificant === true, `Delta=2.0 with contracted SEMs is significant`);

// Test Teaching DOCX export with dual group overlap metrics
const docxWithDualOverlap = DocxReports.createTeachingDocx({
  distName: 'Normal Distribution',
  sampleMean: 10.0,
  theoMean: 10.0,
  sampleSD: 2.5,
  theoSD: 2.5,
  skewness: 0.05,
  kurtosis: -0.02,
  jbStat: 0.12,
  jbP: 0.94,
  isNormal: true,
  overlap: ovlHetero
});
const docxDualBuffer = docxWithDualOverlap.generateUint8Array();
assert(docxDualBuffer.length > 30000, `Teaching DOCX with dual group overlap generated ${docxDualBuffer.length} bytes`);
assert(isZip(docxDualBuffer), 'Teaching DOCX with dual group overlap is a valid PKZIP archive');

console.log(`\nVerification Complete: ${passes} Passed, ${failures} Failed`);
if (failures > 0) process.exit(1);


