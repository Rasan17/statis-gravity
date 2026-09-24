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

console.log('--- Testing Power Analysis ---');
// d = 0.5, alpha = 0.05, power = 0.80 -> N per group is ~64
const pwr = PowerAnalysis.sampleSizeMeans(10, 15, 10, 0.05, 0.80);
assert(pwr.nPerGroup >= 60 && pwr.nPerGroup <= 70, `Sample size per group is ~64, got ${pwr.nPerGroup}`);

console.log(`\nVerification Complete: ${passes} Passed, ${failures} Failed`);
if (failures > 0) process.exit(1);
