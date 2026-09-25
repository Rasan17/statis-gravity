/**
 * Comprehensive Calculation Audit Generator for Statis-Gravity
 * Extracts calculations from all core modules to compare against Python scipy/numpy/statsmodels/sklearn
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

import { Distributions } from '../js/stats/distributions.js';
import { Descriptive } from '../js/stats/descriptive.js';
import { Hypothesis } from '../js/stats/hypothesis.js';
import { Anova } from '../js/stats/anova.js';
import { Categorical } from '../js/stats/categorical.js';
import { Correlation } from '../js/stats/correlation.js';
import { Diagnostic } from '../js/stats/diagnostic.js';
import { PowerAnalysis } from '../js/stats/power.js';
import { Psm } from '../js/stats/psm.js';
import { Multivariate } from '../js/stats/multivariate.js';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const auditData = {
  metadata: {
    timestamp: new Date().toISOString(),
    suite: 'Statis-Gravity Mathematical Calculation Audit'
  },
  distributions: {},
  descriptive: {},
  hypothesis: {},
  anova: {},
  categorical: {},
  correlation: {},
  diagnostic: {},
  power: {},
  multivariate: {},
  psm: {}
};

// 1. Distributions Audit
const zValues = [-2.576, -1.96, -1.645, 0.0, 1.0, 1.645, 1.96, 2.576];
auditData.distributions.normal = zValues.map(z => ({
  z,
  cdf: Distributions.normalCDF(z),
  pTwoTailed: Distributions.normalPValue(z)
}));

const tValues = [
  { t: 2.228, df: 10 },
  { t: -1.753, df: 15 },
  { t: 3.169, df: 10 },
  { t: 2.042, df: 30 }
];
auditData.distributions.studentT = tValues.map(item => ({
  t: item.t,
  df: item.df,
  pTwoTailed: Distributions.tPValue(item.t, item.df)
}));

const chi2Values = [
  { chi2: 3.841, df: 1 },
  { chi2: 5.991, df: 2 },
  { chi2: 9.488, df: 4 },
  { chi2: 18.307, df: 10 }
];
auditData.distributions.chiSquare = chi2Values.map(item => ({
  chi2: item.chi2,
  df: item.df,
  pValue: Distributions.chiSquarePValue(item.chi2, item.df)
}));

const fValues = [
  { f: 3.885, df1: 2, df2: 12 },
  { f: 2.689, df1: 3, df2: 20 },
  { f: 4.459, df1: 1, df2: 18 }
];
auditData.distributions.fisherF = fValues.map(item => ({
  f: item.f,
  df1: item.df1,
  df2: item.df2,
  pValue: Distributions.fPValue(item.f, item.df1, item.df2)
}));

const probValues = [0.005, 0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975, 0.995];
auditData.distributions.inverseNormal = probValues.map(p => ({
  p,
  z: Distributions.invNormalCDF(p)
}));

// 2. Descriptive Statistics Audit
const sample1 = [14.2, 15.1, 13.8, 16.5, 14.9, 15.8, 17.2, 13.5, 15.0, 14.6, 16.1, 14.8, 15.4, 16.0];
const sample2 = [22.4, 28.1, 18.9, 31.5, 24.0, 19.8, 35.2, 26.7, 21.3, 29.4, 33.1, 20.5];
const sampleSkewed = [120, 135, 140, 145, 150, 155, 160, 290, 310, 480];

auditData.descriptive.sample1 = {
  data: sample1,
  results: Descriptive.calculate(sample1)
};
auditData.descriptive.sample2 = {
  data: sample2,
  results: Descriptive.calculate(sample2)
};
auditData.descriptive.sampleSkewed = {
  data: sampleSkewed,
  results: Descriptive.calculate(sampleSkewed)
};

// 3. Hypothesis Testing Audit
auditData.hypothesis.independentT = {
  groupA: sample1,
  groupB: sample2,
  results: Hypothesis.independentTTest(sample1, sample2)
};

auditData.hypothesis.welchT = {
  groupA: sample1,
  groupB: sample2,
  results: Hypothesis.welchTTest(sample1, sample2)
};

const pairedPre = [18.5, 22.1, 19.4, 25.0, 17.8, 24.2, 21.0, 26.5, 20.1, 23.4];
const pairedPost = [12.1, 14.5, 13.0, 16.2, 11.5, 15.8, 14.0, 17.1, 13.2, 15.0];

auditData.hypothesis.pairedT = {
  pre: pairedPre,
  post: pairedPost,
  results: Hypothesis.pairedTTest(pairedPre, pairedPost)
};

auditData.hypothesis.mannWhitney = {
  groupA: sample1,
  groupB: sample2,
  results: Hypothesis.mannWhitneyUTest(sample1, sample2)
};

auditData.hypothesis.wilcoxon = {
  pre: pairedPre,
  post: pairedPost,
  results: Hypothesis.wilcoxonSignedRank(pairedPre, pairedPost)
};

auditData.hypothesis.assumptionsIndep = {
  groupA: sample1,
  groupB: sample2,
  results: Hypothesis.evaluateAssumptions(sample1, sample2, false)
};

auditData.hypothesis.assumptionsPaired = {
  pre: pairedPre,
  post: pairedPost,
  results: Hypothesis.evaluateAssumptions(pairedPre, pairedPost, true)
};

auditData.hypothesis.levene = {
  groups: [
    { name: 'Group 1', data: sample1 },
    { name: 'Group 2', data: sample2 }
  ],
  results: Anova.leveneTest([
    { name: 'Group 1', data: sample1 },
    { name: 'Group 2', data: sample2 }
  ])
};

// 4. ANOVA & Multi-Group Audit
const anovaG1 = [12.1, 14.3, 13.5, 15.0, 14.2, 13.8, 16.1, 14.5, 15.2, 13.9];
const anovaG2 = [18.4, 20.1, 19.2, 22.0, 21.3, 20.5, 23.1, 19.8, 21.0, 22.4];
const anovaG3 = [25.0, 28.2, 27.1, 30.5, 29.0, 26.8, 31.4, 28.5, 29.2, 32.1];

const anovaGroups = [
  { name: 'Group 1 (Control)', data: anovaG1 },
  { name: 'Group 2 (Treatment Low)', data: anovaG2 },
  { name: 'Group 3 (Treatment High)', data: anovaG3 }
];

auditData.anova.oneWay = {
  groups: anovaGroups,
  results: Anova.oneWay(anovaGroups)
};

auditData.anova.welch = {
  groups: anovaGroups,
  results: Anova.welch(anovaGroups)
};

auditData.anova.kruskalWallis = {
  groups: anovaGroups,
  results: Anova.kruskalWallis(anovaGroups)
};

auditData.anova.repeatedMeasures = {
  groups: anovaGroups,
  results: Anova.repeatedMeasures(anovaGroups)
};

auditData.anova.friedman = {
  groups: anovaGroups,
  results: Anova.friedman(anovaGroups)
};

// 5. Categorical 2x2 Audit
const catTables = [
  { name: 'Shunt Infection Trial', a: 14, b: 86, c: 32, d: 68 },
  { name: 'Tumor Recurrence', a: 45, b: 55, c: 20, d: 80 },
  { name: 'Rare Event Study', a: 3, b: 97, c: 12, d: 88 }
];

auditData.categorical.tables = catTables.map(t => ({
  name: t.name,
  a: t.a,
  b: t.b,
  c: t.c,
  d: t.d,
  results: Categorical.twoByTwo(t.a, t.b, t.c, t.d)
}));

// 6. Correlation & Regression Audit
const corrX = [10, 12, 14, 15, 18, 20, 22, 24, 25, 28, 30, 32];
const corrY = [15, 19, 21, 22, 27, 31, 33, 37, 39, 42, 47, 49];

auditData.correlation = {
  x: corrX,
  y: corrY,
  pearson: Correlation.pearson(corrX, corrY),
  spearman: Correlation.spearman(corrX, corrY),
  regression: Correlation.linearRegression(corrX, corrY)
};

// 7. Diagnostic & ROC Audit
const diag2x2 = { tp: 85, fp: 12, fn: 15, tn: 188 };
auditData.diagnostic.eval2x2 = {
  table: diag2x2,
  results: Diagnostic.evaluate2x2(diag2x2.tp, diag2x2.fp, diag2x2.fn, diag2x2.tn)
};

const rocSample = [
  { score: 0.95, status: 1 }, { score: 0.88, status: 1 }, { score: 0.82, status: 1 },
  { score: 0.79, status: 1 }, { score: 0.75, status: 1 }, { score: 0.71, status: 0 },
  { score: 0.68, status: 1 }, { score: 0.62, status: 0 }, { score: 0.58, status: 1 },
  { score: 0.55, status: 0 }, { score: 0.51, status: 0 }, { score: 0.47, status: 1 },
  { score: 0.42, status: 0 }, { score: 0.38, status: 0 }, { score: 0.35, status: 0 },
  { score: 0.29, status: 0 }, { score: 0.24, status: 0 }, { score: 0.18, status: 0 }
];

auditData.diagnostic.roc = {
  samples: rocSample,
  results: Diagnostic.computeROC(rocSample)
};

// 8. Power Analysis Audit
auditData.power.indepMeans = {
  params: { m1: 10, m2: 15, sd: 10, alpha: 0.05, power: 0.80 },
  resultN: PowerAnalysis.twoIndependentMeans({ m1: 10, m2: 15, sd: 10, alpha: 0.05, power: 0.80 }),
  resultPower: PowerAnalysis.twoIndependentMeans({ m1: 10, m2: 15, sd: 10, alpha: 0.05, nPerGroup: 64 })
};

auditData.power.pairedMeans = {
  params: { m1: 120, m2: 112, sdDiff: 10, alpha: 0.05, power: 0.80 },
  resultN: PowerAnalysis.pairedMeans({ m1: 120, m2: 112, sdDiff: 10, alpha: 0.05, power: 0.80 }),
  resultPower: PowerAnalysis.pairedMeans({ m1: 120, m2: 112, sdDiff: 10, alpha: 0.05, sdDiff: 10, nPairs: 15 })
};

// 9. Multivariate Analysis (PCA) Audit
const mixedCohort = Multivariate.getSampleMixedCohort();
const pcaContCols = ['age_months', 'cranial_index', 'operative_time_min', 'blood_loss_ml', 'length_of_stay_days'];
const pcaRes = Multivariate.runPCA(mixedCohort, pcaContCols);

auditData.multivariate.pca = {
  records: mixedCohort.map(r => {
    const obj = {};
    pcaContCols.forEach(c => obj[c] = r[c]);
    return obj;
  }),
  cols: pcaContCols,
  results: {
    eigenvalues: pcaRes.scree.map(s => s.eigenvalue),
    variancePct: pcaRes.scree.map(s => s.variancePct),
    cumulativePct: pcaRes.scree.map(s => s.cumulativePct),
    componentLoadings: pcaRes.loadings,
    individualCoordSample: pcaRes.individuals.slice(0, 5).map(ind => ind.coords)
  }
};

// 10. Propensity Score Matching (PSM) Logistic Regression & Balance
const psmDataset = Psm.getSampleClinicalDataset();
const psmCovars = ['age', 'sex', 'comorbidity_score', 'baseline_severity'];
const psmAnalysis = Psm.executeAnalysis(psmDataset, {
  treatmentCol: 'treatment_col',
  outcomeCol: 'outcome_col',
  covariateCols: psmCovars,
  caliperMultiplier: 0.20,
  enforceCommonSupport: true
});

auditData.psm = {
  dataset: psmDataset,
  treatmentCol: 'treatment_col',
  outcomeCol: 'outcome_col',
  covariateCols: psmCovars,
  results: {
    nTotal: psmAnalysis.prep.nClean,
    nTreated: psmAnalysis.prep.nTreated,
    nControl: psmAnalysis.prep.nControl,
    nMatchedPairs: psmAnalysis.nMatchedPairs,
    coefficients: psmAnalysis.logisticRegression?.coefficients || [],
    pseudoR2: psmAnalysis.logisticRegression?.pseudoR2 || 0,
    balanceTable: psmAnalysis.balance.balanceTable,
    att: psmAnalysis.outcome.att,
    ci95: psmAnalysis.outcome.ci95,
    pValue: psmAnalysis.outcome.pValue
  }
};

const outputPath = path.join(__dirname, 'audit_data.json');
fs.writeFileSync(outputPath, JSON.stringify(auditData, null, 2));
console.log(`Generated JS audit data successfully at ${outputPath}`);
