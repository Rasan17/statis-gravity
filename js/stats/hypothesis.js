/**
 * Statis-Gravity - Hypothesis Testing Module
 * Parametric & Non-Parametric two-group comparisons.
 */

import { Distributions } from './distributions.js';
import { Descriptive } from './descriptive.js';

export const Hypothesis = {
  /**
   * Independent Two-Sample Student's t-test (Equal Variances assumed)
   */
  independentTTest(groupA, groupB) {
    const statsA = Descriptive.calculate(groupA);
    const statsB = Descriptive.calculate(groupB);

    if (statsA.n < 2 || statsB.n < 2) {
      return { error: 'Both groups must contain at least 2 observations.' };
    }

    const n1 = statsA.n;
    const n2 = statsB.n;
    const df = n1 + n2 - 2;

    // Pooled variance
    const sp2 = ((n1 - 1) * statsA.variance + (n2 - 1) * statsB.variance) / df;
    const seDiff = Math.sqrt(sp2 * (1 / n1 + 1 / n2));

    const meanDiff = statsA.mean - statsB.mean;
    const t = seDiff === 0 ? 0 : meanDiff / seDiff;
    const pValue = Distributions.tPValue(t, df);

    // Cohen's d effect size
    const sp = Math.sqrt(sp2);
    const cohensD = sp === 0 ? 0 : meanDiff / sp;

    // 95% CI for difference
    const z95 = Distributions.invNormalCDF(0.975);
    const tFactor = df > 30 ? z95 : z95 * (1 + 1 / (4 * df));
    const ci95 = [meanDiff - tFactor * seDiff, meanDiff + tFactor * seDiff];

    return {
      testName: "Independent Samples Student's t-test",
      groupA: statsA,
      groupB: statsB,
      meanDiff,
      seDiff,
      df,
      statistic: t,
      pValue,
      cohensD,
      ci95,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Welch's t-test (Unequal Variances, Satterthwaite approximation)
   */
  welchTTest(groupA, groupB) {
    const statsA = Descriptive.calculate(groupA);
    const statsB = Descriptive.calculate(groupB);

    if (statsA.n < 2 || statsB.n < 2) {
      return { error: 'Both groups must contain at least 2 observations.' };
    }

    const n1 = statsA.n;
    const n2 = statsB.n;
    const v1 = statsA.variance / n1;
    const v2 = statsB.variance / n2;

    const seDiff = Math.sqrt(v1 + v2);
    const meanDiff = statsA.mean - statsB.mean;
    const t = seDiff === 0 ? 0 : meanDiff / seDiff;

    // Welch-Satterthwaite degrees of freedom
    const num = Math.pow(v1 + v2, 2);
    const den = Math.pow(v1, 2) / (n1 - 1) + Math.pow(v2, 2) / (n2 - 1);
    const df = den === 0 ? 1 : num / den;

    const pValue = Distributions.tPValue(t, df);

    // Cohen's d
    const pooledSD = Math.sqrt((statsA.variance + statsB.variance) / 2);
    const cohensD = pooledSD === 0 ? 0 : meanDiff / pooledSD;

    const z95 = Distributions.invNormalCDF(0.975);
    const tFactor = df > 30 ? z95 : z95 * (1 + 1 / (4 * df));
    const ci95 = [meanDiff - tFactor * seDiff, meanDiff + tFactor * seDiff];

    return {
      testName: "Welch's t-test (Unequal Variances)",
      groupA: statsA,
      groupB: statsB,
      meanDiff,
      seDiff,
      df,
      statistic: t,
      pValue,
      cohensD,
      ci95,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Paired Samples t-test (Dependent samples)
   */
  pairedTTest(groupPre, groupPost) {
    const a = Descriptive.cleanData(groupPre);
    const b = Descriptive.cleanData(groupPost);

    const n = Math.min(a.length, b.length);
    if (n < 2) {
      return { error: 'Paired test requires at least 2 valid paired observations.' };
    }

    const diffs = [];
    for (let i = 0; i < n; i++) {
      diffs.push(b[i] - a[i]); // Post - Pre
    }

    const diffStats = Descriptive.calculate(diffs);
    const df = n - 1;
    const t = diffStats.sem === 0 ? 0 : diffStats.mean / diffStats.sem;
    const pValue = Distributions.tPValue(t, df);
    const cohensD = diffStats.sd === 0 ? 0 : diffStats.mean / diffStats.sd;

    return {
      testName: "Paired Samples t-test",
      n,
      meanDiff: diffStats.mean,
      sdDiff: diffStats.sd,
      seDiff: diffStats.sem,
      df,
      statistic: t,
      pValue,
      cohensD,
      ci95: diffStats.ci95,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Mann-Whitney U test (Wilcoxon Rank-Sum test for independent samples)
   */
  mannWhitneyUTest(groupA, groupB) {
    const a = Descriptive.cleanData(groupA);
    const b = Descriptive.cleanData(groupB);

    const n1 = a.length;
    const n2 = b.length;
    if (n1 === 0 || n2 === 0) {
      return { error: 'Both groups must contain at least 1 observation.' };
    }

    // Combine and rank with mid-ranks for ties
    const combined = [
      ...a.map(val => ({ val, group: 'A' })),
      ...b.map(val => ({ val, group: 'B' }))
    ].sort((x, y) => x.val - y.val);

    const N = combined.length;
    let i = 0;
    while (i < N) {
      let j = i;
      while (j < N - 1 && combined[j + 1].val === combined[i].val) {
        j++;
      }
      const rank = (i + 1 + j + 1) / 2;
      for (let k = i; k <= j; k++) {
        combined[k].rank = rank;
      }
      i = j + 1;
    }

    let rankSumA = 0;
    let rankSumB = 0;
    for (const item of combined) {
      if (item.group === 'A') rankSumA += item.rank;
      else rankSumB += item.rank;
    }

    const u1 = rankSumA - (n1 * (n1 + 1)) / 2;
    const u2 = rankSumB - (n2 * (n2 + 1)) / 2;
    const U = Math.min(u1, u2);

    // Normal approximation for large sample U
    const meanU = (n1 * n2) / 2;
    const sigmaU = Math.sqrt((n1 * n2 * (n1 + n2 + 1)) / 12);
    const z = sigmaU === 0 ? 0 : (U - meanU) / sigmaU;
    const pValue = Distributions.normalPValue(z);

    // Common Language Effect Size (CLES) / Rank-biserial correlation
    const cles = (u1) / (n1 * n2);
    const rankBiserial = 1 - (2 * U) / (n1 * n2);

    return {
      testName: 'Mann-Whitney U Test (Wilcoxon Rank-Sum)',
      n1,
      n2,
      rankSumA,
      rankSumB,
      u1,
      u2,
      statistic: U,
      zScore: z,
      pValue,
      cles,
      rankBiserial,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Wilcoxon Signed-Rank test (Paired non-parametric comparison)
   */
  wilcoxonSignedRank(groupPre, groupPost) {
    const a = Descriptive.cleanData(groupPre);
    const b = Descriptive.cleanData(groupPost);

    const len = Math.min(a.length, b.length);
    if (len < 3) {
      return { error: 'Wilcoxon signed-rank requires at least 3 paired observations.' };
    }

    // Differences discarding zeros
    const diffs = [];
    for (let i = 0; i < len; i++) {
      const d = b[i] - a[i];
      if (d !== 0) {
        diffs.push({ diff: d, absDiff: Math.abs(d) });
      }
    }

    const n = diffs.length;
    if (n < 3) {
      return { error: 'Insufficient non-zero differences for Wilcoxon test.' };
    }

    diffs.sort((x, y) => x.absDiff - y.absDiff);

    // Assign mid-ranks to ties
    let i = 0;
    while (i < n) {
      let j = i;
      while (j < n - 1 && diffs[j + 1].absDiff === diffs[i].absDiff) {
        j++;
      }
      const rank = (i + 1 + j + 1) / 2;
      for (let k = i; k <= j; k++) {
        diffs[k].rank = rank;
      }
      i = j + 1;
    }

    let wPlus = 0;
    let wMinus = 0;
    for (const d of diffs) {
      if (d.diff > 0) wPlus += d.rank;
      else wMinus += d.rank;
    }

    const W = Math.min(wPlus, wMinus);
    const meanW = (n * (n + 1)) / 4;
    const sigmaW = Math.sqrt((n * (n + 1) * (2 * n + 1)) / 24);
    const z = sigmaW === 0 ? 0 : (W - meanW) / sigmaW;
    const pValue = Distributions.normalPValue(z);

    return {
      testName: 'Wilcoxon Signed-Rank Test',
      n,
      wPlus,
      wMinus,
      statistic: W,
      zScore: z,
      pValue,
      isSignificant: pValue < 0.05
    };
  }
};
