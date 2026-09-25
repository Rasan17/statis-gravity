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

    // Matched-pairs rank-biserial correlation effect size: r = (W+ - W-) / (W+ + W-)
    const totalRankSum = wPlus + wMinus;
    const rankBiserial = totalRankSum === 0 ? 0 : (wPlus - wMinus) / totalRankSum;

    const rawDiffs = diffs.map(d => d.diff);
    const diffStats = Descriptive.calculate(rawDiffs);

    return {
      testName: 'Wilcoxon Signed-Rank Test',
      n,
      wPlus,
      wMinus,
      statistic: W,
      zScore: z,
      pValue,
      rankBiserial,
      meanDiff: diffStats.mean,
      medianDiff: diffStats.median,
      sdDiff: diffStats.sd,
      seDiff: diffStats.sem,
      ci95: diffStats.ci95,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Evaluates statistical assumptions (Normality & Equality of Variances)
   * and recommends the optimal test depending on sample design (paired vs independent).
   */
  evaluateAssumptions(groupA, groupB, isPaired = false) {
    const a = Descriptive.cleanData(groupA);
    const b = Descriptive.cleanData(groupB);

    if (isPaired) {
      const len = Math.min(a.length, b.length);
      if (len < 2) {
        return {
          error: 'Paired comparison requires at least 2 valid paired data points.',
          isPaired: true
        };
      }
      const unequalLengths = a.length !== b.length;
      const diffs = [];
      for (let i = 0; i < len; i++) {
        diffs.push(b[i] - a[i]);
      }
      const statsDiff = Descriptive.calculate(diffs);
      const statsA = Descriptive.calculate(a.slice(0, len));
      const statsB = Descriptive.calculate(b.slice(0, len));

      const jbDiff = statsDiff.normality || { isNormal: true, pValue: 1.0, statistic: 0 };
      const skewAlert = Math.abs(statsDiff.skewness) > 1.0;
      const kurtAlert = Math.abs(statsDiff.kurtosis) > 1.5;
      const isNormal = len >= 10 ? jbDiff.isNormal : (!skewAlert && !kurtAlert);

      const recommendedTest = isNormal ? 'paired' : 'wilcoxon';
      const recommendedTestName = isNormal 
        ? "Paired Samples Student's t-test (Parametric)" 
        : "Wilcoxon Signed-Rank Test (Non-parametric)";

      let rationale = `Paired / repeated measures design (${len} paired observations). `;
      if (isNormal) {
        rationale += `Within-subject differences (Δ = Post - Pre) conform to normality (Jarque-Bera p = ${jbDiff.pValue > 0.001 ? jbDiff.pValue.toFixed(3) : '< .001'}, Skewness = ${statsDiff.skewness.toFixed(2)}). The parametric Paired Samples t-test is recommended.`;
      } else {
        rationale += `Within-subject differences (Δ) deviate significantly from normality (Jarque-Bera p = ${jbDiff.pValue > 0.001 ? jbDiff.pValue.toFixed(3) : '< .001'}, Skewness = ${statsDiff.skewness.toFixed(2)}, Excess Kurtosis = ${statsDiff.kurtosis.toFixed(2)}). The non-parametric Wilcoxon Signed-Rank test is recommended.`;
      }

      return {
        isPaired: true,
        unequalLengths,
        n: len,
        nA: a.length,
        nB: b.length,
        statsA,
        statsB,
        statsDiff,
        normality: {
          isNormal,
          testName: 'Jarque-Bera Test of Paired Differences (Δ)',
          statistic: jbDiff.statistic,
          pValue: jbDiff.pValue,
          skewness: statsDiff.skewness,
          kurtosis: statsDiff.kurtosis,
          interpretation: isNormal ? 'Normal Distribution (Parametric suitable)' : 'Non-Normal Distribution (Non-parametric recommended)'
        },
        varianceEquality: {
          applicable: false,
          note: 'Not required for paired repeated measures (within-subject differencing removes inter-subject variance).'
        },
        recommendedTest,
        recommendedTestName,
        rationale
      };
    } else {
      if (a.length < 2 || b.length < 2) {
        return {
          error: 'Independent comparison requires at least 2 observations in each cohort.',
          isPaired: false
        };
      }
      const statsA = Descriptive.calculate(a);
      const statsB = Descriptive.calculate(b);

      const jbA = statsA.normality || { isNormal: true, pValue: 1.0, statistic: 0 };
      const jbB = statsB.normality || { isNormal: true, pValue: 1.0, statistic: 0 };
      const normA = statsA.n >= 10 ? jbA.isNormal : (Math.abs(statsA.skewness) <= 1.0 && Math.abs(statsA.kurtosis) <= 1.5);
      const normB = statsB.n >= 10 ? jbB.isNormal : (Math.abs(statsB.skewness) <= 1.0 && Math.abs(statsB.kurtosis) <= 1.5);
      const isParametric = normA && normB;

      const vA = statsA.variance;
      const vB = statsB.variance;
      let fStat = 1.0;
      let df1 = statsA.n - 1;
      let df2 = statsB.n - 1;
      let varPVal = 1.0;

      if (vA > 0 && vB > 0) {
        if (vA >= vB) {
          fStat = vA / vB;
          df1 = statsA.n - 1;
          df2 = statsB.n - 1;
        } else {
          fStat = vB / vA;
          df1 = statsB.n - 1;
          df2 = statsA.n - 1;
        }
        varPVal = Math.min(1.0, 2 * Distributions.fPValue(fStat, df1, df2));
      }
      const equalVariance = varPVal > 0.05;

      let recommendedTest;
      let recommendedTestName;
      let rationale = `Independent two-cohort design (n₁ = ${statsA.n}, n₂ = ${statsB.n}). `;

      if (isParametric) {
        if (equalVariance) {
          recommendedTest = 'student';
          recommendedTestName = "Student's t-test (Equal Variances)";
          rationale += `Both cohorts conform to normal distributions (Group 1 p = ${jbA.pValue > 0.001 ? jbA.pValue.toFixed(3) : '< .001'}, Group 2 p = ${jbB.pValue > 0.001 ? jbB.pValue.toFixed(3) : '< .001'}) and variance homogeneity is preserved (F(${df1}, ${df2}) = ${fStat.toFixed(2)}, p = ${varPVal > 0.001 ? varPVal.toFixed(3) : '< .001'}). Standard Student's t-test is appropriate.`;
        } else {
          recommendedTest = 'welch';
          recommendedTestName = "Welch's t-test (Unequal Variances - Recommended)";
          rationale += `Both cohorts conform to normal distributions, but variance homogeneity is violated (Heteroscedasticity: F(${df1}, ${df2}) = ${fStat.toFixed(2)}, p = ${varPVal > 0.001 ? varPVal.toFixed(3) : '< .001'}). Welch's t-test with Satterthwaite degrees of freedom is required to control Type I error rates.`;
        }
      } else {
        recommendedTest = 'mannwhitney';
        recommendedTestName = 'Mann-Whitney U Test (Non-parametric)';
        rationale += `Deviation from normality detected (${!normA ? 'Group 1 non-normal (p = ' + (jbA.pValue > 0.001 ? jbA.pValue.toFixed(3) : '< .001') + ')' : ''}${!normA && !normB ? '; ' : ''}${!normB ? 'Group 2 non-normal (p = ' + (jbB.pValue > 0.001 ? jbB.pValue.toFixed(3) : '< .001') + ')' : ''}). Non-parametric rank-sum comparison via Mann-Whitney U test is recommended.`;
      }

      return {
        isPaired: false,
        nA: statsA.n,
        nB: statsB.n,
        statsA,
        statsB,
        normality: {
          isParametric,
          normA,
          normB,
          jbA,
          jbB,
          interpretation: isParametric ? 'Both Cohorts Normal (Parametric suitable)' : 'Normality Violated (Non-parametric recommended)'
        },
        varianceEquality: {
          applicable: true,
          equalVariance,
          fStat,
          df1,
          df2,
          pValue: varPVal,
          varA: vA,
          varB: vB,
          interpretation: equalVariance ? 'Homoscedastic (Equal Variances, p > .05)' : 'Heteroscedastic (Unequal Variances, p ≤ .05)'
        },
        recommendedTest,
        recommendedTestName,
        rationale
      };
    }
  }
};
