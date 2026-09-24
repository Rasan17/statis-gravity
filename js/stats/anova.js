/**
 * Statis-Gravity - ANOVA & Multi-Group Comparison Module
 * One-Way ANOVA, Kruskal-Wallis H Test, Effect Sizes (Eta², Omega²), and Pairwise Tukey HSD.
 */

import { Distributions } from './distributions.js';
import { Descriptive } from './descriptive.js';

export const Anova = {
  /**
   * One-Way Analysis of Variance (ANOVA)
   * @param {Array<{name: string, data: number[]}>} groups
   */
  oneWay(groups) {
    if (!Array.isArray(groups) || groups.length < 2) {
      return { error: 'ANOVA requires at least 2 distinct groups.' };
    }

    const processedGroups = groups.map(g => {
      const stats = Descriptive.calculate(g.data);
      return {
        name: g.name || 'Group',
        stats,
        data: stats.values
      };
    }).filter(g => g.stats.n > 0);

    const k = processedGroups.length;
    if (k < 2) {
      return { error: 'At least 2 groups must have valid numerical data.' };
    }

    const totalN = processedGroups.reduce((acc, g) => acc + g.stats.n, 0);
    const grandSum = processedGroups.reduce((acc, g) => acc + g.stats.sum, 0);
    const grandMean = grandSum / totalN;

    // Degrees of Freedom
    const dfBetween = k - 1;
    const dfWithin = totalN - k;
    const dfTotal = totalN - 1;

    if (dfWithin <= 0) {
      return { error: 'Insufficient degrees of freedom for within-group variance.' };
    }

    // Sum of Squares
    let ssBetween = 0;
    let ssWithin = 0;

    for (const g of processedGroups) {
      ssBetween += g.stats.n * Math.pow(g.stats.mean - grandMean, 2);
      for (const val of g.data) {
        ssWithin += Math.pow(val - g.stats.mean, 2);
      }
    }
    const ssTotal = ssBetween + ssWithin;

    // Mean Squares
    const msBetween = ssBetween / dfBetween;
    const msWithin = ssWithin / dfWithin;

    // F Statistic & p-value
    const F = msWithin === 0 ? 0 : msBetween / msWithin;
    const pValue = Distributions.fPValue(F, dfBetween, dfWithin);

    // Effect Sizes: Eta-squared and Omega-squared
    const etaSquared = ssTotal === 0 ? 0 : ssBetween / ssTotal;
    const omegaSquared = (ssTotal + msWithin) === 0 ? 0 :
      (ssBetween - dfBetween * msWithin) / (ssTotal + msWithin);

    // Pairwise Comparisons (Tukey HSD approximation with Studentized Range Q)
    const pairwise = [];
    const numPairs = (k * (k - 1)) / 2;
    const pooledSD = Math.sqrt(msWithin);

    for (let i = 0; i < k; i++) {
      for (let j = i + 1; j < k; j++) {
        const gA = processedGroups[i];
        const gB = processedGroups[j];
        const meanDiff = gA.stats.mean - gB.stats.mean;
        const seDiff = Math.sqrt(msWithin * (1 / gA.stats.n + 1 / gB.stats.n));
        const seTukey = Math.sqrt((msWithin / 2) * (1 / gA.stats.n + 1 / gB.stats.n));
        const q = seTukey === 0 ? 0 : Math.abs(meanDiff) / seTukey;
        const tEquivalent = q / Math.SQRT2;
        const pPair = Distributions.tPValue(tEquivalent, dfWithin);
        const pAdjusted = Math.min(1.0, pPair * numPairs);

        const z95 = 1.95996;
        const tCrit = dfWithin > 30 ? z95 : z95 * (1 + 1 / (4 * dfWithin));
        const margin = tCrit * seDiff;
        const ci95 = [meanDiff - margin, meanDiff + margin];
        const cohensD = pooledSD === 0 ? 0 : meanDiff / pooledSD;

        pairwise.push({
          groupA: gA.name,
          groupB: gB.name,
          comparison: `${gA.name} vs ${gB.name}`,
          meanA: gA.stats.mean,
          meanB: gB.stats.mean,
          meanDiff,
          seDiff,
          qStatistic: q,
          tStatistic: tEquivalent,
          pValue: pAdjusted,
          pValueRaw: pPair,
          ci95,
          cohensD,
          isSignificant: pAdjusted < 0.05
        });
      }
    }

    return {
      testName: 'One-Way Analysis of Variance (ANOVA)',
      k,
      totalN,
      grandMean,
      groups: processedGroups,
      dfBetween,
      dfWithin,
      dfTotal,
      ssBetween,
      ssWithin,
      ssTotal,
      msBetween,
      msWithin,
      fStatistic: F,
      pValue,
      etaSquared,
      omegaSquared: Math.max(0, omegaSquared),
      pairwise,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Kruskal-Wallis H Test (Non-parametric ANOVA)
   */
  kruskalWallis(groups) {
    if (!Array.isArray(groups) || groups.length < 2) {
      return { error: 'Kruskal-Wallis test requires at least 2 groups.' };
    }

    const allData = [];
    const groupMap = [];

    groups.forEach((g, gIdx) => {
      const clean = Descriptive.cleanData(g.data);
      clean.forEach(val => {
        allData.push({ val, groupIdx: gIdx, name: g.name || `Group ${gIdx + 1}` });
      });
      groupMap.push({ name: g.name || `Group ${gIdx + 1}`, n: clean.length, rankSum: 0 });
    });

    const N = allData.length;
    const k = groups.length;
    if (N < 3) {
      return { error: 'Insufficient total sample size for Kruskal-Wallis test.' };
    }

    allData.sort((a, b) => a.val - b.val);

    // Rank with ties
    let i = 0;
    while (i < N) {
      let j = i;
      while (j < N - 1 && allData[j + 1].val === allData[i].val) {
        j++;
      }
      const rank = (i + 1 + j + 1) / 2;
      for (let m = i; m <= j; m++) {
        allData[m].rank = rank;
      }
      i = j + 1;
    }

    for (const item of allData) {
      groupMap[item.groupIdx].rankSum += item.rank;
    }

    let sumRankSqOverN = 0;
    for (const gm of groupMap) {
      if (gm.n > 0) {
        sumRankSqOverN += (gm.rankSum * gm.rankSum) / gm.n;
      }
    }

    const H = (12 / (N * (N + 1))) * sumRankSqOverN - 3 * (N + 1);
    const df = k - 1;
    const pValue = Distributions.chiSquarePValue(H, df);

    // Epsilon squared effect size
    const epsSq = (H) / ((N * N - 1) / (N + 1));

    return {
      testName: 'Kruskal-Wallis H Test (Non-Parametric ANOVA)',
      k,
      totalN: N,
      df,
      statistic: H,
      pValue,
      epsilonSquared: Math.min(1.0, Math.max(0, epsSq)),
      groups: groupMap,
      isSignificant: pValue < 0.05
    };
  }
};
