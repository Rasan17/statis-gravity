/**
 * Statis-Gravity - ANOVA & Multi-Group Comparison Module
 * One-Way ANOVA, Welch's ANOVA, Kruskal-Wallis H Test, Repeated Measures ANOVA, Friedman Test,
 * Levene's Test of Homoscedasticity, Automated Assumption Diagnostics & Decision Engine.
 */

import { Distributions } from './distributions.js';
import { Descriptive } from './descriptive.js';

export const Anova = {
  /**
   * One-Way Analysis of Variance (Fisher's Standard ANOVA)
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
          statisticLabel: 'q / t',
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
      testKey: 'oneway',
      testName: 'One-Way Analysis of Variance (Fisher\'s ANOVA)',
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
      statistic: F,
      fStatistic: F,
      pValue,
      etaSquared,
      omegaSquared: Math.max(0, omegaSquared),
      effectSizeLabel: 'Eta² (η²) & Omega² (ω²)',
      pairwise,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Welch's ANOVA (Robust ANOVA for Heteroscedastic / Unequal Variances)
   * With Games-Howell Post-Hoc Pairwise Contrasts
   * @param {Array<{name: string, data: number[]}>} groups
   */
  welch(groups) {
    if (!Array.isArray(groups) || groups.length < 2) {
      return { error: 'Welch ANOVA requires at least 2 distinct groups.' };
    }

    const processedGroups = groups.map(g => {
      const stats = Descriptive.calculate(g.data);
      return {
        name: g.name || 'Group',
        stats,
        data: stats.values
      };
    }).filter(g => g.stats.n > 1);

    const k = processedGroups.length;
    if (k < 2) {
      return { error: 'Welch ANOVA requires at least 2 groups with n ≥ 2.' };
    }

    const totalN = processedGroups.reduce((acc, g) => acc + g.stats.n, 0);

    // Weights: w_j = n_j / s_j^2
    const weights = [];
    let sumW = 0;
    for (const g of processedGroups) {
      const varG = Math.max(1e-9, g.stats.variance);
      const w = g.stats.n / varG;
      weights.push(w);
      sumW += w;
    }

    // Weighted Grand Mean: x_prime = sum(w_j * mean_j) / sum(w_j)
    let weightedSumMean = 0;
    for (let j = 0; j < k; j++) {
      weightedSumMean += weights[j] * processedGroups[j].stats.mean;
    }
    const weightedGrandMean = weightedSumMean / sumW;

    // Numerator A = sum(w_j * (mean_j - weightedGrandMean)^2) / (k - 1)
    let sumWeightDevSq = 0;
    for (let j = 0; j < k; j++) {
      sumWeightDevSq += weights[j] * Math.pow(processedGroups[j].stats.mean - weightedGrandMean, 2);
    }
    const numeratorA = sumWeightDevSq / (k - 1);

    // Lambda Term B = sum( (1 - w_j / sumW)^2 / (n_j - 1) )
    let lambdaTermB = 0;
    for (let j = 0; j < k; j++) {
      const g = processedGroups[j];
      const wRatio = 1 - weights[j] / sumW;
      lambdaTermB += Math.pow(wRatio, 2) / (g.stats.n - 1);
    }

    // Denominator adjustment D = 1 + (2 * (k - 2) / (k^2 - 1)) * lambdaTermB
    const denomAdjustment = (k > 1 && (k * k - 1) > 0)
      ? 1 + ((2 * (k - 2)) / (k * k - 1)) * lambdaTermB
      : 1;

    const F_welch = denomAdjustment === 0 ? 0 : numeratorA / denomAdjustment;
    const df1 = k - 1;
    const df2 = lambdaTermB > 0 ? (k * k - 1) / (3 * lambdaTermB) : (totalN - k);
    const pValue = Distributions.fPValue(F_welch, df1, df2);

    // Estimated Omega-squared
    const omegaSquared = Math.max(0, (df1 * (F_welch - 1)) / (df1 * (F_welch - 1) + totalN));
    // Eta-squared approximation
    const etaSquared = Math.max(0, (df1 * F_welch) / (df1 * F_welch + df2));

    // Games-Howell Post-Hoc Pairwise Contrasts (Robust to unequal variances & sample sizes)
    const pairwise = [];
    const numPairs = (k * (k - 1)) / 2;

    for (let i = 0; i < k; i++) {
      for (let j = i + 1; j < k; j++) {
        const gA = processedGroups[i];
        const gB = processedGroups[j];
        const meanDiff = gA.stats.mean - gB.stats.mean;
        const varAOverN = Math.max(1e-9, gA.stats.variance) / gA.stats.n;
        const varBOverN = Math.max(1e-9, gB.stats.variance) / gB.stats.n;
        const seDiff = Math.sqrt(varAOverN + varBOverN);

        const t = seDiff === 0 ? 0 : Math.abs(meanDiff) / seDiff;

        // Welch-Satterthwaite df for this pair
        const dfNumerator = Math.pow(varAOverN + varBOverN, 2);
        const dfDenominator = Math.pow(varAOverN, 2) / (gA.stats.n - 1) + Math.pow(varBOverN, 2) / (gB.stats.n - 1);
        const dfPair = dfDenominator === 0 ? 1 : dfNumerator / dfDenominator;

        // Studentized range q = sqrt(2) * t
        const q = Math.SQRT2 * t;
        const pPairRaw = Distributions.tPValue(t, dfPair);
        const pAdjusted = Math.min(1.0, pPairRaw * numPairs);

        const z95 = 1.95996;
        const tCrit = dfPair > 30 ? z95 : z95 * (1 + 1 / (4 * dfPair));
        const margin = tCrit * seDiff;
        const ci95 = [meanDiff - margin, meanDiff + margin];

        const pooledVar = (gA.stats.variance + gB.stats.variance) / 2;
        const cohensD = pooledVar <= 0 ? 0 : meanDiff / Math.sqrt(pooledVar);

        pairwise.push({
          groupA: gA.name,
          groupB: gB.name,
          comparison: `${gA.name} vs ${gB.name}`,
          meanA: gA.stats.mean,
          meanB: gB.stats.mean,
          meanDiff,
          seDiff,
          statisticLabel: 'Games-Howell t (df\')',
          qStatistic: q,
          tStatistic: t,
          df: dfPair,
          pValue: pAdjusted,
          pValueRaw: pPairRaw,
          ci95,
          cohensD,
          isSignificant: pAdjusted < 0.05
        });
      }
    }

    return {
      testKey: 'welch',
      testName: "Welch's Heteroscedastic ANOVA (Robust)",
      k,
      totalN,
      grandMean: weightedGrandMean,
      groups: processedGroups,
      dfBetween: df1,
      dfWithin: df2,
      df1,
      df2,
      statistic: F_welch,
      fStatistic: F_welch,
      pValue,
      etaSquared,
      omegaSquared,
      effectSizeLabel: 'Robust Omega² (ω²)',
      pairwise,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Kruskal-Wallis H Test (One-Way ANOVA by Ranks for Independent Groups)
   * With Dunn's Post-Hoc Pairwise Contrasts
   * @param {Array<{name: string, data: number[]}>} groups
   */
  kruskalWallis(groups) {
    if (!Array.isArray(groups) || groups.length < 2) {
      return { error: 'Kruskal-Wallis test requires at least 2 groups.' };
    }

    const processedGroups = groups.map((g, idx) => {
      const clean = Descriptive.cleanData(g.data);
      const stats = Descriptive.calculate(clean);
      return {
        name: g.name || `Cohort ${idx + 1}`,
        stats,
        data: clean
      };
    }).filter(g => g.stats.n > 0);

    const k = processedGroups.length;
    if (k < 2) return { error: 'At least 2 groups must contain valid observations.' };

    const allData = [];
    processedGroups.forEach((g, gIdx) => {
      g.data.forEach(val => {
        allData.push({ val, groupIdx: gIdx, name: g.name });
      });
    });

    const N = allData.length;
    if (N < 3) return { error: 'Insufficient total sample size for Kruskal-Wallis test.' };

    allData.sort((a, b) => a.val - b.val);

    // Rank all observations with average ties
    let i = 0;
    const tieCounts = [];
    while (i < N) {
      let j = i;
      while (j < N - 1 && allData[j + 1].val === allData[i].val) {
        j++;
      }
      const tieSize = j - i + 1;
      if (tieSize > 1) tieCounts.push(tieSize);

      const rank = (i + 1 + j + 1) / 2;
      for (let m = i; m <= j; m++) {
        allData[m].rank = rank;
      }
      i = j + 1;
    }

    // Accumulate ranks per group
    const groupRankSums = new Array(k).fill(0);
    for (const item of allData) {
      groupRankSums[item.groupIdx] += item.rank;
    }

    let sumRankSqOverN = 0;
    for (let j = 0; j < k; j++) {
      const n_j = processedGroups[j].stats.n;
      sumRankSqOverN += Math.pow(groupRankSums[j], 2) / n_j;
    }

    let H = (12 / (N * (N + 1))) * sumRankSqOverN - 3 * (N + 1);

    // Tie correction factor
    let tieCorrectionSum = 0;
    for (const t of tieCounts) {
      tieCorrectionSum += (Math.pow(t, 3) - t);
    }
    const tieFactor = 1 - tieCorrectionSum / (Math.pow(N, 3) - N);
    if (tieFactor > 0 && tieFactor < 1) {
      H = H / tieFactor;
    }

    const df = k - 1;
    const pValue = Distributions.chiSquarePValue(H, df);

    // Epsilon-squared effect size: ε² = H / ( (N^2 - 1) / (N + 1) ) = H / (N - 1)
    const epsilonSquared = N > 1 ? Math.min(1.0, Math.max(0, H / (N - 1))) : 0;

    // Dunn's Post-Hoc Pairwise Contrasts (with Bonferroni adjustment)
    const pairwise = [];
    const numPairs = (k * (k - 1)) / 2;
    const tieAdjVariance = (N * (N + 1) / 12) - (tieCorrectionSum / (12 * (N - 1)));

    for (let a = 0; a < k; a++) {
      for (let b = a + 1; b < k; b++) {
        const gA = processedGroups[a];
        const gB = processedGroups[b];
        const meanRankA = groupRankSums[a] / gA.stats.n;
        const meanRankB = groupRankSums[b] / gB.stats.n;
        const rankDiff = meanRankA - meanRankB;
        const seDunn = Math.sqrt(Math.max(1e-9, tieAdjVariance) * (1 / gA.stats.n + 1 / gB.stats.n));
        const z = seDunn === 0 ? 0 : Math.abs(rankDiff) / seDunn;
        const pPairRaw = Distributions.normalPValue(z);
        const pAdjusted = Math.min(1.0, pPairRaw * numPairs);

        // Rank-biserial / effect size r = z / sqrt(N)
        const rEffect = z / Math.sqrt(N);
        const medianDiff = gA.stats.median - gB.stats.median;

        pairwise.push({
          groupA: gA.name,
          groupB: gB.name,
          comparison: `${gA.name} vs ${gB.name}`,
          meanA: gA.stats.mean,
          meanB: gB.stats.mean,
          medianA: gA.stats.median,
          medianB: gB.stats.median,
          meanDiff: medianDiff, // Primary non-parametric difference
          rankDiff,
          seDiff: seDunn,
          statisticLabel: "Dunn's z",
          zStatistic: z,
          tStatistic: z,
          pValue: pAdjusted,
          pValueRaw: pPairRaw,
          ci95: [gA.stats.median - gB.stats.median, gA.stats.median - gB.stats.median],
          cohensD: rEffect,
          isSignificant: pAdjusted < 0.05
        });
      }
    }

    return {
      testKey: 'kruskal',
      testName: 'Kruskal-Wallis H Test (Non-Parametric ANOVA)',
      k,
      totalN: N,
      dfBetween: df,
      dfWithin: N - k,
      df,
      statistic: H,
      fStatistic: H, // for uniform metric display
      pValue,
      etaSquared: epsilonSquared,
      epsilonSquared,
      omegaSquared: epsilonSquared,
      effectSizeLabel: 'Epsilon-Squared (ε²)',
      groups: processedGroups,
      groupRankSums,
      pairwise,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * One-Way Repeated Measures ANOVA (Paired / Within-Subjects Parametric)
   * With Pairwise Bonferroni-Adjusted Paired t-tests & Greenhouse-Geisser Epsilon
   * @param {Array<{name: string, data: number[]}>} groups
   */
  repeatedMeasures(groups) {
    if (!Array.isArray(groups) || groups.length < 2) {
      return { error: 'Repeated Measures ANOVA requires at least 2 conditions/timepoints.' };
    }

    const k = groups.length;
    const cleaned = groups.map((g, idx) => ({
      name: g.name || `Timepoint ${idx + 1}`,
      data: Descriptive.cleanData(g.data)
    }));

    // Find matched sample size N (truncate to minimum length if slight discrepancy)
    const minN = Math.min(...cleaned.map(g => g.data.length));
    if (minN < 2) {
      return { error: 'Repeated Measures ANOVA requires at least 2 complete subjects across all timepoints.' };
    }

    const N = minN;
    const processedGroups = cleaned.map(g => {
      const sliced = g.data.slice(0, N);
      const stats = Descriptive.calculate(sliced);
      return {
        name: g.name,
        stats,
        data: sliced
      };
    });

    // Subject means and Condition means
    const subjectSums = new Array(N).fill(0);
    const conditionSums = new Array(k).fill(0);
    let grandSum = 0;

    for (let i = 0; i < N; i++) {
      for (let j = 0; j < k; j++) {
        const val = processedGroups[j].data[i];
        subjectSums[i] += val;
        conditionSums[j] += val;
        grandSum += val;
      }
    }

    const totalObs = N * k;
    const grandMean = grandSum / totalObs;
    const subjectMeans = subjectSums.map(s => s / k);
    const conditionMeans = conditionSums.map(c => c / N);

    // Sum of Squares
    let ssTotal = 0;
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < k; j++) {
        ssTotal += Math.pow(processedGroups[j].data[i] - grandMean, 2);
      }
    }

    let ssSubjects = 0;
    for (let i = 0; i < N; i++) {
      ssSubjects += k * Math.pow(subjectMeans[i] - grandMean, 2);
    }

    const ssWithin = ssTotal - ssSubjects;

    let ssTreatment = 0;
    for (let j = 0; j < k; j++) {
      ssTreatment += N * Math.pow(conditionMeans[j] - grandMean, 2);
    }

    const ssError = Math.max(0, ssWithin - ssTreatment);

    const dfTreatment = k - 1;
    const dfSubjects = N - 1;
    const dfError = (N - 1) * (k - 1);
    const dfTotal = totalObs - 1;

    const msTreatment = dfTreatment > 0 ? ssTreatment / dfTreatment : 0;
    const msError = dfError > 0 ? ssError / dfError : 1e-9;

    const F = msError === 0 ? 0 : msTreatment / msError;
    const pValue = Distributions.fPValue(F, dfTreatment, dfError);

    // Partial Eta-Squared: η²_p = SS_treatment / (SS_treatment + SS_error)
    const partialEtaSquared = (ssTreatment + ssError) > 0 ? ssTreatment / (ssTreatment + ssError) : 0;

    // Greenhouse-Geisser Epsilon estimate approximation
    const ggEpsilon = Math.max(1 / (k - 1), Math.min(1.0, 1 - 0.5 * (k - 1) / (dfError || 1)));
    const pValueGG = Distributions.fPValue(F, dfTreatment * ggEpsilon, dfError * ggEpsilon);

    // Pairwise Bonferroni-Adjusted Paired t-tests
    const pairwise = [];
    const numPairs = (k * (k - 1)) / 2;

    for (let i = 0; i < k; i++) {
      for (let j = i + 1; j < k; j++) {
        const gA = processedGroups[i];
        const gB = processedGroups[j];
        const diffs = [];
        for (let m = 0; m < N; m++) {
          diffs.push(gA.data[m] - gB.data[m]);
        }
        const diffStats = Descriptive.calculate(diffs);
        const meanDiff = diffStats.mean;
        const seDiff = diffStats.sem;
        const t = seDiff === 0 ? 0 : meanDiff / seDiff;
        const dfPair = N - 1;
        const pPairRaw = Distributions.tPValue(Math.abs(t), dfPair);
        const pAdjusted = Math.min(1.0, pPairRaw * numPairs);

        const z95 = 1.95996;
        const tCrit = dfPair > 30 ? z95 : z95 * (1 + 1 / (4 * dfPair));
        const margin = tCrit * seDiff;
        const ci95 = [meanDiff - margin, meanDiff + margin];
        const cohensD = diffStats.sd === 0 ? 0 : meanDiff / diffStats.sd;

        pairwise.push({
          groupA: gA.name,
          groupB: gB.name,
          comparison: `${gA.name} vs ${gB.name}`,
          meanA: gA.stats.mean,
          meanB: gB.stats.mean,
          meanDiff,
          seDiff,
          statisticLabel: 'Paired t (df)',
          tStatistic: t,
          df: dfPair,
          pValue: pAdjusted,
          pValueRaw: pPairRaw,
          ci95,
          cohensD,
          isSignificant: pAdjusted < 0.05
        });
      }
    }

    return {
      testKey: 'rm_anova',
      testName: 'One-Way Repeated Measures ANOVA (Within-Subjects)',
      k,
      totalN: N,
      matchedN: N,
      grandMean,
      groups: processedGroups,
      dfBetween: dfTreatment,
      dfWithin: dfError,
      dfTreatment,
      dfSubjects,
      dfError,
      dfTotal,
      ssBetween: ssTreatment,
      ssWithin: ssError,
      ssTreatment,
      ssSubjects,
      ssError,
      ssTotal,
      msTreatment,
      msError,
      statistic: F,
      fStatistic: F,
      pValue,
      pValueGG,
      ggEpsilon,
      etaSquared: partialEtaSquared,
      omegaSquared: partialEtaSquared,
      partialEtaSquared,
      effectSizeLabel: 'Partial Eta² (η²_p)',
      pairwise,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Friedman Test (Non-Parametric Two-Way ANOVA by Ranks for Repeated Measures)
   * With Pairwise Wilcoxon Signed-Rank Contrasts & Kendall's W Effect Size
   * @param {Array<{name: string, data: number[]}>} groups
   */
  friedman(groups) {
    if (!Array.isArray(groups) || groups.length < 2) {
      return { error: 'Friedman test requires at least 2 conditions/timepoints.' };
    }

    const k = groups.length;
    const cleaned = groups.map((g, idx) => ({
      name: g.name || `Condition ${idx + 1}`,
      data: Descriptive.cleanData(g.data)
    }));

    const minN = Math.min(...cleaned.map(g => g.data.length));
    if (minN < 2) {
      return { error: 'Friedman test requires at least 2 complete subjects across all conditions.' };
    }

    const N = minN;
    const processedGroups = cleaned.map(g => {
      const sliced = g.data.slice(0, N);
      const stats = Descriptive.calculate(sliced);
      return {
        name: g.name,
        stats,
        data: sliced
      };
    });

    // Rank across conditions for each subject
    const rankMatrix = []; // [subject][condition]
    let totalTieCorrection = 0;

    for (let i = 0; i < N; i++) {
      const row = [];
      for (let j = 0; j < k; j++) {
        row.push({ condIdx: j, val: processedGroups[j].data[i] });
      }
      row.sort((a, b) => a.val - b.val);

      // Rank with ties
      let p = 0;
      while (p < k) {
        let q = p;
        while (q < k - 1 && row[q + 1].val === row[p].val) q++;
        const tieLen = q - p + 1;
        if (tieLen > 1) {
          totalTieCorrection += Math.pow(tieLen, 3) - tieLen;
        }
        const rank = (p + 1 + q + 1) / 2;
        for (let m = p; m <= q; m++) {
          row[m].rank = rank;
        }
        p = q + 1;
      }

      const assignedRanks = new Array(k);
      for (const item of row) {
        assignedRanks[item.condIdx] = item.rank;
      }
      rankMatrix.push(assignedRanks);
    }

    // Sum of ranks per condition
    const conditionRankSums = new Array(k).fill(0);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < k; j++) {
        conditionRankSums[j] += rankMatrix[i][j];
      }
    }

    let sumRankSq = 0;
    for (let j = 0; j < k; j++) {
      sumRankSq += Math.pow(conditionRankSums[j], 2);
    }

    let chiSqF = (12 / (N * k * (k + 1))) * sumRankSq - 3 * N * (k + 1);

    // Tie correction factor: C = 1 - totalTieCorrection / (N * k * (k^2 - 1))
    if (k > 1) {
      const tieDenom = N * k * (k * k - 1);
      if (tieDenom > 0 && totalTieCorrection > 0) {
        const C = 1 - totalTieCorrection / tieDenom;
        if (C > 0 && C < 1) {
          chiSqF = chiSqF / C;
        }
      }
    }

    const df = k - 1;
    const pValue = Distributions.chiSquarePValue(chiSqF, df);

    // Kendall's W Concordance Coefficient: W = chiSqF / (N * (k - 1))
    const kendallsW = (N > 0 && df > 0) ? Math.min(1.0, Math.max(0, chiSqF / (N * df))) : 0;

    // Pairwise Wilcoxon Signed-Rank Contrasts (with Bonferroni adjustment)
    const pairwise = [];
    const numPairs = (k * (k - 1)) / 2;

    for (let a = 0; a < k; a++) {
      for (let b = a + 1; b < k; b++) {
        const gA = processedGroups[a];
        const gB = processedGroups[b];
        const diffs = [];
        for (let m = 0; m < N; m++) {
          const d = gA.data[m] - gB.data[m];
          if (d !== 0) diffs.push({ diff: d, absDiff: Math.abs(d) });
        }

        let z = 0;
        let pPairRaw = 1.0;
        const nDiff = diffs.length;

        if (nDiff >= 2) {
          diffs.sort((x, y) => x.absDiff - y.absDiff);
          let p = 0;
          while (p < nDiff) {
            let q = p;
            while (q < nDiff - 1 && diffs[q + 1].absDiff === diffs[p].absDiff) q++;
            const rank = (p + 1 + q + 1) / 2;
            for (let m = p; m <= q; m++) diffs[m].rank = rank;
            p = q + 1;
          }

          let wPlus = 0;
          let wMinus = 0;
          for (const d of diffs) {
            if (d.diff > 0) wPlus += d.rank;
            else wMinus += d.rank;
          }
          const W_stat = Math.min(wPlus, wMinus);
          const meanW = (nDiff * (nDiff + 1)) / 4;
          const sigmaW = Math.sqrt((nDiff * (nDiff + 1) * (2 * nDiff + 1)) / 24);
          z = sigmaW === 0 ? 0 : (W_stat - meanW) / sigmaW;
          pPairRaw = Distributions.normalPValue(Math.abs(z));
        }

        const pAdjusted = Math.min(1.0, pPairRaw * numPairs);
        const medianDiff = gA.stats.median - gB.stats.median;

        pairwise.push({
          groupA: gA.name,
          groupB: gB.name,
          comparison: `${gA.name} vs ${gB.name}`,
          meanA: gA.stats.mean,
          meanB: gB.stats.mean,
          medianA: gA.stats.median,
          medianB: gB.stats.median,
          meanDiff: medianDiff,
          seDiff: 0,
          statisticLabel: 'Wilcoxon z',
          zStatistic: z,
          tStatistic: z,
          pValue: pAdjusted,
          pValueRaw: pPairRaw,
          ci95: [medianDiff, medianDiff],
          cohensD: Math.abs(z) / Math.sqrt(N),
          isSignificant: pAdjusted < 0.05
        });
      }
    }

    return {
      testKey: 'friedman',
      testName: 'Friedman Test (Non-Parametric Repeated Measures)',
      k,
      totalN: N,
      matchedN: N,
      dfBetween: df,
      dfWithin: N * (k - 1),
      df,
      statistic: chiSqF,
      fStatistic: chiSqF,
      pValue,
      etaSquared: kendallsW,
      omegaSquared: kendallsW,
      kendallsW,
      effectSizeLabel: "Kendall's Concordance W",
      groups: processedGroups,
      conditionRankSums,
      pairwise,
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Levene's Test of Homoscedasticity (Equality of Variances across k groups)
   * Uses Brown-Forsythe median-centered absolute deviations
   * @param {Array<{name: string, data: number[]}>} groups
   */
  leveneTest(groups) {
    if (!Array.isArray(groups) || groups.length < 2) {
      return { error: 'Levene test requires at least 2 groups.' };
    }

    const processed = groups.map(g => {
      const clean = Descriptive.cleanData(g.data);
      const stats = Descriptive.calculate(clean);
      return { name: g.name, data: clean, stats };
    }).filter(g => g.data.length > 1);

    const k = processed.length;
    if (k < 2) return { error: 'Levene test requires at least 2 groups with n ≥ 2.' };

    const totalN = processed.reduce((acc, g) => acc + g.data.length, 0);

    // Compute absolute deviations from group median z_ij = |x_ij - median_j|
    const deviationGroups = processed.map(g => {
      const med = g.stats.median;
      const devData = g.data.map(x => Math.abs(x - med));
      return {
        name: g.name,
        data: devData
      };
    });

    // Run One-Way ANOVA on the deviations
    const anovaOnDeviations = this.oneWay(deviationGroups);
    if (anovaOnDeviations.error) {
      return { error: anovaOnDeviations.error };
    }

    const F = anovaOnDeviations.fStatistic;
    const df1 = anovaOnDeviations.dfBetween;
    const df2 = anovaOnDeviations.dfWithin;
    const pValue = anovaOnDeviations.pValue;

    // Variance ratio max(var) / min(var)
    const variances = processed.map(g => g.stats.variance).filter(v => v > 0);
    const maxVar = variances.length > 0 ? Math.max(...variances) : 1;
    const minVar = variances.length > 0 ? Math.max(1e-9, Math.min(...variances)) : 1;
    const varianceRatio = maxVar / minVar;

    return {
      statistic: F,
      fStatistic: F,
      df1,
      df2,
      pValue,
      varianceRatio,
      equalVariance: pValue >= 0.05,
      interpretation: pValue >= 0.05
        ? 'Homoscedastic (Equal Variances Confirmed, p ≥ .05)'
        : 'Heteroscedastic (Unequal Variances Detected, p < .05)'
    };
  },

  /**
   * Diagnostic Assumption Evaluation Engine for Multi-Group Data
   * Assesses Pairing, Normality across all k cohorts, Homoscedasticity, Sphericity
   * Recommends optimal statistical test & provides clinical decision rationale
   * @param {Array<{name: string, data: number[]}>} groups
   * @param {boolean} isPaired
   */
  evaluateAssumptions(groups, isPaired = false) {
    if (!Array.isArray(groups) || groups.length < 2) {
      return { error: 'Assumption evaluation requires at least 2 cohorts.' };
    }

    const cleaned = groups.map((g, idx) => ({
      name: g.name || `Cohort ${idx + 1}`,
      data: Descriptive.cleanData(g.data)
    })).filter(g => g.data.length > 0);

    const k = cleaned.length;
    if (k < 2) return { error: 'At least 2 cohorts must have numerical observations.' };

    const cohortStats = cleaned.map(g => Descriptive.calculate(g.data));

    // 1. Normality Assessment across all k cohorts
    const normalityDetails = [];
    let allNormal = true;

    for (let j = 0; j < k; j++) {
      const stats = cohortStats[j];
      const jb = stats.normality || { isNormal: true, pValue: 1.0, statistic: 0 };
      const skewAlert = Math.abs(stats.skewness) > 1.0;
      const kurtAlert = Math.abs(stats.kurtosis) > 1.5;
      const isNorm = stats.n >= 10 ? jb.isNormal : (!skewAlert && !kurtAlert);

      if (!isNorm) allNormal = false;

      normalityDetails.push({
        cohortIndex: j + 1,
        name: cleaned[j].name,
        n: stats.n,
        mean: stats.mean,
        median: stats.median,
        sd: stats.sd,
        skewness: stats.skewness,
        kurtosis: stats.kurtosis,
        jbStat: jb.statistic,
        pValue: jb.pValue,
        isNormal: isNorm
      });
    }

    if (isPaired) {
      // PAIRED / REPEATED MEASURES DESIGN
      const lengths = cleaned.map(g => g.data.length);
      const minN = Math.min(...lengths);
      const maxN = Math.max(...lengths);
      const unequalLengths = minN !== maxN;

      if (minN < 2) {
        return {
          error: 'Paired repeated measures design requires at least 2 matched subjects across all cohorts.',
          isPaired: true
        };
      }

      // Check normality of within-subject differences from baseline
      let pairedDiffsNormal = true;
      const baselineData = cleaned[0].data.slice(0, minN);
      for (let j = 1; j < k; j++) {
        const compData = cleaned[j].data.slice(0, minN);
        const diffs = compData.map((val, idx) => val - baselineData[idx]);
        const diffStats = Descriptive.calculate(diffs);
        const jbDiff = diffStats.normality || { isNormal: true, pValue: 1.0, statistic: 0 };
        const skewAlert = Math.abs(diffStats.skewness) > 1.0;
        const kurtAlert = Math.abs(diffStats.kurtosis) > 1.5;
        const diffNorm = minN >= 10 ? jbDiff.isNormal : (!skewAlert && !kurtAlert);
        if (!diffNorm) pairedDiffsNormal = false;
      }

      const isParametric = allNormal || pairedDiffsNormal;
      const recommendedTest = isParametric ? 'rm_anova' : 'friedman';
      const recommendedTestName = isParametric
        ? 'Repeated Measures ANOVA (Paired, Parametric)'
        : 'Friedman Test (Paired, Non-Parametric)';

      let rationale = `Within-subjects repeated measures design with ${k} matched conditions across N = ${minN} participants. `;
      if (unequalLengths) {
        rationale += `Note: Cohort observation counts varied (${lengths.join(', ')}); data was strictly aligned by subject index to the matched sample size N = ${minN}. `;
      }

      if (isParametric) {
        rationale += `Condition measurements and within-subject difference scores conform adequately to normality without severe skewness. One-Way Repeated Measures ANOVA is recommended.`;
      } else {
        rationale += `One or more conditions exhibit significant departure from normality or heavy skewness. Non-parametric ranking via the Friedman test is recommended to prevent Type I error distortion.`;
      }

      return {
        isPaired: true,
        k,
        matchedN: minN,
        unequalLengths,
        lengths,
        cohortStats,
        normality: {
          isParametric,
          allNormal,
          details: normalityDetails,
          interpretation: isParametric ? 'Normal / Symmetrical Distribution' : 'Non-Normal Distribution Detected'
        },
        varianceEquality: {
          applicable: false,
          note: 'Sphericity applies to repeated measures; homoscedasticity evaluated between subjects is not required for purely within-subject contrasts.'
        },
        recommendedTest,
        recommendedTestName,
        rationale
      };
    } else {
      // INDEPENDENT COHORTS DESIGN
      // 2. Homogeneity of Variances via Levene's Test
      const levene = this.leveneTest(cleaned);
      const equalVariance = levene.equalVariance !== undefined ? levene.equalVariance : true;
      const isParametric = allNormal;

      let recommendedTest;
      let recommendedTestName;
      let rationale = `Independent multi-cohort comparison across ${k} distinct groups (total N = ${cohortStats.reduce((a, s) => a + s.n, 0)}). `;

      if (isParametric) {
        if (equalVariance) {
          recommendedTest = 'oneway';
          recommendedTestName = "One-Way ANOVA (Fisher's Standard - Equal Variances)";
          rationale += `All ${k} cohorts conform to normal distributions (Jarque-Bera p ≥ .05) and Levene's test confirms homogeneity of variances (F(${levene.df1}, ${levene.df2}) = ${(levene.fStatistic || 0).toFixed(2)}, p = ${levene.pValue > 0.001 ? levene.pValue.toFixed(3) : '< .001'}). Standard Fisher's One-Way ANOVA with Tukey's HSD post-hoc contrasts is optimal.`;
        } else {
          recommendedTest = 'welch';
          recommendedTestName = "Welch's ANOVA (Robust - Unequal Variances)";
          rationale += `All ${k} cohorts conform to normal distributions, but Levene's test detected significant heteroscedasticity / unequal variances (F(${levene.df1}, ${levene.df2}) = ${(levene.fStatistic || 0).toFixed(2)}, p = ${levene.pValue > 0.001 ? levene.pValue.toFixed(3) : '< .001'}, variance ratio = ${(levene.varianceRatio || 1).toFixed(2)}×). Standard Fisher ANOVA inflates Type I errors under heteroscedasticity; Welch's ANOVA with Games-Howell post-hoc contrasts is strongly recommended.`;
        }
      } else {
        recommendedTest = 'kruskal';
        recommendedTestName = 'Kruskal-Wallis H Test (Non-Parametric)';
        const nonNormalCohorts = normalityDetails.filter(d => !d.isNormal).map(d => `${d.name} (p = ${d.pValue > 0.001 ? d.pValue.toFixed(3) : '< .001'})`);
        rationale += `Departure from normality detected in cohort(s): ${nonNormalCohorts.join(', ')}. The non-parametric Kruskal-Wallis H test by ranks (with Dunn's post-hoc contrasts) is recommended to protect against distribution anomalies.`;
      }

      return {
        isPaired: false,
        k,
        cohortStats,
        normality: {
          isParametric,
          allNormal,
          details: normalityDetails,
          interpretation: isParametric ? 'All Cohorts Normal (Parametric Suitable)' : 'Normality Violated (Non-Parametric Recommended)'
        },
        varianceEquality: {
          applicable: true,
          equalVariance,
          fStat: levene.fStatistic || 0,
          df1: levene.df1 || 0,
          df2: levene.df2 || 0,
          pValue: levene.pValue || 1.0,
          varianceRatio: levene.varianceRatio || 1.0,
          interpretation: levene.interpretation || (equalVariance ? 'Homoscedastic (Equal Variances)' : 'Heteroscedastic (Unequal Variances)')
        },
        recommendedTest,
        recommendedTestName,
        rationale
      };
    }
  },

  /**
   * Main Dispatcher for Multi-Group Analysis
   * @param {Array<{name: string, data: number[]}>} groups
   * @param {string} testType 'auto' | 'oneway' | 'welch' | 'kruskal' | 'rm_anova' | 'friedman'
   * @param {boolean} isPaired
   */
  test(groups, testType = 'auto', isPaired = false) {
    const assumptions = this.evaluateAssumptions(groups, isPaired);
    if (assumptions.error) {
      return { error: assumptions.error };
    }

    let activeTest = testType;
    if (activeTest === 'auto') {
      activeTest = assumptions.recommendedTest;
    }

    let result;
    switch (activeTest) {
      case 'welch':
        result = this.welch(groups);
        break;
      case 'kruskal':
        result = this.kruskalWallis(groups);
        break;
      case 'rm_anova':
        result = this.repeatedMeasures(groups);
        break;
      case 'friedman':
        result = this.friedman(groups);
        break;
      case 'oneway':
      default:
        result = this.oneWay(groups);
        break;
    }

    if (result.error) {
      return { error: result.error, assumptions };
    }

    result.assumptions = assumptions;
    result.requestedMethod = testType;
    result.executedMethod = activeTest;
    result.isPaired = isPaired;
    return result;
  }
};
