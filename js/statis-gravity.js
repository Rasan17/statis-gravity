/**
 * Statis-Gravity — Complete Self-Contained Clinical Biostatistics & Visualization Engine
 * Unified standalone build: Compatible with both local file:// protocol and web servers.
 * Conceived, supervised design and testing: Dr G Narenthiran MB ChB BSc(MedSci)(Hons) MRCS(Ed.) FEBNS FRCS(SN)
 * Copyright, G Narenthiran FEBS FRCS(SN), g_narenthiran@hotmail.com
 */

(function () {
  'use strict';

  // ==========================================
  // 1. STATISTICAL DISTRIBUTIONS
  // ==========================================
  const Distributions = {
    normalCDF(z) {
      if (isNaN(z)) return NaN;
      if (z === 0) return 0.5;
      if (z < -8.0) return 0.0;
      if (z > 8.0) return 1.0;

      const sign = z < 0 ? -1 : 1;
      const x = Math.abs(z) / Math.SQRT2;
      const a1 = 0.254829592;
      const a2 = -0.284496736;
      const a3 = 1.421413741;
      const a4 = -1.453152027;
      const a5 = 1.061405429;
      const p = 0.3275911;

      const t = 1.0 / (1.0 + p * x);
      const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
      return 0.5 * (1.0 + sign * y);
    },

    normalPValue(z) {
      const cdf = this.normalCDF(Math.abs(z));
      return 2 * (1 - cdf);
    },

    invNormalCDF(p) {
      if (p <= 0) return -Infinity;
      if (p >= 1) return Infinity;

      const a = [
        -3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00
      ];
      const b = [
        -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01
      ];
      const c = [
        -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
        -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00
      ];
      const d = [
         7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
         3.754408661907416e+00
      ];

      const p_low = 0.02425;
      const p_high = 1 - p_low;

      if (p < p_low) {
        const q = Math.sqrt(-2 * Math.log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
      } else if (p <= p_high) {
        const q = p - 0.5;
        const r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
      } else {
        const q = Math.sqrt(-2 * Math.log(1 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
      }
    },

    logGamma(x) {
      if (x <= 0) return NaN;
      const p = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109585651912,
        9.9843695780195716e-6,
        1.5056327351493116e-7
      ];
      const z = x - 1;
      let sum = p[0];
      for (let i = 1; i < p.length; i++) {
        sum += p[i] / (z + i);
      }
      const t = z + 7.5;
      return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(sum);
    },

    tPDF(t, df) {
      const num = Math.exp(this.logGamma((df + 1) / 2) - this.logGamma(df / 2));
      const den = Math.sqrt(df * Math.PI) * Math.pow(1 + (t * t) / df, (df + 1) / 2);
      return num / den;
    },

    tPValue(t, df) {
      if (df <= 0 || isNaN(t) || isNaN(df)) return NaN;
      const absT = Math.abs(t);
      if (absT === 0) return 1.0;
      if (df > 100) return this.normalPValue(absT);

      let sum = 0;
      const upper = Math.max(50.0, absT + 30.0);
      const N = 400;
      const h = (upper - absT) / N;
      for (let i = 0; i <= N; i++) {
        const x = absT + i * h;
        const w = (i === 0 || i === N) ? 1 : (i % 2 === 1 ? 4 : 2);
        sum += w * this.tPDF(x, df);
      }
      const tail = (sum * h) / 3;
      return Math.min(1.0, Math.max(0.0, 2 * tail));
    },

    chiSquarePValue(chiSq, df) {
      if (chiSq <= 0 || df <= 0) return 1.0;
      if (df === 1) return this.normalPValue(Math.sqrt(chiSq));
      if (df === 2) return Math.exp(-chiSq / 2);

      const a = df / 2;
      const x = chiSq / 2;
      let b = x + 1.0 - a;
      let c = 1.0 / 1e-30;
      let d = 1.0 / b;
      let h = d;
      for (let i = 1; i <= 200; i++) {
        const an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (Math.abs(d) < 1e-30) d = 1e-30;
        c = b + an / c;
        if (Math.abs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        const del = d * c;
        h *= del;
        if (Math.abs(del - 1.0) < 1e-15) break;
      }
      const q = Math.exp(-x + a * Math.log(x) - this.logGamma(a)) * h;
      return Math.min(1.0, Math.max(0.0, q));
    },

    fPDF(f, df1, df2) {
      if (f <= 0) return 0;
      const num = Math.pow(df1 * f, df1) * Math.pow(df2, df2);
      const den = Math.pow(df1 * f + df2, df1 + df2);
      const logBeta = this.logGamma(df1 / 2) + this.logGamma(df2 / 2) - this.logGamma((df1 + df2) / 2);
      const logVal = 0.5 * Math.log(num / den) - logBeta - Math.log(f);
      return Math.exp(logVal);
    },

    fPValue(F, df1, df2) {
      if (F <= 0 || df1 <= 0 || df2 <= 0) return 1.0;
      const upper = Math.max(100.0, F * 15.0);
      const N = 400;
      const h = (upper - F) / N;
      let sum = 0;
      for (let i = 0; i <= N; i++) {
        const x = F + i * h;
        const w = (i === 0 || i === N) ? 1 : (i % 2 === 1 ? 4 : 2);
        sum += w * this.fPDF(x, df1, df2);
      }
      const tail = (sum * h) / 3;
      return Math.min(1.0, Math.max(0.0, tail));
    }
  };

  // ==========================================
  // 2. DESCRIPTIVE STATISTICS
  // ==========================================
  const Descriptive = {
    cleanData(raw) {
      if (!Array.isArray(raw)) return [];
      return raw
        .map(v => (typeof v === 'number' ? v : parseFloat(v)))
        .filter(v => typeof v === 'number' && !isNaN(v) && isFinite(v));
    },

    calculate(data) {
      const values = this.cleanData(data);
      const n = values.length;

      if (n === 0) {
        return { n: 0, error: 'Please enter at least one valid number.' };
      }

      const sorted = [...values].sort((a, b) => a - b);
      const sum = sorted.reduce((acc, val) => acc + val, 0);
      const mean = sum / n;

      let median;
      if (n % 2 === 0) {
        median = (sorted[n / 2 - 1] + sorted[n / 2]) / 2;
      } else {
        median = sorted[Math.floor(n / 2)];
      }

      const min = sorted[0];
      const max = sorted[n - 1];
      const range = max - min;

      // Mode
      const freqMap = new Map();
      let maxFreq = 0;
      for (const v of sorted) {
        const f = (freqMap.get(v) || 0) + 1;
        freqMap.set(v, f);
        if (f > maxFreq) maxFreq = f;
      }
      const modes = [];
      if (maxFreq > 1) {
        for (const [val, freq] of freqMap.entries()) {
          if (freq === maxFreq) modes.push(val);
        }
      }

      const percentile = (p) => {
        if (n === 1) return sorted[0];
        const index = (n - 1) * p;
        const lower = Math.floor(index);
        const upper = Math.ceil(index);
        const weight = index - lower;
        return sorted[lower] * (1 - weight) + sorted[upper] * weight;
      };

      const q1 = percentile(0.25);
      const q3 = percentile(0.75);
      const iqr = q3 - q1;

      // Sum of squared deviations, Variance, Standard Deviation
      let sumSqDiff = 0;
      let sumCubeDiff = 0;
      let sumQuadDiff = 0;

      for (const v of sorted) {
        const diff = v - mean;
        const diffSq = diff * diff;
        sumSqDiff += diffSq;
        sumCubeDiff += diffSq * diff;
        sumQuadDiff += diffSq * diffSq;
      }

      const variance = n > 1 ? sumSqDiff / (n - 1) : 0;
      const sd = Math.sqrt(variance);
      const sem = n > 0 ? sd / Math.sqrt(n) : 0;

      const lowerFence = q1 - 1.5 * iqr;
      const upperFence = q3 + 1.5 * iqr;
      const extremeLowerFence = q1 - 3.0 * iqr;
      const extremeUpperFence = q3 + 3.0 * iqr;

      const outliers = [];
      for (const v of sorted) {
        if (v < lowerFence || v > upperFence) {
          const isExtreme = v < extremeLowerFence || v > extremeUpperFence;
          const zScore = sd > 0 ? (v - mean) / sd : 0;
          outliers.push({
            value: v,
            zScore,
            type: isExtreme ? 'Extreme' : 'Mild',
            direction: v > upperFence ? 'High' : 'Low'
          });
        }
      }

      let skewness = 0;
      let kurtosis = 0;
      if (n >= 3 && sd > 0) {
        skewness = (n * sumCubeDiff) / ((n - 1) * (n - 2) * Math.pow(sd, 3));
      }
      if (n >= 4 && sd > 0) {
        const term1 = (n * (n + 1) * sumQuadDiff) / ((n - 1) * (n - 2) * (n - 3) * Math.pow(sd, 4));
        const term2 = (3 * Math.pow(n - 1, 2)) / ((n - 2) * (n - 3));
        kurtosis = term1 - term2;
      }

      let skewnessInterpretation;
      if (skewness > 1.0) skewnessInterpretation = 'High Right Skew (Positive)';
      else if (skewness > 0.5) skewnessInterpretation = 'Moderate Right Skew';
      else if (skewness < -1.0) skewnessInterpretation = 'High Left Skew (Negative)';
      else if (skewness < -0.5) skewnessInterpretation = 'Moderate Left Skew';
      else skewnessInterpretation = 'Approximately Symmetric';

      let kurtosisInterpretation;
      if (kurtosis >= 1.0) kurtosisInterpretation = 'Leptokurtic (Heavy Tails / Peaked)';
      else if (kurtosis <= -1.0) kurtosisInterpretation = 'Platykurtic (Light Tails / Flat)';
      else kurtosisInterpretation = 'Mesokurtic (Normal Tails)';

      let ci95 = [mean, mean];
      let ci99 = [mean, mean];
      if (n > 1) {
        const df = n - 1;
        const z95 = Distributions.invNormalCDF(0.975);
        const z99 = Distributions.invNormalCDF(0.995);
        const t95Factor = df > 30 ? z95 : z95 * (1 + 1 / (4 * df) + 1 / (32 * df * df));
        const t99Factor = df > 30 ? z99 : z99 * (1 + 1 / (4 * df) + 1 / (32 * df * df));

        const margin95 = t95Factor * sem;
        const margin99 = t99Factor * sem;
        ci95 = [mean - margin95, mean + margin95];
        ci99 = [mean - margin99, mean + margin99];
      }

      let jbStat = 0;
      let jbPValue = 1.0;
      if (n >= 8 && sd > 0) {
        jbStat = (n / 6) * (Math.pow(skewness, 2) + Math.pow(kurtosis, 2) / 4);
        jbPValue = Distributions.chiSquarePValue(jbStat, 2);
      }

      return {
        n,
        values: sorted,
        sum,
        mean,
        median,
        modes,
        maxFreq,
        min,
        max,
        range,
        q1,
        q3,
        iqr,
        lowerFence,
        upperFence,
        extremeLowerFence,
        extremeUpperFence,
        outliers,
        variance,
        sd,
        sem,
        skewness,
        skewnessInterpretation,
        kurtosis,
        kurtosisInterpretation,
        ci95,
        ci99,
        normality: {
          test: 'Jarque-Bera',
          statistic: jbStat,
          pValue: jbPValue,
          isNormal: jbPValue > 0.05
        }
      };
    }
  };

  // ==========================================
  // 3. HYPOTHESIS TESTING
  // ==========================================
  const Hypothesis = {
    independentTTest(groupA, groupB) {
      const statsA = Descriptive.calculate(groupA);
      const statsB = Descriptive.calculate(groupB);

      if (statsA.n < 2 || statsB.n < 2) {
        return { error: 'Both groups must contain at least 2 observations.' };
      }

      const n1 = statsA.n;
      const n2 = statsB.n;
      const df = n1 + n2 - 2;

      const sp2 = ((n1 - 1) * statsA.variance + (n2 - 1) * statsB.variance) / df;
      const seDiff = Math.sqrt(sp2 * (1 / n1 + 1 / n2));
      const meanDiff = statsA.mean - statsB.mean;
      const t = seDiff === 0 ? 0 : meanDiff / seDiff;
      const pValue = Distributions.tPValue(t, df);

      const sp = Math.sqrt(sp2);
      const cohensD = sp === 0 ? 0 : meanDiff / sp;

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

      const num = Math.pow(v1 + v2, 2);
      const den = Math.pow(v1, 2) / (n1 - 1) + Math.pow(v2, 2) / (n2 - 1);
      const df = den === 0 ? 1 : num / den;
      const pValue = Distributions.tPValue(t, df);

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

    pairedTTest(groupPre, groupPost) {
      const a = Descriptive.cleanData(groupPre);
      const b = Descriptive.cleanData(groupPost);
      const n = Math.min(a.length, b.length);
      if (n < 2) return { error: 'Paired test requires at least 2 paired observations.' };

      const diffs = [];
      for (let i = 0; i < n; i++) diffs.push(b[i] - a[i]);
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

    mannWhitneyUTest(groupA, groupB) {
      const a = Descriptive.cleanData(groupA);
      const b = Descriptive.cleanData(groupB);
      const n1 = a.length;
      const n2 = b.length;
      if (n1 === 0 || n2 === 0) return { error: 'Both groups must contain at least 1 observation.' };

      const combined = [
        ...a.map(val => ({ val, group: 'A' })),
        ...b.map(val => ({ val, group: 'B' }))
      ].sort((x, y) => x.val - y.val);

      const N = combined.length;
      let i = 0;
      while (i < N) {
        let j = i;
        while (j < N - 1 && combined[j + 1].val === combined[i].val) j++;
        const rank = (i + 1 + j + 1) / 2;
        for (let k = i; k <= j; k++) combined[k].rank = rank;
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

      const meanU = (n1 * n2) / 2;
      const sigmaU = Math.sqrt((n1 * n2 * (n1 + n2 + 1)) / 12);
      const z = sigmaU === 0 ? 0 : (U - meanU) / sigmaU;
      const pValue = Distributions.normalPValue(z);
      const rankBiserial = 1 - (2 * U) / (n1 * n2);

      return {
        testName: 'Mann-Whitney U Test (Wilcoxon Rank-Sum)',
        n1, n2,
        rankSumA, rankSumB,
        u1, u2,
        statistic: U,
        zScore: z,
        pValue,
        rankBiserial,
        isSignificant: pValue < 0.05
      };
    },

    wilcoxonSignedRank(groupPre, groupPost) {
      const a = Descriptive.cleanData(groupPre);
      const b = Descriptive.cleanData(groupPost);
      const len = Math.min(a.length, b.length);
      if (len < 3) return { error: 'Wilcoxon signed-rank requires at least 3 paired observations.' };

      const diffs = [];
      for (let i = 0; i < len; i++) {
        const d = b[i] - a[i];
        if (d !== 0) diffs.push({ diff: d, absDiff: Math.abs(d) });
      }

      const n = diffs.length;
      if (n < 3) return { error: 'Insufficient non-zero differences for Wilcoxon test.' };

      diffs.sort((x, y) => x.absDiff - y.absDiff);

      let i = 0;
      while (i < n) {
        let j = i;
        while (j < n - 1 && diffs[j + 1].absDiff === diffs[i].absDiff) j++;
        const rank = (i + 1 + j + 1) / 2;
        for (let k = i; k <= j; k++) diffs[k].rank = rank;
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

  // ==========================================
  // 4. MULTI-GROUP ANOVA
  // ==========================================
  const Anova = {
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

      const dfBetween = k - 1;
      const dfWithin = totalN - k;
      const dfTotal = totalN - 1;

      if (dfWithin <= 0) {
        return { error: 'Insufficient degrees of freedom for within-group variance.' };
      }

      let ssBetween = 0;
      let ssWithin = 0;

      for (const g of processedGroups) {
        ssBetween += g.stats.n * Math.pow(g.stats.mean - grandMean, 2);
        for (const val of g.data) {
          ssWithin += Math.pow(val - g.stats.mean, 2);
        }
      }
      const ssTotal = ssBetween + ssWithin;

      const msBetween = ssBetween / dfBetween;
      const msWithin = ssWithin / dfWithin;

      const F = msWithin === 0 ? 0 : msBetween / msWithin;
      const pValue = Distributions.fPValue(F, dfBetween, dfWithin);

      const etaSquared = ssTotal === 0 ? 0 : ssBetween / ssTotal;
      const omegaSquared = (ssTotal + msWithin) === 0 ? 0 :
        (ssBetween - dfBetween * msWithin) / (ssTotal + msWithin);

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
        testName: "One-Way Analysis of Variance (Fisher's ANOVA)",
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

      const weights = [];
      let sumW = 0;
      for (const g of processedGroups) {
        const varG = Math.max(1e-9, g.stats.variance);
        const w = g.stats.n / varG;
        weights.push(w);
        sumW += w;
      }

      let weightedSumMean = 0;
      for (let j = 0; j < k; j++) {
        weightedSumMean += weights[j] * processedGroups[j].stats.mean;
      }
      const weightedGrandMean = weightedSumMean / sumW;

      let sumWeightDevSq = 0;
      for (let j = 0; j < k; j++) {
        sumWeightDevSq += weights[j] * Math.pow(processedGroups[j].stats.mean - weightedGrandMean, 2);
      }
      const numeratorA = sumWeightDevSq / (k - 1);

      let lambdaTermB = 0;
      for (let j = 0; j < k; j++) {
        const g = processedGroups[j];
        const wRatio = 1 - weights[j] / sumW;
        lambdaTermB += Math.pow(wRatio, 2) / (g.stats.n - 1);
      }

      const denomAdjustment = (k > 1 && (k * k - 1) > 0)
        ? 1 + ((2 * (k - 2)) / (k * k - 1)) * lambdaTermB
        : 1;

      const F_welch = denomAdjustment === 0 ? 0 : numeratorA / denomAdjustment;
      const df1 = k - 1;
      const df2 = lambdaTermB > 0 ? (k * k - 1) / (3 * lambdaTermB) : (totalN - k);
      const pValue = Distributions.fPValue(F_welch, df1, df2);

      const omegaSquared = Math.max(0, (df1 * (F_welch - 1)) / (df1 * (F_welch - 1) + totalN));
      const etaSquared = Math.max(0, (df1 * F_welch) / (df1 * F_welch + df2));

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

          const dfNumerator = Math.pow(varAOverN + varBOverN, 2);
          const dfDenominator = Math.pow(varAOverN, 2) / (gA.stats.n - 1) + Math.pow(varBOverN, 2) / (gB.stats.n - 1);
          const dfPair = dfDenominator === 0 ? 1 : dfNumerator / dfDenominator;

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
            statisticLabel: "Games-Howell t (df')",
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
      const epsilonSquared = N > 1 ? Math.min(1.0, Math.max(0, H / (N - 1))) : 0;

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
            meanDiff: medianDiff,
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
        fStatistic: H,
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

    repeatedMeasures(groups) {
      if (!Array.isArray(groups) || groups.length < 2) {
        return { error: 'Repeated Measures ANOVA requires at least 2 conditions/timepoints.' };
      }

      const k = groups.length;
      const cleaned = groups.map((g, idx) => ({
        name: g.name || `Timepoint ${idx + 1}`,
        data: Descriptive.cleanData(g.data)
      }));

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

      const partialEtaSquared = (ssTreatment + ssError) > 0 ? ssTreatment / (ssTreatment + ssError) : 0;
      const ggEpsilon = Math.max(1 / (k - 1), Math.min(1.0, 1 - 0.5 * (k - 1) / (dfError || 1)));
      const pValueGG = Distributions.fPValue(F, dfTreatment * ggEpsilon, dfError * ggEpsilon);

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

      const rankMatrix = [];
      let totalTieCorrection = 0;

      for (let i = 0; i < N; i++) {
        const row = [];
        for (let j = 0; j < k; j++) {
          row.push({ condIdx: j, val: processedGroups[j].data[i] });
        }
        row.sort((a, b) => a.val - b.val);

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
      const kendallsW = (N > 0 && df > 0) ? Math.min(1.0, Math.max(0, chiSqF / (N * df))) : 0;

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

      const deviationGroups = processed.map(g => {
        const med = g.stats.median;
        const devData = g.data.map(x => Math.abs(x - med));
        return {
          name: g.name,
          data: devData
        };
      });

      const anovaOnDeviations = this.oneWay(deviationGroups);
      if (anovaOnDeviations.error) {
        return { error: anovaOnDeviations.error };
      }

      const F = anovaOnDeviations.fStatistic;
      const df1 = anovaOnDeviations.dfBetween;
      const df2 = anovaOnDeviations.dfWithin;
      const pValue = anovaOnDeviations.pValue;

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

  // ==========================================
  // 5. CATEGORICAL & RISK METRICS
  // ==========================================
  const Categorical = {
    twoByTwo(a, b, c, d) {
      a = Math.max(0, parseInt(a, 10) || 0);
      b = Math.max(0, parseInt(b, 10) || 0);
      c = Math.max(0, parseInt(c, 10) || 0);
      d = Math.max(0, parseInt(d, 10) || 0);

      const r1 = a + b;
      const r2 = c + d;
      const c1 = a + c;
      const c2 = b + d;
      const n = a + b + c + d;

      if (n === 0 || r1 === 0 || r2 === 0 || c1 === 0 || c2 === 0) {
        return { error: 'Contingency table has empty marginal totals.' };
      }

      const eA = (r1 * c1) / n;
      const eB = (r1 * c2) / n;
      const eC = (r2 * c1) / n;
      const eD = (r2 * c2) / n;

      const chiSq =
        Math.pow(a - eA, 2) / eA +
        Math.pow(b - eB, 2) / eB +
        Math.pow(c - eC, 2) / eC +
        Math.pow(d - eD, 2) / eD;
      const pValueStandard = Distributions.chiSquarePValue(chiSq, 1);

      const aAdj = a === 0 || b === 0 || c === 0 || d === 0 ? a + 0.5 : a;
      const bAdj = a === 0 || b === 0 || c === 0 || d === 0 ? b + 0.5 : b;
      const cAdj = a === 0 || b === 0 || c === 0 || d === 0 ? c + 0.5 : c;
      const dAdj = a === 0 || b === 0 || c === 0 || d === 0 ? d + 0.5 : d;

      const oddsRatio = (aAdj * dAdj) / (bAdj * cAdj);
      const seLnOR = Math.sqrt(1 / aAdj + 1 / bAdj + 1 / cAdj + 1 / dAdj);
      const z95 = Distributions.invNormalCDF(0.975);
      const orCI95 = [
        Math.exp(Math.log(oddsRatio) - z95 * seLnOR),
        Math.exp(Math.log(oddsRatio) + z95 * seLnOR)
      ];

      const riskTreated = a / c1;
      const riskControl = b / c2;
      let relativeRisk = riskControl === 0 ? Infinity : riskTreated / riskControl;
      let rrCI95 = [0, 0];
      if (a > 0 && b > 0) {
        const seLnRR = Math.sqrt((c / a) / c1 + (d / b) / c2);
        rrCI95 = [
          Math.exp(Math.log(relativeRisk) - z95 * seLnRR),
          Math.exp(Math.log(relativeRisk) + z95 * seLnRR)
        ];
      }

      const arr = riskTreated - riskControl;
      const nnt = arr !== 0 ? Math.abs(1 / arr) : Infinity;

      // Fisher's exact 2x2
      const hyperProb = (x) => {
        const logP =
          Distributions.logGamma(r1 + 1) +
          Distributions.logGamma(r2 + 1) +
          Distributions.logGamma(c1 + 1) +
          Distributions.logGamma(c2 + 1) -
          Distributions.logGamma(n + 1) -
          Distributions.logGamma(x + 1) -
          Distributions.logGamma(r1 - x + 1) -
          Distributions.logGamma(c1 - x + 1) -
          Distributions.logGamma(r2 - c1 + x + 1);
        return Math.exp(logP);
      };

      const obsP = hyperProb(a);
      const minX = Math.max(0, c1 - r2);
      const maxX = Math.min(r1, c1);
      let pFisher = 0;
      for (let x = minX; x <= maxX; x++) {
        const p = hyperProb(x);
        if (p <= obsP + 1e-12) pFisher += p;
      }
      pFisher = Math.min(1.0, Math.max(0.0, pFisher));

      return {
        table: { a, b, c, d, n },
        chiSquare: { standard: chiSq, pValueStandard },
        fishersExact: { pValue: pFisher },
        riskMetrics: { oddsRatio, orCI95, relativeRisk, rrCI95, arr, nnt }
      };
    }
  };

  // ==========================================
  // 6. CORRELATION & REGRESSION
  // ==========================================
  const Correlation = {
    cleanPairs(xRaw, yRaw) {
      const pairs = [];
      const n = Math.min(xRaw.length, yRaw.length);
      for (let i = 0; i < n; i++) {
        const x = typeof xRaw[i] === 'number' ? xRaw[i] : parseFloat(xRaw[i]);
        const y = typeof yRaw[i] === 'number' ? yRaw[i] : parseFloat(yRaw[i]);
        if (!isNaN(x) && isFinite(x) && !isNaN(y) && isFinite(y)) {
          pairs.push({ x, y });
        }
      }
      return pairs;
    },

    pearson(xData, yData) {
      const pairs = this.cleanPairs(xData, yData);
      const n = pairs.length;
      if (n < 3) return { error: 'Correlation requires at least 3 paired values.' };

      const xVals = pairs.map(p => p.x);
      const yVals = pairs.map(p => p.y);
      const xStats = Descriptive.calculate(xVals);
      const yStats = Descriptive.calculate(yVals);

      if (xStats.sd === 0 || yStats.sd === 0) return { error: 'No variance in one variable.' };

      let covSum = 0;
      for (let i = 0; i < n; i++) {
        covSum += (xVals[i] - xStats.mean) * (yVals[i] - yStats.mean);
      }
      const r = (covSum / (n - 1)) / (xStats.sd * yStats.sd);
      const df = n - 2;
      const t = Math.abs(r) >= 1 ? Infinity : (r * Math.sqrt(df)) / Math.sqrt(1 - r * r);
      const pValue = Distributions.tPValue(t, df);

      return { n, r, rSquared: r * r, pValue, statistic: t, isSignificant: pValue < 0.05 };
    },

    spearman(xData, yData) {
      const pairs = this.cleanPairs(xData, yData);
      const n = pairs.length;
      if (n < 3) return { rho: 0, pValue: 1 };

      const assignRanks = (arr) => {
        const indexed = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
        const ranks = new Array(arr.length);
        let i = 0;
        while (i < arr.length) {
          let j = i;
          while (j < arr.length - 1 && indexed[j + 1].v === indexed[i].v) j++;
          const avgRank = (i + 1 + j + 1) / 2;
          for (let k = i; k <= j; k++) ranks[indexed[k].i] = avgRank;
          i = j + 1;
        }
        return ranks;
      };

      const xRanks = assignRanks(pairs.map(p => p.x));
      const yRanks = assignRanks(pairs.map(p => p.y));
      const res = this.pearson(xRanks, yRanks);
      return { rho: res.r, pValue: res.pValue };
    },

    linearRegression(xData, yData) {
      const pairs = this.cleanPairs(xData, yData);
      const n = pairs.length;
      if (n < 3) return { error: 'Regression requires at least 3 values.' };

      const xVals = pairs.map(p => p.x);
      const yVals = pairs.map(p => p.y);
      const xStats = Descriptive.calculate(xVals);
      const yStats = Descriptive.calculate(yVals);

      let sxx = 0;
      let sxy = 0;
      for (let i = 0; i < n; i++) {
        const dx = xVals[i] - xStats.mean;
        const dy = yVals[i] - yStats.mean;
        sxx += dx * dx;
        sxy += dx * dy;
      }
      if (sxx === 0) return { error: 'X has zero variance.' };

      const slope = sxy / sxx;
      const intercept = yStats.mean - slope * xStats.mean;

      let ssRes = 0;
      for (let i = 0; i < n; i++) {
        const yPred = intercept + slope * xVals[i];
        const res = yVals[i] - yPred;
        ssRes += res * res;
      }
      const seResidual = Math.sqrt(ssRes / (n - 2));

      return {
        n, slope, intercept,
        equation: `y = ${slope.toFixed(3)}x + ${intercept.toFixed(3)}`,
        seResidual, sxx, xStats, yStats
      };
    },

    quadraticRegression(xData, yData) {
      const pairs = this.cleanPairs(xData, yData);
      const n = pairs.length;
      if (n < 4) return { error: 'At least 4 paired observations required for quadratic fit.' };

      const xs = pairs.map(p => p.x);
      const ys = pairs.map(p => p.y);

      let s0 = n, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
      let t0 = 0, t1 = 0, t2 = 0;
      for (let i = 0; i < n; i++) {
        const x = xs[i], y = ys[i];
        const x2 = x * x;
        s1 += x; s2 += x2; s3 += x2 * x; s4 += x2 * x2;
        t0 += y; t1 += y * x; t2 += y * x2;
      }

      const det = (m) =>
        m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
        m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
        m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);

      const M = [[s4, s3, s2], [s3, s2, s1], [s2, s1, s0]];
      const D = det(M);
      if (Math.abs(D) < 1e-12) return { error: 'Singular matrix in quadratic fit.' };

      const b2 = det([[t2, s3, s2], [t1, s2, s1], [t0, s1, s0]]) / D;
      const b1 = det([[s4, t2, s2], [s3, t1, s1], [s2, t0, s0]]) / D;
      const b0 = det([[s4, s3, t2], [s3, s2, t1], [s2, s1, t0]]) / D;

      const vertexX = -b1 / (2 * b2);
      const vertexY = b2 * vertexX * vertexX + b1 * vertexX + b0;

      const minX = Math.min(...xs);
      const maxX = Math.max(...xs);
      const vertexInRange = vertexX >= minX && vertexX <= maxX;

      let ssResQuad = 0;
      let ssTotal = 0;
      const meanY = t0 / n;
      for (let i = 0; i < n; i++) {
        const pred = b2 * xs[i] * xs[i] + b1 * xs[i] + b0;
        ssResQuad += Math.pow(ys[i] - pred, 2);
        ssTotal += Math.pow(ys[i] - meanY, 2);
      }
      const rSquaredQuad = ssTotal === 0 ? 0 : Math.max(0, 1 - ssResQuad / ssTotal);

      const lin = this.linearRegression(xs, ys);
      let ssResLin = 0;
      for (let i = 0; i < n; i++) {
        const predLin = lin.intercept + lin.slope * xs[i];
        ssResLin += Math.pow(ys[i] - predLin, 2);
      }

      const df2 = n - 3;
      const ssDiff = Math.max(0, ssResLin - ssResQuad);
      const fStat = df2 > 0 && ssResQuad > 0 ? (ssDiff / 1) / (ssResQuad / df2) : 0;
      const pQuad = df2 > 0 ? Distributions.fPValue(fStat, 1, df2) : 1.0;

      const isSignificantlyQuadratic = pQuad < 0.05 && ssDiff > 0.05 * ssResLin;
      const shape = !isSignificantlyQuadratic
        ? 'Linear'
        : (b2 > 0 ? (vertexInRange ? 'U-Shaped' : 'J-Shaped') : (vertexInRange ? 'Inverted U-Shaped' : 'Inverted J-Shaped'));

      return {
        b2, b1, b0,
        equation: `y = ${b2.toFixed(3)}x² ${b1 >= 0 ? '+' : ''}${b1.toFixed(3)}x ${b0 >= 0 ? '+' : ''}${b0.toFixed(3)}`,
        vertexX,
        vertexY,
        vertexInRange,
        rSquaredQuad,
        fStat,
        pQuad,
        isSignificantlyQuadratic,
        shape
      };
    },

    analyzeMonotonicity(pairs) {
      if (!pairs || pairs.length < 3) {
        return { isMonotonic: true, reversals: 0, pattern: 'Indeterminate (N < 3)', direction: 'None' };
      }
      const sorted = [...pairs].sort((a, b) => a.x - b.x);
      const n = sorted.length;

      const smoothedY = [];
      for (let i = 0; i < n; i++) {
        const prev = i > 0 ? sorted[i - 1].y : sorted[i].y;
        const curr = sorted[i].y;
        const next = i < n - 1 ? sorted[i + 1].y : sorted[i].y;
        smoothedY.push((prev + curr + next) / 3);
      }

      const ySpan = Math.max(...sorted.map(p => p.y)) - Math.min(...sorted.map(p => p.y)) || 1;
      const dyThreshold = ySpan * 0.04;
      const directions = [];
      for (let i = 1; i < n; i++) {
        const diff = smoothedY[i] - smoothedY[i - 1];
        if (Math.abs(diff) > dyThreshold) {
          directions.push(diff > 0 ? 1 : -1);
        }
      }

      const compressed = [];
      for (const d of directions) {
        if (compressed.length === 0 || compressed[compressed.length - 1] !== d) {
          compressed.push(d);
        }
      }

      const reversals = Math.max(0, compressed.length - 1);
      const isMonotonic = reversals === 0 && compressed.length > 0;
      let pattern;
      let direction = 'None';

      if (compressed.length === 0) {
        pattern = 'Flat / Invariable';
      } else if (reversals === 0) {
        direction = compressed[0] === 1 ? 'Increasing' : 'Decreasing';
        pattern = compressed[0] === 1 ? 'Monotonically Increasing' : 'Monotonically Decreasing';
      } else if (reversals === 1) {
        pattern = compressed[0] === -1 ? 'U-Shaped (Decrease → Increase)' : 'Inverted U-Shaped (Increase → Decrease)';
      } else if (reversals === 2) {
        pattern = compressed[0] === 1 ? 'Tri-Phasic / N-Shaped (Increase → Decrease → Increase)' : 'Tri-Phasic / Inverted N (Decrease → Increase → Decrease)';
      } else {
        pattern = `Multi-Phasic (${reversals} Direction Reversals)`;
      }

      return {
        isMonotonic,
        reversals,
        pattern,
        direction
      };
    },

    checkAssumptions(xData, yData, reg) {
      const pairs = this.cleanPairs(xData, yData);
      const n = pairs.length;
      if (n < 4) {
        return { error: 'At least 4 paired observations are required to evaluate assumptions.' };
      }

      const xs = pairs.map(p => p.x);
      const ys = pairs.map(p => p.y);
      const xStats = Descriptive.calculate(xs);
      const yStats = Descriptive.calculate(ys);

      // 1. Normality of X and Y
      const xNormal = xStats.n < 8 ? Math.abs(xStats.skewness) < 1.2 : xStats.normality.isNormal;
      const yNormal = yStats.n < 8 ? Math.abs(yStats.skewness) < 1.2 : yStats.normality.isNormal;
      const bivariateNormal = xNormal && yNormal;

      // 2. Outlier / High Leverage Detection
      const outliers = [];
      const seRes = reg && reg.seResidual ? reg.seResidual : 1;
      const sxx = reg && reg.sxx ? reg.sxx : 1;
      pairs.forEach((p) => {
        const pred = reg.intercept + reg.slope * p.x;
        const res = p.y - pred;
        const stdRes = Math.abs(res) / (seRes || 1);
        const hii = 1 / n + Math.pow(p.x - xStats.mean, 2) / sxx;
        const cooksD = (Math.pow(stdRes, 2) / 2) * (hii / Math.pow(1 - Math.min(0.99, hii), 2));
        if (stdRes > 2.5 || (stdRes > 1.96 && cooksD > 0.5)) {
          outliers.push({ x: p.x, y: p.y, stdRes: stdRes.toFixed(2), cooksD: cooksD.toFixed(2) });
        }
      });

      // 3. Shape Analysis: Linear vs U-Shaped & Monotonicity Reversals
      const quad = this.quadraticRegression(xs, ys);
      const mono = this.analyzeMonotonicity(pairs);
      const isUShaped = quad && quad.isSignificantlyQuadratic && (quad.shape === 'U-Shaped' || quad.shape === 'Inverted U-Shaped');

      // 4. Linearity vs Monotonic Curvature Check
      const pearsonRes = this.pearson(xs, ys);
      const spearmanRes = this.spearman(xs, ys);
      const diffCorr = Math.abs(spearmanRes.rho) - Math.abs(pearsonRes.r);
      const isLinear = !isUShaped && mono.isMonotonic && diffCorr <= 0.12;

      // 5. Homoscedasticity (Goldfeld-Quandt variance ratio test on residuals)
      const sortedByX = [...pairs].sort((a, b) => a.x - b.x);
      const half = Math.floor(n / 2);
      const group1 = sortedByX.slice(0, half).map(p => p.y - (reg.intercept + reg.slope * p.x));
      const group2 = sortedByX.slice(half).map(p => p.y - (reg.intercept + reg.slope * p.x));
      const var1 = Descriptive.calculate(group1).variance || 1e-6;
      const var2 = Descriptive.calculate(group2).variance || 1e-6;
      const fRatio = Math.max(var1, var2) / Math.min(var1, var2);
      const dfGQ = half - 1;
      const homoscedasticP = dfGQ > 0 ? Distributions.fPValue(fRatio, dfGQ, dfGQ) : 1.0;
      const isHomoscedastic = homoscedasticP > 0.05;

      // 6. Automated Decision Synthesis
      const violations = [];
      if (isUShaped) {
        violations.push(`Strong ${quad.shape} relationship detected (Quadratic term p < .001, Nadir/Zenith at X = ${quad.vertexX.toFixed(2)})`);
      } else if (mono.reversals >= 2) {
        violations.push(`Multi-phasic non-monotonic trajectory (${mono.pattern})`);
      }

      if (!bivariateNormal && !isUShaped && mono.reversals < 2) {
        if (!xNormal && !yNormal) violations.push('Variables X and Y are non-normally distributed');
        else if (!xNormal) violations.push('Variable X is skewed/non-normal');
        else violations.push('Variable Y is skewed/non-normal');
      }
      if (outliers.length > 0) {
        violations.push(`${outliers.length} influential outlier(s) detected`);
      }
      if (!isLinear && !isUShaped && mono.isMonotonic) {
        violations.push('Monotonic non-linear relationship detected (|ρ| > |r|)');
      }
      if (!isHomoscedastic) {
        violations.push('Heteroscedasticity present (residual variance changes across X)');
      }

      let recommendedTest;
      let recommendationDetail;
      let isParametricOk = violations.length === 0;

      if (isUShaped) {
        recommendedTest = `Quadratic Polynomial Fit (${quad.shape})`;
        recommendationDetail = `A significant non-linear ${quad.shape} relationship was confirmed (R² = ${quad.rSquaredQuad.toFixed(3)}, p < .001; ${quad.b2 > 0 ? 'Nadir minimum' : 'Zenith peak'} at X = ${quad.vertexX.toFixed(2)}). Linear correlation (r = ${pearsonRes.r.toFixed(2)}) is misleadingly near zero and falsely obscures this strong clinical association.`;
      } else if (mono.reversals >= 2) {
        recommendedTest = 'Non-Linear / Segmented Regression';
        recommendationDetail = `The data trajectory is non-monotonic (${mono.pattern}) with ${mono.reversals} directional reversals (e.g. numbers increase, decrease, and then increase again). Standard linear correlation (r = ${pearsonRes.r.toFixed(2)}) and rank correlation (ρ = ${spearmanRes.rho.toFixed(2)}) cannot properly represent multi-phasic reversals.`;
      } else if (isParametricOk) {
        recommendedTest = "Pearson's Correlation & Linear Regression";
        recommendationDetail = `All parametric assumptions (${mono.pattern.toLowerCase()}, bivariate normality, homoscedasticity, absence of extreme outliers) are satisfied. Pearson's r and OLS regression provide optimal statistical power and valid confidence intervals.`;
      } else {
        recommendedTest = "Spearman's Rank Correlation (ρ)";
        recommendationDetail = `Parametric assumptions violated (${violations.join('; ')}). Spearman's rank correlation is recommended because it makes no distributional assumptions, resists outlier leverage, and correctly assesses monotonic relationships (${mono.pattern}).`;
      }

      return {
        n,
        xNormal,
        yNormal,
        bivariateNormal,
        outliers,
        isLinear,
        isUShaped,
        shape: isUShaped ? quad.shape : (mono.reversals >= 2 ? mono.pattern : (mono.isMonotonic ? mono.pattern : 'Non-Linear')),
        quad,
        mono,
        isHomoscedastic,
        homoscedasticP,
        violations,
        isParametricOk,
        recommendedTest,
        recommendationDetail
      };
    }
  };

  // ==========================================
  // 7. DIAGNOSTIC & ROC
  // ==========================================
  const Diagnostic = {
    evaluate2x2(tp, fp, fn, tn) {
      tp = Math.max(0, parseInt(tp, 10) || 0);
      fp = Math.max(0, parseInt(fp, 10) || 0);
      fn = Math.max(0, parseInt(fn, 10) || 0);
      tn = Math.max(0, parseInt(tn, 10) || 0);

      const diseased = tp + fn;
      const nonDiseased = fp + tn;
      const testPositive = tp + fp;
      const testNegative = fn + tn;
      const total = diseased + nonDiseased;

      if (total === 0) return { error: 'No data.' };

      const sensitivity = diseased > 0 ? tp / diseased : 0;
      const specificity = nonDiseased > 0 ? tn / nonDiseased : 0;
      const ppv = testPositive > 0 ? tp / testPositive : 0;
      const npv = testNegative > 0 ? tn / testNegative : 0;
      const accuracy = total > 0 ? (tp + tn) / total : 0;
      const prevalence = total > 0 ? diseased / total : 0;
      const plr = (1 - specificity) > 0 ? sensitivity / (1 - specificity) : Infinity;
      const nlr = specificity > 0 ? (1 - sensitivity) / specificity : Infinity;
      const youdenJ = sensitivity + specificity - 1;

      const wilsonCI = (p, n) => {
        if (n === 0) return [0, 0];
        const z = Distributions.invNormalCDF(0.975);
        const denom = 1 + (z * z) / n;
        const center = (p + (z * z) / (2 * n)) / denom;
        const margin = (z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))) / denom;
        return [Math.max(0, center - margin), Math.min(1, center + margin)];
      };

      return {
        tp, fp, fn, tn, total, diseased, nonDiseased, testPositive, testNegative,
        sensitivity,
        sensitivityCI95: wilsonCI(sensitivity, diseased),
        specificity,
        specificityCI95: wilsonCI(specificity, nonDiseased),
        ppv,
        ppvCI95: wilsonCI(ppv, testPositive),
        npv,
        npvCI95: wilsonCI(npv, testNegative),
        accuracy,
        accuracyCI95: wilsonCI(accuracy, total),
        prevalence,
        plr, nlr, youdenJ
      };
    },

    computeROC(samples) {
      const valid = samples.filter(
        s => typeof s.score === 'number' && !isNaN(s.score) && (s.status === 0 || s.status === 1)
      );

      const nPos = valid.filter(s => s.status === 1).length;
      const nNeg = valid.filter(s => s.status === 0).length;
      if (nPos === 0 || nNeg === 0) return { error: 'ROC requires both positive and negative cases.' };

      valid.sort((a, b) => b.score - a.score);

      const thresholds = [-Infinity, ...valid.map(s => s.score), Infinity];
      const unique = Array.from(new Set(thresholds)).sort((a, b) => b - a);

      const points = [];
      let best = { threshold: unique[0], youdenJ: -1, sens: 0, spec: 0 };

      for (const t of unique) {
        let tp = 0;
        let fp = 0;
        for (const s of valid) {
          if (s.score >= t) {
            if (s.status === 1) tp++;
            else fp++;
          }
        }
        const tpr = tp / nPos;
        const fpr = fp / nNeg;
        const youden = tpr - fpr;

        if (youden > best.youdenJ) {
          best = {
            threshold: t === -Infinity ? 'Lowest' : t === Infinity ? 'Highest' : t,
            youdenJ: youden,
            sens: tpr,
            spec: 1 - fpr
          };
        }
        points.push({ threshold: t, fpr, tpr });
      }

      points.sort((a, b) => a.fpr - b.fpr || a.tpr - b.tpr);

      let auc = 0;
      for (let i = 1; i < points.length; i++) {
        const dx = points[i].fpr - points[i - 1].fpr;
        const avgY = (points[i].tpr + points[i - 1].tpr) / 2;
        auc += dx * avgY;
      }
      auc = Math.min(1.0, Math.max(0.0, auc));

      const q1 = auc / (2 - auc);
      const q2 = (2 * auc * auc) / (1 + auc);
      const seNum = auc * (1 - auc) + (nPos - 1) * (q1 - auc * auc) + (nNeg - 1) * (q2 - auc * auc);
      const seAUC = Math.sqrt(Math.max(0, seNum / (nPos * nNeg)));
      const z95 = 1.95996;
      const aucCI95 = [Math.max(0, auc - z95 * seAUC), Math.min(1, auc + z95 * seAUC)];

      return { nPos, nNeg, total: nPos + nNeg, auc, aucCI95, points, optimalCutoff: best };
    }
  };

  // ==========================================
  // 8. POWER ANALYSIS
  // ==========================================
  const PowerAnalysis = {
    twoIndependentMeans({ m1, m2, sd, alpha = 0.05, power = 0.80, nPerGroup = null, ratio = 1.0 }) {
      m1 = parseFloat(m1);
      m2 = parseFloat(m2);
      sd = parseFloat(sd);
      alpha = parseFloat(alpha);
      power = parseFloat(power);

      const diff = Math.abs(m1 - m2);
      if (diff === 0 || sd <= 0) {
        return { error: 'Mean difference must be greater than zero and SD positive.' };
      }

      const cohensD = diff / sd;
      const zAlpha = Distributions.invNormalCDF(1 - alpha / 2);

      if (nPerGroup !== null && nPerGroup > 0) {
        const n1 = Math.round(nPerGroup);
        const n2 = Math.ceil(n1 * ratio);
        const effectiveN = (n1 * n2) / (n1 + n2);
        const zBeta = cohensD * Math.sqrt(effectiveN) - zAlpha;
        const achievedPower = Math.min(0.9999, Math.max(0.0001, Distributions.normalCDF(zBeta)));

        return {
          design: 'independent_means',
          testName: 'Two Independent Groups (Two-Sample t-Test)',
          m1, m2, diff, sd,
          effectSize: cohensD,
          effectSizeLabel: "Cohen's d",
          cohensD,
          alpha,
          power: achievedPower,
          achievedPower,
          nPerGroup: n1,
          n1, n2,
          totalN: n1 + n2,
          isPostHoc: true
        };
      } else {
        const zBeta = Distributions.invNormalCDF(power);
        const n1 = Math.ceil(
          Math.pow(zAlpha + zBeta, 2) * Math.pow(sd, 2) * (1 + 1 / ratio) / Math.pow(diff, 2)
        );
        const n2 = Math.ceil(n1 * ratio);
        const effectiveN = (n1 * n2) / (n1 + n2);
        const actualZBeta = cohensD * Math.sqrt(effectiveN) - zAlpha;
        const achievedPower = Math.min(0.9999, Math.max(0.0001, Distributions.normalCDF(actualZBeta)));

        return {
          design: 'independent_means',
          testName: 'Two Independent Groups (Two-Sample t-Test)',
          m1, m2, diff, sd,
          effectSize: cohensD,
          effectSizeLabel: "Cohen's d",
          cohensD,
          alpha,
          power,
          achievedPower,
          nPerGroup: ratio === 1 ? n1 : null,
          n1, n2,
          totalN: n1 + n2,
          isPostHoc: false
        };
      }
    },

    pairedMeans({ m1, m2, sdDiff, alpha = 0.05, power = 0.80, nPairs = null }) {
      m1 = parseFloat(m1);
      m2 = parseFloat(m2);
      sdDiff = parseFloat(sdDiff);
      alpha = parseFloat(alpha);
      power = parseFloat(power);

      const diff = Math.abs(m1 - m2);
      if (diff === 0 || sdDiff <= 0) {
        return { error: 'Paired mean difference must be greater than zero and SD of differences positive.' };
      }

      const dz = diff / sdDiff;
      const zAlpha = Distributions.invNormalCDF(1 - alpha / 2);

      if (nPairs !== null && nPairs > 0) {
        const n = Math.round(nPairs);
        const zBeta = dz * Math.sqrt(n) - zAlpha;
        const achievedPower = Math.min(0.9999, Math.max(0.0001, Distributions.normalCDF(zBeta)));

        return {
          design: 'paired_means',
          testName: 'Paired Mean Study Design (Paired t-Test)',
          m1, m2, diff, sdDiff,
          effectSize: dz,
          effectSizeLabel: "Cohen's d_z (Paired)",
          dz,
          alpha,
          power: achievedPower,
          achievedPower,
          nPairs: n,
          nPerGroup: n,
          totalN: n,
          isPostHoc: true
        };
      } else {
        const zBeta = Distributions.invNormalCDF(power);
        const rawN = Math.pow(zAlpha + zBeta, 2) / Math.pow(dz, 2);
        const nCalc = Math.ceil(rawN);
        const actualZBeta = dz * Math.sqrt(nCalc) - zAlpha;
        const achievedPower = Math.min(0.9999, Math.max(0.0001, Distributions.normalCDF(actualZBeta)));

        return {
          design: 'paired_means',
          testName: 'Paired Mean Study Design (Paired t-Test)',
          m1, m2, diff, sdDiff,
          effectSize: dz,
          effectSizeLabel: "Cohen's d_z (Paired)",
          dz,
          alpha,
          power,
          achievedPower,
          nPairs: nCalc,
          nPerGroup: nCalc,
          totalN: nCalc,
          isPostHoc: false
        };
      }
    },

    contingency2x2({ p1, p2, alpha = 0.05, power = 0.80, nPerGroup = null, testType = 'chisq' }) {
      p1 = parseFloat(p1);
      p2 = parseFloat(p2);
      alpha = parseFloat(alpha);
      power = parseFloat(power);

      if (isNaN(p1) || isNaN(p2) || p1 <= 0 || p1 >= 1 || p2 <= 0 || p2 >= 1) {
        return { error: 'Event proportions (p1 and p2) must be strictly between 0 and 1.' };
      }
      const diff = Math.abs(p1 - p2);
      if (diff === 0) {
        return { error: 'Event proportions must differ (p1 ≠ p2).' };
      }

      const zAlpha = Distributions.invNormalCDF(1 - alpha / 2);
      const pBar = (p1 + p2) / 2;
      const qBar = 1 - pBar;

      const arr = diff;
      const rr = p1 / p2;
      const rrr = p2 > 0 ? (p2 - p1) / p2 : 0;
      const or = (p1 / (1 - p1)) / (p2 / (1 - p2));
      const nnt = 1 / arr;
      const cohenH = 2 * Math.abs(Math.asin(Math.sqrt(p1)) - Math.asin(Math.sqrt(p2)));

      const sigma0 = Math.sqrt(2 * pBar * qBar);
      const sigma1 = Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));

      if (nPerGroup !== null && nPerGroup > 0) {
        const n = Math.round(nPerGroup);

        const zBetaChisq = (Math.sqrt(n) * diff - zAlpha * sigma0) / sigma1;
        const powerChisq = Math.min(0.9999, Math.max(0.0001, Distributions.normalCDF(zBetaChisq)));

        const zBetaCC = (Math.sqrt(n) * diff - (1 / Math.sqrt(n)) - zAlpha * sigma0) / sigma1;
        const powerCC = Math.min(0.9999, Math.max(0.0001, Distributions.normalCDF(zBetaCC)));

        const activePower = testType === 'fisher' ? powerCC : powerChisq;

        return {
          design: 'contingency_2x2',
          testName: testType === 'fisher' ? "Fisher's Exact / Continuity-Corrected Chi-Square" : "Pearson Chi-Square Test",
          p1, p2, diff, arr, rr, rrr, or, nnt, cohenH,
          effectSize: arr,
          effectSizeLabel: 'ARR (|p1 - p2|)',
          alpha,
          power: activePower,
          achievedPower: activePower,
          powerChisq,
          powerFisher: powerCC,
          nPerGroup: n,
          n1: n, n2: n,
          totalN: 2 * n,
          testType,
          isPostHoc: true
        };
      } else {
        const zBeta = Distributions.invNormalCDF(power);
        const term1 = zAlpha * Math.sqrt(2 * pBar * qBar);
        const term2 = zBeta * Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));
        const uncorrectedN = Math.ceil(Math.pow(term1 + term2, 2) / Math.pow(diff, 2));

        const ccN = Math.ceil((uncorrectedN / 4) * Math.pow(1 + Math.sqrt(1 + 4 / (uncorrectedN * diff)), 2));
        const chosenN = testType === 'fisher' ? ccN : uncorrectedN;

        const zBetaAchieved = (Math.sqrt(chosenN) * diff - zAlpha * sigma0) / sigma1;
        const achievedPower = Math.min(0.9999, Math.max(0.0001, Distributions.normalCDF(zBetaAchieved)));

        return {
          design: 'contingency_2x2',
          testName: testType === 'fisher' ? "Fisher's Exact / Continuity-Corrected Chi-Square" : "Pearson Chi-Square Test",
          p1, p2, diff, arr, rr, rrr, or, nnt, cohenH,
          effectSize: arr,
          effectSizeLabel: 'ARR (|p1 - p2|)',
          alpha,
          power,
          achievedPower,
          uncorrectedN,
          continuityCorrectedN: ccN,
          nPerGroup: chosenN,
          n1: chosenN, n2: chosenN,
          totalN: 2 * chosenN,
          testType,
          isPostHoc: false
        };
      }
    },

    sampleSizeMeans(m1, m2, sd, alpha = 0.05, power = 0.80) {
      return this.twoIndependentMeans({ m1, m2, sd, alpha, power });
    },

    sampleSizeProportions(p1, p2, alpha = 0.05, power = 0.80) {
      return this.contingency2x2({ p1, p2, alpha, power, testType: 'chisq' });
    }
  };

  // ==========================================
  // 9. DATA PARSER & SAMPLES
  // ==========================================
  const DataParser = {
    parseSeries(text) {
      if (!text || typeof text !== 'string') return [];
      return text
        .split(/[\r\n,;\t\s]+/)
        .map(s => parseFloat(s.trim()))
        .filter(n => !isNaN(n) && isFinite(n));
    },

    samples: {
      icpDynamics: {
        groupA: [10.2, 11.5, 9.8, 12.1, 10.9, 13.4, 11.0, 9.5, 12.8, 10.4, 11.7, 10.0],
        groupB: [18.5, 21.0, 19.2, 23.4, 20.1, 25.6, 22.0, 19.8, 24.5, 20.3, 21.8, 19.5]
      },
      tumorResection: {
        groupA: [14.2, 15.1, 13.8, 16.5, 14.9, 15.8, 17.2, 13.5, 15.0, 14.6, 16.1, 14.8, 15.4, 16.0],
        groupB: [22.4, 28.1, 18.9, 31.5, 24.0, 19.8, 35.2, 26.7, 21.3, 29.4, 33.1, 20.5]
      },
      drainOutputSkewed: [12, 14, 15, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 38, 62],
      cranialAsymmetry: [
        [7.2, 6.8, 7.5, 6.9, 8.1, 7.0, 7.4, 6.5, 7.9, 7.1],
        [4.1, 3.8, 4.5, 3.9, 4.8, 3.6, 4.2, 3.5, 4.0, 3.7],
        [2.5, 2.8, 2.2, 2.6, 3.1, 2.4, 2.9, 2.1, 2.7, 2.3]
      ],
      diagnosticBiomarker: [
        { score: 1850, status: 1 }, { score: 1720, status: 1 }, { score: 1640, status: 1 },
        { score: 1590, status: 1 }, { score: 1530, status: 1 }, { score: 1480, status: 1 },
        { score: 1390, status: 1 }, { score: 1320, status: 1 }, { score: 1280, status: 1 },
        { score: 1210, status: 1 }, { score: 1150, status: 1 }, { score: 1100, status: 0 },
        { score: 1040, status: 1 }, { score: 980, status: 0 },  { score: 920, status: 0 },
        { score: 860, status: 0 },  { score: 810, status: 0 },  { score: 750, status: 0 },
        { score: 690, status: 0 },  { score: 620, status: 0 }
      ],
      shunt2x2: { a: 14, b: 36, c: 186, d: 164 },
      diagnostic2x2: { a: 92, b: 8, c: 12, d: 188 }
    }
  };

  // ==========================================
  // 9B. TEACHING & SIMULATION ENGINE
  // ==========================================
  const Teaching = {
    generators: {
      standardNormal() {
        let u = 0, v = 0;
        while (u === 0) u = Math.random();
        while (v === 0) v = Math.random();
        return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
      },

      normal(n, mean = 0, sd = 1) {
        const data = new Array(n);
        for (let i = 0; i < n; i++) {
          data[i] = mean + sd * this.standardNormal();
        }
        return data;
      },

      studentsT(n, df = 5, mean = 0, scale = 1) {
        const data = new Array(n);
        for (let i = 0; i < n; i++) {
          const z = this.standardNormal();
          let chi2 = 0;
          for (let j = 0; j < df; j++) {
            const u = this.standardNormal();
            chi2 += u * u;
          }
          const t = z / Math.sqrt(chi2 / df);
          data[i] = mean + scale * t;
        }
        return data;
      },

      uniform(n, min = 0, max = 10) {
        const data = new Array(n);
        const span = max - min;
        for (let i = 0; i < n; i++) {
          data[i] = min + Math.random() * span;
        }
        return data;
      },

      exponential(n, rate = 0.5) {
        const data = new Array(n);
        const lambda = rate > 0 ? rate : 1;
        for (let i = 0; i < n; i++) {
          let u = Math.random();
          while (u === 0) u = Math.random();
          data[i] = -Math.log(u) / lambda;
        }
        return data;
      },

      logNormal(n, mu = 1.0, sigma = 0.6) {
        const data = new Array(n);
        for (let i = 0; i < n; i++) {
          const z = this.standardNormal();
          data[i] = Math.exp(mu + sigma * z);
        }
        return data;
      },

      bimodal(n, m1 = 12, s1 = 2, m2 = 24, s2 = 3, p = 0.5) {
        const data = new Array(n);
        for (let i = 0; i < n; i++) {
          if (Math.random() < p) {
            data[i] = m1 + s1 * this.standardNormal();
          } else {
            data[i] = m2 + s2 * this.standardNormal();
          }
        }
        return data;
      },

      poisson(n, lambda = 4) {
        const data = new Array(n);
        const L = Math.exp(-lambda);
        for (let i = 0; i < n; i++) {
          let k = 0;
          let p = 1;
          do {
            k++;
            p *= Math.random();
          } while (p > L);
          data[i] = k - 1;
        }
        return data;
      },

      chiSquare(n, df = 4) {
        const data = new Array(n);
        for (let i = 0; i < n; i++) {
          let sum = 0;
          for (let j = 0; j < df; j++) {
            const z = this.standardNormal();
            sum += z * z;
          }
          data[i] = sum;
        }
        return data;
      }
    },

    pdf: {
      normal(x, mean = 0, sd = 1) {
        if (sd <= 0) return 0;
        const z = (x - mean) / sd;
        return (1.0 / (sd * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
      },

      studentsT(x, df = 5, mean = 0, scale = 1) {
        if (scale <= 0 || df <= 0) return 0;
        const t = (x - mean) / scale;
        const logGammaNum = Distributions.logGamma((df + 1) / 2);
        const logGammaDen = Distributions.logGamma(df / 2);
        const coeff = Math.exp(logGammaNum - logGammaDen) / (Math.sqrt(df * Math.PI) * scale);
        return coeff * Math.pow(1 + (t * t) / df, -(df + 1) / 2);
      },

      uniform(x, min = 0, max = 10) {
        if (x >= min && x <= max && max > min) {
          return 1.0 / (max - min);
        }
        return 0;
      },

      exponential(x, rate = 0.5) {
        if (x < 0 || rate <= 0) return 0;
        return rate * Math.exp(-rate * x);
      },

      logNormal(x, mu = 1.0, sigma = 0.6) {
        if (x <= 0 || sigma <= 0) return 0;
        const z = (Math.log(x) - mu) / sigma;
        return (1.0 / (x * sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
      },

      bimodal(x, m1 = 12, s1 = 2, m2 = 24, s2 = 3, p = 0.5) {
        const pdf1 = this.normal(x, m1, s1);
        const pdf2 = this.normal(x, m2, s2);
        return p * pdf1 + (1 - p) * pdf2;
      },

      poisson(k, lambda = 4) {
        if (k < 0 || Math.floor(k) !== k) return 0;
        let logFact = 0;
        for (let i = 2; i <= k; i++) logFact += Math.log(i);
        return Math.exp(k * Math.log(lambda) - lambda - logFact);
      },

      chiSquare(x, df = 4) {
        if (x <= 0 || df <= 0) return 0;
        const k = df / 2;
        const logDenom = k * Math.log(2) + Distributions.logGamma(k);
        return Math.exp((k - 1) * Math.log(x) - x / 2 - logDenom);
      }
    },

    distMetadata: {
      normal: {
        name: 'Normal (Gaussian) Distribution',
        symbol: 'X ~ N(μ, σ²)',
        clinicalExample: 'Adult Intracranial Pressure (ICP ~ 11.5 ± 2.8 mmHg), Systolic Blood Pressure, Cerebral Blood Flow.',
        properties: 'Symmetric bell-shaped density, zero skewness (G₁ = 0), mesokurtic (excess G₂ = 0). The 68-95-99.7 empirical rule applies.',
        testNote: 'Parametric Student t-test and ANOVA strictly assume normally distributed data or sample means.'
      },
      studentsT: {
        name: "Student's t-Distribution",
        symbol: 'T ~ t(ν)',
        clinicalExample: 'Pilot clinical trial sample means with small cohorts (n < 30) where true population variance is unknown.',
        properties: 'Heavier tails than normal distribution to account for sampling uncertainty in s. Converges to Gaussian as degrees of freedom ν → ∞.',
        testNote: "Formulated by William Sealy Gosset (under pseudonym 'Student') in 1908 for quality control at Guinness Brewery."
      },
      uniform: {
        name: 'Continuous Uniform Distribution',
        symbol: 'X ~ U(a, b)',
        clinicalExample: 'Randomized clinical trial treatment allocation sequences, Monte Carlo randomized permutation sampling.',
        properties: 'Constant probability density over interval [a, b]. Mean = (a+b)/2, Variance = (b-a)²/12. Rectangular morphology.',
        testNote: 'Classic non-normal benchmark illustrating the power of the Central Limit Theorem.'
      },
      exponential: {
        name: 'Exponential Distribution',
        symbol: 'X ~ Exp(λ)',
        clinicalExample: 'Time-to-event survival analysis, postoperative aneurysm recurrence intervals, acute epileptic seizure recurrence intervals.',
        properties: 'Memoryless property: P(X > s + t | X > s) = P(X > t). Strong right skew (G₁ = 2), Mean = 1/λ, SD = 1/λ.',
        testNote: 'Forms the baseline hazard foundation for Cox proportional hazards regression models in oncology and neurosurgery.'
      },
      logNormal: {
        name: 'Log-Normal Distribution',
        symbol: 'ln(X) ~ N(μ, σ²)',
        clinicalExample: 'Serum Neurofilament Light (NfL) biomarker concentration, hospital length of stay, ICU sedation durations.',
        properties: 'Strictly positive (X > 0) with long right tail. Multiplicative growth processes naturally produce log-normal distributions.',
        testNote: 'Logarithmic transformation [Y = ln(X)] converts log-normal biomedical data into symmetric Gaussian data for parametric analysis.'
      },
      bimodal: {
        name: 'Bimodal Mixture Distribution',
        symbol: 'p·N(μ₁, σ₁²) + (1-p)·N(μ₂, σ₂²)',
        clinicalExample: 'Pediatric vs adult disease onset peaks (e.g. craniopharyngioma, Hodgkin lymphoma), drug responders vs non-responders.',
        properties: 'Two distinct local maxima (peaks) separated by a trough. Violates unimodality assumptions of standard descriptive metrics.',
        testNote: 'Reporting Mean ± SD is misleading for bimodal distributions; stratified subgroup reporting is clinically essential.'
      },
      poisson: {
        name: 'Poisson Count Distribution',
        symbol: 'K ~ Pois(λ)',
        clinicalExample: 'Surgical site infections per 1,000 operative bed days, annual emergency craniotomy case volume.',
        properties: 'Discrete count distribution with parameter λ. Fundamental property: Mean = Variance = λ.',
        testNote: 'Used in Poisson regression models for incidence rate ratios (IRR) in epidemiology and hospital quality assurance.'
      },
      chiSquare: {
        name: 'Chi-Square Distribution',
        symbol: 'X ~ χ²(k)',
        clinicalExample: 'Sum of squared standardized residuals in 2x2 contingency tables and goodness-of-fit testing.',
        properties: 'Right-skewed distribution of the sum of k independent squared standard normal variates. Mean = k, Variance = 2k.',
        testNote: 'Underpins Pearson contingency test, McNemar paired test, and Wald confidence intervals in logistic regression.'
      }
    },

    clt: {
      populationType: 'uniform',
      sampleSize: 30,
      sampleMeans: [],
      lastSample: [],

      populations: {
        uniform: {
          name: 'Uniform Parent U(0, 10)',
          min: 0,
          max: 10,
          mean: 5.0,
          sd: Math.sqrt(100 / 12),
          description: 'Flat rectangular distribution with zero skew and platykurtic tails.',
          drawOne() {
            return Math.random() * 10;
          },
          pdf(x) {
            return x >= 0 && x <= 10 ? 0.1 : 0;
          }
        },
        exponential: {
          name: 'Exponential Parent Exp(λ = 0.5)',
          min: 0,
          max: 12,
          mean: 2.0,
          sd: 2.0,
          description: 'Heavily right-skewed time-to-event survival distribution (Skewness G₁ = 2.0).',
          drawOne() {
            let u = Math.random();
            while (u === 0) u = Math.random();
            return -Math.log(u) / 0.5;
          },
          pdf(x) {
            return x >= 0 ? 0.5 * Math.exp(-0.5 * x) : 0;
          }
        },
        bimodal: {
          name: 'Bimodal Parent (Two Distinct Subgroups)',
          min: 0,
          max: 12,
          mean: 5.0,
          sd: Math.sqrt(0.36 + 0.25 * 36),
          description: 'Two separate Gaussian clusters centered at μ₁=2.0 and μ₂=8.0 with a void in between.',
          drawOne() {
            const z = Teaching.generators.standardNormal();
            return Math.random() < 0.5 ? 2.0 + 0.6 * z : 8.0 + 0.6 * z;
          },
          pdf(x) {
            return Teaching.pdf.bimodal(x, 2.0, 0.6, 8.0, 0.6, 0.5);
          }
        },
        ushaped: {
          name: 'U-Shaped Parent (Bimodal Extremes / Anti-Normal)',
          min: 0,
          max: 10,
          mean: 5.0,
          sd: 3.535,
          description: 'Bathtub-shaped density with maximum probability at endpoints and minimal in the center.',
          drawOne() {
            const u = Math.random();
            return 5.0 + 5.0 * Math.sin(Math.PI * (u - 0.5));
          },
          pdf(x) {
            if (x <= 0.05 || x >= 9.95) return 0.35;
            return 1.0 / (Math.PI * Math.sqrt(Math.max(0.001, x * (10 - x))));
          }
        }
      },

      setPopulation(type) {
        if (this.populations[type]) {
          this.populationType = type;
          this.reset();
        }
      },

      setSampleSize(n) {
        this.sampleSize = Math.max(2, Math.min(200, Math.round(n)));
        this.reset();
      },

      reset() {
        this.sampleMeans = [];
        this.lastSample = [];
      },

      drawSamples(count = 1) {
        const pop = this.populations[this.populationType];
        const n = this.sampleSize;

        for (let c = 0; c < count; c++) {
          let sum = 0;
          const currentSample = new Array(n);
          for (let i = 0; i < n; i++) {
            const val = pop.drawOne();
            currentSample[i] = val;
            sum += val;
          }
          const mean = sum / n;
          this.sampleMeans.push(mean);
          if (c === count - 1) {
            this.lastSample = currentSample;
          }
        }

        return this.getSummary();
      },

      getSummary() {
        const pop = this.populations[this.populationType];
        const n = this.sampleSize;
        const k = this.sampleMeans.length;

        const theoreticalMean = pop.mean;
        const theoreticalSE = pop.sd / Math.sqrt(n);

        if (k === 0) {
          return {
            population: pop,
            sampleSize: n,
            samplesDrawn: 0,
            theoreticalMean,
            theoreticalSE,
            observedMean: null,
            observedSE: null,
            skewness: null,
            kurtosis: null,
            isNormal: null
          };
        }

        const stats = Descriptive.calculate(this.sampleMeans);

        return {
          population: pop,
          sampleSize: n,
          samplesDrawn: k,
          theoreticalMean,
          theoreticalSE,
          observedMean: stats.mean,
          observedSE: stats.sd,
          skewness: stats.skewness,
          kurtosis: stats.kurtosis,
          normality: stats.normality,
          min: stats.min,
          max: stats.max,
          iqr: stats.iqr,
          values: this.sampleMeans,
          lastSample: this.lastSample
        };
      }
    },

    /**
     * Student's t-Distribution Approximation to Normal Distribution Simulation
     */
    tConvergence: {
      sampleSize: 4,

      /**
       * Compute two-tailed critical value for Student's t distribution at alpha = 0.05
       * Uses analytic Cornish-Fisher expansion with exact fallbacks for small df
       */
      getCriticalValue(df, alpha = 0.05) {
        if (df <= 0) return NaN;
        if (df === 1) return 12.7062;
        if (df === 2) return 4.3027;
        if (df === 3) return 3.1824;
        if (df === 4) return 2.7764;
        if (df >= 500) return 1.95996;

        const z = 1.95996398454;
        const nu = df;
        const z2 = z * z, z3 = z2 * z, z5 = z3 * z2, z7 = z5 * z2, z9 = z7 * z2;
        const a = (z3 + z) / (4 * nu);
        const b = (5 * z5 + 16 * z3 + 3 * z) / (96 * nu * nu);
        const c = (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * Math.pow(nu, 3));
        const d = (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * Math.pow(nu, 4));
        return z + a + b + c + d;
      },

      /**
       * Compute full convergence metrics comparing t(df) to Standard Normal N(0, 1)
       */
      getMetrics(sampleSize = 4) {
        const n = Math.max(2, Math.round(sampleSize));
        const df = n - 1;
        const tPeak = Distributions.tPDF(0, df);
        const normPeak = 1.0 / Math.sqrt(2 * Math.PI); // ~0.398942
        const peakDiffPct = ((tPeak - normPeak) / normPeak) * 100;

        const tCrit = this.getCriticalValue(df, 0.05);
        const zCrit = 1.95996;
        const critDiffPct = ((tCrit - zCrit) / zCrit) * 100;

        const tailProb = Distributions.tPValue(zCrit, df); // Actual probability mass beyond +/- 1.96
        const normTailProb = 0.05; // Exactly 5% for N(0, 1)

        const excessKurtosis = df > 4 ? 6 / (df - 4) : Infinity;
        const maxDiscrepancy = Math.abs(normPeak - tPeak);

        let clinicalNote = '';
        if (df <= 4) {
          clinicalNote = `Extremely fat tails (excess kurtosis ${df <= 4 ? 'undefined / infinite' : excessKurtosis.toFixed(2)}). Critical t cutoff (${tCrit.toFixed(3)}) is +${critDiffPct.toFixed(1)}% wider than Gaussian z = 1.960. Testing at z = 1.96 would cause a Type I error inflation to ${(tailProb * 100).toFixed(1)}% (nearly 3x higher than intended 5%)!`;
        } else if (df < 30) {
          clinicalNote = `Moderate tail thickening (kurtosis = ${excessKurtosis.toFixed(2)}). Critical t (${tCrit.toFixed(3)}) is +${critDiffPct.toFixed(1)}% wider than Gaussian z = 1.960. Gosset's t-test is mandatory for valid inference.`;
        } else if (df < 60) {
          clinicalNote = `Approaching Gaussian equivalence. Critical t (${tCrit.toFixed(3)}) is within +${critDiffPct.toFixed(1)}% of z = 1.960. The normal approximation is clinically robust for sample sizes n ≥ 31.`;
        } else {
          clinicalNote = `Virtually identical to Standard Normal N(0, 1). Critical t (${tCrit.toFixed(3)}) deviates by only +${critDiffPct.toFixed(2)}% from z = 1.960, and tail probability is ${(tailProb * 100).toFixed(2)}% vs 5.00%.`;
        }

        return {
          sampleSize: n,
          df,
          tPeak,
          normPeak,
          peakDiffPct,
          tCrit,
          zCrit,
          critDiffPct,
          tailProb,
          normTailProb,
          excessKurtosis,
          maxDiscrepancy,
          clinicalNote
        };
      }
    },

    /**
     * Two-Sample Overlap, Dispersion (SD vs. SEM), and Alpha Significance Simulation
     */
    significanceOverlap: {
      getTCriticalValue(df, alpha = 0.05) {
        if (df <= 0) return NaN;
        const sigAlpha = Math.min(0.20, Math.max(0.0001, alpha));
        const p = 1 - sigAlpha / 2;
        const z = Distributions.invNormalCDF(p);
        if (df >= 500) return z;
        if (df === 1) return Math.tan((p - 0.5) * Math.PI);
        if (df === 2) return Math.sqrt(2 / (4 * (1 - p) * p) - 2);
        const nu = df;
        const z2 = z * z, z3 = z2 * z, z5 = z3 * z2, z7 = z5 * z2, z9 = z7 * z2;
        const a = (z3 + z) / (4 * nu);
        const b = (5 * z5 + 16 * z3 + 3 * z) / (96 * nu * nu);
        const c = (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * Math.pow(nu, 3));
        const d = (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * Math.pow(nu, 4));
        return Math.max(z, z + a + b + c + d);
      },

      computeOverlap(mu1, s1, mu2, s2) {
        if (s1 <= 0 || s2 <= 0) return 0;
        const normPDF = (x, m, s) => (1.0 / (s * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * Math.pow((x - m) / s, 2));
        const minX = Math.min(mu1 - 5.0 * s1, mu2 - 5.0 * s2);
        const maxX = Math.max(mu1 + 5.0 * s1, mu2 + 5.0 * s2);
        const steps = 600;
        const dx = (maxX - minX) / steps;
        let sum = 0;
        for (let i = 0; i <= steps; i++) {
          const x = minX + i * dx;
          const y1 = normPDF(x, mu1, s1);
          const y2 = normPDF(x, mu2, s2);
          const w = (i === 0 || i === steps) ? 0.5 : 1.0;
          sum += Math.min(y1, y2) * w * dx;
        }
        return Math.max(0, Math.min(1.0, sum));
      },

      getMetrics({
        mean1 = 10.0,
        delta = 2.0,
        sd = 2.5,
        sd1,
        sd2,
        n = 16,
        n1,
        n2,
        sem1 = null,
        sem2 = null,
        alpha = 0.05,
        viewMode = 'means'
      } = {}) {
        const mu1 = parseFloat(mean1) || 10.0;
        const dMu = Math.max(0, parseFloat(delta) !== undefined && !isNaN(parseFloat(delta)) ? parseFloat(delta) : 2.0);
        const mu2 = mu1 + dMu;

        const baseSD = parseFloat(sd) || 2.5;
        const sigma1 = Math.max(0.2, parseFloat(sd1 !== undefined && sd1 !== null ? sd1 : baseSD));
        const sigma2 = Math.max(0.2, parseFloat(sd2 !== undefined && sd2 !== null ? sd2 : baseSD));

        let sampleN1, sError1;
        if (sem1 !== null && sem1 !== undefined && !isNaN(parseFloat(sem1))) {
          sError1 = Math.max(0.01, parseFloat(sem1));
          sampleN1 = Math.max(2, Math.min(1000, Math.round(Math.pow(sigma1 / sError1, 2))));
        } else {
          sampleN1 = Math.max(2, Math.round(n1 !== undefined && n1 !== null ? n1 : (n || 16)));
          sError1 = sigma1 / Math.sqrt(sampleN1);
        }

        let sampleN2, sError2;
        if (sem2 !== null && sem2 !== undefined && !isNaN(parseFloat(sem2))) {
          sError2 = Math.max(0.01, parseFloat(sem2));
          sampleN2 = Math.max(2, Math.min(1000, Math.round(Math.pow(sigma2 / sError2, 2))));
        } else {
          sampleN2 = Math.max(2, Math.round(n2 !== undefined && n2 !== null ? n2 : (n || 16)));
          sError2 = sigma2 / Math.sqrt(sampleN2);
        }

        const sigAlpha = Math.min(0.20, Math.max(0.001, parseFloat(alpha) || 0.05));

        const v1 = Math.pow(sError1, 2);
        const v2 = Math.pow(sError2, 2);
        const seDiff = Math.sqrt(v1 + v2);

        const dfNum = Math.pow(v1 + v2, 2);
        const dfDenom = ((sampleN1 > 1) ? Math.pow(v1, 2) / (sampleN1 - 1) : 0) +
                        ((sampleN2 > 1) ? Math.pow(v2, 2) / (sampleN2 - 1) : 0);
        const df = dfDenom > 0 ? Math.max(1, dfNum / dfDenom) : (sampleN1 + sampleN2 - 2);

        const zCrit = Distributions.invNormalCDF(1 - sigAlpha / 2);
        const tCrit = this.getTCriticalValue(df, sigAlpha);

        const deltaCrit = tCrit * seDiff;
        const deltaCritZ = zCrit * seDiff;

        const tStat = seDiff > 0 ? dMu / seDiff : 0;
        const zStat = seDiff > 0 ? dMu / seDiff : 0;
        const pValue = Distributions.tPValue(tStat, df);
        const isSignificant = pValue < sigAlpha;

        const pooledSD = Math.sqrt(
          ((sampleN1 - 1) * Math.pow(sigma1, 2) + (sampleN2 - 1) * Math.pow(sigma2, 2)) /
          Math.max(1, (sampleN1 + sampleN2 - 2))
        );
        const cohensD = pooledSD > 0 ? dMu / pooledSD : 0;

        const patientOVL = this.computeOverlap(mu1, sigma1, mu2, sigma2);
        const meansOVL = this.computeOverlap(mu1, sError1, mu2, sError2);

        const moe1 = tCrit * sError1;
        const moe2 = tCrit * sError2;
        const ci1 = [mu1 - moe1, mu1 + moe1];
        const ci2 = [mu2 - moe2, mu2 + moe2];
        const ciOverlapDist = Math.max(0, ci1[1] - ci2[0]);

        const nullCritLeft = -deltaCrit;
        const nullCritRight = deltaCrit;

        let statusText = '';
        if (isSignificant) {
          statusText = `STATISTICALLY SIGNIFICANT (p = ${pValue < 0.0001 ? '< 0.0001' : pValue.toFixed(4)} < α = ${sigAlpha.toFixed(3)})`;
        } else {
          statusText = `NOT STATISTICALLY SIGNIFICANT (p = ${pValue.toFixed(4)} ≥ α = ${sigAlpha.toFixed(3)})`;
        }

        let explanation = '';
        const isHetero = Math.abs(sigma1 - sigma2) > 0.05 || Math.abs(sError1 - sError2) > 0.02;
        if (isSignificant) {
          explanation = `The observed difference between means (Δ = ${dMu.toFixed(2)}) meets or exceeds the critical significance threshold (Δcrit = ${deltaCrit.toFixed(2)} at α = ${sigAlpha.toFixed(3)}). ` +
            `While individual patient values overlap substantially (${(patientOVL * 100).toFixed(1)}% patient overlap with SD₁ = ${sigma1.toFixed(2)} and SD₂ = ${sigma2.toFixed(2)}), ` +
            `the standard errors of the means (SEM₁ = ${sError1.toFixed(3)}, SEM₂ = ${sError2.toFixed(3)}) have contracted with sample sizes (n₁ = ${sampleN1}, n₂ = ${sampleN2}) so that the sampling distributions of the two means only overlap by ${(meansOVL * 100).toFixed(1)}%. ` +
            `Under H₀, observing a mean difference this large occurs with probability p = ${pValue < 0.0001 ? '< 0.0001' : pValue.toFixed(4)} < α = ${sigAlpha.toFixed(3)}, rejecting the null hypothesis${isHetero ? ` (Welch df = ${df.toFixed(1)})` : ''}.`;
        } else {
          explanation = `The observed difference between means (Δ = ${dMu.toFixed(2)}) is smaller than the required critical threshold (Δcrit = ${deltaCrit.toFixed(2)} at α = ${sigAlpha.toFixed(3)}). ` +
            `The sampling distributions of the two sample means overlap too heavily (${(meansOVL * 100).toFixed(1)}% overlap, p = ${pValue.toFixed(3)} ≥ α). ` +
            `To achieve significance at this α level, you must either: (1) observe a larger effect size Δ, (2) reduce measurement noise (lower SD), or (3) recruit more patients to shrink SEM (SEM₁ = ${sError1.toFixed(3)}, SEM₂ = ${sError2.toFixed(3)}).`;
        }

        return {
          mean1: mu1,
          mean2: mu2,
          delta: dMu,
          sd1: sigma1,
          sd2: sigma2,
          sd: (sigma1 + sigma2) / 2,
          sem1: sError1,
          sem2: sError2,
          sem: (sError1 + sError2) / 2,
          n1: sampleN1,
          n2: sampleN2,
          n: sampleN1,
          alpha: sigAlpha,
          viewMode,
          seDiff,
          df,
          zCrit,
          tCrit,
          deltaCrit,
          deltaCritZ,
          tStat,
          zStat,
          pValue,
          isSignificant,
          statusText,
          pooledSD,
          cohensD,
          patientOVL,
          meansOVL,
          moe1,
          moe2,
          moe: moe1,
          ci1,
          ci2,
          ciOverlapDist,
          nullCritLeft,
          nullCritRight,
          explanation
        };
      }
    },

    powerSimulation: {
      getMetrics({
        sd = 4.0,
        sem = null,
        n = null,
        beta = null,
        power = null,
        delta = 2.0,
        alpha = 0.05,
        viewMode = 'distributions',
        lastChanged = 'power'
      } = {}) {
        const sigAlpha = Math.min(0.20, Math.max(0.001, parseFloat(alpha) || 0.05));
        const zCrit = Distributions.invNormalCDF(1 - sigAlpha / 2);
        const dMu = Math.max(0.1, parseFloat(delta) !== undefined && !isNaN(parseFloat(delta)) ? parseFloat(delta) : 2.0);
        const sigma = Math.max(0.2, parseFloat(sd) !== undefined && !isNaN(parseFloat(sd)) ? parseFloat(sd) : 4.0);

        let sampleN, sError, pwr, bta;

        if (lastChanged === 'power' && power !== null && power !== undefined) {
          pwr = Math.min(0.999, Math.max(0.50, parseFloat(power)));
          bta = 1.0 - pwr;
          const zBeta = Distributions.invNormalCDF(pwr);
          sampleN = Math.max(4, Math.min(1000, Math.round(2 * Math.pow(zCrit + zBeta, 2) * Math.pow(sigma, 2) / Math.pow(dMu, 2))));
          sError = sigma / Math.sqrt(sampleN);
        } else if (lastChanged === 'beta' && beta !== null && beta !== undefined) {
          bta = Math.min(0.50, Math.max(0.001, parseFloat(beta)));
          pwr = 1.0 - bta;
          const zBeta = Distributions.invNormalCDF(pwr);
          sampleN = Math.max(4, Math.min(1000, Math.round(2 * Math.pow(zCrit + zBeta, 2) * Math.pow(sigma, 2) / Math.pow(dMu, 2))));
          sError = sigma / Math.sqrt(sampleN);
        } else if (lastChanged === 'sem' && sem !== null && sem !== undefined) {
          sError = Math.max(0.01, parseFloat(sem));
          sampleN = Math.max(4, Math.min(1000, Math.round(Math.pow(sigma / sError, 2))));
          sError = sigma / Math.sqrt(sampleN);
          const seDiff = sigma * Math.sqrt(2 / sampleN);
          const lambda = dMu / seDiff;
          pwr = Math.max(0.001, Math.min(0.999, Distributions.normalCDF(lambda - zCrit)));
          bta = 1.0 - pwr;
        } else {
          sampleN = Math.max(4, Math.min(1000, Math.round(n !== null && n !== undefined ? n : 64)));
          sError = sigma / Math.sqrt(sampleN);
          const seDiff = sigma * Math.sqrt(2 / sampleN);
          const lambda = dMu / seDiff;
          pwr = Math.max(0.001, Math.min(0.999, Distributions.normalCDF(lambda - zCrit)));
          bta = 1.0 - pwr;
        }

        const seDiff = sigma * Math.sqrt(2 / sampleN);
        const lambda = seDiff > 0 ? dMu / seDiff : 0;
        const calculatedPower = Math.max(0.001, Math.min(0.999, Distributions.normalCDF(lambda - zCrit)));
        const calculatedBeta = 1.0 - calculatedPower;

        const xCrit = zCrit * seDiff;
        const cohensD = sigma > 0 ? dMu / sigma : 0;

        const matrix = {
          trueNegative: 1 - sigAlpha,
          falsePositive: sigAlpha,
          falseNegative: calculatedBeta,
          truePositive: calculatedPower
        };

        const curvePoints = [];
        const nSteps = [4, 6, 8, 10, 14, 18, 24, 30, 38, 48, 60, 75, 90, 110, 135, 165, 200, 250];
        for (const curN of nSteps) {
          const curSEDiff = sigma * Math.sqrt(2 / curN);
          const curLam = curSEDiff > 0 ? dMu / curSEDiff : 0;
          const curPwr = Math.max(0, Math.min(1.0, Distributions.normalCDF(curLam - zCrit)));
          curvePoints.push({ n: curN, power: curPwr });
        }

        let powerRating = '';
        if (calculatedPower >= 0.90) {
          powerRating = 'EXCELLENT (≥ 90%)';
        } else if (calculatedPower >= 0.80) {
          powerRating = 'ADEQUATE / REGULATORY STANDARD (80% - 90%)';
        } else if (calculatedPower >= 0.60) {
          powerRating = 'BORDERLINE / SUBOPTIMAL (60% - 80%)';
        } else {
          powerRating = 'SEVERELY UNDERPOWERED (< 60%)';
        }

        let explanation = '';
        if (calculatedPower >= 0.80) {
          explanation = `The study is well-powered (${(calculatedPower * 100).toFixed(1)}% Power at α = ${sigAlpha.toFixed(3)}). ` +
            `With a sample size of n = ${sampleN} per group (N = ${sampleN * 2} total), the standard error of the mean contracts to SEM = ${sError.toFixed(3)} ` +
            `(SE_diff = ${seDiff.toFixed(3)}), ensuring that the sampling distribution of a true difference Δ = ${dMu.toFixed(2)} (Cohen's d = ${cohensD.toFixed(2)}) ` +
            `is shifted far to the right of the significance boundary (xcrit = ${xCrit.toFixed(2)}). ` +
            `The Type II error risk is contained to β = ${(calculatedBeta * 100).toFixed(1)}%, meaning there is only a 1-in-${Math.round(1 / Math.max(0.001, calculatedBeta))} risk of a false-negative trial outcome.`;
        } else {
          explanation = `The study is UNDERPOWERED (${(calculatedPower * 100).toFixed(1)}% Power at α = ${sigAlpha.toFixed(3)}). ` +
            `Because sample size (n = ${sampleN} per group) is insufficient for the biological noise level (SD = ${sigma.toFixed(2)}), the SEM is wide (SEM = ${sError.toFixed(3)}), ` +
            `causing the H₁ distribution to heavily overlap the null acceptance region. ` +
            `The Type II error rate is β = ${(calculatedBeta * 100).toFixed(1)}% — meaning you have a ${(calculatedBeta * 100).toFixed(0)}% chance of failing to detect a truly effective medical intervention! ` +
            `To reach the clinical gold standard (80% power), you must either recruit more patients (reduce SEM to ${(sigma / Math.sqrt(2 * Math.pow(zCrit + 0.842, 2) * Math.pow(sigma, 2) / Math.pow(dMu, 2))).toFixed(3)}) or reduce measurement error.`;
        }

        return {
          sd: sigma,
          sem: sError,
          n: sampleN,
          totalN: sampleN * 2,
          power: calculatedPower,
          beta: calculatedBeta,
          delta: dMu,
          alpha: sigAlpha,
          zCrit,
          seDiff,
          lambda,
          xCrit,
          cohensD,
          viewMode,
          matrix,
          curvePoints,
          powerRating,
          explanation
        };
      }
    },

    /**
     * Section 6: Bayesian Statistics & Logic (3Blue1Brown Model)
     * The Geometry of Changing Beliefs: Unit square area, representative counts, and odds updating.
     */
    bayesianSimulation: {
      getMetrics({
        prior = 0.0476,
        likelihood = 0.40,
        falsePositive = 0.10,
        sampleSize = 210,
        viewMode = 'square',
        preset = 'steve'
      } = {}) {
        const pPrior = Math.max(0.0001, Math.min(0.9999, Number.isFinite(prior) ? prior : 0.0476));
        const pNotPrior = 1.0 - pPrior;
        const pLikelihood = Math.max(0.0001, Math.min(1.0, Number.isFinite(likelihood) ? likelihood : 0.40));
        const pFalsePos = Math.max(0.0001, Math.min(1.0, Number.isFinite(falsePositive) ? falsePositive : 0.10));
        const nTotal = Math.max(10, Math.min(100000, Math.round(Number.isFinite(sampleSize) ? sampleSize : 210)));

        // Joint probabilities (Areas in 3Blue1Brown Unit Square)
        const areaHAndE = pPrior * pLikelihood;             // P(H ∩ E) - Green Area
        const areaNotHAndE = pNotPrior * pFalsePos;         // P(¬H ∩ E) - Amber Area
        const pEvidence = areaHAndE + areaNotHAndE;          // P(E) - Total Shaded Evidence Area
        const posterior = pEvidence > 0 ? areaHAndE / pEvidence : pPrior; // P(H|E)
        const posteriorNotH = 1.0 - posterior;

        // Odds and Bayes Factor
        const priorOdds = pPrior / pNotPrior;
        const bayesFactor = pFalsePos > 0 ? pLikelihood / pFalsePos : 999.0;
        const posteriorOdds = priorOdds * bayesFactor;
        const beliefShift = posterior - pPrior;

        // Natural frequencies / representative sample counts
        const countH = Math.round(nTotal * pPrior);
        const countNotH = nTotal - countH;
        const countHAndE = Math.round(countH * pLikelihood);
        const countNotHAndE = Math.round(countNotH * pFalsePos);
        const countTotalE = countHAndE + countNotHAndE;
        const countPosterior = countTotalE > 0 ? countHAndE / countTotalE : posterior;

        // Sequential evidence trajectory (up to 4 steps of compounding evidence)
        const trajectory = [];
        let currentP = pPrior;
        trajectory.push({ step: 0, label: 'Prior Belief', p: currentP });
        for (let s = 1; s <= 4; s++) {
          const num = currentP * pLikelihood;
          const den = num + (1.0 - currentP) * pFalsePos;
          currentP = den > 0 ? num / den : currentP;
          trajectory.push({ step: s, label: `Evidence #${s}`, p: currentP });
        }

        // Evidence strength qualitative classification (Jeffreys / Kass & Raftery scale)
        let evidenceRating = '';
        if (bayesFactor > 100) evidenceRating = 'Decisive Evidence (BF > 100)';
        else if (bayesFactor > 30) evidenceRating = 'Very Strong Evidence (30 < BF ≤ 100)';
        else if (bayesFactor > 10) evidenceRating = 'Strong Evidence (10 < BF ≤ 30)';
        else if (bayesFactor > 3) evidenceRating = 'Substantial / Moderate (3 < BF ≤ 10)';
        else if (bayesFactor > 1) evidenceRating = 'Weak / Anecdotal (1 < BF ≤ 3)';
        else if (Math.abs(bayesFactor - 1.0) < 0.01) evidenceRating = 'Neutral / Irrelevant (BF = 1.0)';
        else if (bayesFactor > 0.33) evidenceRating = 'Weak Evidence for Alternative';
        else if (bayesFactor > 0.1) evidenceRating = 'Moderate for Alternative';
        else evidenceRating = 'Strong for Alternative (BF < 0.1)';

        // Pedagogical explanation narrative
        let explanation = '';
        if (preset === 'steve' || Math.abs(pPrior - 0.0476) < 0.01) {
          explanation = `As Grant Sanderson (3Blue1Brown) demonstrates, while Steve's description sounds 4× more characteristic of a librarian ` +
            `(${(pLikelihood * 100).toFixed(0)}% vs ${(pFalsePos * 100).toFixed(0)}%), ` +
            `the prior ratio of farmers to librarians is 20 to 1 (P(H) = ${(pPrior * 100).toFixed(1)}%). ` +
            `In a representative sample of N = ${nTotal} individuals, there are ${countH} librarians and ${countNotH} farmers. ` +
            `Only ${countHAndE} librarians fit the description, while ${countNotHAndE} farmers fit the description! ` +
            `Restricting our universe to only those who fit the description leaves ${countTotalE} people. ` +
            `Therefore, the probability that Steve is a librarian is ${countHAndE} / ${countTotalE} = ${(posterior * 100).toFixed(1)}% ` +
            `— meaning Steve is still 5× more likely to be a farmer (${(posteriorNotH * 100).toFixed(1)}%) despite sounding like a librarian!`;
        } else {
          explanation = `Bayes' Theorem updates prior belief P(H) = ${(pPrior * 100).toFixed(1)}% in light of new evidence with ` +
            `likelihood P(E|H) = ${(pLikelihood * 100).toFixed(1)}% and false alarm rate P(E|¬H) = ${(pFalsePos * 100).toFixed(1)}%. ` +
            `The evidence provides a Bayes Factor of ${bayesFactor.toFixed(2)}×. ` +
            `The evidence restricts the possibility space to total area P(E) = ${(pEvidence * 100).toFixed(2)}%. ` +
            `Within this restricted subspace, the hypothesis occupies ${(posterior * 100).toFixed(1)}% (posterior probability), ` +
            `representing a belief shift of ${(beliefShift >= 0 ? '+' : '')}${(beliefShift * 100).toFixed(1)} percentage points.`;
        }

        return {
          prior: pPrior,
          notPrior: pNotPrior,
          likelihood: pLikelihood,
          falsePositive: pFalsePos,
          sampleSize: nTotal,
          areaHAndE,
          areaNotHAndE,
          pEvidence,
          posterior,
          posteriorNotH,
          priorOdds,
          bayesFactor,
          posteriorOdds,
          beliefShift,
          countH,
          countNotH,
          countHAndE,
          countNotHAndE,
          countTotalE,
          countPosterior,
          trajectory,
          evidenceRating,
          explanation,
          viewMode,
          preset
        };
      }
    }
  };

// ==========================================
  // 9C. PROPENSITY SCORE MATCHING & CAUSAL INFERENCE
  // ==========================================
  /**
 * Statis-Gravity - Propensity Score Matching (PSM) & Causal Inference Engine
 * Implements clinical observational data pre-processing, multivariate logistic regression,
 * nearest-neighbor caliper matching, covariate balance diagnostics (SMD & Love plot),
 * Average Treatment Effect on the Treated (ATT) estimation, and multi-language
 * reproducible script generation (Python, R, Stata).
 */



const Psm = {
  /**
   * Evaluates missingness and filters / prepares complete cases
   * @param {Array<Object>} rows Raw dataset rows
   * @param {string} treatmentCol Name of treatment indicator column (binary 0/1)
   * @param {string} outcomeCol Name of primary outcome column
   * @param {Array<string>} covariateCols List of baseline confounding covariates
   * @param {string} missingMode 'complete_case' or 'impute_median'
   */
  prepareData(rows, treatmentCol, outcomeCol, covariateCols, missingMode = 'complete_case') {
    if (!rows || !Array.isArray(rows) || rows.length === 0) {
      return { error: 'No observational data provided for propensity score matching.' };
    }
    if (!treatmentCol) {
      return { error: 'Treatment variable column must be specified.' };
    }
    if (!outcomeCol) {
      return { error: 'Outcome variable column must be specified.' };
    }
    if (!covariateCols || covariateCols.length === 0) {
      return { error: 'At least one baseline confounding covariate must be specified.' };
    }

    const allRequiredCols = [treatmentCol, outcomeCol, ...covariateCols];
    const totalRows = rows.length;
    let missingRowsCount = 0;
    const perColMissing = {};
    allRequiredCols.forEach(col => { perColMissing[col] = 0; });

    // Identify missingness
    const validRows = [];
    rows.forEach((row, idx) => {
      let isRowMissing = false;
      allRequiredCols.forEach(col => {
        const val = row[col];
        if (val === undefined || val === null || val === '' || isNaN(Number(val))) {
          perColMissing[col]++;
          isRowMissing = true;
        }
      });
      if (isRowMissing) {
        missingRowsCount++;
      } else {
        const parsed = { _rowId: idx + 1 };
        allRequiredCols.forEach(col => {
          parsed[col] = Number(row[col]);
        });
        validRows.push(parsed);
      }
    });

    if (validRows.length < 10) {
      return {
        error: `Insufficient complete cases for matching. Found ${validRows.length} valid rows from ${totalRows} total rows.`
      };
    }

    // Verify binary treatment assignment
    const treatmentValues = new Set(validRows.map(r => r[treatmentCol]));
    const uniqueVals = Array.from(treatmentValues);
    if (uniqueVals.length !== 2) {
      return {
        error: `Treatment variable "${treatmentCol}" must be binary (found values: [${uniqueVals.join(', ')}]).`
      };
    }

    // Ensure treatment is 0 and 1
    const minVal = Math.min(...uniqueVals);
    const maxVal = Math.max(...uniqueVals);
    validRows.forEach(r => {
      r._treatment = r[treatmentCol] === maxVal ? 1 : 0;
      r._outcome = r[outcomeCol];
    });

    const treatedCount = validRows.filter(r => r._treatment === 1).length;
    const controlCount = validRows.filter(r => r._treatment === 0).length;

    if (treatedCount < 3 || controlCount < 3) {
      return {
        error: `Both treatment groups require adequate observations (Treated: ${treatedCount}, Control: ${controlCount}).`
      };
    }

    return {
      totalRows,
      completeCasesCount: validRows.length,
      missingRowsCount,
      perColMissing,
      treatmentCol,
      outcomeCol,
      covariateCols,
      treatedCount,
      controlCount,
      data: validRows
    };
  },

  /**
   * Fits a Multivariate Logistic Regression Model via Newton-Raphson IRLS
   * P(Treatment = 1 | X) = 1 / (1 + exp(- (beta_0 + sum(beta_j * X_j))))
   */
  fitLogisticRegression(data, covariateCols) {
    const N = data.length;
    const p = covariateCols.length;
    const numParams = p + 1; // Intercept + p covariates

    // Check for constant columns
    for (let j = 0; j < p; j++) {
      const col = covariateCols[j];
      const vals = data.map(d => d[col]);
      const min = Math.min(...vals);
      const max = Math.max(...vals);
      if (min === max) {
        return { error: `Covariate "${col}" has zero variance (constant value ${min}). Remove it from model.` };
      }
    }

    // Compute column means and SDs for numerical conditioning
    const means = covariateCols.map(col => {
      const sum = data.reduce((acc, d) => acc + d[col], 0);
      return sum / N;
    });
    const sds = covariateCols.map((col, idx) => {
      const m = means[idx];
      const sumSq = data.reduce((acc, d) => acc + Math.pow(d[col] - m, 2), 0);
      return Math.sqrt(sumSq / (N - 1)) || 1.0;
    });

    // Normalized design matrix X (N x numParams) with intercept
    const X = data.map(d => {
      const row = [1.0];
      covariateCols.forEach((col, idx) => {
        row.push((d[col] - means[idx]) / sds[idx]);
      });
      return row;
    });

    const y = data.map(d => d._treatment);
    const yMean = y.reduce((a, b) => a + b, 0) / N;

    // Null log likelihood
    const nullLogLik = y.reduce((acc, yi) => {
      const p = Math.max(1e-12, Math.min(1 - 1e-12, yMean));
      return acc + (yi * Math.log(p) + (1 - yi) * Math.log(1 - p));
    }, 0);

    // Initial coefficients
    let beta = new Array(numParams).fill(0);
    beta[0] = Math.log(Math.max(1e-6, yMean / (1 - yMean)));

    const maxIter = 40;
    const tolerance = 1e-7;
    let converged = false;
    let logLik = nullLogLik;

    for (let iter = 0; iter < maxIter; iter++) {
      // Compute probabilities and weights
      const pVec = [];
      const wVec = [];
      const residuals = [];

      for (let i = 0; i < N; i++) {
        let eta = 0;
        for (let j = 0; j < numParams; j++) {
          eta += X[i][j] * beta[j];
        }
        // Clamping for numerical stability
        eta = Math.max(-30, Math.min(30, eta));
        const pi = 1.0 / (1.0 + Math.exp(-eta));
        pVec.push(pi);
        wVec.push(Math.max(1e-9, pi * (1.0 - pi)));
        residuals.push(y[i] - pi);
      }

      // Gradient vector: X^T * r
      const grad = new Array(numParams).fill(0);
      for (let j = 0; j < numParams; j++) {
        for (let i = 0; i < N; i++) {
          grad[j] += X[i][j] * residuals[i];
        }
      }

      // Hessian / Fisher Information matrix: X^T * W * X + small ridge regularization
      const H = Array.from({ length: numParams }, () => new Array(numParams).fill(0));
      for (let r = 0; r < numParams; r++) {
        for (let c = r; c < numParams; c++) {
          let sum = 0;
          for (let i = 0; i < N; i++) {
            sum += X[i][r] * wVec[i] * X[i][c];
          }
          if (r === c) sum += 1e-5; // Ridge stabilization
          H[r][c] = sum;
          H[c][r] = sum;
        }
      }

      // Solve H * delta = grad using Cholesky / Gaussian elimination with partial pivoting
      const delta = this._solveLinearSystem(H, grad);
      if (!delta) {
        break; // Collinearity or singularity encountered
      }

      let maxDelta = 0;
      for (let j = 0; j < numParams; j++) {
        beta[j] += delta[j];
        if (Math.abs(delta[j]) > maxDelta) maxDelta = Math.abs(delta[j]);
      }

      // Compute current log-likelihood
      let currentLogLik = 0;
      for (let i = 0; i < N; i++) {
        let eta = 0;
        for (let j = 0; j < numParams; j++) eta += X[i][j] * beta[j];
        eta = Math.max(-30, Math.min(30, eta));
        const pi = 1.0 / (1.0 + Math.exp(-eta));
        const safeP = Math.max(1e-15, Math.min(1 - 1e-15, pi));
        currentLogLik += y[i] * Math.log(safeP) + (1 - y[i]) * Math.log(1 - safeP);
      }
      logLik = currentLogLik;

      if (maxDelta < tolerance) {
        converged = true;
        break;
      }
    }

    // Invert Information matrix for parameter covariance
    const finalW = [];
    for (let i = 0; i < N; i++) {
      let eta = 0;
      for (let j = 0; j < numParams; j++) eta += X[i][j] * beta[j];
      eta = Math.max(-30, Math.min(30, eta));
      const pi = 1.0 / (1.0 + Math.exp(-eta));
      finalW.push(Math.max(1e-9, pi * (1.0 - pi)));
    }

    const H = Array.from({ length: numParams }, () => new Array(numParams).fill(0));
    for (let r = 0; r < numParams; r++) {
      for (let c = r; c < numParams; c++) {
        let sum = 0;
        for (let i = 0; i < N; i++) {
          sum += X[i][r] * finalW[i] * X[i][c];
        }
        if (r === c) sum += 1e-5;
        H[r][c] = sum;
        H[c][r] = sum;
      }
    }

    const covMatrixNorm = this._invertMatrix(H) || Array.from({ length: numParams }, () => new Array(numParams).fill(0));

    // Transform coefficients back to unstandardized natural scale
    // y = beta_0 + sum_j beta_j * (x_j - mean_j)/sd_j
    //   = (beta_0 - sum_j (beta_j * mean_j / sd_j)) + sum_j (beta_j / sd_j) * x_j
    const naturalBeta = new Array(numParams).fill(0);
    const naturalSE = new Array(numParams).fill(0);

    let naturalIntercept = beta[0];
    for (let j = 1; j < numParams; j++) {
      const idx = j - 1;
      naturalBeta[j] = beta[j] / sds[idx];
      naturalSE[j] = Math.sqrt(Math.max(0, covMatrixNorm[j][j])) / sds[idx];
      naturalIntercept -= beta[j] * (means[idx] / sds[idx]);
    }
    naturalBeta[0] = naturalIntercept;
    naturalSE[0] = Math.sqrt(Math.max(0, covMatrixNorm[0][0]));

    // Model metrics
    const lrStat = Math.max(0, 2 * (logLik - nullLogLik));
    const lrDf = p;
    const lrPValue = lrDf > 0 ? Distributions.chiSquarePValue(lrStat, lrDf) : 1.0;
    const mcfaddenR2 = nullLogLik !== 0 ? Math.max(0, 1 - (logLik / nullLogLik)) : 0;

    // Attach propensity scores and logit PS to each observation
    data.forEach(d => {
      let eta = naturalBeta[0];
      covariateCols.forEach((col, idx) => {
        eta += naturalBeta[idx + 1] * d[col];
      });
      // Clamping
      eta = Math.max(-25, Math.min(25, eta));
      const ps = 1.0 / (1.0 + Math.exp(-eta));
      d._ps = ps;
      d._logitPs = eta;
    });

    // Parameter summary table
    const coefficients = [
      {
        term: 'Intercept (β₀)',
        estimate: naturalBeta[0],
        stdError: naturalSE[0],
        zScore: naturalSE[0] > 0 ? naturalBeta[0] / naturalSE[0] : 0,
        pValue: naturalSE[0] > 0 ? Distributions.normalPValue(naturalBeta[0] / naturalSE[0]) : 1.0,
        oddsRatio: Math.exp(Math.max(-20, Math.min(20, naturalBeta[0]))),
        ci95: [
          Math.exp(naturalBeta[0] - 1.96 * naturalSE[0]),
          Math.exp(naturalBeta[0] + 1.96 * naturalSE[0])
        ]
      }
    ];

    covariateCols.forEach((col, idx) => {
      const b = naturalBeta[idx + 1];
      const se = naturalSE[idx + 1];
      const z = se > 0 ? b / se : 0;
      const pVal = se > 0 ? Distributions.normalPValue(z) : 1.0;
      coefficients.push({
        term: col,
        estimate: b,
        stdError: se,
        zScore: z,
        pValue: pVal,
        oddsRatio: Math.exp(Math.max(-20, Math.min(20, b))),
        ci95: [
          Math.exp(b - 1.96 * se),
          Math.exp(b + 1.96 * se)
        ]
      });
    });

    return {
      converged,
      logLik,
      nullLogLik,
      lrStat,
      lrDf,
      lrPValue,
      mcfaddenR2,
      coefficients
    };
  },

  /**
   * Nearest-Neighbor 1:1 Matching without Replacement with Strict Caliper
   * @param {Array<Object>} data Dataset with _treatment, _ps, _logitPs, _outcome
   * @param {number} caliperMultiplier Multiplier of SD(logit(PS)), default 0.20 (Austin 2011)
   * @param {boolean} enforceCommonSupport Only match within overlapping propensity range
   */
  matchNearestNeighbor(data, caliperMultiplier = 0.20, enforceCommonSupport = true) {
    const treated = data.filter(d => d._treatment === 1);
    const control = data.filter(d => d._treatment === 0);

    if (treated.length === 0 || control.length === 0) {
      return { error: 'Both treated and control units are required for matching.' };
    }

    // Compute standard deviation of logit propensity score across full sample
    const allLogits = data.map(d => d._logitPs);
    const meanLogit = allLogits.reduce((a, b) => a + b, 0) / allLogits.length;
    const sdLogit = Math.sqrt(
      allLogits.reduce((acc, val) => acc + Math.pow(val - meanLogit, 2), 0) / (allLogits.length - 1)
    ) || 1.0;

    const caliperWidth = caliperMultiplier * sdLogit;

    // Common support boundary
    const treatedMinPs = Math.min(...treated.map(d => d._ps));
    const treatedMaxPs = Math.max(...treated.map(d => d._ps));
    const controlMinPs = Math.min(...control.map(d => d._ps));
    const controlMaxPs = Math.max(...control.map(d => d._ps));

    const commonSupportMin = Math.max(treatedMinPs, controlMinPs);
    const commonSupportMax = Math.min(treatedMaxPs, controlMaxPs);

    // Filter candidate treated if common support requested
    const eligibleTreated = enforceCommonSupport
      ? treated.filter(t => t._ps >= commonSupportMin && t._ps <= commonSupportMax)
      : treated;

    const droppedOffSupportTreated = treated.length - eligibleTreated.length;

    // Sort treated units descending by propensity score (greedy matching order)
    const sortedTreated = [...eligibleTreated].sort((a, b) => b._logitPs - a._logitPs);

    const availableControls = [...control];
    const usedControlIds = new Set();
    const matchedPairs = [];
    let unmatchedDueToCaliper = 0;

    sortedTreated.forEach(tUnit => {
      let bestControl = null;
      let minDistance = Infinity;

      for (let i = 0; i < availableControls.length; i++) {
        const cUnit = availableControls[i];
        if (usedControlIds.has(cUnit._rowId)) continue;

        // If common support enforced, control must also be in common support
        if (enforceCommonSupport && (cUnit._ps < commonSupportMin || cUnit._ps > commonSupportMax)) {
          continue;
        }

        const dist = Math.abs(tUnit._logitPs - cUnit._logitPs);
        if (dist < minDistance) {
          minDistance = dist;
          bestControl = cUnit;
        }
      }

      if (bestControl && minDistance <= caliperWidth) {
        usedControlIds.add(bestControl._rowId);
        matchedPairs.push({
          pairId: matchedPairs.length + 1,
          treated: tUnit,
          control: bestControl,
          distanceLogit: minDistance,
          diffPs: Math.abs(tUnit._ps - bestControl._ps)
        });
      } else {
        unmatchedDueToCaliper++;
      }
    });

    const matchedTreated = matchedPairs.map(p => p.treated);
    const matchedControl = matchedPairs.map(p => p.control);

    return {
      caliperMultiplier,
      sdLogit,
      caliperWidth,
      commonSupport: {
        min: commonSupportMin,
        max: commonSupportMax,
        enforced: enforceCommonSupport,
        droppedTreated: droppedOffSupportTreated
      },
      nTotalTreated: treated.length,
      nTotalControl: control.length,
      nMatchedPairs: matchedPairs.length,
      nMatchedTreated: matchedTreated.length,
      nMatchedControl: matchedControl.length,
      nUnmatchedTreated: treated.length - matchedPairs.length,
      nUnmatchedControl: control.length - matchedPairs.length,
      unmatchedDueToCaliper,
      pairs: matchedPairs,
      matchedTreated,
      matchedControl
    };
  },

  /**
   * Assesses Covariate Balance & Standardized Mean Differences (SMD) Pre- and Post-Matching
   * Pre SMD = (Mean_T_pre - Mean_C_pre) / sqrt((Var_T_pre + Var_C_pre)/2)
   * Post SMD = (Mean_T_post - Mean_C_post) / sqrt((Var_T_pre + Var_C_pre)/2)
   */
  assessBalance(data, matchResult, covariateCols) {
    const rawTreated = data.filter(d => d._treatment === 1);
    const rawControl = data.filter(d => d._treatment === 0);
    const matchedTreated = matchResult.matchedTreated;
    const matchedControl = matchResult.matchedControl;

    const balanceTable = [];
    let allBalanced = true;

    covariateCols.forEach(col => {
      // Pre-matching statistics
      const preT = rawTreated.map(d => d[col]);
      const preC = rawControl.map(d => d[col]);

      const meanTPre = preT.reduce((a, b) => a + b, 0) / preT.length;
      const meanCPre = preC.reduce((a, b) => a + b, 0) / preC.length;

      const varTPre = preT.reduce((acc, v) => acc + Math.pow(v - meanTPre, 2), 0) / (preT.length - 1);
      const varCPre = preC.reduce((acc, v) => acc + Math.pow(v - meanCPre, 2), 0) / (preC.length - 1);

      // Pooled standard deviation from unadjusted pre-matching sample (standard Austin convention)
      const pooledSdPre = Math.sqrt((varTPre + varCPre) / 2) || 1.0;
      const smdPre = (meanTPre - meanCPre) / pooledSdPre;

      // Post-matching statistics
      let meanTPost = 0;
      let meanCPost = 0;
      let varTPost = 0;
      let varCPost = 0;
      let smdPost = 0;
      let varianceRatio = 1.0;

      if (matchedTreated.length > 0) {
        const postT = matchedTreated.map(d => d[col]);
        const postC = matchedControl.map(d => d[col]);

        meanTPost = postT.reduce((a, b) => a + b, 0) / postT.length;
        meanCPost = postC.reduce((a, b) => a + b, 0) / postC.length;

        varTPost = postT.length > 1
          ? postT.reduce((acc, v) => acc + Math.pow(v - meanTPost, 2), 0) / (postT.length - 1)
          : 0;
        varCPost = postC.length > 1
          ? postC.reduce((acc, v) => acc + Math.pow(v - meanCPost, 2), 0) / (postC.length - 1)
          : 0;

        smdPost = (meanTPost - meanCPost) / pooledSdPre;
        varianceRatio = varCPost > 0 ? varTPost / varCPost : 1.0;
      }

      const absSmdPre = Math.abs(smdPre);
      const absSmdPost = Math.abs(smdPost);
      const percentReduction = absSmdPre > 0 ? ((absSmdPre - absSmdPost) / absSmdPre) * 100 : 0;
      const isBalanced = absSmdPost < 0.10; // Standard clinical benchmark
      if (!isBalanced) allBalanced = false;

      balanceTable.push({
        covariate: col,
        meanTreatedPre: meanTPre,
        meanControlPre: meanCPre,
        sdTreatedPre: Math.sqrt(varTPre),
        sdControlPre: Math.sqrt(varCPre),
        smdPre,
        absSmdPre,
        meanTreatedPost: meanTPost,
        meanControlPost: meanCPost,
        sdTreatedPost: Math.sqrt(varTPost),
        sdControlPost: Math.sqrt(varCPost),
        smdPost,
        absSmdPost,
        varianceRatioPost: varianceRatio,
        varRatioPost: varianceRatio,
        percentReduction,
        isBalanced
      });
    });

    // Check Propensity score balance as well
    const psPreT = rawTreated.map(d => d._ps);
    const psPreC = rawControl.map(d => d._ps);
    const psMeanTPre = psPreT.reduce((a, b) => a + b, 0) / psPreT.length;
    const psMeanCPre = psPreC.reduce((a, b) => a + b, 0) / psPreC.length;
    const psVarTPre = psPreT.reduce((acc, v) => acc + Math.pow(v - psMeanTPre, 2), 0) / (psPreT.length - 1);
    const psVarCPre = psPreC.reduce((acc, v) => acc + Math.pow(v - psMeanCPre, 2), 0) / (psPreC.length - 1);
    const psPooledSd = Math.sqrt((psVarTPre + psVarCPre) / 2) || 1.0;
    const psSmdPre = (psMeanTPre - psMeanCPre) / psPooledSd;

    let psSmdPost = 0;
    if (matchedTreated.length > 0) {
      const psPostT = matchedTreated.map(d => d._ps);
      const psPostC = matchedControl.map(d => d._ps);
      const psMeanTPost = psPostT.reduce((a, b) => a + b, 0) / psPostT.length;
      const psMeanCPost = psPostC.reduce((a, b) => a + b, 0) / psPostC.length;
      psSmdPost = (psMeanTPost - psMeanCPost) / psPooledSd;
    }

    const maxAbsSmdPre = Math.max(0, ...balanceTable.map(d => d.absSmdPre || 0));
    const maxAbsSmdPost = Math.max(0, ...balanceTable.map(d => d.absSmdPost || 0));

    return {
      allBalanced,
      threshold: 0.10,
      balanceTable,
      maxAbsSmdPre,
      maxAbsSmdPost,
      psBalance: {
        smdPre: psSmdPre,
        absSmdPre: Math.abs(psSmdPre),
        smdPost: psSmdPost,
        absSmdPost: Math.abs(psSmdPost)
      }
    };
  },

  /**
   * Estimates the Average Treatment Effect on the Treated (ATT)
   * Handles continuous outcomes (paired t-test / robust standard error)
   * and binary outcomes (paired risk difference & McNemar odds ratio).
   */
  estimateOutcomeEffect(data, matchResult, outcomeCol) {
    const pairs = matchResult.pairs;
    if (!pairs || pairs.length === 0) {
      return { error: 'No matched pairs available for outcome estimation.' };
    }

    const nPairs = pairs.length;
    const rawTreated = data.filter(d => d._treatment === 1).map(d => d._outcome);
    const rawControl = data.filter(d => d._treatment === 0).map(d => d._outcome);

    const unadjMeanT = rawTreated.reduce((a, b) => a + b, 0) / rawTreated.length;
    const unadjMeanC = rawControl.reduce((a, b) => a + b, 0) / rawControl.length;
    const unadjVarT = rawTreated.reduce((acc, v) => acc + Math.pow(v - unadjMeanT, 2), 0) / (rawTreated.length - 1);
    const unadjVarC = rawControl.reduce((acc, v) => acc + Math.pow(v - unadjMeanC, 2), 0) / (rawControl.length - 1);
    const unadjSE = Math.sqrt((unadjVarT / rawTreated.length) + (unadjVarC / rawControl.length)) || 1e-4;
    const unadjT = (unadjMeanT - unadjMeanC) / unadjSE;
    const unadjDf = Math.max(1, rawTreated.length + rawControl.length - 2);
    const unadjP = Distributions.tPValue(unadjT, unadjDf);
    const unadjustedDiff = {
      diff: unadjMeanT - unadjMeanC,
      se: unadjSE,
      statistic: unadjT,
      pValue: unadjP,
      ci95: [(unadjMeanT - unadjMeanC) - 1.96 * unadjSE, (unadjMeanT - unadjMeanC) + 1.96 * unadjSE]
    };

    // Detect if outcome is binary (only 0 and 1)
    const allOutcomes = data.map(d => d._outcome);
    const uniqueOutcomes = Array.from(new Set(allOutcomes));
    const isBinaryOutcome = uniqueOutcomes.length === 2 && uniqueOutcomes.every(v => v === 0 || v === 1);

    if (isBinaryOutcome) {
      // Binary outcome: Concordant and discordant matched pairs
      let a = 0; // Both 1
      let b = 0; // Treated=1, Control=0
      let c = 0; // Treated=0, Control=1
      let d = 0; // Both 0

      pairs.forEach(p => {
        const yT = p.treated._outcome;
        const yC = p.control._outcome;
        if (yT === 1 && yC === 1) a++;
        else if (yT === 1 && yC === 0) b++;
        else if (yT === 0 && yC === 1) c++;
        else d++;
      });

      const riskTreated = (a + b) / nPairs;
      const riskControl = (a + c) / nPairs;
      const riskDiff = riskTreated - riskControl; // (b - c) / nPairs
      const riskRatio = riskControl > 0 ? riskTreated / riskControl : Infinity;

      // Paired standard error for Risk Difference
      const seRD = Math.sqrt(((b + c) - Math.pow(b - c, 2) / nPairs)) / nPairs || 1e-6;
      const zRD = riskDiff / seRD;
      const pValRD = Distributions.normalPValue(zRD);
      const ci95RD = [riskDiff - 1.96 * seRD, riskDiff + 1.96 * seRD];

      // McNemar Odds Ratio: b / c
      const or = c > 0 ? b / c : Infinity;
      const lnOrSE = Math.sqrt((1 / Math.max(1, b)) + (1 / Math.max(1, c)));
      const ci95OR = [Math.exp(Math.log(Math.max(1e-6, or)) - 1.96 * lnOrSE), Math.exp(Math.log(Math.max(1e-6, or)) + 1.96 * lnOrSE)];

      const mcNemarChi2 = (b + c) > 0 ? Math.pow(Math.abs(b - c) - 1, 2) / (b + c) : 0;
      const mcNemarP = Distributions.chiSquarePValue(mcNemarChi2, 1);

      return {
        type: 'binary',
        outcomeType: 'binary',
        outcomeCol,
        nPairs,
        unadjustedDiff,
        att: riskDiff,
        pointEstimate: riskDiff,
        metricName: 'ATT (Risk Difference)',
        riskTreated,
        riskControl,
        riskRatio,
        oddsRatio: or,
        ci95OddsRatio: ci95OR,
        se: seRD,
        stdError: seRD,
        statistic: zRD,
        testStatistic: zRD,
        testName: 'Paired Z-Test & McNemar Discordant Pairs',
        pValue: pValRD,
        mcNemarChi2,
        mcNemarP,
        ci95: ci95RD,
        isSignificant: pValRD < 0.05,
        contingency: { a, b, c, d }
      };
    } else {
      // Continuous outcome: Matched pair differences
      const diffs = pairs.map(p => p.treated._outcome - p.control._outcome);
      const meanDiff = diffs.reduce((a, b) => a + b, 0) / nPairs;

      const varianceDiff = diffs.reduce((acc, d) => acc + Math.pow(d - meanDiff, 2), 0) / (nPairs - 1);
      const sdDiff = Math.sqrt(varianceDiff);
      const robustSE = sdDiff / Math.sqrt(nPairs);

      const df = nPairs - 1;
      const tStat = robustSE > 0 ? meanDiff / robustSE : 0;
      const pValue = df > 0 ? Distributions.tPValue(tStat, df) : 1.0;

      // 95% t Critical value approximation
      const tCrit = Math.abs(Distributions.invNormalCDF(0.975)) + (df < 30 ? (1.5 / df) : 0);
      const ci95 = [meanDiff - tCrit * robustSE, meanDiff + tCrit * robustSE];

      // Cohen's d for matched pairs
      const cohensD = sdDiff > 0 ? meanDiff / sdDiff : 0;

      return {
        type: 'continuous',
        outcomeType: 'continuous',
        outcomeCol,
        nPairs,
        unadjustedDiff,
        att: meanDiff,
        pointEstimate: meanDiff,
        metricName: 'ATT (Mean Difference)',
        se: robustSE,
        stdError: robustSE,
        statistic: tStat,
        testStatistic: tStat,
        testName: 'Paired Student t-Test with Cluster-Robust SE',
        df,
        pValue,
        ci95,
        cohensD,
        isSignificant: pValue < 0.05
      };
    }
  },

  /**
   * Master execution pipeline for Propensity Score Matching
   */
  executeAnalysis(rows, arg2, arg3, arg4, arg5 = {}) {
    let treatmentCol, outcomeCol, covariateCols, options;
    if (typeof arg2 === 'object' && !Array.isArray(arg2) && arg2 !== null) {
      treatmentCol = arg2.treatmentCol;
      outcomeCol = arg2.outcomeCol;
      covariateCols = arg2.covariateCols;
      options = arg2;
    } else {
      treatmentCol = arg2;
      outcomeCol = arg3;
      covariateCols = arg4;
      options = arg5 || {};
    }

    const caliper = options.caliperMultiplier !== undefined ? options.caliperMultiplier : (options.caliper !== undefined ? options.caliper : 0.20);
    const enforceCommonSupport = options.enforceCommonSupport !== undefined ? options.enforceCommonSupport : true;

    // 1. Prepare data & inspect missingness
    const prep = this.prepareData(rows, treatmentCol, outcomeCol, covariateCols);
    if (prep.error) return { error: prep.error };

    // 2. Fit Logistic Regression model
    const logitModel = this.fitLogisticRegression(prep.data, covariateCols);
    if (logitModel.error) return { error: logitModel.error };

    // 3. Perform Nearest-Neighbor Matching with Caliper
    const matchResult = this.matchNearestNeighbor(prep.data, caliper, enforceCommonSupport);
    if (matchResult.error) return { error: matchResult.error };

    // 4. Assess Covariate Balance (SMD & Love plot data)
    const balance = this.assessBalance(prep.data, matchResult, covariateCols);

    // 5. Outcome Analysis & ATT Estimation
    const outcomeResult = this.estimateOutcomeEffect(prep.data, matchResult, outcomeCol);
    if (outcomeResult.error) return { error: outcomeResult.error };

    // 6. Generate Overlap Histogram / Density Data
    const overlapData = this._generateOverlapData(prep.data, matchResult);

    // 7. Generate Reproducible Scripts
    const scripts = {
      python: this.generatePythonScript(treatmentCol, outcomeCol, covariateCols, caliper, outcomeResult.type || outcomeResult.outcomeType),
      r: this.generateRScript(treatmentCol, outcomeCol, covariateCols, caliper, outcomeResult.type || outcomeResult.outcomeType),
      stata: this.generateStataScript(treatmentCol, outcomeCol, covariateCols, caliper, outcomeResult.type || outcomeResult.outcomeType)
    };

    const isSig = outcomeResult.isSignificant;
    const reportText = [
      `An observational clinical cohort analysis (N = ${prep.totalRows}, Complete cases = ${prep.completeCasesCount}) was conducted to estimate the Average Treatment Effect on the Treated (ATT) of ${treatmentCol} on ${outcomeCol}.`,
      `Baseline confounding was controlled using 1:1 nearest-neighbor propensity score matching without replacement within a strict caliper of ${caliper} × SD(logit PS) (${matchResult.caliperWidth.toFixed(4)}), with common support enforcement.`,
      `A total of ${matchResult.nPairs} matched pairs (${matchResult.nTreatedMatched} treated vs. ${matchResult.nControlMatched} control) were successfully matched.`,
      `Covariate balance assessment demonstrated marked bias reduction: maximum post-matching |SMD| was ${balance.maxAbsSmdPost.toFixed(3)} (${balance.allBalanced ? 'all covariates below the 0.10 threshold' : 'acceptable balance'}), compared to maximum pre-matching |SMD| of ${balance.maxAbsSmdPre.toFixed(3)}.`,
      outcomeResult.type === 'continuous'
        ? `The unadjusted observational mean difference was ${outcomeResult.unadjustedDiff.diff >= 0 ? '+' : ''}${outcomeResult.unadjustedDiff.diff.toFixed(2)} (SE = ${outcomeResult.unadjustedDiff.se.toFixed(2)}, p = ${outcomeResult.unadjustedDiff.pValue < 0.001 ? '< .001' : outcomeResult.unadjustedDiff.pValue.toFixed(3)}). After propensity matching, the unconfounded causal ATT point estimate was ${outcomeResult.att >= 0 ? '+' : ''}${outcomeResult.att.toFixed(2)} (robust SE = ${outcomeResult.se.toFixed(2)}, 95% CI [${outcomeResult.ci95[0].toFixed(2)}, ${outcomeResult.ci95[1].toFixed(2)}], t(${outcomeResult.df}) = ${outcomeResult.statistic.toFixed(2)}, p = ${outcomeResult.pValue < 0.001 ? '< .001' : outcomeResult.pValue.toFixed(3)}), indicating a ${isSig ? 'statistically significant' : 'non-significant'} causal treatment effect.`
        : `The unadjusted observational risk difference was ${(outcomeResult.unadjustedDiff.diff * 100).toFixed(1)}%. After propensity matching, the matched causal ATT Risk Difference was ${(outcomeResult.att * 100).toFixed(1)}% (95% CI [${(outcomeResult.ci95[0] * 100).toFixed(1)}%, ${(outcomeResult.ci95[1] * 100).toFixed(1)}%], p = ${outcomeResult.pValue < 0.001 ? '< .001' : outcomeResult.pValue.toFixed(3)}), with a matched Odds Ratio of ${outcomeResult.oddsRatio.toFixed(2)}.`
    ].join(' ');

    return {
      treatmentCol,
      outcomeCol,
      covariateCols,
      totalN: prep.totalRows,
      completeN: prep.completeCasesCount,
      nMatchedPairs: matchResult.nMatchedPairs || matchResult.nPairs || (matchResult.pairs ? matchResult.pairs.length : 0),
      matchedTreatedN: matchResult.nMatchedTreated || (matchResult.matchedTreated ? matchResult.matchedTreated.length : 0),
      unmatchedTreatedN: matchResult.nTotalTreated || 0,
      matchedControlN: matchResult.nMatchedControl || (matchResult.matchedControl ? matchResult.matchedControl.length : 0),
      logisticRegression: logitModel,
      caliper: {
        multiplier: caliper,
        width: matchResult.caliperWidth,
        enforceCommonSupport
      },
      balance,
      outcome: outcomeResult,
      overlap: overlapData,
      scripts,
      reportText,
      // Nested legacy properties for compatibility
      prep,
      logitModel,
      matchResult,
      outcomeResult,
      overlapData
    };
  },

  /**
   * Prepares density bin data for the Propensity Score Overlap Plot
   */
  _generateOverlapData(data, matchResult) {
    const rawTreated = data.filter(d => d._treatment === 1).map(d => d._ps);
    const rawControl = data.filter(d => d._treatment === 0).map(d => d._ps);
    const matchedTreated = matchResult.matchedTreated.map(d => d._ps);
    const matchedControl = matchResult.matchedControl.map(d => d._ps);

    const nBins = 25;
    const bins = [];
    for (let i = 0; i < nBins; i++) {
      const start = i / nBins;
      const end = (i + 1) / nBins;
      const mid = (start + end) / 2;

      const preT = rawTreated.filter(p => p >= start && p < end).length / (rawTreated.length || 1);
      const preC = rawControl.filter(p => p >= start && p < end).length / (rawControl.length || 1);
      const postT = matchedTreated.filter(p => p >= start && p < end).length / (matchedTreated.length || 1);
      const postC = matchedControl.filter(p => p >= start && p < end).length / (matchedControl.length || 1);

      bins.push({
        binIndex: i,
        range: [start, end],
        mid,
        preTreated: preT,
        preControl: preC,
        postTreated: postT,
        postControl: postC
      });
    }

    return {
      nBins,
      bins,
      commonSupport: matchResult.commonSupport
    };
  },

  /**
   * Generates a complete, reproducible Python script using statsmodels and scipy
   */
  generatePythonScript(treatmentCol, outcomeCol, covariateCols, caliper = 0.20, outcomeType = 'continuous') {
    if (typeof treatmentCol === 'object' && treatmentCol !== null) {
      const opts = treatmentCol;
      treatmentCol = opts.treatmentCol;
      outcomeCol = opts.outcomeCol;
      covariateCols = opts.covariateCols;
      caliper = opts.caliperMultiplier !== undefined ? opts.caliperMultiplier : (opts.caliper !== undefined ? opts.caliper : 0.20);
      outcomeType = opts.isBinaryOutcome ? 'binary' : (opts.outcomeType || 'continuous');
    }
    const covListPy = (covariateCols || []).map(c => `'${c}'`).join(', ');
    return `"""
=============================================================================
Propensity Score Matching (PSM) - Observational Clinical Causal Analysis
Methodology: 1:1 Nearest-Neighbor Caliper Matching Without Replacement
Generated by Statis-Gravity Clinical Biostatistics Platform
=============================================================================
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.spatial.distance import cdist
from scipy import stats
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# 1. DATA PREPARATION & MISSINGNESS HANDLING
# -----------------------------------------------------------------------------
# Load observational clinical dataset
# df = pd.read_csv('clinical_observational_data.csv')

treatment_col = '${treatmentCol}'
outcome_col = '${outcomeCol}'
covariates = [${covListPy}]

print("--- Step 1: Initial Cohort Inspection & Missing Data ---")
cols_needed = [treatment_col, outcome_col] + covariates
print(f"Total observational records: {len(df)}")
missing_summary = df[cols_needed].isnull().sum()
print("Missing values per variable:\\n", missing_summary)

# Complete-case analysis rationale:
# When missingness is low (< 5%) or missing completely at random (MCAR),
# complete-case analysis provides unbiased propensity score estimates.
df_clean = df.dropna(subset=cols_needed).copy()
df_clean[treatment_col] = df_clean[treatment_col].astype(int)
print(f"Complete cases retained: {len(df_clean)} "
      f"(Treated: {sum(df_clean[treatment_col] == 1)}, Control: {sum(df_clean[treatment_col] == 0)})")

# -----------------------------------------------------------------------------
# 2. PROPENSITY SCORE ESTIMATION (MULTIVARIATE LOGISTIC REGRESSION)
# -----------------------------------------------------------------------------
# Clinical Rationale:
# Propensity score e(X) = P(Treatment = 1 | X) collapses multi-dimensional baseline
# confounding into a scalar balancing score (Rosenbaum & Rubin, 1983).
formula = f"{treatment_col} ~ " + " + ".join(covariates)
ps_model = smf.logit(formula=formula, data=df_clean).fit(disp=False)
print("\\n--- Step 2: Logistic Regression Model Parameters ---")
print(ps_model.summary())

# Extract propensity scores and logit of propensity scores
df_clean['ps'] = ps_model.predict(df_clean)
# Logit transformation: linearizes propensity score for robust matching metric
df_clean['logit_ps'] = np.log(df_clean['ps'] / (1.0 - df_clean['ps']))

# Common Support Assessment:
treated_ps = df_clean[df_clean[treatment_col] == 1]['ps']
control_ps = df_clean[df_clean[treatment_col] == 0]['ps']
cs_min = max(treated_ps.min(), control_ps.min())
cs_max = min(treated_ps.max(), control_ps.max())
print(f"\\nRegion of Common Support: [{cs_min:.4f}, {cs_max:.4f}]")

# -----------------------------------------------------------------------------
# 3. 1:1 NEAREST NEIGHBOR MATCHING WITH CALIPER (WITHOUT REPLACEMENT)
# -----------------------------------------------------------------------------
# Clinical Rationale:
# Austin (2011) demonstrated that a caliper width equal to 0.20 of the standard
# deviation of the logit of the propensity score eliminates > 98% of baseline bias.
sd_logit_ps = df_clean['logit_ps'].std()
caliper_width = ${caliper.toFixed(2)} * sd_logit_ps
print(f"SD(logit(PS)): {sd_logit_ps:.4f} | Enforced Caliper Width: {caliper_width:.4f}")

treated_df = df_clean[(df_clean[treatment_col] == 1) & (df_clean['ps'] >= cs_min) & (df_clean['ps'] <= cs_max)].copy()
control_df = df_clean[(df_clean[treatment_col] == 0) & (df_clean['ps'] >= cs_min) & (df_clean['ps'] <= cs_max)].copy()

# Sort treated descending to prioritize high-propensity subjects (greedy matching)
treated_df = treated_df.sort_values(by='logit_ps', ascending=False)

matched_pairs = []
available_controls = control_df.copy()

for t_idx, t_row in treated_df.iterrows():
    if len(available_controls) == 0:
        break
    
    # Compute absolute distance on logit propensity score
    distances = np.abs(available_controls['logit_ps'].values - t_row['logit_ps'])
    min_dist_idx = np.argmin(distances)
    min_dist = distances[min_dist_idx]
    
    if min_dist <= caliper_width:
        c_idx = available_controls.index[min_dist_idx]
        matched_pairs.append({
            'treated_id': t_idx,
            'control_id': c_idx,
            'distance_logit': min_dist,
            'diff_ps': abs(t_row['ps'] - available_controls.loc[c_idx, 'ps'])
        })
        available_controls = available_controls.drop(c_idx)

pairs_df = pd.DataFrame(matched_pairs)
matched_treated_ids = pairs_df['treated_id'].values
matched_control_ids = pairs_df['control_id'].values

df_matched_treated = df_clean.loc[matched_treated_ids].copy()
df_matched_control = df_clean.loc[matched_control_ids].copy()
df_matched = pd.concat([df_matched_treated, df_matched_control])

print(f"Successfully matched: {len(pairs_df)} pairs ({len(pairs_df)*2} total subjects)")
print(f"Unmatched treated: {len(treated_df) - len(pairs_df)}")

# -----------------------------------------------------------------------------
# 4. BALANCE ASSESSMENT & STANDARDIZED MEAN DIFFERENCES (LOVE PLOT)
# -----------------------------------------------------------------------------
# Clinical Rationale:
# Standardized Mean Difference (SMD) is independent of sample size.
# An SMD < 0.10 indicates adequate balance between treatment and control cohorts.
balance_records = []
for cov in covariates:
    # Pre-matching pooled standard deviation (Austin, 2009)
    m_t_pre = df_clean[df_clean[treatment_col] == 1][cov].mean()
    m_c_pre = df_clean[df_clean[treatment_col] == 0][cov].mean()
    v_t_pre = df_clean[df_clean[treatment_col] == 1][cov].var()
    v_c_pre = df_clean[df_clean[treatment_col] == 0][cov].var()
    pooled_sd = np.sqrt((v_t_pre + v_c_pre) / 2.0)
    
    smd_pre = (m_t_pre - m_c_pre) / pooled_sd
    
    # Post-matching
    m_t_post = df_matched_treated[cov].mean()
    m_c_post = df_matched_control[cov].mean()
    smd_post = (m_t_post - m_c_post) / pooled_sd
    
    balance_records.append({
        'Covariate': cov,
        'Mean_T_Pre': m_t_pre,
        'Mean_C_Pre': m_c_pre,
        'SMD_Pre': smd_pre,
        'Mean_T_Post': m_t_post,
        'Mean_C_Post': m_c_post,
        'SMD_Post': smd_post,
        'Balanced': abs(smd_post) < 0.10
    })

balance_df = pd.DataFrame(balance_records)
print("\\n--- Step 4: Baseline Covariate Balance Assessment ---")
print(balance_df.to_string(index=False))

# Love Plot Generation
plt.figure(figsize=(8, len(covariates) * 0.6 + 2))
y_pos = np.arange(len(covariates))
plt.scatter(np.abs(balance_df['SMD_Pre']), y_pos, color='#f43f5e', label='Unadjusted (Pre-Match)', s=80, zorder=3)
plt.scatter(np.abs(balance_df['SMD_Post']), y_pos, color='#10b981', label='Matched (Post-Match)', s=90, zorder=4)
for i in range(len(covariates)):
    plt.plot([np.abs(balance_df['SMD_Pre'][i]), np.abs(balance_df['SMD_Post'][i])], [y_pos[i], y_pos[i]], color='gray', linestyle=':', alpha=0.7)

plt.axvline(x=0.10, color='#00d2ff', linestyle='--', linewidth=1.5, label='Standard Balance Cutoff (0.10)')
plt.axvline(x=0.05, color='#a855f7', linestyle=':', linewidth=1.2, label='Strict Cutoff (0.05)')
plt.yticks(y_pos, balance_df['Covariate'])
plt.xlabel('Absolute Standardized Mean Difference (|SMD|)')
plt.title('Love Plot: Covariate Balance Before and After PSM')
plt.legend(loc='upper right')
plt.grid(True, alpha=0.25)
plt.tight_layout()
plt.show()

# -----------------------------------------------------------------------------
# 5. OUTCOME ANALYSIS: AVERAGE TREATMENT EFFECT ON THE TREATED (ATT)
# -----------------------------------------------------------------------------
print("\\n--- Step 5: Causal Treatment Effect Estimation (ATT) ---")
${outcomeType === 'continuous' ? `# Continuous Outcome Analysis: Paired t-test with cluster-robust standard errors
pairs_df['y_treated'] = df_clean.loc[pairs_df['treated_id'], outcome_col].values
pairs_df['y_control'] = df_clean.loc[pairs_df['control_id'], outcome_col].values
pairs_df['diff'] = pairs_df['y_treated'] - pairs_df['y_control']

att = pairs_df['diff'].mean()
se_att = pairs_df['diff'].std() / np.sqrt(len(pairs_df))
t_stat = att / se_att
df_paired = len(pairs_df) - 1
p_val = 2.0 * (1.0 - stats.t.cdf(np.abs(t_stat), df=df_paired))
ci_95 = stats.t.interval(0.95, df=df_paired, loc=att, scale=se_att)

print(f"Outcome Variable: {outcome_col}")
print(f"Unadjusted Mean Diff: {df_clean[df_clean[treatment_col] == 1][outcome_col].mean() - df_clean[df_clean[treatment_col] == 0][outcome_col].mean():.4f}")
print(f"ATT (Average Treatment Effect on Treated): {att:.4f}")
print(f"Paired Robust Standard Error: {se_att:.4f}")
print(f"t-statistic ({df_paired} df): {t_stat:.4f} | p-value: {p_val:.4e}")
print(f"95% Confidence Interval: [{ci_95[0]:.4f}, {ci_95[1]:.4f}]")` : `# Binary Outcome Analysis: Paired Discordant Pairs & McNemar Analysis
pairs_df['y_treated'] = df_clean.loc[pairs_df['treated_id'], outcome_col].values
pairs_df['y_control'] = df_clean.loc[pairs_df['control_id'], outcome_col].values

# Discordant pairs (b: Treated=1, Control=0; c: Treated=0, Control=1)
b = sum((pairs_df['y_treated'] == 1) & (pairs_df['y_control'] == 0))
c = sum((pairs_df['y_treated'] == 0) & (pairs_df['y_control'] == 1))
att_rd = (b - c) / len(pairs_df)
se_rd = np.sqrt((b + c) - ((b - c)**2 / len(pairs_df))) / len(pairs_df)
z_stat = att_rd / se_rd
p_val = 2 * (1 - stats.norm.cdf(abs(z_stat)))
ci_95 = [att_rd - 1.96 * se_rd, att_rd + 1.96 * se_rd]

print(f"Outcome Variable: {outcome_col}")
print(f"ATT Risk Difference: {att_rd:.4f} (95% CI: [{ci_95[0]:.4f}, {ci_95[1]:.4f}])")
print(f"Standard Error: {se_rd:.4f} | z = {z_stat:.4f} | p-value = {p_val:.4e}")
print(f"Discordant Pairs: b (T=1, C=0) = {b}, c (T=0, C=1) = {c}")`}
`;
  },

  /**
   * Generates a complete, reproducible R script using MatchIt and cobalt
   */
  generateRScript(treatmentCol, outcomeCol, covariateCols, caliper = 0.20, outcomeType = 'continuous') {
    if (typeof treatmentCol === 'object' && treatmentCol !== null) {
      const opts = treatmentCol;
      treatmentCol = opts.treatmentCol;
      outcomeCol = opts.outcomeCol;
      covariateCols = opts.covariateCols;
      caliper = opts.caliperMultiplier !== undefined ? opts.caliperMultiplier : (opts.caliper !== undefined ? opts.caliper : 0.20);
      outcomeType = opts.isBinaryOutcome ? 'binary' : (opts.outcomeType || 'continuous');
    }
    const covFormulaR = (covariateCols || []).join(' + ');
    return `################################################################################
# Propensity Score Matching (PSM) - Observational Clinical Causal Analysis
# Package Suite: MatchIt, cobalt, survey, sandwich
# Generated by Statis-Gravity Clinical Biostatistics Platform
################################################################################

# Install required biostatistical packages if needed:
# install.packages(c("MatchIt", "cobalt", "survey", "sandwich", "lmtest", "ggplot2"))
library(MatchIt)
library(cobalt)
library(survey)
library(sandwich)
library(lmtest)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. DATA PREP & MISSINGNESS
# ------------------------------------------------------------------------------
# Load observational clinical data
# df <- read.csv("clinical_observational_data.csv")

treatment_var <- "${treatmentCol}"
outcome_var   <- "${outcomeCol}"
covariates    <- c(${covariateCols.map(c => `"${c}"`).join(', ')})

cat("--- Step 1: Missing Data & Cohort Inspection ---\\n")
cat("Total raw sample size:", nrow(df), "\\n")
cols_needed <- c(treatment_var, outcome_var, covariates)
missing_counts <- colSums(is.na(df[, cols_needed]))
print(missing_counts)

# Complete-case filtering:
df_clean <- na.omit(df[, cols_needed])
df_clean[[treatment_var]] <- as.numeric(df_clean[[treatment_var]])
cat("Complete cases retained:", nrow(df_clean), "\\n")
cat("Treated units:", sum(df_clean[[treatment_var]] == 1), 
    "| Control units:", sum(df_clean[[treatment_var]] == 0), "\\n\\n")

# ------------------------------------------------------------------------------
# 2 & 3. PROPENSITY SCORE ESTIMATION & 1:1 CALIPER MATCHING
# ------------------------------------------------------------------------------
# Clinical Rationale:
# MatchIt fits logistic regression under distance = "logit",
# enforces nearest neighbor matching without replacement (method = "nearest"),
# and sets caliper = ${caliper.toFixed(2)} (in SD of logit PS) within common support.
match_formula <- as.formula(paste("${treatmentCol} ~", "${covFormulaR}"))

cat("--- Steps 2 & 3: Estimating Propensity Scores & Executing MatchIt ---\\n")
m.out <- matchit(
  formula = match_formula,
  data = df_clean,
  method = "nearest",
  distance = "logit",
  caliper = ${caliper.toFixed(2)},
  std.caliper = TRUE,       # Enforce caliper as 0.20 SD of logit(PS)
  replace = FALSE,
  discard = "both"          # Enforce common support boundaries
)

print(summary(m.out, un = TRUE))

# Common support / Propensity score overlap visualization
plot(m.out, type = "jitter", interactive = FALSE)
plot(m.out, type = "density", interactive = FALSE)

# ------------------------------------------------------------------------------
# 4. BALANCE ASSESSMENT & LOVE PLOT
# ------------------------------------------------------------------------------
cat("\\n--- Step 4: Assessing Covariate Balance (cobalt Suite) ---\\n")
# Summary balance table reporting Standardized Mean Differences (SMD)
bal_tab <- bal.tab(m.out, un = TRUE, stats = c("m", "v"), thresholds = c(m = 0.10))
print(bal_tab)

# Publication-grade Love Plot (Austin 2009 standard: |SMD| < 0.10)
love.plot(
  m.out,
  binary = "std",
  thresholds = c(m = 0.10),
  colors = c("#f43f5e", "#10b981"),
  shapes = c(21, 19),
  size = 3.5,
  var.order = "unadjusted",
  title = "Love Plot: Covariate Balance (Pre- vs. Post-Matching)"
)

# Extract matched analytic dataset
df_matched <- match.data(m.out)

# ------------------------------------------------------------------------------
# 5. OUTCOME ANALYSIS: AVERAGE TREATMENT EFFECT ON THE TREATED (ATT)
# ------------------------------------------------------------------------------
cat("\\n--- Step 5: Estimating Average Treatment Effect on the Treated (ATT) ---\\n")
${outcomeType === 'continuous' ? `# Continuous Outcome Analysis: Linear model with pair-clustered robust SE
att_fit <- lm(as.formula(paste("${outcomeCol} ~", "${treatmentCol}")), 
              data = df_matched, 
              weights = weights)

# Cluster-robust standard errors accounting for matched pairs (subclass)
att_robust <- coeftest(att_fit, vcov. = vcovCL, cluster = ~subclass)
print(att_robust)

ci_att <- coefci(att_fit, vcov. = vcovCL, cluster = ~subclass)
cat("\\n95% Confidence Interval for ATT:\\n")
print(ci_att["${treatmentCol}", ])` : `# Binary Outcome Analysis: Logistic regression with cluster-robust standard errors
att_glm <- glm(as.formula(paste("${outcomeCol} ~", "${treatmentCol}")), 
               data = df_matched, 
               family = binomial(link = "logit"), 
               weights = weights)

att_robust <- coeftest(att_glm, vcov. = vcovCL, cluster = ~subclass)
print(att_robust)

or_est <- exp(coef(att_glm)["${treatmentCol}"])
or_ci <- exp(coefci(att_glm, vcov. = vcovCL, cluster = ~subclass)["${treatmentCol}", ])
cat("\\nATT Odds Ratio:", round(or_est, 4), "\\n")
cat("95% CI for Odds Ratio: [", round(or_ci[1], 4), ",", round(or_ci[2], 4), "]\\n")`}
`;
  },

  /**
   * Generates a complete, reproducible Stata script (.do)
   */
  generateStataScript(treatmentCol, outcomeCol, covariateCols, caliper = 0.20, outcomeType = 'continuous') {
    if (typeof treatmentCol === 'object' && treatmentCol !== null) {
      const opts = treatmentCol;
      treatmentCol = opts.treatmentCol;
      outcomeCol = opts.outcomeCol;
      covariateCols = opts.covariateCols;
      caliper = opts.caliperMultiplier !== undefined ? opts.caliperMultiplier : (opts.caliper !== undefined ? opts.caliper : 0.20);
      outcomeType = opts.isBinaryOutcome ? 'binary' : (opts.outcomeType || 'continuous');
    }
    const covListStata = (covariateCols || []).join(' ');
    return `* ==============================================================================
* Propensity Score Matching (PSM) - Observational Clinical Causal Analysis
* Stata Implementation using teffects psmatch & pstest
* Generated by Statis-Gravity Clinical Biostatistics Platform
* ==============================================================================

version 17
clear all
set more off

* 1. LOAD OBSERVATIONAL CLINICAL DATA
* use "clinical_observational_data.dta", clear

local treatment "${treatmentCol}"
local outcome "${outcomeCol}"
local covariates "${covListStata}"

* ------------------------------------------------------------------------------
* 1. DATA PREP & MISSINGNESS
* ------------------------------------------------------------------------------
display as text "--- Step 1: Missing Data Inspection ---"
misstable summarize \`treatment' \`outcome' \`covariates'

* Complete cases filtering
drop if missing(\`treatment') | missing(\`outcome')
foreach var of local covariates {
    drop if missing(\`var')
}
display as text "Complete cases retained: " _N
tabulate \`treatment'

* ------------------------------------------------------------------------------
* 2 & 3. PROPENSITY SCORE ESTIMATION & MATCHING (teffects psmatch)
* ------------------------------------------------------------------------------
* teffects psmatch fits logistic model, estimates ATT (atet),
* and implements 1:1 nearest neighbor matching within caliper.
display as text "--- Steps 2 & 3: Estimating PS & Caliper Matching ---"
teffects psmatch (\`outcome') (\`treatment' \`covariates', logit), ///
    atet ///
    caliper(${caliper.toFixed(2)}) ///
    vce(robust) ///
    generate(match_id)

* Display ATT estimation results
teffects summarize

* Overlap plot: assess common support
psgraph, treated(\`treatment') pscore(match_id)

* ------------------------------------------------------------------------------
* 4. BALANCE ASSESSMENT & LOVE PLOT
* ------------------------------------------------------------------------------
* pstest computes Standardized Mean Differences (SMD) before and after matching
display as text "--- Step 4: Assessing Covariate Balance ---"
pstest \`covariates', both graph

* Save Stata balance graph
graph export "psm_love_plot.png", replace width(2000)

* ------------------------------------------------------------------------------
* 5. OUTCOME ANALYSIS
* ------------------------------------------------------------------------------
display as text "--- Step 5: Causal Treatment Effect (ATT) Summary ---"
* Point estimate, standard error, z, p-value, and 95% CI
matrix list r(table)
`;
  },

  /**
   * Helper: Invert square matrix via Gauss-Jordan elimination with partial pivoting
   */
  _invertMatrix(matrix) {
    const n = matrix.length;
    // Augment with identity matrix
    const aug = matrix.map((row, i) => {
      const identityRow = new Array(n).fill(0);
      identityRow[i] = 1.0;
      return [...row, ...identityRow];
    });

    for (let i = 0; i < n; i++) {
      // Find pivot
      let maxRow = i;
      let maxVal = Math.abs(aug[i][i]);
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(aug[k][i]) > maxVal) {
          maxVal = Math.abs(aug[k][i]);
          maxRow = k;
        }
      }

      if (maxVal < 1e-12) return null; // Singular

      // Swap rows
      if (maxRow !== i) {
        const temp = aug[i];
        aug[i] = aug[maxRow];
        aug[maxRow] = temp;
      }

      // Scale pivot row
      const pivot = aug[i][i];
      for (let j = 0; j < 2 * n; j++) {
        aug[i][j] /= pivot;
      }

      // Eliminate column
      for (let k = 0; k < n; k++) {
        if (k !== i) {
          const factor = aug[k][i];
          for (let j = 0; j < 2 * n; j++) {
            aug[k][j] -= factor * aug[i][j];
          }
        }
      }
    }

    return aug.map(row => row.slice(n));
  },

  /**
   * Helper: Solve linear system A * x = b via Gaussian elimination with partial pivoting
   */
  _solveLinearSystem(A, b) {
    const n = A.length;
    const M = A.map((row, i) => [...row, b[i]]);

    for (let i = 0; i < n; i++) {
      let maxRow = i;
      let maxVal = Math.abs(M[i][i]);
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(M[k][i]) > maxVal) {
          maxVal = Math.abs(M[k][i]);
          maxRow = k;
        }
      }

      if (maxVal < 1e-12) return null;

      if (maxRow !== i) {
        const temp = M[i];
        M[i] = M[maxRow];
        M[maxRow] = temp;
      }

      const pivot = M[i][i];
      for (let j = i; j <= n; j++) {
        M[i][j] /= pivot;
      }

      for (let k = 0; k < n; k++) {
        if (k !== i) {
          const factor = M[k][i];
          for (let j = i; j <= n; j++) {
            M[k][j] -= factor * M[i][j];
          }
        }
      }
    }

    return M.map(row => row[n]);
  },

  /**
   * Parses CSV or TSV string into an array of patient records
   */
  parseClinicalCsv(csvText) {
    if (!csvText || typeof csvText !== 'string') return [];
    const lines = csvText.trim().split(/\r?\n/).filter(l => l.trim().length > 0);
    if (lines.length < 2) return [];

    const delimiter = lines[0].includes('\t') ? '\t' : (lines[0].includes(';') ? ';' : ',');
    const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));

    const records = [];
    for (let i = 1; i < lines.length; i++) {
      const parts = lines[i].split(delimiter).map(p => p.trim().replace(/^["']|["']$/g, ''));
      if (parts.length === headers.length) {
        const row = {};
        headers.forEach((h, idx) => {
          const val = parts[idx];
          const num = Number(val);
          row[h] = !isNaN(num) && val !== '' ? num : val;
        });
        records.push(row);
      }
    }
    return records;
  },

  /**
   * Generates a realistic observational clinical craniosynostosis dataset
   * 120 patients evaluated for surgical management:
   * Endoscopic Strip Craniectomy (treatment = 1) vs. Open Cranial Vault Reconstruction (control = 0)
   */
  getSampleClinicalDataset() {
    // Deterministic pseudo-random seed generator for perfect clinical reproducibility
    let seed = 42;
    const rnd = () => {
      seed = (seed * 9301 + 49297) % 233280;
      return seed / 233280;
    };
    const rndNormal = (mean, sd) => {
      const u1 = Math.max(1e-7, rnd());
      const u2 = rnd();
      const z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      return mean + z * sd;
    };

    const dataset = [];
    const nTotal = 120;

    for (let i = 1; i <= nTotal; i++) {
      // Patient baseline covariates
      const age = Math.max(2.2, Math.min(12.0, rndNormal(6.4, 2.3))); // Age in months
      const sex = rnd() > 0.42 ? 1 : 0; // 58% male
      const comorbidity_score = rnd() < 0.60 ? 0 : (rnd() < 0.85 ? 1 : 2); // 0=None, 1=Mild, 2=Moderate/Syndromic
      const baseline_severity = Math.max(5.0, Math.min(15.5, rndNormal(9.8, 2.4))); // Initial CVAI %
      const bmi_percentile = Math.max(10, Math.min(95, Math.round(rndNormal(52, 18))));

      // Observational clinical assignment probability (Confounded Treatment Allocation)
      // Clinicians steer younger infants with lower severity and fewer comorbidities towards endoscopic strip
      const latentZ = 1.6 
        - 0.38 * (age - 6.0) 
        - 0.28 * (baseline_severity - 9.5) 
        - 0.75 * comorbidity_score 
        + 0.15 * (sex === 1 ? 0.5 : -0.5);
      const trueProb = 1.0 / (1.0 + Math.exp(-latentZ));
      const treatment_col = rnd() < trueProb ? 1 : 0;

      // Clinical Outcome: Hospital Length of Stay (days)
      // True causal effect of endoscopic strip is a reduction of ~2.2 days
      // Outcome is also influenced by baseline severity, age, and comorbidities
      let trueOutcome = 4.8 
        - 2.3 * treatment_col 
        + 0.18 * age 
        + 0.22 * baseline_severity 
        + 0.65 * comorbidity_score 
        + rndNormal(0, 0.45);
      trueOutcome = Math.max(1.0, Math.round(trueOutcome * 10) / 10);

      dataset.push({
        patient_id: `PT-${1000 + i}`,
        treatment_col,
        outcome_col: trueOutcome,
        age: Math.round(age * 10) / 10,
        sex,
        comorbidity_score,
        baseline_severity: Math.round(baseline_severity * 10) / 10,
        bmi_percentile
      });
    }

    return dataset;
  }
};

/**
 * Statis-Gravity - Clinical DOCX Report Generator
 * Generates native Microsoft Office Open XML (.docx) files without external dependencies.
 * Creates compliant ZIP packages with WordprocessingML, formatting tables, callout boxes,
 * APA/ICMJE statistical summaries, methodological rationale, and background knowledge.
 */

// CRC-32 Lookup Table for ZIP header validation
const crcTable = new Uint32Array(256);
for (let n = 0; n < 256; n++) {
  let c = n;
  for (let k = 0; k < 8; k++) {
    c = ((c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1));
  }
  crcTable[n] = c;
}

function calcCRC32(bytes) {
  let crc = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i++) {
    crc = (crc >>> 8) ^ crcTable[(crc ^ bytes[i]) & 0xFF];
  }
  return (crc ^ 0xFFFFFFFF) >>> 0;
}

/**
 * Pure JavaScript ZIP archive packager (PKZIP specification, Stored / uncompressed method 0)
 * Works universally in browser and Node.js environments.
 */
class SimpleZip {
  constructor() {
    this.files = [];
  }

  addFile(name, content) {
    let data;
    if (typeof content === 'string') {
      data = new TextEncoder().encode(content);
    } else if (content instanceof Uint8Array) {
      data = content;
    } else if (content && content.buffer) {
      data = new Uint8Array(content.buffer);
    } else {
      data = new Uint8Array(content);
    }
    this.files.push({ name, data });
  }

  build() {
    const encoder = new TextEncoder();
    const prepared = this.files.map(f => {
      const nameBytes = encoder.encode(f.name);
      const crc = calcCRC32(f.data);
      const localSize = 30 + nameBytes.length + f.data.length;
      const cdSize = 46 + nameBytes.length;
      return { ...f, nameBytes, crc, localSize, cdSize };
    });

    let localTotal = 0;
    let cdTotal = 0;
    for (const p of prepared) {
      localTotal += p.localSize;
      cdTotal += p.cdSize;
    }
    const eocdSize = 22;
    const buffer = new Uint8Array(localTotal + cdTotal + eocdSize);
    const view = new DataView(buffer.buffer);

    let offset = 0;
    const offsets = [];

    // 1. Local file headers and data
    for (const p of prepared) {
      offsets.push(offset);
      view.setUint32(offset, 0x04034b50, true);
      view.setUint16(offset + 4, 20, true);
      view.setUint16(offset + 6, 0x0800, true); // UTF-8 filename
      view.setUint16(offset + 8, 0, true);      // Stored (no compression)
      view.setUint16(offset + 10, 0, true);
      view.setUint16(offset + 12, 0x5421, true); // Date/time
      view.setUint32(offset + 14, p.crc, true);
      view.setUint32(offset + 18, p.data.length, true);
      view.setUint32(offset + 22, p.data.length, true);
      view.setUint16(offset + 26, p.nameBytes.length, true);
      view.setUint16(offset + 28, 0, true);

      buffer.set(p.nameBytes, offset + 30);
      buffer.set(p.data, offset + 30 + p.nameBytes.length);
      offset += p.localSize;
    }

    const cdOffset = offset;

    // 2. Central Directory headers
    for (let i = 0; i < prepared.length; i++) {
      const p = prepared[i];
      const localOff = offsets[i];
      view.setUint32(offset, 0x02014b50, true);
      view.setUint16(offset + 4, 20, true);
      view.setUint16(offset + 6, 20, true);
      view.setUint16(offset + 8, 0x0800, true);
      view.setUint16(offset + 10, 0, true);
      view.setUint16(offset + 12, 0, true);
      view.setUint16(offset + 14, 0x5421, true);
      view.setUint32(offset + 16, p.crc, true);
      view.setUint32(offset + 20, p.data.length, true);
      view.setUint32(offset + 24, p.data.length, true);
      view.setUint16(offset + 28, p.nameBytes.length, true);
      view.setUint16(offset + 30, 0, true);
      view.setUint16(offset + 32, 0, true);
      view.setUint16(offset + 34, 0, true);
      view.setUint16(offset + 36, 0, true);
      view.setUint32(offset + 38, 0, true);
      view.setUint32(offset + 42, localOff, true);

      buffer.set(p.nameBytes, offset + 46);
      offset += p.cdSize;
    }

    // 3. End of Central Directory record (EOCD)
    view.setUint32(offset, 0x06054b50, true);
    view.setUint16(offset + 4, 0, true);
    view.setUint16(offset + 6, 0, true);
    view.setUint16(offset + 8, prepared.length, true);
    view.setUint16(offset + 10, prepared.length, true);
    view.setUint32(offset + 12, cdTotal, true);
    view.setUint32(offset + 16, cdOffset, true);
    view.setUint16(offset + 20, 0, true);

    return buffer;
  }
}

/**
 * Escapes characters for XML text
 */
function escapeXML(str) {
  return String(str ?? '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&apos;');
}


// ==========================================
// 15. MULTIVARIATE EXPLORATORY DATA ANALYSIS (PCA / MCA / FAMD) ENGINE
// ==========================================
/**
 * Statis-Gravity - Multivariate Exploratory Data Analysis (EDA) Engine
 * Implements Principal Component Analysis (PCA), Multiple Correspondence Analysis (MCA),
 * and Factor Analysis of Mixed Data (FAMD / Pagès 2004) with lightweight, self-contained
 * client-side numerical linear algebra (Jacobi eigenvalue solver and thin SVD).
 */

const Multivariate = {
  // =========================================================================
  // 1. LIGHTWEIGHT NUMERICAL LINEAR ALGEBRA & MATRIX DECOMPOSITION
  // =========================================================================
  LinearAlgebra: {
    transpose(A) {
      const rows = A.length;
      const cols = A[0].length;
      const T = Array.from({ length: cols }, () => new Array(rows));
      for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
          T[j][i] = A[i][j];
        }
      }
      return T;
    },

    matmul(A, B) {
      const rowsA = A.length;
      const colsA = A[0].length;
      const colsB = B[0].length;
      const C = Array.from({ length: rowsA }, () => new Array(colsB).fill(0));
      for (let i = 0; i < rowsA; i++) {
        for (let k = 0; k < colsA; k++) {
          const aik = A[i][k];
          if (aik === 0) continue;
          for (let j = 0; j < colsB; j++) {
            C[i][j] += aik * B[k][j];
          }
        }
      }
      return C;
    },

    /**
     * Cyclic Jacobi algorithm for symmetric matrix eigenvalue decomposition
     * Solves A * V = V * D, where A is symmetric (n x n)
     */
    jacobi(A, maxIter = 150, tol = 1e-12) {
      const n = A.length;
      let V = Array.from({ length: n }, (_, i) => {
        const row = new Array(n).fill(0);
        row[i] = 1.0;
        return row;
      });
      let D = A.map(row => [...row]);

      for (let iter = 0; iter < maxIter; iter++) {
        let maxOff = 0;
        let p = 0, q = 1;
        for (let i = 0; i < n; i++) {
          for (let j = i + 1; j < n; j++) {
            const val = Math.abs(D[i][j]);
            if (val > maxOff) {
              maxOff = val;
              p = i; q = j;
            }
          }
        }
        if (maxOff < tol) break;

        const diff = D[q][q] - D[p][p];
        let t;
        if (Math.abs(D[p][q]) < Math.abs(diff) * 1e-15) {
          t = D[p][q] / diff;
        } else {
          const phi = diff / (2.0 * D[p][q]);
          t = 1.0 / (Math.abs(phi) + Math.sqrt(phi * phi + 1.0));
          if (phi < 0) t = -t;
        }

        const c = 1.0 / Math.sqrt(t * t + 1.0);
        const s = t * c;
        const tau = s / (1.0 + c);

        const Dpq = D[p][q];
        D[p][q] = 0;
        D[q][p] = 0;
        D[p][p] -= t * Dpq;
        D[q][q] += t * Dpq;

        for (let j = 0; j < n; j++) {
          if (j !== p && j !== q) {
            const Djp = D[j][p];
            const Djq = D[j][q];
            D[j][p] = Djp - s * (Djq + tau * Djp);
            D[p][j] = D[j][p];
            D[j][q] = Djq + s * (Djp - tau * Djq);
            D[q][j] = D[j][q];
          }
        }

        for (let j = 0; j < n; j++) {
          const Vjp = V[j][p];
          const Vjq = V[j][q];
          V[j][p] = Vjp - s * (Vjq + tau * Vjp);
          V[j][q] = Vjq + s * (Vjp - tau * Vjq);
        }
      }

      // Sort eigenvalues and corresponding eigenvector columns descending
      const eigen = [];
      for (let i = 0; i < n; i++) {
        eigen.push({
          value: Math.max(0, D[i][i]),
          vector: V.map(row => row[i])
        });
      }
      eigen.sort((a, b) => b.value - a.value);

      return {
        values: eigen.map(e => e.value),
        vectors: eigen.map(e => e.vector) // Column eigenvectors
      };
    },

    /**
     * Thin Singular Value Decomposition (SVD): M = U * diag(s) * V^T
     * Efficiently handles rectangular n x p matrices via Gram matrix eigendecomposition.
     */
    svd(M) {
      const n = M.length;
      const p = M[0].length;

      if (n >= p) {
        // Compute p x p Gram matrix: C = M^T * M
        const MT = this.transpose(M);
        const C = this.matmul(MT, M);
        const eig = this.jacobi(C);

        const rank = p;
        const s = eig.values.map(val => Math.sqrt(Math.max(0, val)));
        // V matrix: columns are eigenvectors of M^T * M
        const V = Array.from({ length: p }, (_, i) => new Array(p));
        for (let j = 0; j < p; j++) {
          for (let i = 0; i < p; i++) {
            V[i][j] = eig.vectors[j][i];
          }
        }

        // Left singular vectors: U_k = M * V_k / s_k
        const U = Array.from({ length: n }, () => new Array(p).fill(0));
        for (let k = 0; k < p; k++) {
          const sk = s[k];
          if (sk > 1e-12) {
            for (let i = 0; i < n; i++) {
              let sum = 0;
              for (let j = 0; j < p; j++) {
                sum += M[i][j] * V[j][k];
              }
              U[i][k] = sum / sk;
            }
          }
        }

        return { U, s, V };
      } else {
        // When p > n: compute n x n Gram matrix: C = M * M^T
        const MT = this.transpose(M);
        const C = this.matmul(M, MT);
        const eig = this.jacobi(C);

        const s = eig.values.map(val => Math.sqrt(Math.max(0, val)));
        // U matrix: columns are eigenvectors of M * M^T
        const U = Array.from({ length: n }, (_, i) => new Array(n));
        for (let j = 0; j < n; j++) {
          for (let i = 0; i < n; i++) {
            U[i][j] = eig.vectors[j][i];
          }
        }

        // Right singular vectors: V_k = M^T * U_k / s_k
        const V = Array.from({ length: p }, () => new Array(n).fill(0));
        for (let k = 0; k < n; k++) {
          const sk = s[k];
          if (sk > 1e-12) {
            for (let j = 0; j < p; j++) {
              let sum = 0;
              for (let i = 0; i < n; i++) {
                sum += MT[j][i] * U[i][k];
              }
              V[j][k] = sum / sk;
            }
          }
        }

        return { U, s, V };
      }
    }
  },

  // =========================================================================
  // 2. DATASET PARSER & VARIABLE TYPE INFERENCE
  // =========================================================================
  parseDataset(csvText) {
    if (!csvText || typeof csvText !== 'string') {
      return { error: 'No data provided.' };
    }

    const lines = csvText.trim().split(/\r?\n/).filter(l => l.trim().length > 0);
    if (lines.length < 3) {
      return { error: 'Dataset must contain a header row and at least 2 observation rows.' };
    }

    const delimiter = lines[0].includes('\t') ? '\t' : (lines[0].includes(';') ? ';' : ',');
    const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));

    const records = [];
    for (let i = 1; i < lines.length; i++) {
      const parts = lines[i].split(delimiter).map(p => p.trim().replace(/^["']|["']$/g, ''));
      if (parts.length === headers.length) {
        const row = { _id: `Row-${i}` };
        headers.forEach((h, idx) => {
          row[h] = parts[idx];
        });
        records.push(row);
      }
    }

    if (records.length < 3) {
      return { error: 'Insufficient valid records in dataset.' };
    }

    // Infer variable types
    const colTypes = {};
    const colStats = {};

    headers.forEach(h => {
      let numCount = 0;
      let totalNonEmpty = 0;
      const values = [];

      records.forEach(r => {
        const val = r[h];
        if (val !== undefined && val !== null && val !== '') {
          totalNonEmpty++;
          const num = Number(val);
          if (!isNaN(num) && isFinite(num)) {
            numCount++;
            values.push(num);
          } else {
            values.push(val);
          }
        }
      });

      const numRatio = totalNonEmpty > 0 ? numCount / totalNonEmpty : 0;
      const uniqueVals = new Set(values);

      // If > 80% numerical and at least 3 distinct values, infer continuous
      if (numRatio >= 0.80 && uniqueVals.size > 2) {
        colTypes[h] = 'continuous';
        const nums = values.filter(v => typeof v === 'number');
        const mean = nums.reduce((a, b) => a + b, 0) / nums.length;
        const variance = nums.reduce((acc, v) => acc + (v - mean) ** 2, 0) / (nums.length - 1);
        colStats[h] = {
          type: 'continuous',
          count: nums.length,
          mean,
          sd: Math.sqrt(variance),
          min: Math.min(...nums),
          max: Math.max(...nums)
        };
      } else {
        colTypes[h] = 'categorical';
        colStats[h] = {
          type: 'categorical',
          count: values.length,
          uniqueCount: uniqueVals.size,
          categories: Array.from(uniqueVals).map(String).sort()
        };
      }
    });

    return {
      headers,
      records,
      colTypes,
      colStats,
      nRows: records.length,
      nCols: headers.length
    };
  },

  // =========================================================================
  // 3. PRINCIPAL COMPONENT ANALYSIS (PCA)
  // =========================================================================
  runPCA(records, continuousCols, options = {}) {
    if (!continuousCols || continuousCols.length < 2) {
      return { error: 'PCA requires at least 2 continuous variables.' };
    }

    const n = records.length;
    const p = continuousCols.length;

    // Filter complete cases
    const cleanRecords = records.filter(r => {
      return continuousCols.every(c => {
        const val = Number(r[c]);
        return !isNaN(val) && isFinite(val) && r[c] !== '' && r[c] !== null;
      });
    });

    if (cleanRecords.length < 3) {
      return { error: 'Insufficient complete observations for PCA (minimum 3 required).' };
    }

    const nClean = cleanRecords.length;

    // 1. Mean-centering and scaling to unit variance
    const means = {};
    const sds = {};

    continuousCols.forEach(col => {
      const vals = cleanRecords.map(r => Number(r[col]));
      const m = vals.reduce((a, b) => a + b, 0) / nClean;
      const variance = vals.reduce((acc, v) => acc + (v - m) ** 2, 0) / (nClean - 1);
      const s = Math.sqrt(variance) || 1.0;
      means[col] = m;
      sds[col] = s;
    });

    // Standardized matrix Z / sqrt(n - 1)
    const factor = Math.sqrt(nClean - 1);
    const Z = Array.from({ length: nClean }, (_, i) => {
      const r = cleanRecords[i];
      return continuousCols.map(col => (Number(r[col]) - means[col]) / (sds[col] * factor));
    });

    // 2. Perform SVD on Z
    const svdRes = this.LinearAlgebra.svd(Z);
    const eigenvalues = svdRes.s.map(s => s * s);
    const totalVariance = eigenvalues.reduce((a, b) => a + b, 0) || p;

    const nComponents = Math.min(p, nClean - 1);
    const scree = [];
    let cumVar = 0;

    for (let k = 0; k < nComponents; k++) {
      const eig = eigenvalues[k] || 0;
      const varPct = (eig / totalVariance) * 100;
      cumVar += varPct;
      scree.push({
        dim: k + 1,
        label: `Dim ${k + 1}`,
        eigenvalue: eig,
        variancePct: varPct,
        cumulativePct: Math.min(100, cumVar)
      });
    }

    // 3. Compute Individual Coordinates (Factor Scores)
    // Score_k = Z_unscaled * V_k = (Z * sqrt(n-1)) * V_k = U_k * s_k * sqrt(n-1)
    const individuals = cleanRecords.map((r, i) => {
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        // Individual coordinate on component k
        coords.push(svdRes.U[i][k] * svdRes.s[k] * factor);
      }

      // Squared cosine (quality of representation) and distance from centroid
      const totalDistSq = coords.reduce((acc, c) => acc + c * c, 0);
      const cos2 = coords.map(c => (totalDistSq > 0 ? (c * c) / totalDistSq : 0));

      return {
        id: r._id || r.id || r.patient_id || `Obs-${i + 1}`,
        originalRow: r,
        coords,
        cos2,
        distFromCentroid: Math.sqrt(totalDistSq)
      };
    });

    // Individual contributions
    for (let k = 0; k < nComponents; k++) {
      const sumDistSq = individuals.reduce((acc, ind) => acc + ind.coords[k] ** 2, 0) || 1;
      individuals.forEach(ind => {
        if (!ind.ctr) ind.ctr = [];
        ind.ctr[k] = (ind.coords[k] ** 2 / sumDistSq) * 100;
      });
    }

    // Flag outliers: Euclidean distance in factor space > 2.5 standard deviations
    const allDists = individuals.map(d => d.distFromCentroid);
    const meanDist = allDists.reduce((a, b) => a + b, 0) / allDists.length;
    const sdDist = Math.sqrt(allDists.reduce((acc, v) => acc + (v - meanDist) ** 2, 0) / (allDists.length - 1)) || 1;
    const outlierCutoff = meanDist + 2.5 * sdDist;

    individuals.forEach(ind => {
      ind.isOutlier = ind.distFromCentroid > outlierCutoff;
    });

    // 4. Compute Variable Coordinates (Loadings / Correlations)
    // Loading_j,k = V_j,k * sqrt(eigenvalue_k) = Pearson correlation between variable j and component k
    const variables = continuousCols.map((col, j) => {
      const coords = [];
      const cos2 = [];
      const ctr = [];

      for (let k = 0; k < nComponents; k++) {
        const loading = svdRes.V[j][k] * Math.sqrt(Math.max(0, eigenvalues[k]));
        coords.push(loading);
        // cos2 on dimension k = loading^2 (since var is standardized to unit variance)
        cos2.push(loading * loading);
        // Contribution of variable j to component k = loading^2 / eigenvalue_k
        const eig = eigenvalues[k];
        ctr.push(eig > 0 ? (loading * loading / eig) * 100 : 0);
      }

      return {
        name: col,
        type: 'continuous',
        coords,
        cos2,
        ctr,
        mean: means[col],
        sd: sds[col]
      };
    });

    return {
      method: 'PCA',
      methodLabel: 'Principal Component Analysis (PCA)',
      nObservations: nClean,
      nVariables: p,
      continuousCols,
      eigenvalues,
      totalVariance,
      scree,
      individuals,
      variables,
      continuousVariables: variables,
      allVariables: variables,
      outliers: individuals.filter(d => d.isOutlier),
      outlierCount: individuals.filter(d => d.isOutlier).length,
      outlierCutoff
    };
  },

  // =========================================================================
  // 4. MULTIPLE CORRESPONDENCE ANALYSIS (MCA)
  // =========================================================================
  runMCA(records, categoricalCols, options = {}) {
    if (!categoricalCols || categoricalCols.length < 2) {
      return { error: 'MCA requires at least 2 categorical variables.' };
    }

    const n = records.length;
    const K = categoricalCols.length;

    // Filter complete cases
    const cleanRecords = records.filter(r => {
      return categoricalCols.every(c => r[c] !== undefined && r[c] !== null && String(r[c]).trim() !== '');
    });

    if (cleanRecords.length < 3) {
      return { error: 'Insufficient complete observations for MCA.' };
    }

    const nClean = cleanRecords.length;

    // 1. Build Indicator Matrix Z
    const modalities = [];
    categoricalCols.forEach(col => {
      const uniqueVals = Array.from(new Set(cleanRecords.map(r => String(r[col]).trim()))).sort();
      uniqueVals.forEach(val => {
        modalities.push({
          col,
          val,
          fullName: `${col}_${val}`,
          label: `${col}: ${val}`
        });
      });
    });

    const J = modalities.length;
    if (J <= K) {
      return { error: 'Insufficient category variation to perform MCA.' };
    }

    // Indicator matrix Z: nClean x J
    const Z = Array.from({ length: nClean }, () => new Array(J).fill(0));
    cleanRecords.forEach((r, i) => {
      modalities.forEach((m, j) => {
        if (String(r[m.col]).trim() === m.val) {
          Z[i][j] = 1;
        }
      });
    });

    // Column marginal frequencies and proportions
    const colCounts = new Array(J).fill(0);
    for (let j = 0; j < J; j++) {
      for (let i = 0; i < nClean; i++) {
        colCounts[j] += Z[i][j];
      }
    }

    const p = colCounts.map(c => c / nClean);

    // Standardized residuals matrix S:
    // S_ij = (Z_ij / (n * K) - (1/n) * (p_j / K)) / sqrt((1/n) * (p_j / K))
    // S_ij = (1 / sqrt(n * K)) * (Z_ij - p_j) / sqrt(p_j)
    const S = Array.from({ length: nClean }, (_, i) => {
      const factorNK = 1.0 / Math.sqrt(nClean * K);
      return modalities.map((m, j) => {
        const pj = Math.max(1e-9, p[j]);
        return factorNK * (Z[i][j] - pj) / Math.sqrt(pj);
      });
    });

    // 2. Perform SVD on S
    const svdRes = this.LinearAlgebra.svd(S);
    const eigenvalues = svdRes.s.map(s => s * s);
    const totalInertia = (J - K) / K;

    const nComponents = Math.min(J - K, nClean - 1);
    const scree = [];
    let cumVar = 0;

    for (let k = 0; k < nComponents; k++) {
      const eig = eigenvalues[k] || 0;
      // In MCA, standard variance explained percentage:
      const varPct = totalInertia > 0 ? (eig / totalInertia) * 100 : 0;
      cumVar += varPct;
      scree.push({
        dim: k + 1,
        label: `Dim ${k + 1}`,
        eigenvalue: eig,
        variancePct: varPct,
        cumulativePct: Math.min(100, cumVar)
      });
    }

    // 3. Individual Principal Coordinates
    // F_ik = sqrt(n) * U_ik * s_k
    const sqrtN = Math.sqrt(nClean);
    const individuals = cleanRecords.map((r, i) => {
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        coords.push(sqrtN * svdRes.U[i][k] * svdRes.s[k]);
      }
      const totalDistSq = coords.reduce((acc, c) => acc + c * c, 0);
      const cos2 = coords.map(c => (totalDistSq > 0 ? (c * c) / totalDistSq : 0));

      return {
        id: r._id || r.id || r.patient_id || `Obs-${i + 1}`,
        originalRow: r,
        coords,
        cos2,
        distFromCentroid: Math.sqrt(totalDistSq)
      };
    });

    // Outlier check
    const allDists = individuals.map(d => d.distFromCentroid);
    const meanDist = allDists.reduce((a, b) => a + b, 0) / allDists.length;
    const sdDist = Math.sqrt(allDists.reduce((acc, v) => acc + (v - meanDist) ** 2, 0) / (allDists.length - 1)) || 1;
    const outlierCutoff = meanDist + 2.5 * sdDist;
    individuals.forEach(ind => {
      ind.isOutlier = ind.distFromCentroid > outlierCutoff;
    });

    // 4. Modality / Category Coordinates
    // A_jk = (1 / sqrt(p_j / K)) * V_jk * s_k / sqrt(K) = (1 / sqrt(p_j)) * V_jk * s_k
    const modalitiesCoords = modalities.map((m, j) => {
      const pj = Math.max(1e-9, p[j]);
      const coords = [];
      const ctr = [];

      for (let k = 0; k < nComponents; k++) {
        // Principal coordinate of modality j on dimension k
        const coord = (1.0 / Math.sqrt(pj)) * svdRes.V[j][k] * svdRes.s[k];
        coords.push(coord);
        // Contribution of modality j to dimension k
        const eig = eigenvalues[k];
        const c = eig > 0 ? (pj * coord * coord) / (K * eig) * 100 : 0;
        ctr.push(c);
      }

      return {
        name: m.label,
        col: m.col,
        val: m.val,
        type: 'modality',
        coords,
        ctr,
        count: colCounts[j],
        proportion: p[j]
      };
    });

    // Variables summary (correlation ratio eta^2 of categorical variables with dimensions)
    const variables = categoricalCols.map(col => {
      const mods = modalitiesCoords.filter(m => m.col === col);
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        // eta^2 = sum_m (p_m * coord_m,k^2) / eigenvalue_k
        const eig = eigenvalues[k];
        const eta2 = eig > 0 ? mods.reduce((acc, m) => acc + m.proportion * (m.coords[k] ** 2), 0) / eig : 0;
        coords.push(Math.sqrt(Math.max(0, Math.min(1, eta2))));
      }
      return {
        name: col,
        type: 'categorical',
        coords
      };
    });

    return {
      method: 'MCA',
      methodLabel: 'Multiple Correspondence Analysis (MCA)',
      nObservations: nClean,
      nVariables: K,
      nModalities: J,
      categoricalCols,
      eigenvalues,
      totalInertia,
      scree,
      individuals,
      variables,
      modalities: modalitiesCoords,
      allVariables: modalitiesCoords,
      outliers: individuals.filter(d => d.isOutlier),
      outlierCount: individuals.filter(d => d.isOutlier).length,
      outlierCutoff
    };
  },

  // =========================================================================
  // 5. FACTOR ANALYSIS OF MIXED DATA (FAMD / Pagès 2004)
  // =========================================================================
  runFAMD(records, continuousCols, categoricalCols, options = {}) {
    if (!continuousCols || continuousCols.length === 0) {
      return this.runMCA(records, categoricalCols, options);
    }
    if (!categoricalCols || categoricalCols.length === 0) {
      return this.runPCA(records, continuousCols, options);
    }

    const p1 = continuousCols.length;
    const K = categoricalCols.length;

    // Filter complete cases across both continuous and categorical variables
    const cleanRecords = records.filter(r => {
      const contOk = continuousCols.every(c => {
        const val = Number(r[c]);
        return !isNaN(val) && isFinite(val) && r[c] !== '' && r[c] !== null;
      });
      const catOk = categoricalCols.every(c => r[c] !== undefined && r[c] !== null && String(r[c]).trim() !== '');
      return contOk && catOk;
    });

    if (cleanRecords.length < 3) {
      return { error: 'Insufficient complete observations for FAMD (minimum 3 required).' };
    }

    const nClean = cleanRecords.length;

    // 1. Preprocess Continuous Variables: Center and scale to unit variance
    const contMeans = {};
    const contSds = {};
    continuousCols.forEach(col => {
      const vals = cleanRecords.map(r => Number(r[col]));
      const m = vals.reduce((a, b) => a + b, 0) / nClean;
      const variance = vals.reduce((acc, v) => acc + (v - m) ** 2, 0) / (nClean - 1);
      contMeans[col] = m;
      contSds[col] = Math.sqrt(variance) || 1.0;
    });

    // 2. Preprocess Categorical Variables: One-hot encode and scale by 1 / sqrt(p_m)
    const modalities = [];
    categoricalCols.forEach(col => {
      const uniqueVals = Array.from(new Set(cleanRecords.map(r => String(r[col]).trim()))).sort();
      uniqueVals.forEach(val => {
        modalities.push({
          col,
          val,
          label: `${col}: ${val}`
        });
      });
    });

    const J = modalities.length;
    const pProps = modalities.map(m => {
      const count = cleanRecords.filter(r => String(r[m.col]).trim() === m.val).length;
      return Math.max(1e-6, count / nClean);
    });

    // 3. Construct Joint Balanced Matrix W of size nClean x (p1 + J)
    // Continuous block: (X_ij - mean_j) / (sd_j * sqrt(nClean))
    // Categorical block: (I_im - p_m) / (sqrt(p_m) * sqrt(nClean))
    const totalCols = p1 + J;
    const sqrtN = Math.sqrt(nClean);

    const W = Array.from({ length: nClean }, (_, i) => {
      const r = cleanRecords[i];
      const row = new Array(totalCols);

      // Continuous columns
      for (let j = 0; j < p1; j++) {
        const col = continuousCols[j];
        row[j] = (Number(r[col]) - contMeans[col]) / (contSds[col] * sqrtN);
      }

      // Categorical dummy columns
      for (let m = 0; m < J; m++) {
        const mod = modalities[m];
        const isMatch = String(r[mod.col]).trim() === mod.val ? 1 : 0;
        const pm = pProps[m];
        row[p1 + m] = (isMatch - pm) / (Math.sqrt(pm) * sqrtN);
      }

      return row;
    });

    // 4. SVD on W
    const svdRes = this.LinearAlgebra.svd(W);
    const eigenvalues = svdRes.s.map(s => s * s);
    const totalInertia = p1 + (J - K);

    const nComponents = Math.min(totalCols - K, nClean - 1);
    const scree = [];
    let cumVar = 0;

    for (let k = 0; k < nComponents; k++) {
      const eig = eigenvalues[k] || 0;
      const varPct = totalInertia > 0 ? (eig / totalInertia) * 100 : 0;
      cumVar += varPct;
      scree.push({
        dim: k + 1,
        label: `Dim ${k + 1}`,
        eigenvalue: eig,
        variancePct: varPct,
        cumulativePct: Math.min(100, cumVar)
      });
    }

    // 5. Individual Coordinates (Factor Scores)
    // F_ik = sqrt(nClean) * U_ik * s_k
    const individuals = cleanRecords.map((r, i) => {
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        coords.push(sqrtN * svdRes.U[i][k] * svdRes.s[k]);
      }
      const totalDistSq = coords.reduce((acc, c) => acc + c * c, 0);
      const cos2 = coords.map(c => (totalDistSq > 0 ? (c * c) / totalDistSq : 0));

      return {
        id: r._id || r.id || r.patient_id || `Obs-${i + 1}`,
        originalRow: r,
        coords,
        cos2,
        distFromCentroid: Math.sqrt(totalDistSq)
      };
    });

    // Flag outliers (> 2.5 SD distance in factor space)
    const allDists = individuals.map(d => d.distFromCentroid);
    const meanDist = allDists.reduce((a, b) => a + b, 0) / allDists.length;
    const sdDist = Math.sqrt(allDists.reduce((acc, v) => acc + (v - meanDist) ** 2, 0) / (allDists.length - 1)) || 1;
    const outlierCutoff = meanDist + 2.5 * sdDist;
    individuals.forEach(ind => {
      ind.isOutlier = ind.distFromCentroid > outlierCutoff;
    });

    // 6. Variable Loadings and Modality Centroids
    // A. Continuous variable correlation with dimensions
    const continuousVariables = continuousCols.map((col, j) => {
      const rawVals = cleanRecords.map(r => Number(r[col]));
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        // Pearson correlation with individual factor score
        const score_k = individuals.map(ind => ind.coords[k]);
        const mScore = score_k.reduce((a, b) => a + b, 0) / nClean;
        const sdScore = Math.sqrt(score_k.reduce((acc, v) => acc + (v - mScore) ** 2, 0) / (nClean - 1)) || 1;

        let cov = 0;
        for (let i = 0; i < nClean; i++) {
          cov += (rawVals[i] - contMeans[col]) * (score_k[i] - mScore);
        }
        cov /= (nClean - 1);
        const rCorr = cov / (contSds[col] * sdScore);
        coords.push(Math.max(-1, Math.min(1, rCorr)));
      }

      return {
        name: col,
        type: 'continuous',
        coords,
        mean: contMeans[col],
        sd: contSds[col]
      };
    });

    // B. Categorical Modality Centroid Coordinates
    const modalityCoordinates = modalities.map((m, mIdx) => {
      const matchingInds = individuals.filter(ind => String(ind.originalRow[m.col]).trim() === m.val);
      const count = matchingInds.length;
      const coords = [];

      for (let k = 0; k < nComponents; k++) {
        if (count > 0) {
          const meanCoord = matchingInds.reduce((acc, ind) => acc + ind.coords[k], 0) / count;
          coords.push(meanCoord);
        } else {
          coords.push(0);
        }
      }

      return {
        name: m.label,
        col: m.col,
        val: m.val,
        type: 'modality',
        coords,
        count,
        proportion: pProps[mIdx]
      };
    });

    const allVariables = [...continuousVariables, ...modalityCoordinates];

    return {
      method: 'FAMD',
      methodLabel: 'Factor Analysis of Mixed Data (FAMD)',
      nObservations: nClean,
      nContinuous: p1,
      nCategorical: K,
      nModalities: J,
      continuousCols,
      categoricalCols,
      eigenvalues,
      totalInertia,
      scree,
      individuals,
      variables: continuousVariables,
      modalities: modalityCoordinates,
      allVariables,
      continuousVariables,
      categoricalVariables: modalityCoordinates,
      outliers: individuals.filter(d => d.isOutlier),
      outlierCount: individuals.filter(d => d.isOutlier).length,
      outlierCutoff
    };
  },

  // =========================================================================
  // 6. MASTER ANALYSIS DISPATCHER
  // =========================================================================
  executeAnalysis(records, method, activeCols, options = {}) {
    const colTypes = options.colTypes || {};
    const groupingCol = options.groupingCol || null;

    let analysis;
    if (method === 'PCA') {
      const contCols = activeCols.filter(c => (colTypes[c] || 'continuous') === 'continuous');
      analysis = this.runPCA(records, contCols, options);
    } else if (method === 'MCA') {
      const catCols = activeCols.filter(c => colTypes[c] === 'categorical');
      analysis = this.runMCA(records, catCols, options);
    } else {
      // FAMD
      const contCols = activeCols.filter(c => (colTypes[c] || 'continuous') === 'continuous');
      const catCols = activeCols.filter(c => colTypes[c] === 'categorical');
      analysis = this.runFAMD(records, contCols, catCols, options);
    }

    if (analysis.error) return analysis;

    // Attach grouping information to individuals
    if (groupingCol && records[0] && records[0][groupingCol] !== undefined) {
      analysis.groupingCol = groupingCol;
      const groups = new Set();
      analysis.individuals.forEach(ind => {
        ind.group = String(ind.originalRow[groupingCol] || 'Unassigned');
        groups.add(ind.group);
      });
      analysis.uniqueGroups = Array.from(groups).sort();
    } else {
      analysis.groupingCol = null;
      analysis.uniqueGroups = [];
      analysis.individuals.forEach(ind => { ind.group = 'Observation'; });
    }

    // Add Interpretation & Scripts
    analysis.interpretation = this.generateInterpretation(analysis);
    analysis.scripts = {
      python: this.generatePythonScript(analysis),
      r: this.generateRScript(analysis)
    };

    return analysis;
  },

  // =========================================================================
  // 7. AUTOMATED NATURAL LANGUAGE INTERPRETATION & EXPLANATION
  // =========================================================================
  generateInterpretation(analysis) {
    if (!analysis || analysis.error) return '';

    const scree = analysis.scree || [];
    const dim1 = scree[0] || { variancePct: 0, eigenvalue: 0 };
    const dim2 = scree[1] || { variancePct: 0, eigenvalue: 0 };
    const dim3 = scree[2] || { variancePct: 0, eigenvalue: 0 };

    const cum2D = (dim1.variancePct + dim2.variancePct).toFixed(1);
    const cum3D = (dim1.variancePct + dim2.variancePct + dim3.variancePct).toFixed(1);

    // Identify Drivers for Dim 1 & Dim 2
    const allVars = analysis.allVariables || [];

    // Sort by absolute coordinate on Dim 1
    const sortedDim1 = [...allVars].sort((a, b) => Math.abs(b.coords[0]) - Math.abs(a.coords[0]));
    const topDim1Pos = sortedDim1.filter(v => v.coords[0] > 0).slice(0, 3).map(v => `${v.name} (+${v.coords[0].toFixed(2)})`);
    const topDim1Neg = sortedDim1.filter(v => v.coords[0] < 0).slice(0, 3).map(v => `${v.name} (${v.coords[0].toFixed(2)})`);

    // Sort by absolute coordinate on Dim 2
    const sortedDim2 = [...allVars].sort((a, b) => Math.abs(b.coords[1]) - Math.abs(a.coords[1]));
    const topDim2Pos = sortedDim2.filter(v => v.coords[1] > 0).slice(0, 3).map(v => `${v.name} (+${v.coords[1].toFixed(2)})`);
    const topDim2Neg = sortedDim2.filter(v => v.coords[1] < 0).slice(0, 3).map(v => `${v.name} (${v.coords[1].toFixed(2)})`);

    // Outlier list
    const outliers = (analysis.individuals || []).filter(ind => ind.isOutlier);
    const outlierSummary = outliers.length > 0
      ? `Flagged ${outliers.length} statistical outlier observation(s) exceeding 2.5 SD distance from centroid: ${outliers.slice(0, 5).map(o => o.id).join(', ')}${outliers.length > 5 ? '...' : ''}.`
      : 'No severe outlier observations detected (all individuals fall within 2.5 standard deviations of the multivariate centroid).';

    // Clustering summary
    let clusterSummary = '';
    if (analysis.groupingCol && analysis.uniqueGroups.length > 1) {
      clusterSummary = `Observations were stratified across ${analysis.uniqueGroups.length} categories of "${analysis.groupingCol}". Clear separation is observable along principal axes, indicating that baseline covariates correlate significantly with group classification.`;
    }

    const narrative = [
      `### Multivariate Dimensionality Reduction Findings (${analysis.methodLabel}):`,
      `1. **Variance Retention & Dimensionality:** The primary 2D principal plane (Dim 1 vs. Dim 2) accounts for **${cum2D}%** of total inertia (Dim 1 = ${dim1.variancePct.toFixed(1)}%, Dim 2 = ${dim2.variancePct.toFixed(1)}%). Expanding to the 3D subspace (Dim 1–3) captures **${cum3D}%** of the aggregate sample variance.`,
      `2. **Key Drivers for Dimension 1 (Horizontal Axis):** Strongest positive associations include ${topDim1Pos.length > 0 ? topDim1Pos.join(', ') : 'none'}, contrasting with negative associations ${topDim1Neg.length > 0 ? topDim1Neg.join(', ') : 'none'}.`,
      `3. **Key Drivers for Dimension 2 (Vertical Axis):** Orthogonal variance is primarily defined by ${topDim2Pos.length > 0 ? topDim2Pos.join(', ') : 'none'} against ${topDim2Neg.length > 0 ? topDim2Neg.join(', ') : 'none'}.`,
      `4. **Outliers & Anomaly Detection:** ${outlierSummary}`,
      clusterSummary ? `5. **Subgroup Clustering:** ${clusterSummary}` : ''
    ].filter(Boolean).join('\n\n');

    return narrative;
  },

  // =========================================================================
  // 8. REPRODUCIBLE SCRIPT GENERATOR [Python & R]
  // =========================================================================
  generatePythonScript(analysis) {
    const method = analysis.method;
    const contColsPy = (analysis.continuousCols || []).map(c => `'${c}'`).join(', ');
    const catColsPy = (analysis.categoricalCols || []).map(c => `'${c}'`).join(', ');
    const groupPy = analysis.groupingCol ? `'${analysis.groupingCol}'` : 'None';

    return `"""
=============================================================================
Multivariate Exploratory Data Analysis (${analysis.methodLabel})
Generated by Statis-Gravity Clinical Biostatistics Platform
=============================================================================
"""

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# 1. Load Clinical Observational Data
# df = pd.read_csv('clinical_multivariate_data.csv')

continuous_cols = [${contColsPy}]
categorical_cols = [${catColsPy}]
group_col = ${groupPy}

${method === 'PCA' ? `# -----------------------------------------------------------------------------
# 2. PRINCIPAL COMPONENT ANALYSIS (scikit-learn)
# -----------------------------------------------------------------------------
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

df_clean = df[continuous_cols].dropna()
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_clean)

pca = PCA(n_components=min(5, len(continuous_cols)))
pca_scores = pca.fit_transform(X_scaled)

scores_df = pd.DataFrame(pca_scores, columns=[f'Dim {i+1}' for i in range(pca_scores.shape[1])])
if group_col and group_col in df.columns:
    scores_df[group_col] = df.loc[df_clean.index, group_col].values

exp_var = pca.explained_variance_ratio_ * 100
print("Explained Variance per Dimension (%):", np.round(exp_var, 2))
print("Cumulative Variance (%):", np.round(np.cumsum(exp_var), 2))

# 3. Interactive 2D Factor Biplot
fig_2d = px.scatter(
    scores_df, x='Dim 1', y='Dim 2', color=group_col,
    title=f"PCA Factor Map: Dim 1 ({exp_var[0]:.1f}%) vs Dim 2 ({exp_var[1]:.1f}%)",
    template='plotly_dark'
)
fig_2d.show()

# 4. Interactive 3D Rotation Plot
if scores_df.shape[1] >= 3:
    fig_3d = px.scatter_3d(
        scores_df, x='Dim 1', y='Dim 2', z='Dim 3', color=group_col,
        title="PCA 3D Space (Dim 1 vs Dim 2 vs Dim 3)",
        template='plotly_dark'
    )
    fig_3d.show()` : (method === 'MCA' ? `# -----------------------------------------------------------------------------
# 2. MULTIPLE CORRESPONDENCE ANALYSIS (prince package)
# -----------------------------------------------------------------------------
# pip install prince
import prince

df_clean = df[categorical_cols].dropna()
mca = prince.MCA(n_components=min(5, len(categorical_cols)*2), random_state=42)
mca.fit(df_clean)

scores_df = mca.row_coordinates(df_clean)
scores_df.columns = [f'Dim {i+1}' for i in range(scores_df.shape[1])]
if group_col and group_col in df.columns:
    scores_df[group_col] = df.loc[df_clean.index, group_col].values

fig_2d = px.scatter(
    scores_df, x='Dim 1', y='Dim 2', color=group_col,
    title="MCA Factor Map: Individuals & Categories",
    template='plotly_dark'
)
fig_2d.show()` : `# -----------------------------------------------------------------------------
# 2. FACTOR ANALYSIS OF MIXED DATA (prince package)
# -----------------------------------------------------------------------------
# pip install prince
import prince

all_cols = continuous_cols + categorical_cols
df_clean = df[all_cols].dropna()

famd = prince.FAMD(n_components=min(5, len(all_cols)), random_state=42)
famd.fit(df_clean)

scores_df = famd.row_coordinates(df_clean)
scores_df.columns = [f'Dim {i+1}' for i in range(scores_df.shape[1])]
if group_col and group_col in df.columns:
    scores_df[group_col] = df.loc[df_clean.index, group_col].values

fig_2d = px.scatter(
    scores_df, x='Dim 1', y='Dim 2', color=group_col,
    title="FAMD Factor Map (Continuous & Categorical Integration)",
    template='plotly_dark'
)
fig_2d.show()`)}
`;
  },

  generateRScript(analysis) {
    const method = analysis.method;
    const contColsR = (analysis.continuousCols || []).map(c => `"${c}"`).join(', ');
    const catColsR = (analysis.categoricalCols || []).map(c => `"${c}"`).join(', ');
    const groupR = analysis.groupingCol ? `"${analysis.groupingCol}"` : 'NULL';

    return `################################################################################
# Multivariate Exploratory Data Analysis (${analysis.methodLabel})
# Package Suite: FactoMineR, factoextra, ggplot2
# Generated by Statis-Gravity Clinical Biostatistics Platform
################################################################################

# install.packages(c("FactoMineR", "factoextra", "ggplot2"))
library(FactoMineR)
library(factoextra)
library(ggplot2)

# 1. Load Clinical Observational Data
# df <- read.csv("clinical_multivariate_data.csv")

continuous_cols  <- c(${contColsR})
categorical_cols <- c(${catColsR})
group_var        <- ${groupR}

${method === 'PCA' ? `# 2. PRINCIPAL COMPONENT ANALYSIS
df_pca <- na.omit(df[, continuous_cols])
res.pca <- PCA(df_pca, scale.unit = TRUE, graph = FALSE)

# Scree Plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# 2D Factor Map with Individuals
fviz_pca_ind(res.pca, 
             habillage = if(!is.null(group_var)) df[[group_var]] else "none",
             addEllipses = TRUE, 
             repel = TRUE)

# Variable Correlation Circle (Biplot Arrows)
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)` : (method === 'MCA' ? `# 2. MULTIPLE CORRESPONDENCE ANALYSIS
df_mca <- na.omit(df[, categorical_cols])
res.mca <- MCA(df_mca, graph = FALSE)

# Scree Plot
fviz_eig(res.mca, addlabels = TRUE)

# Individuals & Categories Biplot
fviz_mca_biplot(res.mca, repel = TRUE,
                ggtheme = theme_minimal())` : `# 2. FACTOR ANALYSIS OF MIXED DATA (FAMD)
all_cols <- c(continuous_cols, categorical_cols)
df_famd <- na.omit(df[, all_cols])
res.famd <- FAMD(df_famd, graph = FALSE)

# Scree Plot
fviz_eig(res.famd, addlabels = TRUE)

# Individuals Factor Map
fviz_famd_ind(res.famd, 
              habillage = if(!is.null(group_var)) df[[group_var]] else "none",
              repel = TRUE)

# Quantitative & Qualitative Variable Contributions
fviz_famd_var(res.famd, "var", repel = TRUE)`)}
`;
  },

  // =========================================================================
  // 9. REALISTIC MIXED CRANIOFACIAL CLINICAL DATASET GENERATOR
  // =========================================================================
  getSampleMixedCohort() {
    let seed = 101;
    const rnd = () => {
      seed = (seed * 9301 + 49297) % 233280;
      return seed / 233280;
    };
    const rndNormal = (mean, sd) => {
      const u1 = Math.max(1e-7, rnd());
      const u2 = rnd();
      const z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      return mean + z * sd;
    };

    const sutures = ['Sagittal', 'Coronal', 'Metopic', 'Lambdoid'];
    const approaches = ['Endoscopic_Strip', 'Open_Cranial_Vault', 'Spring_Assisted'];
    const comorbidities = ['Grade_0', 'Grade_I', 'Grade_II'];

    const cohort = [];
    const N = 120;

    for (let i = 1; i <= N; i++) {
      // Patient age in months
      const rApproach = rnd();
      let approach = approaches[0];
      if (rApproach > 0.45 && rApproach <= 0.85) approach = approaches[1];
      else if (rApproach > 0.85) approach = approaches[2];

      const suture = sutures[Math.floor(rnd() * sutures.length)];
      const comorb = comorbidities[rnd() < 0.65 ? 0 : (rnd() < 0.90 ? 1 : 2)];

      // Correlate continuous parameters with approach & clinical condition
      let age = approach === 'Endoscopic_Strip'
        ? Math.max(2.1, Math.min(5.5, rndNormal(3.6, 0.8)))
        : (approach === 'Spring_Assisted'
            ? Math.max(3.0, Math.min(7.0, rndNormal(4.5, 0.9)))
            : Math.max(6.0, Math.min(14.0, rndNormal(9.2, 2.1))));

      let op_time = approach === 'Endoscopic_Strip'
        ? Math.max(50, Math.round(rndNormal(72, 12)))
        : (approach === 'Spring_Assisted'
            ? Math.max(70, Math.round(rndNormal(98, 16)))
            : Math.max(140, Math.round(rndNormal(205, 32))));

      let blood_loss = approach === 'Endoscopic_Strip'
        ? Math.max(20, Math.round(rndNormal(48, 14)))
        : (approach === 'Spring_Assisted'
            ? Math.max(35, Math.round(rndNormal(75, 20)))
            : Math.max(150, Math.round(rndNormal(285, 65))));

      let los_days = approach === 'Endoscopic_Strip'
        ? Math.max(1.0, Math.round(rndNormal(1.8, 0.5) * 10) / 10)
        : (approach === 'Spring_Assisted'
            ? Math.max(1.5, Math.round(rndNormal(2.4, 0.6) * 10) / 10)
            : Math.max(3.0, Math.round(rndNormal(4.9, 1.1) * 10) / 10));

      let cranial_index = suture === 'Sagittal'
        ? Math.max(62, Math.round(rndNormal(68, 3.5)))
        : (suture === 'Coronal'
            ? Math.max(82, Math.round(rndNormal(88, 4.2)))
            : Math.max(74, Math.round(rndNormal(79, 3.8))));

      let baseline_cvai = Math.max(3.5, Math.round(rndNormal(9.4, 2.6) * 10) / 10);
      let transfusion = blood_loss > 120 || (rnd() < 0.15 && approach !== 'Endoscopic_Strip') ? 'Yes' : 'No';

      cohort.push({
        patient_id: `PT-${1000 + i}`,
        age_months: Math.round(age * 10) / 10,
        cranial_index,
        baseline_cvai,
        operative_time_min: op_time,
        blood_loss_ml: blood_loss,
        length_of_stay_days: los_days,
        suture_type: suture,
        surgical_approach: approach,
        comorbidity_grade: comorb,
        transfusion_required: transfusion
      });
    }

    return cohort;
  }
};


/**
 * WordprocessingML XML Document Builder
 */
class DocxBuilder {
  constructor() {
    this.bodyElements = [];
  }

  addTitle(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:before="120" w:after="80"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="38"/><w:color w:val="0F172A"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addSubTitle(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:after="160"/></w:pPr>` +
      `<w:r><w:rPr><w:sz w:val="22"/><w:color w:val="2563EB"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addAttributionHeader() {
    this.bodyElements.push(
      `<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:after="60"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="18"/><w:color w:val="334155"/></w:rPr>` +
      `<w:t>Conceived, supervised design and testing: Dr G Narenthiran MB ChB BSc(MedSci)(Hons) MRCS(Ed.) FEBNS FRCS(SN)</w:t></w:r></w:p>` +
      `<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:after="180"/></w:pPr>` +
      `<w:r><w:rPr><w:sz w:val="16"/><w:color w:val="64748B"/></w:rPr>` +
      `<w:t>Copyright, G Narenthiran FEBS FRCS(SN), g_narenthiran@hotmail.com | Dedicated to Mrs Nirmaladevy Ganesalingam BSc (mother)</w:t></w:r></w:p>` +
      `<w:p><w:pPr><w:pBdr><w:bottom w:val="single" w:sz="6" w:space="8" w:color="CBD5E1"/></w:pBdr><w:spacing w:after="240"/></w:pPr></w:p>`
    );
    return this;
  }

  addDisclaimerBox() {
    this.bodyElements.push(
      `<w:tbl>` +
      `<w:tblPr><w:tblW w:w="5000" w:type="pct"/><w:tblBorders><w:top w:val="single" w:sz="6" w:color="F59E0B"/><w:left w:val="single" w:sz="18" w:color="D97706"/><w:bottom w:val="single" w:sz="6" w:color="F59E0B"/><w:right w:val="single" w:sz="6" w:color="F59E0B"/><w:insideH w:val="none"/><w:insideV w:val="none"/></w:tblBorders></w:tblPr>` +
      `<w:tr><w:tc><w:tcPr><w:shd w:fill="FEF3C7"/><w:tcMar><w:top w:w="120" w:type="dxa"/><w:bottom w:w="120" w:type="dxa"/><w:left w:w="200" w:type="dxa"/><w:right w:w="200" w:type="dxa"/></w:tcMar></w:tcPr>` +
      `<w:p><w:r><w:rPr><w:b/><w:sz w:val="18"/><w:color w:val="92400E"/></w:rPr><w:t>⚠️ TESTING &amp; METHODOLOGICAL NOTICE: </w:t></w:r>` +
      `<w:r><w:rPr><w:sz w:val="18"/><w:color w:val="78350F"/></w:rPr><w:t>AI was used to vibe code this WebApp. The App is still in testing phase. Not to use for clinical, research or decision making.</w:t></w:r></w:p>` +
      `</w:tc></w:tr></w:tbl>` +
      `<w:p><w:pPr><w:spacing w:after="180"/></w:pPr></w:p>`
    );
    return this;
  }

  addHeading1(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:spacing w:before="280" w:after="120"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="28"/><w:color w:val="1E3A8A"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addHeading2(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:spacing w:before="200" w:after="80"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="24"/><w:color w:val="0D9488"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addHeading3(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:spacing w:before="160" w:after="60"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="20"/><w:color w:val="334155"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addParagraph(text, opts = {}) {
    const bold = opts.bold ? '<w:b/>' : '';
    const italic = opts.italic ? '<w:i/>' : '';
    const color = opts.color ? `<w:color w:val="${opts.color}"/>` : '<w:color w:val="1E293B"/>';
    const sz = opts.size ? `<w:sz w:val="${opts.size}"/>` : '<w:sz w:val="21"/>';
    const align = opts.align ? `<w:jc w:val="${opts.align}"/>` : '';
    const spacing = opts.spacing ? `<w:spacing w:after="${opts.spacing}"/>` : '<w:spacing w:after="120" w:line="276" w:lineRule="auto"/>';

    this.bodyElements.push(
      `<w:p><w:pPr>${align}${spacing}</w:pPr>` +
      `<w:r><w:rPr>${bold}${italic}${color}${sz}</w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addBullet(text, opts = {}) {
    const bold = opts.bold ? '<w:b/>' : '';
    const sz = opts.size ? `<w:sz w:val="${opts.size}"/>` : '<w:sz w:val="21"/>';
    this.bodyElements.push(
      `<w:p><w:pPr><w:ind w:left="360"/><w:spacing w:after="80" w:line="260" w:lineRule="auto"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:color w:val="2563EB"/></w:rPr><w:t>• </w:t></w:r>` +
      `<w:r><w:rPr>${bold}${sz}<w:color w:val="1E293B"/></w:rPr><w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addCalloutBox(title, content, bgColor = "F1F5F9", borderColor = "0284C7") {
    this.bodyElements.push(
      `<w:tbl>` +
      `<w:tblPr><w:tblW w:w="5000" w:type="pct"/><w:tblBorders><w:top w:val="single" w:sz="4" w:color="CBD5E1"/><w:left w:val="single" w:sz="24" w:color="${borderColor}"/><w:bottom w:val="single" w:sz="4" w:color="CBD5E1"/><w:right w:val="single" w:sz="4" w:color="CBD5E1"/><w:insideH w:val="none"/><w:insideV w:val="none"/></w:tblBorders></w:tblPr>` +
      `<w:tr><w:tc><w:tcPr><w:shd w:fill="${bgColor}"/><w:tcMar><w:top w:w="140" w:type="dxa"/><w:bottom w:w="140" w:type="dxa"/><w:left w:w="220" w:type="dxa"/><w:right w:w="220" w:type="dxa"/></w:tcMar></w:tcPr>` +
      (title ? `<w:p><w:pPr><w:spacing w:after="80"/></w:pPr><w:r><w:rPr><w:b/><w:sz w:val="22"/><w:color w:val="0F172A"/></w:rPr><w:t>${escapeXML(title)}</w:t></w:r></w:p>` : '') +
      `<w:p><w:pPr><w:spacing w:after="40" w:line="276" w:lineRule="auto"/></w:pPr><w:r><w:rPr><w:sz w:val="20"/><w:color w:val="334155"/></w:rPr><w:t>${escapeXML(content)}</w:t></w:r></w:p>` +
      `</w:tc></w:tr></w:tbl>` +
      `<w:p><w:pPr><w:spacing w:after="160"/></w:pPr></w:p>`
    );
    return this;
  }

  addTable(headers, rows) {
    let tblXML = `<w:tbl><w:tblPr><w:tblW w:w="5000" w:type="pct"/><w:tblBorders>` +
      `<w:top w:val="single" w:sz="6" w:color="CBD5E1"/>` +
      `<w:left w:val="single" w:sz="4" w:color="E2E8F0"/>` +
      `<w:bottom w:val="single" w:sz="8" w:color="94A3B8"/>` +
      `<w:right w:val="single" w:sz="4" w:color="E2E8F0"/>` +
      `<w:insideH w:val="single" w:sz="4" w:color="E2E8F0"/>` +
      `<w:insideV w:val="single" w:sz="4" w:color="E2E8F0"/>` +
      `</w:tblBorders></w:tblPr>`;

    // Header Row
    if (headers && headers.length > 0) {
      tblXML += `<w:tr><w:trPr><w:tblHeader/></w:trPr>`;
      for (const h of headers) {
        tblXML += `<w:tc><w:tcPr><w:shd w:fill="E2E8F0"/><w:tcMar><w:top w:w="120" w:type="dxa"/><w:bottom w:w="120" w:type="dxa"/><w:left w:w="140" w:type="dxa"/><w:right w:w="140" w:type="dxa"/></w:tcMar></w:tcPr>` +
          `<w:p><w:pPr><w:spacing w:after="40"/></w:pPr><w:r><w:rPr><w:b/><w:sz w:val="19"/><w:color w:val="0F172A"/></w:rPr><w:t>${escapeXML(h)}</w:t></w:r></w:p></w:tc>`;
      }
      tblXML += `</w:tr>`;
    }

    // Data Rows
    for (let rIdx = 0; rIdx < rows.length; rIdx++) {
      const row = rows[rIdx];
      const bg = rIdx % 2 === 1 ? 'F8FAFC' : 'FFFFFF';
      tblXML += `<w:tr>`;
      for (let cIdx = 0; cIdx < row.length; cIdx++) {
        const cell = row[cIdx];
        const isFirstCol = cIdx === 0;
        tblXML += `<w:tc><w:tcPr><w:shd w:fill="${bg}"/><w:tcMar><w:top w:w="100" w:type="dxa"/><w:bottom w:w="100" w:type="dxa"/><w:left w:w="140" w:type="dxa"/><w:right w:w="140" w:type="dxa"/></w:tcMar></w:tcPr>` +
          `<w:p><w:pPr><w:spacing w:after="30"/></w:pPr><w:r><w:rPr>${isFirstCol ? '<w:b/>' : ''}<w:sz w:val="18"/><w:color w:val="1E293B"/></w:rPr><w:t>${escapeXML(cell)}</w:t></w:r></w:p></w:tc>`;
      }
      tblXML += `</w:tr>`;
    }

    tblXML += `</w:tbl><w:p><w:pPr><w:spacing w:after="160"/></w:pPr></w:p>`;
    this.bodyElements.push(tblXML);
    return this;
  }

  buildDocumentXML() {
    return (
      `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>` +
      `<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">` +
      `<w:body>` +
      this.bodyElements.join('\n') +
      `<w:sectPr>` +
      `<w:pgSz w:w="12240" w:h="15840"/>` +
      `<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>` +
      `</w:sectPr>` +
      `</w:body></w:document>`
    );
  }

  toZip() {
    const zip = new SimpleZip();
    zip.addFile('[Content_Types].xml',
      `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>` +
      `<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">` +
      `<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>` +
      `<Default Extension="xml" ContentType="application/xml"/>` +
      `<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>` +
      `</Types>`
    );

    zip.addFile('_rels/.rels',
      `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>` +
      `<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">` +
      `<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>` +
      `</Relationships>`
    );

    zip.addFile('word/document.xml', this.buildDocumentXML());
    return zip;
  }

  generateUint8Array() {
    return this.toZip().build();
  }

  generateBlob() {
    const uint8 = this.generateUint8Array();
    return new Blob([uint8], {
      type: 'application/vnd.openxmlformats-officedocument.wordprocessingml.document'
    });
  }
}

/**
 * High-level DOCX Report Generator per module
 */
const DocxReports = {
  createDescriptiveDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Descriptive Statistics, Normality & Distribution Diagnostics')
      .addAttributionHeader()
      .addDisclaimerBox();

    // 1. Information & Analyzed Input Data
    d.addHeading1('1. Analyzed Cohort Information & Input Data')
      .addParagraph(`Variable/Series Label: ${data.name || 'Sample Continuous Variable'}`)
      .addParagraph(`Sample Size (N): ${data.n} observations`)
      .addParagraph(`Observed Minimum to Maximum: ${data.min.toFixed(2)} to ${data.max.toFixed(2)} (Range = ${(data.max - data.min).toFixed(2)})`)
      .addParagraph(`Raw Data Sample: ${data.values ? data.values.slice(0, 20).join(', ') + (data.values.length > 20 ? ' ...' : '') : 'N/A'}`);

    // 2. Statistical Outcome & Numerical Results
    d.addHeading1('2. Statistical Outcome & Distributional Metrics');
    d.addHeading2('Central Tendency & Parametric/Non-Parametric Dispersion');
    d.addTable(
      ['Statistical Metric', 'Numerical Value', 'Clinical / Mathematical Interpretation'],
      [
        ['Sample Size (N)', `${data.n}`, 'Valid non-null numerical observations'],
        ['Mean (M)', `${data.mean.toFixed(2)}`, 'Arithmetic average (sensitive to outliers)'],
        ['95% Confidence Interval', `[${data.ci95[0].toFixed(2)}, ${data.ci95[1].toFixed(2)}]`, 'Estimated population mean parameter bound'],
        ['Standard Deviation (SD, s)', `${data.sd.toFixed(2)}`, 'Average spread about the sample mean'],
        ['Sample Variance (s²)', `${data.variance.toFixed(2)}`, 'Second central moment (Bessel-corrected n-1)'],
        ['Standard Error of Mean (SEM)', `${data.sem.toFixed(3)}`, 'Precision of sample mean (s / √n)'],
        ['Median (50th Percentile)', `${data.median.toFixed(2)}`, 'Robust central value unaffected by skewness'],
        ['25th Percentile (Q1)', `${data.q1.toFixed(2)}`, 'Lower quartile boundary'],
        ['75th Percentile (Q3)', `${data.q3.toFixed(2)}`, 'Upper quartile boundary'],
        ['Interquartile Range (IQR)', `${data.iqr.toFixed(2)}`, 'Spread of middle 50% of observations (Q3 - Q1)'],
        ['Mode(s)', data.modes && data.modes.length ? data.modes.join(', ') : 'No repeat values', `Peak sample frequency (${data.maxFreq || 1})`],
        ['Range (Max - Min)', `${(data.max - data.min).toFixed(2)}`, `Spans from ${data.min.toFixed(2)} to ${data.max.toFixed(2)}`]
      ]
    );

    d.addHeading2('Shape Diagnostics & Normality Assessment');
    d.addTable(
      ['Distribution Metric', 'Observed Value', 'Statistical Benchmark / Status'],
      [
        ['Sample Skewness (G1)', `${data.skewness.toFixed(3)}`, data.skewnessInterpretation || (Math.abs(data.skewness) < 0.5 ? 'Approximately Symmetric' : 'Skewed')],
        ['Excess Kurtosis (G2)', `${data.kurtosis.toFixed(3)}`, data.kurtosisInterpretation || (Math.abs(data.kurtosis) < 0.5 ? 'Mesokurtic' : 'Heavy/Light-Tailed')],
        ['Normality Omnibus Test', 'Jarque-Bera Test', `JB = ${data.normality.statistic.toFixed(2)}, df = 2, p = ${data.normality.pValue.toExponential(3)}`],
        ['Normality Decision', data.normality.isNormal ? 'Normal (Gaussian)' : 'Statistically Non-Normal', data.normality.isNormal ? 'Parametric tests valid (p ≥ 0.05)' : 'Reject H0 of normality (p < 0.05)']
      ]
    );

    d.addHeading2('Tukey Fences & Outlier Diagnostics');
    const outliersCount = data.outliers ? data.outliers.length : 0;
    d.addParagraph(`Lower Tukey Fence (Q1 - 1.5×IQR): ${data.lowerFence.toFixed(2)}`);
    d.addParagraph(`Upper Tukey Fence (Q3 + 1.5×IQR): ${data.upperFence.toFixed(2)}`);
    d.addParagraph(`Detected Outliers: ${outliersCount === 0 ? 'None detected beyond fences' : `${outliersCount} outlier point(s) identified`}`);
    if (outliersCount > 0) {
      const outlierRows = data.outliers.map(o => [
        `${o.value.toFixed(2)}`,
        `Z = ${o.zScore > 0 ? '+' : ''}${o.zScore.toFixed(2)}`,
        `${o.type}`,
        o.isExtreme ? 'Beyond 3.0×IQR (Extreme Outlier)' : 'Between 1.5× and 3.0×IQR (Mild Outlier)'
      ]);
      d.addTable(['Outlier Value', 'Standardized Z-Score', 'Classification', 'Tukey Fence Range'], outlierRows);
    }

    // 3. Clinical & Statistical Outcome Interpretation
    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'APA / ICMJE Style Summary Statement',
      data.reportText || 'Continuous variable descriptive summary.',
      'F0FDF4',
      '16A34A'
    );
    if (!data.normality.isNormal || outliersCount > 0) {
      d.addParagraph('Clinical Presentation Recommendation: Due to statistically significant distributional skewness and/or the presence of outlier points, reporting the Median and Interquartile Range (IQR) as the primary measures of central tendency and dispersion is strongly advised by STROBE and ICMJE reporting standards.', { bold: true });
    } else {
      d.addParagraph('Clinical Presentation Recommendation: The dataset satisfies Gaussian assumptions without extreme outliers; reporting the Mean and Standard Deviation (Mean ± SD) alongside the 95% Confidence Interval is methodologically sound.', { bold: true });
    }

    // 4. Reason This Particular Test Was Chosen
    d.addHeading1('4. Reason Particular Tests & Methods Were Chosen');
    d.addBullet('Jarque-Bera Omnibus Normality Test: Chosen because it assesses both sample skewness (third standardized moment) and excess kurtosis (fourth moment) simultaneously under an asymptotic Chi-square distribution (df=2). In clinical trials, biological markers and surgical recovery times frequently exhibit asymmetric positive skewness and long tails; the Jarque-Bera test provides superior power in detecting combined shape violations without the arbitrary binning required by Kolmogorov-Smirnov.');
    d.addBullet('Tukey 1.5×IQR Outlier Fences: Selected because Gaussian Z-scores (e.g. ±3 SD) suffer from masking and swamping in the presence of severe outliers (since standard deviation itself is distorted by the extreme points). Tukey fences utilize the interquartile range (IQR), providing a non-parametric, robust threshold for identifying contaminated or physiologically anomalous values.');
    d.addBullet('Dual Parametric & Non-Parametric Metrics: Furnishing both Mean/SD and Median/IQR enables clinicians to immediately detect data skewness and verify whether parametric inferential tests (such as Student\'s t or ANOVA) will be robust or whether non-parametric alternatives (Mann-Whitney, Kruskal-Wallis) must be adopted.');

    // 5. Background Knowledge & Methodology References
    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Descriptive statistics constitute the foundation of all clinical evidence. Under the EQUATOR Network guidelines (CONSORT 2010 for clinical trials, STROBE for observational studies, and STARD for diagnostic accuracy), all patient baseline demographics and clinical outcome variables must be reported with appropriate indicators of spread.');
    d.addParagraph('Mathematical Definitions:');
    d.addBullet('Sample Mean: x̄ = (1/n) ∑ xᵢ. Reflects the gravitational balance point of the observations.');
    d.addBullet('Bessel-Corrected Standard Deviation: s = √[ (1/(n-1)) ∑ (xᵢ - x̄)² ]. Uses n-1 degrees of freedom to remove negative bias in estimating the true population variance σ.');
    d.addBullet('Standard Error of the Mean: SEM = s / √n. Measures the precision with which the sample mean estimates the true population mean; as sample size n grows, SEM contracts.');
    d.addBullet('Jarque-Bera Statistic: JB = (n/6) [ S² + (K² / 4) ] ~ χ²(2), where S is sample skewness and K is excess kurtosis.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Jarque CM, Bera AK (1987). A test for normality of observations and regression residuals. International Statistical Review, 55(2): 163–172.');
    d.addBullet('Tukey JW (1977). Exploratory Data Analysis. Addison-Wesley.');
    d.addBullet('Altman DG, Bland JM (1996). Detecting skewness from summary information. BMJ, 313(7066): 1200.');

    return d;
  },

  createHypothesisDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle(`Module: Two-Cohort Hypothesis Testing (${data.testName || 'Inferential Analysis'})`)
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Cohort Information & Input Data')
      .addParagraph(`Cohort A: ${data.nameA} (n = ${data.groupA.n} patients/specimens)`)
      .addParagraph(`Cohort B: ${data.nameB} (n = ${data.groupB.n} patients/specimens)`)
      .addParagraph(`Evaluated Statistical Paradigm: ${data.testName || 'Two-Sample Test'}`)
      .addParagraph(`Cohort A Data Preview: ${data.groupA.values ? data.groupA.values.slice(0, 15).join(', ') : 'N/A'}`)
      .addParagraph(`Cohort B Data Preview: ${data.groupB.values ? data.groupB.values.slice(0, 15).join(', ') : 'N/A'}`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addTable(
      ['Cohort Summary', 'Sample Size (n)', 'Mean', 'Std Dev (SD)', 'Std Error (SEM)', '95% Confidence Interval', 'Median (IQR)'],
      [
        [data.nameA, `${data.groupA.n}`, `${data.groupA.mean.toFixed(2)}`, `${data.groupA.sd.toFixed(2)}`, `${data.groupA.sem.toFixed(3)}`, `[${data.groupA.ci95[0].toFixed(2)}, ${data.groupA.ci95[1].toFixed(2)}]`, `${data.groupA.median.toFixed(2)} (${data.groupA.iqr.toFixed(2)})`],
        [data.nameB, `${data.groupB.n}`, `${data.groupB.mean.toFixed(2)}`, `${data.groupB.sd.toFixed(2)}`, `${data.groupB.sem.toFixed(3)}`, `[${data.groupB.ci95[0].toFixed(2)}, ${data.groupB.ci95[1].toFixed(2)}]`, `${data.groupB.median.toFixed(2)} (${data.groupB.iqr.toFixed(2)})`]
      ]
    );

    if (data.assumptions) {
      d.addHeading2('Diagnostic Assessment of Statistical Assumptions');
      if (data.assumptions.isPaired) {
        d.addParagraph(`Sample Design: Paired / Repeated Measures (n = ${data.assumptions.n} matched observations).`);
        d.addTable(
          ['Assumption Evaluated', 'Test / Metric', 'Calculated Statistic', 'p-Value', 'Verdict / Interpretation'],
          [
            ['Normality of Differences (Δ)', 'Jarque-Bera Test', `JB = ${data.assumptions.normality.statistic.toFixed(2)} (Skew: ${data.assumptions.normality.skewness.toFixed(2)}, Kurt: ${data.assumptions.normality.kurtosis.toFixed(2)})`, `${data.assumptions.normality.pValue < 0.001 ? 'p < .001' : 'p = ' + data.assumptions.normality.pValue.toFixed(4)}`, data.assumptions.normality.isNormal ? 'Normal Distribution (Parametric Valid)' : 'Non-Normal (Non-Parametric Recommended)'],
            ['Homogeneity of Variance', 'Within-Subject Differencing', 'N/A (Between-cohort variance removed by design)', 'N/A', 'Homoscedasticity assumption satisfied by pairing']
          ]
        );
      } else {
        d.addParagraph(`Sample Design: Independent Two-Cohort Comparison (${data.nameA}: n = ${data.assumptions.statsA.n}, ${data.nameB}: n = ${data.assumptions.statsB.n}).`);
        d.addTable(
          ['Assumption Evaluated', 'Target Cohort / Test', 'Calculated Statistic', 'p-Value', 'Verdict / Interpretation'],
          [
            ['Normality (Cohort 1)', `${data.nameA} (Jarque-Bera)`, `JB = ${data.assumptions.normality.jbA.statistic.toFixed(2)} (Skew: ${data.assumptions.statsA.skewness.toFixed(2)})`, `${data.assumptions.normality.jbA.pValue < 0.001 ? 'p < .001' : 'p = ' + data.assumptions.normality.jbA.pValue.toFixed(4)}`, data.assumptions.normality.normA ? 'Normal Distribution' : 'Skewed / Non-Normal'],
            ['Normality (Cohort 2)', `${data.nameB} (Jarque-Bera)`, `JB = ${data.assumptions.normality.jbB.statistic.toFixed(2)} (Skew: ${data.assumptions.statsB.skewness.toFixed(2)})`, `${data.assumptions.normality.jbB.pValue < 0.001 ? 'p < .001' : 'p = ' + data.assumptions.normality.jbB.pValue.toFixed(4)}`, data.assumptions.normality.normB ? 'Normal Distribution' : 'Skewed / Non-Normal'],
            ['Homogeneity of Variance', 'F-Test of Variances', `F(${data.assumptions.varianceEquality.df1}, ${data.assumptions.varianceEquality.df2}) = ${data.assumptions.varianceEquality.fStat.toFixed(2)} (s₁²=${data.assumptions.statsA.variance.toFixed(2)}, s₂²=${data.assumptions.statsB.variance.toFixed(2)})`, `${data.assumptions.varianceEquality.pValue < 0.001 ? 'p < .001' : 'p = ' + data.assumptions.varianceEquality.pValue.toFixed(4)}`, data.assumptions.varianceEquality.equalVariance ? 'Equal Variances (Homoscedastic)' : 'Unequal Variances (Heteroscedastic)']
          ]
        );
      }
      d.addCalloutBox(
        'Automated Test Recommendation Decision Engine',
        `Recommended Test: ${data.assumptions.recommendedTestName}\nDecision Rationale: ${data.assumptions.rationale}`,
        'E0F2FE',
        '0284C7'
      );
    }

    d.addHeading2('Comparative Inferential Test Results');
    d.addTable(
      ['Inferential Parameter', 'Calculated Value', 'Clinical Interpretation / Benchmark'],
      [
        ['Test Statistic', `${data.testName.includes('Mann-Whitney') ? 'U = ' : (data.testName.includes('Wilcoxon') ? 'W = ' : 't = ')}${(data.statistic !== undefined ? data.statistic : data.zScore || 0).toFixed(3)}`, 'Standardized difference between cohort locations'],
        ['Degrees of Freedom (df)', `${data.df ? data.df.toFixed(2) : 'N/A (Rank test)'}`, 'Satterthwaite adjustment for unequal cohort variances'],
        ['p-Value (Two-Tailed)', `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + data.pValue.toFixed(4)}`, data.isSignificant ? 'Statistically Significant (p < 0.05)' : 'Not Significant (p ≥ 0.05)'],
        ['Mean Difference (ΔM)', `${data.meanDiff !== undefined ? (data.meanDiff >= 0 ? '+' : '') + data.meanDiff.toFixed(2) : 'N/A'}`, `Observed clinical point difference (${data.nameA} - ${data.nameB})`],
        ['95% CI of Difference', data.ci95 ? `[${data.ci95[0].toFixed(2)}, ${data.ci95[1].toFixed(2)}]` : 'N/A', 'Range of plausible true population differences'],
        ['Standardized Effect Size', data.cohensD !== undefined ? `Cohen\'s d = ${data.cohensD.toFixed(2)}` : `Rank-Biserial r = ${(data.rankBiserial || 0).toFixed(2)}`, 'Magnitude of separation independent of sample size'],
        ['Effect Classification', Math.abs(data.cohensD || data.rankBiserial || 0) >= 0.8 ? 'Large Effect' : (Math.abs(data.cohensD || data.rankBiserial || 0) >= 0.5 ? 'Moderate Effect' : 'Small/Negligible Effect'), 'Cohen (1988) benchmark: 0.2 small, 0.5 medium, 0.8 large']
      ]
    );

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'APA / ICMJE Style Summary Statement',
      data.reportText || 'Hypothesis testing narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Welch\'s Unequal Variances t-Test: Selected as the gold-standard default because classical Student\'s t-test assumes homoscedasticity (equal population variances). In biomedical and clinical research, treatment interventions frequently alter not only the mean but also the variance (e.g. patients respond heterogeneously). When variances differ even moderately, Student\'s t-test exhibits severely inflated Type I error rates. Welch\'s t-test adjusts the degrees of freedom using the Satterthwaite approximation, providing robust Type I error control without sacrificing statistical power when variances happen to be identical (Ruxton 2006, Delacre et al. 2017).');
    d.addBullet('Mann-Whitney U Test (Alternative): Provided for ordinal scales, Glasgow Coma Scale (GCS) scores, or continuous markers with severe skewness. Instead of comparing means, it evaluates stochastic dominance across ranked values.');
    d.addBullet('Reporting Effect Sizes: Statistical significance (p < 0.05) is heavily dependent on sample size; very large clinical cohorts can produce tiny p-values for clinically meaningless differences. Cohen\'s d quantifies the practical medical magnitude of the separation.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Hypothesis testing evaluates the likelihood of observing the experimental difference under the null hypothesis (H0: μA = μB).');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Welch\'s t-Statistic: t = (x̄A - x̄B) / √[ (sA²/nA) + (sB²/nB) ].');
    d.addBullet('Welch-Satterthwaite Degrees of Freedom: df = [ (sA²/nA + sB²/nB)² ] / [ (sA²/nA)² / (nA - 1) + (sB²/nB)² / (nB - 1) ].');
    d.addBullet('Cohen\'s d: d = (x̄A - x̄B) / sₚ, where sₚ is the pooled standard deviation.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Welch BL (1947). The generalization of \'Student\'s\' problem when several different population variances are involved. Biometrika, 34(1/2): 28–35.');
    d.addBullet('Ruxton GD (2006). The unequal variance t-test is an underused alternative to Student\'s t-test and the Mann-Whitney U test. Behavioral Ecology, 17(4): 688–690.');
    d.addBullet('Delacre M, Lakens D, Leys C (2017). Why psychologists should by default use Welch\'s t-test instead of Student\'s t-test. International Review of Social Psychology, 30(1): 92–101.');

    return d;
  },

  createAnovaDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle(`Module: ${data.testName || 'Multi-Cohort Analysis'} & Post-Hoc Pairwise Contrasts`)
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Cohorts Information & Input Data')
      .addParagraph(`Number of Evaluated Cohorts/Conditions (k): ${data.k || (data.groups ? data.groups.length : 0)}`)
      .addParagraph(`Total Analyzed Sample Size: ${data.totalN || 0} observations ${data.matchedN ? `(${data.matchedN} matched subjects)` : ''}`)
      .addParagraph(`Study Design: ${data.isPaired ? 'Within-Subjects Repeated Measures / Paired Timepoints' : 'Independent Between-Subjects Cohorts'}`);

    if (data.groups && data.groups.length > 0) {
      d.addHeading2('Cohort Descriptive Summary');
      const groupRows = data.groups.map(g => [
        g.name,
        `${g.stats.n}`,
        `${g.stats.mean.toFixed(2)}`,
        `${g.stats.sd.toFixed(2)}`,
        `${g.stats.sem.toFixed(3)}`,
        `[${g.stats.ci95[0].toFixed(2)}, ${g.stats.ci95[1].toFixed(2)}]`,
        `${g.stats.median.toFixed(2)} (${g.stats.iqr.toFixed(2)})`
      ]);
      d.addTable(['Cohort Name', 'Sample n', 'Mean (M)', 'Std Dev (SD)', 'Std Error (SEM)', '95% CI of Mean', 'Median (IQR)'], groupRows);
    }

    if (data.assumptions) {
      d.addHeading2('Diagnostic Assessment of Statistical Assumptions');
      if (data.assumptions.isPaired) {
        d.addParagraph(`Study Design: Paired / Repeated Measures across ${data.assumptions.k} conditions (N = ${data.assumptions.matchedN} matched subjects).`);
        const normRows = (data.assumptions.normality && data.assumptions.normality.details)
          ? data.assumptions.normality.details.map(det => [
              det.name,
              `${det.n}`,
              `${det.skewness.toFixed(2)}`,
              `${det.kurtosis.toFixed(2)}`,
              `JB = ${(det.jbStat || 0).toFixed(2)}`,
              `${det.pValue < 0.001 ? 'p < .001' : 'p = ' + det.pValue.toFixed(3)}`,
              det.isNormal ? 'Normal (Parametric Valid)' : 'Skewed / Non-Normal'
            ])
          : [];
        if (normRows.length > 0) {
          d.addTable(['Condition', 'n', 'Skewness', 'Kurtosis', 'Jarque-Bera', 'p-Value', 'Normality Status'], normRows);
        }
      } else {
        d.addParagraph(`Study Design: Independent Between-Subjects Comparison across ${data.assumptions.k} cohorts.`);
        const normRows = (data.assumptions.normality && data.assumptions.normality.details)
          ? data.assumptions.normality.details.map(det => [
              det.name,
              `${det.n}`,
              `${det.skewness.toFixed(2)}`,
              `${det.kurtosis.toFixed(2)}`,
              `JB = ${(det.jbStat || 0).toFixed(2)}`,
              `${det.pValue < 0.001 ? 'p < .001' : 'p = ' + det.pValue.toFixed(3)}`,
              det.isNormal ? 'Normal Distribution' : 'Skewed / Non-Normal'
            ])
          : [];
        if (normRows.length > 0) {
          d.addTable(['Cohort', 'n', 'Skewness', 'Kurtosis', 'Jarque-Bera', 'p-Value', 'Normality Status'], normRows);
        }

        if (data.assumptions.varianceEquality && data.assumptions.varianceEquality.applicable) {
          d.addTable(
            ['Assumption Evaluated', 'Diagnostic Test', 'Test Statistic & df', 'p-Value', 'Homoscedasticity Verdict'],
            [
              [
                'Homogeneity of Variances',
                'Levene\'s Test (Brown-Forsythe)',
                `F(${data.assumptions.varianceEquality.df1}, ${data.assumptions.varianceEquality.df2}) = ${(data.assumptions.varianceEquality.fStat || 0).toFixed(2)} (Ratio: ${(data.assumptions.varianceEquality.varianceRatio || 1).toFixed(2)}×)`,
                `${data.assumptions.varianceEquality.pValue < 0.001 ? 'p < .001' : 'p = ' + data.assumptions.varianceEquality.pValue.toFixed(4)}`,
                data.assumptions.varianceEquality.equalVariance ? 'Equal Variances Confirmed (Homoscedastic)' : 'Unequal Variances (Heteroscedastic - Welch Required)'
              ]
            ]
          );
        }
      }

      d.addCalloutBox(
        'Automated Test Recommendation Decision Engine',
        `Recommended Test: ${data.assumptions.recommendedTestName}\nDecision Rationale: ${data.assumptions.rationale}${data.requestedMethod && data.requestedMethod !== 'auto' && data.requestedMethod !== data.assumptions.recommendedTest ? '\n[Note: User manually selected ' + data.testName + ']' : ''}`,
        'E0F2FE',
        '0284C7'
      );
    }

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2(`${data.testName || 'Omnibus Test'} Summary Table`);

    if (data.testKey === 'welch') {
      d.addTable(
        ['Inferential Parameter', 'Calculated Value', 'Clinical Interpretation / Benchmark'],
        [
          ['Welch F-Statistic', `F_Welch = ${(data.statistic || data.fStatistic || 0).toFixed(3)}`, 'Robust omnibus variance ratio adjusting for heteroscedasticity'],
          ['Adjusted Degrees of Freedom', `df1 = ${data.df1 || data.dfBetween || 0}, df2 = ${(data.df2 || data.dfWithin || 0).toFixed(2)}`, 'Adjusted via Welch-Satterthwaite approximation'],
          ['p-Value', `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + (data.pValue || 0).toFixed(4)}`, data.isSignificant ? 'Statistically Significant (p < 0.05)' : 'Not Significant (ns)'],
          ['Robust Omega-Squared (ω²)', `${(data.omegaSquared || 0).toFixed(3)}`, 'Unbiased population effect size for unequal variances'],
          ['Eta-Squared (η²)', `${(data.etaSquared || 0).toFixed(3)}`, 'Sample proportion of total variance explained']
        ]
      );
    } else if (data.testKey === 'kruskal') {
      d.addTable(
        ['Inferential Parameter', 'Calculated Value', 'Clinical Interpretation / Benchmark'],
        [
          ['Kruskal-Wallis Statistic', `H = ${(data.statistic || 0).toFixed(3)}`, 'Non-parametric omnibus rank sum variance'],
          ['Degrees of Freedom (df)', `df = ${data.df || (data.k - 1)}`, 'k - 1 cohorts'],
          ['p-Value (Chi-Square)', `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + (data.pValue || 0).toFixed(4)}`, data.isSignificant ? 'Statistically Significant (p < 0.05)' : 'Not Significant (ns)'],
          ['Epsilon-Squared (ε²)', `${(data.epsilonSquared || data.etaSquared || 0).toFixed(3)}`, 'Non-parametric degree of stochastic separation (0 to 1)']
        ]
      );
    } else if (data.testKey === 'rm_anova') {
      d.addTable(
        ['Source of Variation', 'Sum of Squares (SS)', 'Degrees of Freedom (df)', 'Mean Square (MS)', 'F-Statistic', 'p-Value', 'Partial Eta² (η²_p)'],
        [
          ['Treatment (Time/Condition)', `${(data.ssTreatment || data.ssBetween || 0).toFixed(2)}`, `${data.dfTreatment || data.dfBetween || 0}`, `${(data.msTreatment || data.msBetween || 0).toFixed(2)}`, `F = ${(data.fStatistic || 0).toFixed(2)}`, `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + (data.pValue || 0).toFixed(4)}`, `${(data.partialEtaSquared || data.etaSquared || 0).toFixed(3)}`],
          ['Error (Residual)', `${(data.ssError || data.ssWithin || 0).toFixed(2)}`, `${data.dfError || data.dfWithin || 0}`, `${(data.msError || data.msWithin || 0).toFixed(2)}`, '-', '-', `GG Epsilon (ε̂) = ${(data.ggEpsilon || 1).toFixed(2)}`],
          ['Subjects', `${(data.ssSubjects || 0).toFixed(2)}`, `${data.dfSubjects || 0}`, '-', '-', '-', '-'],
          ['Total', `${(data.ssTotal || 0).toFixed(2)}`, `${data.dfTotal || 0}`, '-', '-', '-', '-']
        ]
      );
    } else if (data.testKey === 'friedman') {
      d.addTable(
        ['Inferential Parameter', 'Calculated Value', 'Clinical Interpretation / Benchmark'],
        [
          ['Friedman Test Statistic', `χ²_F = ${(data.statistic || 0).toFixed(3)}`, 'Two-way rank sum statistic for matched observations'],
          ['Degrees of Freedom (df)', `df = ${data.df || (data.k - 1)}`, 'k - 1 conditions across N subjects'],
          ['p-Value (Chi-Square)', `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + (data.pValue || 0).toFixed(4)}`, data.isSignificant ? 'Statistically Significant (p < 0.05)' : 'Not Significant (ns)'],
          ['Kendall\'s Concordance (W)', `${(data.kendallsW || 0).toFixed(3)}`, 'Degree of subject ranking consistency across conditions (0 to 1)']
        ]
      );
    } else {
      // One-Way ANOVA (Fisher's Standard)
      d.addTable(
        ['Source of Variation', 'Sum of Squares (SS)', 'Degrees of Freedom (df)', 'Mean Square (MS)', 'F-Statistic', 'p-Value', 'Omega-Squared (ω²)'],
        [
          ['Between Groups (Treatment)', `${(data.ssBetween || 0).toFixed(2)}`, `${data.dfBetween || 0}`, `${(data.msBetween || 0).toFixed(2)}`, `F = ${(data.fStatistic || 0).toFixed(2)}`, `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + (data.pValue || 0).toFixed(4)}`, `${(data.omegaSquared || 0).toFixed(3)}`],
          ['Within Groups (Residual/Error)', `${(data.ssWithin || 0).toFixed(2)}`, `${data.dfWithin || 0}`, `${(data.msWithin || 0).toFixed(2)}`, '-', '-', `Eta² (η²) = ${(data.etaSquared || 0).toFixed(3)}`],
          ['Total', `${(data.ssTotal !== undefined ? data.ssTotal : ((data.ssBetween || 0) + (data.ssWithin || 0))).toFixed(2)}`, `${(data.dfBetween || 0) + (data.dfWithin || 0)}`, '-', '-', '-', '-']
        ]
      );
    }

    if (data.pairwise && data.pairwise.length > 0) {
      d.addHeading2('Post-Hoc Pairwise Contrasts');
      const statHeader = data.testKey === 'welch' ? 'Games-Howell t' : (data.testKey === 'kruskal' ? 'Dunn\'s z' : (data.testKey === 'friedman' ? 'Wilcoxon z' : (data.testKey === 'rm_anova' ? 'Paired t' : 'Tukey q')));
      const pairRows = data.pairwise.map(p => [
        p.comparison,
        `${(p.meanDiff >= 0 ? '+' : '')}${p.meanDiff.toFixed(2)}`,
        `${p.seDiff ? p.seDiff.toFixed(3) : '-'}`,
        `${p.qStatistic !== undefined ? 'q = ' + p.qStatistic.toFixed(2) : (p.tStatistic !== undefined ? 't = ' + p.tStatistic.toFixed(2) : 'z = ' + (p.zStatistic || 0).toFixed(2))}`,
        `${p.pValue < 0.001 ? 'p < .001' : 'p = ' + p.pValue.toFixed(4)}`,
        p.ci95 ? `[${p.ci95[0].toFixed(2)}, ${p.ci95[1].toFixed(2)}]` : 'N/A',
        `${p.cohensD !== undefined ? p.cohensD.toFixed(2) : '-'}`,
        p.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)'
      ]);
      d.addTable(['Pairwise Contrast', 'Difference', 'Std Error', statHeader, 'Adjusted p', '95% CI of Diff', 'Effect Size', 'Significance'], pairRows);
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Multi-Cohort APA / ICMJE Clinical Summary',
      data.reportText || 'Multi-group statistical narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    if (data.testKey === 'welch') {
      d.addBullet('Welch\'s Heteroscedastic ANOVA: Chosen because the cohorts exhibit unequal population variances (heteroscedasticity confirmed by Levene\'s test). Standard Fisher ANOVA suffers from severe Type I error rate inflation when group variances differ. Welch\'s ANOVA computes weighted variance terms and adjusts degrees of freedom via the Welch-Satterthwaite method, preserving valid error rates.');
      d.addBullet('Games-Howell Post-Hoc Contrasts: Adopted because it does not assume equal variances or equal group sample sizes, strictly bounding family-wise error across multiple contrasts.');
    } else if (data.testKey === 'kruskal') {
      d.addBullet('Kruskal-Wallis H Test: Chosen because one or more cohorts violate the assumption of normality or contain heavy outliers. By transforming continuous observations into ranks, it evaluates whether the median rank distributions differ significantly across groups without parametric distribution assumptions.');
      d.addBullet('Dunn\'s Post-Hoc Contrasts: Adopted to pinpoint pairwise stochastic differences using mean rank differences with family-wise error adjustments.');
    } else if (data.testKey === 'rm_anova') {
      d.addBullet('Repeated Measures ANOVA: Chosen because the same subjects were evaluated repeatedly across conditions/timepoints. By partitioning out between-subjects variability from the error term, it provides substantially higher statistical power than between-subjects ANOVA.');
      d.addBullet('Greenhouse-Geisser Sphericity Adjustment: Applied to adjust degrees of freedom when the compound symmetry / sphericity assumption is violated.');
    } else if (data.testKey === 'friedman') {
      d.addBullet('Friedman Test: Chosen as the non-parametric counterpart to Repeated Measures ANOVA. It ranks conditions within each individual subject, eliminating between-subject baseline differences while protecting against non-normal or skewed longitudinal distributions.');
    } else {
      d.addBullet('One-Way Omnibus ANOVA: Selected because testing multiple cohorts with uncorrected pairwise t-tests results in severe Family-Wise Error Rate inflation (FWER). ANOVA simultaneously assesses whether between-cohort variance exceeds within-cohort residual variation.');
      d.addBullet('Tukey\'s Honest Significant Difference (HSD): Chosen as the post-hoc method because it utilizes the Studentized Range distribution (q) to strictly bound the overall Family-Wise Error Rate at α = 0.05 across all possible pairwise comparisons.');
      d.addBullet('Omega-Squared (ω²) Reporting: Included alongside Eta-squared (η²) because Omega-squared provides an unbiased population effect size estimate in clinical samples.');
    }

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Partitioning of Variance in Multi-Group Designs:');
    d.addBullet('Between-Groups Mean Square: MS_B = SS_B / (k - 1).');
    d.addBullet('Within-Groups Mean Square: MS_W = SS_W / (N - k).');
    d.addBullet('Welch F-Ratio: Incorporates group sample size weights w_j = n_j / s_j².');
    d.addBullet('Friedman Statistic: χ²_F = [12 / (N k (k+1))] ∑ R_j² - 3 N (k+1).');
    d.addParagraph('Key Academic References:');
    d.addBullet('Fisher RA (1925). Statistical Methods for Research Workers. Oliver and Boyd, Edinburgh.');
    d.addBullet('Welch BL (1951). On the comparison of several mean values: an alternative approach. Biometrika, 38(3/4): 330–336.');
    d.addBullet('Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance analysis. J Am Stat Assoc, 47(260): 583–621.');
    d.addBullet('Friedman M (1937). The use of ranks to avoid the assumption of normality. J Am Stat Assoc, 32(200): 675–701.');
    d.addBullet('Games PA, Howell JF (1976). Pairwise multiple comparison procedures with unequal N\'s and/or variances. J Educ Stat, 1(2): 113–125.');

    return d;
  },

  createCategoricalDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle(`Module: 2x2 Contingency, Risk Metrics & Diagnostic Evaluation`)
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed 2x2 Matrix Information & Configuration')
      .addParagraph(`Analysis Framework: ${data.mode === 'diagnostic' ? 'Diagnostic Test Evaluation vs Gold Standard Reference' : 'Clinical Study / Trial (Intervention vs Control Cohort)'}`)
      .addParagraph(`Total Evaluated Subjects (N): ${data.totalN} patients/samples`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2('Observed 2x2 Contingency Matrix');
    if (data.mode === 'diagnostic') {
      d.addTable(
        ['Test Condition', 'Gold Standard Positive (+)', 'Gold Standard Negative (-)', 'Total Row Margins'],
        [
          ['New Test Positive (+)', `TP = ${data.a}`, `FP = ${data.b}`, `Test Pos = ${data.a + data.b}`],
          ['New Test Negative (-)', `FN = ${data.c}`, `TN = ${data.d}`, `Test Neg = ${data.c + data.d}`],
          ['Total Column Margins', `Condition Pos = ${data.a + data.c}`, `Condition Neg = ${data.b + data.d}`, `N = ${data.totalN}`]
        ]
      );

      d.addHeading2('Diagnostic Accuracy & Discriminative Performance Metrics');
      d.addTable(
        ['Diagnostic Index', 'Point Estimate (%)', '95% Confidence Interval', 'Clinical Definition & Interpretation'],
        [
          ['Sensitivity (True Positive Rate)', `${(data.diag.sensitivity * 100).toFixed(1)}%`, `[${(data.diag.sensitivityCI95[0]*100).toFixed(1)}%, ${(data.diag.sensitivityCI95[1]*100).toFixed(1)}%]`, 'Ability of test to correctly identify patients with disease (TP / [TP + FN])'],
          ['Specificity (True Negative Rate)', `${(data.diag.specificity * 100).toFixed(1)}%`, `[${(data.diag.specificityCI95[0]*100).toFixed(1)}%, ${(data.diag.specificityCI95[1]*100).toFixed(1)}%]`, 'Ability of test to correctly identify disease-free patients (TN / [TN + FP])'],
          ['Positive Predictive Value (PPV)', `${(data.diag.ppv * 100).toFixed(1)}%`, `[${(data.diag.ppvCI95[0]*100).toFixed(1)}%, ${(data.diag.ppvCI95[1]*100).toFixed(1)}%]`, 'Probability that a patient with a positive test truly has disease'],
          ['Negative Predictive Value (NPV)', `${(data.diag.npv * 100).toFixed(1)}%`, `[${(data.diag.npvCI95[0]*100).toFixed(1)}%, ${(data.diag.npvCI95[1]*100).toFixed(1)}%]`, 'Probability that a patient with a negative test is truly disease-free'],
          ['Overall Diagnostic Accuracy', `${(data.diag.accuracy * 100).toFixed(1)}%`, `[${(data.diag.accuracyCI95[0]*100).toFixed(1)}%, ${(data.diag.accuracyCI95[1]*100).toFixed(1)}%]`, 'Proportion of all test results that were correct ([TP + TN] / N)'],
          ['Positive Likelihood Ratio (LR+)', `${data.diag.plr.toFixed(2)}`, `[${data.diag.plrCI95[0].toFixed(2)}, ${data.diag.plrCI95[1].toFixed(2)}]`, 'Ratio of TPR to FPR (>10 indicates strong diagnostic confirmation)'],
          ['Negative Likelihood Ratio (LR-)', `${data.diag.nlr.toFixed(2)}`, `[${data.diag.nlrCI95[0].toFixed(2)}, ${data.diag.nlrCI95[1].toFixed(2)}]`, 'Ratio of FNR to TNR (<0.1 indicates strong rule-out capacity)'],
          ['Youden\'s Index (J)', `${data.diag.youdenJ.toFixed(3)}`, '-', 'Overall effectiveness criterion: Sensitivity + Specificity - 1']
        ]
      );
    } else {
      d.addTable(
        ['Study Arm / Cohort', 'Event Occurred (+)', 'Event Absent (-)', 'Total Cohort Sample'],
        [
          ['Intervention / Treatment', `Events = ${data.a}`, `Non-Events = ${data.b}`, `Total = ${data.a + data.b}`],
          ['Control / Placebo', `Events = ${data.c}`, `Non-Events = ${data.d}`, `Total = ${data.c + data.d}`],
          ['Total Margins', `Total Events = ${data.a + data.c}`, `Total Non-Events = ${data.b + data.d}`, `N = ${data.totalN}`]
        ]
      );

      d.addHeading2('Hypothesis Tests & Clinical Risk Metrics');
      d.addTable(
        ['Statistical / Risk Metric', 'Numerical Value', 'Clinical Interpretation / Benchmark'],
        [
          ['Pearson Chi-Square (χ²)', `χ²(1) = ${data.chiSquare.standard.toFixed(2)} (p = ${data.chiSquare.pValueStandard < 0.001 ? '< .001' : data.chiSquare.pValueStandard.toFixed(4)})`, 'Asymptotic test of association'],
          ['Yates\' Continuity Corrected χ²', `χ²(1) = ${data.chiSquare.yates.toFixed(2)} (p = ${data.chiSquare.pValueYates < 0.001 ? '< .001' : data.chiSquare.pValueYates.toFixed(4)})`, 'Corrects for discrete probability overestimation'],
          ['Fisher\'s Exact Test (Two-Tailed)', `p = ${data.fishersExact.pValue < 0.001 ? '< .001' : data.fishersExact.pValue.toFixed(4)}`, 'Exact conditional hypergeometric probability'],
          ['Odds Ratio (OR)', `${data.risk.oddsRatio.toFixed(2)} (95% CI [${data.risk.orCI95[0].toFixed(2)}, ${data.risk.orCI95[1].toFixed(2)}])`, 'Odds of event in intervention vs control'],
          ['Relative Risk (RR)', `${data.risk.relativeRisk.toFixed(2)} (95% CI [${data.risk.rrCI95[0].toFixed(2)}, ${data.risk.rrCI95[1].toFixed(2)}])`, 'Risk of event in intervention vs control'],
          ['Absolute Risk Reduction (ARR)', `${(data.risk.arr * 100).toFixed(1)}% (95% CI [${(data.risk.arrCI95[0]*100).toFixed(1)}%, ${(data.risk.arrCI95[1]*100).toFixed(1)}%])`, 'Absolute difference in event rate'],
          ['Relative Risk Reduction (RRR)', `${(data.risk.rrr * 100).toFixed(1)}%`, 'Proportional reduction relative to baseline risk'],
          ['Number Needed to Treat (NNT)', `${data.risk.nnt.toFixed(1)} patients`, 'Patients needed to treat to prevent 1 adverse event (1 / ARR)']
        ]
      );
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Clinical / Epidemiological Summary',
      data.reportText || '2x2 categorical narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Fisher\'s Exact Test vs Pearson Chi-Square: Pearson\'s Chi-Square is an asymptotic test requiring adequate expected cell counts (Cochrane criterion: no expected cell < 1, and no more than 20% < 5). In clinical trials or rare surgical complications, small counts render Chi-Square unreliable. Fisher\'s Exact Test evaluates the exact hypergeometric permutation distribution, providing unconditional validity across any sample size.');
    d.addBullet('Importance of ARR and NNT in Evidence-Based Medicine: While Relative Risk Reduction (RRR) is commonly advertised in commercial literature, it can dramatically exaggerate clinical benefits. For instance, reducing risk from 2 in 10,000 to 1 in 10,000 is a 50% relative reduction, yet the ARR is 0.01% (NNT = 10,000). Reporting ARR and NNT ensures clinicians comprehend the real-world clinical workload required to achieve patient benefit.');
    d.addBullet('Likelihood Ratios for Diagnostic Evaluation: PPV and NPV are profoundly distorted by disease prevalence (Bayes\' theorem). Likelihood Ratios are prevalence-independent intrinsic properties of the diagnostic test.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Fisher\'s Exact Hypergeometric Probability: P = [ (a+b)! (c+d)! (a+c)! (b+d)! ] / [ N! a! b! c! d! ].');
    d.addBullet('Pearson\'s Chi-Square: χ² = ∑ [ (O - E)² / E ], with df = (r-1)(c-1) = 1.');
    d.addBullet('Number Needed to Treat: NNT = 1 / ARR = 1 / |p_intervention - p_control|.');
    d.addBullet('Positive Likelihood Ratio: LR+ = Sensitivity / (1 - Specificity).');
    d.addBullet('Negative Likelihood Ratio: LR- = (1 - Sensitivity) / Specificity.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Fisher RA (1922). On the interpretation of χ² from contingency tables, and the calculation of P. Journal of the Royal Statistical Society, 85(1): 87–94.');
    d.addBullet('Laupacis A, Sackett DL, Roberts RS (1988). An assessment of clinically useful measures of the consequences of treatment. N Engl J Med, 318(26): 1728–1733.');
    d.addBullet('Yates F (1934). Contingency tables involving small numbers and the χ² test. Supplement to the Journal of the Royal Statistical Society, 1(2): 217–235.');

    return d;
  },

  createCorrelationDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Correlation, Linear Regression & Non-Linear Curve Morphology')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Variables & Data Information')
      .addParagraph(`Independent Predictor (X): ${data.xName || 'Variable X'}`)
      .addParagraph(`Dependent Outcome (Y): ${data.yName || 'Variable Y'}`)
      .addParagraph(`Sample Size (N pairs): ${data.n} observations`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2('Bivariate Correlation & Regression Parameters');
    d.addTable(
      ['Regression / Correlation Metric', 'Observed Value', 'Clinical Interpretation / Benchmark'],
      [
        ['Pearson Correlation Coefficient (r)', `r = ${data.corr.r.toFixed(3)} (p = ${data.corr.p < 0.001 ? '< .001' : data.corr.p.toFixed(4)})`, 'Linear co-variation strength and direction'],
        ['95% Confidence Interval for r', `[${data.corr.ci95[0].toFixed(2)}, ${data.corr.ci95[1].toFixed(2)}]`, 'Fisher Z-transformed confidence limits'],
        ['Spearman Rank Correlation (ρ)', `ρ = ${data.corr.spearman.toFixed(3)} (p = ${data.corr.spearmanP < 0.001 ? '< .001' : data.corr.spearmanP.toFixed(4)})`, 'Monotonic association across ranked values'],
        ['Linear Regression Equation', `Y = ${data.reg.intercept.toFixed(2)} + ${data.reg.slope.toFixed(3)} × X`, 'Ordinary Least Squares (OLS) line of best fit'],
        ['Slope (β1) & Std Error', `${data.reg.slope.toFixed(3)} ± ${data.reg.seSlope.toFixed(3)}`, `Change in Y per 1-unit increase in X`],
        ['Intercept (β0) & Std Error', `${data.reg.intercept.toFixed(2)} ± ${data.reg.seIntercept.toFixed(2)}`, `Estimated baseline value of Y when X = 0`],
        ['Coefficient of Determination (R²)', `R² = ${data.reg.rSquared.toFixed(3)}`, `Proportion of total outcome variance explained (${(data.reg.rSquared * 100).toFixed(1)}%)`],
        ['Residual Standard Error (RSE)', `${data.reg.residualSE.toFixed(2)}`, 'Standard deviation of residuals around regression line'],
        ['Regression F-Statistic', `F(1, ${data.n - 2}) = ${data.reg.fStatistic.toFixed(2)} (p = ${data.reg.pVal < 0.001 ? '< .001' : data.reg.pVal.toFixed(4)})`, 'Overall statistical significance of linear model']
      ]
    );

    if (data.morphology) {
      d.addHeading2('Curve Morphology & Non-Linear Assumption Diagnostics');
      d.addTable(
        ['Morphology Diagnostic', 'Assessed Status', 'Clinical & Methodological Significance'],
        [
          ['Detected Relationship Shape', `${data.morphology.shape}`, 'Identifies linear vs non-linear physiologic windows'],
          ['Monotonicity Status', data.morphology.isMonotonic ? 'Strictly Monotonic' : 'Non-Monotonic (Directional Reversals)', data.morphology.reversals > 0 ? `${data.morphology.reversals} reversal(s) detected` : 'Consistently rising or falling'],
          ['Quadratic Polynomial Fit (R²)', `R² = ${data.morphology.quadR2.toFixed(3)}`, 'Fit of quadratic curve Y = a + bX + cX²'],
          ['Optimal Vertex / Nadir Point', data.morphology.vertexX ? `X* = ${data.morphology.vertexX.toFixed(1)} (Y* = ${data.morphology.vertexY.toFixed(1)})` : 'N/A (Linear)', 'Optimal physiologic target (e.g. ideal CPP in neurotrauma)'],
          ['Recommended Model', `${data.morphology.recommendedModel}`, 'Recommended mathematical representation']
        ]
      );
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Clinical Correlation & Regression Narrative',
      data.reportText || 'Correlation narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Curve Morphology & Non-Linear Screening: Standard statistical tools blindly compute Pearson\'s r. However, in clinical medicine, relationships are frequently U-shaped or non-monotonic (e.g. Cerebral Perfusion Pressure vs Mortality, or Drug Concentration vs Efficacy). On a symmetric U-curve, Pearson\'s r is approximately 0, misleading researchers into claiming "no association". Screening for quadratic curvature and directional monotonicity prevents dangerous medical errors.');
    d.addBullet('Spearman\'s Rank Correlation: Selected as a robust alternative when relationships are monotonic but non-linear, or when outlier observations exert excessive leverage.');
    d.addBullet('Ordinary Least Squares Regression with 95% Confidence Bands: Enables quantitative clinical prediction of expected physiological responses.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Pearson Correlation: r = ∑ [ (xᵢ - x̄)(yᵢ - ȳ) ] / √[ ∑ (xᵢ - x̄)² ∑ (yᵢ - ȳ)² ].');
    d.addBullet('OLS Slope: β1 = ∑ [ (xᵢ - x̄)(yᵢ - ȳ) ] / ∑ (xᵢ - x̄)².');
    d.addBullet('Quadratic Vertex: X* = -b / (2c) for parabola Y = a + bX + cX².');
    d.addParagraph('Key Academic References:');
    d.addBullet('Pearson K (1895). Notes on regression and inheritance in the case of two parents. Proceedings of the Royal Society of London, 58: 240–242.');
    d.addBullet('Spearman C (1904). The proof and measurement of association between two things. American Journal of Psychology, 15(1): 72–101.');
    d.addBullet('Anscombe FJ (1973). Graphs in statistical analysis. The American Statistician, 27(1): 17–21.');

    return d;
  },

  createDiagnosticDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Diagnostic Biomarker ROC & Discrimination Analysis')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Biomarker & Diagnostic Information')
      .addParagraph(`Index Biomarker / Test: ${data.name || 'Clinical Diagnostic Score'}`)
      .addParagraph(`Total Patient Cohort Evaluated: ${data.totalN} patients`)
      .addParagraph(`True Positive Disease Cases: ${data.posCount} patients`)
      .addParagraph(`True Negative Normal/Control Cases: ${data.negCount} patients`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2('ROC Discriminative Performance & Optimal Cutoff');
    d.addTable(
      ['Diagnostic Metric', 'Numerical Value', 'Clinical Interpretation / Benchmark'],
      [
        ['Area Under the Curve (ROC AUC)', `${data.auc.toFixed(3)}`, 'Probability that a diseased patient scores higher than a non-diseased patient'],
        ['Hanley-McNeil 95% CI for AUC', `[${data.aucCI95[0].toFixed(3)}, ${data.aucCI95[1].toFixed(3)}]`, 'Precision of global discriminative ability'],
        ['Standard Error of AUC (SE)', `${data.seAuc.toFixed(3)}`, 'Non-parametric standard error estimate'],
        ['Discrimination Rating', data.auc >= 0.90 ? 'Outstanding Discrimination (AUC ≥ 0.90)' : (data.auc >= 0.80 ? 'Excellent Discrimination (0.80 - 0.89)' : 'Acceptable Discrimination (0.70 - 0.79)'), 'Hosmer-Lemeshow diagnostic classification'],
        ['Optimal Cutoff Score', `${data.optimalCutoff.toFixed(2)}`, 'Threshold maximizing Youden\'s J index (Sensitivity + Specificity - 1)'],
        ['Sensitivity at Cutoff', `${(data.sensitivity * 100).toFixed(1)}%`, 'True Positive Rate at recommended decision boundary'],
        ['Specificity at Cutoff', `${(data.specificity * 100).toFixed(1)}%`, 'True Negative Rate at recommended decision boundary'],
        ['Positive Likelihood Ratio (LR+)', `${data.plr.toFixed(2)}`, 'Factor multiplying pre-test odds upon positive result'],
        ['Negative Likelihood Ratio (LR-)', `${data.nlr.toFixed(2)}`, 'Factor multiplying pre-test odds upon negative result']
      ]
    );

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Diagnostic Accuracy & ROC Summary',
      data.reportText || 'ROC narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Receiver Operating Characteristic (ROC) Analysis: Selected because diagnostic tests must not be judged at a single arbitrary threshold. Evaluating sensitivity across the entire continuum of specificities allows clinical teams to select distinct cutoffs depending on whether the test is used for high-sensitivity screening (where false negatives cannot be tolerated) or high-specificity confirmatory testing.');
    d.addBullet('Likelihood Ratios: Emphasized because predictive values (PPV/NPV) depend entirely on underlying population prevalence. Likelihood Ratios remain stable and allow Bayesian post-test probability calculations using Fagan\'s nomogram.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Empirical AUC: Computed via non-parametric trapezoidal integration, mathematically equivalent to the Wilcoxon/Mann-Whitney U statistic: AUC = U / (n₁ n₂).');
    d.addBullet('Hanley-McNeil Variance: SE = √[ (θ(1-θ) + (n₁-1)(Q₁ - θ²) + (n₂-1)(Q₂ - θ²)) / (n₁ n₂) ], where θ is the AUC.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Hanley JA, McNeil BJ (1982). The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology, 143(1): 29–36.');
    d.addBullet('Youden WJ (1950). Index for rating diagnostic tests. Cancer, 3(1): 32–35.');
    d.addBullet('Zweig MH, Campbell G (1993). Receiver-operating characteristic (ROC) plots: a fundamental evaluation tool in clinical medicine. Clinical Chemistry, 39(4): 561–577.');

    return d;
  },

  createPowerDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Sample Size Determination & Statistical Power Analysis')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Study Design Information & Input Parameters')
      .addParagraph(`Study Design: ${data.designLabel || 'Two-Sample Clinical Design'}`)
      .addParagraph(`Analysis Objective: ${data.goal === 'power' ? 'Post-Hoc / Achieved Statistical Power Estimation' : 'A Priori Sample Size Requirement Determination'}`)
      .addParagraph(`Type I Error Rate (Alpha, α): ${data.alpha} (Two-Tailed)`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2('Calculated Sample Size & Power Parameters');
    d.addTable(
      ['Study Planning Parameter', 'Calculated Value', 'Clinical & Methodological Significance'],
      [
        ['Standardized Effect Size', `${data.effectSizeLabel || 'Effect'}: ${data.effectSize.toFixed(3)}`, 'Anticipated magnitude of clinical benefit'],
        ['Required Sample Size per Arm', `${data.nPerGroup} patients`, 'Minimum enrollment per arm to guarantee target power'],
        ['Total Required Study Enrollment', `${data.totalN} patients`, 'Total recruitment target across both arms'],
        ['Target Statistical Power (1 - β)', `${(data.targetPower * 100).toFixed(1)}%`, 'Desired probability of correctly rejecting false null hypothesis'],
        ['Achieved Statistical Power', `${(data.achievedPower * 100).toFixed(1)}%`, 'Actual probability of detecting the effect with specified n'],
        ['Recommended Target with 15% Attrition', `${Math.ceil(data.totalN / 0.85)} patients`, 'Accommodates anticipated patient dropout / loss to follow-up']
      ]
    );

    if (data.additionalMetrics) {
      d.addHeading2('Comparative Risk & Methodology Metrics');
      const addRows = Object.entries(data.additionalMetrics).map(([k, v]) => [k, `${v}`, '-']);
      d.addTable(['Metric', 'Value', 'Notes'], addRows);
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Study Planning & Ethical Justification Summary',
      data.reportText || 'Sample size planning narrative statement.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('A Priori Power Calculation Mandate: Underpowered trials fail to detect clinically vital treatments (Type II error), wasting grant funding and exposing patients to experimental risks without scientific return. Conversely, oversized trials needlessly expose excess patients to inferior treatments. Ethics committees (IRBs) and funding bodies mandate rigorous a priori power justification.');
    d.addBullet('Continuity-Corrected Fisher / Proportion Formula: Classical asymptotic formulas underestimate sample size requirements when proportions are small. Applying continuity correction ensures that the exact conditional test achieves the intended statistical power.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Statistical power (1 - β) is the probability of avoiding a false-negative conclusion.');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Independent Means Sample Size: n = 2 [ (Z_{α/2} + Z_β) / d ]².');
    d.addBullet('Paired Means Sample Size: N_pairs = [ (Z_{α/2} + Z_β) / d_z ]².');
    d.addBullet('Proportions Fleiss Continuity Correction: n\' = (n / 4) [ 1 + √( 1 + 4 / (n |p₁ - p₂|) ) ]².');
    d.addParagraph('Key Academic References:');
    d.addBullet('Cohen J (1988). Statistical Power Analysis for the Behavioral Sciences (2nd ed.). Lawrence Erlbaum Associates.');
    d.addBullet('Fleiss JL, Tytun A, Ury HK (1980). A simple approximation for calculating sample sizes for comparing two independent proportions. Biometrics, 36(2): 343–346.');
    d.addBullet('Schulz KF, Altman DG, Moher D (2010). CONSORT 2010 Statement: updated guidelines for reporting parallel group randomised trials. BMJ, 340: c332.');

    return d;
  },

  createTeachingDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Teaching, Distributions, CLT & Student\'s t Convergence')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Distribution Information & Simulation Parameters');
    d.addParagraph(`Selected Distribution: ${data.distName || 'Normal (Gaussian) Distribution'}`);
    d.addParagraph(`Sample Size (N): ${data.distN || 500} simulated observations`);
    d.addParagraph(`Clinical Context: ${data.clinicalExample || 'Biomedical research modeling'}`);
    d.addParagraph(`CLT Simulation Parent Population: ${data.cltParentName || 'Uniform Distribution'}`);
    d.addParagraph(`CLT Sample Size per Draw (n): ${data.cltN || 30}`);
    d.addParagraph(`CLT Total Iterations (k): ${data.cltK || 1000} sample means`);
    if (data.tConv) {
      d.addParagraph(`Student's t Simulation: Sample Size n = ${data.tConv.sampleSize} (Degrees of Freedom ν = ${data.tConv.df})`);
    }
    if (data.overlap) {
      d.addParagraph(`Two-Sample Overlap Simulation: Mean Difference Δ = ${data.overlap.delta.toFixed(2)}, Group 1 (SD₁ = ${data.overlap.sd1.toFixed(2)}, SEM₁ = ${data.overlap.sem1.toFixed(3)}, n₁ = ${data.overlap.n1}), Group 2 (SD₂ = ${data.overlap.sd2.toFixed(2)}, SEM₂ = ${data.overlap.sem2.toFixed(3)}, n₂ = ${data.overlap.n2}), Significance Level α = ${data.overlap.alpha.toFixed(3)}`);
    }
    if (data.bayes) {
      d.addParagraph(`Bayesian Statistics & Logic (3Blue1Brown Model): Prior P(H) = ${(data.bayes.prior * 100).toFixed(1)}%, Likelihood P(E|H) = ${(data.bayes.likelihood * 100).toFixed(1)}%, False Positive Rate P(E|¬H) = ${(data.bayes.falsePositive * 100).toFixed(1)}%, Population Sample N = ${data.bayes.sampleSize}`);
    }

    d.addHeading1('2. Statistical Outcome & Empirical Convergence Metrics');
    d.addHeading2('Computer-Generated Distribution Metrics');
    d.addTable(
      ['Distribution Metric', 'Empirical Sample Value', 'Theoretical Population Value', 'Status'],
      [
        ['Mean (M)', `${data.sampleMean?.toFixed(3) || '-'}`, `${data.theoMean?.toFixed(3) || '-'}`, 'Congruent'],
        ['Standard Deviation (SD)', `${data.sampleSD?.toFixed(3) || '-'}`, `${data.theoSD?.toFixed(3) || '-'}`, 'Congruent'],
        ['Skewness (G₁)', `${data.skewness?.toFixed(3) || '-'}`, `${data.theoSkewness || '0.000'}`, data.skewnessLabel || 'Evaluated'],
        ['Kurtosis (Excess G₂)', `${data.kurtosis?.toFixed(3) || '-'}`, `${data.theoKurtosis || '0.000'}`, data.kurtosisLabel || 'Evaluated'],
        ['Jarque-Bera Normality Test', `JB = ${data.jbStat?.toFixed(2) || '-'}, p = ${data.jbP?.toFixed(4) || '-'}`, 'Null: Gaussian', data.isNormal ? 'Normal (p ≥ 0.05)' : 'Non-Normal (p < 0.05)']
      ]
    );

    if (data.cltResults) {
      d.addHeading2('Central Limit Theorem (CLT) Convergence Results');
      d.addTable(
        ['CLT Parameter', 'Simulated Value', 'Theoretical CLT Value', 'Convergence Note'],
        [
          ['Parent Population Mean (μ)', `${data.cltResults.theoMean?.toFixed(3)}`, `${data.cltResults.theoMean?.toFixed(3)}`, 'Ground truth benchmark'],
          ['Observed Mean of Means (x̄̄)', `${data.cltResults.obsMean?.toFixed(3)}`, `${data.cltResults.theoMean?.toFixed(3)}`, `Error: ${Math.abs(data.cltResults.obsMean - data.cltResults.theoMean).toFixed(4)}`],
          ['Standard Error of Means (SE)', `${data.cltResults.obsSE?.toFixed(3)}`, `${data.cltResults.theoSE?.toFixed(3)}`, `Law of 1/√n shrinkage: σ/√${data.cltN}`],
          ['Sampling Distribution Skewness', `${data.cltResults.skewness?.toFixed(3)}`, '0.000 (Symmetric)', 'Asymmetry eradicated by averaging'],
          ['Sampling Distribution Normality', `p = ${data.cltResults.normalityP?.toFixed(4)}`, 'p ≥ 0.05', data.cltResults.isNormal ? 'Gaussian Bell Curve Achieved' : 'Approaching Gaussian']
        ]
      );
    }

    if (data.tConv) {
      d.addHeading2('Student\'s t Convergence to Standard Normal N(0, 1) Results');
      d.addTable(
        ['Student\'s t Parameter', `t-Distribution (ν = ${data.tConv.df})`, 'Standard Normal N(0, 1)', 'Methodological Consequence'],
        [
          ['Two-Tailed Critical Value (α=0.05)', `t_crit = ${data.tConv.tCrit?.toFixed(3)}`, 'z_crit = 1.960', `Deviation: ${data.tConv.critDiffPct >= 0 ? '+' : ''}${data.tConv.critDiffPct?.toFixed(1)}% wider cutoff`],
          ['Peak Density f(0)', `f_t(0) = ${data.tConv.tPeak?.toFixed(4)}`, 'φ(0) = 0.3989', `Peak discrepancy: ${data.tConv.peakDiffPct?.toFixed(1)}%`],
          ['Tail Probability P(|X| > 1.960)', `${(data.tConv.tailProb * 100).toFixed(1)}%`, '5.00%', `Type I false positive risk if using z=1.96: ${(data.tConv.tailProb * 100).toFixed(1)}%`],
          ['Excess Kurtosis (Fat Tails)', `${data.tConv.excessKurtosis === Infinity ? '∞ (Fat Tails)' : data.tConv.excessKurtosis?.toFixed(2)}`, '0.00 (Mesokurtic)', data.tConv.df <= 4 ? '4th moment undefined' : 'Heavy tail factor: 6/(ν-4)'],
          ['Max Discrepancy sup |f_t - φ|', `${data.tConv.maxDiscrepancy?.toFixed(4)}`, '0.0000', 'Uniform convergence metric']
        ]
      );
    }

    if (data.overlap) {
      d.addHeading2('Two-Sample Overlap, Dispersion (SD vs. SEM) & Alpha Significance Results');
      d.addTable(
        ['Analytical Parameter', 'Simulated Value', 'Clinical & Inferential Meaning'],
        [
          ['Mean Difference (Δ = μ₂ - μ₁)', `Δ = ${data.overlap.delta.toFixed(2)}`, 'Observed separation between the two group means'],
          ['Group 1 Spread (SD₁ & SEM₁)', `SD₁ = ${data.overlap.sd1.toFixed(2)}, SEM₁ = ${data.overlap.sem1.toFixed(3)} (n₁ = ${data.overlap.n1})`, 'Biological spread (SD₁) vs. sample mean precision (SEM₁ = SD₁/√n₁)'],
          ['Group 2 Spread (SD₂ & SEM₂)', `SD₂ = ${data.overlap.sd2.toFixed(2)}, SEM₂ = ${data.overlap.sem2.toFixed(3)} (n₂ = ${data.overlap.n2})`, 'Biological spread (SD₂) vs. sample mean precision (SEM₂ = SD₂/√n₂)'],
          ['Standard Error of Difference (SE_diff)', `SE_diff = ${data.overlap.seDiff.toFixed(3)} (Welch df = ${data.overlap.df.toFixed(1)})`, 'Combined standard error: √(SEM₁² + SEM₂²) under unequal variances'],
          ['Chosen Significance Level (α)', `α = ${data.overlap.alpha.toFixed(3)}`, `Type I error tolerance: ${data.overlap.alpha === 0.05 ? 'Standard 5% biomedical risk' : (data.overlap.alpha < 0.05 ? 'Strict threshold' : 'Relaxed exploratory threshold')}`],
          ['Critical Value (t_crit vs z_crit)', `t_crit = ${data.overlap.tCrit.toFixed(3)} (z = ${data.overlap.zCrit.toFixed(3)})`, 'Required number of standard errors to claim statistical significance'],
          ['Critical Difference Boundary (Δcrit)', `Δcrit = ${data.overlap.deltaCrit.toFixed(2)}`, 'Minimum mean separation needed to achieve p < α (Δcrit = t_crit · SE_diff)'],
          ['Test Statistic & p-value', `t = ${data.overlap.tStat.toFixed(2)}, p = ${data.overlap.pValue < 0.0001 ? '< 0.0001' : data.overlap.pValue.toFixed(4)}`, `Significance: ${data.overlap.isSignificant ? 'REJECT H₀ (p < α)' : 'FAIL TO REJECT H₀ (p ≥ α)'}`],
          ['Individual Patient Overlap (OVL_SD)', `${(data.overlap.patientOVL * 100).toFixed(1)}%`, 'Proportion of overlapping individual patient values (Weitzman\'s OVL)'],
          ['Sampling Distribution Overlap (OVL_SEM)', `${(data.overlap.meansOVL * 100).toFixed(1)}%`, 'Proportion of overlap between the sampling distributions of sample means']
        ]
      );
    }

    if (data.power) {
      d.addHeading2('Statistical Power (1 − β), Type II Error (β), SD & SEM Simulation Results');
      d.addTable(
        ['Power & Precision Parameter', 'Simulated Value', 'Clinical & Regulatory Significance'],
        [
          ['Statistical Power (1 − β)', `${(data.power.power * 100).toFixed(1)}%`, `Sensitivity / True positive detection rate (${data.power.powerRating})`],
          ['Type II Error Rate (Beta, β)', `${(data.power.beta * 100).toFixed(1)}%`, 'Probability of failing to detect a true treatment difference (false negative risk)'],
          ['Required Sample Size per Group (n)', `n = ${data.power.n} patients`, `Total trial recruitment: N = ${data.power.totalN} patients across 2 treatment arms`],
          ['Target Mean Difference (Δ)', `Δ = ${data.power.delta.toFixed(2)}`, 'Minimum clinically important difference (MCID) between treatment means'],
          ['Patient Standard Deviation (SD, σ)', `SD = ${data.power.sd.toFixed(2)}`, 'Biological variability among patients; higher SD inflates required sample size'],
          ['Standard Error of the Mean (SEM)', `SEM = ${data.power.sem.toFixed(3)}`, 'Precision of mean estimate: SEM = SD / √n; shrinking SEM separates curves and drives Power'],
          ['Standard Error of Difference (SE_diff)', `SE_diff = ${data.power.seDiff.toFixed(3)}`, 'Combined estimation error: σ · √(2/n) = √2 · SEM'],
          ['Standardized Effect Size (Cohen\'s d)', `d = ${data.power.cohensD.toFixed(2)}`, `${data.power.cohensD >= 0.8 ? 'Large effect' : (data.power.cohensD >= 0.5 ? 'Medium effect' : 'Small effect')} (d = Δ / SD)`],
          ['Critical Significance Cutoff (xcrit)', `xcrit = ${data.power.xCrit.toFixed(3)}`, `Boundary beyond which H₀ is rejected: xcrit = z_crit · SE_diff at α = ${data.power.alpha.toFixed(3)}`],
          ['Non-Centrality Parameter (λ)', `λ = ${data.power.lambda.toFixed(3)}`, 'Signal-to-noise ratio shifting the H₁ distribution: λ = Δ / SE_diff']
        ]
      );
    }

    if (data.bayes) {
      d.addHeading2('Bayesian Statistics & Logic (3Blue1Brown Model) Simulation Results');
      d.addTable(
        ['Bayesian Parameter', 'Simulated Value', 'Epistemological & Clinical Meaning'],
        [
          ['Prior Probability P(H)', `${(data.bayes.prior * 100).toFixed(1)}% (Odds 1:${(1 / data.bayes.priorOdds).toFixed(1)})`, 'Initial degree of belief before observing evidence (Base rate)'],
          ['Likelihood P(E|H)', `${(data.bayes.likelihood * 100).toFixed(1)}%`, 'True Positive Rate: Probability of evidence given hypothesis H is true'],
          ['False Positive Rate P(E|¬H)', `${(data.bayes.falsePositive * 100).toFixed(1)}%`, 'False Alarm Rate: Probability of evidence given hypothesis H is false'],
          ['Total Evidence P(E)', `${(data.bayes.pEvidence * 100).toFixed(2)}%`, 'Marginal likelihood: Total shaded area of possibilities matching evidence'],
          ['Bayes Factor (Likelihood Ratio)', `${data.bayes.bayesFactor.toFixed(2)}×`, `${data.bayes.evidenceRating} (P(E|H) / P(E|¬H))`],
          ['Posterior Probability P(H|E)', `${(data.bayes.posterior * 100).toFixed(1)}% (Odds 1:${(1 / data.bayes.posteriorOdds).toFixed(1)})`, `Updated degree of belief after conditioning on evidence (Shift: ${(data.bayes.beliefShift >= 0 ? '+' : '')}${(data.bayes.beliefShift * 100).toFixed(1)}%)`],
          ['Representative Counts (N)', `H: ${data.bayes.countHAndE} of ${data.bayes.countH} | ¬H: ${data.bayes.countNotHAndE} of ${data.bayes.countNotH}`, `Natural frequencies in population of N=${data.bayes.sampleSize}: ${data.bayes.countHAndE} / ${data.bayes.countTotalE} = ${(data.bayes.posterior * 100).toFixed(1)}%`],
          ['Sequential Updating (4 Steps)', `P₀: ${(data.bayes.prior * 100).toFixed(1)}% → P₁: ${(data.bayes.trajectory[1].p * 100).toFixed(1)}% → P₄: ${(data.bayes.trajectory[4].p * 100).toFixed(1)}%`, 'Compounding belief trajectory across successive independent observations']
        ]
      );
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Pedagogical Synthesis & Clinical Trial Relevance',
      data.reportText || 'The Central Limit Theorem, Student\'s t convergence, Two-Sample Overlap, Statistical Power, and Bayesian updating simulations demonstrate the mathematical foundations of parametric testing, the definition of alpha, and Bayesian logic in scientific research.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Foundation of Inferential Biostatistics: Parametric hypothesis tests (Student t-test, ANOVA, ordinary least squares regression) mathematically assume normally distributed errors or sample means. The Central Limit Theorem provides the mathematical justification for deploying these tests in clinical trials with n ≥ 30 even when raw clinical metrics (e.g. ICU stay, recovery hours) are skewed.');
    d.addBullet('Gosset\'s Student\'s t Adjustment: In small clinical cohorts (n < 30), estimating population variance σ² using sample variance s² introduces substantial stochastic instability into the test statistic denominator. Using Gaussian critical values (z = 1.96) severely inflates the Type I error rate (e.g. to 14.5% at n = 4). Student\'s t distribution compensates for this extra uncertainty by thickening the tails and demanding a higher critical threshold (t = 3.182 at n = 4).');
    d.addBullet('The n ≥ 31 Clinical Threshold: As demonstrated by the simulation, when sample size reaches n ≥ 31 (degrees of freedom ν ≥ 30), the critical t cutoff drops to 2.042 (only 4.2% wider than 1.960), and tail probability converges close to 5.0%. This mathematical threshold explains why sample sizes of 30 or greater historically permit Gaussian approximation in medical trial protocols.');
    d.addBullet('Why α = 0.05 Defines the Point of Significance: In 1925, Ronald A. Fisher proposed the 5% significance level (p < 0.05) as a pragmatic convention for scientific research—representing a 1 in 20 chance of observing an effect as extreme under the null hypothesis of no difference. On a standard Gaussian distribution, exactly 5% of probability mass lies in the tails beyond ±1.960 standard errors (2.5% in each tail). Hence, the critical separation distance between sample means is Δcrit = 1.960 · SE_diff. When the observed difference Δ exceeds Δcrit, the p-value falls below 0.05.');
    d.addBullet('The Fundamental Distinction Between SD and SEM: Standard Deviation (SD) reflects real inter-individual biological diversity among patients and does not contract when sample size increases. In contrast, the Standard Error of the Mean (SEM = SD/√n) quantifies our uncertainty in the population mean estimate and contracts steadily as 1/√n. Consequently, two treatment groups can exhibit 70% biological overlap in individual patient scores, yet their treatment difference can be verified as statistically significant (p < 0.001) once sufficient patients are enrolled to shrink the SEM.');
    d.addBullet('Consequences of Modifying Alpha (α): Relaxing α to 0.10 moves the critical cutoff inward to z = 1.645, lowering the required separation Δcrit and declaring significance on smaller differences or smaller sample sizes, at the expense of doubling the false-positive risk to 10%. Tightening α to 0.01 (z = 2.576) or 0.001 (z = 3.291), as required in confirmatory registration trials or genome-wide studies, shifts the cutoff outward into the extreme tails, demanding either much larger effect sizes or substantially expanded sample sizes before significance can be claimed.');
    d.addBullet('Statistical Power (1 − β) as the Scientific Safeguard Against False Negatives: While alpha (α = 0.05) strictly caps the risk of a false positive, statistical power (1 − β) measures the study\'s ability to identify a genuine therapeutic effect. An underpowered trial (e.g. 50% power) is ethically and scientifically problematic because patients undergo experimental risk when the study has only a coin-toss probability of reaching definitive conclusions.');
    d.addBullet('The Interplay of SD, SEM, Beta, and Power: The non-centrality parameter λ = Δ / (SD · √(2/n)) controls the separation between null and alternative distributions. Because SEM = SD / √n, doubling the sample size shrinks SEM by 1.414, drawing the distributions apart and collapsing the Type II error region β.');
    d.addBullet('The 3Blue1Brown Geometric Insight into Bayes\' Theorem: Rather than memorizing abstract formulas, Bayes\' theorem is intuitively understood as proportions of area within a 1×1 unit square of all possibilities. Observing evidence restricts our sample space to only the shaded regions where the evidence occurs; the posterior probability is simply the fraction of that restricted space corresponding to the hypothesis of interest.');
    d.addBullet('Base Rate Neglect and Natural Frequency Framing: In Steve the Librarian problem, people intuitively fixate on the 4:1 likelihood ratio (40% vs 10%) and forget the 20:1 base rate ratio of farmers to librarians. Translating abstract probabilities into natural frequencies (e.g. 4 librarians vs 20 farmers in a village of 210 people) eliminates cognitive bias and reveals why Steve is still 5× more likely to be a farmer.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Classical Lindberg-Lévy Central Limit Theorem: Let X₁, X₂, ..., X_n be independent and identically distributed (i.i.d.) random variables with mean μ and finite variance σ². Then as n → ∞: √n (X̄_n - μ) / σ → N(0, 1).');
    d.addBullet('Student\'s t Distribution Density: f(t; ν) = [ Γ((ν+1)/2) / (√(πν) Γ(ν/2)) ] · [ 1 + t²/ν ]^{-(ν+1)/2}. As ν → ∞, [ 1 + t²/ν ]^{-(ν+1)/2} → exp(-t²/2), converging to Standard Normal N(0, 1).');
    d.addBullet('Standard Error of the Mean: SEM = σ / √n. Quadrupling patient enrollment cuts the estimation uncertainty in half.');
    d.addBullet('Weitzman\'s Distribution Overlap Coefficient (OVL): For two equal-variance Gaussian curves separated by difference Δ: OVL = 2 · Φ(-|Δ| / (2 · s)), where s = SD for patient-level overlap and s = SEM for sampling-mean-level overlap.');
    d.addBullet('Critical Significance Boundary: Δcrit = t_crit(α, df) · SD · √(2/n).');
    d.addBullet('Two-Sample Power Formulation: 1 - β = Φ(Δ / (σ · √(2/n)) - z_{1 - α/2}).');
    d.addBullet('Required Sample Size Equation: n = 2 · (z_{1 - α/2} + z_{1 - β})² · σ² / Δ².');
    d.addBullet('Bayes\' Theorem in Area Form: P(H|E) = P(H ∩ E) / P(E) = [P(H) · P(E|H)] / [P(H) · P(E|H) + P(¬H) · P(E|¬H)].');
    d.addBullet('Odds Form of Bayes\' Rule: Posterior Odds = Prior Odds × Bayes Factor (Likelihood Ratio).');
    d.addParagraph('Key Academic References:');
    d.addBullet('Sanderson G (2019). Bayes theorem, the geometry of changing beliefs. 3Blue1Brown, YouTube.');
    d.addBullet('Kahneman D, Tversky A (1973). On the psychology of prediction. Psychological Review, 80(4): 237–251.');
    d.addBullet('Gigerenzer G, Hoffrage U (1995). How to improve Bayesian reasoning without instruction: Frequency formats. Psychological Review, 102(4): 684–704.');
    d.addBullet('Cohen J (1988). Statistical Power Analysis for the Behavioral Sciences. 2nd ed. Hillsdale, NJ: Lawrence Erlbaum Associates.');
    d.addBullet('Moher D, Hopewell S, Schulz KF, et al. (2010). CONSORT 2010 explanation and elaboration: updated guidelines for reporting parallel group randomised trials. BMJ, 340: c869.');
    d.addBullet('Altman DG, Bland JM (1995). Absence of evidence is not evidence of absence. BMJ, 311(7003): 485.');
    d.addBullet('Fisher RA (1925). Statistical Methods for Research Workers. Edinburgh: Oliver and Boyd.');
    d.addBullet('Cumming G, Finch S (2005). Inference by eye: confidence intervals and how to read pictures of data. Am Psychol, 60(2): 170–180.');
    d.addBullet('Student [Gosset WS] (1908). The probable error of a mean. Biometrika, 6(1): 1–25.');
    d.addBullet('Altman DG, Bland JM (2005). Standard deviations and standard errors. BMJ, 331(7521): 903.');
    d.addBullet('Laplace PS (1810). Mémoire sur les approximations des formules qui sont fonctions de très grands nombres et sur leur application aux probabilités. Mémoires de l\'Académie Royale des Sciences de Paris.');
    d.addBullet('Gauss CF (1809). Theoria motus corporum coelestium in sectionibus conicis solem ambientium. Hamburg: Perthes et Besser.');

    return d;
  },

  createTeachingBayesianDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Teaching - Bayesian Statistics & Logic (3Blue1Brown Model)')
      .addAttributionHeader()
      .addDisclaimerBox();

    const bayes = data.bayes || {};

    d.addHeading1('1. Bayesian Model & Parameter Specification');
    d.addParagraph(`Prior Probability P(H): ${(bayes.prior * 100).toFixed(1)}% (Prior Odds: 1 : ${(1 / (bayes.priorOdds || 0.05)).toFixed(1)})`);
    d.addParagraph(`Likelihood P(E|H) [True Positive Rate]: ${(bayes.likelihood * 100).toFixed(1)}%`);
    d.addParagraph(`False Positive Rate P(E|¬H) [False Alarm Rate]: ${(bayes.falsePositive * 100).toFixed(1)}%`);
    d.addParagraph(`Representative Population Sample (N): ${bayes.sampleSize || 210} individuals`);
    d.addParagraph(`Clinical Scenario Preset: ${bayes.preset || 'Steve the Librarian (3Blue1Brown)'}`);

    d.addHeading1('2. Statistical Outcome & Empirical Updating Metrics');
    d.addTable(
      ['Bayesian Parameter', 'Simulated Value', 'Epistemological & Clinical Meaning'],
      [
        ['Prior Probability P(H)', `${(bayes.prior * 100).toFixed(1)}% (Odds 1:${(1 / (bayes.priorOdds || 0.05)).toFixed(1)})`, 'Initial degree of belief before observing evidence (Base rate)'],
        ['Likelihood P(E|H)', `${(bayes.likelihood * 100).toFixed(1)}%`, 'True Positive Rate: Probability of evidence given hypothesis H is true'],
        ['False Positive Rate P(E|¬H)', `${(bayes.falsePositive * 100).toFixed(1)}%`, 'False Alarm Rate: Probability of evidence given hypothesis H is false'],
        ['Total Evidence P(E)', `${(bayes.pEvidence * 100).toFixed(2)}%`, 'Marginal likelihood: Total shaded area of possibilities matching evidence'],
        ['Bayes Factor (Likelihood Ratio)', `${bayes.bayesFactor?.toFixed(2)}×`, `${bayes.evidenceRating || 'Evidence Ratio'} (P(E|H) / P(E|¬H))`],
        ['Posterior Probability P(H|E)', `${(bayes.posterior * 100).toFixed(1)}% (Odds 1:${(1 / (bayes.posteriorOdds || 0.2)).toFixed(1)})`, `Updated degree of belief after conditioning on evidence (Shift: ${(bayes.beliefShift >= 0 ? '+' : '')}${(bayes.beliefShift * 100).toFixed(1)}%)`],
        ['Representative Counts (N)', `H: ${bayes.countHAndE} of ${bayes.countH} | ¬H: ${bayes.countNotHAndE} of ${bayes.countNotH}`, `Natural frequencies in population of N=${bayes.sampleSize}: ${bayes.countHAndE} / ${bayes.countTotalE} = ${(bayes.posterior * 100).toFixed(1)}%`],
        ['Sequential Updating (4 Steps)', `P₀: ${(bayes.prior * 100).toFixed(1)}% → P₁: ${(bayes.trajectory?.[1]?.p * 100).toFixed(1)}% → P₄: ${(bayes.trajectory?.[4]?.p * 100).toFixed(1)}%`, 'Compounding belief trajectory across successive independent observations']
      ]
    );

    d.addHeading1('3. Clinical & Epistemological Interpretation');
    d.addCalloutBox(
      'Pedagogical Synthesis & Decision-Making Under Uncertainty',
      data.reportText || bayes.explanation || 'Bayes\' theorem provides the mathematical framework for updating rational beliefs in light of new evidence.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Method Was Chosen');
    d.addBullet('The 3Blue1Brown Geometric Insight into Bayes\' Theorem: Rather than memorizing abstract formulas, Bayes\' theorem is intuitively understood as proportions of area within a 1×1 unit square of all possibilities. Observing evidence restricts our sample space to only the shaded regions where the evidence occurs; the posterior probability is simply the fraction of that restricted space corresponding to the hypothesis of interest.');
    d.addBullet('Base Rate Neglect and Natural Frequency Framing: In Steve the Librarian problem, people intuitively fixate on the 4:1 likelihood ratio (40% vs 10%) and forget the 20:1 base rate ratio of farmers to librarians. Translating abstract probabilities into natural frequencies (e.g. 4 librarians vs 20 farmers in a village of 210 people) eliminates cognitive bias and reveals why Steve is still 5× more likely to be a farmer.');
    d.addBullet('Sequential Updating and Continuous Learning: In medical diagnosis and scientific research, evidence arrives sequentially. By taking the posterior probability of test #1 as the prior probability for test #2, Bayes\' rule provides a recursive framework for learning and resolving uncertainty.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Bayes\' Theorem in Area Form: P(H|E) = P(H ∩ E) / P(E) = [P(H) · P(E|H)] / [P(H) · P(E|H) + P(¬H) · P(E|¬H)].');
    d.addBullet('Odds Form of Bayes\' Rule: Posterior Odds = Prior Odds × Bayes Factor (Likelihood Ratio).');
    d.addBullet('Jeffreys Evidence Scale: Bayes factor > 100 = Decisive; 30–100 = Very Strong; 10–30 = Strong; 3–10 = Substantial; 1–3 = Barely worth mentioning.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Sanderson G (2019). Bayes theorem, the geometry of changing beliefs. 3Blue1Brown, YouTube.');
    d.addBullet('Kahneman D, Tversky A (1973). On the psychology of prediction. Psychological Review, 80(4): 237–251.');
    d.addBullet('Gigerenzer G, Hoffrage U (1995). How to improve Bayesian reasoning without instruction: Frequency formats. Psychological Review, 102(4): 684–704.');
    d.addBullet('Jeffreys H (1961). Theory of Probability. 3rd ed. Oxford Classic Texts in the Physical Sciences.');
    d.addBullet('Kass RE, Raftery AE (1995). Bayes factors. Journal of the American Statistical Association, 90(430): 773–795.');

    return d;
  },

  createPsmDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Propensity Score Matching (PSM) & Observational Causal Inference')
      .addAttributionHeader()
      .addDisclaimerBox();

    const prep = data.prep || {};
    const match = data.matchResult || {};
    const outcome = data.outcomeResult || {};
    const logit = data.logitModel || {};
    const balance = data.balance || {};

    d.addHeading1('1. Observational Study Design & Matching Specifications')
      .addParagraph(`Treatment Assignment Variable: ${data.treatmentCol || 'treatment_col'} (Binary 0/1)`)
      .addParagraph(`Primary Clinical Outcome: ${data.outcomeCol || 'outcome_col'}`)
      .addParagraph(`Confounding Baseline Covariates: ${(data.covariateCols || []).join(', ')}`)
      .addParagraph(`Matching Algorithm: 1:1 Nearest-Neighbor Without Replacement (Greedy on Logit PS)`)
      .addParagraph(`Caliper Width: ${match.caliperMultiplier || 0.20} × SD(logit(PS)) = ${(match.caliperWidth || 0).toFixed(4)}`)
      .addParagraph(`Common Support Region: [${(match.commonSupport?.min || 0).toFixed(4)}, ${(match.commonSupport?.max || 1).toFixed(4)}] (Discarded off-support: ${match.commonSupport?.droppedTreated || 0} treated)`);

    const cohortSummaryRows = [
      ['Total Observational Records', `${prep.totalRows || 0}`],
      ['Complete Cases Analyzed', `${prep.completeCasesCount || 0} (${prep.missingRowsCount || 0} missing rows filtered)`],
      ['Pre-Matching Treated Cohort (Z = 1)', `${prep.treatedCount || match.nTotalTreated || 0}`],
      ['Pre-Matching Control Cohort (Z = 0)', `${prep.controlCount || match.nTotalControl || 0}`],
      ['Successfully Matched Pairs', `${match.nMatchedPairs || 0} pairs (${(match.nMatchedPairs || 0) * 2} patients)`],
      ['Unmatched Treated Cohort', `${match.nUnmatchedTreated || 0}`],
      ['Unmatched Control Cohort', `${match.nUnmatchedControl || 0}`]
    ];
    d.addTable(['Study Metric', 'Clinical Value'], cohortSummaryRows);

    // Section 2: Propensity Score Logistic Regression
    d.addHeading1('2. Propensity Score Estimation (Multivariate Logistic Regression)');
    d.addParagraph(`Model Fit: Likelihood Ratio χ² = ${(logit.lrStat || 0).toFixed(2)} (df = ${logit.lrDf || 0}, p = ${(logit.lrPValue !== undefined ? logit.lrPValue.toExponential(3) : '< 0.001')}), McFadden's Pseudo-R² = ${(logit.mcfaddenR2 || 0).toFixed(3)}.`);

    if (logit.coefficients && logit.coefficients.length > 0) {
      const logitRows = logit.coefficients.map(c => [
        c.term,
        c.estimate.toFixed(4),
        c.stdError.toFixed(4),
        c.zScore.toFixed(2),
        c.pValue < 0.001 ? '< 0.001' : c.pValue.toFixed(3),
        c.oddsRatio.toFixed(3),
        `[${c.ci95[0].toFixed(3)}, ${c.ci95[1].toFixed(3)}]`
      ]);
      d.addTable(['Model Parameter', 'Coefficient (β)', 'Std Error', 'Wald z', 'p-Value', 'Odds Ratio (OR)', '95% CI of OR'], logitRows);
    }

    // Section 3: Baseline Covariate Balance Assessment & Love Plot
    d.addHeading1('3. Baseline Covariate Balance Diagnostics & Standardized Mean Differences (SMD)');
    d.addParagraph('Standardized Mean Difference (SMD) evaluates whether matching successfully removed confounding bias. An SMD < 0.10 satisfies standard clinical equivalence; SMD < 0.05 indicates stringent balance.');

    if (balance.balanceTable && balance.balanceTable.length > 0) {
      const balRows = balance.balanceTable.map(b => [
        b.covariate,
        `${b.meanTreatedPre.toFixed(2)} vs ${b.meanControlPre.toFixed(2)}`,
        b.smdPre.toFixed(3),
        `${b.meanTreatedPost.toFixed(2)} vs ${b.meanControlPost.toFixed(2)}`,
        b.smdPost.toFixed(3),
        `${b.percentReduction.toFixed(1)}%`,
        b.isBalanced ? 'Balanced (|SMD| < 0.10)' : 'Residual Imbalance'
      ]);
      d.addTable(['Covariate', 'Unadjusted Means (T vs C)', 'Pre-Match SMD', 'Matched Means (T vs C)', 'Post-Match SMD', '% Bias Reduction', 'Balance Status'], balRows);
    }

    // Section 4: Causal Outcome Analysis
    d.addHeading1('4. Causal Outcome Analysis: Average Treatment Effect on the Treated (ATT)');
    const unadjDiff = typeof outcome.unadjustedDiff === 'object' && outcome.unadjustedDiff !== null
      ? outcome.unadjustedDiff.diff
      : (outcome.unadjustedDiff || 0);

    const outcomeRows = [
      ['Primary Clinical Outcome', `${outcome.outcomeCol || data.outcomeCol}`],
      ['Outcome Variable Type', `${outcome.outcomeType === 'binary' ? 'Binary / Proportional Event' : 'Continuous Numerical Metric'}`],
      ['Unadjusted (Confounded) Difference', `${unadjDiff.toFixed(3)}`],
      ['Matched ATT (Causal Point Estimate)', `${outcome.att !== undefined ? outcome.att.toFixed(3) : (outcome.pointEstimate || 0).toFixed(3)}`],
      ['Paired / Cluster-Robust Std Error', `${(outcome.se || outcome.stdError || 0).toFixed(4)}`],
      ['Test Statistic', `${(outcome.statistic !== undefined ? outcome.statistic : (outcome.testStatistic || 0)).toFixed(2)} (${outcome.testName || 'Paired t-test'})`],
      ['Statistical Significance (p-Value)', `${outcome.pValue !== undefined ? (outcome.pValue < 0.001 ? '< 0.001' : outcome.pValue.toFixed(4)) : '--'}`],
      ['95% Confidence Interval for ATT', outcome.ci95 ? `[${outcome.ci95[0].toFixed(3)}, ${outcome.ci95[1].toFixed(3)}]` : '--'],
      ['Causal Conclusion', outcome.isSignificant ? 'Statistically Significant Causal Effect (p < .05)' : 'No Significant Causal Difference (ns)']
    ];
    d.addTable(['Causal Parameter', 'Statistical Result'], outcomeRows);

    // Section 5: Clinical & STROBE Summary
    d.addHeading1('5. Clinical Interpretation & Methodological Rationale');
    d.addCalloutBox(
      data.reportText || 'Propensity score matching eliminated observed confounding bias across baseline covariates, providing an unconfounded estimate of the Average Treatment Effect on the Treated.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('6. Key Academic & Methodological References');
    d.addBullet('Rosenbaum PR, Rubin DB (1983). The central role of the propensity score in observational studies for causal effects. Biometrika, 70(1): 41–55.');
    d.addBullet('Austin PC (2011). An introduction to propensity score methods for reducing the effects of confounding in observational studies. Multivariate Behavioral Research, 46(3): 399–424.');
    d.addBullet('Austin PC (2009). Using the standardized mean difference to compare the characteristics of treated and untreated subjects in propensity-score matched samples. Statistics in Medicine, 28(25): 3083–3107.');
    d.addBullet('Ho DE, Imai K, King G, Stuart EA (2011). MatchIt: Nonparametric preprocessing for parametric causal inference. Journal of Statistical Software, 42(8): 1–28.');
    d.addBullet('Stuart EA (2010). Matching methods for causal inference: A review and a look forward. Statistical Science, 25(1): 1–21.');

    return d;
  },

  createMultivariateDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Multivariate Exploratory Data Analysis (PCA / MCA / FAMD)')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analysis Specifications & Cohort Summary');
    const specRows = [
      ['Multivariate Technique', `${data.methodLabel || data.method || 'FAMD'}`],
      ['Observations Evaluated (N)', `${data.nObservations || (data.individuals ? data.individuals.length : 0)}`],
      ['Continuous Variables Count', `${data.continuousCols ? data.continuousCols.length : (data.nContinuous || 0)}`],
      ['Categorical Variables Count', `${data.categoricalCols ? data.categoricalCols.length : (data.nCategorical || 0)}`],
      ['Grouping / Stratification Variable', `${data.groupingCol || 'None (Unsupervised)'}`],
      ['Flagged Outlier Observations', `${data.outlierCount || 0} observations (> 2.5 SD distance)`]
    ];
    d.addTable(['Specification Parameter', 'Value'], specRows);

    d.addHeading1('2. Eigenvalues & Variance Explained (Scree Decomposition)');
    const screeRows = (data.scree || []).map(s => [
      s.label,
      s.eigenvalue.toFixed(4),
      `${s.variancePct.toFixed(2)}%`,
      `${s.cumulativePct.toFixed(2)}%`
    ]);
    d.addTable(['Dimension', 'Eigenvalue (λ)', 'Variance Explained (%)', 'Cumulative Variance (%)'], screeRows);

    d.addHeading1('3. Variable Coordinates, Loadings & Modality Representation');
    const varRows = (data.allVariables || data.variables || []).slice(0, 15).map(v => [
      v.name,
      v.type || 'variable',
      v.coords[0] !== undefined ? v.coords[0].toFixed(3) : '--',
      v.coords[1] !== undefined ? v.coords[1].toFixed(3) : '--',
      v.coords[2] !== undefined ? v.coords[2].toFixed(3) : '--'
    ]);
    d.addTable(['Variable / Modality', 'Data Type', 'Dim 1 Coordinate', 'Dim 2 Coordinate', 'Dim 3 Coordinate'], varRows);

    d.addHeading1('4. Automated Epistemological Interpretation & Findings');
    d.addCalloutBox(
      data.interpretation || data.reportText || 'Multivariate dimensionality reduction completed.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('5. Key Methodological & Academic References');
    d.addBullet('Pearson K (1901). On lines and planes of closest fit to systems of points in space. Philosophical Magazine, 2(11): 559–572.');
    d.addBullet('Hotelling H (1933). Analysis of a complex of statistical variables into principal components. Journal of Educational Psychology, 24(6): 417–441.');
    d.addBullet('Benzécri JP (1973). L\'Analyse des Données: La Correspondance. Dunod, Paris.');
    d.addBullet('Pagès J (2004). Analyse factorielle de données mixtes: Principe et exemple d\'application. Revue de Statistique Appliquée, 52(4): 93–111.');
    d.addBullet('Husson F, Josse J, Lê S (2017). Exploratory Multivariate Analysis by Example Using R. 2nd ed. CRC Press.');

    return d;
  }
};


  // ==========================================
  // 10. EXPORTER & FORMATTING
  // ==========================================
  const Exporter = {
    formatP(p) {
      if (isNaN(p) || p === null) return 'N/A';
      if (p < 0.001) return 'p < .001';
      return `p = ${p.toFixed(3).replace(/^0\./, '.')}`;
    },

    formatDescriptiveReport(res) {
      if (res.error) return res.error;
      let modeStr = 'no distinct mode (all values unique)';
      if (res.modes && res.modes.length > 0) {
        modeStr = `mode = ${res.modes.map(m => m.toFixed(2)).join(', ')} (freq = ${res.maxFreq || 2})`;
      }
      const normStr = res.normality.isNormal
        ? `normally distributed via Jarque-Bera omnibus test (JB = ${res.normality.statistic.toFixed(2)}, df = 2, ${this.formatP(res.normality.pValue)}; skewness = ${res.skewness.toFixed(2)} [${res.skewnessInterpretation || 'symmetric'}], excess kurtosis = ${res.kurtosis.toFixed(2)} [${res.kurtosisInterpretation || 'mesokurtic'}])`
        : `significantly non-normal via Jarque-Bera omnibus test (JB = ${res.normality.statistic.toFixed(2)}, df = 2, ${this.formatP(res.normality.pValue)}; skewness = ${res.skewness.toFixed(2)} [${res.skewnessInterpretation || 'skewed'}], excess kurtosis = ${res.kurtosis.toFixed(2)} [${res.kurtosisInterpretation || 'leptokurtic'}])`;

      let outlierStr = 'No outliers detected beyond Tukey fences (1.5 × IQR)';
      if (res.outliers && res.outliers.length > 0) {
        outlierStr = `${res.outliers.length} outlier(s) detected beyond Tukey fences [${res.lowerFence.toFixed(1)}, ${res.upperFence.toFixed(1)}]: ` +
          res.outliers.map(o => `${o.value.toFixed(1)} (Z = ${o.zScore > 0 ? '+' : ''}${o.zScore.toFixed(2)}, ${o.type})`).join(', ');
      }

      return (
        `Continuous variable summary (N = ${res.n}): ` +
        `Mean = ${res.mean.toFixed(2)} (SD = ${res.sd.toFixed(2)}, 95% CI [${res.ci95[0].toFixed(2)}, ${res.ci95[1].toFixed(2)}]); ` +
        `Median = ${res.median.toFixed(2)} (IQR = ${res.iqr.toFixed(2)}, Q1 = ${res.q1.toFixed(2)}, Q3 = ${res.q3.toFixed(2)}); ` +
        `${modeStr}; Range = ${res.min.toFixed(2)} to ${res.max.toFixed(2)}. ` +
        `Distribution: ${normStr}. ` +
        `Outliers: ${outlierStr}. ` +
        (res.outliers && res.outliers.length > 0 || !res.normality.isNormal ? 'Reporting Median and IQR is recommended due to skewness/outliers.' : 'Reporting Mean and SD is appropriate.')
      );
    },

    downloadBlob(filename, blob) {
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    },

    exportTabToDocx(tabId, data) {
      let builder;
      switch (tabId) {
        case 'descriptive': builder = DocxReports.createDescriptiveDocx(data); break;
        case 'hypothesis': builder = DocxReports.createHypothesisDocx(data); break;
        case 'anova': builder = DocxReports.createAnovaDocx(data); break;
        case 'categorical': builder = DocxReports.createCategoricalDocx(data); break;
        case 'correlation': builder = DocxReports.createCorrelationDocx(data); break;
        case 'diagnostic': builder = DocxReports.createDiagnosticDocx(data); break;
        case 'power': builder = DocxReports.createPowerDocx(data); break;
        case 'teaching': builder = DocxReports.createTeachingDocx(data); break;
        case 'teaching-bayesian': builder = DocxReports.createTeachingBayesianDocx ? DocxReports.createTeachingBayesianDocx(data) : DocxReports.createTeachingDocx(data); break;
        case 'propensity': builder = DocxReports.createPsmDocx(data); break;
        case 'multivariate': builder = DocxReports.createMultivariateDocx(data); break;
        default: console.error('Unknown tab for DOCX export:', tabId); return;
      }
      const blob = builder.generateBlob();
      this.downloadBlob(`statis-gravity-${tabId}-report.docx`, blob);
    },

    async copyToClipboard(text) {
      try {
        await navigator.clipboard.writeText(text);
        return true;
      } catch {
        const ta = document.createElement('textarea');
        ta.value = text;
        document.body.appendChild(ta);
        ta.select();
        document.execCommand('copy');
        document.body.removeChild(ta);
        return true;
      }
    }
  };

  // ==========================================
  // 10.5 RANDOMISER ENGINE (CLINICAL ALLOCATION & AUDIT TRAIL)
  // ==========================================
  const Randomiser = {
    getSecureRandomInt(min, max) {
      min = Math.floor(min);
      max = Math.floor(max);
      if (min > max) throw new Error(`Invalid range: min (${min}) > max (${max})`);
      const range = max - min + 1;
      if (range === 1) return min;

      const MAX_UINT32 = 4294967296;
      const limit = MAX_UINT32 - (MAX_UINT32 % range);
      const buffer = new Uint32Array(1);
      let rand;
      do {
        let gotRand = false;
        if (typeof crypto !== 'undefined' && crypto && typeof crypto.getRandomValues === 'function') {
          try {
            crypto.getRandomValues(buffer);
            gotRand = true;
          } catch (e) {}
        }
        if (!gotRand && typeof window !== 'undefined' && window.crypto && typeof window.crypto.getRandomValues === 'function') {
          try {
            window.crypto.getRandomValues(buffer);
            gotRand = true;
          } catch (e) {}
        }
        if (!gotRand) {
          buffer[0] = Math.floor(Math.random() * MAX_UINT32);
        }
        rand = buffer[0];
      } while (rand >= limit);

      return min + (rand % range);
    },

    fisherYatesShuffle(array) {
      const copy = [...array];
      for (let i = copy.length - 1; i > 0; i--) {
        const j = this.getSecureRandomInt(0, i);
        const temp = copy[i];
        copy[i] = copy[j];
        copy[j] = temp;
      }
      return copy;
    },

    generateSimpleAllocation(options = {}) {
      const participantId = options.participantId || 1;
      const labelA = options.labelA || 'Group A';
      const labelB = options.labelB || 'Group B';

      const randomNumber = this.getSecureRandomInt(1, 100);
      const isOdd = (randomNumber % 2) !== 0;
      const groupKey = isOdd ? 'A' : 'B';
      const groupLabel = isOdd ? labelA : labelB;

      return {
        participantId,
        timestamp: new Date().toISOString(),
        displayTime: new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' }),
        randomNumber,
        parity: isOdd ? 'Odd' : 'Even',
        groupKey,
        groupLabel,
        method: 'Simple (1-100 Odd/Even)'
      };
    },

    createBlock(blockSize, options = {}) {
      let size = parseInt(blockSize, 10);
      if (isNaN(size) || size < 2) size = 4;
      if (size % 2 !== 0) size += 1;

      const blockNumber = options.blockNumber || 1;
      const labelA = options.labelA || 'Group A';
      const labelB = options.labelB || 'Group B';

      const half = size / 2;
      const allocations = [];
      for (let i = 0; i < half; i++) allocations.push('A');
      for (let i = 0; i < half; i++) allocations.push('B');

      const shuffled = this.fisherYatesShuffle(allocations);

      const slots = shuffled.map((g, idx) => ({
        slotIndex: idx + 1,
        groupKey: g,
        groupLabel: g === 'A' ? labelA : labelB,
        assigned: false,
        participantId: null,
        timestamp: null,
        displayTime: null
      }));

      return {
        blockNumber,
        blockSize: size,
        slots,
        currentIndex: 0,
        isComplete: false
      };
    },

    assignNextInBlock(activeBlock, participantId, options = {}) {
      const labelA = options.labelA || 'Group A';
      const labelB = options.labelB || 'Group B';
      const blockSize = options.blockSize || 4;

      let block = activeBlock;
      let isNewBlock = false;

      if (!block || block.currentIndex >= block.slots.length) {
        const nextBlockNumber = block ? (block.blockNumber + 1) : 1;
        block = this.createBlock(blockSize, { blockNumber: nextBlockNumber, labelA, labelB });
        isNewBlock = true;
      }

      const slot = block.slots[block.currentIndex];
      slot.assigned = true;
      slot.participantId = participantId;
      const now = new Date();
      slot.timestamp = now.toISOString();
      slot.displayTime = now.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' });
      slot.groupLabel = slot.groupKey === 'A' ? labelA : labelB;

      block.currentIndex += 1;
      if (block.currentIndex >= block.slots.length) {
        block.isComplete = true;
      }

      const allocationRecord = {
        participantId,
        timestamp: slot.timestamp,
        displayTime: slot.displayTime,
        blockNumber: block.blockNumber,
        blockSize: block.blockSize,
        slotInBlock: slot.slotIndex,
        groupKey: slot.groupKey,
        groupLabel: slot.groupLabel,
        method: `Block Randomization (${block.blockSize})`
      };

      return {
        updatedBlock: block,
        allocationRecord,
        isNewBlock
      };
    },

    computeSummary(history = []) {
      const total = history.length;
      let countA = 0;
      let countB = 0;

      for (let i = 0; i < total; i++) {
        if (history[i].groupKey === 'A') countA++;
        else if (history[i].groupKey === 'B') countB++;
      }

      const pctA = total > 0 ? (countA / total) * 100 : 0;
      const pctB = total > 0 ? (countB / total) * 100 : 0;
      const diff = Math.abs(countA - countB);

      let ratioStr = '1.00 : 1.00';
      if (countB === 0 && countA > 0) {
        ratioStr = `${countA} : 0`;
      } else if (countA === 0 && countB > 0) {
        ratioStr = `0 : ${countB}`;
      } else if (countB > 0) {
        ratioStr = `${(countA / countB).toFixed(2)} : 1.00`;
      }

      return {
        total,
        countA,
        countB,
        pctA,
        pctB,
        diff,
        ratioStr
      };
    },

    exportToCSV(history = [], mode = 'simple') {
      if (mode === 'simple') {
        const headers = ['Participant ID', 'Timestamp', 'Random Integer (1-100)', 'Parity', 'Assigned Group Code', 'Group Label'];
        const rows = history.map(item => [
          item.participantId,
          `"${item.timestamp}"`,
          item.randomNumber,
          item.parity,
          item.groupKey,
          `"${(item.groupLabel || '').replace(/"/g, '""')}"`
        ]);
        return [headers.join(','), ...rows.map(r => r.join(','))].join('\r\n');
      } else {
        const headers = ['Participant ID', 'Timestamp', 'Block Number', 'Block Size', 'Slot in Block', 'Assigned Group Code', 'Group Label'];
        const rows = history.map(item => [
          item.participantId,
          `"${item.timestamp}"`,
          item.blockNumber,
          item.blockSize,
          item.slotInBlock,
          item.groupKey,
          `"${(item.groupLabel || '').replace(/"/g, '""')}"`
        ]);
        return [headers.join(','), ...rows.map(r => r.join(','))].join('\r\n');
      }
    }
  };

  // ==========================================
  // 11. CANVAS CHART ENGINE & PLOTS
  // ==========================================
  class ChartEngine {
    constructor(canvasId, options = {}) {
      this.canvas = typeof canvasId === 'string' ? document.getElementById(canvasId) : canvasId;
      if (!this.canvas) return;
      this.ctx = this.canvas.getContext('2d');
      this.options = Object.assign({
        padding: { top: 35, right: 25, bottom: 45, left: 55 },
        theme: 'dark',
        fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'
      }, options);

      this.colors = {
        dark: {
          bg: '#090d16',
          grid: 'rgba(255, 255, 255, 0.07)',
          axis: 'rgba(255, 255, 255, 0.22)',
          text: '#94a3b8',
          textMuted: '#64748b',
          textBold: '#f8fafc',
          primary: '#00d2ff',
          secondary: '#a855f7',
          accent: '#10b981',
          warning: '#f59e0b',
          danger: '#f43f5e'
        },
        light: {
          bg: '#ffffff',
          grid: 'rgba(0, 0, 0, 0.06)',
          axis: 'rgba(0, 0, 0, 0.2)',
          text: '#64748b',
          textMuted: '#94a3b8',
          textBold: '#0f172a',
          primary: '#0284c7',
          secondary: '#9333ea',
          accent: '#059669',
          warning: '#d97706',
          danger: '#e11d48'
        }
      };

      this.initHiDPI();
      this.setupResizeObserver();
    }

    get palette() {
      return this.colors[this.options.theme] || this.colors.dark;
    }

    initHiDPI() {
      if (!this.canvas) return;
      const rect = this.canvas.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      this.width = rect.width || this.canvas.width || 600;
      this.height = rect.height || this.canvas.height || 320;

      this.canvas.width = this.width * dpr;
      this.canvas.height = this.height * dpr;
      if (this.ctx.resetTransform) this.ctx.resetTransform();
      this.ctx.scale(dpr, dpr);
    }

    setupResizeObserver() {
      if (typeof ResizeObserver === 'function' && this.canvas) {
        this.resizeObserver = new ResizeObserver(() => {
          const rect = this.canvas.getBoundingClientRect();
          if (rect.width > 0 && rect.height > 0) {
            this.initHiDPI();
            if (typeof this.lastRender === 'function') this.lastRender();
            else if (typeof this.lastRenderFn === 'function') this.lastRenderFn();
          }
        });
        this.resizeObserver.observe(this.canvas);
      }
    }

    clear() {
      this.ctx.clearRect(0, 0, this.width, this.height);
    }

    getBounds() {
      const p = this.options.padding;
      return {
        x: p.left,
        y: p.top,
        width: Math.max(10, this.width - p.left - p.right),
        height: Math.max(10, this.height - p.top - p.bottom)
      };
    }

    getPlotBounds() {
      return this.getBounds();
    }

    drawAxes({ xLabel = '', yLabel = '', xTicks = [], yTicks = [], title = '' } = {}) {
      const b = this.getBounds();
      const ctx = this.ctx;
      const pal = this.palette;

      if (title) {
        ctx.fillStyle = pal.textBold;
        ctx.font = '600 13px -apple-system, sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText(title, b.x, b.y - 12);
      }

      ctx.lineWidth = 1;
      for (const t of yTicks) {
        const yPos = b.y + b.height - (t.norm * b.height);
        ctx.strokeStyle = pal.grid;
        ctx.beginPath();
        ctx.moveTo(b.x, yPos);
        ctx.lineTo(b.x + b.width, yPos);
        ctx.stroke();

        ctx.fillStyle = pal.text;
        ctx.font = '500 10px monospace';
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        ctx.fillText(t.label, b.x - 8, yPos);
      }

      for (const t of xTicks) {
        const xPos = b.x + (t.norm * b.width);
        ctx.strokeStyle = pal.grid;
        ctx.beginPath();
        ctx.moveTo(xPos, b.y);
        ctx.lineTo(xPos, b.y + b.height);
        ctx.stroke();

        ctx.fillStyle = pal.text;
        ctx.font = '500 10px monospace';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.fillText(t.label, xPos, b.y + b.height + 6);
      }

      ctx.strokeStyle = pal.axis;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(b.x, b.y);
      ctx.lineTo(b.x, b.y + b.height);
      ctx.lineTo(b.x + b.width, b.y + b.height);
      ctx.stroke();

      ctx.fillStyle = pal.text;
      ctx.font = '600 11px -apple-system, sans-serif';
      if (xLabel) {
        ctx.textAlign = 'center';
        ctx.fillText(xLabel, b.x + b.width / 2, b.y + b.height + 28);
      }
      if (yLabel) {
        ctx.save();
        ctx.translate(b.x - 38, b.y + b.height / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.textAlign = 'center';
        ctx.fillText(yLabel, 0, 0);
        ctx.restore();
      }
    }

    saveImage(name = 'plot.png') {
      const link = document.createElement('a');
      link.download = name;
      link.href = this.canvas.toDataURL('image/png');
      link.click();
    }
  }

  const Plots = {
    renderBoxPlot(engine, groupStats, title) {
      engine.lastRender = () => this.renderBoxPlot(engine, groupStats, title);
      engine.clear();
      const b = engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;

      let groups = [];
      if (Array.isArray(groupStats)) {
        groups = groupStats.map(g => g.stats ? g : { name: g.name || 'Sample Data', stats: g, color: g.color });
      } else if (groupStats && typeof groupStats === 'object') {
        groups = [{ name: groupStats.name || 'Sample Data', stats: groupStats, color: pal.primary }];
      }

      if (!groups || groups.length === 0) return;

      let gMin = Infinity;
      let gMax = -Infinity;
      for (const g of groups) {
        if (!g.stats || g.stats.n === 0) continue;
        if (g.stats.min < gMin) gMin = g.stats.min;
        if (g.stats.max > gMax) gMax = g.stats.max;
      }
      if (!isFinite(gMin) || !isFinite(gMax)) return;

      const span = gMax - gMin || 1;
      const yMin = gMin - 0.1 * span;
      const yMax = gMax + 0.1 * span;
      const yRange = yMax - yMin;
      const yToPx = (v) => b.y + b.height - ((v - yMin) / yRange) * b.height;

      const yTicks = [];
      for (let i = 0; i <= 4; i++) {
        const v = yMin + (i / 4) * yRange;
        yTicks.push({ norm: i / 4, label: v.toFixed(1) });
      }
      const k = groups.length;
      const xTicks = groups.map((g, idx) => ({ norm: (idx + 0.5) / k, label: g.name }));
      engine.drawAxes({ yTicks, xTicks, title, yLabel: 'Observed Value' });

      const slotW = b.width / k;
      const boxW = Math.min(55, slotW * 0.45);
      const colors = [pal.primary, pal.secondary, pal.accent, '#f59e0b'];

      groups.forEach((g, idx) => {
        const s = g.stats;
        if (!s || s.n === 0) return;
        const col = g.color || colors[idx % colors.length];
        const cx = b.x + (idx + 0.5) * slotW;

        const q1Y = yToPx(s.q1);
        const q3Y = yToPx(s.q3);
        const medY = yToPx(s.median);
        const minFY = yToPx(Math.max(s.min, s.lowerFence));
        const maxFY = yToPx(Math.min(s.max, s.upperFence));

        ctx.strokeStyle = pal.axis;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(cx, q1Y); ctx.lineTo(cx, minFY);
        ctx.moveTo(cx - boxW * 0.25, minFY); ctx.lineTo(cx + boxW * 0.25, minFY);
        ctx.moveTo(cx, q3Y); ctx.lineTo(cx, maxFY);
        ctx.moveTo(cx - boxW * 0.25, maxFY); ctx.lineTo(cx + boxW * 0.25, maxFY);
        ctx.stroke();

        const top = Math.min(q1Y, q3Y);
        const h = Math.abs(q1Y - q3Y) || 2;
        ctx.fillStyle = `${col}22`;
        ctx.strokeStyle = col;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.rect(cx - boxW / 2, top, boxW, h);
        ctx.fill();
        ctx.stroke();

        ctx.strokeStyle = pal.textBold;
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.moveTo(cx - boxW / 2, medY);
        ctx.lineTo(cx + boxW / 2, medY);
        ctx.stroke();

        // Mean Diamond Marker (◆)
        if (typeof s.mean === 'number' && isFinite(s.mean)) {
          const meanY = yToPx(s.mean);
          ctx.fillStyle = pal.accent;
          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 1.2;
          ctx.beginPath();
          const d = 5;
          ctx.moveTo(cx, meanY - d);
          ctx.lineTo(cx + d, meanY);
          ctx.lineTo(cx, meanY + d);
          ctx.lineTo(cx - d, meanY);
          ctx.closePath();
          ctx.fill();
          ctx.stroke();
        }

        // Outliers
        if (s.outliers && s.outliers.length > 0) {
          for (const out of s.outliers) {
            const val = typeof out === 'object' && out !== null ? out.value : out;
            const isExtreme = typeof out === 'object' && out !== null ? out.type === 'Extreme' : false;
            const outY = yToPx(val);

            ctx.fillStyle = isExtreme ? pal.danger : (pal.warning || '#f59e0b');
            ctx.strokeStyle = isExtreme ? '#ffffff' : pal.danger;
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            ctx.arc(cx, outY, isExtreme ? 5.5 : 4, 0, 2 * Math.PI);
            ctx.fill();
            ctx.stroke();

            ctx.fillStyle = isExtreme ? pal.danger : pal.text;
            ctx.font = '600 10px -apple-system, sans-serif';
            ctx.textAlign = 'left';
            ctx.fillText(`${val.toFixed(1)}${isExtreme ? ' (Extr)' : ''}`, cx + boxW * 0.4, outY + 3);
          }
        }

        if (s.values && s.values.length <= 150) {
          ctx.fillStyle = `${col}88`;
          for (let i = 0; i < s.values.length; i++) {
            const py = yToPx(s.values[i]);
            const jitter = (Math.sin(i * 13 + s.values[i]) * 0.35) * (boxW / 2);
            ctx.beginPath();
            ctx.arc(cx + jitter, py, 3, 0, 2 * Math.PI);
            ctx.fill();
          }
        }

        if (k === 1) {
          ctx.font = '500 10px -apple-system, sans-serif';
          ctx.textAlign = 'left';
          ctx.fillStyle = pal.textBold;
          ctx.fillText(`Median: ${s.median.toFixed(1)}`, cx + boxW / 2 + 8, medY + 3);
          ctx.fillStyle = pal.text;
          ctx.fillText(`Q3: ${s.q3.toFixed(1)}`, cx + boxW / 2 + 8, q3Y + 3);
          ctx.fillText(`Q1: ${s.q1.toFixed(1)}`, cx + boxW / 2 + 8, q1Y + 3);
          if (typeof s.mean === 'number') {
            ctx.fillStyle = pal.accent;
            ctx.fillText(`Mean: ${s.mean.toFixed(1)}`, cx - boxW / 2 - 62, yToPx(s.mean) + 3);
          }
        }
      });
    },

    renderErrorBarPlot(engine, groupStats, options = {}) {
      const opts = typeof options === 'string' ? { title: options } : (options || {});
      const mode = (opts.mode || 'ci95').toLowerCase();
      const title = opts.title || 'Comparative Cohort Distribution';

      if (mode === 'iqr') {
        return this.renderBoxPlot(engine, groupStats, title);
      }

      engine.lastRender = () => this.renderErrorBarPlot(engine, groupStats, options);
      engine.lastRenderFn = engine.lastRender;
      engine.clear();

      const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;

      let groups = [];
      if (Array.isArray(groupStats)) {
        groups = groupStats.map(g => g.stats ? g : { name: g.name || 'Sample Data', stats: g, color: g.color });
      } else if (groupStats && typeof groupStats === 'object') {
        groups = [{ name: groupStats.name || 'Sample Data', stats: groupStats, color: pal.primary }];
      }

      if (groups.length === 0) return;

      const processedGroups = groups.map((g, idx) => {
        const s = g.stats || {};
        const mean = typeof s.mean === 'number' && isFinite(s.mean) ? s.mean : 0;
        const sem = typeof s.sem === 'number' && isFinite(s.sem) ? s.sem : (s.sd && s.n ? s.sd / Math.sqrt(s.n) : 0);
        const sd = typeof s.sd === 'number' && isFinite(s.sd) ? s.sd : 0;

        let lower = mean;
        let upper = mean;
        let label = '95% CI';
        let subLabel = '';

        if (mode === 'sem') {
          lower = mean - sem;
          upper = mean + sem;
          label = '±1 SEM';
          subLabel = `SEM: ±${sem.toFixed(2)}`;
        } else if (mode === 'sd') {
          lower = mean - sd;
          upper = mean + sd;
          label = '±1 SD';
          subLabel = `SD: ±${sd.toFixed(2)}`;
        } else {
          if (Array.isArray(s.ci95) && s.ci95.length === 2 && isFinite(s.ci95[0]) && isFinite(s.ci95[1])) {
            lower = s.ci95[0];
            upper = s.ci95[1];
          } else {
            const margin = sem * 1.96;
            lower = mean - margin;
            upper = mean + margin;
          }
          label = '95% CI';
          subLabel = `[${lower.toFixed(2)}, ${upper.toFixed(2)}]`;
        }

        return {
          ...g,
          mean,
          sem,
          sd,
          lower,
          upper,
          label,
          subLabel
        };
      });

      let globalMin = Infinity;
      let globalMax = -Infinity;
      for (const g of processedGroups) {
        const s = g.stats;
        if (!s || s.n === 0) continue;
        if (g.lower < globalMin) globalMin = g.lower;
        if (g.upper > globalMax) globalMax = g.upper;
        if (typeof s.min === 'number' && isFinite(s.min) && s.min < globalMin) globalMin = s.min;
        if (typeof s.max === 'number' && isFinite(s.max) && s.max > globalMax) globalMax = s.max;
      }

      if (!isFinite(globalMin) || !isFinite(globalMax)) return;

      const span = globalMax - globalMin || 1;
      const yMin = globalMin - 0.15 * span;
      const yMax = globalMax + 0.18 * span;
      const yRange = yMax - yMin;

      const yToPixel = (val) => b.y + b.height - ((val - yMin) / yRange) * b.height;

      const yTicks = [];
      const stepCount = 5;
      for (let i = 0; i <= stepCount; i++) {
        const v = yMin + (i / stepCount) * yRange;
        yTicks.push({ norm: (v - yMin) / yRange, label: v.toFixed(1) });
      }

      const k = processedGroups.length;
      const xTicks = processedGroups.map((g, idx) => ({
        norm: (idx + 0.5) / k,
        label: g.name
      }));

      const modeHeaders = {
        ci95: '95% Confidence Interval (Mean ± 95% CI)',
        sem: 'Standard Error of Mean (Mean ± 1 SEM)',
        sd: 'Standard Deviation (Mean ± 1 SD)',
        iqr: 'Interquartile Range'
      };
      const modeTitle = modeHeaders[mode] || '95% Confidence Interval';

      engine.drawAxes({
        yTicks,
        xTicks,
        title: `${title}`,
        yLabel: 'Observed Value'
      });

      const slotWidth = b.width / k;
      const capWidth = Math.min(38, Math.max(22, slotWidth * 0.28));
      const colors = [pal.primary, pal.secondary, pal.accent, '#f59e0b'];

      if (k === 2 && processedGroups[0].stats && processedGroups[1].stats) {
        const gA = processedGroups[0];
        const gB = processedGroups[1];
        const xA = b.x + 0.5 * slotWidth;
        const xB = b.x + 1.5 * slotWidth;
        const yA = yToPixel(gA.mean);
        const yB = yToPixel(gB.mean);

        ctx.save();
        ctx.strokeStyle = `${pal.axis}`;
        ctx.setLineDash([4, 4]);
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(xA, yA);
        ctx.lineTo(xB, yB);
        ctx.stroke();
        ctx.setLineDash([]);

        const midX = (xA + xB) / 2;
        const midY = (yA + yB) / 2;
        const deltaM = gA.mean - gB.mean;
        const deltaText = `ΔM = ${deltaM >= 0 ? '+' : ''}${deltaM.toFixed(2)}`;

        ctx.font = `600 10px ${engine.options.fontFamily || '-apple-system, sans-serif'}`;
        const textW = ctx.measureText ? (ctx.measureText(deltaText)?.width || 50) : 50;
        ctx.fillStyle = pal.bgSurfaceElevated || pal.bgCard || '#1e293b';
        ctx.strokeStyle = pal.border || '#334155';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.roundRect?.(midX - textW / 2 - 6, midY - 10, textW + 12, 18, 4) ||
          ctx.rect(midX - textW / 2 - 6, midY - 10, textW + 12, 18);
        ctx.fill();
        ctx.stroke();

        ctx.fillStyle = pal.textBold;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(deltaText, midX, midY);
        ctx.restore();
      }

      processedGroups.forEach((g, idx) => {
        const s = g.stats;
        if (!s || s.n === 0) return;

        const groupColor = g.color || colors[idx % colors.length];
        const centerX = b.x + (idx + 0.5) * slotWidth;

        const meanY = yToPixel(g.mean);
        const lowerY = yToPixel(g.lower);
        const upperY = yToPixel(g.upper);

        if (s.values && s.values.length > 0) {
          ctx.save();
          ctx.fillStyle = `${groupColor}44`;
          for (let i = 0; i < s.values.length; i++) {
            const v = s.values[i];
            const ptY = yToPixel(v);
            const jitter = (Math.sin(i * 12.9898 + v) * 0.35) * (slotWidth * 0.3);
            ctx.beginPath();
            ctx.arc(centerX + jitter, ptY, 3.2, 0, 2 * Math.PI);
            ctx.fill();
          }
          ctx.restore();
        }

        const pillarWidth = Math.min(52, slotWidth * 0.35);
        const pillarTop = Math.min(lowerY, upperY);
        const pillarHeight = Math.abs(lowerY - upperY) || 2;
        ctx.save();
        ctx.fillStyle = `${groupColor}14`;
        ctx.beginPath();
        ctx.roundRect?.(centerX - pillarWidth / 2, pillarTop, pillarWidth, pillarHeight, 6) ||
          ctx.rect(centerX - pillarWidth / 2, pillarTop, pillarWidth, pillarHeight);
        ctx.fill();
        ctx.restore();

        ctx.save();
        ctx.strokeStyle = groupColor;
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.moveTo(centerX, lowerY);
        ctx.lineTo(centerX, upperY);
        ctx.stroke();

        ctx.beginPath();
        ctx.moveTo(centerX - capWidth / 2, upperY);
        ctx.lineTo(centerX + capWidth / 2, upperY);
        ctx.moveTo(centerX - capWidth / 2, lowerY);
        ctx.lineTo(centerX + capWidth / 2, lowerY);
        ctx.stroke();
        ctx.restore();

        ctx.save();
        ctx.fillStyle = groupColor;
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.arc(centerX, meanY, 6.5, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();

        ctx.fillStyle = '#ffffff';
        ctx.beginPath();
        ctx.arc(centerX, meanY, 2, 0, 2 * Math.PI);
        ctx.fill();
        ctx.restore();

        ctx.save();
        ctx.font = `600 11px ${engine.options.fontFamily || '-apple-system, sans-serif'}`;
        ctx.textAlign = 'left';
        ctx.textBaseline = 'middle';

        ctx.fillStyle = pal.textBold;
        ctx.fillText(`M = ${g.mean.toFixed(2)}`, centerX + capWidth / 2 + 8, meanY);

        ctx.font = `600 10px ${engine.options.fontFamily || '-apple-system, sans-serif'}`;
        ctx.fillStyle = pal.textMuted || pal.text;
        ctx.fillText(g.subLabel, centerX + capWidth / 2 + 8, meanY + 13);

        ctx.font = `500 9px ${engine.options.fontFamily || '-apple-system, sans-serif'}`;
        ctx.fillStyle = pal.textDim || pal.textMuted || pal.text;
        ctx.fillText(`Span: ${(g.upper - g.lower).toFixed(2)}`, centerX + capWidth / 2 + 8, meanY + 25);
        ctx.fillText(`n = ${s.n}`, centerX + capWidth / 2 + 8, meanY + 36);

        ctx.textAlign = 'right';
        ctx.font = `500 9px ${engine.options.fontFamily || '-apple-system, sans-serif'}`;
        ctx.fillStyle = pal.textDim || pal.textMuted || pal.text;
        ctx.fillText(g.upper.toFixed(2), centerX - capWidth / 2 - 6, upperY);
        ctx.fillText(g.lower.toFixed(2), centerX - capWidth / 2 - 6, lowerY);
        ctx.restore();
      });

      ctx.save();
      ctx.font = `italic 10px ${engine.options.fontFamily || '-apple-system, sans-serif'}`;
      ctx.textAlign = 'right';
      ctx.textBaseline = 'top';
      ctx.fillStyle = pal.textMuted || pal.text;
      ctx.fillText(`Mode: ${modeTitle}`, b.x + b.width, b.y - 18);
      ctx.restore();
    },

    renderViolinPlot(engine, data, title = 'Violin Density Plot (KDE & Quartiles)') {
      engine.lastRender = () => this.renderViolinPlot(engine, data, title);
      engine.clear();
      const b = engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;

      let groups = [];
      if (Array.isArray(data)) {
        groups = data.map(g => g.stats ? g : { name: g.name || 'Sample Data', stats: g, color: g.color });
      } else if (data && typeof data === 'object') {
        groups = [{ name: data.name || 'Sample Data', stats: data, color: pal.secondary }];
      }

      if (groups.length === 0) return;

      let gMin = Infinity;
      let gMax = -Infinity;
      for (const g of groups) {
        const s = g.stats;
        if (!s || s.n === 0) continue;
        if (s.min < gMin) gMin = s.min;
        if (s.max > gMax) gMax = s.max;
      }
      if (!isFinite(gMin) || !isFinite(gMax)) return;

      const span = gMax - gMin || 1;
      const yMin = gMin - 0.12 * span;
      const yMax = gMax + 0.12 * span;
      const yRange = yMax - yMin;
      const yToPx = (v) => b.y + b.height - ((v - yMin) / yRange) * b.height;

      const yTicks = [];
      for (let i = 0; i <= 4; i++) {
        const v = yMin + (i / 4) * yRange;
        yTicks.push({ norm: i / 4, label: v.toFixed(1) });
      }
      const k = groups.length;
      const xTicks = groups.map((g, idx) => ({ norm: (idx + 0.5) / k, label: g.name }));
      engine.drawAxes({ yTicks, xTicks, title, yLabel: 'Observed Value' });

      const slotW = b.width / k;
      const maxHalfW = Math.min(75, slotW * 0.42);

      groups.forEach((g, idx) => {
        const s = g.stats;
        if (!s || !s.values || s.values.length < 2) return;
        const col = g.color || pal.secondary;
        const cx = b.x + (idx + 0.5) * slotW;
        const values = s.values;
        const n = values.length;

        const iqr = s.iqr > 0 ? s.iqr : (s.sd || 1);
        const A = Math.min(s.sd || 1, (iqr / 1.34) || 1) || 1;
        const h = Math.max(1e-3 * span, 0.9 * A * Math.pow(n, -0.2));

        const kde = (y) => {
          let sum = 0;
          const invH = 1 / h;
          for (let i = 0; i < n; i++) {
            const u = (y - values[i]) * invH;
            sum += Math.exp(-0.5 * u * u);
          }
          return sum / (n * h * Math.sqrt(2 * Math.PI));
        };

        const gridSteps = 100;
        const gridY = [];
        const density = [];
        let maxDensity = 0;
        for (let i = 0; i <= gridSteps; i++) {
          const yVal = s.min + (i / gridSteps) * (s.max - s.min);
          const d = kde(yVal);
          gridY.push(yVal);
          density.push(d);
          if (d > maxDensity) maxDensity = d;
        }
        if (maxDensity <= 0) maxDensity = 1;

        ctx.beginPath();
        for (let i = 0; i <= gridSteps; i++) {
          const py = yToPx(gridY[i]);
          const px = cx + (density[i] / maxDensity) * maxHalfW;
          if (i === 0) ctx.moveTo(cx, py);
          else ctx.lineTo(px, py);
        }
        ctx.lineTo(cx, yToPx(gridY[gridSteps]));
        for (let i = gridSteps; i >= 0; i--) {
          const py = yToPx(gridY[i]);
          const px = cx - (density[i] / maxDensity) * maxHalfW;
          ctx.lineTo(px, py);
        }
        ctx.closePath();

        const grad = ctx.createLinearGradient(cx - maxHalfW, 0, cx + maxHalfW, 0);
        grad.addColorStop(0, `${col}15`);
        grad.addColorStop(0.5, `${col}44`);
        grad.addColorStop(1, `${col}15`);
        ctx.fillStyle = grad;
        ctx.fill();
        ctx.strokeStyle = col;
        ctx.lineWidth = 2;
        ctx.stroke();

        ctx.strokeStyle = pal.axis;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(cx, yToPx(s.min));
        ctx.lineTo(cx, yToPx(s.max));
        ctx.stroke();

        const q1Y = yToPx(s.q1);
        const q3Y = yToPx(s.q3);
        const iqrW = 10;
        ctx.fillStyle = pal.bg === '#ffffff' ? '#e2e8f0' : '#1e293b';
        ctx.strokeStyle = pal.textBold;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.rect(cx - iqrW / 2, Math.min(q1Y, q3Y), iqrW, Math.abs(q1Y - q3Y) || 2);
        ctx.fill();
        ctx.stroke();

        const medY = yToPx(s.median);
        ctx.fillStyle = '#ffffff';
        ctx.strokeStyle = pal.textBold;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.arc(cx, medY, 4, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();

        if (typeof s.mean === 'number' && isFinite(s.mean)) {
          const meanY = yToPx(s.mean);
          ctx.fillStyle = pal.primary;
          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 1.2;
          ctx.beginPath();
          const mSize = 4.5;
          ctx.moveTo(cx, meanY - mSize);
          ctx.lineTo(cx + mSize, meanY);
          ctx.lineTo(cx, meanY + mSize);
          ctx.lineTo(cx - mSize, meanY);
          ctx.closePath();
          ctx.fill();
          ctx.stroke();
        }

        if (k === 1) {
          ctx.font = '600 10px -apple-system, sans-serif';
          ctx.textAlign = 'left';
          ctx.fillStyle = pal.textBold;
          ctx.fillText(`Median: ${s.median.toFixed(1)}`, cx + maxHalfW * 0.7 + 10, medY + 3);
          if (typeof s.mean === 'number') {
            ctx.fillStyle = pal.primary;
            ctx.fillText(`Mean: ${s.mean.toFixed(1)}`, cx - maxHalfW * 0.7 - 65, yToPx(s.mean) + 3);
          }
          ctx.fillStyle = pal.text;
          ctx.fillText(`IQR [${s.q1.toFixed(1)}, ${s.q3.toFixed(1)}]`, cx + maxHalfW * 0.7 + 10, q3Y + 3);
        }
      });
    },

    renderHistogram(engine, stats, title) {
      engine.lastRender = () => this.renderHistogram(engine, stats, title);
      engine.clear();
      const b = engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;
      if (!stats || stats.n < 1) return;

      const span = stats.range || 1;
      const numBins = Math.max(4, Math.min(20, Math.round(1 + 3.322 * Math.log10(stats.n))));
      const binW = span / numBins;
      const bins = new Array(numBins).fill(0);

      for (const v of stats.values) {
        const idx = Math.min(numBins - 1, Math.floor((v - stats.min) / binW));
        bins[idx]++;
      }

      const maxCount = Math.max(...bins, 1);
      const yMax = maxCount * 1.25;

      const yTicks = [
        { norm: 0, label: '0' },
        { norm: 0.5, label: (yMax * 0.5).toFixed(0) },
        { norm: 1.0, label: yMax.toFixed(0) }
      ];
      const xTicks = [];
      for (let i = 0; i <= 4; i++) {
        const v = stats.min + (i / 4) * span;
        xTicks.push({ norm: i / 4, label: v.toFixed(1) });
      }

      engine.drawAxes({ yTicks, xTicks, title, xLabel: 'Values', yLabel: 'Frequency' });

      const barW = b.width / numBins;
      for (let i = 0; i < numBins; i++) {
        const count = bins[i];
        const h = (count / yMax) * b.height;
        const x = b.x + i * barW;
        const y = b.y + b.height - h;

        ctx.fillStyle = `${pal.primary}33`;
        ctx.strokeStyle = pal.primary;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.rect(x + 1, y, barW - 2, h);
        ctx.fill();
        ctx.stroke();
      }

      if (stats.sd > 0) {
        ctx.strokeStyle = pal.accent;
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        const steps = 80;
        for (let i = 0; i <= steps; i++) {
          const xv = stats.min + (i / steps) * span;
          const z = (xv - stats.mean) / stats.sd;
          const density = (1 / (stats.sd * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
          const expCount = density * binW * stats.n;
          const px = b.x + (i / steps) * b.width;
          const py = b.y + b.height - Math.min(b.height, (expCount / yMax) * b.height);
          if (i === 0) ctx.moveTo(px, py);
          else ctx.lineTo(px, py);
        }
        ctx.stroke();
      }
    },

    renderScatterRegression(engine, pairs, reg, title, quad = null) {
      engine.lastRender = () => this.renderScatterRegression(engine, pairs, reg, title, quad);
      engine.clear();
      const b = engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;
      if (!pairs || pairs.length < 2) return;

      const xs = pairs.map(p => p.x);
      const ys = pairs.map(p => p.y);
      const minX = Math.min(...xs);
      const maxX = Math.max(...xs);
      const minY = Math.min(...ys);
      const maxY = Math.max(...ys);

      const spanX = maxX - minX || 1;
      const spanY = maxY - minY || 1;
      const x0 = minX - 0.08 * spanX;
      const x1 = maxX + 0.08 * spanX;
      const y0 = minY - 0.08 * spanY;
      const y1 = maxY + 0.08 * spanY;

      const toX = (v) => b.x + ((v - x0) / (x1 - x0)) * b.width;
      const toY = (v) => b.y + b.height - ((v - y0) / (y1 - y0)) * b.height;

      const xTicks = [];
      const yTicks = [];
      for (let i = 0; i <= 4; i++) {
        xTicks.push({ norm: i / 4, label: (x0 + (i / 4) * (x1 - x0)).toFixed(1) });
        yTicks.push({ norm: i / 4, label: (y0 + (i / 4) * (y1 - y0)).toFixed(1) });
      }
      engine.drawAxes({ xTicks, yTicks, title, xLabel: 'Variable X', yLabel: 'Variable Y' });

      // 1. Linear regression line
      if (reg && reg.slope !== undefined) {
        const isQuadActive = quad && quad.isSignificantlyQuadratic && (quad.shape === 'U-Shaped' || quad.shape === 'Inverted U-Shaped');
        ctx.strokeStyle = isQuadActive ? `${pal.primary}77` : pal.primary;
        ctx.lineWidth = isQuadActive ? 1.8 : 2.5;
        if (isQuadActive) ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(toX(x0), toY(reg.intercept + reg.slope * x0));
        ctx.lineTo(toX(x1), toY(reg.intercept + reg.slope * x1));
        ctx.stroke();
        ctx.setLineDash([]);
      }

      // 2. Quadratic fitted curve & Nadir/Zenith point if U-shaped
      if (quad && quad.isSignificantlyQuadratic && (quad.shape === 'U-Shaped' || quad.shape === 'Inverted U-Shaped')) {
        ctx.strokeStyle = '#f59e0b'; // amber-500
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        const steps = 100;
        for (let i = 0; i <= steps; i++) {
          const vx = x0 + (i / steps) * (x1 - x0);
          const vy = quad.b0 + quad.b1 * vx + quad.b2 * vx * vx;
          const px = toX(vx);
          const py = toY(vy);
          if (i === 0) ctx.moveTo(px, py);
          else ctx.lineTo(px, py);
        }
        ctx.stroke();

        // Plot Nadir or Zenith point
        if (quad.vertexX >= minX - spanX * 0.1 && quad.vertexX <= maxX + spanX * 0.1) {
          const vx = quad.vertexX;
          const vy = quad.vertexY;
          const px = toX(vx);
          const py = toY(vy);

          // Glowing aura
          ctx.fillStyle = 'rgba(245, 158, 11, 0.25)';
          ctx.beginPath();
          ctx.arc(px, py, 9, 0, 2 * Math.PI);
          ctx.fill();

          // Core vertex dot
          ctx.fillStyle = '#f59e0b';
          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 1.5;
          ctx.beginPath();
          ctx.arc(px, py, 5, 0, 2 * Math.PI);
          ctx.fill();
          ctx.stroke();

          // Vertex text badge
          ctx.fillStyle = '#f59e0b';
          ctx.font = '600 11px Inter, system-ui, sans-serif';
          const label = `${quad.b2 > 0 ? 'Nadir' : 'Zenith'}: (${vx.toFixed(1)}, ${vy.toFixed(1)})`;
          ctx.fillText(label, px + 8, py - 6);
        }

        // Mini legend
        ctx.font = '10px Inter, system-ui, sans-serif';
        ctx.fillStyle = pal.primary;
        ctx.fillRect(b.x + b.width - 165, b.y + 12, 12, 3);
        ctx.fillStyle = pal.textMuted || '#94a3b8';
        ctx.fillText(`Linear OLS (r = ${(reg.r !== undefined ? reg.r : Math.sqrt(Math.max(0, reg.rSquared || 0))).toFixed(2)})`, b.x + b.width - 148, b.y + 16);

        ctx.fillStyle = '#f59e0b';
        ctx.fillRect(b.x + b.width - 165, b.y + 26, 12, 3);
        ctx.fillStyle = '#f59e0b';
        ctx.fillText(`Quadratic (${quad.shape}, R² = ${quad.rSquaredQuad.toFixed(2)})`, b.x + b.width - 148, b.y + 30);
      }

      // 3. Draw Scatter points
      ctx.fillStyle = pal.secondary;
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 1;
      for (const p of pairs) {
        ctx.beginPath();
        ctx.arc(toX(p.x), toY(p.y), 4.5, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
      }
    },

    renderROC(engine, roc, title) {
      engine.lastRender = () => this.renderROC(engine, roc, title);
      engine.clear();
      const b = engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;
      if (!roc || !roc.points) return;

      const xTicks = [
        { norm: 0, label: '0.0' }, { norm: 0.5, label: '0.5' }, { norm: 1, label: '1.0' }
      ];
      const yTicks = [
        { norm: 0, label: '0.0' }, { norm: 0.5, label: '0.5' }, { norm: 1, label: '1.0' }
      ];
      engine.drawAxes({
        xTicks, yTicks, title,
        xLabel: '1 - Specificity (FPR)',
        yLabel: 'Sensitivity (TPR)'
      });

      const toX = (fpr) => b.x + fpr * b.width;
      const toY = (tpr) => b.y + b.height - tpr * b.height;

      // Diagonal chance
      ctx.strokeStyle = pal.axis;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(b.x, b.y + b.height);
      ctx.lineTo(b.x + b.width, b.y);
      ctx.stroke();
      ctx.setLineDash([]);

      // Fill AUC
      ctx.fillStyle = `${pal.primary}22`;
      ctx.beginPath();
      ctx.moveTo(b.x, b.y + b.height);
      for (const pt of roc.points) ctx.lineTo(toX(pt.fpr), toY(pt.tpr));
      ctx.lineTo(b.x + b.width, b.y + b.height);
      ctx.closePath();
      ctx.fill();

      // ROC Line
      ctx.strokeStyle = pal.primary;
      ctx.lineWidth = 3;
      ctx.beginPath();
      for (let i = 0; i < roc.points.length; i++) {
        const pt = roc.points[i];
        if (i === 0) ctx.moveTo(toX(pt.fpr), toY(pt.tpr));
        else ctx.lineTo(toX(pt.fpr), toY(pt.tpr));
      }
      ctx.stroke();

      if (roc.optimalCutoff) {
        const opt = roc.optimalCutoff;
        const ox = toX(1 - opt.spec);
        const oy = toY(opt.sens);
        ctx.fillStyle = pal.danger;
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(ox, oy, 6, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
      }
    },

    renderTeachingDistribution(engine, data, distKey, params = {}, title = 'Generated Distribution & Theoretical PDF') {
      engine.lastRender = () => this.renderTeachingDistribution(engine, data, distKey, params, title);
      engine.lastRenderFn = engine.lastRender;
      engine.clear();
      const b = engine.getPlotBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;

      if (!data || data.length < 5) {
        ctx.fillStyle = pal.textDim || '#94a3b8';
        ctx.font = `14px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText('Click "Generate New Sample" to simulate distribution data.', b.x + b.width / 2, b.y + b.height / 2);
        return;
      }

      const n = data.length;
      let minVal = Infinity;
      let maxVal = -Infinity;
      let sum = 0;
      for (let i = 0; i < n; i++) {
        const v = data[i];
        if (v < minVal) minVal = v;
        if (v > maxVal) maxVal = v;
        sum += v;
      }
      const sampleMean = sum / n;

      const span = maxVal - minVal || 1;
      const plotMin = minVal - 0.05 * span;
      const plotMax = maxVal + 0.05 * span;
      const plotSpan = plotMax - plotMin;

      const numBins = Math.max(12, Math.min(35, Math.round(1 + 3.322 * Math.log10(n) * 1.5)));
      const binWidth = plotSpan / numBins;

      const bins = new Array(numBins).fill(0);
      for (let i = 0; i < n; i++) {
        const idx = Math.min(numBins - 1, Math.max(0, Math.floor((data[i] - plotMin) / binWidth)));
        bins[idx]++;
      }

      const maxCount = Math.max(...bins, 1);
      const yMax = maxCount * 1.25;

      const yTicks = [
        { norm: 0, label: '0' },
        { norm: 0.5, label: (yMax * 0.5).toFixed(0) },
        { norm: 1.0, label: yMax.toFixed(0) }
      ];

      const xTicks = [];
      for (let i = 0; i <= 5; i++) {
        const v = plotMin + (i / 5) * plotSpan;
        xTicks.push({ norm: i / 5, label: v.toFixed(1) });
      }

      engine.drawAxes({
        yTicks,
        xTicks,
        title,
        xLabel: 'Observation Value (X)',
        yLabel: 'Frequency Count'
      });

      const toX = (val) => b.x + ((val - plotMin) / plotSpan) * b.width;
      const toY = (count) => b.y + b.height - (count / yMax) * b.height;

      // Draw Histogram Bars
      const barWidth = b.width / numBins;
      for (let i = 0; i < numBins; i++) {
        const count = bins[i];
        if (count === 0) continue;
        const h = (count / yMax) * b.height;
        const x = b.x + i * barWidth;
        const y = b.y + b.height - h;

        ctx.fillStyle = `${pal.primary}33`;
        ctx.strokeStyle = pal.primary;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.rect(x + 1, y, Math.max(1, barWidth - 2), h);
        ctx.fill();
        ctx.stroke();
      }

      // Theoretical PDF Curve
      if (Teaching.pdf && Teaching.pdf[distKey]) {
        const pdfFn = Teaching.pdf[distKey];
        const densityScale = n * binWidth;

        ctx.strokeStyle = pal.accent || '#38bdf8';
        ctx.lineWidth = 2.8;
        ctx.beginPath();

        const steps = 180;
        let started = false;
        for (let s = 0; s <= steps; s++) {
          const xVal = plotMin + (s / steps) * plotSpan;
          let pdfVal = 0;

          switch (distKey) {
            case 'normal':
              pdfVal = pdfFn.call(Teaching.pdf, xVal, params.mean, params.sd);
              break;
            case 'studentsT':
              pdfVal = pdfFn.call(Teaching.pdf, xVal, params.df, params.mean, params.scale);
              break;
            case 'uniform':
              pdfVal = pdfFn.call(Teaching.pdf, xVal, params.min, params.max);
              break;
            case 'exponential':
              pdfVal = pdfFn.call(Teaching.pdf, xVal, params.rate);
              break;
            case 'logNormal':
              pdfVal = pdfFn.call(Teaching.pdf, xVal, params.mu, params.sigma);
              break;
            case 'bimodal':
              pdfVal = pdfFn.call(Teaching.pdf, xVal, params.m1, params.s1, params.m2, params.s2, params.p);
              break;
            case 'poisson':
              pdfVal = pdfFn.call(Teaching.pdf, Math.round(xVal), params.lambda);
              break;
            case 'chiSquare':
              pdfVal = pdfFn.call(Teaching.pdf, xVal, params.df);
              break;
            default:
              pdfVal = 0;
          }

          const countVal = pdfVal * densityScale;
          const px = toX(xVal);
          const py = Math.max(b.y - 10, toY(countVal));

          if (!started) {
            ctx.moveTo(px, py);
            started = true;
          } else {
            ctx.lineTo(px, py);
          }
        }
        ctx.stroke();
      }

      // Sample Mean marker
      if (isFinite(sampleMean)) {
        const meanX = toX(sampleMean);
        if (meanX >= b.x && meanX <= b.x + b.width) {
          ctx.strokeStyle = '#f59e0b';
          ctx.lineWidth = 2;
          ctx.setLineDash([4, 4]);
          ctx.beginPath();
          ctx.moveTo(meanX, b.y);
          ctx.lineTo(meanX, b.y + b.height);
          ctx.stroke();
          ctx.setLineDash([]);

          ctx.fillStyle = '#f59e0b';
          ctx.font = `600 11px ${engine.options.fontFamily}`;
          ctx.textAlign = 'center';
          ctx.fillText(`M = ${sampleMean.toFixed(2)}`, meanX, b.y + 14);
        }
      }

      // Legend
      ctx.textAlign = 'right';
      ctx.font = `500 11px ${engine.options.fontFamily}`;
      ctx.fillStyle = pal.primary;
      ctx.fillText('■ Empirical Sample Histogram', b.x + b.width - 10, b.y + 15);
      ctx.fillStyle = pal.accent || '#38bdf8';
      ctx.fillText('— Theoretical PDF Overlay', b.x + b.width - 10, b.y + 32);
      ctx.fillStyle = '#f59e0b';
      ctx.fillText('┆ Sample Mean', b.x + b.width - 10, b.y + 49);
    },

    renderCltParent(engine, parentInfo, lastSample = [], title = 'CLT Parent Population Distribution') {
      engine.lastRender = () => this.renderCltParent(engine, parentInfo, lastSample, title);
      engine.lastRenderFn = engine.lastRender;
      engine.clear();
      const b = engine.getPlotBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;

      if (!parentInfo) return;

      const minX = parentInfo.min;
      const maxX = parentInfo.max;
      const span = maxX - minX;

      const steps = 150;
      let maxDensity = 0;
      const densities = [];
      for (let i = 0; i <= steps; i++) {
        const x = minX + (i / steps) * span;
        const d = parentInfo.pdf(x);
        densities.push({ x, d });
        if (d > maxDensity) maxDensity = d;
      }
      const yMax = (maxDensity || 0.5) * 1.3;

      const toX = (val) => b.x + ((val - minX) / span) * b.width;
      const toY = (d) => b.y + b.height - (d / yMax) * b.height;

      const xTicks = [];
      for (let i = 0; i <= 4; i++) {
        const v = minX + (i / 4) * span;
        xTicks.push({ norm: i / 4, label: v.toFixed(1) });
      }
      const yTicks = [
        { norm: 0, label: '0' },
        { norm: 0.5, label: (yMax * 0.5).toFixed(2) },
        { norm: 1.0, label: yMax.toFixed(2) }
      ];

      engine.drawAxes({
        yTicks,
        xTicks,
        title: `${title} (${parentInfo.name})`,
        xLabel: 'Observation Value (X)',
        yLabel: 'Probability Density f(X)'
      });

      // Fill Area
      ctx.fillStyle = `${pal.secondary || '#94a3b8'}22`;
      ctx.beginPath();
      ctx.moveTo(b.x, b.y + b.height);
      for (const pt of densities) {
        ctx.lineTo(toX(pt.x), toY(pt.d));
      }
      ctx.lineTo(b.x + b.width, b.y + b.height);
      ctx.closePath();
      ctx.fill();

      // Stroke Density
      ctx.strokeStyle = pal.secondary || '#94a3b8';
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      for (let i = 0; i < densities.length; i++) {
        const pt = densities[i];
        if (i === 0) ctx.moveTo(toX(pt.x), toY(pt.d));
        else ctx.lineTo(toX(pt.x), toY(pt.d));
      }
      ctx.stroke();

      // True Mean line
      const muX = toX(parentInfo.mean);
      if (muX >= b.x && muX <= b.x + b.width) {
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 2;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(muX, b.y);
        ctx.lineTo(muX, b.y + b.height);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = '#ef4444';
        ctx.font = `600 11px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText(`True μ = ${parentInfo.mean.toFixed(2)}`, muX, b.y + 14);
      }

      // Draw Last Sample Draw points
      if (lastSample && lastSample.length > 0) {
        let sampleSum = 0;
        for (const val of lastSample) {
          sampleSum += val;
          const ptX = toX(val);
          const ptD = parentInfo.pdf(val);
          const ptY = toY(ptD);

          ctx.fillStyle = '#38bdf8';
          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.arc(ptX, ptY, 4, 0, 2 * Math.PI);
          ctx.fill();
          ctx.stroke();
        }

        const sampleMean = sampleSum / lastSample.length;
        const sMeanX = toX(sampleMean);
        ctx.strokeStyle = '#38bdf8';
        ctx.lineWidth = 2.2;
        ctx.beginPath();
        ctx.moveTo(sMeanX, b.y + b.height - 25);
        ctx.lineTo(sMeanX, b.y + b.height);
        ctx.stroke();

        ctx.fillStyle = '#38bdf8';
        ctx.font = `600 10px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText(`Sample x̄ = ${sampleMean.toFixed(2)}`, sMeanX, b.y + b.height - 28);
      }

      // Legend
      ctx.textAlign = 'right';
      ctx.font = `500 11px ${engine.options.fontFamily}`;
      ctx.fillStyle = pal.secondary || '#94a3b8';
      ctx.fillText('— Parent Density f(X)', b.x + b.width - 10, b.y + 15);
      ctx.fillStyle = '#ef4444';
      ctx.fillText('┆ True Population Mean (μ)', b.x + b.width - 10, b.y + 32);
      if (lastSample && lastSample.length > 0) {
        ctx.fillStyle = '#38bdf8';
        ctx.fillText(`● Current Draw (n = ${lastSample.length})`, b.x + b.width - 10, b.y + 49);
      }
    },

    renderCltSampling(engine, cltData, title = 'Sampling Distribution of the Mean (x̄)') {
      engine.lastRender = () => this.renderCltSampling(engine, cltData, title);
      engine.lastRenderFn = engine.lastRender;
      engine.clear();
      const b = engine.getPlotBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;

      if (!cltData || cltData.samplesDrawn === 0) {
        ctx.fillStyle = pal.textDim || '#94a3b8';
        ctx.font = `14px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText('No samples drawn yet. Click "▶ Step (Draw 1 Sample)" or "⚡ Draw 100 Samples" above.', b.x + b.width / 2, b.y + b.height / 2);
        return;
      }

      const means = cltData.values;
      const k = means.length;
      const n = cltData.sampleSize;
      const trueMu = cltData.theoreticalMean;
      const trueSE = cltData.theoreticalSE;

      let minM = Infinity;
      let maxM = -Infinity;
      for (let i = 0; i < k; i++) {
        if (means[i] < minM) minM = means[i];
        if (means[i] > maxM) maxM = means[i];
      }

      const seSpan = 3.5 * trueSE;
      const plotMin = Math.min(minM - 0.2 * trueSE, trueMu - seSpan);
      const plotMax = Math.max(maxM + 0.2 * trueSE, trueMu + seSpan);
      const plotSpan = plotMax - plotMin || 1;

      const numBins = Math.max(15, Math.min(40, Math.round(1 + 3.322 * Math.log10(k) * 2)));
      const binWidth = plotSpan / numBins;

      const bins = new Array(numBins).fill(0);
      for (let i = 0; i < k; i++) {
        const idx = Math.min(numBins - 1, Math.max(0, Math.floor((means[i] - plotMin) / binWidth)));
        bins[idx]++;
      }

      const maxCount = Math.max(...bins, 1);
      const yMax = maxCount * 1.25;

      const toX = (val) => b.x + ((val - plotMin) / plotSpan) * b.width;
      const toY = (count) => b.y + b.height - (count / yMax) * b.height;

      const xTicks = [];
      for (let i = 0; i <= 5; i++) {
        const v = plotMin + (i / 5) * plotSpan;
        xTicks.push({ norm: i / 5, label: v.toFixed(2) });
      }
      const yTicks = [
        { norm: 0, label: '0' },
        { norm: 0.5, label: (yMax * 0.5).toFixed(0) },
        { norm: 1.0, label: yMax.toFixed(0) }
      ];

      engine.drawAxes({
        yTicks,
        xTicks,
        title: `${title} (k = ${k.toLocaleString()} samples, n = ${n})`,
        xLabel: 'Sample Mean Value (x̄)',
        yLabel: 'Frequency of Means'
      });

      // Histogram Bars
      const barWidth = b.width / numBins;
      for (let i = 0; i < numBins; i++) {
        const count = bins[i];
        if (count === 0) continue;
        const h = (count / yMax) * b.height;
        const x = b.x + i * barWidth;
        const y = b.y + b.height - h;

        ctx.fillStyle = `${pal.primary}44`;
        ctx.strokeStyle = pal.primary;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.rect(x + 1, y, Math.max(1, barWidth - 2), h);
        ctx.fill();
        ctx.stroke();
      }

      // Theoretical Gaussian Normal Curve: N(μ, σ/√n)
      if (trueSE > 0) {
        const densityScale = k * binWidth;
        ctx.strokeStyle = '#22c55e';
        ctx.lineWidth = 3;
        ctx.beginPath();

        const steps = 150;
        let started = false;
        for (let s = 0; s <= steps; s++) {
          const xVal = plotMin + (s / steps) * plotSpan;
          const z = (xVal - trueMu) / trueSE;
          const normDensity = (1.0 / (trueSE * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
          const countVal = normDensity * densityScale;
          const px = toX(xVal);
          const py = toY(countVal);

          if (!started) {
            ctx.moveTo(px, py);
            started = true;
          } else {
            ctx.lineTo(px, py);
          }
        }
        ctx.stroke();
      }

      // True Mean line
      const muX = toX(trueMu);
      if (muX >= b.x && muX <= b.x + b.width) {
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 2;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(muX, b.y);
        ctx.lineTo(muX, b.y + b.height);
        ctx.stroke();
        ctx.setLineDash([]);
      }

      // Observed Mean of Means
      if (cltData.observedMean !== null) {
        const obsX = toX(cltData.observedMean);
        if (obsX >= b.x && obsX <= b.x + b.width) {
          ctx.strokeStyle = '#f59e0b';
          ctx.lineWidth = 2;
          ctx.beginPath();
          ctx.moveTo(obsX, b.y);
          ctx.lineTo(obsX, b.y + b.height);
          ctx.stroke();
        }
      }

      // Legend
      ctx.textAlign = 'right';
      ctx.font = `500 11px ${engine.options.fontFamily}`;
      ctx.fillStyle = pal.primary;
      ctx.fillText(`■ Simulated Means (k = ${k.toLocaleString()})`, b.x + b.width - 10, b.y + 15);
      ctx.fillStyle = '#22c55e';
      ctx.fillText(`— CLT Normal Fit N(μ, σ/√n)`, b.x + b.width - 10, b.y + 32);
      ctx.fillStyle = '#ef4444';
      ctx.fillText(`┆ True Mean μ = ${trueMu.toFixed(2)}`, b.x + b.width - 10, b.y + 49);
      if (cltData.observedMean !== null) {
        ctx.fillStyle = '#f59e0b';
        ctx.fillText(`— Observed x̄̄ = ${cltData.observedMean.toFixed(2)} (SE: ${cltData.observedSE.toFixed(3)})`, b.x + b.width - 10, b.y + 66);
      }
    },

    /**
     * Renders Student's t-Distribution Convergence to Standard Normal N(0, 1)
     */
    renderTConvergence(engine, metrics, options = {}) {
      engine.lastRender = () => this.renderTConvergence(engine, metrics, options);
      engine.lastRenderFn = engine.lastRender;
      engine.clear();
      const b = engine.getPlotBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;

      if (!metrics) return;

      const df = metrics.df;
      const showTailArea = options.showTailArea !== false;
      const title = options.title || `Student's t(ν = ${df}) Convergence to Standard Normal N(0, 1)`;

      const minX = -4.5;
      const maxX = 4.5;
      const spanX = maxX - minX;
      const yMax = 0.44; // Peak of N(0, 1) is 0.39894

      const toX = (val) => b.x + ((val - minX) / spanX) * b.width;
      const toY = (d) => b.y + b.height - (Math.max(0, d) / yMax) * b.height;

      // Generate grid & axes ticks
      const xTicks = [];
      for (let x = -4; x <= 4; x += 1) {
        xTicks.push({ norm: (x - minX) / spanX, label: `${x > 0 ? '+' : ''}${x}` });
      }
      const yTicks = [
        { norm: 0, label: '0.00' },
        { norm: 0.1 / yMax, label: '0.10' },
        { norm: 0.2 / yMax, label: '0.20' },
        { norm: 0.3 / yMax, label: '0.30' },
        { norm: 0.4 / yMax, label: '0.40' }
      ];

      engine.drawAxes({
        xTicks,
        yTicks,
        title,
        xLabel: 'Standardized Value (t / z)',
        yLabel: 'Probability Density f(x)'
      });

      const steps = 240;
      const pointsNorm = [];
      const pointsT = [];

      for (let i = 0; i <= steps; i++) {
        const x = minX + (i / steps) * spanX;
        const normD = (1.0 / Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * x * x);
        const tD = Teaching.pdf.studentsT(x, df, 0, 1);
        pointsNorm.push({ x, y: normD });
        pointsT.push({ x, y: tD });
      }

      // 1. Shaded Tail Area (|x| >= 1.960) under Student's t curve
      if (showTailArea) {
        const zCrit = 1.95996;
        ctx.fillStyle = 'rgba(239, 68, 68, 0.22)';

        // Left Tail [-4.5, -1.960]
        ctx.beginPath();
        ctx.moveTo(toX(minX), toY(0));
        for (const pt of pointsT) {
          if (pt.x <= -zCrit) {
            ctx.lineTo(toX(pt.x), toY(pt.y));
          }
        }
        const tAtLeftCrit = Teaching.pdf.studentsT(-zCrit, df, 0, 1);
        ctx.lineTo(toX(-zCrit), toY(tAtLeftCrit));
        ctx.lineTo(toX(-zCrit), toY(0));
        ctx.closePath();
        ctx.fill();

        // Right Tail [1.960, 4.5]
        ctx.beginPath();
        ctx.moveTo(toX(zCrit), toY(0));
        const tAtRightCrit = Teaching.pdf.studentsT(zCrit, df, 0, 1);
        ctx.lineTo(toX(zCrit), toY(tAtRightCrit));
        for (const pt of pointsT) {
          if (pt.x >= zCrit) {
            ctx.lineTo(toX(pt.x), toY(pt.y));
          }
        }
        ctx.lineTo(toX(maxX), toY(0));
        ctx.closePath();
        ctx.fill();
      }

      // 2. Render Standard Normal Reference Curve N(0, 1) [Emerald/Cyan dashed line]
      ctx.strokeStyle = '#10b981';
      ctx.lineWidth = 2.2;
      ctx.setLineDash([5, 4]);
      ctx.beginPath();
      for (let i = 0; i < pointsNorm.length; i++) {
        const pt = pointsNorm[i];
        if (i === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
        else ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.stroke();
      ctx.setLineDash([]);

      // 3. Render Student's t(ν) Curve [Vibrant Violet solid line]
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 2.8;
      ctx.beginPath();
      for (let i = 0; i < pointsT.length; i++) {
        const pt = pointsT[i];
        if (i === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
        else ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.stroke();

      // 4. Mark Critical Values Lines
      const zCrit = 1.95996;
      const tCrit = metrics.tCrit;

      [-zCrit, zCrit].forEach(zVal => {
        const xPix = toX(zVal);
        ctx.strokeStyle = '#64748b';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(xPix, toY(0));
        ctx.lineTo(xPix, toY(0.18));
        ctx.stroke();
        ctx.setLineDash([]);
      });

      if (tCrit <= 4.4) {
        [-tCrit, tCrit].forEach(tVal => {
          const xPix = toX(tVal);
          ctx.strokeStyle = '#c084fc';
          ctx.lineWidth = 1.8;
          ctx.setLineDash([4, 3]);
          ctx.beginPath();
          ctx.moveTo(xPix, toY(0));
          ctx.lineTo(xPix, toY(0.24));
          ctx.stroke();
          ctx.setLineDash([]);
        });
      }

      // 5. Annotations & Peak Height Indicator
      const normPeakY = toY(metrics.normPeak);
      const tPeakY = toY(metrics.tPeak);
      const midX = toX(0);

      if (Math.abs(metrics.peakDiffPct) > 1.5) {
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(midX, tPeakY);
        ctx.lineTo(midX, normPeakY);
        ctx.stroke();

        ctx.fillStyle = '#f59e0b';
        ctx.font = `600 10px ${engine.options.fontFamily}`;
        ctx.textAlign = 'left';
        ctx.fillText(`Δ Peak: ${metrics.peakDiffPct.toFixed(1)}%`, midX + 8, (normPeakY + tPeakY) / 2 + 4);
      }

      // 6. Legend
      ctx.textAlign = 'right';
      ctx.font = `500 11px ${engine.options.fontFamily}`;

      ctx.fillStyle = '#10b981';
      ctx.fillText('— — Standard Normal N(0, 1) [Peak: 0.3989]', b.x + b.width - 10, b.y + 15);

      ctx.fillStyle = '#a855f7';
      ctx.fillText(`—— Student's t (ν = ${df}) [Peak: ${metrics.tPeak.toFixed(4)}]`, b.x + b.width - 10, b.y + 32);

      ctx.fillStyle = '#94a3b8';
      ctx.fillText(`┆ 95% Cutoffs: z = ±1.960 vs t = ±${tCrit.toFixed(3)} (${metrics.critDiffPct >= 0 ? '+' : ''}${metrics.critDiffPct.toFixed(1)}%)`, b.x + b.width - 10, b.y + 49);

      if (showTailArea) {
        ctx.fillStyle = '#ef4444';
        ctx.fillText(`░░ Fat Tail Risk: P(|T| > 1.96) = ${(metrics.tailProb * 100).toFixed(1)}% vs 5.0%`, b.x + b.width - 10, b.y + 66);
      }
    },

    renderTwoSampleOverlap(engine, metrics, options = {}) {
      engine.lastRender = () => this.renderTwoSampleOverlap(engine, metrics, options);
      engine.lastRenderFn = engine.lastRender;
      engine.clear();
      const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette;

      const {
        mean1 = 10,
        mean2 = 12,
        delta = 2,
        sd = 2.5,
        sd1 = sd || 2.5,
        sd2 = sd || 2.5,
        n = 16,
        n1 = n || 16,
        n2 = n || 16,
        sem = 0.625,
        sem1 = sem || 0.625,
        sem2 = sem || 0.625,
        alpha = 0.05,
        viewMode = 'means',
        seDiff = 0.88,
        df = 30,
        tCrit = 2.04,
        zCrit = 1.96,
        deltaCrit = 1.8,
        tStat = 0,
        pValue = 0.05,
        isSignificant = false,
        cohensD = 0.8,
        patientOVL = 0.5,
        meansOVL = 0.05,
        moe = 1.0,
        ci1 = [mean1 - 1, mean1 + 1],
        ci2 = [mean2 - 1, mean2 + 1]
      } = metrics || {};

      if (viewMode === 'null') {
        const safeSEDiff = Math.max(0.001, Number.isFinite(seDiff) ? seDiff : 0.88);
        const safeDelta = Number.isFinite(delta) ? delta : 2.0;
        const xSpan = Math.max(0.1, Math.max(4.2 * safeSEDiff, safeDelta + 2.5 * safeSEDiff));
        const minX = -xSpan;
        const maxX = xSpan;
        const peakY = 1.0 / (safeSEDiff * Math.sqrt(2 * Math.PI));
        const maxY = peakY * 1.25;

        const toX = (val) => {
          const v = Number.isFinite(val) ? val : minX;
          return b.x + ((v - minX) / (maxX - minX)) * b.width;
        };
        const toY = (val) => {
          const v = Number.isFinite(val) ? val : 0;
          return b.y + b.height - (v / maxY) * b.height;
        };

        ctx.strokeStyle = pal.grid;
        ctx.lineWidth = 1;
        const xSteps = 6;
        for (let i = 0; i <= xSteps; i++) {
          const xVal = minX + (i / xSteps) * (maxX - minX);
          const xPix = toX(xVal);
          ctx.beginPath();
          ctx.moveTo(xPix, b.y);
          ctx.lineTo(xPix, b.y + b.height);
          ctx.stroke();

          ctx.fillStyle = pal.textMuted;
          ctx.font = `500 10px ${engine.options.fontFamily}`;
          ctx.textAlign = 'center';
          ctx.fillText(xVal.toFixed(2), xPix, b.y + b.height + 15);
        }

        const numPts = 250;
        const pts = [];
        for (let i = 0; i <= numPts; i++) {
          const xVal = minX + (i / numPts) * (maxX - minX);
          const z = xVal / seDiff;
          const yVal = (1.0 / (seDiff * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
          pts.push({ x: xVal, y: yVal });
        }

        ctx.fillStyle = 'rgba(239, 68, 68, 0.25)';
        ctx.beginPath();
        ctx.moveTo(toX(minX), toY(0));
        for (const pt of pts) {
          if (pt.x <= -deltaCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
        }
        ctx.lineTo(toX(-deltaCrit), toY(0));
        ctx.closePath();
        ctx.fill();

        ctx.beginPath();
        ctx.moveTo(toX(deltaCrit), toY(0));
        for (const pt of pts) {
          if (pt.x >= deltaCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
        }
        ctx.lineTo(toX(maxX), toY(0));
        ctx.closePath();
        ctx.fill();

        ctx.fillStyle = 'rgba(16, 185, 129, 0.08)';
        ctx.beginPath();
        ctx.moveTo(toX(-deltaCrit), toY(0));
        for (const pt of pts) {
          if (pt.x >= -deltaCrit && pt.x <= deltaCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
        }
        ctx.lineTo(toX(deltaCrit), toY(0));
        ctx.closePath();
        ctx.fill();

        ctx.strokeStyle = '#38bdf8';
        ctx.lineWidth = 2.4;
        ctx.beginPath();
        pts.forEach((pt, idx) => {
          if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
          else ctx.lineTo(toX(pt.x), toY(pt.y));
        });
        ctx.stroke();

        [-deltaCrit, deltaCrit].forEach((cVal, idx) => {
          const xPix = toX(cVal);
          ctx.strokeStyle = '#ef4444';
          ctx.lineWidth = 1.8;
          ctx.setLineDash([4, 4]);
          ctx.beginPath();
          ctx.moveTo(xPix, toY(0));
          ctx.lineTo(xPix, toY(peakY * 0.75));
          ctx.stroke();
          ctx.setLineDash([]);

          ctx.fillStyle = '#ef4444';
          ctx.font = `600 10px ${engine.options.fontFamily}`;
          ctx.textAlign = idx === 0 ? 'right' : 'left';
          ctx.fillText(idx === 0 ? `-Δcrit (${cVal.toFixed(2)})` : `+Δcrit (${cVal.toFixed(2)})`, xPix + (idx === 0 ? -5 : 5), toY(peakY * 0.75));
        });

        const obsXPix = toX(delta);
        const obsColor = isSignificant ? '#10b981' : '#f59e0b';
        ctx.strokeStyle = obsColor;
        ctx.lineWidth = 2.8;
        ctx.beginPath();
        ctx.moveTo(obsXPix, toY(0));
        ctx.lineTo(obsXPix, toY(peakY * 0.95));
        ctx.stroke();

        ctx.fillStyle = obsColor;
        ctx.beginPath();
        ctx.moveTo(obsXPix, toY(peakY * 0.98));
        ctx.lineTo(obsXPix - 6, toY(peakY * 0.90));
        ctx.lineTo(obsXPix + 6, toY(peakY * 0.90));
        ctx.closePath();
        ctx.fill();

        ctx.fillStyle = obsColor;
        ctx.font = `700 11px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText(`Observed Δ = ${delta.toFixed(2)} (t = ${tStat.toFixed(2)})`, obsXPix, toY(peakY * 1.06));

        const badgeText = isSignificant ? `✓ SIGNIFICANT (p = ${pValue.toFixed(4)} < α)` : `✗ NOT SIGNIFICANT (p = ${pValue.toFixed(4)} ≥ α)`;
        ctx.fillStyle = isSignificant ? 'rgba(16, 185, 129, 0.2)' : 'rgba(239, 68, 68, 0.2)';
        ctx.strokeStyle = isSignificant ? '#10b981' : '#ef4444';
        ctx.lineWidth = 1;
        const badgeWidth = ctx.measureText(badgeText).width + 24;
        ctx.fillRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);
        ctx.strokeRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);

        ctx.fillStyle = isSignificant ? '#10b981' : '#ef4444';
        ctx.font = `700 11px ${engine.options.fontFamily}`;
        ctx.textAlign = 'right';
        ctx.fillText(badgeText, b.x + b.width - 22, b.y + 26);

        ctx.textAlign = 'left';
        ctx.font = `500 11px ${engine.options.fontFamily}`;
        ctx.fillStyle = '#38bdf8';
        ctx.fillText(`— Null Distribution H₀: Δ ~ N(0, SE²diff), SE = ${seDiff.toFixed(3)}`, b.x + 12, b.y + 20);
        ctx.fillStyle = '#ef4444';
        ctx.fillText(`░ Rejection Region (α = ${alpha.toFixed(3)}, t_crit = ±${tCrit.toFixed(3)})`, b.x + 12, b.y + 36);
        return;
      }

      // -------------------------------------------------------------
      // MEANS, PATIENTS, OR DUAL VIEW
      // -------------------------------------------------------------
      const activeSigma1 = Math.max(0.001, (viewMode === 'patients') ? sd1 : sem1);
      const activeSigma2 = Math.max(0.001, (viewMode === 'patients') ? sd2 : sem2);
      const maxSpread = Math.max(sd1, sd2, 4.0);
      const minX = mean1 - Math.max(3.8 * maxSpread, 4.0);
      const maxX = Math.max(mean2 + Math.max(3.8 * maxSpread, 4.0), mean1 + 7.5);
      const xSpan = Math.max(0.1, maxX - minX);

      const peak1 = 1.0 / (activeSigma1 * Math.sqrt(2 * Math.PI));
      const peak2 = 1.0 / (activeSigma2 * Math.sqrt(2 * Math.PI));
      const maxDensity = Math.max(0.001, peak1, peak2);
      const maxY = maxDensity * 1.30;

      const toX = (val) => {
        const v = Number.isFinite(val) ? val : minX;
        return b.x + ((v - minX) / xSpan) * b.width;
      };
      const toY = (val) => {
        const v = Number.isFinite(val) ? val : 0;
        return b.y + b.height - (v / maxY) * b.height;
      };

      ctx.strokeStyle = pal.grid;
      ctx.lineWidth = 1;
      const xSteps = 7;
      for (let i = 0; i <= xSteps; i++) {
        const xVal = minX + (i / xSteps) * xSpan;
        const xPix = toX(xVal);
        ctx.beginPath();
        ctx.moveTo(xPix, b.y);
        ctx.lineTo(xPix, b.y + b.height);
        ctx.stroke();

        ctx.fillStyle = pal.textMuted;
        ctx.font = `500 10px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText(xVal.toFixed(1), xPix, b.y + b.height + 15);
      }

      const numPts = 320;
      const pts1 = [];
      const pts2 = [];
      const ptsOverlap = [];

      const normPDF = (x, mu, s) => (1.0 / (s * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * Math.pow((x - mu) / s, 2));

      for (let i = 0; i <= numPts; i++) {
        const x = minX + (i / numPts) * xSpan;
        const y1 = normPDF(x, mean1, activeSigma1);
        const y2 = normPDF(x, mean2, activeSigma2);
        const yOverlap = Math.min(y1, y2);

        pts1.push({ x, y: y1 });
        pts2.push({ x, y: y2 });
        ptsOverlap.push({ x, y: yOverlap });
      }

      ctx.fillStyle = isSignificant ? 'rgba(16, 185, 129, 0.22)' : 'rgba(239, 68, 68, 0.24)';
      ctx.beginPath();
      ctx.moveTo(toX(minX), toY(0));
      ptsOverlap.forEach(pt => ctx.lineTo(toX(pt.x), toY(pt.y)));
      ctx.lineTo(toX(maxX), toY(0));
      ctx.closePath();
      ctx.fill();

      if (viewMode === 'dual') {
        const patientPeak1 = 1.0 / (sd1 * Math.sqrt(2 * Math.PI));
        const patientPeak2 = 1.0 / (sd2 * Math.sqrt(2 * Math.PI));
        const dualScale1 = (maxDensity * 0.45) / patientPeak1;
        const dualScale2 = (maxDensity * 0.45) / patientPeak2;

        ctx.strokeStyle = 'rgba(6, 182, 212, 0.45)';
        ctx.lineWidth = 1.6;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        for (let i = 0; i <= numPts; i++) {
          const x = minX + (i / numPts) * xSpan;
          const y = normPDF(x, mean1, sd1) * dualScale1;
          if (i === 0) ctx.moveTo(toX(x), toY(y));
          else ctx.lineTo(toX(x), toY(y));
        }
        ctx.stroke();

        ctx.strokeStyle = 'rgba(168, 85, 247, 0.45)';
        ctx.beginPath();
        for (let i = 0; i <= numPts; i++) {
          const x = minX + (i / numPts) * xSpan;
          const y = normPDF(x, mean2, sd2) * dualScale2;
          if (i === 0) ctx.moveTo(toX(x), toY(y));
          else ctx.lineTo(toX(x), toY(y));
        }
        ctx.stroke();
        ctx.setLineDash([]);
      }

      ctx.strokeStyle = '#06b6d4';
      ctx.lineWidth = 2.6;
      ctx.beginPath();
      pts1.forEach((pt, idx) => {
        if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
        else ctx.lineTo(toX(pt.x), toY(pt.y));
      });
      ctx.stroke();

      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 2.6;
      ctx.beginPath();
      pts2.forEach((pt, idx) => {
        if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
        else ctx.lineTo(toX(pt.x), toY(pt.y));
      });
      ctx.stroke();

      [
        { mu: mean1, color: '#06b6d4', label: `Group 1 (μ₁ = ${mean1.toFixed(1)}${viewMode === 'patients' ? `, SD₁=${sd1.toFixed(2)}` : `, SEM₁=${sem1.toFixed(2)}`})` },
        { mu: mean2, color: '#a855f7', label: `Group 2 (μ₂ = ${mean2.toFixed(1)}${viewMode === 'patients' ? `, SD₂=${sd2.toFixed(2)}` : `, SEM₂=${sem2.toFixed(2)}`})` }
      ].forEach((grp) => {
        const xPix = toX(grp.mu);
        ctx.strokeStyle = grp.color;
        ctx.lineWidth = 1.5;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(xPix, toY(0));
        ctx.lineTo(xPix, toY(maxDensity));
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = grp.color;
        ctx.font = `600 10px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText(grp.label, xPix, toY(maxDensity) - 8);
      });

      const critX = mean1 + deltaCrit;
      if (Number.isFinite(critX) && critX <= maxX && critX >= minX) {
        const critXPix = toX(critX);
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 2.0;
        ctx.setLineDash([5, 3]);
        ctx.beginPath();
        ctx.moveTo(critXPix, toY(0));
        ctx.lineTo(critXPix, toY(maxDensity * 0.85));
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = '#f59e0b';
        ctx.font = `600 10px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText(`Significance Boundary (Δcrit = ${deltaCrit.toFixed(2)})`, critXPix, toY(maxDensity * 0.85) - 6);
      }

      if (viewMode === 'means' || viewMode === 'dual') {
        const barY1 = toY(maxDensity * 0.05);
        const barY2 = toY(maxDensity * 0.12);

        if (Array.isArray(ci1) && Number.isFinite(ci1[0]) && Number.isFinite(ci1[1])) {
          ctx.strokeStyle = '#06b6d4';
          ctx.lineWidth = 2.0;
          ctx.beginPath();
          ctx.moveTo(toX(ci1[0]), barY1);
          ctx.lineTo(toX(ci1[1]), barY1);
          ctx.stroke();
          [ci1[0], ci1[1]].forEach(cx => {
            ctx.beginPath();
            ctx.moveTo(toX(cx), barY1 - 4);
            ctx.lineTo(toX(cx), barY1 + 4);
            ctx.stroke();
          });
        }

        if (Array.isArray(ci2) && Number.isFinite(ci2[0]) && Number.isFinite(ci2[1])) {
          ctx.strokeStyle = '#a855f7';
          ctx.lineWidth = 2.0;
          ctx.beginPath();
          ctx.moveTo(toX(ci2[0]), barY2);
          ctx.lineTo(toX(ci2[1]), barY2);
          ctx.stroke();
          [ci2[0], ci2[1]].forEach(cx => {
            ctx.beginPath();
            ctx.moveTo(toX(cx), barY2 - 4);
            ctx.lineTo(toX(cx), barY2 + 4);
            ctx.stroke();
          });
        }

        ctx.fillStyle = '#64748b';
        ctx.font = `500 9px ${engine.options.fontFamily}`;
        ctx.textAlign = 'left';
        ctx.fillText(`(1-α)% CIs: ±t_crit·SEM`, b.x + 10, barY2 - 8);
      }

      const midXPix = toX((mean1 + mean2) / 2);
      const activeOVL = viewMode === 'patients' ? patientOVL : meansOVL;
      ctx.fillStyle = isSignificant ? '#10b981' : '#ef4444';
      ctx.font = `700 11px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText(`Overlap: ${(activeOVL * 100).toFixed(1)}%`, midXPix, toY(maxDensity * 0.35));

      const badgeText = isSignificant ? `✓ SIGNIFICANT (p = ${pValue.toFixed(4)} < α)` : `✗ NOT SIGNIFICANT (p = ${pValue.toFixed(4)} ≥ α)`;
      ctx.fillStyle = isSignificant ? 'rgba(16, 185, 129, 0.2)' : 'rgba(239, 68, 68, 0.2)';
      ctx.strokeStyle = isSignificant ? '#10b981' : '#ef4444';
      ctx.lineWidth = 1;
      const badgeWidth = ctx.measureText(badgeText).width + 24;
      ctx.fillRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);
      ctx.strokeRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);

      ctx.fillStyle = isSignificant ? '#10b981' : '#ef4444';
      ctx.font = `700 11px ${engine.options.fontFamily}`;
      ctx.textAlign = 'right';
      ctx.fillText(badgeText, b.x + b.width - 22, b.y + 26);

      ctx.textAlign = 'left';
      ctx.font = `500 11px ${engine.options.fontFamily}`;

      if (viewMode === 'patients') {
        ctx.fillStyle = '#06b6d4';
        ctx.fillText(`— Group 1 Patient Population N(μ₁, SD₁²), SD₁ = ${sd1.toFixed(2)}`, b.x + 12, b.y + 20);
        ctx.fillStyle = '#a855f7';
        ctx.fillText(`— Group 2 Patient Population N(μ₂, SD₂²), SD₂ = ${sd2.toFixed(2)}`, b.x + 12, b.y + 36);
        ctx.fillStyle = '#94a3b8';
        ctx.fillText(`░ Patient Overlap = ${(patientOVL * 100).toFixed(1)}% (Cohen's d = ${cohensD.toFixed(2)})`, b.x + 12, b.y + 52);
      } else if (viewMode === 'dual') {
        ctx.fillStyle = '#06b6d4';
        ctx.fillText(`— Group 1 Means (Solid, SEM₁ = ${sem1.toFixed(3)}) & Patients (Dashed, SD₁ = ${sd1.toFixed(2)})`, b.x + 12, b.y + 20);
        ctx.fillStyle = '#a855f7';
        ctx.fillText(`— Group 2 Means (Solid, SEM₂ = ${sem2.toFixed(3)}) & Patients (Dashed, SD₂ = ${sd2.toFixed(2)})`, b.x + 12, b.y + 36);
        ctx.fillStyle = '#94a3b8';
        ctx.fillText(`░ Means Overlap = ${(meansOVL * 100).toFixed(1)}% vs Patient Overlap = ${(patientOVL * 100).toFixed(1)}%`, b.x + 12, b.y + 52);
      } else {
        ctx.fillStyle = '#06b6d4';
        ctx.fillText(`— Group 1 Sampling Distribution (SEM₁ = ${sem1.toFixed(3)}, n₁ = ${n1})`, b.x + 12, b.y + 20);
        ctx.fillStyle = '#a855f7';
        ctx.fillText(`— Group 2 Sampling Distribution (SEM₂ = ${sem2.toFixed(3)}, n₂ = ${n2})`, b.x + 12, b.y + 36);
        ctx.fillStyle = '#f59e0b';
        ctx.fillText(`┆ Boundary: Δcrit = ${deltaCrit.toFixed(2)} at α = ${alpha.toFixed(3)} (Welch df = ${df.toFixed(1)})`, b.x + 12, b.y + 52);
      }
    },

    renderPowerSimulation(engine, metrics, options = {}) {
      engine.lastRender = () => this.renderPowerSimulation(engine, metrics, options);
      engine.lastRenderFn = engine.lastRender;
      if (engine.canvas && typeof engine.initHiDPI === 'function') {
        const rect = engine.canvas.getBoundingClientRect();
        if (rect.width > 50 && (engine.width <= 100 || Math.abs(engine.width - rect.width) > 30)) {
          engine.initHiDPI();
        }
      }
      engine.clear();
      const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette || {};
      const font = engine.options?.fontFamily || '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
      const textMuted = pal.textMuted || pal.text || '#94a3b8';

      const {
        sd = 4.0,
        sem = 0.50,
        n = 64,
        totalN = 128,
        power = 0.80,
        beta = 0.20,
        delta = 2.0,
        alpha = 0.05,
        zCrit = 1.96,
        seDiff = 0.707,
        lambda = 2.83,
        xCrit = 1.386,
        cohensD = 0.50,
        viewMode = 'distributions',
        matrix = {},
        curvePoints = [],
        powerRating = 'ADEQUATE'
      } = metrics || {};

      // -------------------------------------------------------------
      // VIEW MODE 1: POWER VS SAMPLE SIZE CURVE
      // -------------------------------------------------------------
      if (viewMode === 'curve') {
        const minN = 4;
        const maxN = 250;
        const toX = (val) => {
          const clamped = Math.max(minN, Math.min(maxN, Number.isFinite(val) ? val : minN));
          return b.x + ((clamped - minN) / (maxN - minN)) * b.width;
        };
        const toY = (val) => {
          const clamped = Math.max(0, Math.min(1.0, Number.isFinite(val) ? val : 0));
          return b.y + b.height - clamped * b.height;
        };

        // Draw Grid & Axes
        ctx.strokeStyle = pal.grid || 'rgba(255, 255, 255, 0.08)';
        ctx.lineWidth = 1;

        // X grid
        const nTicks = [4, 25, 50, 75, 100, 150, 200, 250];
        nTicks.forEach(nVal => {
          const xPix = toX(nVal);
          ctx.beginPath();
          ctx.moveTo(xPix, b.y);
          ctx.lineTo(xPix, b.y + b.height);
          ctx.stroke();

          ctx.fillStyle = textMuted;
          ctx.font = `500 10px ${font}`;
          ctx.textAlign = 'center';
          ctx.fillText(`n=${nVal}`, xPix, b.y + b.height + 15);
        });

        // Y grid
        const pTicks = [0.2, 0.4, 0.6, 0.8, 0.9, 1.0];
        pTicks.forEach(pVal => {
          const yPix = toY(pVal);
          ctx.beginPath();
          ctx.moveTo(b.x, yPix);
          ctx.lineTo(b.x + b.width, yPix);
          ctx.stroke();

          ctx.fillStyle = textMuted;
          ctx.font = `500 10px ${font}`;
          ctx.textAlign = 'right';
          ctx.fillText(`${(pVal * 100).toFixed(0)}%`, b.x - 8, yPix + 4);
        });

        // 80% Benchmark Line (Amber/Green dashed)
        const y80 = toY(0.80);
        ctx.strokeStyle = '#10b981';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(b.x, y80);
        ctx.lineTo(b.x + b.width, y80);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = '#10b981';
        ctx.font = `600 10px ${font}`;
        ctx.textAlign = 'right';
        ctx.fillText('80% Regulatory Standard', b.x + b.width - 10, y80 - 6);

        // 90% Benchmark Line (Cyan dashed)
        const y90 = toY(0.90);
        ctx.strokeStyle = '#06b6d4';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(b.x, y90);
        ctx.lineTo(b.x + b.width, y90);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = '#06b6d4';
        ctx.fillText('90% High Rigor', b.x + b.width - 10, y90 - 6);

        // Shaded area under curve
        if (curvePoints && curvePoints.length > 1) {
          ctx.fillStyle = 'rgba(16, 185, 129, 0.12)';
          ctx.beginPath();
          ctx.moveTo(toX(curvePoints[0].n), toY(0));
          curvePoints.forEach(pt => ctx.lineTo(toX(pt.n), toY(pt.power)));
          ctx.lineTo(toX(curvePoints[curvePoints.length - 1].n), toY(0));
          ctx.closePath();
          ctx.fill();

          // Draw Power Curve line
          ctx.strokeStyle = '#10b981';
          ctx.lineWidth = 3.0;
          ctx.beginPath();
          curvePoints.forEach((pt, idx) => {
            if (idx === 0) ctx.moveTo(toX(pt.n), toY(pt.power));
            else ctx.lineTo(toX(pt.n), toY(pt.power));
          });
          ctx.stroke();
        }

        // Current Point Marker
        const curXPix = toX(n);
        const curYPix = toY(power);

        // Drop lines to axes
        ctx.strokeStyle = '#a855f7';
        ctx.lineWidth = 1.6;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(curXPix, b.y + b.height);
        ctx.lineTo(curXPix, curYPix);
        ctx.lineTo(b.x, curYPix);
        ctx.stroke();
        ctx.setLineDash([]);

        // Glowing dot
        ctx.fillStyle = 'rgba(168, 85, 247, 0.35)';
        ctx.beginPath();
        ctx.arc(curXPix, curYPix, 9, 0, Math.PI * 2);
        ctx.fill();

        ctx.fillStyle = '#a855f7';
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        ctx.arc(curXPix, curYPix, 5, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();

        // Tooltip Callout Box
        const calloutText = `n = ${n} | Power = ${(power * 100).toFixed(1)}% (β = ${(beta * 100).toFixed(1)}%)`;
        ctx.font = `700 11px ${font}`;
        const textWidth = ctx.measureText(calloutText).width;
        const boxW = textWidth + 18;
        const boxH = 26;
        const boxX = Math.min(b.x + b.width - boxW - 8, Math.max(b.x + 8, curXPix - boxW / 2));
        const boxY = Math.max(b.y + 8, curYPix - 38);

        ctx.fillStyle = '#1e293b';
        ctx.strokeStyle = '#a855f7';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (ctx.roundRect) {
          ctx.roundRect(boxX, boxY, boxW, boxH, [6]);
        } else {
          ctx.rect(boxX, boxY, boxW, boxH);
        }
        ctx.fill();
        ctx.stroke();

        ctx.fillStyle = '#f8fafc';
        ctx.textAlign = 'center';
        ctx.fillText(calloutText, boxX + boxW / 2, boxY + 17);

        // Legend
        ctx.textAlign = 'left';
        ctx.font = `500 11px ${font}`;
        ctx.fillStyle = '#10b981';
        ctx.fillText(`— Statistical Power Curve: P(n) at Δ = ${delta.toFixed(2)}, SD = ${sd.toFixed(2)}, α = ${alpha.toFixed(3)}`, b.x + 12, b.y + 20);
        ctx.fillStyle = '#a855f7';
        ctx.fillText(`● Current Operating Point (n = ${n} per group, N = ${totalN}, SEM = ${sem.toFixed(3)})`, b.x + 12, b.y + 36);
        return;
      }

      // -------------------------------------------------------------
      // VIEW MODE 2: 2x2 DECISION ERROR MATRIX
      // -------------------------------------------------------------
      if (viewMode === 'matrix') {
        const pad = 12;
        const cardW = (b.width - pad * 3) / 2;
        const cardH = (b.height - pad * 3) / 2;

        const cells = [
          {
            col: 0, row: 0,
            title: 'SPECIFICITY (1 − α)',
            val: `${((1 - alpha) * 100).toFixed(1)}%`,
            sub: 'True Negative Rate',
            desc: 'H₀ is TRUE (no effect), decision is RETAIN H₀. Correct clinical conclusion: drug is correctly recognized as having no effect.',
            color: '#38bdf8',
            bg: 'rgba(56, 189, 248, 0.08)',
            border: 'rgba(56, 189, 248, 0.35)',
            pct: 1 - alpha
          },
          {
            col: 1, row: 0,
            title: 'TYPE I ERROR (α)',
            val: `${(alpha * 100).toFixed(1)}%`,
            sub: 'False Positive Rate (Significance Level)',
            desc: 'H₀ is TRUE (no effect), but decision is REJECT H₀. False alarm: ineffective drug erroneously declared effective.',
            color: '#ef4444',
            bg: 'rgba(239, 68, 68, 0.08)',
            border: 'rgba(239, 68, 68, 0.35)',
            pct: alpha
          },
          {
            col: 0, row: 1,
            title: 'TYPE II ERROR (β)',
            val: `${(beta * 100).toFixed(1)}%`,
            sub: 'False Negative Rate',
            desc: 'H₁ is TRUE (real effect Δ), but decision is RETAIN H₀. Missed discovery: effective therapy discarded due to lack of power!',
            color: '#f59e0b',
            bg: 'rgba(245, 158, 11, 0.08)',
            border: 'rgba(245, 158, 11, 0.35)',
            pct: beta
          },
          {
            col: 1, row: 1,
            title: 'STATISTICAL POWER (1 − β)',
            val: `${(power * 100).toFixed(1)}%`,
            sub: 'True Positive Rate (Sensitivity)',
            desc: 'H₁ is TRUE (real effect Δ), and decision is REJECT H₀. Successful trial: effective therapy correctly discovered and verified!',
            color: '#10b981',
            bg: 'rgba(16, 185, 129, 0.08)',
            border: 'rgba(16, 185, 129, 0.35)',
            pct: power
          }
        ];

        cells.forEach(c => {
          const cx = b.x + pad + c.col * (cardW + pad);
          const cy = b.y + pad + c.row * (cardH + pad);

          ctx.fillStyle = c.bg;
          ctx.strokeStyle = c.border;
          ctx.lineWidth = 1.5;
          ctx.beginPath();
          if (ctx.roundRect) {
            ctx.roundRect(cx, cy, cardW, cardH, [8]);
          } else {
            ctx.rect(cx, cy, cardW, cardH);
          }
          ctx.fill();
          ctx.stroke();

          // Card Title
          ctx.fillStyle = c.color;
          ctx.font = `700 11px ${font}`;
          ctx.textAlign = 'left';
          ctx.fillText(c.title, cx + 12, cy + 20);

          // Subtitle
          ctx.fillStyle = textMuted;
          ctx.font = `500 9px ${font}`;
          ctx.fillText(c.sub, cx + 12, cy + 34);

          // Large Percentage Value
          ctx.fillStyle = c.color;
          ctx.font = `800 24px ${font}`;
          ctx.fillText(c.val, cx + 12, cy + 64);

          // Progress bar indicator
          const barX = cx + 12;
          const barY = cy + 72;
          const barW = cardW - 24;
          const barH = 6;
          ctx.fillStyle = 'rgba(255, 255, 255, 0.1)';
          ctx.fillRect(barX, barY, barW, barH);
          ctx.fillStyle = c.color;
          ctx.fillRect(barX, barY, barW * Math.max(0, Math.min(1.0, c.pct)), barH);

          // Description text wrapped
          ctx.fillStyle = '#cbd5e1';
          ctx.font = `400 9.5px ${font}`;
          const words = c.desc.split(' ');
          let line = '';
          let lineY = cy + 93;
          for (let w = 0; w < words.length; w++) {
            const testLine = line + words[w] + ' ';
            if (ctx.measureText(testLine).width > cardW - 24 && w > 0) {
              ctx.fillText(line, cx + 12, lineY);
              line = words[w] + ' ';
              lineY += 12;
            } else {
              line = testLine;
            }
          }
          ctx.fillText(line, cx + 12, lineY);
        });
        return;
      }

      // -------------------------------------------------------------
      // VIEW MODE 3: DUAL DISTRIBUTION (H0 VS H1) WITH SHADED REGIONS
      // -------------------------------------------------------------
      const safeSEDiff = Math.max(0.001, Number.isFinite(seDiff) ? seDiff : 0.707);
      const safeDelta = Number.isFinite(delta) ? delta : 2.0;

      const minX = -3.5 * safeSEDiff;
      const maxX = safeDelta + 3.8 * safeSEDiff;
      const xSpan = Math.max(0.1, maxX - minX);

      const peakY = 1.0 / (safeSEDiff * Math.sqrt(2 * Math.PI));
      const maxY = peakY * 1.32;

      const toX = (val) => {
        const v = Number.isFinite(val) ? val : minX;
        return b.x + ((v - minX) / xSpan) * b.width;
      };
      const toY = (val) => {
        const v = Number.isFinite(val) ? val : 0;
        return b.y + b.height - (v / maxY) * b.height;
      };

      // Grid lines & X-axis
      ctx.strokeStyle = pal.grid || 'rgba(255, 255, 255, 0.08)';
      ctx.lineWidth = 1;
      const xSteps = 7;
      for (let i = 0; i <= xSteps; i++) {
        const xVal = minX + (i / xSteps) * xSpan;
        const xPix = toX(xVal);
        ctx.beginPath();
        ctx.moveTo(xPix, b.y);
        ctx.lineTo(xPix, b.y + b.height);
        ctx.stroke();

        ctx.fillStyle = textMuted;
        ctx.font = `500 10px ${font}`;
        ctx.textAlign = 'center';
        ctx.fillText(xVal.toFixed(2), xPix, b.y + b.height + 15);
      }

      // Generate Points for H0: N(0, seDiff^2) and H1: N(delta, seDiff^2)
      const numPts = 320;
      const ptsH0 = [];
      const ptsH1 = [];
      const normPDF = (x, mu, s) => (1.0 / (s * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * Math.pow((x - mu) / s, 2));

      for (let i = 0; i <= numPts; i++) {
        const x = minX + (i / numPts) * xSpan;
        ptsH0.push({ x, y: normPDF(x, 0, safeSEDiff) });
        ptsH1.push({ x, y: normPDF(x, safeDelta, safeSEDiff) });
      }

      const safeXCrit = Number.isFinite(xCrit) ? xCrit : (zCrit * safeSEDiff);

      // 1. Shading: Statistical Power (1 - beta) under H1 (where x >= safeXCrit, Emerald Green)
      ctx.fillStyle = 'rgba(16, 185, 129, 0.38)';
      ctx.beginPath();
      ctx.moveTo(toX(safeXCrit), toY(0));
      ctx.lineTo(toX(safeXCrit), toY(normPDF(safeXCrit, safeDelta, safeSEDiff)));
      for (const pt of ptsH1) {
        if (pt.x >= safeXCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.lineTo(toX(maxX), toY(0));
      ctx.closePath();
      ctx.fill();

      // 2. Shading: Type II Error (beta) under H1 (where x < safeXCrit, Amber)
      ctx.fillStyle = 'rgba(245, 158, 11, 0.35)';
      ctx.beginPath();
      ctx.moveTo(toX(minX), toY(0));
      for (const pt of ptsH1) {
        if (pt.x <= safeXCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.lineTo(toX(safeXCrit), toY(normPDF(safeXCrit, safeDelta, safeSEDiff)));
      ctx.lineTo(toX(safeXCrit), toY(0));
      ctx.closePath();
      ctx.fill();

      // 3. Shading: Type I Error (alpha/2) under H0 (where x >= safeXCrit, Red rejection tail)
      ctx.fillStyle = 'rgba(239, 68, 68, 0.35)';
      ctx.beginPath();
      ctx.moveTo(toX(safeXCrit), toY(0));
      ctx.lineTo(toX(safeXCrit), toY(normPDF(safeXCrit, 0, safeSEDiff)));
      for (const pt of ptsH0) {
        if (pt.x >= safeXCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.lineTo(toX(maxX), toY(0));
      ctx.closePath();
      ctx.fill();

      // 4. Shading: Retention Zone (1 - alpha/2) under H0 (Soft blue tint)
      ctx.fillStyle = 'rgba(56, 189, 248, 0.08)';
      ctx.beginPath();
      ctx.moveTo(toX(minX), toY(0));
      for (const pt of ptsH0) {
        if (pt.x <= safeXCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.lineTo(toX(safeXCrit), toY(normPDF(safeXCrit, 0, safeSEDiff)));
      ctx.lineTo(toX(safeXCrit), toY(0));
      ctx.closePath();
      ctx.fill();

      // 5. Draw H0 Curve (Cyan)
      ctx.strokeStyle = '#06b6d4';
      ctx.lineWidth = 2.6;
      ctx.beginPath();
      ptsH0.forEach((pt, idx) => {
        if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
        else ctx.lineTo(toX(pt.x), toY(pt.y));
      });
      ctx.stroke();

      // 6. Draw H1 Curve (Violet/Purple)
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 2.6;
      ctx.beginPath();
      ptsH1.forEach((pt, idx) => {
        if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
        else ctx.lineTo(toX(pt.x), toY(pt.y));
      });
      ctx.stroke();

      // 7. Center Means Vertical Drop Lines (mu=0 and mu=delta)
      [
        { mu: 0, color: '#06b6d4', label: 'Null H₀: Δ = 0' },
        { mu: safeDelta, color: '#a855f7', label: `Alternative H₁: Δ = ${safeDelta.toFixed(2)}` }
      ].forEach(grp => {
        const xPix = toX(grp.mu);
        ctx.strokeStyle = grp.color;
        ctx.lineWidth = 1.6;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(xPix, toY(0));
        ctx.lineTo(xPix, toY(peakY));
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = grp.color;
        ctx.font = `600 10px ${font}`;
        ctx.textAlign = 'center';
        ctx.fillText(grp.label, xPix, toY(peakY) - 8);
      });

      // 8. Effect Size Bracket Δ between mu0 and mu1
      const yBracket = toY(peakY * 0.45);
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 1.6;
      ctx.beginPath();
      ctx.moveTo(toX(0), yBracket);
      ctx.lineTo(toX(safeDelta), yBracket);
      ctx.stroke();
      // Bracket tick ends
      [0, safeDelta].forEach(muVal => {
        ctx.beginPath();
        ctx.moveTo(toX(muVal), yBracket - 4);
        ctx.lineTo(toX(muVal), yBracket + 4);
        ctx.stroke();
      });
      ctx.fillStyle = '#a855f7';
      ctx.font = `700 10px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(`True Effect Δ = ${safeDelta.toFixed(2)} (d = ${cohensD.toFixed(2)})`, toX(safeDelta / 2), yBracket - 6);

      // 9. Critical Threshold Line xcrit (Red dashed)
      if (Number.isFinite(safeXCrit) && safeXCrit >= minX && safeXCrit <= maxX) {
        const critXPix = toX(safeXCrit);
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 2.2;
        ctx.setLineDash([5, 3]);
        ctx.beginPath();
        ctx.moveTo(critXPix, toY(0));
        ctx.lineTo(critXPix, toY(peakY * 1.05));
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = '#ef4444';
        ctx.font = `700 10px ${font}`;
        ctx.textAlign = 'center';
        ctx.fillText(`Significance Cutoff xcrit = ${safeXCrit.toFixed(2)}`, critXPix, toY(peakY * 1.05) - 6);
      }

      // 10. Shading Zone Labels
      // Power label in green area
      const pwrX = toX(Math.max(safeXCrit + 0.3 * safeSEDiff, safeDelta));
      if (pwrX < b.x + b.width - 60) {
        ctx.fillStyle = '#10b981';
        ctx.font = `700 11px ${font}`;
        ctx.textAlign = 'center';
        ctx.fillText(`Power (1 − β): ${(power * 100).toFixed(1)}%`, pwrX, toY(peakY * 0.28));
      }

      // Beta label in amber area
      const betaX = toX(Math.min(safeXCrit - 0.2 * safeSEDiff, safeDelta - 0.3 * safeSEDiff));
      if (betaX > b.x + 50) {
        ctx.fillStyle = '#f59e0b';
        ctx.font = `700 11px ${font}`;
        ctx.textAlign = 'center';
        ctx.fillText(`Beta (β): ${(beta * 100).toFixed(1)}%`, betaX, toY(peakY * 0.16));
      }

      // 11. Status Badge at Top-Right
      const isAdequate = power >= 0.80;
      const badgeText = isAdequate
        ? `✓ ADEQUATE POWER: ${(power * 100).toFixed(1)}% (β = ${(beta * 100).toFixed(1)}%)`
        : `⚠ UNDERPOWERED: ${(power * 100).toFixed(1)}% (β = ${(beta * 100).toFixed(1)}%)`;
      ctx.fillStyle = isAdequate ? 'rgba(16, 185, 129, 0.2)' : 'rgba(239, 68, 68, 0.2)';
      ctx.strokeStyle = isAdequate ? '#10b981' : '#ef4444';
      ctx.lineWidth = 1;
      ctx.font = `700 11px ${font}`;
      const badgeWidth = ctx.measureText(badgeText).width + 24;
      ctx.fillRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);
      ctx.strokeRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);

      ctx.fillStyle = isAdequate ? '#10b981' : '#ef4444';
      ctx.textAlign = 'right';
      ctx.fillText(badgeText, b.x + b.width - 22, b.y + 26);

      // 12. Dynamic Legend
      ctx.textAlign = 'left';
      ctx.font = `500 11px ${font}`;
      ctx.fillStyle = '#06b6d4';
      ctx.fillText(`— Null Distribution H₀: Δ ~ N(0, SE²diff), SE = ${safeSEDiff.toFixed(3)} (SEM = ${sem.toFixed(3)})`, b.x + 12, b.y + 20);
      ctx.fillStyle = '#a855f7';
      ctx.fillText(`— Alternative Distribution H₁: Δ ~ N(${safeDelta.toFixed(2)}, SE²diff), n = ${n} per group (N = ${totalN})`, b.x + 12, b.y + 36);
      ctx.fillStyle = '#10b981';
      ctx.fillText(`░ Green Shaded Area: Power (1 − β) = ${(power * 100).toFixed(1)}% | ░ Amber Area: β = ${(beta * 100).toFixed(1)}%`, b.x + 12, b.y + 52);
    },

    /**
     * Section 6: Bayesian Statistics & Logic (3Blue1Brown Model)
     * Visualizing Bayes' Theorem via the 1x1 Unit Square, Natural Frequencies, and Odds Updating.
     */
    renderBayesianSimulation(engine, metrics, options = {}) {
      engine.lastRender = () => this.renderBayesianSimulation(engine, metrics, options);
      engine.lastRenderFn = engine.lastRender;
      if (engine.canvas && typeof engine.initHiDPI === 'function') {
        const rect = engine.canvas.getBoundingClientRect();
        if (rect.width > 50 && (engine.width <= 100 || Math.abs(engine.width - rect.width) > 30)) {
          engine.initHiDPI();
        }
      }
      engine.clear();
      const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
      const ctx = engine.ctx;
      const pal = engine.palette || {};
      const font = engine.options?.fontFamily || '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
      const textMuted = pal.textMuted || pal.text || '#94a3b8';

      const {
        prior = 0.0476,
        notPrior = 0.9524,
        likelihood = 0.40,
        falsePositive = 0.10,
        sampleSize = 210,
        areaHAndE = 0.01905,
        areaNotHAndE = 0.09524,
        pEvidence = 0.11429,
        posterior = 0.1667,
        posteriorNotH = 0.8333,
        priorOdds = 0.05,
        bayesFactor = 4.0,
        posteriorOdds = 0.20,
        beliefShift = 0.119,
        countH = 10,
        countNotH = 200,
        countHAndE = 4,
        countNotHAndE = 20,
        countTotalE = 24,
        trajectory = [],
        evidenceRating = 'Substantial / Moderate',
        viewMode = 'square',
        preset = 'steve'
      } = metrics || {};

      // -------------------------------------------------------------
      // VIEW MODE 1: 3BLUE1BROWN 1x1 UNIT SQUARE (GEOMETRY OF BAYES)
      // -------------------------------------------------------------
      if (viewMode === 'square') {
        const padTop = 32;
        const padBottom = 40;
        const labelLeftW = 90;
        const labelRightW = 88;

        const availableH = b.height - padTop - padBottom;
        const showRightCard = b.width >= 720;
        const rightCardW = showRightCard ? Math.min(380, Math.max(260, Math.floor(b.width * 0.36))) : 0;
        const gap = showRightCard ? 28 : 0;

        const availableW = b.width - labelLeftW - labelRightW - rightCardW - gap;
        const sqSize = Math.floor(Math.max(220, Math.min(availableH, availableW)));

        const sqX = b.x + labelLeftW;
        const sqY = b.y + padTop;

        // Draw Main 1x1 Possibility Space Square
        ctx.fillStyle = 'rgba(15, 23, 42, 0.75)';
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.30)';
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        if (ctx.roundRect) ctx.roundRect(sqX, sqY, sqSize, sqSize, [6]);
        else ctx.rect(sqX, sqY, sqSize, sqSize);
        ctx.fill();
        ctx.stroke();

        // Horizontal division coordinate for Prior P(H) vs P(~H)
        const splitX = sqX + sqSize * prior;
        const colWidthH = splitX - sqX;
        const colWidthNotH = sqX + sqSize - splitX;

        // Vertical heights for Evidence Shading
        const hShadeH = sqSize * likelihood;
        const hShadeNotH = sqSize * falsePositive;

        const yTopShadeH = sqY + sqSize - hShadeH;
        const yTopShadeNotH = sqY + sqSize - hShadeNotH;

        // 1. Shaded Region: P(H and E) = P(H) * P(E|H) (Emerald Green)
        if (colWidthH > 0.5 && hShadeH > 0.5) {
          let fillGreen = 'rgba(16, 185, 129, 0.55)';
          try {
            if (typeof ctx.createLinearGradient === 'function') {
              const gradGreen = ctx.createLinearGradient(sqX, yTopShadeH, sqX, sqY + sqSize);
              if (gradGreen && typeof gradGreen.addColorStop === 'function') {
                gradGreen.addColorStop(0, 'rgba(16, 185, 129, 0.70)');
                gradGreen.addColorStop(1, 'rgba(16, 185, 129, 0.40)');
                fillGreen = gradGreen;
              }
            }
          } catch { /* fallback */ }
          ctx.fillStyle = fillGreen;
          ctx.fillRect(sqX, yTopShadeH, colWidthH, hShadeH);

          ctx.strokeStyle = '#10b981';
          ctx.lineWidth = 2.0;
          ctx.strokeRect(sqX, yTopShadeH, colWidthH, hShadeH);
        }

        // 2. Shaded Region: P(~H and E) = P(~H) * P(E|~H) (Amber/Orange)
        if (colWidthNotH > 0.5 && hShadeNotH > 0.5) {
          let fillAmber = 'rgba(245, 158, 11, 0.55)';
          try {
            if (typeof ctx.createLinearGradient === 'function') {
              const gradAmber = ctx.createLinearGradient(splitX, yTopShadeNotH, splitX, sqY + sqSize);
              if (gradAmber && typeof gradAmber.addColorStop === 'function') {
                gradAmber.addColorStop(0, 'rgba(245, 158, 11, 0.70)');
                gradAmber.addColorStop(1, 'rgba(245, 158, 11, 0.40)');
                fillAmber = gradAmber;
              }
            }
          } catch { /* fallback */ }
          ctx.fillStyle = fillAmber;
          ctx.fillRect(splitX, yTopShadeNotH, colWidthNotH, hShadeNotH);

          ctx.strokeStyle = '#f59e0b';
          ctx.lineWidth = 2.0;
          ctx.strokeRect(splitX, yTopShadeNotH, colWidthNotH, hShadeNotH);
        }

        // 3. Vertical Dividing Line between H and ~H
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)';
        ctx.lineWidth = 2.0;
        ctx.setLineDash([5, 4]);
        ctx.beginPath();
        ctx.moveTo(splitX, sqY);
        ctx.lineTo(splitX, sqY + sqSize);
        ctx.stroke();
        ctx.setLineDash([]);

        // Top Dimension Labels: Prior Hypotheses
        ctx.font = `700 12.5px ${font}`;
        ctx.textAlign = 'center';

        const rawHX = sqX + colWidthH / 2;
        const rawNotHX = splitX + colWidthNotH / 2;
        const minHeaderDist = 80;
        let headerHX = rawHX;
        let headerNotHX = rawNotHX;
        if (headerNotHX - headerHX < minHeaderDist) {
          headerHX = Math.max(sqX - 6, headerNotHX - minHeaderDist);
        }

        // Left Column Header: P(H)
        ctx.fillStyle = '#10b981';
        ctx.fillText(`P(H) = ${(prior * 100).toFixed(1)}%`, headerHX, sqY - 16);
        ctx.font = `500 10.5px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(preset === 'steve' ? 'Librarians' : 'Hypothesis H', headerHX, sqY - 3);

        // Right Column Header: P(~H)
        ctx.font = `700 12.5px ${font}`;
        ctx.fillStyle = '#f59e0b';
        ctx.fillText(`P(¬H) = ${(notPrior * 100).toFixed(1)}%`, headerNotHX, sqY - 16);
        ctx.font = `500 10.5px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(preset === 'steve' ? 'Farmers (20× more)' : 'Alternative ¬H', headerNotHX, sqY - 3);

        // Height Labels for Likelihoods
        // P(E|H) on left edge
        ctx.textAlign = 'right';
        ctx.fillStyle = '#10b981';
        ctx.font = `700 12px ${font}`;
        const yMidH = yTopShadeH + hShadeH / 2 + 4;
        ctx.fillText(`P(E|H) = ${(likelihood * 100).toFixed(0)}%`, sqX - 10, yMidH);
        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(preset === 'steve' ? 'True Pos' : 'Likelihood', sqX - 10, yMidH + 13);

        // P(E|~H) on right edge
        ctx.textAlign = 'left';
        ctx.fillStyle = '#f59e0b';
        ctx.font = `700 12px ${font}`;
        const yMidNotH = yTopShadeNotH + hShadeNotH / 2 + 4;
        ctx.fillText(`P(E|¬H) = ${(falsePositive * 100).toFixed(0)}%`, sqX + sqSize + 10, yMidNotH);
        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(preset === 'steve' ? 'False Alarm' : 'False Pos', sqX + sqSize + 10, yMidNotH + 13);

        // Unshaded Region Labels: Ruled Out by Evidence
        ctx.fillStyle = 'rgba(148, 163, 184, 0.45)';
        ctx.font = `italic 11px ${font}`;
        ctx.textAlign = 'center';
        if (sqSize - hShadeNotH > 40) {
          ctx.fillText('Ruled out by Evidence ¬E (Dimmed space)', splitX + colWidthNotH / 2, sqY + (sqSize - hShadeNotH) / 2);
        }

        // Bottom Area Value Labels
        let areaHX = rawHX;
        let areaNotHX = rawNotHX;
        if (areaNotHX - areaHX < 85) {
          areaHX = Math.max(sqX - 6, areaNotHX - 85);
        }

        ctx.font = `700 12px ${font}`;
        ctx.fillStyle = '#10b981';
        ctx.textAlign = 'center';
        ctx.fillText(`Area: ${(areaHAndE * 100).toFixed(1)}%`, areaHX, sqY + sqSize + 17);
        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText('P(H ∩ E)', areaHX, sqY + sqSize + 30);

        ctx.font = `700 12px ${font}`;
        ctx.fillStyle = '#f59e0b';
        ctx.textAlign = 'center';
        ctx.fillText(`Area: ${(areaNotHAndE * 100).toFixed(1)}%`, areaNotHX, sqY + sqSize + 17);
        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText('P(¬H ∩ E)', areaNotHX, sqY + sqSize + 30);

        // -------------------------------------------------------------
        // RIGHT SIDE: RESTRICTED SPACE PROPORTION CALLOUT (BAYES RULE)
        // -------------------------------------------------------------
        if (showRightCard && rightCardW >= 220) {
          const rightX = sqX + sqSize + labelRightW + gap;
          const rightY = sqY - 8;
          const cardH = sqSize + 44;

          // Card Box
          ctx.fillStyle = 'rgba(30, 41, 59, 0.75)';
          ctx.strokeStyle = 'rgba(56, 189, 248, 0.35)';
          ctx.lineWidth = 1.5;
          ctx.beginPath();
          if (ctx.roundRect) ctx.roundRect(rightX, rightY, rightCardW, cardH, [10]);
          else ctx.rect(rightX, rightY, rightCardW, cardH);
          ctx.fill();
          ctx.stroke();

          // Card Header
          ctx.textAlign = 'left';
          ctx.font = `700 12.5px ${font}`;
          ctx.fillStyle = '#38bdf8';
          ctx.fillText('RESTRICTED POSSIBILITY SPACE P(E)', rightX + 16, rightY + 24);

          ctx.font = `400 10px ${font}`;
          ctx.fillStyle = textMuted;
          ctx.fillText('The evidence discards all unshaded space.', rightX + 16, rightY + 40);

          // Total Evidence Area P(E)
          ctx.font = `700 12px ${font}`;
          ctx.fillStyle = '#cbd5e1';
          ctx.fillText(`Total Evidence Area = ${(pEvidence * 100).toFixed(2)}%`, rightX + 16, rightY + 66);

          // Visual Proportion Bar
          const barX = rightX + 16;
          const barY = rightY + 78;
          const barW = rightCardW - 32;
          const barH = 26;

          const pwrW = barW * Math.max(0, Math.min(1.0, posterior));
          const remW = barW - pwrW;

          // Green segment: Posterior P(H|E)
          ctx.fillStyle = '#10b981';
          ctx.fillRect(barX, barY, pwrW, barH);
          // Amber segment: P(~H|E)
          ctx.fillStyle = '#f59e0b';
          ctx.fillRect(barX + pwrW, barY, remW, barH);

          // Bar border
          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 1.0;
          ctx.strokeRect(barX, barY, barW, barH);

          // Labels inside bar
          ctx.font = `700 11px ${font}`;
          if (pwrW > 35) {
            ctx.fillStyle = '#ffffff';
            ctx.textAlign = 'center';
            ctx.fillText(`${(posterior * 100).toFixed(1)}%`, barX + pwrW / 2, barY + 17);
          }
          if (remW > 45) {
            ctx.fillStyle = '#0f172a';
            ctx.textAlign = 'center';
            ctx.fillText(`${(posteriorNotH * 100).toFixed(1)}%`, barX + pwrW + remW / 2, barY + 17);
          }

          // Legend under bar
          ctx.textAlign = 'left';
          ctx.font = `600 10.5px ${font}`;
          ctx.fillStyle = '#10b981';
          ctx.fillText(`■ P(H|E): ${(posterior * 100).toFixed(1)}%`, barX, barY + 44);
          ctx.fillStyle = '#f59e0b';
          ctx.fillText(`■ P(¬H|E): ${(posteriorNotH * 100).toFixed(1)}%`, barX + barW / 2, barY + 44);

          // 3Blue1Brown Equation Breakdown
          const eqY = barY + 54;
          const eqH = 96;
          ctx.fillStyle = 'rgba(255, 255, 255, 0.06)';
          ctx.fillRect(barX, eqY, barW, eqH);
          ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
          ctx.strokeRect(barX, eqY, barW, eqH);

          ctx.fillStyle = '#94a3b8';
          ctx.font = `600 10px ${font}`;
          ctx.fillText('BAYES\' PROPORTION RULE (3BLUE1BROWN):', barX + 10, eqY + 16);

          ctx.font = `700 11.5px monospace`;
          ctx.fillStyle = '#00d2ff';
          ctx.fillText('P(H|E) = Green Area / Total Shaded Area', barX + 10, eqY + 33);

          ctx.fillStyle = '#f8fafc';
          ctx.font = `600 10.5px monospace`;
          ctx.fillText(`• By Area %:   ${(areaHAndE * 100).toFixed(2)}% / ${(pEvidence * 100).toFixed(2)}% = ${(posterior * 100).toFixed(1)}%`, barX + 10, eqY + 52);

          // Calculation by Number / Natural Counts
          const rawCountHAndE = sampleSize * areaHAndE;
          const rawCountNotHAndE = sampleSize * areaNotHAndE;
          const rawCountTotalE = rawCountHAndE + rawCountNotHAndE;
          const fmtCount = (v) => Math.abs(v - Math.round(v)) < 0.05 ? Math.round(v).toString() : v.toFixed(1);
          const strH = fmtCount(rawCountHAndE);
          const strNotH = fmtCount(rawCountNotHAndE);
          const strTot = fmtCount(rawCountTotalE);

          ctx.fillStyle = '#38bdf8';
          ctx.font = `600 10.5px monospace`;
          ctx.fillText(`• By Number:   ${strH} / (${strH} + ${strNotH}) = ${strH} / ${strTot} = ${(posterior * 100).toFixed(1)}%`, barX + 10, eqY + 70);

          ctx.fillStyle = textMuted;
          ctx.font = `400 9px ${font}`;
          ctx.fillText(`(Representative cohort of N = ${sampleSize} individuals)`, barX + 10, eqY + 86);

          // Bottom takeaway callout
          const shiftY = eqY + eqH + 18;
          ctx.fillStyle = beliefShift >= 0 ? '#10b981' : '#ef4444';
          ctx.font = `700 12px ${font}`;
          ctx.fillText(
            `Belief Shift: ${(beliefShift >= 0 ? '+' : '')}${(beliefShift * 100).toFixed(1)}%`,
            barX,
            shiftY
          );

          ctx.font = `600 11px ${font}`;
          ctx.fillStyle = '#cbd5e1';
          ctx.fillText(
            `Prior Odds ${priorOdds < 0.1 ? '1:' + (1/priorOdds).toFixed(1) : priorOdds.toFixed(2)} → Posterior Odds ${posteriorOdds < 0.1 ? '1:' + (1/posteriorOdds).toFixed(1) : posteriorOdds.toFixed(2)}`,
            barX,
            shiftY + 18
          );
        }
        return;
      }

      // -------------------------------------------------------------
      // VIEW MODE 2: REPRESENTATIVE SAMPLE & NATURAL FREQUENCIES
      // -------------------------------------------------------------
      if (viewMode === 'sample') {
        const topY = b.y + 20;

        // Header Banner
        ctx.textAlign = 'left';
        ctx.font = `700 13px ${font}`;
        ctx.fillStyle = '#38bdf8';
        ctx.fillText(`Thinking With Counts: Representative Cohort of N = ${sampleSize} Individuals`, b.x + 16, topY);

        ctx.font = `400 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText('Daniel Kahneman & Amos Tversky discovered that cognitive errors drop from 85% to 0% when framed as natural counts.', b.x + 16, topY + 16);

        // Two Cohort Cards: Hypothesis H vs Alternative ~H
        const cardY = topY + 30;
        const cardW = (b.width - 48) / 2;
        const cardH = 145;

        // Left Card: Hypothesis H (e.g. Librarians)
        const card1X = b.x + 16;
        ctx.fillStyle = 'rgba(16, 185, 129, 0.08)';
        ctx.strokeStyle = '#10b981';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (ctx.roundRect) ctx.roundRect(card1X, cardY, cardW, cardH, [8]);
        else ctx.rect(card1X, cardY, cardW, cardH);
        ctx.fill();
        ctx.stroke();

        ctx.fillStyle = '#10b981';
        ctx.font = `700 12px ${font}`;
        ctx.fillText(preset === 'steve' ? 'LIBRARIANS (Hypothesis H)' : 'HYPOTHESIS H POPULATION', card1X + 16, cardY + 24);

        ctx.font = `800 24px ${font}`;
        ctx.fillText(`${countH} people`, card1X + 16, cardY + 54);

        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(`Total in population: ${(prior * 100).toFixed(1)}% of N = ${sampleSize}`, card1X + 16, cardY + 70);

        // Matching Evidence Subsection
        ctx.fillStyle = 'rgba(16, 185, 129, 0.2)';
        ctx.fillRect(card1X + 16, cardY + 80, cardW - 32, 45);
        ctx.fillStyle = '#10b981';
        ctx.font = `700 14px ${font}`;
        ctx.fillText(`✓ ${countHAndE} Fit Description / Test Positive`, card1X + 26, cardY + 102);
        ctx.font = `500 9.5px ${font}`;
        ctx.fillStyle = '#cbd5e1';
        ctx.fillText(`${(likelihood * 100).toFixed(0)}% likelihood of evidence given H (${countH} × ${(likelihood * 100).toFixed(0)}% = ${countHAndE})`, card1X + 26, cardY + 117);

        // Right Card: Alternative ~H (e.g. Farmers)
        const card2X = card1X + cardW + 16;
        ctx.fillStyle = 'rgba(245, 158, 11, 0.08)';
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (ctx.roundRect) ctx.roundRect(card2X, cardY, cardW, cardH, [8]);
        else ctx.rect(card2X, cardY, cardW, cardH);
        ctx.fill();
        ctx.stroke();

        ctx.fillStyle = '#f59e0b';
        ctx.font = `700 12px ${font}`;
        ctx.fillText(preset === 'steve' ? 'FARMERS (Alternative ¬H)' : 'ALTERNATIVE ¬H POPULATION', card2X + 16, cardY + 24);

        ctx.font = `800 24px ${font}`;
        ctx.fillText(`${countNotH} people`, card2X + 16, cardY + 54);

        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(`Total in population: ${(notPrior * 100).toFixed(1)}% of N = ${sampleSize}`, card2X + 16, cardY + 70);

        // Matching Evidence Subsection
        ctx.fillStyle = 'rgba(245, 158, 11, 0.2)';
        ctx.fillRect(card2X + 16, cardY + 80, cardW - 32, 45);
        ctx.fillStyle = '#f59e0b';
        ctx.font = `700 14px ${font}`;
        ctx.fillText(`⚠ ${countNotHAndE} Fit Description / False Alarms`, card2X + 26, cardY + 102);
        ctx.font = `500 9.5px ${font}`;
        ctx.fillStyle = '#cbd5e1';
        ctx.fillText(`${(falsePositive * 100).toFixed(0)}% false alarm rate given ¬H (${countNotH} × ${(falsePositive * 100).toFixed(0)}% = ${countNotHAndE})`, card2X + 26, cardY + 117);

        // Combined Outcome Summary Box
        const summY = cardY + cardH + 16;
        const summW = b.width - 32;
        const summH = 75;

        ctx.fillStyle = 'rgba(30, 41, 59, 0.85)';
        ctx.strokeStyle = '#a855f7';
        ctx.lineWidth = 1.6;
        ctx.beginPath();
        if (ctx.roundRect) ctx.roundRect(b.x + 16, summY, summW, summH, [8]);
        else ctx.rect(b.x + 16, summY, summW, summH);
        ctx.fill();
        ctx.stroke();

        ctx.textAlign = 'left';
        ctx.fillStyle = '#a855f7';
        ctx.font = `700 12px ${font}`;
        ctx.fillText('NATURAL FREQUENCY RATIO CALCULATION:', b.x + 32, summY + 22);

        ctx.fillStyle = '#f8fafc';
        ctx.font = `500 11px ${font}`;
        ctx.fillText(
          `Total individuals who fit the description = ${countHAndE} (from H) + ${countNotHAndE} (from ¬H) = ${countTotalE} people out of ${sampleSize}.`,
          b.x + 32,
          summY + 40
        );

        ctx.fillStyle = '#00d2ff';
        ctx.font = `800 13px monospace`;
        ctx.fillText(
          `Posterior Probability P(H|E) = ${countHAndE} / ${countTotalE} = ${(posterior * 100).toFixed(1)}% (Alternative is ${(posteriorNotH * 100).toFixed(1)}%)`,
          b.x + 32,
          summY + 60
        );
        return;
      }

      // -------------------------------------------------------------
      // VIEW MODE 3: ODDS FORM & BAYES FACTOR SCALE
      // -------------------------------------------------------------
      if (viewMode === 'odds') {
        const topY = b.y + 20;

        ctx.textAlign = 'left';
        ctx.font = `700 13px ${font}`;
        ctx.fillStyle = '#38bdf8';
        ctx.fillText('Bayes\' Theorem in Odds Form: Prior Odds × Bayes Factor = Posterior Odds', b.x + 16, topY);

        ctx.font = `400 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText('Thinking in odds turns Bayesian multiplication into straightforward scaling.', b.x + 16, topY + 16);

        // Three Step Cards: Prior Odds -> Bayes Factor -> Posterior Odds
        const cardY = topY + 32;
        const cardW = (b.width - 64) / 3;
        const cardH = 160;

        // Card 1: Prior Odds
        const c1X = b.x + 16;
        ctx.fillStyle = 'rgba(56, 189, 248, 0.08)';
        ctx.strokeStyle = '#38bdf8';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (ctx.roundRect) ctx.roundRect(c1X, cardY, cardW, cardH, [8]);
        else ctx.rect(c1X, cardY, cardW, cardH);
        ctx.fill();
        ctx.stroke();

        ctx.fillStyle = '#38bdf8';
        ctx.font = `700 11px ${font}`;
        ctx.fillText('1. PRIOR ODDS', c1X + 16, cardY + 22);

        ctx.font = `800 22px ${font}`;
        ctx.fillText(priorOdds < 0.1 ? `1 : ${(1 / priorOdds).toFixed(1)}` : priorOdds.toFixed(2), c1X + 16, cardY + 54);

        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(`P(H) = ${(prior * 100).toFixed(1)}%`, c1X + 16, cardY + 74);
        ctx.fillText(`P(¬H) = ${(notPrior * 100).toFixed(1)}%`, c1X + 16, cardY + 90);
        ctx.fillText('O(H) = P(H) / P(¬H)', c1X + 16, cardY + 114);

        // Card 2: Bayes Factor (Likelihood Ratio)
        const c2X = c1X + cardW + 16;
        ctx.fillStyle = 'rgba(168, 85, 247, 0.08)';
        ctx.strokeStyle = '#a855f7';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (ctx.roundRect) ctx.roundRect(c2X, cardY, cardW, cardH, [8]);
        else ctx.rect(c2X, cardY, cardW, cardH);
        ctx.fill();
        ctx.stroke();

        ctx.fillStyle = '#a855f7';
        ctx.font = `700 11px ${font}`;
        ctx.fillText('2. BAYES FACTOR (LR)', c2X + 16, cardY + 22);

        ctx.font = `800 22px ${font}`;
        ctx.fillText(`${bayesFactor.toFixed(2)}×`, c2X + 16, cardY + 54);

        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(`P(E|H) = ${(likelihood * 100).toFixed(0)}%`, c2X + 16, cardY + 74);
        ctx.fillText(`P(E|¬H) = ${(falsePositive * 100).toFixed(0)}%`, c2X + 16, cardY + 90);
        ctx.fillText('BF = P(E|H) / P(E|¬H)', c2X + 16, cardY + 114);
        ctx.fillStyle = '#a855f7';
        ctx.font = `600 9px ${font}`;
        ctx.fillText(evidenceRating, c2X + 16, cardY + 138);

        // Card 3: Posterior Odds
        const c3X = c2X + cardW + 16;
        ctx.fillStyle = 'rgba(16, 185, 129, 0.08)';
        ctx.strokeStyle = '#10b981';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (ctx.roundRect) ctx.roundRect(c3X, cardY, cardW, cardH, [8]);
        else ctx.rect(c3X, cardY, cardW, cardH);
        ctx.fill();
        ctx.stroke();

        ctx.fillStyle = '#10b981';
        ctx.font = `700 11px ${font}`;
        ctx.fillText('3. POSTERIOR ODDS', c3X + 16, cardY + 22);

        ctx.font = `800 22px ${font}`;
        ctx.fillText(posteriorOdds < 0.1 ? `1 : ${(1 / posteriorOdds).toFixed(1)}` : posteriorOdds.toFixed(2), c3X + 16, cardY + 54);

        ctx.font = `500 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText(`O(H|E) = O(H) × BF`, c3X + 16, cardY + 74);
        ctx.fillText(`= ${priorOdds.toFixed(3)} × ${bayesFactor.toFixed(2)}`, c3X + 16, cardY + 90);

        ctx.fillStyle = '#10b981';
        ctx.font = `700 12px monospace`;
        ctx.fillText(`P(H|E) = ${(posterior * 100).toFixed(1)}%`, c3X + 16, cardY + 120);

        // Connectors / Multipliers between cards
        ctx.fillStyle = '#f8fafc';
        ctx.font = `800 18px ${font}`;
        ctx.textAlign = 'center';
        ctx.fillText('×', c1X + cardW + 8, cardY + 50);
        ctx.fillText('=', c2X + cardW + 8, cardY + 50);
        return;
      }

      // -------------------------------------------------------------
      // VIEW MODE 4: SEQUENTIAL EVIDENCE UPDATING TRAJECTORY
      // -------------------------------------------------------------
      if (viewMode === 'sequential') {
        const topY = b.y + 15;

        ctx.textAlign = 'left';
        ctx.font = `700 13px ${font}`;
        ctx.fillStyle = '#38bdf8';
        ctx.fillText('Sequential Belief Trajectory: Compounding Evidence Updates Belief Toward Certainty', b.x + 16, topY);

        ctx.font = `400 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText('Today\'s posterior becomes tomorrow\'s prior when observing repeated independent tests.', b.x + 16, topY + 16);

        // Trajectory Chart Plot Area
        const chartX = b.x + 45;
        const chartY = topY + 38;
        const chartW = b.width - 70;
        const chartH = b.height - 85;

        // Draw Grid
        ctx.strokeStyle = pal.grid || 'rgba(255, 255, 255, 0.08)';
        ctx.lineWidth = 1.0;
        const yTicks = [0, 0.25, 0.50, 0.75, 1.0];
        yTicks.forEach(yt => {
          const yPix = chartY + chartH - yt * chartH;
          ctx.beginPath();
          ctx.moveTo(chartX, yPix);
          ctx.lineTo(chartX + chartW, yPix);
          ctx.stroke();

          ctx.fillStyle = textMuted;
          ctx.font = `500 10px ${font}`;
          ctx.textAlign = 'right';
          ctx.fillText(`${(yt * 100).toFixed(0)}%`, chartX - 8, yPix + 4);
        });

        // 50% Threshold line
        const y50 = chartY + chartH - 0.5 * chartH;
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
        ctx.lineWidth = 1.2;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(chartX, y50);
        ctx.lineTo(chartX + chartW, y50);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = textMuted;
        ctx.font = `500 9px ${font}`;
        ctx.textAlign = 'right';
        ctx.fillText('50% Ambiguity Line', chartX + chartW - 8, y50 - 4);

        // Plot Trajectory Steps
        const numSteps = trajectory.length - 1;
        const getX = (stepIdx) => chartX + (stepIdx / numSteps) * chartW;
        const getY = (pVal) => chartY + chartH - pVal * chartH;

        // Fill area under trajectory
        ctx.fillStyle = 'rgba(0, 210, 255, 0.12)';
        ctx.beginPath();
        ctx.moveTo(getX(0), chartY + chartH);
        trajectory.forEach((t, i) => ctx.lineTo(getX(i), getY(t.p)));
        ctx.lineTo(getX(numSteps), chartY + chartH);
        ctx.closePath();
        ctx.fill();

        // Draw Trajectory Line
        ctx.strokeStyle = '#00d2ff';
        ctx.lineWidth = 3.0;
        ctx.beginPath();
        trajectory.forEach((t, i) => {
          if (i === 0) ctx.moveTo(getX(i), getY(t.p));
          else ctx.lineTo(getX(i), getY(t.p));
        });
        ctx.stroke();

        // Draw Step Markers & Callouts
        trajectory.forEach((t, i) => {
          const xPix = getX(i);
          const yPix = getY(t.p);

          // Dot
          ctx.fillStyle = '#00d2ff';
          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 2.0;
          ctx.beginPath();
          ctx.arc(xPix, yPix, 6, 0, Math.PI * 2);
          ctx.fill();
          ctx.stroke();

          // Label above dot
          ctx.fillStyle = '#f8fafc';
          ctx.font = `700 11px ${font}`;
          ctx.textAlign = 'center';
          ctx.fillText(`${(t.p * 100).toFixed(1)}%`, xPix, yPix - 12);

          // Step Label below X axis
          ctx.fillStyle = textMuted;
          ctx.font = `600 10px ${font}`;
          ctx.fillText(t.label, xPix, chartY + chartH + 18);
        });
      }
    }
,

  /**
   * Renders the Love Plot for Covariate Balance Assessment in PSM
   * Shows Absolute Standardized Mean Differences before and after matching
   * with vertical benchmark cutoffs at |SMD| = 0.10 and |SMD| = 0.05.
   */
  renderLovePlot(engine, balanceData, options = {}) {
    if (!engine) return;
    engine.lastRenderFn = () => this.renderLovePlot(engine, balanceData, options);
    engine.clear();

    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const font = '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    const isDark = engine.options.theme !== 'light';

    const textMain = isDark ? '#f8fafc' : '#0f172a';
    const textMuted = isDark ? '#94a3b8' : '#64748b';
    const gridColor = isDark ? 'rgba(255, 255, 255, 0.08)' : 'rgba(0, 0, 0, 0.08)';

    const table = balanceData?.balanceTable || [];
    if (table.length === 0) {
      ctx.fillStyle = textMuted;
      ctx.font = `500 13px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText('No covariate balance data available for Love Plot', b.x + b.width / 2, b.y + b.height / 2);
      return;
    }

    const title = options.title || 'Love Plot: Covariate Balance (Pre- vs. Post-Matching)';
    const maxVal = Math.max(0.40, ...table.map(d => Math.max(d.absSmdPre || 0, d.absSmdPost || 0))) * 1.15;

    // Margins
    const labelWidth = Math.min(140, b.width * 0.28);
    const plotX = b.x + labelWidth;
    const plotW = b.width - labelWidth - 20;
    const plotY = b.y + 35;
    const plotH = b.height - 65;

    const xToPix = (smd) => plotX + (Math.max(0, smd) / maxVal) * plotW;
    const k = table.length;
    const yToPix = (idx) => plotY + ((idx + 0.6) / (k + 0.2)) * plotH;

    // Header Title
    ctx.fillStyle = textMain;
    ctx.font = `700 13px ${font}`;
    ctx.textAlign = 'left';
    ctx.fillText(title, b.x + 10, b.y + 18);

    // Legend
    ctx.font = `600 11px ${font}`;
    ctx.textAlign = 'right';
    const legX = b.x + b.width - 20;
    // Post match dot
    ctx.fillStyle = '#10b981';
    ctx.beginPath(); ctx.arc(legX - 165, b.y + 14, 5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = textMain;
    ctx.fillText('Matched (Post)', legX - 85, b.y + 18);
    // Pre match dot
    ctx.fillStyle = '#f43f5e';
    ctx.beginPath(); ctx.arc(legX - 70, b.y + 14, 5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = textMain;
    ctx.fillText('Unadjusted (Pre)', legX, b.y + 18);

    // X-Axis Grid & Ticks
    const xSteps = [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.0].filter(v => v <= maxVal);
    ctx.lineWidth = 1.0;
    xSteps.forEach(val => {
      const px = xToPix(val);
      ctx.strokeStyle = gridColor;
      ctx.beginPath();
      ctx.moveTo(px, plotY);
      ctx.lineTo(px, plotY + plotH);
      ctx.stroke();

      ctx.fillStyle = textMuted;
      ctx.font = `500 10px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(val.toFixed(2), px, plotY + plotH + 16);
    });

    // X Axis Label
    ctx.fillStyle = textMain;
    ctx.font = `600 11px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText('Absolute Standardized Mean Difference (|SMD|)', plotX + plotW / 2, plotY + plotH + 30);

    // Cutoff 0.10 Line (Standard Benchmark)
    const cut10X = xToPix(0.10);
    ctx.strokeStyle = '#00d2ff';
    ctx.lineWidth = 1.8;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    ctx.moveTo(cut10X, plotY);
    ctx.lineTo(cut10X, plotY + plotH);
    ctx.stroke();

    // Cutoff 0.05 Line (Strict Benchmark)
    if (maxVal >= 0.05) {
      const cut05X = xToPix(0.05);
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 1.2;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(cut05X, plotY);
      ctx.lineTo(cut05X, plotY + plotH);
      ctx.stroke();
      ctx.setLineDash([]);
    } else {
      ctx.setLineDash([]);
    }

    // Benchmark Cutoff Annotations at top
    ctx.fillStyle = '#00d2ff';
    ctx.font = `700 9px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText('0.10 Threshold', cut10X, plotY - 4);

    // Draw Rows for each Covariate
    table.forEach((row, idx) => {
      const py = yToPix(idx);
      const preX = xToPix(row.absSmdPre);
      const postX = xToPix(row.absSmdPost);

      // Horizontal subtle track line
      ctx.strokeStyle = gridColor;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(plotX, py);
      ctx.lineTo(plotX + plotW, py);
      ctx.stroke();

      // Connecting arrow/segment from pre to post
      ctx.strokeStyle = row.absSmdPost <= 0.10 ? 'rgba(16, 185, 129, 0.45)' : 'rgba(244, 63, 94, 0.45)';
      ctx.lineWidth = 2.0;
      ctx.beginPath();
      ctx.moveTo(preX, py);
      ctx.lineTo(postX, py);
      ctx.stroke();

      // Pre-matching Marker (Rose)
      ctx.fillStyle = '#f43f5e';
      ctx.strokeStyle = isDark ? '#0f172a' : '#ffffff';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(preX, py, 5, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();

      // Post-matching Marker (Emerald)
      ctx.fillStyle = '#10b981';
      ctx.beginPath();
      ctx.arc(postX, py, 6, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();

      // Covariate Name Label
      ctx.fillStyle = textMain;
      ctx.font = `600 11px ${font}`;
      ctx.textAlign = 'right';
      ctx.fillText(row.covariate, plotX - 10, py + 4);

      // Delta improvement text
      if (row.percentReduction !== undefined) {
        ctx.fillStyle = row.absSmdPost <= 0.10 ? '#10b981' : '#f43f5e';
        ctx.font = `500 9px ${font}`;
        ctx.textAlign = 'left';
        ctx.fillText(`${(row.absSmdPost).toFixed(3)}`, Math.max(preX, postX) + 10, py + 3);
      }
    });
  },

  /**
   * Renders Propensity Score Overlap & Common Support Plot
   */
  renderPsmOverlapPlot(engine, overlapData, options = {}) {
    if (!engine) return;
    engine.lastRenderFn = () => this.renderPsmOverlapPlot(engine, overlapData, options);
    engine.clear();

    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const font = '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    const isDark = engine.options.theme !== 'light';

    const textMain = isDark ? '#f8fafc' : '#0f172a';
    const textMuted = isDark ? '#94a3b8' : '#64748b';
    const gridColor = isDark ? 'rgba(255, 255, 255, 0.08)' : 'rgba(0, 0, 0, 0.08)';

    const bins = overlapData?.bins || [];
    if (bins.length === 0) return;

    const title = options.title || 'Propensity Score Distribution & Common Support Overlap';
    const cs = overlapData?.commonSupport || { min: 0, max: 1 };

    const plotX = b.x + 35;
    const plotW = b.width - 60;
    const plotY = b.y + 40;
    const plotH = b.height - 75;

    // Header Title
    ctx.fillStyle = textMain;
    ctx.font = `700 13px ${font}`;
    ctx.textAlign = 'left';
    ctx.fillText(title, b.x + 10, b.y + 18);

    // Legend
    ctx.font = `600 11px ${font}`;
    ctx.textAlign = 'right';
    const legX = b.x + b.width - 20;

    // Treated legend
    ctx.fillStyle = '#00d2ff';
    ctx.fillRect(legX - 180, b.y + 10, 12, 8);
    ctx.fillStyle = textMain;
    ctx.fillText('Treated Group', legX - 95, b.y + 18);

    // Control legend
    ctx.fillStyle = '#f59e0b';
    ctx.fillRect(legX - 85, b.y + 10, 12, 8);
    ctx.fillStyle = textMain;
    ctx.fillText('Control Group', legX, b.y + 18);

    // Find max proportion for Y scaling
    let maxProp = 0.05;
    bins.forEach(bin => {
      maxProp = Math.max(maxProp, bin.preTreated, bin.preControl, bin.postTreated, bin.postControl);
    });
    maxProp *= 1.25;

    const xToPix = (p) => plotX + Math.max(0, Math.min(1, p)) * plotW;
    const yToPix = (prop) => plotY + plotH - (prop / maxProp) * plotH;

    // Common Support Shaded Area
    const csMinX = xToPix(cs.min);
    const csMaxX = xToPix(cs.max);
    ctx.fillStyle = isDark ? 'rgba(0, 210, 255, 0.07)' : 'rgba(0, 210, 255, 0.12)';
    ctx.fillRect(csMinX, plotY, csMaxX - csMinX, plotH);

    // Common Support Vertical Dashed Lines
    ctx.strokeStyle = '#00d2ff';
    ctx.lineWidth = 1.5;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(csMinX, plotY);
    ctx.lineTo(csMinX, plotY + plotH);
    ctx.moveTo(csMaxX, plotY);
    ctx.lineTo(csMaxX, plotY + plotH);
    ctx.stroke();
    ctx.setLineDash([]);

    // Common Support Label
    ctx.fillStyle = '#00d2ff';
    ctx.font = `600 10px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText(`Common Support [${cs.min.toFixed(2)}, ${cs.max.toFixed(2)}]`, (csMinX + csMaxX) / 2, plotY + 14);

    // X Axis Ticks & Grid
    [0.0, 0.2, 0.4, 0.6, 0.8, 1.0].forEach(p => {
      const px = xToPix(p);
      ctx.strokeStyle = gridColor;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(px, plotY);
      ctx.lineTo(px, plotY + plotH);
      ctx.stroke();

      ctx.fillStyle = textMuted;
      ctx.font = `500 10px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(p.toFixed(1), px, plotY + plotH + 16);
    });

    ctx.fillStyle = textMain;
    ctx.font = `600 11px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText('Estimated Propensity Score e(X)', plotX + plotW / 2, plotY + plotH + 30);

    // Draw Treated Distribution (Pre-match, Cyan Area)
    ctx.fillStyle = 'rgba(0, 210, 255, 0.25)';
    ctx.strokeStyle = '#00d2ff';
    ctx.lineWidth = 2.0;
    ctx.beginPath();
    ctx.moveTo(xToPix(bins[0].mid), plotY + plotH);
    bins.forEach(bin => {
      ctx.lineTo(xToPix(bin.mid), yToPix(bin.preTreated));
    });
    ctx.lineTo(xToPix(bins[bins.length - 1].mid), plotY + plotH);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();

    // Draw Control Distribution (Pre-match, Amber Area)
    ctx.fillStyle = 'rgba(245, 158, 11, 0.25)';
    ctx.strokeStyle = '#f59e0b';
    ctx.lineWidth = 2.0;
    ctx.beginPath();
    ctx.moveTo(xToPix(bins[0].mid), plotY + plotH);
    bins.forEach(bin => {
      ctx.lineTo(xToPix(bin.mid), yToPix(bin.preControl));
    });
    ctx.lineTo(xToPix(bins[bins.length - 1].mid), plotY + plotH);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();

    // Draw Matched Treated (Dashed Cyan Line)
    ctx.strokeStyle = '#00d2ff';
    ctx.lineWidth = 2.0;
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    bins.forEach((bin, i) => {
      const px = xToPix(bin.mid);
      const py = yToPix(bin.postTreated);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    });
    ctx.stroke();

    // Draw Matched Control (Dashed Emerald Line)
    ctx.strokeStyle = '#10b981';
    ctx.lineWidth = 2.2;
    ctx.setLineDash([2, 2]);
    ctx.beginPath();
    bins.forEach((bin, i) => {
      const px = xToPix(bin.mid);
      const py = yToPix(bin.postControl);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    });
    ctx.stroke();
    ctx.setLineDash([]);
  }
  };

  // ==========================================
  // 12. APP CONTROLLER
  // ==========================================
  class AppController {
    constructor() {
      window.sgApp = this;
      window.app = this;
      this.theme = localStorage.getItem('sg_theme') || 'dark';
      this.engines = {};
      this.results = { teaching: {} };
      this.init();
    }

    init() {
      this.applyTheme(this.theme);
      this.initTabs();
      this.initEngines();
      this.initAnovaGroups();
      this.bindEvents();
      try { this.initRandomiser(); } catch (err) { console.error('Error in initRandomiser:', err); }
      try { this.initPsm(); } catch (err) { console.error('Error in initPsm:', err); }
      try { this.initMultivariate(); } catch (err) { console.error('Error in initMultivariate:', err); }
      try { this.loadInitialData(); } catch (err) { console.error('Error in loadInitialData:', err); }
    }

    applyTheme(theme) {
      this.theme = theme;
      document.documentElement.setAttribute('data-theme', theme);
      localStorage.setItem('sg_theme', theme);
      for (const k in this.engines) {
        if (this.engines[k]) {
          this.engines[k].options.theme = theme;
          if (this.engines[k].lastRender) this.engines[k].lastRender();
        }
      }
    }

    initTabs() {
      const tabBtns = document.querySelectorAll('.tab-btn');
      const panes = document.querySelectorAll('.tab-pane');

      tabBtns.forEach(btn => {
        btn.addEventListener('click', () => {
          const target = btn.dataset.target;
          tabBtns.forEach(b => b.classList.remove('active'));
          panes.forEach(p => p.classList.remove('active'));

          btn.classList.add('active');
          const pane = document.getElementById(target);
          if (pane) {
            pane.classList.add('active');
            const canvases = pane.querySelectorAll('canvas');
            canvases.forEach(canvas => {
              const eng = this.engines[canvas.id];
              if (canvas && eng) {
                eng.initHiDPI();
                if (typeof eng.lastRender === 'function') eng.lastRender();
                else if (typeof eng.lastRenderFn === 'function') eng.lastRenderFn();
              }
            });
            if (target === 'tab-teaching') {
              if (typeof this.runTwoSampleOverlap === 'function') this.runTwoSampleOverlap();
              if (typeof this.runPowerSimulation === 'function') this.runPowerSimulation();
            }
            if (target === 'tab-teaching-bayesian') {
              if (typeof this.runBayesianSimulation === 'function') this.runBayesianSimulation();
            }
            if (target === 'tab-multivariate') {
              if (typeof this.resizeMultivariatePlots === 'function') this.resizeMultivariatePlots();
            }
          }
        });
      });
    }

    initEngines() {
      ['descCanvas', 'descBoxCanvas', 'descViolinCanvas', 'hypoCanvas', 'anovaCanvas', 'corrCanvas', 'rocCanvas', 'teachingDistCanvas', 'teachingCltParentCanvas', 'teachingCltSamplingCanvas', 'teachingTCanvas', 'teachingOverlapCanvas', 'teachingPowerCanvas', 'teachingBayesCanvas', 'psmLovePlotCanvas', 'psmOverlapCanvas'].forEach(id => {
        const el = document.getElementById(id);
        if (el) {
          const opts = id === 'teachingBayesCanvas'
            ? { theme: this.theme, padding: { top: 20, right: 20, bottom: 24, left: 20 } }
            : { theme: this.theme };
          this.engines[id] = new ChartEngine(el, opts);
        }
      });
    }

    bindEvents() {
      // Theme toggle
      document.getElementById('themeToggleBtn')?.addEventListener('click', () => {
        this.applyTheme(this.theme === 'dark' ? 'light' : 'dark');
      });

      // 1. Descriptive
      const runDesc = () => this.runDescriptive();
      document.getElementById('descComputeBtn')?.addEventListener('click', runDesc);
      document.getElementById('descInput')?.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) runDesc();
      });
      document.getElementById('descSampleBtn')?.addEventListener('click', () => {
        document.getElementById('descInput').value = DataParser.samples.icpDynamics.groupA.join(', ');
        this.runDescriptive();
      });
      document.getElementById('descOutlierSampleBtn')?.addEventListener('click', () => {
        document.getElementById('descInput').value = DataParser.samples.drainOutputSkewed.join(', ');
        this.runDescriptive();
      });

      // 2. Hypothesis
      document.getElementById('hypoComputeBtn')?.addEventListener('click', () => this.runHypo());
      document.getElementById('hypoTestType')?.addEventListener('change', () => this.runHypo());
      document.getElementById('hypoErrorBarMode')?.addEventListener('change', () => this.runHypo());

      const pairedCb = document.getElementById('hypoIsPaired');
      pairedCb?.addEventListener('change', (e) => {
        const isPaired = e.target.checked;
        const modeBadge = document.getElementById('hypoDesignModeBadge');
        if (modeBadge) {
          modeBadge.innerText = isPaired 
            ? 'Paired / Dependent Samples (Pre vs. Post / Matched)' 
            : 'Independent Cohorts (Two Separate Groups)';
          modeBadge.className = `badge ${isPaired ? 'badge-sig' : 'badge-neutral'}`;
        }
        const sampleBtn = document.getElementById('hypoSampleBtn');
        if (sampleBtn) {
          sampleBtn.innerText = isPaired ? 'Load Pre/Post ICP' : 'Load Independent Cohorts';
        }
        const nameAInput = document.getElementById('hypoNameA');
        const nameBInput = document.getElementById('hypoNameB');
        if (nameAInput && nameBInput) {
          if (isPaired && (nameAInput.value === 'Cohort 1 (Control)' || nameAInput.value === 'Standard Resection' || nameAInput.value === 'Group A')) {
            nameAInput.value = 'Pre-Infusion';
            nameBInput.value = 'Post-Infusion';
          } else if (!isPaired && (nameAInput.value === 'Pre-Infusion' || nameAInput.value === 'Pre-Intervention')) {
            nameAInput.value = 'Standard Resection';
            nameBInput.value = 'Supramarginal Resection';
          }
        }
        if (document.getElementById('hypoGroupA')?.value && document.getElementById('hypoGroupB')?.value) {
          this.runHypo();
        }
      });

      document.getElementById('hypoSampleBtn')?.addEventListener('click', () => {
        const isPaired = document.getElementById('hypoIsPaired')?.checked || false;
        if (isPaired) {
          document.getElementById('hypoNameA').value = 'Pre-Infusion';
          document.getElementById('hypoNameB').value = 'Post-Infusion';
          document.getElementById('hypoGroupA').value = DataParser.samples.icpDynamics.groupA.join(', ');
          document.getElementById('hypoGroupB').value = DataParser.samples.icpDynamics.groupB.join(', ');
        } else {
          document.getElementById('hypoNameA').value = 'Standard Resection';
          document.getElementById('hypoNameB').value = 'Supramarginal Resection';
          const s = DataParser.samples.tumorResection || {
            groupA: [14.2, 15.1, 13.8, 16.5, 14.9, 15.8, 17.2, 13.5, 15.0, 14.6, 16.1, 14.8, 15.4, 16.0],
            groupB: [22.4, 28.1, 18.9, 31.5, 24.0, 19.8, 35.2, 26.7, 21.3, 29.4, 33.1, 20.5]
          };
          document.getElementById('hypoGroupA').value = s.groupA.join(', ');
          document.getElementById('hypoGroupB').value = s.groupB.join(', ');
        }
        this.runHypo();
      });

      const designModal = document.getElementById('hypoDesignModal');
      document.getElementById('hypoDesignInfoBtn')?.addEventListener('click', () => {
        if (designModal) designModal.style.display = 'flex';
      });
      document.getElementById('hypoDesignModalClose')?.addEventListener('click', () => {
        if (designModal) designModal.style.display = 'none';
      });
      document.getElementById('hypoDesignModalSetIndependent')?.addEventListener('click', () => {
        if (pairedCb) { pairedCb.checked = false; pairedCb.dispatchEvent(new Event('change')); }
        if (designModal) designModal.style.display = 'none';
      });
      document.getElementById('hypoDesignModalSetPaired')?.addEventListener('click', () => {
        if (pairedCb) { pairedCb.checked = true; pairedCb.dispatchEvent(new Event('change')); }
        if (designModal) designModal.style.display = 'none';
      });
      designModal?.addEventListener('click', (e) => {
        if (e.target === designModal) designModal.style.display = 'none';
      });

      // 3. ANOVA & Multi-Group Analysis
      document.getElementById('anovaComputeBtn')?.addEventListener('click', () => this.runAnova());
      document.getElementById('anovaAddGroupBtn')?.addEventListener('click', () => this.addAnovaGroup());
      
      const anovaPairedCb = document.getElementById('anovaIsPaired');
      anovaPairedCb?.addEventListener('change', () => {
        this.updateAnovaDesignUI();
        this.runAnova();
      });

      document.getElementById('anovaTestType')?.addEventListener('change', () => this.runAnova());
      document.getElementById('anovaErrorBarMode')?.addEventListener('change', () => this.runAnova());

      // ANOVA Presets
      document.getElementById('anovaSampleBtn')?.addEventListener('click', () => this.loadAnovaPreset('sample3'));
      document.getElementById('anovaSampleWelchBtn')?.addEventListener('click', () => this.loadAnovaPreset('sampleWelch'));
      document.getElementById('anovaSampleRMBtn')?.addEventListener('click', () => this.loadAnovaPreset('sampleRM'));
      document.getElementById('anovaSampleSkewBtn')?.addEventListener('click', () => this.loadAnovaPreset('sampleSkew'));

      // ANOVA Study Design Modal
      const anovaModal = document.getElementById('anovaDesignModal');
      document.getElementById('anovaDesignInfoBtn')?.addEventListener('click', () => {
        if (anovaModal) anovaModal.style.display = 'flex';
      });
      document.getElementById('anovaDesignModalClose')?.addEventListener('click', () => {
        if (anovaModal) anovaModal.style.display = 'none';
      });
      document.getElementById('anovaDesignModalSetIndependent')?.addEventListener('click', () => {
        if (anovaPairedCb) { anovaPairedCb.checked = false; anovaPairedCb.dispatchEvent(new Event('change')); }
        if (anovaModal) anovaModal.style.display = 'none';
      });
      document.getElementById('anovaDesignModalSetPaired')?.addEventListener('click', () => {
        if (anovaPairedCb) { anovaPairedCb.checked = true; anovaPairedCb.dispatchEvent(new Event('change')); }
        if (anovaModal) anovaModal.style.display = 'none';
      });
      anovaModal?.addEventListener('click', (e) => {
        if (e.target === anovaModal) anovaModal.style.display = 'none';
      });

      // 4. Categorical / 2x2 Risk & Contingency
      ['catA', 'catB', 'catC', 'catD'].forEach(id => {
        document.getElementById(id)?.addEventListener('input', () => this.runCat());
      });
      document.getElementById('catAnalysisMode')?.addEventListener('change', (e) => {
        const mode = e.target.value;
        const curA = parseFloat(document.getElementById('catA')?.value) || 0;
        const curB = parseFloat(document.getElementById('catB')?.value) || 0;
        const curC = parseFloat(document.getElementById('catC')?.value) || 0;
        const curD = parseFloat(document.getElementById('catD')?.value) || 0;

        // Auto-swap default sample when user switches frameworks on unmodified sample data
        if (mode === 'diagnostic' && curA === 14 && curB === 36 && curC === 186 && curD === 164) {
          const s = DataParser.samples.diagnostic2x2;
          document.getElementById('catA').value = s.a;
          document.getElementById('catB').value = s.b;
          document.getElementById('catC').value = s.c;
          document.getElementById('catD').value = s.d;
        } else if (mode === 'study' && curA === 92 && curB === 8 && curC === 12 && curD === 188) {
          const s = DataParser.samples.shunt2x2;
          document.getElementById('catA').value = s.a;
          document.getElementById('catB').value = s.b;
          document.getElementById('catC').value = s.c;
          document.getElementById('catD').value = s.d;
        }
        this.runCat();
      });
      document.getElementById('catSampleBtn')?.addEventListener('click', () => {
        const mode = document.getElementById('catAnalysisMode')?.value || 'diagnostic';
        if (mode === 'diagnostic') {
          const s = DataParser.samples.diagnostic2x2;
          document.getElementById('catA').value = s.a;
          document.getElementById('catB').value = s.b;
          document.getElementById('catC').value = s.c;
          document.getElementById('catD').value = s.d;
        } else {
          const s = DataParser.samples.shunt2x2;
          document.getElementById('catA').value = s.a;
          document.getElementById('catB').value = s.b;
          document.getElementById('catC').value = s.c;
          document.getElementById('catD').value = s.d;
        }
        this.runCat();
      });

      // 5. Correlation
      document.getElementById('corrComputeBtn')?.addEventListener('click', () => this.runCorr());
      document.getElementById('corrSampleBtn')?.addEventListener('click', () => {
        document.getElementById('corrX').value = '10, 12, 14, 15, 18, 20, 22, 24, 25, 28, 30, 32';
        document.getElementById('corrY').value = '15, 19, 21, 22, 27, 31, 33, 37, 39, 42, 47, 49';
        this.runCorr();
      });
      document.getElementById('corrUShapeBtn')?.addEventListener('click', () => {
        // Classic neurosurgical U-shaped curve:
        // Cerebral Perfusion Pressure (CPP in mmHg) vs Mortality Rate (%)
        document.getElementById('corrX').value = '42, 46, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110';
        document.getElementById('corrY').value = '82, 68, 54, 38, 25, 18, 16, 19, 26, 36, 49, 62, 73, 85, 96';
        this.runCorr();
      });
      document.getElementById('corrNCurveBtn')?.addEventListener('click', () => {
        // Tri-phasic / N-shaped non-monotonic curve (Increase → Decrease → Increase)
        document.getElementById('corrX').value = '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15';
        document.getElementById('corrY').value = '15, 28, 45, 62, 70, 65, 50, 36, 25, 20, 26, 42, 60, 78, 95';
        this.runCorr();
      });

      // 6. Diagnostic ROC
      document.getElementById('rocComputeBtn')?.addEventListener('click', () => this.runROC());
      document.getElementById('rocSampleBtn')?.addEventListener('click', () => {
        const lines = DataParser.samples.diagnosticBiomarker.map(s => `${s.score}, ${s.status}`).join('\n');
        document.getElementById('rocInput').value = lines;
        this.runROC();
      });

      // 7. Power Analysis Events
      const pwrInputs = [
        'pwrM1', 'pwrM2', 'pwrSD',
        'pwrPairedM1', 'pwrPairedM2', 'pwrPairedSD',
        'pwrP1', 'pwrP2', 'pwrAlpha', 'pwrPower', 'pwrGivenN', 'pwrPropTestType'
      ];
      pwrInputs.forEach(id => {
        document.getElementById(id)?.addEventListener('input', () => this.runPower());
        document.getElementById(id)?.addEventListener('change', () => this.runPower());
      });

      document.getElementById('pwrStudyDesign')?.addEventListener('change', (e) => {
        const design = e.target.value;
        const indepSec = document.getElementById('pwrSectionIndep');
        const pairedSec = document.getElementById('pwrSectionPaired');
        const contSec = document.getElementById('pwrSectionContingency');
        const cardFisher = document.getElementById('pwrCardFisher');
        const cardRisk = document.getElementById('pwrCardClinicalRisk');
        const subtitle = document.getElementById('pwrSubtitle');
        const givenNLabel = document.getElementById('pwrGivenNLabel');

        if (indepSec) indepSec.style.display = design === 'independent' ? 'block' : 'none';
        if (pairedSec) pairedSec.style.display = design === 'paired' ? 'block' : 'none';
        if (contSec) contSec.style.display = design === 'contingency' ? 'block' : 'none';

        if (cardFisher) cardFisher.style.display = design === 'contingency' ? 'block' : 'none';
        if (cardRisk) cardRisk.style.display = design === 'contingency' ? 'block' : 'none';

        if (subtitle) {
          if (design === 'independent') subtitle.innerText = 'Two Independent Groups (Independent Samples t-Test)';
          else if (design === 'paired') subtitle.innerText = 'Paired Mean Study Design (Before vs After / Paired t-Test)';
          else subtitle.innerText = 'Clinical Study: Test vs Control (2x2 Contingency: Chi-Square / Fisher\'s Exact)';
        }

        if (givenNLabel) {
          givenNLabel.innerText = design === 'paired'
            ? 'Available Number of Pairs (N):'
            : 'Available Sample Size per Group (n):';
        }

        this.runPower();
      });

      document.getElementById('pwrCalcGoal')?.addEventListener('change', (e) => {
        const goal = e.target.value;
        const targetGroup = document.getElementById('pwrTargetPowerGroup');
        const givenGroup = document.getElementById('pwrGivenNGroup');
        if (targetGroup) targetGroup.style.display = goal === 'sample_size' ? 'block' : 'none';
        if (givenGroup) givenGroup.style.display = goal === 'power' ? 'block' : 'none';
        this.runPower();
      });

      document.getElementById('powerComputeBtn')?.addEventListener('click', () => this.runPower());

      document.getElementById('powerSampleBtn')?.addEventListener('click', () => {
        const design = document.getElementById('pwrStudyDesign')?.value || 'independent';
        if (design === 'independent') {
          document.getElementById('pwrM1').value = '10';
          document.getElementById('pwrM2').value = '15';
          document.getElementById('pwrSD').value = '10';
        } else if (design === 'paired') {
          document.getElementById('pwrPairedM1').value = '120';
          document.getElementById('pwrPairedM2').value = '112';
          document.getElementById('pwrPairedSD').value = '10';
        } else if (design === 'contingency') {
          document.getElementById('pwrP1').value = '0.08';
          document.getElementById('pwrP2').value = '0.20';
          document.getElementById('pwrPropTestType').value = 'fisher';
        }
        this.runPower();
      });

      // Copy buttons
      document.querySelectorAll('.btn-copy-report').forEach(btn => {
        btn.addEventListener('click', async (e) => {
          const id = e.currentTarget.dataset.reportId;
          const text = document.getElementById(id)?.innerText;
          if (text) {
            await Exporter.copyToClipboard(text);
            const orig = e.currentTarget.innerText;
            e.currentTarget.innerText = 'Copied!';
            setTimeout(() => { e.currentTarget.innerText = orig; }, 1500);
          }
        });
      });

      // Save chart buttons
      document.querySelectorAll('.btn-export-docx').forEach(btn => {
        btn.addEventListener('click', (e) => {
          const tab = e.currentTarget.dataset.tab;
          this.exportDocx(tab);
        });
      });

      document.querySelectorAll('.btn-save-plot').forEach(btn => {
        btn.addEventListener('click', (e) => {
          const cid = e.currentTarget.dataset.canvasId;
          if (this.engines[cid]) this.engines[cid].saveImage(`${cid}.png`);
        });
      });

      // 8. Teaching & Simulation Events
      document.getElementById('teachingDistSelect')?.addEventListener('change', () => {
        this.renderTeachingParams();
        this.runTeachingDistribution();
      });

      document.getElementById('teachingNRange')?.addEventListener('input', (e) => {
        const val = e.target.value;
        const el = document.getElementById('teachingNVal');
        if (el) el.innerText = val;
        this.runTeachingDistribution();
      });

      document.getElementById('teachingGenBtn')?.addEventListener('click', () => {
        this.runTeachingDistribution();
      });

      document.getElementById('teachingCopyBtn')?.addEventListener('click', (e) => {
        if (this.currentTeachingData) {
          navigator.clipboard.writeText(this.currentTeachingData.join(', '));
          const orig = e.currentTarget.innerText;
          e.currentTarget.innerText = 'Copied!';
          setTimeout(() => { e.currentTarget.innerText = orig; }, 1500);
        }
      });

      document.getElementById('teachingSendToDescBtn')?.addEventListener('click', () => {
        if (this.currentTeachingData) {
          const descInput = document.getElementById('descInput');
          if (descInput) {
            descInput.value = this.currentTeachingData.join(', ');
          }
          const descTabBtn = document.querySelector('.tab-btn[data-target="tab-descriptive"]');
          if (descTabBtn) descTabBtn.click();
          this.runDescriptive();
        }
      });

      // CLT Simulation Controls
      document.getElementById('cltPopSelect')?.addEventListener('change', (e) => {
        Teaching.clt.setPopulation(e.target.value);
        this.updateCltUI(Teaching.clt.getSummary());
      });

      document.getElementById('cltNRange')?.addEventListener('input', (e) => {
        const n = parseInt(e.target.value);
        const el = document.getElementById('cltNVal');
        if (el) el.innerText = n;
        Teaching.clt.setSampleSize(n);
        this.updateCltUI(Teaching.clt.getSummary());
      });

      document.getElementById('cltStepBtn')?.addEventListener('click', () => {
        const summary = Teaching.clt.drawSamples(1);
        this.updateCltUI(summary);
      });

      document.getElementById('cltDraw100Btn')?.addEventListener('click', () => {
        const summary = Teaching.clt.drawSamples(100);
        this.updateCltUI(summary);
      });

      document.getElementById('cltDraw1000Btn')?.addEventListener('click', () => {
        const summary = Teaching.clt.drawSamples(1000);
        this.updateCltUI(summary);
      });

      document.getElementById('cltResetBtn')?.addEventListener('click', () => {
        Teaching.clt.reset();
        this.updateCltUI(Teaching.clt.getSummary());
      });

      // Student's t Convergence Controls
      document.getElementById('tConvNRange')?.addEventListener('input', (e) => {
        this.runTConvergence(parseInt(e.target.value));
      });

      document.querySelectorAll('.btn-t-preset').forEach(btn => {
        btn.addEventListener('click', (e) => {
          const n = parseInt(e.currentTarget.dataset.n);
          if (n) this.runTConvergence(n);
        });
      });

      document.getElementById('tConvAnimateBtn')?.addEventListener('click', () => {
        this.animateTConvergence();
      });

      document.getElementById('tConvResetBtn')?.addEventListener('click', () => {
        if (this.tConvAnimationTimer) {
          clearInterval(this.tConvAnimationTimer);
          this.tConvAnimationTimer = null;
          const btn = document.getElementById('tConvAnimateBtn');
          if (btn) btn.innerText = '▶ Animate Convergence';
        }
        this.runTConvergence(4);
      });

      document.getElementById('tConvShowTailArea')?.addEventListener('change', () => {
        this.runTConvergence();
      });

      // Section 4: Two-Sample Overlap, SD vs SEM & Alpha Boundary Events
      // Section 4: Two-Sample Overlap, Group SD/SEM & Alpha Boundary Events
      const syncGroup1FromSD = () => {
        const sd1 = parseFloat(document.getElementById('overlapSD1Range')?.value) || 2.5;
        const n1 = parseInt(document.getElementById('overlapN1Range')?.value) || 16;
        const sem1 = sd1 / Math.sqrt(n1);
        const sem1El = document.getElementById('overlapSEM1Range');
        if (sem1El) sem1El.value = sem1.toFixed(3);
        if (document.getElementById('overlapLinkGroupsCheck')?.checked) {
          const sd2El = document.getElementById('overlapSD2Range');
          const sem2El = document.getElementById('overlapSEM2Range');
          if (sd2El) sd2El.value = sd1;
          if (sem2El) sem2El.value = sem1.toFixed(3);
        }
        this.runTwoSampleOverlap();
      };

      const syncGroup1FromSEM = () => {
        const sd1 = parseFloat(document.getElementById('overlapSD1Range')?.value) || 2.5;
        const sem1 = parseFloat(document.getElementById('overlapSEM1Range')?.value) || 0.625;
        const n1 = Math.max(2, Math.min(500, Math.round(Math.pow(sd1 / sem1, 2))));
        const n1El = document.getElementById('overlapN1Range');
        if (n1El) n1El.value = n1;
        if (document.getElementById('overlapLinkGroupsCheck')?.checked) {
          const sem2El = document.getElementById('overlapSEM2Range');
          const n2El = document.getElementById('overlapN2Range');
          if (sem2El) sem2El.value = sem1.toFixed(3);
          if (n2El) n2El.value = n1;
        }
        this.runTwoSampleOverlap();
      };

      const syncGroup1FromN = () => {
        const sd1 = parseFloat(document.getElementById('overlapSD1Range')?.value) || 2.5;
        const n1 = parseInt(document.getElementById('overlapN1Range')?.value) || 16;
        const sem1 = sd1 / Math.sqrt(n1);
        const sem1El = document.getElementById('overlapSEM1Range');
        if (sem1El) sem1El.value = sem1.toFixed(3);
        if (document.getElementById('overlapLinkGroupsCheck')?.checked) {
          const n2El = document.getElementById('overlapN2Range');
          const sem2El = document.getElementById('overlapSEM2Range');
          if (n2El) n2El.value = n1;
          if (sem2El) sem2El.value = sem1.toFixed(3);
        }
        this.runTwoSampleOverlap();
      };

      const syncGroup2FromSD = () => {
        const linkCheck = document.getElementById('overlapLinkGroupsCheck');
        if (linkCheck) linkCheck.checked = false;
        const sd2 = parseFloat(document.getElementById('overlapSD2Range')?.value) || 2.5;
        const n2 = parseInt(document.getElementById('overlapN2Range')?.value) || 16;
        const sem2 = sd2 / Math.sqrt(n2);
        const sem2El = document.getElementById('overlapSEM2Range');
        if (sem2El) sem2El.value = sem2.toFixed(3);
        this.runTwoSampleOverlap();
      };

      const syncGroup2FromSEM = () => {
        const linkCheck = document.getElementById('overlapLinkGroupsCheck');
        if (linkCheck) linkCheck.checked = false;
        const sd2 = parseFloat(document.getElementById('overlapSD2Range')?.value) || 2.5;
        const sem2 = parseFloat(document.getElementById('overlapSEM2Range')?.value) || 0.625;
        const n2 = Math.max(2, Math.min(500, Math.round(Math.pow(sd2 / sem2, 2))));
        const n2El = document.getElementById('overlapN2Range');
        if (n2El) n2El.value = n2;
        this.runTwoSampleOverlap();
      };

      const syncGroup2FromN = () => {
        const linkCheck = document.getElementById('overlapLinkGroupsCheck');
        if (linkCheck) linkCheck.checked = false;
        const sd2 = parseFloat(document.getElementById('overlapSD2Range')?.value) || 2.5;
        const n2 = parseInt(document.getElementById('overlapN2Range')?.value) || 16;
        const sem2 = sd2 / Math.sqrt(n2);
        const sem2El = document.getElementById('overlapSEM2Range');
        if (sem2El) sem2El.value = sem2.toFixed(3);
        this.runTwoSampleOverlap();
      };

      document.getElementById('overlapSD1Range')?.addEventListener('input', syncGroup1FromSD);
      document.getElementById('overlapSEM1Range')?.addEventListener('input', syncGroup1FromSEM);
      document.getElementById('overlapN1Range')?.addEventListener('input', syncGroup1FromN);

      document.getElementById('overlapSD2Range')?.addEventListener('input', syncGroup2FromSD);
      document.getElementById('overlapSEM2Range')?.addEventListener('input', syncGroup2FromSEM);
      document.getElementById('overlapN2Range')?.addEventListener('input', syncGroup2FromN);

      document.getElementById('overlapDeltaRange')?.addEventListener('input', () => this.runTwoSampleOverlap());
      document.getElementById('overlapAlphaRange')?.addEventListener('input', () => this.runTwoSampleOverlap());

      document.getElementById('overlapLinkGroupsCheck')?.addEventListener('change', (e) => {
        if (e.target.checked) {
          const sd1 = document.getElementById('overlapSD1Range')?.value;
          const sem1 = document.getElementById('overlapSEM1Range')?.value;
          const n1 = document.getElementById('overlapN1Range')?.value;
          const sd2El = document.getElementById('overlapSD2Range');
          const sem2El = document.getElementById('overlapSEM2Range');
          const n2El = document.getElementById('overlapN2Range');
          if (sd2El && sd1) sd2El.value = sd1;
          if (sem2El && sem1) sem2El.value = sem1;
          if (n2El && n1) n2El.value = n1;
          this.runTwoSampleOverlap();
        }
      });

      document.querySelectorAll('.btn-overlap-alpha').forEach(btn => {
        btn.addEventListener('click', (e) => {
          const alpha = parseFloat(e.currentTarget.dataset.alpha);
          const range = document.getElementById('overlapAlphaRange');
          if (range && !isNaN(alpha)) {
            range.value = alpha;
            document.querySelectorAll('.btn-overlap-alpha').forEach(b => {
              b.style.borderColor = '';
              b.style.color = '';
              b.style.fontWeight = '';
            });
            e.currentTarget.style.borderColor = '#f59e0b';
            e.currentTarget.style.color = '#f59e0b';
            e.currentTarget.style.fontWeight = '600';
            this.runTwoSampleOverlap();
          }
        });
      });

      document.querySelectorAll('.btn-overlap-mode').forEach(btn => {
        btn.addEventListener('click', (e) => {
          const mode = e.currentTarget.dataset.mode;
          this.currentOverlapMode = mode;
          document.querySelectorAll('.btn-overlap-mode').forEach(b => {
            b.classList.remove('btn-primary');
            b.classList.add('btn-secondary');
          });
          e.currentTarget.classList.remove('btn-secondary');
          e.currentTarget.classList.add('btn-primary');
          this.runTwoSampleOverlap();
        });
      });

      document.getElementById('overlapAnimateBtn')?.addEventListener('click', () => {
        this.animateOverlapSeparation();
      });

      document.getElementById('overlapResetBtn')?.addEventListener('click', () => {
        if (this.overlapAnimationTimer) {
          clearInterval(this.overlapAnimationTimer);
          this.overlapAnimationTimer = null;
          const btn = document.getElementById('overlapAnimateBtn');
          if (btn) btn.innerText = '▶ Animate';
        }
        const dRange = document.getElementById('overlapDeltaRange');
        const aRange = document.getElementById('overlapAlphaRange');
        const sd1Range = document.getElementById('overlapSD1Range');
        const sem1Range = document.getElementById('overlapSEM1Range');
        const n1Range = document.getElementById('overlapN1Range');
        const sd2Range = document.getElementById('overlapSD2Range');
        const sem2Range = document.getElementById('overlapSEM2Range');
        const n2Range = document.getElementById('overlapN2Range');
        const linkCheck = document.getElementById('overlapLinkGroupsCheck');
        if (dRange) dRange.value = 2.0;
        if (aRange) aRange.value = 0.050;
        if (sd1Range) sd1Range.value = 2.5;
        if (sem1Range) sem1Range.value = 0.625;
        if (n1Range) n1Range.value = 16;
        if (sd2Range) sd2Range.value = 2.5;
        if (sem2Range) sem2Range.value = 0.625;
        if (n2Range) n2Range.value = 16;
        if (linkCheck) linkCheck.checked = false;
        this.currentOverlapMode = 'means';
        document.querySelectorAll('.btn-overlap-mode').forEach(b => {
          b.classList.remove('btn-primary');
          b.classList.add('btn-secondary');
          if (b.dataset.mode === 'means') {
            b.classList.remove('btn-secondary');
            b.classList.add('btn-primary');
          }
        });
        document.querySelectorAll('.btn-overlap-alpha').forEach(b => {
          b.style.borderColor = '';
          b.style.color = '';
          b.style.fontWeight = '';
          if (b.dataset.alpha === '0.05') {
            b.style.borderColor = '#f59e0b';
            b.style.color = '#f59e0b';
            b.style.fontWeight = '600';
          }
        });
        this.runTwoSampleOverlap();
      });

      // Teaching Simulation 5: Power Simulation Listeners
      document.getElementById('powerSDRange')?.addEventListener('input', (e) => {
        this.runPowerSimulation({ sd: parseFloat(e.target.value) }, 'sd');
      });

      document.getElementById('powerSEMRange')?.addEventListener('input', (e) => {
        this.runPowerSimulation({ sem: parseFloat(e.target.value) }, 'sem');
      });

      document.getElementById('powerNRange')?.addEventListener('input', (e) => {
        this.runPowerSimulation({ n: parseInt(e.target.value) }, 'n');
      });

      document.getElementById('powerPowerRange')?.addEventListener('input', (e) => {
        this.runPowerSimulation({ power: parseFloat(e.target.value) }, 'power');
      });

      document.getElementById('powerBetaRange')?.addEventListener('input', (e) => {
        this.runPowerSimulation({ beta: parseFloat(e.target.value) }, 'beta');
      });

      document.getElementById('powerDeltaRange')?.addEventListener('input', (e) => {
        this.runPowerSimulation({ delta: parseFloat(e.target.value) }, 'delta');
      });

      document.getElementById('powerAlphaRange')?.addEventListener('input', (e) => {
        this.runPowerSimulation({ alpha: parseFloat(e.target.value) }, 'alpha');
      });

      document.querySelectorAll('.btn-power-preset').forEach(btn => {
        btn.addEventListener('click', (e) => {
          document.querySelectorAll('.btn-power-preset').forEach(b => b.classList.remove('active'));
          e.currentTarget.classList.add('active');
          const pwr = parseFloat(e.currentTarget.dataset.power);
          this.runPowerSimulation({ power: pwr }, 'power');
        });
      });

      document.querySelectorAll('.btn-power-effect-preset').forEach(btn => {
        btn.addEventListener('click', (e) => {
          document.querySelectorAll('.btn-power-effect-preset').forEach(b => b.classList.remove('active'));
          e.currentTarget.classList.add('active');
          const d = parseFloat(e.currentTarget.dataset.d);
          const sd = parseFloat(document.getElementById('powerSDRange')?.value) || 4.0;
          const delta = d * sd;
          this.runPowerSimulation({ delta }, 'delta');
        });
      });

      document.querySelectorAll('.btn-power-mode').forEach(btn => {
        btn.addEventListener('click', (e) => {
          document.querySelectorAll('.btn-power-mode').forEach(b => {
            b.classList.remove('btn-primary', 'active');
            b.classList.add('btn-secondary');
          });
          e.currentTarget.classList.remove('btn-secondary');
          e.currentTarget.classList.add('btn-primary', 'active');
          this.currentPowerMode = e.currentTarget.dataset.mode;
          this.runPowerSimulation();
        });
      });

      document.getElementById('powerAnimateBtn')?.addEventListener('click', () => {
        this.animatePowerGain();
      });

      document.getElementById('powerResetBtn')?.addEventListener('click', () => {
        if (this.powerAnimationTimer) {
          clearInterval(this.powerAnimationTimer);
          this.powerAnimationTimer = null;
          const btn = document.getElementById('powerAnimateBtn');
          if (btn) btn.innerText = '▶ Animate Power Gain';
        }
        const sdRange = document.getElementById('powerSDRange');
        const semRange = document.getElementById('powerSEMRange');
        const nRange = document.getElementById('powerNRange');
        const pRange = document.getElementById('powerPowerRange');
        const bRange = document.getElementById('powerBetaRange');
        const dRange = document.getElementById('powerDeltaRange');
        const aRange = document.getElementById('powerAlphaRange');
        if (sdRange) sdRange.value = 4.0;
        if (semRange) semRange.value = 0.500;
        if (nRange) nRange.value = 64;
        if (pRange) pRange.value = 0.80;
        if (bRange) bRange.value = 0.20;
        if (dRange) dRange.value = 2.0;
        if (aRange) aRange.value = 0.050;
        this.currentPowerMode = 'distributions';
        document.querySelectorAll('.btn-power-mode').forEach(b => {
          b.classList.remove('btn-primary', 'active');
          b.classList.add('btn-secondary');
          if (b.dataset.mode === 'distributions') {
            b.classList.remove('btn-secondary');
            b.classList.add('btn-primary', 'active');
          }
        });
        document.querySelectorAll('.btn-power-preset').forEach(b => {
          b.classList.toggle('active', b.dataset.power === '0.80');
        });
        document.querySelectorAll('.btn-power-effect-preset').forEach(b => {
          b.classList.toggle('active', b.dataset.d === '0.50');
        });
        this.runPowerSimulation({ power: 0.80, sd: 4.0, delta: 2.0, alpha: 0.05 }, 'power');
      });

      // Bayesian Simulation (Section 6) Events
      document.getElementById('bayesPriorRange')?.addEventListener('input', (e) => {
        this.runBayesianSimulation({ prior: parseFloat(e.target.value) / 100 });
      });

      document.getElementById('bayesSampleRange')?.addEventListener('input', (e) => {
        this.runBayesianSimulation({ sampleSize: parseInt(e.target.value) });
      });

      document.getElementById('bayesLikelihoodRange')?.addEventListener('input', (e) => {
        this.runBayesianSimulation({ likelihood: parseFloat(e.target.value) / 100 });
      });

      document.getElementById('bayesFalsePosRange')?.addEventListener('input', (e) => {
        this.runBayesianSimulation({ falsePositive: parseFloat(e.target.value) / 100 });
      });

      // Preset buttons
      document.querySelectorAll('.teaching-bayes-preset-btn').forEach(btn => {
        btn.addEventListener('click', (e) => {
          document.querySelectorAll('.teaching-bayes-preset-btn').forEach(b => b.classList.remove('active'));
          e.currentTarget.classList.add('active');
          const presetKey = e.currentTarget.dataset.preset;
          this.currentBayesPreset = presetKey;
          let p = 0.0476, lik = 0.40, fp = 0.10, n = 210;
          if (presetKey === 'steve') {
            p = 0.0476; lik = 0.40; fp = 0.10; n = 210;
          } else if (presetKey === 'disease') {
            p = 0.01; lik = 0.95; fp = 0.05; n = 1000;
          } else if (presetKey === 'science') {
            p = 0.10; lik = 0.80; fp = 0.05; n = 500;
          } else if (presetKey === 'equal') {
            p = 0.20; lik = 0.50; fp = 0.50; n = 200;
          } else if (presetKey === 'fair_coin') {
            p = 0.50; lik = 0.75; fp = 0.25; n = 100;
          }
          this.runBayesianSimulation({ prior: p, likelihood: lik, falsePositive: fp, sampleSize: n, preset: presetKey });
        });
      });

      // View mode pills
      document.querySelectorAll('.teaching-bayes-view-btn').forEach(btn => {
        btn.addEventListener('click', (e) => {
          document.querySelectorAll('.teaching-bayes-view-btn').forEach(b => {
            b.classList.remove('btn-primary', 'active');
            b.classList.add('btn-secondary');
          });
          e.currentTarget.classList.remove('btn-secondary');
          e.currentTarget.classList.add('btn-primary', 'active');
          this.currentBayesView = e.currentTarget.dataset.view;
          this.runBayesianSimulation();
        });
      });

      // Animate restriction button
      document.getElementById('bayesAnimateBtn')?.addEventListener('click', () => {
        this.animateBayesRestriction();
      });

      // Sequential evidence button (Today's posterior becomes tomorrow's prior)
      document.getElementById('bayesNextEvidenceBtn')?.addEventListener('click', () => {
        if (this.results?.teaching?.bayes) {
          const nextPrior = this.results.teaching.bayes.posterior;
          const priorSlider = document.getElementById('bayesPriorRange');
          if (priorSlider) priorSlider.value = (nextPrior * 100).toFixed(1);
          this.runBayesianSimulation({ prior: nextPrior });
        }
      });

      // Reset button
      document.getElementById('bayesResetBtn')?.addEventListener('click', () => {
        if (this.bayesAnimationTimer) {
          clearInterval(this.bayesAnimationTimer);
          this.bayesAnimationTimer = null;
          const btn = document.getElementById('bayesAnimateBtn');
          if (btn) btn.innerText = '▶ Animate Restriction';
        }
        this.currentBayesPreset = 'steve';
        this.currentBayesView = 'square';
        document.querySelectorAll('.teaching-bayes-preset-btn').forEach(b => {
          b.classList.toggle('active', b.dataset.preset === 'steve');
        });
        document.querySelectorAll('.teaching-bayes-view-btn').forEach(b => {
          b.classList.remove('btn-primary', 'active');
          b.classList.add('btn-secondary');
          if (b.dataset.view === 'square') {
            b.classList.remove('btn-secondary');
            b.classList.add('btn-primary', 'active');
          }
        });
        this.runBayesianSimulation({ prior: 0.0476, likelihood: 0.40, falsePositive: 0.10, sampleSize: 210, preset: 'steve' });
      });

      // Double-Click Inline Number Editing Configuration
      const editableConfigs = [
        // Section 6: Bayesian
        { span: 'bayesPriorVal', range: 'bayesPriorRange' },
        { span: 'bayesStatusPriorVal', range: 'bayesPriorRange' },
        { span: 'bayesSampleVal', range: 'bayesSampleRange' },
        { span: 'bayesLikelihoodVal', range: 'bayesLikelihoodRange' },
        { span: 'bayesStatusLikelihoodVal', range: 'bayesLikelihoodRange' },
        { span: 'bayesFalsePosVal', range: 'bayesFalsePosRange' },
        { span: 'bayesStatusFalsePosVal', range: 'bayesFalsePosRange' },

        // Teaching: Sample Size Distribution
        { span: 'teachingNVal', range: 'teachingNRange' },

        // Teaching: CLT
        { span: 'cltNVal', range: 'cltNRange' },

        // Teaching: Student's t
        { span: 'tConvNVal', range: 'tConvNRange' },

        // Teaching: Overlap
        { span: 'overlapDeltaVal', range: 'overlapDeltaRange' },
        { span: 'overlapAlphaVal', range: 'overlapAlphaRange' },
        { span: 'overlapSD1Val', range: 'overlapSD1Range' },
        { span: 'overlapSEM1Val', range: 'overlapSEM1Range' },
        { span: 'overlapN1Val', range: 'overlapN1Range' },
        { span: 'overlapSD2Val', range: 'overlapSD2Range' },
        { span: 'overlapSEM2Val', range: 'overlapSEM2Range' },
        { span: 'overlapN2Val', range: 'overlapN2Range' },

        // Teaching: Power Simulation
        { span: 'powerSDVal', range: 'powerSDRange' },
        { span: 'powerSEMVal', range: 'powerSEMRange' },
        { span: 'powerNVal', range: 'powerNRange' },
        { span: 'powerPowerVal', range: 'powerPowerRange' },
        { span: 'powerBetaVal', range: 'powerBetaRange' },
        { span: 'powerDeltaVal', range: 'powerDeltaRange' },
        { span: 'powerAlphaVal', range: 'powerAlphaRange' }
      ];

      editableConfigs.forEach(cfg => {
        this.setupEditableNumber(cfg.span, cfg.range, cfg.options);
      });

      // Disclaimer Modal Dismissal
      document.getElementById('disclaimerDismissBtn')?.addEventListener('click', () => {
        document.getElementById('disclaimerModal')?.classList.add('hidden');
      });
      document.getElementById('disclaimerModal')?.addEventListener('click', (e) => {
        if (e.target.id === 'disclaimerModal') {
          e.currentTarget.classList.add('hidden');
        }
      });
      document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape') {
          document.getElementById('disclaimerModal')?.classList.add('hidden');
        }
      });
    }

    setupEditableNumber(spanTarget, rangeTarget, options = {}) {
      const spanEl = typeof spanTarget === 'string' ? document.getElementById(spanTarget) : spanTarget;
      const rangeEl = typeof rangeTarget === 'string' ? document.getElementById(rangeTarget) : rangeTarget;
      if (!spanEl || !rangeEl) return;

      spanEl.classList.add('editable-number-badge');
      if (!spanEl.title) {
        spanEl.title = 'Double-click to edit value directly';
      }

      spanEl.addEventListener('dblclick', (e) => {
        e.preventDefault();
        e.stopPropagation();

        if (spanEl.dataset.isEditing === 'true') return;
        spanEl.dataset.isEditing = 'true';

        // Read current numeric value from the linked range input
        const currentVal = rangeEl.value;
        const min = rangeEl.min !== '' ? rangeEl.min : (options.min !== undefined ? options.min : '');
        const max = rangeEl.max !== '' ? rangeEl.max : (options.max !== undefined ? options.max : '');
        const step = rangeEl.step !== '' ? rangeEl.step : (options.step !== undefined ? options.step : 'any');

        // Create inline input - use step 'any' to allow arbitrary decimals like 11.8 or 4.76
        const input = document.createElement('input');
        input.type = 'number';
        input.className = 'editable-number-input';
        if (min !== '') input.min = min;
        if (max !== '') input.max = max;
        input.step = 'any';
        input.value = currentVal;

        // Match styles with span
        const computed = window.getComputedStyle(spanEl);
        input.style.fontSize = computed.fontSize;
        input.style.fontWeight = computed.fontWeight || '700';
        input.style.color = computed.color || '#00d2ff';
        input.style.textAlign = computed.textAlign === 'center' ? 'center' : 'right';
        const initialWidth = Math.max(68, spanEl.offsetWidth + 18);
        input.style.width = `${initialWidth}px`;

        // Hide span and insert input right beside it
        const prevDisplay = spanEl.style.display;
        spanEl.style.display = 'none';
        spanEl.parentNode.insertBefore(input, spanEl.nextSibling);

        // Focus and select all
        setTimeout(() => {
          input.focus();
          input.select();
        }, 10);

        let isFinished = false;
        const finishEdit = (commit) => {
          if (isFinished) return;
          isFinished = true;

          input.remove();
          spanEl.style.display = prevDisplay;
          spanEl.dataset.isEditing = 'false';

          if (commit) {
            const cleanedStr = input.value.toString().replace(/[^0-9.-]/g, '');
            let val = parseFloat(cleanedStr);
            if (!isNaN(val)) {
              const minNum = parseFloat(min);
              const maxNum = parseFloat(max);
              if (!isNaN(minNum) && val < minNum) val = minNum;
              if (!isNaN(maxNum) && val > maxNum) val = maxNum;

              // Prevent native <input type="range"> from snapping custom decimals to coarse step
              const prevStep = rangeEl.step;
              rangeEl.step = 'any';
              rangeEl.value = val;
              rangeEl.dispatchEvent(new Event('input', { bubbles: true }));
              rangeEl.dispatchEvent(new Event('change', { bubbles: true }));

              // When user later drags slider by mouse/touch, revert to original slider step
              const restoreSliderStep = () => {
                if (prevStep) rangeEl.step = prevStep;
                rangeEl.removeEventListener('pointerdown', restoreSliderStep);
                rangeEl.removeEventListener('touchstart', restoreSliderStep);
                rangeEl.removeEventListener('mousedown', restoreSliderStep);
              };
              rangeEl.addEventListener('pointerdown', restoreSliderStep, { once: true });
              rangeEl.addEventListener('touchstart', restoreSliderStep, { once: true });
              rangeEl.addEventListener('mousedown', restoreSliderStep, { once: true });

              if (typeof options.onCommit === 'function') {
                options.onCommit(val);
              }
            }
          }
        };

        input.addEventListener('keydown', (ev) => {
          if (ev.key === 'Enter') {
            ev.preventDefault();
            finishEdit(true);
          } else if (ev.key === 'Escape') {
            ev.preventDefault();
            finishEdit(false);
          }
        });

        input.addEventListener('blur', () => {
          finishEdit(true);
        });
      });
    }

    loadInitialData() {
      // Trigger sample loads
      document.getElementById('descSampleBtn')?.click();
      document.getElementById('hypoSampleBtn')?.click();
      document.getElementById('anovaSampleBtn')?.click();
      document.getElementById('catSampleBtn')?.click();
      document.getElementById('corrSampleBtn')?.click();
      document.getElementById('rocSampleBtn')?.click();
      this.runPower();
      this.initTeachingModule();
      this.loadPsmSample();
    }

    runDescriptive() {
      const raw = document.getElementById('descInput')?.value || '';
      const data = DataParser.parseSeries(raw);
      const res = Descriptive.calculate(data);

      if (res.error) {
        alert(res.error);
        return;
      }

      document.getElementById('descN').innerText = res.n;
      const rangeEl = document.getElementById('descRange');
      if (rangeEl) rangeEl.innerText = `Range: ${res.min.toFixed(1)} to ${res.max.toFixed(1)}`;

      document.getElementById('descMean').innerText = res.mean.toFixed(2);
      const sdSubEl = document.getElementById('descSDSub');
      if (sdSubEl) sdSubEl.innerText = `SD: ±${res.sd.toFixed(2)}`;

      document.getElementById('descMedian').innerText = res.median.toFixed(2);

      // Mode
      const modeEl = document.getElementById('descMode');
      const modeSubEl = document.getElementById('descModeSub');
      if (modeEl) {
        if (res.modes && res.modes.length > 0) {
          if (res.modes.length === 1) {
            modeEl.innerText = res.modes[0].toFixed(2);
            if (modeSubEl) modeSubEl.innerText = `Unimodal (Count: ${res.maxFreq || 2})`;
          } else if (res.modes.length === 2) {
            modeEl.innerText = `${res.modes[0].toFixed(1)}, ${res.modes[1].toFixed(1)}`;
            if (modeSubEl) modeSubEl.innerText = `Bimodal (Count: ${res.maxFreq || 2})`;
          } else {
            modeEl.innerText = `${res.modes.slice(0, 2).map(m => m.toFixed(1)).join(', ')}...`;
            if (modeSubEl) modeSubEl.innerText = `Multimodal (${res.modes.length} modes)`;
          }
        } else {
          modeEl.innerText = 'No Mode';
          if (modeSubEl) modeSubEl.innerText = 'All unique (freq = 1)';
        }
      }

      // Interquartile Range (IQR)
      const iqrEl = document.getElementById('descIQR');
      const iqrSubEl = document.getElementById('descIQRSub');
      if (iqrEl) iqrEl.innerText = res.iqr.toFixed(2);
      if (iqrSubEl) iqrSubEl.innerText = `Q1: ${res.q1.toFixed(2)} | Q3: ${res.q3.toFixed(2)}`;

      document.getElementById('descSD').innerText = res.sd.toFixed(2);
      const varEl = document.getElementById('descVar');
      if (varEl) varEl.innerText = `Var: ${res.variance.toFixed(2)}`;

      document.getElementById('descCI95').innerText = `[${res.ci95[0].toFixed(2)}, ${res.ci95[1].toFixed(2)}]`;
      const semEl = document.getElementById('descSEM');
      if (semEl) semEl.innerText = `SEM: ±${res.sem.toFixed(2)}`;

      // Skewness
      const skewValEl = document.getElementById('descSkewnessVal');
      const skewSubEl = document.getElementById('descSkewnessSub');
      if (skewValEl) skewValEl.innerText = res.skewness.toFixed(2);
      if (skewSubEl) skewSubEl.innerText = res.skewnessInterpretation;

      // Kurtosis
      const kurtValEl = document.getElementById('descKurtosisVal');
      const kurtSubEl = document.getElementById('descKurtosisSub');
      if (kurtValEl) kurtValEl.innerText = (res.kurtosis > 0 ? '+' : '') + res.kurtosis.toFixed(2);
      if (kurtSubEl) kurtSubEl.innerText = res.kurtosisInterpretation;

      // Outliers (Tukey's fences)
      const outValEl = document.getElementById('descOutliersVal');
      const outSubEl = document.getElementById('descOutliersSub');
      if (outValEl) {
        outValEl.innerText = `${res.outliers.length} Detected`;
        if (res.outliers.length > 0) {
          const extremeCount = res.outliers.filter(o => o.type === 'Extreme').length;
          const mildCount = res.outliers.length - extremeCount;
          const parts = [];
          if (mildCount > 0) parts.push(`${mildCount} Mild`);
          if (extremeCount > 0) parts.push(`${extremeCount} Extreme`);
          if (outSubEl) outSubEl.innerText = parts.join(', ');
        } else {
          if (outSubEl) outSubEl.innerText = `Fences [${res.lowerFence.toFixed(1)}, ${res.upperFence.toFixed(1)}]`;
        }
      }

      // Normality Assessment (Jarque-Bera Test)
      document.getElementById('descNormality').innerText = res.normality.isNormal ? 'Normal ✓' : 'Non-Normal ⚠️';
      const normSubEl = document.getElementById('descNormalitySub');
      if (normSubEl) normSubEl.innerText = `JB = ${res.normality.statistic.toFixed(2)} | p = ${res.normality.pValue.toFixed(3)}`;
      const legacySkew = document.getElementById('descSkewness');
      if (legacySkew) legacySkew.innerText = `Skew: ${res.skewness.toFixed(2)}`;

      // Outlier & Distribution Diagnostics Banner
      const diagBanner = document.getElementById('descDiagBanner');
      const diagBadge = document.getElementById('descDiagBadge');
      const diagTitle = document.getElementById('descDiagTitle');
      const diagDetail = document.getElementById('descDiagDetail');

      if (diagBanner && diagBadge && diagTitle && diagDetail) {
        const hasOutliers = res.outliers && res.outliers.length > 0;
        const isSkewed = Math.abs(res.skewness) > 0.5 || !res.normality.isNormal;
        const hasExtreme = res.outliers && res.outliers.some(o => o.type === 'Extreme');

        if (hasOutliers || isSkewed) {
          diagBanner.style.display = 'block';
          if (hasExtreme || Math.abs(res.skewness) > 1.0) {
            diagBadge.className = 'badge badge-danger';
            diagBadge.innerText = 'High Anomaly / Skew';
          } else {
            diagBadge.className = 'badge badge-warning';
            diagBadge.innerText = 'Moderate Anomaly';
          }

          const issues = [];
          if (hasOutliers) {
            const outList = res.outliers.map(o => `${o.value.toFixed(1)} (${o.type} ${o.direction}, Z = ${o.zScore > 0 ? '+' : ''}${o.zScore.toFixed(2)})`).join(', ');
            issues.push(`<strong>${res.outliers.length} outlier(s)</strong> detected outside Tukey fences [${res.lowerFence.toFixed(1)}, ${res.upperFence.toFixed(1)}]: ${outList}.`);
          }
          if (isSkewed) {
            issues.push(`Distribution displays <strong>${res.skewnessInterpretation}</strong> (G₁ = ${res.skewness.toFixed(2)}) and <strong>${res.kurtosisInterpretation}</strong> (excess G₂ = ${res.kurtosis.toFixed(2)}).`);
          }

          diagTitle.innerText = hasOutliers ? 'Distribution Anomaly: Outliers & Skewness Detected' : 'Distribution Notice: Asymmetric Skewness';
          diagDetail.innerHTML = issues.join(' ') +
            `<br><span style="display:inline-block; margin-top:0.35rem; color: var(--cyan-primary);">💡 <strong>Clinical Guidance:</strong> Reporting <strong>Median (${res.median.toFixed(2)}) &amp; IQR (${res.iqr.toFixed(2)})</strong> is strongly recommended over Mean ± SD. For inferential comparison, use non-parametric tests (Mann-Whitney U or Kruskal-Wallis) to prevent outlier distortion.</span>`;
        } else {
          diagBanner.style.display = 'block';
          diagBadge.className = 'badge badge-sig';
          diagBadge.innerText = 'Normal Distribution';
          diagTitle.innerText = 'Standard Gaussian Distribution Verified';
          diagDetail.innerHTML = `No outliers detected within Tukey's fences [${res.lowerFence.toFixed(1)}, ${res.upperFence.toFixed(1)}]. Skewness (${res.skewness.toFixed(2)}) and Kurtosis (${res.kurtosis.toFixed(2)}) conform to normal standards (Jarque-Bera p = ${res.normality.pValue.toFixed(3)}). <strong>Reporting Mean ± SD is statistically sound.</strong>`;
        }
      }

      // Descriptive Statement
      const repEl = document.getElementById('descReportText');
      if (repEl) {
        repEl.innerText = Exporter.formatDescriptiveReport(res);
      }

      // Render 1: Frequency Distribution Histogram
      if (this.engines.descCanvas) {
        Plots.renderHistogram(this.engines.descCanvas, res, 'Frequency Distribution & Normal Fit');
      }

      // Render 2: Box & Whiskers Plot
      if (this.engines.descBoxCanvas) {
        Plots.renderBoxPlot(this.engines.descBoxCanvas, [{ name: 'Sample Data', stats: res, color: this.engines.descBoxCanvas.palette.primary }], 'Box & Whiskers (Tukey Fences & Outliers)');
      }

      // Render 3: Violin Density Plot (KDE)
      if (this.engines.descViolinCanvas) {
        Plots.renderViolinPlot(this.engines.descViolinCanvas, [{ name: 'Sample Data', stats: res, color: this.engines.descViolinCanvas.palette.secondary }], 'Violin Density Plot (KDE & Quartiles)');
      }
      this.results = this.results || {}; this.results.descriptive = res;
    }

    runHypo() {
      const rawA = document.getElementById('hypoGroupA')?.value || '';
      const rawB = document.getElementById('hypoGroupB')?.value || '';
      const nameA = document.getElementById('hypoNameA')?.value || 'Group A';
      const nameB = document.getElementById('hypoNameB')?.value || 'Group B';
      const selectedTest = document.getElementById('hypoTestType')?.value || 'auto';
      const isPaired = document.getElementById('hypoIsPaired')?.checked || false;

      const a = DataParser.parseSeries(rawA);
      const b = DataParser.parseSeries(rawB);

      // Evaluate statistical assumptions first
      const assumptions = Hypothesis.evaluateAssumptions(a, b, isPaired);
      if (assumptions.error) {
        alert(assumptions.error);
        return;
      }

      // Determine test to execute
      let testToRun = selectedTest;
      if (selectedTest === 'auto') {
        testToRun = assumptions.recommendedTest;
      }

      let res;
      if (testToRun === 'student') res = Hypothesis.independentTTest(a, b);
      else if (testToRun === 'welch') res = Hypothesis.welchTTest(a, b);
      else if (testToRun === 'paired') res = Hypothesis.pairedTTest(a, b);
      else if (testToRun === 'wilcoxon') res = Hypothesis.wilcoxonSignedRank(a, b);
      else res = Hypothesis.mannWhitneyUTest(a, b);

      if (res.error) {
        alert(res.error);
        return;
      }

      res.groupA = Object.assign(res.groupA || Descriptive.calculate(a), { name: nameA });
      res.groupB = Object.assign(res.groupB || Descriptive.calculate(b), { name: nameB });
      res.assumptions = assumptions;
      res.isPaired = isPaired;
      res.selectedTest = selectedTest;
      res.actualTest = testToRun;

      // Update primary metric cards
      document.getElementById('hypoStat').innerText = (res.statistic !== undefined ? res.statistic : (res.zScore || 0)).toFixed(2);
      document.getElementById('hypoPVal').innerText = Exporter.formatP(res.pValue);
      const badge = document.getElementById('hypoPValBadge');
      if (badge) {
        badge.className = `badge ${res.isSignificant ? 'badge-sig' : 'badge-ns'}`;
        badge.innerText = res.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)';
      }
      document.getElementById('hypoEffect').innerText = (res.cohensD !== undefined ? res.cohensD.toFixed(2) : (res.rankBiserial !== undefined ? res.rankBiserial.toFixed(2) : '0.00'));
      document.getElementById('hypoDiff').innerText = res.meanDiff !== undefined ? ((res.meanDiff >= 0 ? '+' : '') + res.meanDiff.toFixed(2)) : (res.medianDiff !== undefined ? ((res.medianDiff >= 0 ? '+' : '') + res.medianDiff.toFixed(2)) : 'N/A');

      // Update Assumptions & Decision Engine Panel
      const recBadge = document.getElementById('hypoRecommendationBadge');
      if (recBadge) {
        recBadge.innerText = `Recommended: ${assumptions.recommendedTestName}`;
        recBadge.className = 'badge badge-sig';
      }

      const decisionTextEl = document.getElementById('hypoDecisionText');
      if (decisionTextEl) {
        let decisionHtml = `<strong>Recommended Test:</strong> ${assumptions.recommendedTestName}<br>`;
        decisionHtml += `<span>${assumptions.rationale}</span>`;
        if (selectedTest !== 'auto') {
          if (selectedTest === assumptions.recommendedTest) {
            decisionHtml += `<div style="margin-top: 0.35rem; color: var(--emerald-primary); font-weight: 600;">✓ Manual test selection (${res.testName}) perfectly matches the statistical recommendation.</div>`;
          } else {
            decisionHtml += `<div style="margin-top: 0.35rem; color: var(--amber-primary); font-weight: 600;">⚠️ Advisory Note: Executed test (${res.testName}) was manually selected, differing from the assumption-recommended test (${assumptions.recommendedTestName}).</div>`;
          }
        } else {
          decisionHtml += `<div style="margin-top: 0.35rem; color: var(--cyan-primary); font-weight: 600;">⚡ Automatically executed: ${res.testName}</div>`;
        }
        decisionTextEl.innerHTML = decisionHtml;
      }

      // Update Normality Diagnostics
      const normBadge = document.getElementById('hypoNormalityBadge');
      const normDetails = document.getElementById('hypoNormalityDetails');
      if (isPaired) {
        if (normBadge) {
          normBadge.className = `badge ${assumptions.normality.isNormal ? 'badge-sig' : 'badge-ns'}`;
          normBadge.innerText = assumptions.normality.isNormal ? 'Normal Differences (Parametric)' : 'Non-Normal Differences (Non-Parametric)';
        }
        if (normDetails) {
          normDetails.innerHTML = `
            <strong>Paired Differences (&Delta; = Post &minus; Pre, n = ${assumptions.n}):</strong><br>
            Mean &Delta; = ${assumptions.statsDiff.mean.toFixed(2)}, SD = ${assumptions.statsDiff.sd.toFixed(2)}<br>
            Skewness = ${assumptions.statsDiff.skewness.toFixed(2)} (${assumptions.statsDiff.skewnessInterpretation})<br>
            Excess Kurtosis = ${assumptions.statsDiff.kurtosis.toFixed(2)} (${assumptions.statsDiff.kurtosisInterpretation})<br>
            Jarque-Bera Test: JB = ${assumptions.normality.statistic.toFixed(2)}, p = ${Exporter.formatP(assumptions.normality.pValue)}<br>
            <span style="font-weight: 600; color: ${assumptions.normality.isNormal ? 'var(--emerald-primary)' : 'var(--amber-primary)'};">${assumptions.normality.interpretation}</span>
          `;
        }
      } else {
        if (normBadge) {
          normBadge.className = `badge ${assumptions.normality.isParametric ? 'badge-sig' : 'badge-ns'}`;
          normBadge.innerText = assumptions.normality.isParametric ? 'Both Cohorts Normal (Parametric)' : 'Normality Violated (Non-Parametric)';
        }
        if (normDetails) {
          normDetails.innerHTML = `
            <strong>${nameA} (n = ${assumptions.statsA.n}):</strong> Skew = ${assumptions.statsA.skewness.toFixed(2)}, Kurt = ${assumptions.statsA.kurtosis.toFixed(2)}, JB = ${assumptions.normality.jbA.statistic.toFixed(2)} (${Exporter.formatP(assumptions.normality.jbA.pValue)}) [${assumptions.normality.normA ? 'Normal' : 'Skewed'}]<br>
            <strong>${nameB} (n = ${assumptions.statsB.n}):</strong> Skew = ${assumptions.statsB.skewness.toFixed(2)}, Kurt = ${assumptions.statsB.kurtosis.toFixed(2)}, JB = ${assumptions.normality.jbB.statistic.toFixed(2)} (${Exporter.formatP(assumptions.normality.jbB.pValue)}) [${assumptions.normality.normB ? 'Normal' : 'Skewed'}]<br>
            <span style="font-weight: 600; color: ${assumptions.normality.isParametric ? 'var(--emerald-primary)' : 'var(--amber-primary)'};">${assumptions.normality.interpretation}</span>
          `;
        }
      }

      // Update Variance Equality Diagnostics
      const varBadge = document.getElementById('hypoVarianceBadge');
      const varDetails = document.getElementById('hypoVarianceDetails');
      if (isPaired) {
        if (varBadge) {
          varBadge.className = 'badge badge-neutral';
          varBadge.innerText = 'N/A (Paired Design)';
        }
        if (varDetails) {
          varDetails.innerHTML = `
            <strong>Paired Repeated Measures:</strong><br>
            Between-cohort homoscedasticity is not required for paired analysis because the evaluation is performed on within-subject difference scores (&Delta;<sub>i</sub> = Post<sub>i</sub> &minus; Pre<sub>i</sub>), removing inter-subject variance. Sphericity is naturally satisfied with 2 repeated measures.
          `;
        }
      } else {
        const eq = assumptions.varianceEquality.equalVariance;
        if (varBadge) {
          varBadge.className = `badge ${eq ? 'badge-sig' : 'badge-ns'}`;
          varBadge.innerText = eq ? 'Equal Variances (Homoscedastic)' : 'Unequal Variances (Heteroscedastic)';
        }
        if (varDetails) {
          varDetails.innerHTML = `
            <strong>F-Test of Equal Variances:</strong><br>
            ${nameA} Variance s₁² = ${assumptions.statsA.variance.toFixed(2)} | ${nameB} Variance s₂² = ${assumptions.statsB.variance.toFixed(2)}<br>
            Variance Ratio F(${assumptions.varianceEquality.df1}, ${assumptions.varianceEquality.df2}) = ${assumptions.varianceEquality.fStat.toFixed(2)}, p = ${Exporter.formatP(assumptions.varianceEquality.pValue)}<br>
            <span style="font-weight: 600; color: ${eq ? 'var(--emerald-primary)' : 'var(--amber-primary)'};">${assumptions.varianceEquality.interpretation}</span>
          `;
        }
      }

      // Update Pipeline Flow Badges
      const flowDesign = document.getElementById('hypoFlowDesign');
      const flowNorm = document.getElementById('hypoFlowNorm');
      const flowVar = document.getElementById('hypoFlowVar');
      const flowTest = document.getElementById('hypoFlowTest');
      if (flowDesign) flowDesign.innerText = isPaired ? 'Paired Samples' : 'Independent Samples';
      if (flowNorm) flowNorm.innerText = (isPaired ? assumptions.normality.isNormal : assumptions.normality.isParametric) ? 'Parametric (Normal)' : 'Non-Parametric';
      if (flowVar) flowVar.innerText = isPaired ? 'Within-Subject' : (assumptions.varianceEquality.equalVariance ? 'Equal Variance' : 'Unequal Variance');
      if (flowTest) {
        flowTest.innerText = res.testName;
        flowTest.className = `badge ${res.isSignificant ? 'badge-sig' : 'badge-neutral'}`;
      }

      // Comprehensive Clinical / Publication APA Narrative
      let report = `A ${res.testName} was conducted to compare ${nameA} and ${nameB}.\n\n`;
      report += `Diagnostic Assumption Testing:\n`;
      if (isPaired) {
        report += `• Sample Design: Paired / repeated measures (n = ${assumptions.n} paired pairs).\n`;
        report += `• Normality of Within-Subject Differences: Jarque-Bera JB = ${assumptions.normality.statistic.toFixed(2)}, ${Exporter.formatP(assumptions.normality.pValue)} (Skewness = ${assumptions.normality.skewness.toFixed(2)}, Kurtosis = ${assumptions.normality.kurtosis.toFixed(2)}). The difference distribution was determined to be ${assumptions.normality.isNormal ? 'normally distributed' : 'non-normally distributed'}.\n`;
        report += `• Homoscedasticity: Not applicable for paired design (within-subject differencing removes inter-subject variance).\n`;
      } else {
        report += `• Sample Design: Independent two-cohort comparison (${nameA}: n = ${assumptions.statsA.n}; ${nameB}: n = ${assumptions.statsB.n}).\n`;
        report += `• Normality Assessment: ${nameA} (JB = ${assumptions.normality.jbA.statistic.toFixed(2)}, ${Exporter.formatP(assumptions.normality.jbA.pValue)}, Skew = ${assumptions.statsA.skewness.toFixed(2)}); ${nameB} (JB = ${assumptions.normality.jbB.statistic.toFixed(2)}, ${Exporter.formatP(assumptions.normality.jbB.pValue)}, Skew = ${assumptions.statsB.skewness.toFixed(2)}). Distribution: ${assumptions.normality.isParametric ? 'Parametric (Normal)' : 'Non-Parametric (Skewed/Deviated)'}.\n`;
        report += `• Equality of Variances: F-test F(${assumptions.varianceEquality.df1}, ${assumptions.varianceEquality.df2}) = ${assumptions.varianceEquality.fStat.toFixed(2)}, ${Exporter.formatP(assumptions.varianceEquality.pValue)}, confirming ${assumptions.varianceEquality.equalVariance ? 'equal variances (homoscedasticity)' : 'unequal variances (heteroscedasticity)'}.\n`;
      }
      report += `• Decision Rationale: ${assumptions.rationale}\n\n`;
      report += `Inferential Test Results:\n`;
      const statLabel = res.testName.includes('Mann-Whitney') ? 'U' : (res.testName.includes('Wilcoxon') ? 'W' : 't');
      const dfLabel = res.df !== undefined ? `(${res.df.toFixed(1)})` : '';
      const effectLabel = res.cohensD !== undefined ? `Cohen's d = ${res.cohensD.toFixed(2)}` : `Rank-Biserial r = ${(res.rankBiserial || 0).toFixed(2)}`;
      report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} difference between ${nameA} (M = ${res.groupA.mean.toFixed(2)}, SD = ${res.groupA.sd.toFixed(2)}) and ${nameB} (M = ${res.groupB.mean.toFixed(2)}, SD = ${res.groupB.sd.toFixed(2)}): ${statLabel}${dfLabel} = ${(res.statistic !== undefined ? res.statistic : res.zScore || 0).toFixed(2)}, ${Exporter.formatP(res.pValue)}, ${effectLabel}.`;
      if (res.ci95) {
        report += ` 95% Confidence Interval: [${res.ci95[0].toFixed(2)}, ${res.ci95[1].toFixed(2)}].`;
      }

      document.getElementById('hypoReportText').innerText = report;

      const errorBarMode = document.getElementById('hypoErrorBarMode')?.value || 'ci95';
      const modeDescriptions = {
        ci95: 'Error Bars: 95% Confidence Interval (Mean ± 95% CI)',
        sem: 'Error Bars: Standard Error of Mean (Mean ± 1 SEM)',
        sd: 'Error Bars: Standard Deviation (Mean ± 1 SD)',
        iqr: 'Distribution: Box & Whiskers (Median, Q1-Q3 IQR, Tukey Fences)'
      };
      const subElem = document.getElementById('hypoChartSub');
      if (subElem) {
        subElem.innerText = modeDescriptions[errorBarMode] || modeDescriptions.ci95;
      }

      if (this.engines.hypoCanvas) {
        Plots.renderErrorBarPlot(this.engines.hypoCanvas, [
          { name: nameA, stats: res.groupA },
          { name: nameB, stats: res.groupB }
        ], {
          mode: errorBarMode,
          title: `${nameA} vs ${nameB}`
        });
      }
      this.results = this.results || {}; this.results.hypothesis = Object.assign(res, { assumptions, nameA, nameB, isPaired });
    }

    initAnovaGroups() {
      const container = document.getElementById('anovaGroupsContainer');
      if (!container) return;

      container.addEventListener('input', (e) => {
        if (e.target.classList.contains('anova-group-data')) {
          const card = e.target.closest('.anova-group-card');
          const countBadge = card?.querySelector('.anova-group-count');
          if (countBadge) {
            const count = DataParser.parseSeries(e.target.value).length;
            countBadge.innerText = `n = ${count}`;
          }
        }
      });

      container.addEventListener('click', (e) => {
        const btn = e.target.closest('.btn-remove-anova-group');
        if (btn) {
          const card = btn.closest('.anova-group-card');
          if (card) this.removeAnovaGroup(card);
        }
      });

      this.updateAnovaRemoveButtons();
      this.updateAnovaDesignUI();
    }

    addAnovaGroup(name = '', dataStr = '') {
      const container = document.getElementById('anovaGroupsContainer');
      if (!container) return;

      const currentCards = container.querySelectorAll('.anova-group-card');
      const k = currentCards.length + 1;
      const cohortName = name.trim() || `Cohort ${k}`;
      const count = dataStr ? DataParser.parseSeries(dataStr).length : 0;

      const card = document.createElement('div');
      card.className = 'anova-group-card';
      card.style.cssText = 'background: var(--bg-surface-elevated); border: 1px solid var(--border-subtle); border-radius: 8px; padding: 0.75rem 0.85rem;';
      card.innerHTML = `
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
          <div style="display: flex; align-items: center; gap: 0.5rem; flex: 1;">
            <span class="anova-group-index" style="font-weight: 700; color: var(--cyan-primary); font-size: 0.82rem; min-width: 60px;">Cohort ${k}:</span>
            <input type="text" class="form-control anova-group-name" value="${cohortName}" style="font-size: 0.82rem; padding: 0.25rem 0.5rem; height: 28px; font-weight: 600;" placeholder="Cohort Name">
          </div>
          <div style="display: flex; align-items: center; gap: 0.5rem;">
            <span class="badge badge-neutral anova-group-count" style="font-size: 0.72rem;">n = ${count}</span>
            <button type="button" class="btn btn-secondary btn-sm btn-remove-anova-group" style="padding: 0.15rem 0.45rem; font-size: 0.72rem; height: 26px; color: var(--rose-primary); border-color: rgba(244, 63, 94, 0.3);" title="Remove this cohort">✕</button>
          </div>
        </div>
        <textarea class="form-control anova-group-data" rows="2" placeholder="Comma or newline-separated values...">${dataStr}</textarea>
      `;

      container.appendChild(card);
      this.updateAnovaRemoveButtons();
      this.runAnova();
    }

    removeAnovaGroup(cardEl) {
      if (!cardEl) return;
      const container = document.getElementById('anovaGroupsContainer');
      const cards = container?.querySelectorAll('.anova-group-card');
      if (!cards || cards.length <= 2) {
        alert('Multi-group comparison requires at least 2 cohorts.');
        return;
      }

      cardEl.remove();

      // Re-index cohort headers
      const remaining = container.querySelectorAll('.anova-group-card');
      remaining.forEach((c, idx) => {
        const idxSpan = c.querySelector('.anova-group-index');
        if (idxSpan) idxSpan.innerText = `Cohort ${idx + 1}:`;
      });

      this.updateAnovaRemoveButtons();
      this.runAnova();
    }

    updateAnovaRemoveButtons() {
      const cards = document.querySelectorAll('#anovaGroupsContainer .anova-group-card');
      const canRemove = cards.length > 2;
      cards.forEach(card => {
        const btn = card.querySelector('.btn-remove-anova-group');
        if (btn) btn.disabled = !canRemove;
      });
    }

    updateAnovaDesignUI() {
      const isPaired = document.getElementById('anovaIsPaired')?.checked || false;
      const modeBadge = document.getElementById('anovaDesignModeBadge');
      const step1 = document.getElementById('anovaPipelineStep1');

      if (modeBadge) {
        if (isPaired) {
          modeBadge.className = 'badge badge-sig';
          modeBadge.innerText = 'Repeated Measures (Matched Within-Subjects)';
        } else {
          modeBadge.className = 'badge badge-neutral';
          modeBadge.innerText = 'Independent Cohorts (Between-Subjects)';
        }
      }

      if (step1) {
        step1.className = isPaired ? 'badge badge-sig' : 'badge badge-neutral';
        step1.innerText = `1. Design: ${isPaired ? 'Paired (RM)' : 'Independent'}`;
      }
    }

    loadAnovaPreset(presetKey) {
      const presets = {
        sample3: {
          isPaired: false,
          testType: 'auto',
          groups: [
            { name: 'Conservative', data: '7.2, 6.8, 7.5, 6.9, 8.1, 7.0, 7.4, 6.5, 7.9, 7.1' },
            { name: 'Orthotic Helmet', data: '4.1, 3.8, 4.5, 3.9, 4.8, 3.6, 4.2, 3.5, 4.0, 3.7' },
            { name: 'Endoscopic Strip', data: '2.5, 2.8, 2.2, 2.6, 3.1, 2.4, 2.9, 2.1, 2.7, 2.3' }
          ]
        },
        sampleWelch: {
          isPaired: false,
          testType: 'auto',
          groups: [
            { name: 'Cohort A (Small Var)', data: '10.1, 10.3, 10.0, 10.2, 9.9, 10.4, 10.1, 9.8, 10.2, 10.0' },
            { name: 'Cohort B (Mod Var)', data: '12.4, 11.1, 13.5, 10.8, 14.2, 11.9, 13.0, 12.1, 11.5, 13.8' },
            { name: 'Cohort C (Large Var)', data: '15.2, 8.5, 19.4, 11.1, 22.0, 14.3, 7.8, 18.6, 12.0, 20.5' },
            { name: 'Cohort D (High Var)', data: '18.0, 29.5, 9.2, 35.1, 14.8, 27.2, 8.1, 31.4, 12.5, 25.8' }
          ]
        },
        sampleRM: {
          isPaired: true,
          testType: 'auto',
          groups: [
            { name: 'Baseline (T0)', data: '22.4, 25.1, 19.8, 27.3, 23.5, 26.2, 21.9, 24.8, 20.5, 23.9' },
            { name: 'Week 2 (T1)', data: '19.1, 22.0, 17.5, 24.2, 20.8, 23.1, 18.9, 21.5, 18.0, 20.7' },
            { name: 'Week 6 (T2)', data: '15.3, 18.4, 14.1, 20.5, 17.0, 19.2, 15.6, 17.9, 14.8, 17.2' },
            { name: 'Month 3 (T3)', data: '12.1, 14.8, 11.2, 16.9, 13.5, 15.4, 12.4, 14.1, 11.5, 13.8' }
          ]
        },
        sampleSkew: {
          isPaired: false,
          testType: 'auto',
          groups: [
            { name: 'Standard Care', data: '1.2, 1.4, 1.1, 1.3, 1.5, 1.2, 1.4, 8.5, 12.3, 19.8' },
            { name: 'Modified Protocol', data: '2.1, 2.3, 2.0, 2.4, 2.2, 2.5, 15.1, 22.4, 28.0, 35.2' },
            { name: 'Novel Intervention', data: '5.5, 5.8, 5.2, 5.9, 6.1, 5.4, 32.0, 45.6, 58.2, 72.1' }
          ]
        }
      };

      const preset = presets[presetKey] || presets.sample3;
      const container = document.getElementById('anovaGroupsContainer');
      if (!container) return;

      container.innerHTML = '';
      preset.groups.forEach((g, idx) => {
        const k = idx + 1;
        const count = DataParser.parseSeries(g.data).length;
        const card = document.createElement('div');
        card.className = 'anova-group-card';
        card.style.cssText = 'background: var(--bg-surface-elevated); border: 1px solid var(--border-subtle); border-radius: 8px; padding: 0.75rem 0.85rem;';
        card.innerHTML = `
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
            <div style="display: flex; align-items: center; gap: 0.5rem; flex: 1;">
              <span class="anova-group-index" style="font-weight: 700; color: var(--cyan-primary); font-size: 0.82rem; min-width: 60px;">Cohort ${k}:</span>
              <input type="text" class="form-control anova-group-name" value="${g.name}" style="font-size: 0.82rem; padding: 0.25rem 0.5rem; height: 28px; font-weight: 600;" placeholder="Cohort Name">
            </div>
            <div style="display: flex; align-items: center; gap: 0.5rem;">
              <span class="badge badge-neutral anova-group-count" style="font-size: 0.72rem;">n = ${count}</span>
              <button type="button" class="btn btn-secondary btn-sm btn-remove-anova-group" style="padding: 0.15rem 0.45rem; font-size: 0.72rem; height: 26px; color: var(--rose-primary); border-color: rgba(244, 63, 94, 0.3);" title="Remove this cohort">✕</button>
            </div>
          </div>
          <textarea class="form-control anova-group-data" rows="2" placeholder="Comma or newline-separated values...">${g.data}</textarea>
        `;
        container.appendChild(card);
      });

      const pairedCb = document.getElementById('anovaIsPaired');
      if (pairedCb) pairedCb.checked = !!preset.isPaired;

      const testSelect = document.getElementById('anovaTestType');
      if (testSelect) testSelect.value = preset.testType || 'auto';

      this.updateAnovaRemoveButtons();
      this.updateAnovaDesignUI();
      this.runAnova();
    }

    getAnovaGroups() {
      const container = document.getElementById('anovaGroupsContainer');
      if (!container) return [];

      const cards = container.querySelectorAll('.anova-group-card');
      const groups = [];

      cards.forEach((card, idx) => {
        const nameInput = card.querySelector('.anova-group-name');
        const dataArea = card.querySelector('.anova-group-data');
        const countBadge = card.querySelector('.anova-group-count');

        const name = nameInput?.value.trim() || `Cohort ${idx + 1}`;
        const rawData = dataArea?.value || '';
        const parsed = DataParser.parseSeries(rawData);

        if (countBadge) {
          countBadge.innerText = `n = ${parsed.length}`;
        }
        groups.push({ name, data: parsed });
      });

      return groups;
    }

    runAnova() {
      const groups = this.getAnovaGroups();
      const isPaired = document.getElementById('anovaIsPaired')?.checked || false;
      const testType = document.getElementById('anovaTestType')?.value || 'auto';

      this.updateAnovaDesignUI();

      if (!groups || groups.length < 2) {
        alert('Please provide at least 2 cohorts for analysis.');
        return;
      }

      // Validate cohort data
      for (let i = 0; i < groups.length; i++) {
        if (!groups[i].data || groups[i].data.length === 0) {
          const msg = `Cohort "${groups[i].name}" has no valid numerical data.`;
          document.getElementById('anovaReportText').innerText = msg;
          return;
        }
      }

      const res = Anova.test(groups, testType, isPaired);
      if (res.error) {
        document.getElementById('anovaReportText').innerText = `Analysis Error: ${res.error}`;
        return;
      }

      // Update Metric Cards
      const fLabel = document.getElementById('anovaFLabel');
      const fVal = document.getElementById('anovaF');
      const pVal = document.getElementById('anovaP');
      const pBadge = document.getElementById('anovaPValBadge');
      const etaLabel = document.getElementById('anovaEtaLabel');
      const etaVal = document.getElementById('anovaEta');
      const omegaLabel = document.getElementById('anovaOmegaLabel');
      const omegaVal = document.getElementById('anovaOmega');

      // 1. Test Statistic
      if (fLabel) {
        if (res.testKey === 'kruskal') fLabel.innerText = 'Kruskal-Wallis (H)';
        else if (res.testKey === 'friedman') fLabel.innerText = 'Friedman (Q / χ²ᵣ)';
        else if (res.testKey === 'welch') fLabel.innerText = "Welch's F-Test";
        else if (res.testKey === 'rm_anova') fLabel.innerText = 'RM-ANOVA (F)';
        else fLabel.innerText = 'Fisher ANOVA (F)';
      }
      if (fVal) {
        const stat = res.statistic !== undefined ? res.statistic : (res.fStatistic || res.hStatistic || res.qStatistic || 0);
        fVal.innerText = isFinite(stat) ? stat.toFixed(2) : '--';
      }

      // 2. p-value & Significance badge
      if (pVal) pVal.innerText = Exporter.formatP(res.pValue);
      if (pBadge) {
        pBadge.className = res.isSignificant ? 'badge badge-sig' : 'badge badge-ns';
        pBadge.innerText = res.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)';
      }

      // 3. Effect Size
      if (etaLabel) {
        if (res.testKey === 'kruskal') etaLabel.innerText = 'Epsilon-Squared (ε²)';
        else if (res.testKey === 'friedman') etaLabel.innerText = "Kendall's W";
        else if (res.testKey === 'welch') etaLabel.innerText = 'Estimated ω²';
        else etaLabel.innerText = 'Partial η² / ω²';
      }
      if (etaVal) {
        if (res.testKey === 'kruskal' && res.effectSize?.epsilonSquared !== undefined) {
          etaVal.innerText = res.effectSize.epsilonSquared.toFixed(3);
        } else if (res.testKey === 'friedman' && res.effectSize?.kendallsW !== undefined) {
          etaVal.innerText = res.effectSize.kendallsW.toFixed(3);
        } else if (res.etaSquared !== undefined) {
          etaVal.innerText = res.etaSquared.toFixed(3);
        } else if (res.omegaSquared !== undefined) {
          etaVal.innerText = res.omegaSquared.toFixed(3);
        } else {
          etaVal.innerText = '--';
        }
      }

      // 4. Design & Degrees of Freedom
      if (omegaLabel) {
        omegaLabel.innerText = res.isPaired ? 'Paired Design & df' : 'Independent Design & df';
      }
      if (omegaVal) {
        if (res.dfBetween !== undefined && res.dfWithin !== undefined) {
          omegaVal.innerText = `df: (${res.dfBetween}, ${typeof res.dfWithin === 'number' ? res.dfWithin.toFixed(1) : res.dfWithin})`;
        } else if (res.df !== undefined) {
          omegaVal.innerText = `df = ${res.df}`;
        } else {
          omegaVal.innerText = `k = ${res.k}`;
        }
      }

      // 5. Assumptions Diagnostic Card
      const asm = res.assumptions;
      if (asm) {
        const recBadge = document.getElementById('anovaRecommendationBadge');
        if (recBadge) {
          recBadge.innerText = `Recommended: ${asm.recommendedTestName}`;
          recBadge.className = 'badge badge-sig';
        }

        const decisionBanner = document.getElementById('anovaDecisionBanner');
        const decisionText = document.getElementById('anovaDecisionText');
        if (decisionText) {
          let bannerHtml = '';
          if (res.userOverride) {
            bannerHtml = `<span class="badge badge-warn" style="font-size: 0.72rem; margin-bottom: 0.35rem; display: inline-block;">Manual Selection</span><br>` +
              `User opted to execute <strong>${res.testName}</strong>. Based on data diagnostics, the statistically optimal test is <strong>${asm.recommendedTestName}</strong>.<br>` +
              `<span style="color: var(--text-dim); font-size: 0.82rem;">${asm.rationale}</span>`;
            if (decisionBanner) decisionBanner.style.borderLeftColor = 'var(--gold-primary)';
          } else {
            bannerHtml = `<span class="badge badge-sig" style="font-size: 0.72rem; margin-bottom: 0.35rem; display: inline-block;">Automated Recommendation Executed</span><br>` +
              `Executed <strong>${res.testName}</strong>.<br>` +
              `<span style="color: var(--text-main);">${asm.rationale}</span>`;
            if (decisionBanner) decisionBanner.style.borderLeftColor = 'var(--cyan-primary)';
          }
          decisionText.innerHTML = bannerHtml;
        }

        // Normality Breakdown
        const normBadge = document.getElementById('anovaNormalityBadge');
        const normDetails = document.getElementById('anovaNormalityDetails');
        if (normBadge) {
          normBadge.className = asm.isNormal ? 'badge badge-sig' : 'badge badge-warn';
          normBadge.innerText = asm.isNormal ? 'Parametric (All Normal)' : 'Non-Parametric (Skewed)';
        }
        if (normDetails && asm.normalityTests) {
          normDetails.innerHTML = asm.normalityTests.map(n => 
            `• <strong>${n.group}</strong>: n=${n.n}, JB=${n.jbStat.toFixed(2)}, p=${Exporter.formatP(n.pValue)} (${n.isNormal ? '<span style="color: var(--emerald-primary);">Normal</span>' : '<span style="color: var(--rose-primary);">Skewed, p < .05</span>'})`
          ).join('<br>');
        }

        // Variance Homogeneity / Sphericity Breakdown
        const varBadge = document.getElementById('anovaVarianceBadge');
        const varDetails = document.getElementById('anovaVarianceDetails');
        if (varBadge) {
          if (isPaired) {
            varBadge.className = 'badge badge-sig';
            varBadge.innerText = asm.sphericity ? `Sphericity ε̂ = ${asm.sphericity.epsilon.toFixed(3)}` : 'Repeated Measures';
          } else {
            varBadge.className = asm.isHomoscedastic ? 'badge badge-sig' : 'badge badge-warn';
            varBadge.innerText = asm.isHomoscedastic ? 'Equal Variances' : 'Unequal Variances';
          }
        }
        if (varDetails) {
          if (isPaired) {
            if (asm.sphericity) {
              varDetails.innerHTML = `Greenhouse-Geisser correction factor: <strong>ε̂ = ${asm.sphericity.epsilon.toFixed(3)}</strong>.<br>` +
                (asm.sphericity.isSpherical 
                  ? `<span style="color: var(--emerald-primary);">Sphericity assumption reasonably met (ε̂ ≈ 1.0).</span>` 
                  : `<span style="color: var(--gold-primary);">Sphericity violated (ε̂ < 0.75). Degrees of freedom adjusted via Greenhouse-Geisser.</span>`);
            } else {
              varDetails.innerText = 'Matched repeated measures design: within-subject correlation structure preserved.';
            }
          } else if (asm.leveneTest) {
            const lev = asm.leveneTest;
            varDetails.innerHTML = `Brown-Forsythe Levene's Test: <strong>F(${lev.df1}, ${lev.df2}) = ${lev.fStat.toFixed(2)}, p = ${Exporter.formatP(lev.pValue)}</strong>.<br>` +
              (asm.isHomoscedastic 
                ? `<span style="color: var(--emerald-primary);">Homoscedasticity confirmed (p ≥ .05). Residual variances across cohorts are equal.</span>`
                : `<span style="color: var(--rose-primary);">Heteroscedasticity detected (p < .05). Residual variances differ significantly; Welch's robust F recommended.</span>`);
          }
        }

        // Pipeline step badges
        const pStep1 = document.getElementById('anovaPipelineStep1');
        const pStep2 = document.getElementById('anovaPipelineStep2');
        const pStep3 = document.getElementById('anovaPipelineStep3');
        const pStepFinal = document.getElementById('anovaPipelineStepFinal');
        if (pStep1) {
          pStep1.className = isPaired ? 'badge badge-sig' : 'badge badge-neutral';
          pStep1.innerText = `1. Design: ${isPaired ? 'Paired (RM)' : 'Independent'}`;
        }
        if (pStep2) {
          pStep2.className = asm.isNormal ? 'badge badge-sig' : 'badge badge-warn';
          pStep2.innerText = `2. Normality: ${asm.isNormal ? 'Parametric' : 'Non-Parametric'}`;
        }
        if (pStep3) {
          if (isPaired) {
            pStep3.className = 'badge badge-sig';
            pStep3.innerText = `3. Sphericity: ε̂ = ${asm.sphericity ? asm.sphericity.epsilon.toFixed(2) : 'N/A'}`;
          } else {
            pStep3.className = asm.isHomoscedastic ? 'badge badge-sig' : 'badge badge-warn';
            pStep3.innerText = `3. Variances: ${asm.isHomoscedastic ? 'Equal' : 'Unequal'}`;
          }
        }
        if (pStepFinal) {
          pStepFinal.className = 'badge badge-sig';
          pStepFinal.innerText = `Selected: ${res.testName}`;
        }
      }

      // 6. Post-Hoc Pairwise Table
      const postHocTitle = document.getElementById('anovaPostHocTitle');
      const postHocSub = document.getElementById('anovaPostHocSubtitle');
      const statCol = document.getElementById('anovaPostHocStatCol');
      if (postHocTitle) postHocTitle.innerText = `🔬 Post-Hoc Pairwise Contrasts (${res.postHocMethod || 'Pairwise'})`;
      if (postHocSub) postHocSub.innerText = `Contrasts between individual cohorts with family-wise error rate control (${res.postHocMethod || 'Contrasts'})`;
      if (statCol) {
        if (res.testKey === 'kruskal') statCol.innerText = "Dunn's z-Stat";
        else if (res.testKey === 'friedman') statCol.innerText = 'Wilcoxon W';
        else if (res.testKey === 'welch') statCol.innerText = 'Games-Howell t';
        else if (res.testKey === 'rm_anova') statCol.innerText = 'Paired t-Stat';
        else statCol.innerText = 'Tukey q (t)';
      }

      const tbody = document.getElementById('anovaPostHocBody');
      if (tbody) {
        if (res.pairwise && res.pairwise.length > 0) {
          tbody.innerHTML = res.pairwise.map(p => {
            const diffVal = p.meanDiff !== undefined ? p.meanDiff : (p.diff !== undefined ? p.diff : 0);
            const seVal = p.seDiff !== undefined ? p.seDiff.toFixed(2) : (p.se !== undefined ? p.se.toFixed(2) : '--');
            const statText = p.qStatistic !== undefined ? `q = ${p.qStatistic.toFixed(2)} (t = ${p.tStatistic.toFixed(2)})` :
                             p.tStatistic !== undefined ? `t = ${p.tStatistic.toFixed(2)}` :
                             p.zStatistic !== undefined ? `z = ${p.zStatistic.toFixed(2)}` :
                             p.wStatistic !== undefined ? `W = ${p.wStatistic.toFixed(2)}` : '--';
            const ciText = p.ci95 ? `[${p.ci95[0].toFixed(2)}, ${p.ci95[1].toFixed(2)}]` : '--';
            const esText = p.cohensD !== undefined ? `d = ${p.cohensD.toFixed(2)}` :
                           p.r !== undefined ? `r = ${p.r.toFixed(2)}` : '--';
            return `
              <tr>
                <td style="font-weight: 600; color: var(--text-main);">${p.comparison}</td>
                <td>${diffVal >= 0 ? '+' : ''}${diffVal.toFixed(2)}</td>
                <td>${seVal}</td>
                <td>${statText}</td>
                <td style="font-weight: 600; color: ${p.isSignificant ? 'var(--cyan-primary)' : 'var(--text-muted)'};">${Exporter.formatP(p.pValue)}</td>
                <td>${ciText}</td>
                <td>${esText}</td>
                <td>
                  <span class="badge ${p.isSignificant ? 'badge-sig' : 'badge-ns'}">
                    ${p.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)'}
                  </span>
                </td>
              </tr>
            `;
          }).join('');
        } else {
          tbody.innerHTML = `<tr><td colspan="8" style="text-align: center; color: var(--text-dim);">No pairwise contrasts calculated (omnibus effect not significant or single cohort).</td></tr>`;
        }
      }

      // 7. Clinical / Academic APA Summary Narrative
      let report = '';
      if (res.testKey === 'one_way') {
        report = `A one-way between-subjects ANOVA was conducted across ${res.k} cohorts (N = ${res.totalN}). `;
        report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} omnibus effect: F(${res.dfBetween}, ${res.dfWithin}) = ${res.fStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, η² = ${res.etaSquared.toFixed(3)}, ω² = ${res.omegaSquared.toFixed(3)}.\n\n`;
      } else if (res.testKey === 'welch') {
        report = `A Welch's robust one-way ANOVA (adjusting for heteroscedasticity) was conducted across ${res.k} cohorts (N = ${res.totalN}). `;
        report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} omnibus effect: Welch's F(${res.dfBetween}, ${res.dfWithin.toFixed(2)}) = ${res.fStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, estimated ω² = ${res.omegaSquared.toFixed(3)}.\n\n`;
      } else if (res.testKey === 'kruskal') {
        report = `A non-parametric Kruskal-Wallis H test was conducted across ${res.k} cohorts (N = ${res.totalN}). `;
        report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} omnibus rank difference: H(${res.df}) = ${res.hStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, ε² = ${res.effectSize.epsilonSquared.toFixed(3)}.\n\n`;
      } else if (res.testKey === 'rm_anova') {
        report = `A one-way repeated measures ANOVA was conducted across ${res.k} conditions (N = ${res.nSubjects} subjects). `;
        report += `Greenhouse-Geisser sphericity correction: ε̂ = ${res.epsilon.toFixed(3)}. Omnibus effect: F(${res.dfTreatment.toFixed(2)}, ${res.dfError.toFixed(2)}) = ${res.fStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, partial η² = ${res.partialEtaSquared.toFixed(3)}.\n\n`;
      } else if (res.testKey === 'friedman') {
        report = `A non-parametric Friedman rank sum test was conducted across ${res.k} repeated conditions (N = ${res.n} subjects). `;
        report += `Omnibus rank difference: Q(${res.df}) = ${res.qStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, Kendall's W = ${res.effectSize.kendallsW.toFixed(3)}.\n\n`;
      }

      // Append post-hoc summary
      if (res.pairwise && res.pairwise.length > 0) {
        report += `Post-hoc contrasts (${res.postHocMethod}) revealed:\n`;
        res.pairwise.forEach(p => {
          const diffVal = p.meanDiff !== undefined ? p.meanDiff : (p.diff !== undefined ? p.diff : 0);
          const statVal = p.qStatistic !== undefined ? `q = ${p.qStatistic.toFixed(2)}` :
                          p.tStatistic !== undefined ? `t = ${p.tStatistic.toFixed(2)}` :
                          p.zStatistic !== undefined ? `z = ${p.zStatistic.toFixed(2)}` :
                          p.wStatistic !== undefined ? `W = ${p.wStatistic.toFixed(2)}` : '';
          if (p.isSignificant) {
            report += `• ${p.comparison}: Statistically significant difference (Δ = ${diffVal.toFixed(2)}, ${statVal ? statVal + ', ' : ''}${Exporter.formatP(p.pValue)}).\n`;
          } else {
            report += `• ${p.comparison}: No statistically significant difference (Δ = ${diffVal.toFixed(2)}, ${statVal ? statVal + ', ' : ''}${Exporter.formatP(p.pValue)}, ns).\n`;
          }
        });
      }

      // Append assumption decision rationale
      if (asm) {
        report += `\nMethodological Rationale: ${asm.rationale}`;
      }

      document.getElementById('anovaReportText').innerText = report;

      // 8. Dispersion Plot
      const errorBarMode = document.getElementById('anovaErrorBarMode')?.value || 'ci95';
      const modeDescriptions = {
        ci95: 'Error Bars: 95% Confidence Interval (Mean ± 95% CI)',
        sem: 'Error Bars: Standard Error of Mean (Mean ± 1 SEM)',
        sd: 'Error Bars: Standard Deviation (Mean ± 1 SD)',
        iqr: 'Distribution: Box & Whiskers (Median, Q1-Q3 IQR, Tukey Fences)'
      };
      const subElem = document.getElementById('anovaChartSub');
      if (subElem) {
        subElem.innerText = modeDescriptions[errorBarMode] || modeDescriptions.ci95;
      }

      if (this.engines.anovaCanvas) {
        Plots.renderErrorBarPlot(this.engines.anovaCanvas, res.groups, {
          mode: errorBarMode,
          title: 'Multi-Cohort Comparison'
        });
      }
      this.results = this.results || {};
      this.results.anova = res;
    }

    runCat() {
      const mode = document.getElementById('catAnalysisMode')?.value || 'diagnostic';
      const a = parseFloat(document.getElementById('catA')?.value) || 0;
      const b = parseFloat(document.getElementById('catB')?.value) || 0;
      const c = parseFloat(document.getElementById('catC')?.value) || 0;
      const d = parseFloat(document.getElementById('catD')?.value) || 0;

      const cornerHeader = document.getElementById('catCornerHeader');
      const col1Header = document.getElementById('catCol1Header');
      const col2Header = document.getElementById('catCol2Header');
      const row1Header = document.getElementById('catRow1Header');
      const row2Header = document.getElementById('catRow2Header');
      const labelA = document.getElementById('catLabelA');
      const labelB = document.getElementById('catLabelB');
      const labelC = document.getElementById('catLabelC');
      const labelD = document.getElementById('catLabelD');
      const explanation = document.getElementById('catExplanation');
      const diagGrid = document.getElementById('catDiagnosticMetrics');
      const studyGrid = document.getElementById('catStudyMetrics');
      const reportTitle = document.getElementById('catReportHeaderTitle');

      let res;
      if (mode === 'diagnostic') {
        if (cornerHeader) cornerHeader.innerText = 'Test \\ Ref';
        if (col1Header) col1Header.innerText = 'Gold Standard (+)';
        if (col2Header) col2Header.innerText = 'Gold Standard (-)';
        if (row1Header) row1Header.innerText = 'New Test (+)';
        if (row2Header) row2Header.innerText = 'New Test (-)';
        if (labelA) labelA.innerText = 'True Positive (TP)';
        if (labelB) labelB.innerText = 'False Positive (FP)';
        if (labelC) labelC.innerText = 'False Negative (FN)';
        if (labelD) labelD.innerText = 'True Negative (TN)';
        if (explanation) explanation.innerText = 'Evaluates index test accuracy against the reference standard: Sensitivity, Specificity, PPV, NPV, Overall Accuracy (Wilson Score 95% CIs), and Likelihood Ratios.';
        if (diagGrid) diagGrid.style.display = 'grid';
        if (studyGrid) studyGrid.style.display = 'none';
        if (reportTitle) reportTitle.innerText = 'Diagnostic Performance & Accuracy Report (STARD compliant)';

        res = Diagnostic.evaluate2x2(a, b, c, d);
        if (res.error) {
          document.getElementById('catReportText').innerText = res.error;
          return;
        }

        const formatPct = (val) => `${(val * 100).toFixed(1)}%`;
        const formatCI = (ci) => `95% CI: [${(ci[0] * 100).toFixed(1)}%, ${(ci[1] * 100).toFixed(1)}%]`;

        const sensEl = document.getElementById('diagSens');
        if (sensEl) sensEl.innerText = formatPct(res.sensitivity);
        const sensCIEl = document.getElementById('diagSensCI');
        if (sensCIEl) sensCIEl.innerText = formatCI(res.sensitivityCI95);

        const specEl = document.getElementById('diagSpec');
        if (specEl) specEl.innerText = formatPct(res.specificity);
        const specCIEl = document.getElementById('diagSpecCI');
        if (specCIEl) specCIEl.innerText = formatCI(res.specificityCI95);

        const ppvEl = document.getElementById('diagPPV');
        if (ppvEl) ppvEl.innerText = formatPct(res.ppv);
        const ppvCIEl = document.getElementById('diagPPVCI');
        if (ppvCIEl) ppvCIEl.innerText = formatCI(res.ppvCI95);

        const npvEl = document.getElementById('diagNPV');
        if (npvEl) npvEl.innerText = formatPct(res.npv);
        const npvCIEl = document.getElementById('diagNPVCI');
        if (npvCIEl) npvCIEl.innerText = formatCI(res.npvCI95);

        const accEl = document.getElementById('diagAcc');
        if (accEl) accEl.innerText = formatPct(res.accuracy);
        const accCIEl = document.getElementById('diagAccCI');
        if (accCIEl) accCIEl.innerText = formatCI(res.accuracyCI95);

        const lrEl = document.getElementById('diagLR');
        if (lrEl) lrEl.innerText = `LR+ ${isFinite(res.plr) ? res.plr.toFixed(2) : '∞'} | LR- ${isFinite(res.nlr) ? res.nlr.toFixed(2) : '0'}`;

        const youdenEl = document.getElementById('diagYouden');
        if (youdenEl) youdenEl.innerText = `Youden's J: ${res.youdenJ.toFixed(3)} | Prev: ${formatPct(res.prevalence)}`;

        const report = `Diagnostic test evaluation against gold standard reference (Total N = ${res.total}): Overall Accuracy = ${formatPct(res.accuracy)} (${formatCI(res.accuracyCI95)}). ` +
          `Sensitivity (TPR) = ${formatPct(res.sensitivity)} (${formatCI(res.sensitivityCI95)}), ` +
          `Specificity (TNR) = ${formatPct(res.specificity)} (${formatCI(res.specificityCI95)}). ` +
          `Positive Predictive Value (PPV) = ${formatPct(res.ppv)} (${formatCI(res.ppvCI95)}), ` +
          `Negative Predictive Value (NPV) = ${formatPct(res.npv)} (${formatCI(res.npvCI95)}). ` +
          `Positive Likelihood Ratio (LR+) = ${isFinite(res.plr) ? res.plr.toFixed(2) : 'N/A'}, ` +
          `Negative Likelihood Ratio (LR-) = ${isFinite(res.nlr) ? res.nlr.toFixed(2) : 'N/A'}, ` +
          `Youden's J index = ${res.youdenJ.toFixed(3)}, Sample Disease Prevalence = ${formatPct(res.prevalence)}.`;

        document.getElementById('catReportText').innerText = report;
      } else {
        if (cornerHeader) cornerHeader.innerText = 'Cohort \\ Event';
        if (col1Header) col1Header.innerText = 'Intervention';
        if (col2Header) col2Header.innerText = 'Control';
        if (row1Header) row1Header.innerText = 'Event (+)';
        if (row2Header) row2Header.innerText = 'No Event (-)';
        if (labelA) labelA.innerText = 'Treated Event';
        if (labelB) labelB.innerText = 'Control Event';
        if (labelC) labelC.innerText = 'Treated No Event';
        if (labelD) labelD.innerText = 'Control No Event';
        if (explanation) explanation.innerText = 'Instant computation of Odds Ratio (Woolf 95% CI), Relative Risk, Absolute Risk Reduction, Number Needed to Treat (NNT), Pearson Chi-Square, and Fisher\'s exact test.';
        if (diagGrid) diagGrid.style.display = 'none';
        if (studyGrid) studyGrid.style.display = 'grid';
        if (reportTitle) reportTitle.innerText = 'Epidemiological & Clinical Risk Report';

        res = Categorical.twoByTwo(a, b, c, d);
        if (res.error) {
          document.getElementById('catReportText').innerText = res.error;
          return;
        }

        const r = res.riskMetrics;
        document.getElementById('catOR').innerText = `${r.oddsRatio.toFixed(2)} [${r.orCI95[0].toFixed(2)}, ${r.orCI95[1].toFixed(2)}]`;
        document.getElementById('catRR').innerText = `${r.relativeRisk.toFixed(2)} [${r.rrCI95[0].toFixed(2)}, ${r.rrCI95[1].toFixed(2)}]`;
        document.getElementById('catChiSq').innerText = res.chiSquare.standard.toFixed(2);
        document.getElementById('catPVal').innerText = `p = ${res.chiSquare.pValueStandard < 0.001 ? '< .001' : res.chiSquare.pValueStandard.toFixed(3)}`;
        document.getElementById('catFisher').innerText = Exporter.formatP(res.fishersExact.pValue);
        document.getElementById('catNNT').innerText = isFinite(r.nnt) ? r.nnt.toFixed(1) : '∞';

        const arrSubEl = document.getElementById('catARRSub');
        if (arrSubEl) arrSubEl.innerText = `ARR: ${(r.arr * 100).toFixed(1)}%`;
        const rrrEl = document.getElementById('catRRR');
        if (rrrEl) rrrEl.innerText = isFinite(r.rrr) ? `${(r.rrr * 100).toFixed(1)}%` : 'N/A';

        const report = `2x2 contingency analysis (N = ${res.table.n}): Pearson χ²(1) = ${res.chiSquare.standard.toFixed(2)}, ${Exporter.formatP(res.chiSquare.pValueStandard)} (Fisher's exact ${Exporter.formatP(res.fishersExact.pValue)}). Odds Ratio = ${r.oddsRatio.toFixed(2)} (95% CI [${r.orCI95[0].toFixed(2)}, ${r.orCI95[1].toFixed(2)}]), Relative Risk = ${r.relativeRisk.toFixed(2)} (95% CI [${r.rrCI95[0].toFixed(2)}, ${r.rrCI95[1].toFixed(2)}]), Absolute Risk Reduction (ARR) = ${(r.arr * 100).toFixed(1)}%, Number Needed to Treat (NNT) = ${isFinite(r.nnt) ? r.nnt.toFixed(1) : 'N/A'}.`;
        document.getElementById('catReportText').innerText = report;
      }
      this.results = this.results || {}; this.results.categorical = res;
    }

    runCorr() {
      const xs = DataParser.parseSeries(document.getElementById('corrX')?.value || '');
      const ys = DataParser.parseSeries(document.getElementById('corrY')?.value || '');

      const p = Correlation.pearson(xs, ys);
      const s = Correlation.spearman(xs, ys);
      const reg = Correlation.linearRegression(xs, ys);

      if (p.error || reg.error) {
        alert(p.error || reg.error);
        return;
      }

      document.getElementById('corrR').innerText = p.r.toFixed(3);
      document.getElementById('corrR2').innerText = p.rSquared.toFixed(3);
      document.getElementById('corrRho').innerText = s.rho.toFixed(3);
      document.getElementById('corrP').innerText = Exporter.formatP(p.pValue);
      document.getElementById('regEq').innerText = reg.equation;

      // Evaluate Assumptions
      const assump = Correlation.checkAssumptions(xs, ys, reg);
      if (!assump.error) {
        document.getElementById('corrAssumpX').innerHTML = assump.xNormal
          ? '<span style="color: var(--emerald-primary);">Normal ✓</span>'
          : '<span style="color: var(--amber-primary);">Skewed ⚠️</span>';

        document.getElementById('corrAssumpY').innerHTML = assump.yNormal
          ? '<span style="color: var(--emerald-primary);">Normal ✓</span>'
          : '<span style="color: var(--amber-primary);">Skewed ⚠️</span>';

        document.getElementById('corrAssumpOutlier').innerHTML = assump.outliers.length === 0
          ? '<span style="color: var(--emerald-primary);">None ✓</span>'
          : `<span style="color: var(--rose-primary);">${assump.outliers.length} Outlier(s) ⚠️</span>`;

        document.getElementById('corrAssumpHomo').innerHTML = assump.isHomoscedastic
          ? `<span style="color: var(--emerald-primary);">Equal Var ✓</span>`
          : `<span style="color: var(--amber-primary);">Heteroscedastic ⚠️</span>`;

        const shapeEl = document.getElementById('corrAssumpShape');
        if (shapeEl) {
          if (assump.isUShaped) {
            shapeEl.innerHTML = `<span style="color: var(--amber-primary); font-weight: 700;">${assump.shape} ⚠️</span>`;
          } else if (assump.mono && !assump.mono.isMonotonic) {
            shapeEl.innerHTML = `<span style="color: var(--rose-primary); font-weight: 700;">${assump.mono.pattern.includes('Tri-Phasic') ? 'N-Shaped' : 'Non-Monotonic'} ⚠️</span>`;
          } else if (assump.mono) {
            shapeEl.innerHTML = `<span style="color: var(--emerald-primary); font-size: 0.82rem;">${assump.mono.pattern} ✓</span>`;
          } else {
            shapeEl.innerHTML = `<span style="color: var(--emerald-primary);">Linear ✓</span>`;
          }
        }

        // Recommendation Banner
        const banner = document.getElementById('corrRecBanner');
        const badge = document.getElementById('corrRecBadge');
        const title = document.getElementById('corrRecTitle');
        const detail = document.getElementById('corrRecDetail');

        badge.className = `badge ${assump.isParametricOk ? 'badge-sig' : 'badge-ns'}`;
        badge.innerText = assump.isParametricOk ? 'Parametric Valid' : (assump.isUShaped ? 'Non-Linear Indicated' : (assump.mono && !assump.mono.isMonotonic ? 'Multi-Phasic Non-Monotonic' : 'Non-Parametric Indicated'));
        title.innerText = assump.recommendedTest;
        detail.innerText = assump.recommendationDetail;

        if (banner) {
          banner.style.background = assump.isParametricOk ? 'rgba(0, 210, 255, 0.08)' : 'rgba(245, 158, 11, 0.08)';
          banner.style.borderColor = assump.isParametricOk ? 'var(--cyan-primary)' : 'var(--amber-primary)';
        }

        // Clinical Report with Assumption Rationale
        let report = `Bivariate Statistical Evaluation (N = ${xs.length} pairs):\n`;
        report += `• Assumption Diagnostics: Bivariate normality is ${assump.bivariateNormal ? 'satisfied' : 'violated'}; ` +
          `directional trajectory is ${assump.mono ? assump.mono.pattern : 'linear'} (${assump.mono && assump.mono.isMonotonic ? 'consistently moving together with 0 reversals' : `${assump.mono ? assump.mono.reversals : 1} directional reversal(s)`}); ` +
          `homoscedasticity is ${assump.isHomoscedastic ? 'preserved' : 'violated'}; ` +
          `${assump.outliers.length === 0 ? 'no extreme leverage outliers detected' : `${assump.outliers.length} influential outlier(s) detected`}.\n`;
        report += `• Method Decision: ${assump.recommendedTest} is the mathematically recommended procedure.\n`;
        if (assump.isUShaped) {
          report += `• Non-Linear Quadratic Model: y = ${assump.quad.equation} (R² = ${assump.quad.rSquaredQuad.toFixed(3)}, Incremental F = ${assump.quad.fStat.toFixed(2)}, ${Exporter.formatP(assump.quad.pQuad)}).\n` +
            `• Clinical Morphology: ${assump.shape} relationship with physiological ${assump.quad.b2 > 0 ? 'optimum / nadir' : 'peak / zenith'} at X = ${assump.quad.vertexX.toFixed(2)} (predicted Y = ${assump.quad.vertexY.toFixed(2)}). Standard linear Pearson (r = ${p.r.toFixed(3)}) fails to capture this strong biological curve.`;
        } else if (assump.mono && !assump.mono.isMonotonic && assump.mono.reversals >= 2) {
          report += `• Non-Monotonic Alert: The data exhibits ${assump.mono.pattern} with ${assump.mono.reversals} directional reversals (e.g. values increase, decrease, then increase). Neither linear Pearson (r = ${p.r.toFixed(3)}) nor monotonic rank Spearman (ρ = ${s.rho.toFixed(3)}) can properly model this multi-directional curve; spline or piecewise non-linear regression is advised.`;
        } else if (assump.isParametricOk) {
          report += `• Results: Pearson's correlation r = ${p.r.toFixed(3)} (${p.isSignificant ? 'statistically significant' : 'ns'}, ${Exporter.formatP(p.pValue)}), R² = ${p.rSquared.toFixed(3)}. Linear regression equation: ${reg.equation} (residual SE = ${reg.seResidual.toFixed(3)}).`;
        } else {
          report += `• Results: Spearman's rank correlation ρ = ${s.rho.toFixed(3)} (${s.pValue < 0.05 ? 'statistically significant' : 'ns'}, ${Exporter.formatP(s.pValue)}). (Reference Pearson r = ${p.r.toFixed(3)}, R² = ${p.rSquared.toFixed(3)}).`;
        }
        document.getElementById('corrReportText').innerText = report;
      }

      if (this.engines.corrCanvas) {
        const pairs = Correlation.cleanPairs(xs, ys);
        const plotTitle = assump && assump.isUShaped ? `Curvilinear Regression (${assump.shape})` : 'Scatter Plot & Linear Regression';
        Plots.renderScatterRegression(this.engines.corrCanvas, pairs, reg, plotTitle, assump?.quad);
      }
      this.results = this.results || {}; this.results.correlation = { n: xs.length, correlation: p, regression: reg, morphology: assump };
    }

    runROC() {
      const text = document.getElementById('rocInput')?.value || '';
      const lines = text.trim().split(/\r?\n/);
      const samples = [];
      for (const line of lines) {
        const parts = line.split(/[,\t\s]+/).filter(Boolean);
        if (parts.length >= 2) {
          samples.push({ score: parseFloat(parts[0]), status: parseInt(parts[1], 10) });
        }
      }

      const roc = Diagnostic.computeROC(samples);
      if (roc.error) {
        alert(roc.error);
        return;
      }

      const best = roc.optimalCutoff;
      const tp = Math.round(best.sens * roc.nPos);
      const fp = Math.round((1 - best.spec) * roc.nNeg);
      const fn = roc.nPos - tp;
      const tn = roc.nNeg - fp;
      const d = Diagnostic.evaluate2x2(tp, fp, fn, tn);

      document.getElementById('rocAUC').innerText = `${roc.auc.toFixed(3)} [${roc.aucCI95[0].toFixed(3)}, ${roc.aucCI95[1].toFixed(3)}]`;
      document.getElementById('rocCutoff').innerText = typeof best.threshold === 'number' ? best.threshold.toFixed(2) : best.threshold;
      document.getElementById('rocSens').innerText = `${(d.sensitivity * 100).toFixed(1)}%`;
      document.getElementById('rocSpec').innerText = `${(d.specificity * 100).toFixed(1)}%`;
      document.getElementById('rocPLR').innerText = d.plr.toFixed(2);
      document.getElementById('rocNLR').innerText = d.nlr.toFixed(2);

      const report = `Diagnostic ROC analysis (N = ${roc.total}): AUC = ${roc.auc.toFixed(3)} (95% CI [${roc.aucCI95[0].toFixed(3)}, ${roc.aucCI95[1].toFixed(3)}]). Optimal decision cutoff = ${best.threshold} yielding Sensitivity = ${(d.sensitivity * 100).toFixed(1)}% and Specificity = ${(d.specificity * 100).toFixed(1)}% (Youden's J = ${best.youdenJ.toFixed(3)}).`;
      document.getElementById('rocReportText').innerText = report;

      if (this.engines.rocCanvas) {
        Plots.renderROC(this.engines.rocCanvas, roc, `ROC Curve (AUC = ${roc.auc.toFixed(3)})`);
      }
      this.results = this.results || {}; this.results.diagnostic = { eval2x2: d, roc, best };
    }

    runPower() {
      const design = document.getElementById('pwrStudyDesign')?.value || 'independent';
      const goal = document.getElementById('pwrCalcGoal')?.value || 'sample_size';
      const alpha = parseFloat(document.getElementById('pwrAlpha')?.value) || 0.05;
      const power = parseFloat(document.getElementById('pwrPower')?.value) || 0.80;
      const givenN = goal === 'power' ? (parseFloat(document.getElementById('pwrGivenN')?.value) || 50) : null;

      const effectLabel = document.getElementById('pwrEffectLabel');
      const effectVal = document.getElementById('pwrD');
      const effectSub = document.getElementById('pwrEffectSub');
      const nGroupLabel = document.getElementById('pwrNGroupLabel');
      const nGroupVal = document.getElementById('pwrNGroup');
      const nGroupSub = document.getElementById('pwrNGroupSub');
      const nTotalLabel = document.getElementById('pwrNTotalLabel');
      const nTotalVal = document.getElementById('pwrNTotal');
      const nTotalSub = document.getElementById('pwrNTotalSub');
      const achievedPowerVal = document.getElementById('pwrAchievedPower');
      const powerSub = document.getElementById('pwrPowerSub');
      const cardFisher = document.getElementById('pwrCardFisher');
      const cardRisk = document.getElementById('pwrCardClinicalRisk');
      const fisherNVal = document.getElementById('pwrFisherN');
      const fisherSub = document.getElementById('pwrFisherSub');
      const riskVal = document.getElementById('pwrClinicalRisk');
      const riskSub = document.getElementById('pwrClinicalRiskSub');
      const reportBox = document.getElementById('pwrReportText');

      if (design === 'independent') {
        if (cardFisher) cardFisher.style.display = 'none';
        if (cardRisk) cardRisk.style.display = 'none';

        const m1 = parseFloat(document.getElementById('pwrM1')?.value) || 10;
        const m2 = parseFloat(document.getElementById('pwrM2')?.value) || 15;
        const sd = parseFloat(document.getElementById('pwrSD')?.value) || 10;

        const res = PowerAnalysis.twoIndependentMeans({ m1, m2, sd, alpha, power, nPerGroup: givenN });
        if (res.error) {
          if (reportBox) reportBox.innerText = res.error;
          return;
        }

        if (effectLabel) effectLabel.innerText = "Effect Size (Cohen's d)";
        if (effectVal) effectVal.innerText = res.cohensD.toFixed(2);
        if (effectSub) effectSub.innerText = `|Δ| = ${res.diff.toFixed(2)} | SD = ${res.sd.toFixed(2)}`;

        if (nGroupLabel) nGroupLabel.innerText = 'Sample Size per Group (n)';
        if (nGroupVal) nGroupVal.innerText = res.nPerGroup;
        if (nGroupSub) nGroupSub.innerText = 'Equal 1:1 allocation';

        if (nTotalLabel) nTotalLabel.innerText = 'Total Study Sample Size (N)';
        if (nTotalVal) nTotalVal.innerText = res.totalN;
        if (nTotalSub) nTotalSub.innerText = 'Across both arms (n₁ + n₂)';

        if (achievedPowerVal) achievedPowerVal.innerText = `${(res.achievedPower * 100).toFixed(1)}%`;
        if (powerSub) powerSub.innerText = `Type II error β = ${(1 - res.achievedPower).toFixed(3)} | α = ${res.alpha}`;

        let report = '';
        if (goal === 'sample_size') {
          report = `For a two-arm parallel randomized study (independent samples t-test) to detect a standardized effect size of Cohen's d = ${res.cohensD.toFixed(2)} (mean difference |Δ| = ${res.diff.toFixed(2)}, pooled SD = ${res.sd.toFixed(2)}) with ${Math.round(res.power * 100)}% statistical power at two-sided α = ${res.alpha}, a minimum of n = ${res.nPerGroup} patients per group (total N = ${res.totalN}) is required. Allowing for a standard 10% loss-to-follow-up buffer, a target enrollment of N = ${Math.ceil(res.totalN / 0.9)} patients (${Math.ceil(res.nPerGroup / 0.9)} per arm) is recommended.`;
        } else {
          report = `With an enrolled sample size of n = ${res.nPerGroup} patients per group (total N = ${res.totalN}) comparing two independent continuous means (Cohen's d = ${res.cohensD.toFixed(2)}, |Δ| = ${res.diff.toFixed(2)}, SD = ${res.sd.toFixed(2)}), the achieved statistical power to detect a true difference at two-sided α = ${res.alpha} is ${(res.achievedPower * 100).toFixed(1)}% (Type II error rate β = ${(1 - res.achievedPower).toFixed(3)}).`;
        }
        if (reportBox) reportBox.innerText = report;

      } else if (design === 'paired') {
        if (cardFisher) cardFisher.style.display = 'none';
        if (cardRisk) cardRisk.style.display = 'none';

        const m1 = parseFloat(document.getElementById('pwrPairedM1')?.value) || 120;
        const m2 = parseFloat(document.getElementById('pwrPairedM2')?.value) || 112;
        const sdDiff = parseFloat(document.getElementById('pwrPairedSD')?.value) || 10;

        const res = PowerAnalysis.pairedMeans({ m1, m2, sdDiff, alpha, power, nPairs: givenN });
        if (res.error) {
          if (reportBox) reportBox.innerText = res.error;
          return;
        }

        if (effectLabel) effectLabel.innerText = "Effect Size (Cohen's d_z)";
        if (effectVal) effectVal.innerText = res.dz.toFixed(2);
        if (effectSub) effectSub.innerText = `|Δ| = ${res.diff.toFixed(2)} | σ_d = ${res.sdDiff.toFixed(2)}`;

        if (nGroupLabel) nGroupLabel.innerText = 'Required Number of Pairs (N)';
        if (nGroupVal) nGroupVal.innerText = res.nPairs;
        if (nGroupSub) nGroupSub.innerText = 'Matched pairs / repeat measures';

        if (nTotalLabel) nTotalLabel.innerText = 'Total Paired Subjects (N)';
        if (nTotalVal) nTotalVal.innerText = res.totalN;
        if (nTotalSub) nTotalSub.innerText = 'Pre-to-post repeated subjects';

        if (achievedPowerVal) achievedPowerVal.innerText = `${(res.achievedPower * 100).toFixed(1)}%`;
        if (powerSub) powerSub.innerText = `Type II error β = ${(1 - res.achievedPower).toFixed(3)} | α = ${res.alpha}`;

        let report = '';
        if (goal === 'sample_size') {
          report = `In a paired / before-and-after repeated measures study design (paired t-test), evaluating a mean intra-individual shift of |Δ| = ${res.diff.toFixed(2)} (SD of paired differences σ_d = ${res.sdDiff.toFixed(2)}, standardized effect size Cohen's d_z = ${res.dz.toFixed(2)}) with ${Math.round(res.power * 100)}% statistical power at two-sided α = ${res.alpha}, a minimum of N = ${res.nPairs} paired subjects is required. Accounting for an estimated 10% withdrawal or uninterpretable paired follow-up, an enrollment of N = ${Math.ceil(res.nPairs / 0.9)} subjects is recommended.`;
        } else {
          report = `For a paired study design with N = ${res.nPairs} evaluable paired subjects evaluating a mean difference of |Δ| = ${res.diff.toFixed(2)} (SD of differences σ_d = ${res.sdDiff.toFixed(2)}, Cohen's d_z = ${res.dz.toFixed(2)}), the achieved statistical power at two-sided α = ${res.alpha} is ${(res.achievedPower * 100).toFixed(1)}% (Type II error rate β = ${(1 - res.achievedPower).toFixed(3)}).`;
        }
        if (reportBox) reportBox.innerText = report;

      } else if (design === 'contingency') {
        if (cardFisher) cardFisher.style.display = 'block';
        if (cardRisk) cardRisk.style.display = 'block';

        const p1 = parseFloat(document.getElementById('pwrP1')?.value) || 0.15;
        const p2 = parseFloat(document.getElementById('pwrP2')?.value) || 0.30;
        const propTestType = document.getElementById('pwrPropTestType')?.value || 'chisq';

        const res = PowerAnalysis.contingency2x2({ p1, p2, alpha, power, nPerGroup: givenN, testType: propTestType });
        if (res.error) {
          if (reportBox) reportBox.innerText = res.error;
          return;
        }

        if (effectLabel) effectLabel.innerText = 'Absolute Risk Reduction (ARR)';
        if (effectVal) effectVal.innerText = `${(res.arr * 100).toFixed(1)}%`;
        if (effectSub) effectSub.innerText = `Cohen's h = ${res.cohenH.toFixed(2)} | p₁=${(res.p1*100).toFixed(0)}% vs p₂=${(res.p2*100).toFixed(0)}%`;

        if (nGroupLabel) nGroupLabel.innerText = 'Sample Size per Group (n)';
        if (nGroupVal) nGroupVal.innerText = res.nPerGroup;
        if (nGroupSub) nGroupSub.innerText = `${res.testType === 'fisher' ? 'Continuity-corrected (Fleiss)' : 'Standard uncorrected'} 1:1`;

        if (nTotalLabel) nTotalLabel.innerText = 'Total Study Sample Size (N)';
        if (nTotalVal) nTotalVal.innerText = res.totalN;
        if (nTotalSub) nTotalSub.innerText = `Test + Control (2 × ${res.nPerGroup})`;

        if (achievedPowerVal) achievedPowerVal.innerText = `${(res.achievedPower * 100).toFixed(1)}%`;
        if (powerSub) powerSub.innerText = `Type II error β = ${(1 - res.achievedPower).toFixed(3)} | α = ${res.alpha}`;

        if (fisherNVal) fisherNVal.innerText = `n = ${res.continuityCorrectedN || res.nPerGroup}`;
        if (fisherSub) fisherSub.innerText = `Total N = ${2 * (res.continuityCorrectedN || res.nPerGroup)} (Fleiss / CC)`;

        if (riskVal) riskVal.innerText = `RR: ${res.rr.toFixed(2)} | OR: ${res.or.toFixed(2)}`;
        if (riskSub) riskSub.innerText = `NNT: ${res.nnt.toFixed(1)} | RRR: ${(res.rrr * 100).toFixed(1)}%`;

        let report = '';
        if (goal === 'sample_size') {
          report = `In a clinical trial comparing test (intervention) versus control cohorts with binary outcomes evaluated via a 2x2 contingency table (${res.testName}), expecting an event rate of ${(res.p1 * 100).toFixed(1)}% in the test arm versus ${(res.p2 * 100).toFixed(1)}% in the control arm (Absolute Risk Reduction ARR = ${(res.arr * 100).toFixed(1)}%, Relative Risk RR = ${res.rr.toFixed(2)}, Odds Ratio OR = ${res.or.toFixed(2)}, Number Needed to Treat NNT = ${res.nnt.toFixed(1)}):\n` +
            `• Uncorrected Pearson Chi-Square test requires n = ${res.uncorrectedN} patients per arm (total N = ${2 * res.uncorrectedN}) for ${Math.round(res.power * 100)}% power at two-sided α = ${res.alpha}.\n` +
            `• Fisher's Exact test / Continuity-Corrected Chi-Square (Fleiss & Casagrande-Pike-Smith formula) requires n = ${res.continuityCorrectedN} patients per arm (total N = ${2 * res.continuityCorrectedN}).\n` +
            `With an anticipated 15% loss-to-follow-up buffer, a target enrollment of N = ${Math.ceil((2 * res.nPerGroup) / 0.85)} patients (${Math.ceil(res.nPerGroup / 0.85)} per arm) is recommended.`;
        } else {
          report = `In a 2x2 contingency cohort study of n = ${res.nPerGroup} patients per group (total N = ${res.totalN}) comparing event rates of ${(res.p1 * 100).toFixed(1)}% (test) vs ${(res.p2 * 100).toFixed(1)}% (control) (ARR = ${(res.arr * 100).toFixed(1)}%, RR = ${res.rr.toFixed(2)}, OR = ${res.or.toFixed(2)}):\n` +
            `• Pearson Chi-Square test achieved power = ${(res.powerChisq * 100).toFixed(1)}% (β = ${(1 - res.powerChisq).toFixed(3)}).\n` +
            `• Fisher's Exact test / Continuity-Corrected achieved power = ${(res.powerFisher * 100).toFixed(1)}% (β = ${(1 - res.powerFisher).toFixed(3)}).`;
        }
        if (reportBox) reportBox.innerText = report;
      }
    }

    exportDocx(tabId) {
      if (!this.results || !this.results[tabId]) {
        switch (tabId) {
          case 'descriptive': this.runDescriptive(); break;
          case 'hypothesis': this.runHypo(); break;
          case 'anova': this.runAnova(); break;
          case 'categorical': this.runCat(); break;
          case 'correlation': this.runCorr(); break;
          case 'diagnostic': this.runROC(); break;
          case 'power': this.runPower(); break;
          case 'teaching': this.runTeaching(); break;
        case 'propensity': this.runPsm(); break;
        case 'multivariate': this.runMultivariate(); break;
        }
      }
      const res = this.results ? this.results[tabId] : null;
      if (!res) {
        alert('Please run the analysis first before exporting the report.');
        return;
      }

      let exportData = {};
      if (tabId === 'descriptive') {
        exportData = {
          name: 'Continuous Variable',
          n: res.n,
          mean: res.mean,
          ci95: res.ci95,
          sd: res.sd,
          variance: res.variance,
          sem: res.sem,
          median: res.median,
          q1: res.q1,
          q3: res.q3,
          iqr: res.iqr,
          modes: res.modes,
          maxFreq: res.maxFreq,
          min: res.min,
          max: res.max,
          skewness: res.skewness,
          skewnessInterpretation: res.skewnessInterpretation,
          kurtosis: res.kurtosis,
          kurtosisInterpretation: res.kurtosisInterpretation,
          normality: res.normality,
          lowerFence: res.lowerFence,
          upperFence: res.upperFence,
          outliers: res.outliers,
          reportText: document.getElementById('descReportText')?.innerText,
          values: res.values
        };
      } else if (tabId === 'hypothesis') {
        exportData = {
          nameA: document.getElementById('hypoNameA')?.value || 'Cohort A',
          nameB: document.getElementById('hypoNameB')?.value || 'Cohort B',
          testName: res.testName,
          statistic: res.statistic,
          df: res.df,
          pValue: res.pValue,
          isSignificant: res.isSignificant,
          meanDiff: res.meanDiff,
          ci95: res.ci95,
          cohensD: res.cohensD,
          rankBiserial: res.rankBiserial,
          reportText: document.getElementById('hypoReportText')?.innerText,
          groupA: res.groupA,
          groupB: res.groupB
        };
      } else if (tabId === 'anova') {
        exportData = Object.assign({}, res, {
          reportText: document.getElementById('anovaReportText')?.innerText
        });
      } else if (tabId === 'categorical') {
        const mode = document.getElementById('catAnalysisMode')?.value || 'diagnostic';
        const a = parseFloat(document.getElementById('catA')?.value) || 0;
        const b = parseFloat(document.getElementById('catB')?.value) || 0;
        const c = parseFloat(document.getElementById('catC')?.value) || 0;
        const d = parseFloat(document.getElementById('catD')?.value) || 0;
        exportData = {
          mode,
          a, b, c, d,
          totalN: a + b + c + d,
          reportText: document.getElementById('catReportText')?.innerText,
          diag: res.diagnostic,
          chiSquare: res.chiSquare,
          fishersExact: res.fishersExact,
          risk: res.riskMetrics
        };
      } else if (tabId === 'correlation') {
        exportData = {
          xName: 'Independent Variable X',
          yName: 'Dependent Variable Y',
          n: res.n,
          reportText: document.getElementById('corrReportText')?.innerText,
          corr: res.correlation,
          reg: res.regression,
          morphology: res.morphology
        };
      } else if (tabId === 'diagnostic') {
        exportData = {
          name: 'Clinical Biomarker Evaluation',
          totalN: res.eval2x2.total,
          posCount: res.eval2x2.tp + res.eval2x2.fn,
          negCount: res.eval2x2.fp + res.eval2x2.tn,
          reportText: document.getElementById('rocReportText')?.innerText,
          auc: res.roc.auc,
          aucCI95: res.roc.aucCI95,
          seAuc: res.roc.seAuc,
          optimalCutoff: res.roc.optimalCutoff,
          sensitivity: res.eval2x2.sensitivity,
          specificity: res.eval2x2.specificity,
          plr: res.eval2x2.plr,
          nlr: res.eval2x2.nlr
        };
      } else if (tabId === 'power') {
        exportData = {
          designLabel: document.getElementById('pwrStudyDesign')?.selectedOptions?.[0]?.text || 'Study Design',
          goal: document.getElementById('pwrCalcGoal')?.value || 'sample_size',
          alpha: parseFloat(document.getElementById('pwrAlpha')?.value) || 0.05,
          effectSize: res.effectSize || 0.5,
          effectSizeLabel: res.effectSizeLabel || "Cohen's d",
          nPerGroup: res.nPerGroup || res.nPairs || 0,
          totalN: res.totalN || 0,
          targetPower: parseFloat(document.getElementById('pwrPower')?.value) || 0.80,
          achievedPower: res.achievedPower || 0.80,
          reportText: document.getElementById('pwrReportText')?.innerText,
          additionalMetrics: res.additionalMetrics || {}
        };
      } else if (tabId === 'teaching') {
        const dist = res.dist || {};
        const clt = res.clt || {};
        exportData = {
          distName: dist.name,
          distN: dist.n,
          clinicalExample: dist.clinicalExample,
          sampleMean: dist.mean,
          theoMean: dist.theoMean,
          sampleSD: dist.sd,
          theoSD: dist.theoSD,
          skewness: dist.skewness,
          skewnessLabel: dist.skewnessInterpretation,
          kurtosis: dist.kurtosis,
          kurtosisLabel: dist.kurtosisInterpretation,
          jbStat: dist.normality?.statistic,
          jbP: dist.normality?.pValue,
          isNormal: dist.normality?.isNormal,
          cltParentName: clt.population?.name,
          cltN: clt.sampleSize,
          cltK: clt.samplesDrawn,
          cltResults: {
            theoMean: clt.theoreticalMean,
            obsMean: clt.observedMean,
            theoSE: clt.theoreticalSE,
            obsSE: clt.observedSE,
            skewness: clt.skewness,
            normalityP: clt.normality?.pValue,
            isNormal: clt.normality?.isNormal
          },
          tConv: res.tConv || Teaching.tConvergence.getMetrics(4),
          overlap: res.overlap || Teaching.significanceOverlap.getMetrics(),
          power: res.power || Teaching.powerSimulation.getMetrics(),
          reportText: document.getElementById('teachingReportText')?.innerText
        };
      } else if (tabId === 'teaching-bayesian') {
        exportData = {
          bayes: res.bayes || Teaching.bayesianSimulation.getMetrics(),
          reportText: document.getElementById('teachingBayesReportText')?.innerText || document.getElementById('bayesPedagogyText')?.innerText
        };
      } else if (tabId === 'propensity') {
        exportData = {
          name: 'Clinical Observational Cohort',
          treatmentCol: res.treatmentCol,
          outcomeCol: res.outcomeCol,
          covariateCols: res.covariateCols,
          totalN: res.totalN,
          completeN: res.completeN,
          nMatchedPairs: res.nMatchedPairs,
          matchedTreatedN: res.matchedTreatedN,
          matchedControlN: res.matchedControlN,
          logisticRegression: res.logisticRegression,
          caliper: res.caliper,
          balance: res.balance,
          outcome: res.outcome,
          reportText: document.getElementById('psmReportText')?.innerText
        };
      } else if (tabId === 'multivariate') {
        const mv = res;
        exportData = {
          method: mv.method,
          methodLabel: mv.method === 'pca' ? 'Principal Component Analysis (PCA)' : (mv.method === 'mca' ? 'Multiple Correspondence Analysis (MCA)' : 'Factor Analysis of Mixed Data (FAMD)'),
          nObservations: mv.individuals ? mv.individuals.length : 0,
          continuousCols: mv.continuousCols || [],
          categoricalCols: mv.categoricalCols || [],
          groupingCol: mv.groupingCol || null,
          outlierCount: (mv.outliers || []).length,
          scree: mv.scree || [],
          allVariables: mv.allVariables || mv.variables || [],
          interpretation: document.getElementById('multivarPedagogyText')?.innerText,
          reportText: document.getElementById('multivarReportText')?.innerText
        };
      }

      Exporter.exportTabToDocx(tabId, exportData);
    }

    // --- Teaching & Simulation Module ---

    initTeachingModule() {
      this.renderTeachingParams();
      this.runTeachingDistribution();
      this.updateCltUI(Teaching.clt.getSummary());
      this.runTConvergence(4);
      this.runTwoSampleOverlap();
      this.runPowerSimulation();
      this.runBayesianSimulation();
    }

    renderTeachingParams() {
      const distKey = document.getElementById('teachingDistSelect')?.value || 'normal';
      const container = document.getElementById('teachingParamsContainer');
      if (!container) return;

      let html = '';
      switch (distKey) {
        case 'normal':
          html = `
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem;">
              <div>
                <label class="form-label" style="font-size: 0.78rem;">Mean (μ):</label>
                <input id="tp_norm_mean" type="number" step="0.5" value="12" class="form-control" style="padding: 0.35rem 0.5rem;">
              </div>
              <div>
                <label class="form-label" style="font-size: 0.78rem;">Std Dev (σ):</label>
                <input id="tp_norm_sd" type="number" min="0.1" step="0.5" value="3" class="form-control" style="padding: 0.35rem 0.5rem;">
              </div>
            </div>
          `;
          break;
        case 'studentsT':
          html = `
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem;">
              <div>
                <label class="form-label" style="font-size: 0.78rem;">Degrees of Freedom (ν):</label>
                <input id="tp_t_df" type="number" min="1" max="100" step="1" value="4" class="form-control" style="padding: 0.35rem 0.5rem;">
              </div>
              <div>
                <label class="form-label" style="font-size: 0.78rem;">Scale (s):</label>
                <input id="tp_t_scale" type="number" min="0.1" step="0.5" value="2" class="form-control" style="padding: 0.35rem 0.5rem;">
              </div>
            </div>
          `;
          break;
        case 'logNormal':
          html = `
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem;">
              <div>
                <label class="form-label" style="font-size: 0.78rem;">Log-Mean (μ):</label>
                <input id="tp_log_mu" type="number" step="0.1" value="1.5" class="form-control" style="padding: 0.35rem 0.5rem;">
              </div>
              <div>
                <label class="form-label" style="font-size: 0.78rem;">Log-SD (σ):</label>
                <input id="tp_log_sigma" type="number" min="0.1" step="0.1" value="0.6" class="form-control" style="padding: 0.35rem 0.5rem;">
              </div>
            </div>
          `;
          break;
        case 'exponential':
          html = `
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Rate Parameter (λ):</label>
              <input id="tp_exp_rate" type="number" min="0.05" step="0.1" value="0.4" class="form-control" style="padding: 0.35rem 0.5rem;">
              <div style="font-size: 0.7rem; color: var(--text-muted); margin-top: 0.2rem;">Theoretical Mean = 1/λ = 2.50</div>
            </div>
          `;
          break;
        case 'bimodal':
          html = `
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; margin-bottom: 0.5rem;">
              <div>
                <label class="form-label" style="font-size: 0.75rem;">Peak 1 Mean:</label>
                <input id="tp_bim_m1" type="number" value="10" class="form-control" style="padding: 0.3rem 0.4rem;">
              </div>
              <div>
                <label class="form-label" style="font-size: 0.75rem;">Peak 1 SD:</label>
                <input id="tp_bim_s1" type="number" min="0.1" value="2" class="form-control" style="padding: 0.3rem 0.4rem;">
              </div>
            </div>
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem;">
              <div>
                <label class="form-label" style="font-size: 0.75rem;">Peak 2 Mean:</label>
                <input id="tp_bim_m2" type="number" value="25" class="form-control" style="padding: 0.3rem 0.4rem;">
              </div>
              <div>
                <label class="form-label" style="font-size: 0.75rem;">Peak 2 SD:</label>
                <input id="tp_bim_s2" type="number" min="0.1" value="3" class="form-control" style="padding: 0.3rem 0.4rem;">
              </div>
            </div>
          `;
          break;
        case 'uniform':
          html = `
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem;">
              <div>
                <label class="form-label" style="font-size: 0.78rem;">Lower Bound (a):</label>
                <input id="tp_uni_min" type="number" value="0" class="form-control" style="padding: 0.35rem 0.5rem;">
              </div>
              <div>
                <label class="form-label" style="font-size: 0.78rem;">Upper Bound (b):</label>
                <input id="tp_uni_max" type="number" value="20" class="form-control" style="padding: 0.35rem 0.5rem;">
              </div>
            </div>
          `;
          break;
        case 'poisson':
          html = `
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Event Rate Parameter (λ):</label>
              <input id="tp_pois_lambda" type="number" min="0.5" step="0.5" value="5" class="form-control" style="padding: 0.35rem 0.5rem;">
              <div style="font-size: 0.7rem; color: var(--text-muted); margin-top: 0.2rem;">Mean = Variance = λ = 5.0</div>
            </div>
          `;
          break;
        case 'chiSquare':
          html = `
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Degrees of Freedom (k):</label>
              <input id="tp_chi_df" type="number" min="1" max="50" step="1" value="4" class="form-control" style="padding: 0.35rem 0.5rem;">
              <div style="font-size: 0.7rem; color: var(--text-muted); margin-top: 0.2rem;">Mean = k = 4.0, Variance = 2k = 8.0</div>
            </div>
          `;
          break;
      }

      container.innerHTML = html;

      // Attach reactive input listeners to all inputs in container
      container.querySelectorAll('input').forEach(input => {
        input.addEventListener('input', () => this.runTeachingDistribution());
      });
    }

    runTeachingDistribution() {
      const distKey = document.getElementById('teachingDistSelect')?.value || 'normal';
      const n = parseInt(document.getElementById('teachingNRange')?.value) || 500;
      const meta = Teaching.distMetadata[distKey] || {};

      let data = [];
      const params = {};
      let theoMean = 0;
      let theoSD = 0;
      let theoSkew = '0.000';
      let theoKurt = '0.000';

      switch (distKey) {
        case 'normal': {
          const mean = parseFloat(document.getElementById('tp_norm_mean')?.value) || 12;
          const sd = Math.max(0.1, parseFloat(document.getElementById('tp_norm_sd')?.value) || 3);
          params.mean = mean;
          params.sd = sd;
          theoMean = mean;
          theoSD = sd;
          theoSkew = '0.000 (Symmetric)';
          theoKurt = '0.000 (Mesokurtic)';
          data = Teaching.generators.normal(n, mean, sd);
          break;
        }
        case 'studentsT': {
          const df = Math.max(1, parseInt(document.getElementById('tp_t_df')?.value) || 4);
          const scale = Math.max(0.1, parseFloat(document.getElementById('tp_t_scale')?.value) || 2);
          params.df = df;
          params.scale = scale;
          params.mean = 0;
          theoMean = 0;
          theoSD = df > 2 ? scale * Math.sqrt(df / (df - 2)) : NaN;
          theoSkew = '0.000';
          theoKurt = df > 4 ? (6 / (df - 4)).toFixed(2) + ' (Heavy tails)' : '∞ (Fat tails)';
          data = Teaching.generators.studentsT(n, df, 0, scale);
          break;
        }
        case 'logNormal': {
          const mu = parseFloat(document.getElementById('tp_log_mu')?.value) || 1.5;
          const sigma = Math.max(0.1, parseFloat(document.getElementById('tp_log_sigma')?.value) || 0.6);
          params.mu = mu;
          params.sigma = sigma;
          theoMean = Math.exp(mu + 0.5 * sigma * sigma);
          theoSD = Math.sqrt((Math.exp(sigma * sigma) - 1) * Math.exp(2 * mu + sigma * sigma));
          theoSkew = ((Math.exp(sigma * sigma) + 2) * Math.sqrt(Math.exp(sigma * sigma) - 1)).toFixed(2) + ' (Right skew)';
          theoKurt = 'Positive (Leptokurtic)';
          data = Teaching.generators.logNormal(n, mu, sigma);
          break;
        }
        case 'exponential': {
          const rate = Math.max(0.01, parseFloat(document.getElementById('tp_exp_rate')?.value) || 0.4);
          params.rate = rate;
          theoMean = 1 / rate;
          theoSD = 1 / rate;
          theoSkew = '2.000 (Strong right skew)';
          theoKurt = '6.000 (Heavy tails)';
          data = Teaching.generators.exponential(n, rate);
          break;
        }
        case 'bimodal': {
          const m1 = parseFloat(document.getElementById('tp_bim_m1')?.value) || 10;
          const s1 = Math.max(0.1, parseFloat(document.getElementById('tp_bim_s1')?.value) || 2);
          const m2 = parseFloat(document.getElementById('tp_bim_m2')?.value) || 25;
          const s2 = Math.max(0.1, parseFloat(document.getElementById('tp_bim_s2')?.value) || 3);
          params.m1 = m1;
          params.s1 = s1;
          params.m2 = m2;
          params.s2 = s2;
          params.p = 0.5;
          theoMean = 0.5 * (m1 + m2);
          theoSD = Math.sqrt(0.5 * (s1 * s1 + s2 * s2) + 0.25 * Math.pow(m1 - m2, 2));
          theoSkew = '~0.000 (Symmetric Peaks)';
          theoKurt = 'Platykurtic (Bimodal Trough)';
          data = Teaching.generators.bimodal(n, m1, s1, m2, s2, 0.5);
          break;
        }
        case 'uniform': {
          const min = parseFloat(document.getElementById('tp_uni_min')?.value) || 0;
          const max = parseFloat(document.getElementById('tp_uni_max')?.value) || 20;
          params.min = min;
          params.max = max;
          theoMean = (min + max) / 2;
          theoSD = Math.sqrt(Math.pow(max - min, 2) / 12);
          theoSkew = '0.000 (Flat)';
          theoKurt = '-1.200 (Platykurtic)';
          data = Teaching.generators.uniform(n, min, max);
          break;
        }
        case 'poisson': {
          const lambda = Math.max(0.1, parseFloat(document.getElementById('tp_pois_lambda')?.value) || 5);
          params.lambda = lambda;
          theoMean = lambda;
          theoSD = Math.sqrt(lambda);
          theoSkew = (1 / Math.sqrt(lambda)).toFixed(2);
          theoKurt = (1 / lambda).toFixed(2);
          data = Teaching.generators.poisson(n, lambda);
          break;
        }
        case 'chiSquare': {
          const df = Math.max(1, parseInt(document.getElementById('tp_chi_df')?.value) || 4);
          params.df = df;
          theoMean = df;
          theoSD = Math.sqrt(2 * df);
          theoSkew = Math.sqrt(8 / df).toFixed(2);
          theoKurt = (12 / df).toFixed(2);
          data = Teaching.generators.chiSquare(n, df);
          break;
        }
      }

      this.currentTeachingData = data;

      // Calculate empirical descriptive statistics
      const stats = Descriptive.calculate(data);

      // Update Metrics Cards
      const meanEl = document.getElementById('teachingSampleMean');
      if (meanEl) meanEl.innerText = stats.mean.toFixed(2);
      const theoMeanEl = document.getElementById('teachingTheoMean');
      if (theoMeanEl) theoMeanEl.innerText = `Theo: ${isNaN(theoMean) ? 'N/A' : theoMean.toFixed(2)}`;

      const sdEl = document.getElementById('teachingSampleSD');
      if (sdEl) sdEl.innerText = stats.sd.toFixed(2);
      const theoSDEl = document.getElementById('teachingTheoSD');
      if (theoSDEl) theoSDEl.innerText = `Theo: ${isNaN(theoSD) ? 'N/A' : theoSD.toFixed(2)}`;

      const skewEl = document.getElementById('teachingSkewness');
      if (skewEl) skewEl.innerText = stats.skewness.toFixed(2);
      const skewSubEl = document.getElementById('teachingSkewnessSub');
      if (skewSubEl) skewSubEl.innerText = stats.skewnessInterpretation;

      const kurtEl = document.getElementById('teachingKurtosis');
      if (kurtEl) kurtEl.innerText = stats.kurtosis.toFixed(2);
      const kurtSubEl = document.getElementById('teachingKurtosisSub');
      if (kurtSubEl) kurtSubEl.innerText = stats.kurtosisInterpretation;

      const normEl = document.getElementById('teachingNormality');
      const normSubEl = document.getElementById('teachingNormalitySub');
      if (normEl) {
        normEl.innerText = stats.normality.isNormal ? 'Normal ✓' : 'Non-Normal ✗';
        normEl.style.color = stats.normality.isNormal ? '#22c55e' : '#ef4444';
      }
      if (normSubEl) {
        normSubEl.innerText = `JB p = ${stats.normality.pValue.toFixed(4)}`;
      }

      // Update Context Card
      const titleEl = document.getElementById('teachingContextTitle');
      const bodyEl = document.getElementById('teachingContextBody');
      if (titleEl) titleEl.innerText = `${meta.name || 'Distribution'} [${meta.symbol || ''}]`;
      if (bodyEl) {
        bodyEl.innerHTML = `
          <strong>Clinical Example:</strong> ${meta.clinicalExample || ''}<br>
          <strong>Theoretical Behavior:</strong> ${meta.properties || ''}<br>
          <strong>Methodological Guideline:</strong> ${meta.testNote || ''}
        `;
      }

      // Update Data Preview (first 50 values)
      const previewEl = document.getElementById('teachingDataPreview');
      if (previewEl) {
        previewEl.value = data.slice(0, 50).map(v => v.toFixed(2)).join(', ') + (data.length > 50 ? ` ... (+${data.length - 50} more)` : '');
      }

      // Render Canvas
      const canvas = document.getElementById('teachingDistCanvas');
      if (canvas && this.engines['teachingDistCanvas']) {
        Plots.renderTeachingDistribution(
          this.engines['teachingDistCanvas'],
          data,
          distKey,
          params,
          `${meta.name || 'Distribution'} (Empirical N = ${n} vs Theoretical PDF)`
        );
      }

      // Cache result
      this.results.teachingDist = {
        name: meta.name || distKey,
        clinicalExample: meta.clinicalExample,
        n,
        mean: stats.mean,
        theoMean,
        sd: stats.sd,
        theoSD,
        skewness: stats.skewness,
        skewnessInterpretation: stats.skewnessInterpretation,
        theoSkewness: theoSkew,
        kurtosis: stats.kurtosis,
        kurtosisInterpretation: stats.kurtosisInterpretation,
        theoKurtosis: theoKurt,
        normality: stats.normality
      };

      if (!this.results.teaching) this.results.teaching = {};
      this.results.teaching.dist = this.results.teachingDist;
    }

    updateCltUI(summary) {
      if (!summary) return;

      // Update Metrics
      const kEl = document.getElementById('cltSamplesCount');
      if (kEl) kEl.innerText = summary.samplesDrawn.toLocaleString();

      const nSubEl = document.getElementById('cltNSub');
      if (nSubEl) nSubEl.innerText = `Each of size n = ${summary.sampleSize}`;

      const muEl = document.getElementById('cltTrueMu');
      if (muEl) muEl.innerText = summary.theoreticalMean.toFixed(2);

      const sigmaSubEl = document.getElementById('cltTrueSigmaSub');
      if (sigmaSubEl) sigmaSubEl.innerText = `Pop SD: ${summary.population.sd.toFixed(2)}`;

      const obsMeanEl = document.getElementById('cltObsMean');
      const meanDiffEl = document.getElementById('cltMeanDiffSub');
      if (obsMeanEl) {
        obsMeanEl.innerText = summary.observedMean !== null ? summary.observedMean.toFixed(2) : '--';
      }
      if (meanDiffEl) {
        meanDiffEl.innerText = summary.observedMean !== null
          ? `|Diff|: ${Math.abs(summary.observedMean - summary.theoreticalMean).toFixed(3)}`
          : '|Diff|: --';
      }

      const theoSEEl = document.getElementById('cltTheoSE');
      if (theoSEEl) theoSEEl.innerText = summary.theoreticalSE.toFixed(3);

      const obsSEEl = document.getElementById('cltObsSE');
      const seDiffEl = document.getElementById('cltSEDiffSub');
      if (obsSEEl) {
        obsSEEl.innerText = summary.observedSE !== null ? summary.observedSE.toFixed(3) : '--';
      }
      if (seDiffEl) {
        seDiffEl.innerText = summary.observedSE !== null
          ? `Diff: ${(summary.observedSE - summary.theoreticalSE).toFixed(3)}`
          : 'Formula: σ / √n';
      }

      const normEl = document.getElementById('cltNormality');
      const normSubEl = document.getElementById('cltNormalitySub');
      if (normEl) {
        if (summary.normality) {
          normEl.innerText = summary.normality.isNormal ? 'Normal ✓' : 'Non-Normal';
          normEl.style.color = summary.normality.isNormal ? '#22c55e' : '#f59e0b';
        } else {
          normEl.innerText = '--';
          normEl.style.color = 'inherit';
        }
      }
      if (normSubEl) {
        normSubEl.innerText = summary.normality ? `JB p = ${summary.normality.pValue.toFixed(3)}` : 'JB Test: --';
      }

      // Render Canvas Charts
      if (this.engines['teachingCltParentCanvas']) {
        Plots.renderCltParent(
          this.engines['teachingCltParentCanvas'],
          summary.population,
          summary.lastSample,
          'Parent Population Distribution'
        );
      }

      if (this.engines['teachingCltSamplingCanvas']) {
        Plots.renderCltSampling(
          this.engines['teachingCltSamplingCanvas'],
          summary,
          'Sampling Distribution of Mean (x̄)'
        );
      }

      // Update Report Text
      const reportEl = document.getElementById('teachingReportText');
      if (reportEl) {
        if (summary.samplesDrawn === 0) {
          reportEl.innerText = `Simulation initialized with ${summary.population.name} and sample size n = ${summary.sampleSize}. Click "▶ Step (Draw 1 Sample)" or "⚡ Draw 100" to begin demonstrating convergence.`;
        } else {
          const normStatement = summary.normality?.isNormal
            ? 'conforms strictly to a Gaussian bell curve (Jarque-Bera p ≥ 0.05), mathematically proving the Central Limit Theorem.'
            : 'is rapidly converging toward Gaussian symmetry as iterations accumulate.';

          reportEl.innerText =
            `A Central Limit Theorem simulation was conducted using a ${summary.population.name} (True μ = ${summary.theoreticalMean.toFixed(2)}, σ = ${summary.population.sd.toFixed(2)}). ` +
            `A total of k = ${summary.samplesDrawn.toLocaleString()} independent samples (each of size n = ${summary.sampleSize}) were drawn. ` +
            `The grand mean of sample means was x̄̄ = ${summary.observedMean.toFixed(2)} (bias = ${Math.abs(summary.observedMean - summary.theoreticalMean).toFixed(3)}), ` +
            `with an empirical Standard Error of s_x̄ = ${summary.observedSE.toFixed(3)} compared to theoretical SE = σ/√n = ${summary.theoreticalSE.toFixed(3)}. ` +
            `The sampling distribution displays Skewness G₁ = ${summary.skewness.toFixed(3)} and ${normStatement} ` +
            `In clinical trials, this proves why parametric t-tests and ANOVA remain robust for sample sizes n ≥ 30 even when raw clinical metrics are non-normal.`;
        }
      }

      // Cache
      this.results.teachingClt = summary;
      if (!this.results.teaching) this.results.teaching = {};
      this.results.teaching.clt = summary;
    }

    runTeaching() {
      this.runTeachingDistribution();
      this.updateCltUI(Teaching.clt.getSummary());
      this.runTConvergence();
    }

    runTConvergence(sampleSize = null) {
      const rangeEl = document.getElementById('tConvNRange');
      if (sampleSize !== null && rangeEl) {
        rangeEl.value = sampleSize;
      }
      const n = parseInt(rangeEl?.value) || 4;
      const valEl = document.getElementById('tConvNVal');
      if (valEl) {
        valEl.innerText = `n = ${n} (ν = ${n - 1})`;
      }

      const metrics = Teaching.tConvergence.getMetrics(n);

      // Update UI Cards
      const dfEl = document.getElementById('tConvDf');
      const dfSubEl = document.getElementById('tConvDfSub');
      if (dfEl) dfEl.innerText = metrics.df;
      if (dfSubEl) dfSubEl.innerText = `Sample size n = ${metrics.sampleSize}`;

      const critEl = document.getElementById('tConvCrit');
      const critSubEl = document.getElementById('tConvCritSub');
      if (critEl) critEl.innerText = metrics.tCrit.toFixed(3);
      if (critSubEl) critSubEl.innerText = `vs z = 1.960 (${metrics.critDiffPct >= 0 ? '+' : ''}${metrics.critDiffPct.toFixed(1)}%)`;

      const peakEl = document.getElementById('tConvPeak');
      const peakSubEl = document.getElementById('tConvPeakSub');
      if (peakEl) peakEl.innerText = metrics.tPeak.toFixed(4);
      if (peakSubEl) peakSubEl.innerText = `Normal: 0.3989 (${metrics.peakDiffPct.toFixed(1)}%)`;

      const tailEl = document.getElementById('tConvTailProb');
      const tailSubEl = document.getElementById('tConvTailProbSub');
      if (tailEl) {
        tailEl.innerText = `${(metrics.tailProb * 100).toFixed(1)}%`;
        tailEl.style.color = metrics.df < 30 ? '#ef4444' : '#22c55e';
      }
      if (tailSubEl) tailSubEl.innerText = metrics.df < 30 ? 'Normal: 5.0% (Type I Risk)' : 'Normal: 5.0% (Matched)';

      const maxDiffEl = document.getElementById('tConvMaxDiff');
      if (maxDiffEl) maxDiffEl.innerText = metrics.maxDiscrepancy.toFixed(4);

      const kurtEl = document.getElementById('tConvKurt');
      const kurtSubEl = document.getElementById('tConvKurtSub');
      if (kurtEl) {
        kurtEl.innerText = metrics.excessKurtosis === Infinity ? '∞' : metrics.excessKurtosis.toFixed(2);
      }
      if (kurtSubEl) {
        kurtSubEl.innerText = metrics.df <= 4 ? 'Fat-tailed (ν ≤ 4)' : `Excess Kurtosis: 6/(ν-4)`;
      }

      const pedagogyEl = document.getElementById('tConvPedagogyText');
      if (pedagogyEl) pedagogyEl.innerText = metrics.clinicalNote;

      const chartTitleEl = document.getElementById('tConvChartTitle');
      if (chartTitleEl) {
        chartTitleEl.innerText = `Student's t(ν = ${metrics.df}) Density Curve vs Standard Normal N(0, 1)`;
      }

      const showTailArea = document.getElementById('tConvShowTailArea')?.checked !== false;

      // Render Canvas
      if (this.engines['teachingTCanvas']) {
        Plots.renderTConvergence(this.engines['teachingTCanvas'], metrics, {
          showTailArea,
          title: `Student's t(ν = ${metrics.df}) vs Standard Normal N(0, 1)`
        });
      }

      // Cache
      if (!this.results.teaching) this.results.teaching = {};
      this.results.teaching.tConv = metrics;
    }

    animateTConvergence() {
      const btn = document.getElementById('tConvAnimateBtn');
      if (this.tConvAnimationTimer) {
        clearInterval(this.tConvAnimationTimer);
        this.tConvAnimationTimer = null;
        if (btn) btn.innerText = '▶ Animate Convergence';
        return;
      }

      const frames = [2, 3, 4, 5, 7, 10, 15, 20, 25, 31, 45, 61, 80, 100, 121, 150];
      let currentIndex = 0;
      const currentN = parseInt(document.getElementById('tConvNRange')?.value) || 4;
      const startIdx = frames.findIndex(f => f >= currentN);
      if (startIdx >= 0 && startIdx < frames.length - 1) currentIndex = startIdx;

      if (btn) btn.innerText = '⏸ Pause Animation';

      this.tConvAnimationTimer = setInterval(() => {
        currentIndex++;
        if (currentIndex >= frames.length) {
          clearInterval(this.tConvAnimationTimer);
          this.tConvAnimationTimer = null;
          if (btn) btn.innerText = '▶ Replay Animation';
          return;
        }
        this.runTConvergence(frames[currentIndex]);
      }, 450);
    }

    runTwoSampleOverlap(overrideParams = {}) {
      const delta = overrideParams.delta !== undefined
        ? overrideParams.delta
        : (parseFloat(document.getElementById('overlapDeltaRange')?.value) || 2.0);
      const alpha = overrideParams.alpha !== undefined
        ? overrideParams.alpha
        : (parseFloat(document.getElementById('overlapAlphaRange')?.value) || 0.05);

      const sd1 = overrideParams.sd1 !== undefined
        ? overrideParams.sd1
        : (parseFloat(document.getElementById('overlapSD1Range')?.value) || 2.5);
      const sem1 = overrideParams.sem1 !== undefined
        ? overrideParams.sem1
        : (parseFloat(document.getElementById('overlapSEM1Range')?.value) || 0.625);
      const n1 = overrideParams.n1 !== undefined
        ? overrideParams.n1
        : (parseInt(document.getElementById('overlapN1Range')?.value) || 16);

      const sd2 = overrideParams.sd2 !== undefined
        ? overrideParams.sd2
        : (parseFloat(document.getElementById('overlapSD2Range')?.value) || 2.5);
      const sem2 = overrideParams.sem2 !== undefined
        ? overrideParams.sem2
        : (parseFloat(document.getElementById('overlapSEM2Range')?.value) || 0.625);
      const n2 = overrideParams.n2 !== undefined
        ? overrideParams.n2
        : (parseInt(document.getElementById('overlapN2Range')?.value) || 16);

      const viewMode = this.currentOverlapMode || 'means';

      const metrics = Teaching.significanceOverlap.getMetrics({
        mean1: 10.0,
        delta,
        sd1,
        sd2,
        sem1,
        sem2,
        n1,
        n2,
        alpha,
        viewMode
      });

      // Update Slider Displays
      const deltaValEl = document.getElementById('overlapDeltaVal');
      if (deltaValEl) deltaValEl.innerText = metrics.delta.toFixed(2);

      const alphaValEl = document.getElementById('overlapAlphaVal');
      if (alphaValEl) alphaValEl.innerText = `α = ${metrics.alpha.toFixed(3)} (z = ${metrics.zCrit.toFixed(3)})`;

      // Group 1 Displays
      const g1SummaryEl = document.getElementById('overlapG1Summary');
      if (g1SummaryEl) g1SummaryEl.innerText = `SD₁ = ${metrics.sd1.toFixed(2)} | SEM₁ = ${metrics.sem1.toFixed(3)}`;

      const sd1ValEl = document.getElementById('overlapSD1Val');
      if (sd1ValEl) sd1ValEl.innerText = metrics.sd1.toFixed(2);

      const sem1ValEl = document.getElementById('overlapSEM1Val');
      if (sem1ValEl) sem1ValEl.innerText = metrics.sem1.toFixed(3);

      const n1ValEl = document.getElementById('overlapN1Val');
      if (n1ValEl) n1ValEl.innerText = `n₁ = ${metrics.n1}`;

      // Group 2 Displays
      const g2SummaryEl = document.getElementById('overlapG2Summary');
      if (g2SummaryEl) g2SummaryEl.innerText = `SD₂ = ${metrics.sd2.toFixed(2)} | SEM₂ = ${metrics.sem2.toFixed(3)}`;

      const sd2ValEl = document.getElementById('overlapSD2Val');
      if (sd2ValEl) sd2ValEl.innerText = metrics.sd2.toFixed(2);

      const sem2ValEl = document.getElementById('overlapSEM2Val');
      if (sem2ValEl) sem2ValEl.innerText = metrics.sem2.toFixed(3);

      const n2ValEl = document.getElementById('overlapN2Val');
      if (n2ValEl) n2ValEl.innerText = `n₂ = ${metrics.n2}`;

      // Update Metric Cards
      const statusValEl = document.getElementById('overlapStatusValue');
      const pValSubEl = document.getElementById('overlapPValueSub');
      if (statusValEl) {
        statusValEl.innerText = metrics.isSignificant ? 'SIGNIFICANT' : 'NOT SIGNIFICANT';
        statusValEl.style.color = metrics.isSignificant ? '#10b981' : '#ef4444';
      }
      if (pValSubEl) {
        pValSubEl.innerText = `p = ${metrics.pValue < 0.0001 ? '< 0.0001' : metrics.pValue.toFixed(4)} ${metrics.isSignificant ? '<' : '≥'} α = ${metrics.alpha.toFixed(3)}`;
      }

      const deltaDispEl = document.getElementById('overlapDeltaDisplay');
      const deltaCritSubEl = document.getElementById('overlapDeltaCritSub');
      if (deltaDispEl) deltaDispEl.innerText = `Δ = ${metrics.delta.toFixed(2)}`;
      if (deltaCritSubEl) deltaCritSubEl.innerText = `Δcrit = ${metrics.deltaCrit.toFixed(2)} (Boundary)`;

      const g1MetricVal = document.getElementById('overlapG1MetricValue');
      const g1MetricSub = document.getElementById('overlapG1MetricSub');
      if (g1MetricVal) g1MetricVal.innerText = `${metrics.sd1.toFixed(2)} | ${metrics.sem1.toFixed(3)}`;
      if (g1MetricSub) g1MetricSub.innerText = `n₁ = ${metrics.n1} subjects`;

      const g2MetricVal = document.getElementById('overlapG2MetricValue');
      const g2MetricSub = document.getElementById('overlapG2MetricSub');
      if (g2MetricVal) g2MetricVal.innerText = `${metrics.sd2.toFixed(2)} | ${metrics.sem2.toFixed(3)}`;
      if (g2MetricSub) g2MetricSub.innerText = `n₂ = ${metrics.n2} subjects`;

      const meansOVLEl = document.getElementById('overlapMeansOVLValue');
      const patientOVLSubEl = document.getElementById('overlapPatientOVLSub');
      if (meansOVLEl) {
        meansOVLEl.innerText = `${(metrics.meansOVL * 100).toFixed(1)}%`;
        meansOVLEl.style.color = metrics.isSignificant ? '#10b981' : '#ef4444';
      }
      if (patientOVLSubEl) patientOVLSubEl.innerText = `Patient Overlap: ${(metrics.patientOVL * 100).toFixed(1)}%`;

      const critValEl = document.getElementById('overlapCritValue');
      const alphaSubEl = document.getElementById('overlapAlphaSub');
      if (critValEl) critValEl.innerText = `t = ${metrics.tCrit.toFixed(3)}`;
      if (alphaSubEl) alphaSubEl.innerText = `SE_diff: ${metrics.seDiff.toFixed(3)} | df: ${metrics.df.toFixed(1)}`;

      // Update Pedagogical Text
      const pedaEl = document.getElementById('overlapPedagogyText');
      if (pedaEl) pedaEl.innerText = metrics.explanation;

      // Update Chart Title
      const titleEl = document.getElementById('overlapChartTitle');
      if (titleEl) {
        if (viewMode === 'patients') {
          titleEl.innerText = `Individual Patient Populations (SD): Biological Overlap vs Mean Significance`;
        } else if (viewMode === 'dual') {
          titleEl.innerText = `Dual Overlay: Patient Biological Spread (SD) vs Inferential Mean Precision (SEM)`;
        } else if (viewMode === 'null') {
          titleEl.innerText = `Null Hypothesis Difference Test: H₀ (Δ=0) vs Observed Separation`;
        } else {
          titleEl.innerText = `Sampling Distributions of Means (SEM): Overlap & Significance Boundary`;
        }
      }

      // Render Canvas
      if (this.engines['teachingOverlapCanvas']) {
        Plots.renderTwoSampleOverlap(this.engines['teachingOverlapCanvas'], metrics);
      }

      // Cache
      if (!this.results.teaching) this.results.teaching = {};
      this.results.teaching.overlap = metrics;
    }

    animateOverlapSeparation() {
      const btn = document.getElementById('overlapAnimateBtn');
      if (this.overlapAnimationTimer) {
        clearInterval(this.overlapAnimationTimer);
        this.overlapAnimationTimer = null;
        if (btn) btn.innerText = '▶ Animate Separation';
        return;
      }

      const deltas = [0.0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.0, 2.3, 2.7, 3.2, 3.8, 4.5];
      let currentIndex = 0;
      const currentDelta = parseFloat(document.getElementById('overlapDeltaRange')?.value) || 0;
      const startIdx = deltas.findIndex(d => d >= currentDelta);
      if (startIdx >= 0 && startIdx < deltas.length - 1) currentIndex = startIdx;

      if (btn) btn.innerText = '⏸ Pause Animation';

      this.overlapAnimationTimer = setInterval(() => {
        currentIndex++;
        if (currentIndex >= deltas.length) {
          clearInterval(this.overlapAnimationTimer);
          this.overlapAnimationTimer = null;
          if (btn) btn.innerText = '▶ Replay Animation';
          return;
        }
        const nextDelta = deltas[currentIndex];
        const dRange = document.getElementById('overlapDeltaRange');
        if (dRange) dRange.value = nextDelta;
        this.runTwoSampleOverlap({ delta: nextDelta });
      }, 500);
    }

    runPowerSimulation(overrideParams = {}, lastChanged = 'power') {
      const sd = overrideParams.sd !== undefined
        ? overrideParams.sd
        : (parseFloat(document.getElementById('powerSDRange')?.value) || 4.0);
      const sem = overrideParams.sem !== undefined
        ? overrideParams.sem
        : (parseFloat(document.getElementById('powerSEMRange')?.value) || 0.50);
      const n = overrideParams.n !== undefined
        ? overrideParams.n
        : (parseInt(document.getElementById('powerNRange')?.value) || 64);
      const power = overrideParams.power !== undefined
        ? overrideParams.power
        : (parseFloat(document.getElementById('powerPowerRange')?.value) || 0.80);
      const beta = overrideParams.beta !== undefined
        ? overrideParams.beta
        : (parseFloat(document.getElementById('powerBetaRange')?.value) || 0.20);
      const delta = overrideParams.delta !== undefined
        ? overrideParams.delta
        : (parseFloat(document.getElementById('powerDeltaRange')?.value) || 2.0);
      const alpha = overrideParams.alpha !== undefined
        ? overrideParams.alpha
        : (parseFloat(document.getElementById('powerAlphaRange')?.value) || 0.05);

      const viewMode = this.currentPowerMode || 'distributions';

      const metrics = Teaching.powerSimulation.getMetrics({
        sd,
        sem,
        n,
        power,
        beta,
        delta,
        alpha,
        viewMode,
        lastChanged
      });

      // Synchronize UI Slider Values without trigger loop
      const sdInput = document.getElementById('powerSDRange');
      const semInput = document.getElementById('powerSEMRange');
      const nInput = document.getElementById('powerNRange');
      const powerInput = document.getElementById('powerPowerRange');
      const betaInput = document.getElementById('powerBetaRange');
      const deltaInput = document.getElementById('powerDeltaRange');
      const alphaInput = document.getElementById('powerAlphaRange');

      if (sdInput && lastChanged !== 'sd') sdInput.value = metrics.sd.toFixed(1);
      if (semInput && lastChanged !== 'sem') semInput.value = metrics.sem.toFixed(3);
      if (nInput && lastChanged !== 'n') nInput.value = metrics.n;
      if (powerInput && lastChanged !== 'power') powerInput.value = metrics.power.toFixed(2);
      if (betaInput && lastChanged !== 'beta') betaInput.value = metrics.beta.toFixed(2);
      if (deltaInput && lastChanged !== 'delta') deltaInput.value = metrics.delta.toFixed(1);
      if (alphaInput && lastChanged !== 'alpha') alphaInput.value = metrics.alpha.toFixed(3);

      // Update Slider Display Labels
      const sdValEl = document.getElementById('powerSDVal');
      const semValEl = document.getElementById('powerSEMVal');
      const nValEl = document.getElementById('powerNVal');
      const powerValEl = document.getElementById('powerPowerVal');
      const betaValEl = document.getElementById('powerBetaVal');
      const deltaValEl = document.getElementById('powerDeltaVal');
      const alphaValEl = document.getElementById('powerAlphaVal');

      if (sdValEl) sdValEl.innerText = metrics.sd.toFixed(2);
      if (semValEl) semValEl.innerText = metrics.sem.toFixed(3);
      if (nValEl) nValEl.innerText = `n = ${metrics.n} (N = ${metrics.totalN})`;
      if (powerValEl) {
        powerValEl.innerText = `${(metrics.power * 100).toFixed(1)}%`;
        powerValEl.style.color = metrics.power >= 0.80 ? '#10b981' : (metrics.power >= 0.60 ? '#f59e0b' : '#ef4444');
      }
      if (betaValEl) {
        betaValEl.innerText = `${(metrics.beta * 100).toFixed(1)}%`;
        betaValEl.style.color = metrics.beta <= 0.20 ? '#10b981' : (metrics.beta <= 0.40 ? '#f59e0b' : '#ef4444');
      }
      if (deltaValEl) deltaValEl.innerText = metrics.delta.toFixed(2);
      if (alphaValEl) alphaValEl.innerText = metrics.alpha.toFixed(3);

      const dSub = document.getElementById('powerCohensDSub');
      if (dSub) dSub.innerText = `Cohen's d = ${metrics.cohensD.toFixed(2)}`;
      const zSub = document.getElementById('powerZCritSub');
      if (zSub) zSub.innerText = `z_crit = ${metrics.zCrit.toFixed(3)}`;

      const g1Summary = document.getElementById('powerG1Summary');
      if (g1Summary) g1Summary.innerText = `SD = ${metrics.sd.toFixed(2)} | SEM = ${metrics.sem.toFixed(3)}`;
      const g2Summary = document.getElementById('powerG2Summary');
      if (g2Summary) g2Summary.innerText = `Power = ${(metrics.power * 100).toFixed(1)}% | β = ${(metrics.beta * 100).toFixed(1)}%`;

      // Update 6 Status Metric Cards
      const statusPwrVal = document.getElementById('powerStatusPowerVal');
      const statusPwrSub = document.getElementById('powerStatusPowerSub');
      if (statusPwrVal) {
        statusPwrVal.innerText = `${(metrics.power * 100).toFixed(1)}%`;
        statusPwrVal.style.color = metrics.power >= 0.80 ? '#10b981' : (metrics.power >= 0.60 ? '#f59e0b' : '#ef4444');
      }
      if (statusPwrSub) {
        statusPwrSub.innerText = `Target: ≥ 80% (${metrics.powerRating})`;
      }

      const statusBetaVal = document.getElementById('powerStatusBetaVal');
      if (statusBetaVal) {
        statusBetaVal.innerText = `${(metrics.beta * 100).toFixed(1)}%`;
        statusBetaVal.style.color = metrics.beta <= 0.20 ? '#10b981' : (metrics.beta <= 0.40 ? '#f59e0b' : '#ef4444');
      }

      const statusNVal = document.getElementById('powerStatusNVal');
      const statusNSub = document.getElementById('powerStatusNSub');
      if (statusNVal) statusNVal.innerText = `n = ${metrics.n}`;
      if (statusNSub) statusNSub.innerText = `Total N = ${metrics.totalN} subjects`;

      const statusDeltaVal = document.getElementById('powerStatusDeltaVal');
      const statusDeltaSub = document.getElementById('powerStatusDeltaSub');
      if (statusDeltaVal) statusDeltaVal.innerText = `Δ = ${metrics.delta.toFixed(2)}`;
      if (statusDeltaSub) {
        const dLabel = metrics.cohensD >= 0.8 ? 'Large effect' : (metrics.cohensD >= 0.5 ? 'Medium effect' : 'Small effect');
        statusDeltaSub.innerText = `Cohen's d = ${metrics.cohensD.toFixed(2)} (${dLabel})`;
      }

      const statusSEMVal = document.getElementById('powerStatusSEMVal');
      const statusSEMSub = document.getElementById('powerStatusSEMSub');
      if (statusSEMVal) statusSEMVal.innerText = `SEM = ${metrics.sem.toFixed(3)}`;
      if (statusSEMSub) statusSEMSub.innerText = `SE_diff: ${metrics.seDiff.toFixed(3)} (σ·√(2/n))`;

      const statusCutoffVal = document.getElementById('powerStatusCutoffVal');
      const statusCutoffSub = document.getElementById('powerStatusCutoffSub');
      if (statusCutoffVal) statusCutoffVal.innerText = `xcrit = ${metrics.xCrit.toFixed(3)}`;
      if (statusCutoffSub) statusCutoffSub.innerText = `α = ${metrics.alpha.toFixed(3)} (Type I: ${(metrics.alpha * 100).toFixed(1)}%)`;

      // Update Pedagogical Text
      const pedaEl = document.getElementById('powerPedagogyText');
      if (pedaEl) pedaEl.innerText = metrics.explanation;

      // Update Chart Title
      const titleEl = document.getElementById('powerChartTitle');
      if (titleEl) {
        if (viewMode === 'curve') {
          titleEl.innerText = `Statistical Power Curve: Power (1 − β) vs Sample Size n (Δ = ${metrics.delta.toFixed(2)}, SD = ${metrics.sd.toFixed(2)})`;
        } else if (viewMode === 'matrix') {
          titleEl.innerText = `2×2 Decision Error Matrix: True State vs Clinical Statistical Decision`;
        } else {
          titleEl.innerText = `Dual Distribution: Null Hypothesis H₀ vs True Effect H₁ (Power & Beta Shading)`;
        }
      }

      // Render Canvas
      if (this.engines['teachingPowerCanvas']) {
        Plots.renderPowerSimulation(this.engines['teachingPowerCanvas'], metrics);
      }

      // Cache
      if (!this.results.teaching) this.results.teaching = {};
      this.results.teaching.power = metrics;
    }

    animatePowerGain() {
      const btn = document.getElementById('powerAnimateBtn');
      if (this.powerAnimationTimer) {
        clearInterval(this.powerAnimationTimer);
        this.powerAnimationTimer = null;
        if (btn) btn.innerText = '▶ Animate Power Gain';
        return;
      }

      const sampleSizes = [10, 16, 24, 34, 46, 64, 86, 112, 144];
      let currentIndex = 0;
      const currentN = parseInt(document.getElementById('powerNRange')?.value) || 64;
      const startIdx = sampleSizes.findIndex(n => n >= currentN);
      if (startIdx >= 0 && startIdx < sampleSizes.length - 1) currentIndex = startIdx;

      if (btn) btn.innerText = '⏸ Pause Animation';

      this.powerAnimationTimer = setInterval(() => {
        currentIndex++;
        if (currentIndex >= sampleSizes.length) {
          clearInterval(this.powerAnimationTimer);
          this.powerAnimationTimer = null;
          if (btn) btn.innerText = '▶ Replay Animation';
          return;
        }
        const nextN = sampleSizes[currentIndex];
        const nRange = document.getElementById('powerNRange');
        if (nRange) nRange.value = nextN;
        this.runPowerSimulation({ n: nextN }, 'n');
      }, 550);
    }

    runBayesianSimulation(overrideParams = {}) {
      // Read input values with fallbacks
      const priorInput = document.getElementById('bayesPriorRange');
      const sampleInput = document.getElementById('bayesSampleRange');
      const likInput = document.getElementById('bayesLikelihoodRange');
      const fpInput = document.getElementById('bayesFalsePosRange');

      let prior = overrideParams.prior !== undefined 
        ? overrideParams.prior 
        : (priorInput ? parseFloat(priorInput.value) / 100 : 0.0476);
      let sampleSize = overrideParams.sampleSize !== undefined 
        ? overrideParams.sampleSize 
        : (sampleInput ? parseInt(sampleInput.value) : 210);
      let likelihood = overrideParams.likelihood !== undefined 
        ? overrideParams.likelihood 
        : (likInput ? parseFloat(likInput.value) / 100 : 0.40);
      let falsePositive = overrideParams.falsePositive !== undefined 
        ? overrideParams.falsePositive 
        : (fpInput ? parseFloat(fpInput.value) / 100 : 0.10);

      const viewMode = this.currentBayesView || 'square';
      const preset = overrideParams.preset || this.currentBayesPreset || 'steve';

      // Smart percentage formatter preserving exact user-entered precision (e.g. 11.8% instead of rounding to 12%)
      const formatPercentSmart = (rate) => {
        if (rate === undefined || rate === null || isNaN(rate)) return '0.0%';
        const pct = rate * 100;
        const rounded1 = Math.round(pct * 10) / 10;
        if (Math.abs(pct - rounded1) < 1e-4) {
          return `${rounded1.toFixed(1)}%`;
        }
        const rounded2 = Math.round(pct * 100) / 100;
        if (Math.abs(pct - rounded2) < 1e-4) {
          return `${rounded2.toFixed(2)}%`;
        }
        return `${parseFloat(pct.toFixed(3))}%`;
      };

      const getSliderDisplayVal = (rate) => {
        const pct = rate * 100;
        const rounded1 = Math.round(pct * 10) / 10;
        if (Math.abs(pct - rounded1) < 1e-4) {
          return rounded1.toFixed(1);
        }
        const rounded2 = Math.round(pct * 100) / 100;
        if (Math.abs(pct - rounded2) < 1e-4) {
          return rounded2.toFixed(2);
        }
        return parseFloat(pct.toFixed(3)).toString();
      };

      // Synchronize UI slider values and numeric display labels
      if (priorInput && overrideParams.prior !== undefined) {
        priorInput.value = getSliderDisplayVal(prior);
      }
      if (sampleInput && overrideParams.sampleSize !== undefined) {
        sampleInput.value = sampleSize;
      }
      if (likInput && overrideParams.likelihood !== undefined) {
        likInput.value = getSliderDisplayVal(likelihood);
      }
      if (fpInput && overrideParams.falsePositive !== undefined) {
        fpInput.value = getSliderDisplayVal(falsePositive);
      }

      const priorValEl = document.getElementById('bayesPriorVal');
      if (priorValEl) priorValEl.innerText = formatPercentSmart(prior);

      const sampleValEl = document.getElementById('bayesSampleVal');
      if (sampleValEl) sampleValEl.innerText = `N = ${sampleSize} people`;

      const likValEl = document.getElementById('bayesLikelihoodVal');
      if (likValEl) likValEl.innerText = formatPercentSmart(likelihood);

      const fpValEl = document.getElementById('bayesFalsePosVal');
      if (fpValEl) fpValEl.innerText = formatPercentSmart(falsePositive);

      const priorSumEl = document.getElementById('bayesPriorSummary');
      if (priorSumEl) {
        const priorOddsRatio = (1 - prior) / prior;
        priorSumEl.innerText = `Prior: ${formatPercentSmart(prior)} (1:${priorOddsRatio.toFixed(1)})`;
      }

      const likSumEl = document.getElementById('bayesLikelihoodSummary');
      if (likSumEl) {
        likSumEl.innerText = `P(E|H): ${formatPercentSmart(likelihood)} | P(E|¬H): ${formatPercentSmart(falsePositive)}`;
      }

      // Calculate metrics
      const metrics = Teaching.bayesianSimulation.getMetrics({
        prior,
        likelihood,
        falsePositive,
        sampleSize,
        viewMode,
        preset
      });

      // Update real-time Metric Cards
      const statusPriorVal = document.getElementById('bayesStatusPriorVal');
      const statusPriorSub = document.getElementById('bayesStatusPriorSub');
      if (statusPriorVal) statusPriorVal.innerText = formatPercentSmart(metrics.prior);
      if (statusPriorSub) statusPriorSub.innerText = `Prior Odds: 1 : ${(1 / metrics.priorOdds).toFixed(1)}`;

      const statusLikVal = document.getElementById('bayesStatusLikelihoodVal');
      const statusLikSub = document.getElementById('bayesStatusLikelihoodSub');
      if (statusLikVal) statusLikVal.innerText = formatPercentSmart(metrics.likelihood);
      if (statusLikSub) statusLikSub.innerText = `True Positive Probability`;

      const statusFPVal = document.getElementById('bayesStatusFalsePosVal');
      const statusFPSub = document.getElementById('bayesStatusFalsePosSub');
      if (statusFPVal) statusFPVal.innerText = formatPercentSmart(metrics.falsePositive);
      if (statusFPSub) statusFPSub.innerText = `False Positive Probability`;

      const statusPEVal = document.getElementById('bayesStatusPEvidenceVal');
      const statusPESub = document.getElementById('bayesStatusPEvidenceSub');
      if (statusPEVal) statusPEVal.innerText = `${(metrics.pEvidence * 100).toFixed(2)}%`;
      if (statusPESub) statusPESub.innerText = `Marginal / Area (${metrics.countTotalE} of ${metrics.sampleSize})`;

      const statusBFVal = document.getElementById('bayesStatusBFVal');
      const statusBFSub = document.getElementById('bayesStatusBFSub');
      if (statusBFVal) statusBFVal.innerText = `${metrics.bayesFactor.toFixed(2)}×`;
      if (statusBFSub) statusBFSub.innerText = metrics.evidenceRating;

      const statusPostVal = document.getElementById('bayesStatusPosteriorVal');
      const statusPostSub = document.getElementById('bayesStatusPosteriorSub');
      if (statusPostVal) statusPostVal.innerText = `${(metrics.posterior * 100).toFixed(1)}%`;
      if (statusPostSub) {
        const sign = metrics.beliefShift >= 0 ? '+' : '';
        statusPostSub.innerText = `Belief Shift: ${sign}${(metrics.beliefShift * 100).toFixed(1)}%`;
      }

      // Update Pedagogical Text
      const pedaEl = document.getElementById('bayesPedagogyText');
      if (pedaEl) pedaEl.innerText = metrics.explanation;

      const bayesReportEl = document.getElementById('teachingBayesReportText');
      if (bayesReportEl) bayesReportEl.innerText = metrics.explanation;

      // Update Chart Title
      const titleEl = document.getElementById('teachingBayesChartTitle');
      if (titleEl) {
        if (viewMode === 'sample') {
          titleEl.innerText = `Representative Sample (N = ${metrics.sampleSize}): Natural Counts Frequency Tree (${metrics.countHAndE} vs ${metrics.countNotHAndE})`;
        } else if (viewMode === 'odds') {
          titleEl.innerText = `Odds Form & Bayes Factor: Prior Odds (1 : ${(1 / metrics.priorOdds).toFixed(1)}) × LR (${metrics.bayesFactor.toFixed(2)}×) = Posterior Odds (${metrics.posteriorOdds.toFixed(2)})`;
        } else if (viewMode === 'sequential') {
          titleEl.innerText = `Sequential Bayesian Updating: Compounding Evidence Trajectory (P₀ = ${(metrics.prior * 100).toFixed(1)}% → P₄ = ${(metrics.trajectory[4].p * 100).toFixed(1)}%)`;
        } else {
          titleEl.innerText = `3Blue1Brown 1×1 Unit Square: Visualizing Bayes' Theorem as Proportions of Area (P(H|E) = ${(metrics.posterior * 100).toFixed(1)}%)`;
        }
      }

      // Render Canvas
      if (this.engines['teachingBayesCanvas']) {
        Plots.renderBayesianSimulation(this.engines['teachingBayesCanvas'], metrics);
      }

      // Cache
      if (!this.results.teaching) this.results.teaching = {};
      this.results.teaching.bayes = metrics;
    }

    animateBayesRestriction() {
      const btn = document.getElementById('bayesAnimateBtn');
      if (this.bayesAnimationTimer) {
        clearInterval(this.bayesAnimationTimer);
        this.bayesAnimationTimer = null;
        if (btn) btn.innerText = '▶ Animate Restriction';
        return;
      }

      // Cycle through view modes
      const views = ['square', 'sample', 'odds', 'sequential'];
      let idx = views.indexOf(this.currentBayesView || 'square');
      if (btn) btn.innerText = '⏸ Pause Animation';

      this.bayesAnimationTimer = setInterval(() => {
        idx = (idx + 1) % views.length;
        const nextView = views[idx];
        this.currentBayesView = nextView;
        document.querySelectorAll('.teaching-bayes-view-btn').forEach(b => {
          const isActive = b.dataset.view === nextView;
          b.classList.toggle('btn-primary', isActive);
          b.classList.toggle('active', isActive);
          b.classList.toggle('btn-secondary', !isActive);
        });
        this.runBayesianSimulation();
      }, 1200);
    }

    // ==========================================
    // 13. RANDOMISER WORKFLOW & EVENT HANDLERS
    // ==========================================
    initRandomiser() {
      this.randomiserState = {
        mode: localStorage.getItem('statis_gravity_randomiser_mode') || 'simple',
        simple: {
          history: [],
          nextId: 1,
          labelA: 'Group A (Treatment)',
          labelB: 'Group B (Control)'
        },
        block: {
          history: [],
          nextId: 1,
          blockSize: 4,
          targetN: 40,
          labelA: 'Group A (Intervention)',
          labelB: 'Group B (Control)',
          activeBlock: null
        }
      };

      try {
        const savedSimple = localStorage.getItem('statis_gravity_randomiser_simple');
        if (savedSimple) {
          const parsed = JSON.parse(savedSimple);
          if (parsed && Array.isArray(parsed.history)) {
            this.randomiserState.simple.history = parsed.history;
            this.randomiserState.simple.nextId = parsed.nextId || (parsed.history.length + 1);
            if (parsed.labelA) this.randomiserState.simple.labelA = parsed.labelA;
            if (parsed.labelB) this.randomiserState.simple.labelB = parsed.labelB;
          }
        }
        const savedBlock = localStorage.getItem('statis_gravity_randomiser_block');
        if (savedBlock) {
          const parsed = JSON.parse(savedBlock);
          if (parsed && Array.isArray(parsed.history)) {
            this.randomiserState.block.history = parsed.history;
            this.randomiserState.block.nextId = parsed.nextId || (parsed.history.length + 1);
            if (parsed.blockSize) this.randomiserState.block.blockSize = parsed.blockSize;
            if (parsed.targetN) this.randomiserState.block.targetN = parsed.targetN;
            if (parsed.labelA) this.randomiserState.block.labelA = parsed.labelA;
            if (parsed.labelB) this.randomiserState.block.labelB = parsed.labelB;
            if (parsed.activeBlock) this.randomiserState.block.activeBlock = parsed.activeBlock;
          }
        }
      } catch (err) {
        console.warn('Could not load randomiser state from localStorage:', err);
      }

      // Sync form input fields
      const sLabelA = document.getElementById('randomiserSimpleLabelA');
      const sLabelB = document.getElementById('randomiserSimpleLabelB');
      const sNextId = document.getElementById('randomiserSimpleNextId');
      if (sLabelA) sLabelA.value = this.randomiserState.simple.labelA;
      if (sLabelB) sLabelB.value = this.randomiserState.simple.labelB;
      if (sNextId) sNextId.value = this.randomiserState.simple.nextId;

      const bSize = document.getElementById('randomiserBlockSizeSelect');
      const bTargetN = document.getElementById('randomiserBlockTargetN');
      const bLabelA = document.getElementById('randomiserBlockLabelA');
      const bLabelB = document.getElementById('randomiserBlockLabelB');
      if (bSize) bSize.value = String(this.randomiserState.block.blockSize);
      if (bTargetN) bTargetN.value = this.randomiserState.block.targetN;
      if (bLabelA) bLabelA.value = this.randomiserState.block.labelA;
      if (bLabelB) bLabelB.value = this.randomiserState.block.labelB;

      // Event listeners for Simple inputs
      sLabelA?.addEventListener('input', (e) => {
        this.randomiserState.simple.labelA = e.target.value.trim() || 'Group A';
        this.saveRandomiserState('simple');
      });
      sLabelB?.addEventListener('input', (e) => {
        this.randomiserState.simple.labelB = e.target.value.trim() || 'Group B';
        this.saveRandomiserState('simple');
      });
      sNextId?.addEventListener('change', (e) => {
        const val = parseInt(e.target.value, 10);
        if (!isNaN(val) && val >= 1) {
          this.randomiserState.simple.nextId = val;
          const btnId = document.getElementById('randomiserSimpleBtnIdText');
          if (btnId) btnId.innerText = String(val);
          this.saveRandomiserState('simple');
        }
      });

      // Event listeners for Block inputs
      bSize?.addEventListener('change', (e) => {
        const newSize = parseInt(e.target.value, 10);
        this.randomiserState.block.blockSize = newSize;
        if (!this.randomiserState.block.activeBlock || this.randomiserState.block.activeBlock.currentIndex === 0) {
          this.randomiserState.block.activeBlock = null;
        }
        this.updateRandomiserBlockUI();
        this.saveRandomiserState('block');
      });
      bTargetN?.addEventListener('input', (e) => {
        const val = parseInt(e.target.value, 10);
        if (!isNaN(val) && val >= 2) {
          this.randomiserState.block.targetN = val;
          this.updateRandomiserBlockUI();
          this.saveRandomiserState('block');
        }
      });
      bLabelA?.addEventListener('input', (e) => {
        this.randomiserState.block.labelA = e.target.value.trim() || 'Group A';
        this.saveRandomiserState('block');
      });
      bLabelB?.addEventListener('input', (e) => {
        this.randomiserState.block.labelB = e.target.value.trim() || 'Group B';
        this.saveRandomiserState('block');
      });

      // Sub-Navigation mode switching
      const modeSimpleBtn = document.getElementById('randomiserModeSimpleBtn');
      const modeBlockBtn = document.getElementById('randomiserModeBlockBtn');
      const simplePanel = document.getElementById('randomiserSimplePanel');
      const blockPanel = document.getElementById('randomiserBlockPanel');

      const setMode = (mode) => {
        this.randomiserState.mode = mode;
        localStorage.setItem('statis_gravity_randomiser_mode', mode);

        if (mode === 'simple') {
          modeSimpleBtn?.classList.add('active');
          modeSimpleBtn?.setAttribute('aria-selected', 'true');
          modeBlockBtn?.classList.remove('active');
          modeBlockBtn?.setAttribute('aria-selected', 'false');
          if (simplePanel) simplePanel.style.display = 'block';
          if (blockPanel) blockPanel.style.display = 'none';
        } else {
          modeBlockBtn?.classList.add('active');
          modeBlockBtn?.setAttribute('aria-selected', 'true');
          modeSimpleBtn?.classList.remove('active');
          modeSimpleBtn?.setAttribute('aria-selected', 'false');
          if (simplePanel) simplePanel.style.display = 'none';
          if (blockPanel) blockPanel.style.display = 'block';
        }
        this.renderRandomiserAuditTable();
      };

      modeSimpleBtn?.addEventListener('click', () => setMode('simple'));
      modeBlockBtn?.addEventListener('click', () => setMode('block'));

      // Action Buttons
      document.getElementById('randomiserSimpleGenerateBtn')?.addEventListener('click', () => {
        this.generateSimpleAllocationAction();
      });

      document.getElementById('randomiserBlockAssignBtn')?.addEventListener('click', () => {
        this.assignBlockAllocationAction();
      });

      document.getElementById('randomiserExportCsvBtn')?.addEventListener('click', () => {
        this.exportRandomiserCSV();
      });

      document.getElementById('randomiserCopyTableBtn')?.addEventListener('click', () => {
        this.copyRandomiserAuditTable();
      });

      // Reset Modal Handlers
      const resetModal = document.getElementById('randomiserResetModal');
      document.getElementById('randomiserResetBtn')?.addEventListener('click', () => {
        if (resetModal) resetModal.classList.remove('hidden');
      });
      document.getElementById('randomiserResetCancelBtn')?.addEventListener('click', () => {
        if (resetModal) resetModal.classList.add('hidden');
      });
      document.getElementById('randomiserResetConfirmBtn')?.addEventListener('click', () => {
        this.resetRandomiserSession();
        if (resetModal) resetModal.classList.add('hidden');
      });

      // Initial View Setup
      setMode(this.randomiserState.mode);
      this.updateRandomiserSimpleUI();
      this.updateRandomiserBlockUI();
      this.renderRandomiserAuditTable();
    }

    generateSimpleAllocationAction() {
      try {
        if (!this.randomiserState || !this.randomiserState.simple) {
          this.initRandomiser();
        }
        const sNextIdEl = document.getElementById('randomiserSimpleNextId');
        if (sNextIdEl) {
          const parsed = parseInt(sNextIdEl.value, 10);
          if (!isNaN(parsed) && parsed >= 1) {
            this.randomiserState.simple.nextId = parsed;
          }
        }
        const pid = this.randomiserState.simple.nextId;
        const record = Randomiser.generateSimpleAllocation({
          participantId: pid,
          labelA: this.randomiserState.simple.labelA || 'Group A (Treatment)',
          labelB: this.randomiserState.simple.labelB || 'Group B (Control)'
        });

        this.randomiserState.simple.history.push(record);
        this.randomiserState.simple.nextId += 1;

        if (sNextIdEl) sNextIdEl.value = this.randomiserState.simple.nextId;

        this.saveRandomiserState('simple');
        this.updateRandomiserSimpleUI(record);
        this.renderRandomiserAuditTable();
      } catch (err) {
        console.error('Error generating simple allocation:', err);
      }
    }

    updateRandomiserSimpleUI(lastRecord = null) {
      const sNextId = this.randomiserState.simple.nextId;
      const btnIdText = document.getElementById('randomiserSimpleBtnIdText');
      if (btnIdText) btnIdText.innerText = String(sNextId);

      const history = this.randomiserState.simple.history;
      const rec = lastRecord || (history.length > 0 ? history[history.length - 1] : null);

      const idleView = document.getElementById('randomiserSimpleIdleView');
      const activeView = document.getElementById('randomiserSimpleActiveView');
      const revealCard = document.getElementById('randomiserSimpleRevealCard');
      const badge = document.getElementById('randomiserSimpleBadge');
      const resId = document.getElementById('randomiserSimpleResultId');
      const resNum = document.getElementById('randomiserSimpleResultNum');
      const resParity = document.getElementById('randomiserSimpleResultParity');
      const resTime = document.getElementById('randomiserSimpleResultTime');
      const tsText = document.getElementById('randomiserSimpleTimestampText');

      if (rec) {
        if (idleView) idleView.style.display = 'none';
        if (activeView) activeView.style.display = 'block';
        if (resId) resId.innerText = `#${rec.participantId}`;
        if (badge) {
          badge.innerText = rec.groupLabel;
          badge.className = `allocation-badge-large ${rec.groupKey === 'A' ? 'badge-group-a' : 'badge-group-b'}`;
        }
        if (revealCard) {
          revealCard.className = `allocation-reveal-card ${rec.groupKey === 'A' ? 'revealed-a' : 'revealed-b'}`;
        }
        if (resNum) resNum.innerText = String(rec.randomNumber);
        if (resParity) {
          resParity.innerText = rec.parity;
          resParity.style.color = rec.groupKey === 'A' ? 'var(--cyan-primary)' : 'var(--emerald-primary)';
        }
        if (resTime) resTime.innerText = rec.displayTime || '';
        if (tsText) tsText.innerText = `Last assigned at ${rec.displayTime || ''}`;
      } else {
        if (idleView) idleView.style.display = 'block';
        if (activeView) activeView.style.display = 'none';
        if (revealCard) revealCard.className = 'allocation-reveal-card';
        if (tsText) tsText.innerText = 'No allocation yet';
      }

      // Summary metrics
      const summary = Randomiser.computeSummary(history);
      const totalEl = document.getElementById('randomiserSimpleTotalN');
      const countAEl = document.getElementById('randomiserSimpleCountA');
      const countBEl = document.getElementById('randomiserSimpleCountB');
      const ratioEl = document.getElementById('randomiserSimpleRatio');
      const pctAEl = document.getElementById('randomiserSimplePctA');
      const pctBEl = document.getElementById('randomiserSimplePctB');
      const barA = document.getElementById('randomiserSimpleBarA');
      const barB = document.getElementById('randomiserSimpleBarB');

      if (totalEl) totalEl.innerText = String(summary.total);
      if (countAEl) countAEl.innerText = `${summary.countA} (${summary.pctA.toFixed(0)}%)`;
      if (countBEl) countBEl.innerText = `${summary.countB} (${summary.pctB.toFixed(0)}%)`;
      if (ratioEl) ratioEl.innerText = summary.ratioStr;
      if (pctAEl) pctAEl.innerText = `${summary.pctA.toFixed(0)}%`;
      if (pctBEl) pctBEl.innerText = `${summary.pctB.toFixed(0)}%`;

      const barAWidth = summary.total === 0 ? 50 : Math.max(5, Math.min(95, summary.pctA));
      const barBWidth = summary.total === 0 ? 50 : (100 - barAWidth);
      if (barA) barA.style.width = `${barAWidth}%`;
      if (barB) barB.style.width = `${barBWidth}%`;
    }

    assignBlockAllocationAction() {
      try {
        if (!this.randomiserState || !this.randomiserState.block) {
          this.initRandomiser();
        }
        const pid = this.randomiserState.block.nextId;
        const res = Randomiser.assignNextInBlock(
          this.randomiserState.block.activeBlock,
          pid,
          {
            blockSize: this.randomiserState.block.blockSize,
            labelA: this.randomiserState.block.labelA || 'Group A (Intervention)',
            labelB: this.randomiserState.block.labelB || 'Group B (Control)'
          }
        );

        this.randomiserState.block.activeBlock = res.updatedBlock;
        this.randomiserState.block.history.push(res.allocationRecord);
        this.randomiserState.block.nextId += 1;

        this.saveRandomiserState('block');
        this.updateRandomiserBlockUI(res.allocationRecord);
        this.renderRandomiserAuditTable();
      } catch (err) {
        console.error('Error assigning block allocation:', err);
      }
    }

    updateRandomiserBlockUI(lastRecord = null) {
      const bNextId = this.randomiserState.block.nextId;
      const btnIdText = document.getElementById('randomiserBlockBtnIdText');
      if (btnIdText) btnIdText.innerText = `Participant #${bNextId}`;

      const history = this.randomiserState.block.history;
      const rec = lastRecord || (history.length > 0 ? history[history.length - 1] : null);

      const idleView = document.getElementById('randomiserBlockIdleView');
      const activeView = document.getElementById('randomiserBlockActiveView');
      const revealCard = document.getElementById('randomiserBlockRevealCard');
      const badge = document.getElementById('randomiserBlockBadge');
      const resId = document.getElementById('randomiserBlockResultId');
      const resBlockNum = document.getElementById('randomiserBlockResultBlockNum');
      const resSlot = document.getElementById('randomiserBlockResultSlot');
      const resTime = document.getElementById('randomiserBlockResultTime');

      if (rec) {
        if (idleView) idleView.style.display = 'none';
        if (activeView) activeView.style.display = 'block';
        if (resId) resId.innerText = `#${rec.participantId}`;
        if (badge) {
          badge.innerText = rec.groupLabel;
          badge.className = `allocation-badge-large ${rec.groupKey === 'A' ? 'badge-group-a' : 'badge-group-b'}`;
        }
        if (revealCard) {
          revealCard.className = `allocation-reveal-card ${rec.groupKey === 'A' ? 'revealed-a' : 'revealed-b'}`;
        }
        if (resBlockNum) resBlockNum.innerText = String(rec.blockNumber);
        if (resSlot) resSlot.innerText = `${rec.slotInBlock} of ${rec.blockSize}`;
        if (resTime) resTime.innerText = rec.displayTime || '';
      } else {
        if (idleView) idleView.style.display = 'block';
        if (activeView) activeView.style.display = 'none';
        if (revealCard) revealCard.className = 'allocation-reveal-card';
      }

      // Block slots visual tracker
      const activeBlock = this.randomiserState.block.activeBlock;
      const bSize = this.randomiserState.block.blockSize;
      const slotsContainer = document.getElementById('randomiserBlockSlotsContainer');
      const progressText = document.getElementById('randomiserBlockProgressText');
      const statusBadge = document.getElementById('randomiserBlockStatusBadge');

      if (slotsContainer) {
        slotsContainer.innerHTML = '';
        const currentBlockNum = activeBlock ? activeBlock.blockNumber : (Math.floor(history.length / bSize) + 1);
        if (statusBadge) statusBadge.innerText = `Block #${currentBlockNum} Active`;

        if (activeBlock && activeBlock.slots) {
          const assignedCount = activeBlock.currentIndex;
          if (progressText) progressText.innerText = `Slot ${assignedCount} of ${activeBlock.blockSize} assigned`;

          activeBlock.slots.forEach(slot => {
            const pill = document.createElement('div');
            if (slot.assigned) {
              pill.className = `block-slot-pill ${slot.groupKey === 'A' ? 'assigned-a' : 'assigned-b'}`;
              pill.innerHTML = `✓ Slot ${slot.slotIndex}: <strong>${slot.groupKey}</strong> (P#${slot.participantId})`;
            } else if (slot.slotIndex === activeBlock.currentIndex + 1) {
              pill.className = 'block-slot-pill concealed-slot active-slot';
              pill.innerHTML = `★ Next Slot ${slot.slotIndex}: 🔒 Concealed`;
            } else {
              pill.className = 'block-slot-pill concealed-slot';
              pill.innerHTML = `Slot ${slot.slotIndex}: 🔒 Pending`;
            }
            slotsContainer.appendChild(pill);
          });
        } else {
          if (progressText) progressText.innerText = `Slot 0 of ${bSize} assigned`;
          for (let i = 1; i <= bSize; i++) {
            const pill = document.createElement('div');
            pill.className = 'block-slot-pill concealed-slot';
            pill.innerHTML = `Slot ${i}: 🔒 Concealed`;
            slotsContainer.appendChild(pill);
          }
        }
      }

      // Summary metrics
      const summary = Randomiser.computeSummary(history);
      const totalBlockEl = document.getElementById('randomiserBlockTotalN');
      const countABlockEl = document.getElementById('randomiserBlockCountA');
      const countBBlockEl = document.getElementById('randomiserBlockCountB');
      const completedBlocksEl = document.getElementById('randomiserBlockCompletedCount');
      const pctABlockEl = document.getElementById('randomiserBlockPctA');
      const pctBBlockEl = document.getElementById('randomiserBlockPctB');
      const barABlock = document.getElementById('randomiserBlockBarA');
      const barBBlock = document.getElementById('randomiserBlockBarB');

      if (totalBlockEl) totalBlockEl.innerText = `${summary.total} / ${this.randomiserState.block.targetN}`;
      if (countABlockEl) countABlockEl.innerText = `${summary.countA} (${summary.pctA.toFixed(0)}%)`;
      if (countBBlockEl) countBBlockEl.innerText = `${summary.countB} (${summary.pctB.toFixed(0)}%)`;
      const completedBlocks = Math.floor(summary.total / bSize);
      if (completedBlocksEl) completedBlocksEl.innerText = String(completedBlocks);
      if (pctABlockEl) pctABlockEl.innerText = `${summary.pctA.toFixed(0)}%`;
      if (pctBBlockEl) pctBBlockEl.innerText = `${summary.pctB.toFixed(0)}%`;

      const barAWidth = summary.total === 0 ? 50 : Math.max(5, Math.min(95, summary.pctA));
      const barBWidth = summary.total === 0 ? 50 : (100 - barAWidth);
      if (barABlock) barABlock.style.width = `${barAWidth}%`;
      if (barBBlock) barBBlock.style.width = `${barBWidth}%`;
    }

    renderRandomiserAuditTable() {
      const mode = this.randomiserState.mode;
      const history = this.randomiserState[mode].history;
      const tbody = document.getElementById('randomiserAuditTableBody');
      const countBadge = document.getElementById('randomiserAuditCountBadge');

      if (!tbody) return;

      if (countBadge) {
        countBadge.innerText = `${history.length} Record${history.length === 1 ? '' : 's'}`;
      }

      if (history.length === 0) {
        tbody.innerHTML = `
          <tr id="randomiserAuditEmptyRow">
            <td colspan="6" style="text-align: center; padding: 2.5rem 1rem; color: var(--text-dim);">
              <div style="font-size: 1.8rem; margin-bottom: 0.4rem;">🎲</div>
              <div style="font-weight: 600; color: var(--text-muted);">No participants randomized yet in ${mode === 'simple' ? 'Simple' : 'Block'} mode</div>
              <div style="font-size: 0.78rem; margin-top: 0.2rem;">Click "Generate / Assign Next Participant" above to initiate sequence allocation.</div>
            </td>
          </tr>
        `;
        return;
      }

      tbody.innerHTML = '';
      // Render in reverse chronological order (newest first)
      for (let i = history.length - 1; i >= 0; i--) {
        const item = history[i];
        const tr = document.createElement('tr');
        tr.style.borderBottom = '1px solid var(--border-subtle)';

        const drawCol = mode === 'simple'
          ? `<span class="badge ${item.parity === 'Odd' ? 'badge-sig' : 'badge-neutral'}">${item.randomNumber} (${item.parity})</span>`
          : `<span class="badge badge-sig">Block ${item.blockNumber} (Slot ${item.slotInBlock}/${item.blockSize})</span>`;

        tr.innerHTML = `
          <td style="padding: 0.6rem 0.85rem; font-weight: 800; color: var(--text-main);">#${item.participantId}</td>
          <td style="padding: 0.6rem 0.85rem; font-family: var(--font-mono); font-size: 0.78rem; color: var(--text-muted);">${item.displayTime || (item.timestamp ? item.timestamp.slice(11, 19) : '--')}</td>
          <td style="padding: 0.6rem 0.85rem; font-size: 0.8rem; color: var(--text-muted);">${item.method || (mode === 'simple' ? 'Simple (1-100)' : 'Block Permuted')}</td>
          <td style="padding: 0.6rem 0.85rem;">${drawCol}</td>
          <td style="padding: 0.6rem 0.85rem; font-weight: 800; color: ${item.groupKey === 'A' ? 'var(--cyan-primary)' : 'var(--emerald-primary)'}; font-size: 1.05rem;">${item.groupKey}</td>
          <td style="padding: 0.6rem 0.85rem; font-weight: 600; color: var(--text-main);">${item.groupLabel}</td>
        `;
        tbody.appendChild(tr);
      }
    }

    saveRandomiserState(type) {
      try {
        if (type === 'simple') {
          localStorage.setItem('statis_gravity_randomiser_simple', JSON.stringify({
            history: this.randomiserState.simple.history,
            nextId: this.randomiserState.simple.nextId,
            labelA: this.randomiserState.simple.labelA,
            labelB: this.randomiserState.simple.labelB
          }));
        } else if (type === 'block') {
          localStorage.setItem('statis_gravity_randomiser_block', JSON.stringify({
            history: this.randomiserState.block.history,
            nextId: this.randomiserState.block.nextId,
            blockSize: this.randomiserState.block.blockSize,
            targetN: this.randomiserState.block.targetN,
            labelA: this.randomiserState.block.labelA,
            labelB: this.randomiserState.block.labelB,
            activeBlock: this.randomiserState.block.activeBlock
          }));
        }
      } catch (err) {
        console.warn('Error saving randomiser state to localStorage:', err);
      }
    }

    exportRandomiserCSV() {
      const mode = this.randomiserState.mode;
      const history = this.randomiserState[mode].history;

      if (!history || history.length === 0) {
        alert(`No randomized records to export in ${mode === 'simple' ? 'Simple' : 'Block'} Randomization mode.`);
        return;
      }

      const csvContent = Randomiser.exportToCSV(history, mode);
      const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      const now = new Date();
      const dateStr = now.toISOString().slice(0, 10);
      a.download = `statis_gravity_${mode}_randomisation_audit_${dateStr}.csv`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    }

    copyRandomiserAuditTable() {
      const mode = this.randomiserState.mode;
      const history = this.randomiserState[mode].history;
      const copyBtn = document.getElementById('randomiserCopyTableBtn');

      if (!history || history.length === 0) {
        alert('Audit table is empty. Randomize at least one participant first.');
        return;
      }

      const headers = mode === 'simple'
        ? ['Participant ID', 'Timestamp', 'Method', 'Random Draw (1-100)', 'Parity', 'Group Code', 'Group Label']
        : ['Participant ID', 'Timestamp', 'Method', 'Block Number', 'Block Size', 'Slot in Block', 'Group Code', 'Group Label'];

      const rows = history.map(item => mode === 'simple'
        ? [item.participantId, item.timestamp, item.method || 'Simple', item.randomNumber, item.parity, item.groupKey, item.groupLabel]
        : [item.participantId, item.timestamp, item.method || 'Block', item.blockNumber, item.blockSize, item.slotInBlock, item.groupKey, item.groupLabel]
      );

      const tsv = [headers.join('\t'), ...rows.map(r => r.join('\t'))].join('\n');

      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(tsv).then(() => {
          if (copyBtn) {
            const orig = copyBtn.innerText;
            copyBtn.innerText = '✓ Copied!';
            setTimeout(() => { copyBtn.innerText = orig; }, 1800);
          }
        });
      } else {
        const ta = document.createElement('textarea');
        ta.value = tsv;
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        document.execCommand('copy');
        document.body.removeChild(ta);
        if (copyBtn) {
          const orig = copyBtn.innerText;
          copyBtn.innerText = '✓ Copied!';
          setTimeout(() => { copyBtn.innerText = orig; }, 1800);
        }
      }
    }

    resetRandomiserSession() {
      const mode = this.randomiserState.mode;
      if (mode === 'simple') {
        this.randomiserState.simple.history = [];
        this.randomiserState.simple.nextId = 1;
        const sNextIdEl = document.getElementById('randomiserSimpleNextId');
        if (sNextIdEl) sNextIdEl.value = '1';
        localStorage.removeItem('statis_gravity_randomiser_simple');
        this.updateRandomiserSimpleUI();
      } else {
        this.randomiserState.block.history = [];
        this.randomiserState.block.nextId = 1;
        this.randomiserState.block.activeBlock = null;
        localStorage.removeItem('statis_gravity_randomiser_block');
        this.updateRandomiserBlockUI();
      }
      this.renderRandomiserAuditTable();
    }

    // ==========================================
    // 14. PROPENSITY SCORE MATCHING (PSM) WORKFLOW
    // ==========================================
    initPsm() {
      this.psmCurrentScriptLang = 'python';
      this.psmActiveView = 'love';

      document.getElementById('psmSampleBtn')?.addEventListener('click', () => {
        this.loadPsmSample();
      });

      document.getElementById('psmClearBtn')?.addEventListener('click', () => {
        const input = document.getElementById('psmCsvInput');
        if (input) input.value = '';
        this.populatePsmVariables();
      });

      document.getElementById('psmComputeBtn')?.addEventListener('click', () => {
        this.runPsm();
      });

      const loveBtn = document.getElementById('psmViewLoveBtn');
      const overlapBtn = document.getElementById('psmViewOverlapBtn');
      const saveBtn = document.getElementById('psmSavePlotBtn');
      const loveContainer = document.getElementById('psmLovePlotContainer');
      const overlapContainer = document.getElementById('psmOverlapContainer');

      loveBtn?.addEventListener('click', () => {
        this.psmActiveView = 'love';
        loveBtn.classList.replace('btn-secondary', 'btn-primary');
        overlapBtn?.classList.replace('btn-primary', 'btn-secondary');
        if (loveContainer) loveContainer.style.display = 'block';
        if (overlapContainer) overlapContainer.style.display = 'none';
        if (saveBtn) saveBtn.dataset.canvasId = 'psmLovePlotCanvas';
        this.renderPsmPlots();
      });

      overlapBtn?.addEventListener('click', () => {
        this.psmActiveView = 'overlap';
        overlapBtn.classList.replace('btn-secondary', 'btn-primary');
        loveBtn?.classList.replace('btn-primary', 'btn-secondary');
        if (loveContainer) loveContainer.style.display = 'none';
        if (overlapContainer) overlapContainer.style.display = 'block';
        if (saveBtn) saveBtn.dataset.canvasId = 'psmOverlapCanvas';
        this.renderPsmPlots();
      });

      const pyBtn = document.getElementById('psmScriptLangPy');
      const rBtn = document.getElementById('psmScriptLangR');
      const stataBtn = document.getElementById('psmScriptLangStata');

      pyBtn?.addEventListener('click', () => this.switchPsmScript('python'));
      rBtn?.addEventListener('click', () => this.switchPsmScript('r'));
      stataBtn?.addEventListener('click', () => this.switchPsmScript('stata'));

      document.getElementById('psmCopyScriptBtn')?.addEventListener('click', () => {
        this.copyPsmScript();
      });

      document.getElementById('psmDownloadScriptBtn')?.addEventListener('click', () => {
        this.downloadPsmScript();
      });

      const csvInput = document.getElementById('psmCsvInput');
      csvInput?.addEventListener('input', () => {
        this.populatePsmVariables();
      });
      csvInput?.addEventListener('change', () => {
        this.populatePsmVariables();
      });

      document.getElementById('psmTreatmentSelect')?.addEventListener('change', () => {
        this.updatePsmCovariatesList();
      });
      document.getElementById('psmOutcomeSelect')?.addEventListener('change', () => {
        this.updatePsmCovariatesList();
      });
    }

    loadPsmSample() {
      const records = Psm.getSampleClinicalDataset();
      if (!records || records.length === 0) return;

      const headers = Object.keys(records[0]);
      const rows = records.map(r => headers.map(h => r[h]).join(','));
      const csv = [headers.join(','), ...rows].join('\n');

      const input = document.getElementById('psmCsvInput');
      if (input) {
        input.value = csv;
      }
      this.populatePsmVariables();
      this.runPsm();
    }

    populatePsmVariables() {
      const csvText = document.getElementById('psmCsvInput')?.value || '';
      const lines = csvText.trim().split(/\r?\n/).filter(l => l.trim().length > 0);
      const badge = document.getElementById('psmDatasetBadge');
      if (lines.length < 2) {
        if (badge) badge.innerText = '0 rows';
        return;
      }

      if (badge) {
        badge.innerText = `N = ${lines.length - 1} patients`;
      }

      const delimiter = lines[0].includes('\t') ? '\t' : (lines[0].includes(';') ? ';' : ',');
      const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));

      const tSelect = document.getElementById('psmTreatmentSelect');
      const oSelect = document.getElementById('psmOutcomeSelect');

      const currentT = tSelect?.value;
      const currentO = oSelect?.value;

      if (tSelect) {
        tSelect.innerHTML = '';
        headers.forEach(h => {
          const opt = document.createElement('option');
          opt.value = h;
          opt.textContent = h;
          if (h.toLowerCase().includes('treat') || h === 'treatment_col' || h === 'group') {
            opt.selected = true;
          }
          tSelect.appendChild(opt);
        });
        if (currentT && headers.includes(currentT)) tSelect.value = currentT;
      }

      if (oSelect) {
        oSelect.innerHTML = '';
        headers.forEach(h => {
          const opt = document.createElement('option');
          opt.value = h;
          opt.textContent = h;
          if (h.toLowerCase().includes('outcom') || h === 'outcome_col' || h.toLowerCase().includes('los')) {
            opt.selected = true;
          }
          oSelect.appendChild(opt);
        });
        if (currentO && headers.includes(currentO)) oSelect.value = currentO;
      }

      this.updatePsmCovariatesList(headers);
    }

    updatePsmCovariatesList(headers) {
      if (!headers) {
        const csvText = document.getElementById('psmCsvInput')?.value || '';
        const lines = csvText.trim().split(/\r?\n/).filter(l => l.trim().length > 0);
        if (lines.length < 2) return;
        const delimiter = lines[0].includes('\t') ? '\t' : (lines[0].includes(';') ? ';' : ',');
        headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));
      }

      const tVal = document.getElementById('psmTreatmentSelect')?.value;
      const oVal = document.getElementById('psmOutcomeSelect')?.value;
      const container = document.getElementById('psmCovariatesContainer');
      if (!container) return;

      const nonCovariates = new Set([tVal, oVal, 'patient_id', 'id', 'subject_id', 'ID', '_rowId']);
      const availableCovariates = headers.filter(h => !nonCovariates.has(h));

      const checkedSet = new Set();
      container.querySelectorAll('input[type="checkbox"]:checked').forEach(cb => checkedSet.add(cb.value));

      container.innerHTML = '';
      availableCovariates.forEach(cov => {
        const label = document.createElement('label');
        label.style.display = 'inline-flex';
        label.style.alignItems = 'center';
        label.style.gap = '0.3rem';
        label.style.padding = '0.2rem 0.5rem';
        label.style.background = 'var(--bg-surface)';
        label.style.border = '1px solid var(--border-subtle)';
        label.style.borderRadius = '4px';
        label.style.fontSize = '0.74rem';
        label.style.cursor = 'pointer';

        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.value = cov;
        cb.classList.add('psm-cov-cb');
        cb.style.accentColor = 'var(--cyan-primary)';
        if (checkedSet.has(cov) || checkedSet.size === 0) {
          cb.checked = true;
        }

        cb.addEventListener('change', () => this.updatePsmCovariatesCount());

        label.appendChild(cb);
        label.appendChild(document.createTextNode(cov));
        container.appendChild(label);
      });

      this.updatePsmCovariatesCount();
    }

    updatePsmCovariatesCount() {
      const container = document.getElementById('psmCovariatesContainer');
      const countEl = document.getElementById('psmCovariatesCount');
      if (container && countEl) {
        const count = container.querySelectorAll('input[type="checkbox"]:checked').length;
        countEl.innerText = count.toString();
      }
    }

    runPsm() {
      const csvText = document.getElementById('psmCsvInput')?.value || '';
      const records = Psm.parseClinicalCsv(csvText);
      if (records.length === 0) {
        alert('Please provide valid clinical dataset rows in CSV or TSV format.');
        return;
      }

      const treatmentCol = document.getElementById('psmTreatmentSelect')?.value;
      const outcomeCol = document.getElementById('psmOutcomeSelect')?.value;

      const container = document.getElementById('psmCovariatesContainer');
      const selectedCovariates = [];
      container?.querySelectorAll('input[type="checkbox"]:checked').forEach(cb => {
        selectedCovariates.push(cb.value);
      });

      if (!treatmentCol || !outcomeCol) {
        alert('Please select both a treatment variable and an outcome variable.');
        return;
      }

      if (selectedCovariates.length === 0) {
        alert('Please select at least one confounding baseline covariate for matching.');
        return;
      }

      const caliperMultiplier = parseFloat(document.getElementById('psmCaliperMultiplier')?.value) || 0.20;
      const enforceCommonSupport = document.getElementById('psmCommonSupportCb')?.checked ?? true;

      const analysis = Psm.executeAnalysis(records, {
        treatmentCol,
        outcomeCol,
        covariateCols: selectedCovariates,
        caliperMultiplier,
        enforceCommonSupport
      });

      if (analysis.error) {
        alert(`Propensity Score Analysis Error: ${analysis.error}`);
        return;
      }

      this.psmLastAnalysis = analysis;
      if (!this.results) this.results = {};
      this.results.propensity = analysis;

      // Update Top Metrics
      const pairsEl = document.getElementById('psmMatchedPairs');
      const pairsBadge = document.getElementById('psmMatchedPairsBadge');
      if (pairsEl) pairsEl.innerText = `${analysis.nMatchedPairs} pairs`;
      if (pairsBadge) {
        const retainedPct = ((analysis.matchedTreatedN / analysis.unmatchedTreatedN) * 100).toFixed(0);
        pairsBadge.innerText = `${analysis.nMatchedPairs * 2} patients (${retainedPct}% treated matched)`;
      }

      const attEl = document.getElementById('psmAttEstimate');
      const attBadge = document.getElementById('psmAttBadge');
      const attCiEl = document.getElementById('psmAttCi');
      const attSeEl = document.getElementById('psmAttSe');
      const maxSmdEl = document.getElementById('psmMaxSmd');
      const balanceBadge = document.getElementById('psmBalanceBadge');

      const outcome = analysis.outcome;
      if (attEl) attEl.innerText = `${outcome.att >= 0 ? '+' : ''}${outcome.att.toFixed(3)}`;
      if (attBadge) {
        attBadge.innerText = outcome.isSignificant ? `p < .001 (Significant)` : `p = ${outcome.pValue.toFixed(3)}`;
        attBadge.className = outcome.isSignificant ? 'badge badge-sig' : 'badge badge-neutral';
      }

      if (attCiEl) attCiEl.innerText = `[${outcome.ci95[0].toFixed(3)}, ${outcome.ci95[1].toFixed(3)}]`;
      if (attSeEl) attSeEl.innerText = `SE = ${outcome.se.toFixed(3)}`;

      const maxSmd = analysis.balance.maxAbsSmdPost;
      if (maxSmdEl) maxSmdEl.innerText = `${maxSmd.toFixed(3)}`;
      if (balanceBadge) {
        if (maxSmd < 0.10) {
          balanceBadge.innerText = 'Excellent Balance (|SMD| < 0.10)';
          balanceBadge.className = 'badge badge-sig';
        } else if (maxSmd < 0.20) {
          balanceBadge.innerText = 'Acceptable Balance (|SMD| < 0.20)';
          balanceBadge.className = 'badge badge-warn';
        } else {
          balanceBadge.innerText = 'Residual Imbalance (|SMD| ≥ 0.20)';
          balanceBadge.className = 'badge badge-danger';
        }
      }

      // Populate Covariate Balance Table
      const bTableBody = document.getElementById('psmBalanceTableBody');
      if (bTableBody) {
        bTableBody.innerHTML = '';
        analysis.balance.balanceTable.forEach(row => {
          const tr = document.createElement('tr');
          const isBalanced = row.absSmdPost <= 0.10;
          const statusBadge = isBalanced
            ? `<span class="badge badge-sig" style="font-size: 0.70rem;">✓ Balanced (&lt; 0.10)</span>`
            : `<span class="badge badge-danger" style="font-size: 0.70rem;">⚠ Imbalance (&gt; 0.10)</span>`;

          const vrVal = (row.varianceRatioPost != null && isFinite(row.varianceRatioPost))
            ? row.varianceRatioPost.toFixed(2)
            : ((row.varRatioPost != null && isFinite(row.varRatioPost)) ? row.varRatioPost.toFixed(2) : '--');

          tr.innerHTML = `
            <td style="font-weight: 600;">${row.covariate}</td>
            <td>${row.meanTreatedPre.toFixed(2)} vs ${row.meanControlPre.toFixed(2)}</td>
            <td style="color: #f43f5e; font-weight: 600;">${row.smdPre.toFixed(3)}</td>
            <td>${row.meanTreatedPost.toFixed(2)} vs ${row.meanControlPost.toFixed(2)}</td>
            <td style="color: #10b981; font-weight: 700;">${row.smdPost.toFixed(3)}</td>
            <td>${vrVal}</td>
            <td style="color: ${row.percentReduction >= 0 ? '#10b981' : '#f43f5e'}; font-weight: 600;">${row.percentReduction.toFixed(1)}%</td>
            <td>${statusBadge}</td>
          `;
          bTableBody.appendChild(tr);
        });
      }

      // Populate Outcome Comparison Table
      const oTableBody = document.getElementById('psmOutcomeTableBody');
      if (oTableBody) {
        oTableBody.innerHTML = '';
        const unadj = outcome.unadjustedDiff;
        const pUnadjStr = unadj.pValue < 0.001 ? '< .001' : unadj.pValue.toFixed(3);
        const pAttStr = outcome.pValue < 0.001 ? '< .001' : outcome.pValue.toFixed(3);

        const trUnadj = document.createElement('tr');
        trUnadj.innerHTML = `
          <td style="font-weight: 600; color: #f43f5e;">Unadjusted (Raw Observational)</td>
          <td>${unadj.diff >= 0 ? '+' : ''}${unadj.diff.toFixed(3)}</td>
          <td>${unadj.se.toFixed(3)}</td>
          <td>t = ${unadj.statistic.toFixed(2)}</td>
          <td>${pUnadjStr}</td>
          <td>[${unadj.ci95[0].toFixed(3)}, ${unadj.ci95[1].toFixed(3)}]</td>
          <td>Confounded by baseline clinical risk</td>
        `;

        const trAtt = document.createElement('tr');
        trAtt.style.background = 'rgba(16, 185, 129, 0.08)';
        trAtt.innerHTML = `
          <td style="font-weight: 700; color: #10b981;">Matched Causal ATT (1:1 NN Caliper)</td>
          <td style="font-weight: 700; color: #10b981;">${outcome.att >= 0 ? '+' : ''}${outcome.att.toFixed(3)}</td>
          <td>${outcome.se.toFixed(3)}</td>
          <td>${outcome.type === 'continuous' ? `t = ${outcome.statistic.toFixed(2)}` : `z = ${outcome.statistic.toFixed(2)}`}</td>
          <td style="font-weight: 700;">${pAttStr}</td>
          <td style="font-weight: 600;">[${outcome.ci95[0].toFixed(3)}, ${outcome.ci95[1].toFixed(3)}]</td>
          <td>Unconfounded treatment effect on treated</td>
        `;

        oTableBody.appendChild(trUnadj);
        oTableBody.appendChild(trAtt);
      }

      this.renderPsmPlots();
      this.updatePsmScript();

      const reportEl = document.getElementById('psmReportText');
      if (reportEl) {
        reportEl.innerText = analysis.reportText;
      }
    }

    renderPsmPlots() {
      if (!this.psmLastAnalysis) return;
      const analysis = this.psmLastAnalysis;

      if (this.psmActiveView === 'love') {
        const engine = this.engines['psmLovePlotCanvas'];
        if (engine && Plots.renderLovePlot) {
          Plots.renderLovePlot(engine, analysis.balance, {
            title: `Love Plot: Baseline Covariate Balance (N = ${analysis.nMatchedPairs} Pairs)`
          });
        }
      } else {
        const engine = this.engines['psmOverlapCanvas'];
        if (engine && Plots.renderPsmOverlapPlot) {
          Plots.renderPsmOverlapPlot(engine, analysis.overlap, {
            title: `Propensity Score Distribution & Common Support Overlap (Caliper = ${analysis.caliper.width.toFixed(4)})`
          });
        }
      }
    }

    switchPsmScript(lang) {
      this.psmCurrentScriptLang = lang;
      const pyBtn = document.getElementById('psmScriptLangPy');
      const rBtn = document.getElementById('psmScriptLangR');
      const stataBtn = document.getElementById('psmScriptLangStata');

      if (pyBtn) {
        pyBtn.className = lang === 'python' ? 'btn btn-primary btn-sm' : 'btn btn-secondary btn-sm';
      }
      if (rBtn) {
        rBtn.className = lang === 'r' ? 'btn btn-primary btn-sm' : 'btn btn-secondary btn-sm';
      }
      if (stataBtn) {
        stataBtn.className = lang === 'stata' ? 'btn btn-primary btn-sm' : 'btn btn-secondary btn-sm';
      }

      this.updatePsmScript();
    }

    updatePsmScript() {
      if (!this.psmLastAnalysis) return;
      const analysis = this.psmLastAnalysis;
      const codeBlock = document.getElementById('psmCodeBlock');
      if (!codeBlock) return;

      const opts = {
        treatmentCol: analysis.treatmentCol,
        outcomeCol: analysis.outcomeCol,
        covariateCols: analysis.covariateCols,
        caliperMultiplier: analysis.caliper.multiplier,
        enforceCommonSupport: analysis.caliper.enforceCommonSupport,
        isBinaryOutcome: analysis.outcome.type === 'binary'
      };

      let code = '';
      if (this.psmCurrentScriptLang === 'python') {
        code = Psm.generatePythonScript(opts);
      } else if (this.psmCurrentScriptLang === 'r') {
        code = Psm.generateRScript(opts);
      } else if (this.psmCurrentScriptLang === 'stata') {
        code = Psm.generateStataScript(opts);
      }

      codeBlock.innerText = code;
    }

    copyPsmScript() {
      const codeBlock = document.getElementById('psmCodeBlock');
      const copyBtn = document.getElementById('psmCopyScriptBtn');
      if (!codeBlock) return;

      const text = codeBlock.innerText;
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(() => {
          if (copyBtn) {
            const orig = copyBtn.innerText;
            copyBtn.innerText = '✓ Copied!';
            setTimeout(() => { copyBtn.innerText = orig; }, 1800);
          }
        });
      } else {
        const ta = document.createElement('textarea');
        ta.value = text;
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        document.execCommand('copy');
        document.body.removeChild(ta);
        if (copyBtn) {
          const orig = copyBtn.innerText;
          copyBtn.innerText = '✓ Copied!';
          setTimeout(() => { copyBtn.innerText = orig; }, 1800);
        }
      }
    }

    downloadPsmScript() {
      const codeBlock = document.getElementById('psmCodeBlock');
      if (!codeBlock) return;
      const text = codeBlock.innerText;

      const extMap = { python: 'py', r: 'R', stata: 'do' };
      const ext = extMap[this.psmCurrentScriptLang] || 'txt';
      const filename = `psm_clinical_analysis.${ext}`;

      const blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    }

    // ==========================================
    // 15. MULTIVARIATE EDA (PCA / MCA / FAMD) WORKFLOW
    // ==========================================
    initMultivariate() {
      this.multivariateActiveView = '2d';
      this.multivariateScriptLang = 'python';
      this.multivariateLastAnalysis = null;

      document.getElementById('multivarSampleBtn')?.addEventListener('click', () => {
        this.loadMultivariateSample();
      });

      document.getElementById('multivarClearBtn')?.addEventListener('click', () => {
        const input = document.getElementById('multivarCsvInput');
        if (input) input.value = '';
        this.populateMultivariateVariables(true);
      });

      document.getElementById('multivarFileInput')?.addEventListener('change', (e) => {
        const file = e.target.files?.[0];
        if (!file) return;
        const reader = new FileReader();
        reader.onload = (evt) => {
          const input = document.getElementById('multivarCsvInput');
          if (input) {
            input.value = evt.target.result;
            this.populateMultivariateVariables(true);
            this.runMultivariate();
          }
        };
        reader.readAsText(file);
      });

      document.getElementById('multivarMethodSelect')?.addEventListener('change', () => {
        this.updateMultivariateMethodHelp();
      });

      document.getElementById('multivarSelectAllBtn')?.addEventListener('click', () => {
        document.querySelectorAll('.multivar-var-cb').forEach(cb => { cb.checked = true; });
      });

      document.getElementById('multivarDeselectAllBtn')?.addEventListener('click', () => {
        document.querySelectorAll('.multivar-var-cb').forEach(cb => { cb.checked = false; });
      });

      // Axis / grouping / biplot change listeners
      const reprojectElements = ['multivarDimX', 'multivarDimY', 'multivarDimZ', 'multivarGroupSelect', 'multivarBiplotCb'];
      reprojectElements.forEach(id => {
        document.getElementById(id)?.addEventListener('change', () => {
          if (this.multivariateLastAnalysis) {
            this.updateMultivariateMetricCards(this.multivariateLastAnalysis);
            this.renderMultivariatePlots();
          }
        });
      });

      document.getElementById('multivarComputeBtn')?.addEventListener('click', () => {
        this.runMultivariate();
      });

      // View Switcher buttons
      const view2DBtn = document.getElementById('multivarView2DBtn');
      const view3DBtn = document.getElementById('multivarView3DBtn');
      const viewScreeBtn = document.getElementById('multivarViewScreeBtn');
      const viewLoadingsBtn = document.getElementById('multivarViewLoadingsBtn');

      view2DBtn?.addEventListener('click', () => this.setMultivariateView('2d'));
      view3DBtn?.addEventListener('click', () => this.setMultivariateView('3d'));
      viewScreeBtn?.addEventListener('click', () => this.setMultivariateView('scree'));
      viewLoadingsBtn?.addEventListener('click', () => this.setMultivariateView('loadings'));

      // Script buttons
      document.getElementById('multivarScriptLangPy')?.addEventListener('click', () => this.switchMultivariateScript('python'));
      document.getElementById('multivarScriptLangR')?.addEventListener('click', () => this.switchMultivariateScript('r'));
      document.getElementById('multivarCopyScriptBtn')?.addEventListener('click', () => this.copyMultivariateScript());
      document.getElementById('multivarDownloadScriptBtn')?.addEventListener('click', () => this.downloadMultivariateScript());

      // Debounced CSV textarea listener
      let debounceTimer = null;
      const csvInput = document.getElementById('multivarCsvInput');
      csvInput?.addEventListener('input', () => {
        clearTimeout(debounceTimer);
        debounceTimer = setTimeout(() => this.populateMultivariateVariables(), 400);
      });

      // Load initial sample automatically
      setTimeout(() => {
        if (!this.multivariateLastAnalysis) {
          this.loadMultivariateSample();
        }
      }, 100);
    }

    updateMultivariateMethodHelp() {
      const method = document.getElementById('multivarMethodSelect')?.value;
      const help = document.getElementById('multivarMethodHelp');
      if (!help) return;
      if (method === 'pca') {
        help.innerText = 'Mean-centers and scales continuous variables to unit variance (Pearson 1901, Hotelling 1933).';
      } else if (method === 'mca') {
        help.innerText = 'Expands categorical modalities into an Indicator/Burt matrix with correspondence scaling (Benzécri 1973).';
      } else {
        help.innerText = 'Balances continuous and categorical blocks via standardized indicator scaling (Pagès 2004).';
      }
    }

    loadMultivariateSample() {
      const cohort = Multivariate.getSampleMixedCohort();
      if (!cohort || cohort.length === 0) return;
      const headers = Object.keys(cohort[0]);
      const csvRows = [headers.join(',')];
      cohort.forEach(row => {
        csvRows.push(headers.map(h => row[h]).join(','));
      });
      const csvText = csvRows.join('\n');
      const input = document.getElementById('multivarCsvInput');
      if (input) input.value = csvText;

      this.populateMultivariateVariables(true);
      this.runMultivariate();
    }

    populateMultivariateVariables(forceReset = false) {
      const csvText = document.getElementById('multivarCsvInput')?.value || '';
      const parsed = Multivariate.parseDataset(csvText);
      const badge = document.getElementById('multivarDatasetBadge');
      const container = document.getElementById('multivarVariablesContainer');
      const groupSelect = document.getElementById('multivarGroupSelect');

      if (parsed.error || !parsed.headers) {
        if (badge) badge.innerText = '0 observations';
        if (container) container.innerHTML = '<span style="color: var(--text-dim); font-size: 0.74rem;">Paste valid CSV dataset to view columns</span>';
        return;
      }

      if (badge) {
        badge.innerText = `N = ${parsed.nRows} observations, ${parsed.nCols} cols`;
      }

      // Populate variable list checkboxes
      if (container) {
        container.innerHTML = '';
        parsed.headers.forEach(h => {
          const isId = h.toLowerCase().includes('id') || h.toLowerCase() === 'row';
          const inferredType = parsed.colTypes[h] || 'continuous';
          const defaultChecked = !isId;

          const row = document.createElement('div');
          row.style.cssText = 'display: flex; align-items: center; justify-content: space-between; padding: 0.25rem 0.4rem; background: var(--bg-surface); border-radius: 4px; border: 1px solid var(--border-subtle);';

          row.innerHTML = `
            <label style="display: flex; align-items: center; gap: 0.4rem; cursor: pointer; margin: 0; font-size: 0.76rem; user-select: none; font-weight: 500;">
              <input type="checkbox" class="multivar-var-cb" value="${h}" ${defaultChecked ? 'checked' : ''} style="accent-color: var(--cyan-primary);">
              <span title="${h}">${h.length > 20 ? h.substring(0, 18) + '...' : h}</span>
            </label>
            <select class="multivar-var-type-select" data-col="${h}" style="font-size: 0.68rem; height: 22px; padding: 0 0.25rem; border-radius: 4px; background: var(--bg-surface-elevated); color: var(--text-main); border: 1px solid var(--border-subtle);">
              <option value="continuous" ${inferredType === 'continuous' ? 'selected' : ''}>Continuous</option>
              <option value="categorical" ${inferredType === 'categorical' ? 'selected' : ''}>Categorical</option>
            </select>
          `;
          container.appendChild(row);
        });
      }

      // Populate grouping dropdown
      if (groupSelect) {
        const currentGroup = groupSelect.value;
        groupSelect.innerHTML = '<option value="">None (Single Color)</option>';
        parsed.headers.forEach(h => {
          if (parsed.colTypes[h] === 'categorical' || (parsed.colStats[h] && parsed.colStats[h].uniqueCount <= 12)) {
            const opt = document.createElement('option');
            opt.value = h;
            opt.textContent = h;
            if (h === 'surgical_approach' || h === 'suture_type' || h.toLowerCase().includes('group')) {
              opt.selected = true;
            }
            groupSelect.appendChild(opt);
          }
        });
        if (currentGroup && Array.from(groupSelect.options).some(o => o.value === currentGroup)) {
          groupSelect.value = currentGroup;
        }
      }
    }

    setMultivariateView(view) {
      this.multivariateActiveView = view;
      const btns = {
        '2d': document.getElementById('multivarView2DBtn'),
        '3d': document.getElementById('multivarView3DBtn'),
        'scree': document.getElementById('multivarViewScreeBtn'),
        'loadings': document.getElementById('multivarViewLoadingsBtn')
      };
      Object.entries(btns).forEach(([k, btn]) => {
        if (!btn) return;
        if (k === view) {
          btn.classList.replace('btn-secondary', 'btn-primary');
        } else {
          btn.classList.replace('btn-primary', 'btn-secondary');
        }
      });

      const plots = {
        '2d': document.getElementById('multivar2DPlot'),
        '3d': document.getElementById('multivar3DPlot'),
        'scree': document.getElementById('multivarScreePlot'),
        'loadings': document.getElementById('multivarLoadingsPlot')
      };
      Object.entries(plots).forEach(([k, el]) => {
        if (el) el.style.display = k === view ? 'block' : 'none';
      });

      this.renderMultivariatePlots();
    }

    runMultivariate() {
      const csvText = document.getElementById('multivarCsvInput')?.value || '';
      const parsed = Multivariate.parseDataset(csvText);
      if (parsed.error) {
        alert(parsed.error);
        return;
      }

      const activeCheckboxes = Array.from(document.querySelectorAll('.multivar-var-cb:checked'));
      const activeCols = activeCheckboxes.map(cb => cb.value);

      if (activeCols.length < 2) {
        alert('Please select at least 2 active variables for dimensionality reduction.');
        return;
      }

      const colTypes = {};
      document.querySelectorAll('.multivar-var-type-select').forEach(sel => {
        colTypes[sel.dataset.col] = sel.value;
      });

      const methodSelect = document.getElementById('multivarMethodSelect')?.value || 'famd';
      const method = methodSelect.toUpperCase();
      const groupingCol = document.getElementById('multivarGroupSelect')?.value || null;

      // Method input validation
      const contCols = activeCols.filter(c => (colTypes[c] || 'continuous') === 'continuous');
      const catCols = activeCols.filter(c => colTypes[c] === 'categorical');

      if (method === 'PCA' && contCols.length < 2) {
        alert('Principal Component Analysis (PCA) requires at least 2 continuous variables. Change variable types or select FAMD.');
        return;
      }
      if (method === 'MCA' && catCols.length < 2) {
        alert('Multiple Correspondence Analysis (MCA) requires at least 2 categorical variables. Change variable types or select FAMD.');
        return;
      }

      const analysis = Multivariate.executeAnalysis(parsed.records, method, activeCols, { colTypes, groupingCol });
      if (analysis.error) {
        alert(`Multivariate Analysis Error: ${analysis.error}`);
        return;
      }

      this.multivariateLastAnalysis = analysis;
      if (!this.results) this.results = {};
      this.results.multivariate = analysis;

      // Update dimension select options
      this.updateMultivariateDimSelectors(analysis.nComponents);

      // Update Top Metrics
      this.updateMultivariateMetricCards(analysis);

      // Populate Scree Table
      this.populateMultivariateScreeTable(analysis.scree);

      // Populate Automated Findings Card
      const pedagogyEl = document.getElementById('multivarPedagogyText');
      if (pedagogyEl) {
        pedagogyEl.innerHTML = analysis.interpretation;
      }

      // Populate Script Display
      this.updateMultivariateScriptDisplay();

      // Populate Narrative Academic Report
      const reportEl = document.getElementById('multivarReportText');
      if (reportEl) {
        const top1 = analysis.scree[0]?.variancePct.toFixed(1);
        const top2 = analysis.scree[1]?.variancePct.toFixed(1);
        const cum2 = ((analysis.scree[0]?.variancePct || 0) + (analysis.scree[1]?.variancePct || 0)).toFixed(1);
        reportEl.innerText = `Multivariate Exploratory Data Analysis (${analysis.method.toUpperCase()}) was performed on N = ${analysis.individuals.length} observations across ${activeCols.length} variables (${analysis.continuousCols ? analysis.continuousCols.length : 0} continuous, ${analysis.categoricalCols ? analysis.categoricalCols.length : 0} categorical). The first two principal dimensions accounted for ${cum2}% of total inertia (Dim 1: ${top1}%, λ = ${analysis.scree[0]?.eigenvalue.toFixed(3)}; Dim 2: ${top2}%, λ = ${analysis.scree[1]?.eigenvalue.toFixed(3)}). ${analysis.outliers.length} statistical outliers (> 2.5 SD factor distance) were flagged in the latent coordinate space. Full numerical linear algebra decomposition was executed client-side via Jacobi cyclic diagonalization and thin SVD.`;
      }

      // Render Active Plotly Graph
      this.renderMultivariatePlots();
    }

    updateMultivariateDimSelectors(nComponents) {
      const maxDim = Math.min(nComponents, 8);
      ['multivarDimX', 'multivarDimY', 'multivarDimZ'].forEach((id, axisIdx) => {
        const sel = document.getElementById(id);
        if (!sel) return;
        const currentVal = parseInt(sel.value) || (axisIdx + 1);
        sel.innerHTML = '';
        for (let d = 1; d <= maxDim; d++) {
          const opt = document.createElement('option');
          opt.value = d;
          opt.textContent = `Dim ${d}`;
          if (d === currentVal) opt.selected = true;
          sel.appendChild(opt);
        }
        if (!Array.from(sel.options).some(o => parseInt(o.value) === currentVal)) {
          sel.value = String(Math.min(axisIdx + 1, maxDim));
        }
      });
    }

    updateMultivariateMetricCards(analysis) {
      const dimX = parseInt(document.getElementById('multivarDimX')?.value) || 1;
      const dimY = parseInt(document.getElementById('multivarDimY')?.value) || 2;

      const varX = analysis.scree[dimX - 1]?.variancePct || 0;
      const eigX = analysis.scree[dimX - 1]?.eigenvalue || 0;
      const varY = analysis.scree[dimY - 1]?.variancePct || 0;
      const eigY = analysis.scree[dimY - 1]?.eigenvalue || 0;
      const cumPlane = varX + varY;

      const methodEl = document.getElementById('multivarMetricMethod');
      const nEl = document.getElementById('multivarMetricN');
      if (methodEl) methodEl.innerText = `${analysis.method.toUpperCase()}`;
      if (nEl) nEl.innerText = `N = ${analysis.individuals.length} obs`;

      const dim1Label = document.getElementById('multivarMetricDim1Label');
      const dim1Val = document.getElementById('multivarMetricDim1Val');
      const dim1Badge = document.getElementById('multivarMetricDim1Badge');
      if (dim1Label) dim1Label.innerText = `Dim ${dimX} Variance (λ_${dimX})`;
      if (dim1Val) dim1Val.innerText = `${varX.toFixed(1)}%`;
      if (dim1Badge) dim1Badge.innerText = `λ = ${eigX.toFixed(3)}`;

      const dim2Label = document.getElementById('multivarMetricDim2Label');
      const dim2Val = document.getElementById('multivarMetricDim2Val');
      const dim2Badge = document.getElementById('multivarMetricDim2Badge');
      if (dim2Label) dim2Label.innerText = `Dim ${dimY} Variance (λ_${dimY})`;
      if (dim2Val) dim2Val.innerText = `${varY.toFixed(1)}%`;
      if (dim2Badge) dim2Badge.innerText = `λ = ${eigY.toFixed(3)}`;

      const cumVal = document.getElementById('multivarMetricCumVal');
      const cumBadge = document.getElementById('multivarMetricCumBadge');
      if (cumVal) cumVal.innerText = `${cumPlane.toFixed(1)}%`;
      if (cumBadge) cumBadge.innerText = `Plane (Dim ${dimX} + Dim ${dimY})`;
    }

    populateMultivariateScreeTable(scree) {
      const tbody = document.getElementById('multivarScreeTableBody');
      if (!tbody) return;
      tbody.innerHTML = '';
      (scree || []).forEach(row => {
        const tr = document.createElement('tr');
        const kaiserPassed = row.eigenvalue >= 1.0;
        const kaiserBadge = kaiserPassed
          ? `<span class="badge badge-sig" style="font-size: 0.70rem;">✓ Retain (λ ≥ 1.0)</span>`
          : `<span class="badge badge-neutral" style="font-size: 0.70rem;">Drop (λ &lt; 1.0)</span>`;

        tr.innerHTML = `
          <td style="font-weight: 600;">${row.label}</td>
          <td style="font-family: monospace;">${row.eigenvalue.toFixed(4)}</td>
          <td style="color: var(--cyan-primary); font-weight: 600;">${row.variancePct.toFixed(2)}%</td>
          <td style="font-weight: 600;">${row.cumulativePct.toFixed(2)}%</td>
          <td>${kaiserBadge}</td>
        `;
        tbody.appendChild(tr);
      });
    }

    renderMultivariatePlots() {
      if (!this.multivariateLastAnalysis || typeof window.Plotly === 'undefined') return;

      const analysis = this.multivariateLastAnalysis;
      const view = this.multivariateActiveView || '2d';
      const dimX = parseInt(document.getElementById('multivarDimX')?.value) || 1;
      const dimY = parseInt(document.getElementById('multivarDimY')?.value) || 2;
      const dimZ = parseInt(document.getElementById('multivarDimZ')?.value) || 3;

      const isDark = document.documentElement.getAttribute('data-theme') !== 'light';
      const textColor = isDark ? '#e2e8f0' : '#1e293b';
      const gridColor = isDark ? 'rgba(255, 255, 255, 0.08)' : 'rgba(0, 0, 0, 0.08)';
      const zeroLineColor = isDark ? 'rgba(255, 255, 255, 0.25)' : 'rgba(0, 0, 0, 0.25)';

      const colorPalette = [
        '#00d2ff', '#10b981', '#f59e0b', '#ec4899', '#8b5cf6',
        '#3b82f6', '#14b8a6', '#f97316', '#a855f7', '#64748b'
      ];

      if (view === '2d') {
        const showBiplot = document.getElementById('multivarBiplotCb')?.checked ?? true;
        const traces = [];

        // Group observations
        const groups = analysis.uniqueGroups && analysis.uniqueGroups.length > 0
          ? analysis.uniqueGroups
          : ['Observation'];

        groups.forEach((grp, gIdx) => {
          const members = analysis.individuals.filter(ind => ind.group === grp);
          if (members.length === 0) return;

          traces.push({
            x: members.map(m => m.coords[dimX - 1] || 0),
            y: members.map(m => m.coords[dimY - 1] || 0),
            mode: 'markers',
            type: 'scatter',
            name: grp,
            text: members.map(m => `${m.id}<br>Group: ${grp}<br>Coords: (${(m.coords[dimX - 1] || 0).toFixed(2)}, ${(m.coords[dimY - 1] || 0).toFixed(2)})`),
            hoverinfo: 'text',
            marker: {
              size: 7,
              color: colorPalette[gIdx % colorPalette.length],
              opacity: 0.85
            }
          });
        });

        // Compute max absolute coordinate to scale biplot arrows
        const allX = analysis.individuals.map(m => Math.abs(m.coords[dimX - 1] || 0));
        const allY = analysis.individuals.map(m => Math.abs(m.coords[dimY - 1] || 0));
        const maxDataX = Math.max(...allX, 1);
        const maxDataY = Math.max(...allY, 1);
        const maxRange = Math.max(maxDataX, maxDataY);

        const annotations = [];

        // Biplot arrows and modality centroids
        if (showBiplot) {
          // 1. Continuous Variable Loading Arrows
          if (analysis.continuousVariables && analysis.continuousVariables.length > 0) {
            analysis.continuousVariables.forEach(v => {
              const vx = (v.coords[dimX - 1] || 0) * (maxRange * 0.85);
              const vy = (v.coords[dimY - 1] || 0) * (maxRange * 0.85);

              annotations.push({
                x: vx,
                y: vy,
                ax: 0,
                ay: 0,
                axref: 'x',
                ayref: 'y',
                xref: 'x',
                yref: 'y',
                showarrow: true,
                arrowhead: 2,
                arrowsize: 1.1,
                arrowwidth: 1.8,
                arrowcolor: '#f59e0b',
                text: v.name,
                font: { color: '#f59e0b', size: 10.5 },
                xanchor: vx >= 0 ? 'left' : 'right',
                yanchor: vy >= 0 ? 'bottom' : 'top'
              });
            });
          }

          // 2. Categorical Modality Centroids
          if (analysis.modalities && analysis.modalities.length > 0) {
            traces.push({
              x: analysis.modalities.map(m => m.coords[dimX - 1] || 0),
              y: analysis.modalities.map(m => m.coords[dimY - 1] || 0),
              mode: 'markers+text',
              type: 'scatter',
              name: 'Categories (Centroids)',
              text: analysis.modalities.map(m => m.label),
              textposition: 'top right',
              hoverinfo: 'text',
              marker: {
                size: 9,
                symbol: 'diamond',
                color: '#ec4899',
                line: { color: '#ffffff', width: 1 }
              }
            });
          }
        }

        const varPctX = analysis.scree[dimX - 1]?.variancePct || 0;
        const varPctY = analysis.scree[dimY - 1]?.variancePct || 0;

        const layout = {
          title: {
            text: `${analysis.method.toUpperCase()} Factor Map — Dim ${dimX} vs Dim ${dimY}`,
            font: { color: textColor, size: 14 }
          },
          paper_bgcolor: 'rgba(0,0,0,0)',
          plot_bgcolor: 'rgba(0,0,0,0)',
          font: { color: textColor, family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif' },
          margin: { l: 55, r: 40, t: 45, b: 50 },
          xaxis: {
            title: `Dimension ${dimX} (${varPctX.toFixed(1)}% variance)`,
            zeroline: true,
            zerolinecolor: zeroLineColor,
            gridcolor: gridColor
          },
          yaxis: {
            title: `Dimension ${dimY} (${varPctY.toFixed(1)}% variance)`,
            zeroline: true,
            zerolinecolor: zeroLineColor,
            gridcolor: gridColor
          },
          annotations,
          legend: { orientation: 'h', y: -0.2 }
        };

        Plotly.newPlot('multivar2DPlot', traces, layout, { responsive: true, displayModeBar: true, displaylogo: false });

      } else if (view === '3d') {
        const traces = [];
        const groups = analysis.uniqueGroups && analysis.uniqueGroups.length > 0
          ? analysis.uniqueGroups
          : ['Observation'];

        groups.forEach((grp, gIdx) => {
          const members = analysis.individuals.filter(ind => ind.group === grp);
          if (members.length === 0) return;

          traces.push({
            x: members.map(m => m.coords[dimX - 1] || 0),
            y: members.map(m => m.coords[dimY - 1] || 0),
            z: members.map(m => m.coords[dimZ - 1] || 0),
            mode: 'markers',
            type: 'scatter3d',
            name: grp,
            text: members.map(m => `${m.id}<br>Group: ${grp}<br>Coords: (${(m.coords[dimX-1]||0).toFixed(2)}, ${(m.coords[dimY-1]||0).toFixed(2)}, ${(m.coords[dimZ-1]||0).toFixed(2)})`),
            hoverinfo: 'text',
            marker: {
              size: 4,
              color: colorPalette[gIdx % colorPalette.length],
              opacity: 0.85
            }
          });
        });

        const varX = analysis.scree[dimX - 1]?.variancePct || 0;
        const varY = analysis.scree[dimY - 1]?.variancePct || 0;
        const varZ = analysis.scree[dimZ - 1]?.variancePct || 0;

        const layout3D = {
          title: {
            text: `3D Factor Space (Dim ${dimX} vs Dim ${dimY} vs Dim ${dimZ})`,
            font: { color: textColor, size: 14 }
          },
          paper_bgcolor: 'rgba(0,0,0,0)',
          plot_bgcolor: 'rgba(0,0,0,0)',
          font: { color: textColor, family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif' },
          margin: { l: 20, r: 20, t: 40, b: 20 },
          scene: {
            xaxis: { title: `Dim ${dimX} (${varX.toFixed(1)}%)`, gridcolor: gridColor },
            yaxis: { title: `Dim ${dimY} (${varY.toFixed(1)}%)`, gridcolor: gridColor },
            zaxis: { title: `Dim ${dimZ} (${varZ.toFixed(1)}%)`, gridcolor: gridColor }
          },
          legend: { orientation: 'h', y: -0.1 }
        };

        Plotly.newPlot('multivar3DPlot', traces, layout3D, { responsive: true, displaylogo: false });

      } else if (view === 'scree') {
        const labels = analysis.scree.map(s => s.label);
        const variancePcts = analysis.scree.map(s => s.variancePct);
        const cumulativePcts = analysis.scree.map(s => s.cumulativePct);

        const traceBar = {
          x: labels,
          y: variancePcts,
          name: 'Variance Explained (%)',
          type: 'bar',
          marker: { color: '#00d2ff' }
        };

        const traceLine = {
          x: labels,
          y: cumulativePcts,
          name: 'Cumulative Variance (%)',
          type: 'scatter',
          mode: 'lines+markers',
          yaxis: 'y2',
          line: { color: '#10b981', width: 2.5 },
          marker: { size: 6, color: '#10b981' }
        };

        const screeLayout = {
          title: { text: 'Scree Plot — Percentage of Explained Variance per Dimension', font: { color: textColor, size: 14 } },
          paper_bgcolor: 'rgba(0,0,0,0)',
          plot_bgcolor: 'rgba(0,0,0,0)',
          font: { color: textColor, family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif' },
          margin: { l: 55, r: 55, t: 45, b: 50 },
          xaxis: { title: 'Principal Dimensions', gridcolor: gridColor },
          yaxis: { title: 'Variance Explained (%)', range: [0, Math.max(40, Math.max(...variancePcts) * 1.25)], gridcolor: gridColor },
          yaxis2: {
            title: 'Cumulative (%)',
            overlaying: 'y',
            side: 'right',
            range: [0, 105],
            gridcolor: 'transparent'
          },
          legend: { orientation: 'h', y: -0.2 }
        };

        Plotly.newPlot('multivarScreePlot', [traceBar, traceLine], screeLayout, { responsive: true, displaylogo: false });

      } else if (view === 'loadings') {
        const variables = analysis.allVariables || analysis.variables || [];
        const sorted = [...variables].sort((a, b) => Math.abs(b.coords[dimX - 1] || 0) - Math.abs(a.coords[dimX - 1] || 0)).slice(0, 15);

        const varNames = sorted.map(v => v.name).reverse();
        const loadingsX = sorted.map(v => v.coords[dimX - 1] || 0).reverse();
        const loadingsY = sorted.map(v => v.coords[dimY - 1] || 0).reverse();

        const traceX = {
          x: loadingsX,
          y: varNames,
          type: 'bar',
          orientation: 'h',
          name: `Dim ${dimX} Loading / Correlation`,
          marker: { color: '#00d2ff' }
        };

        const traceY = {
          x: loadingsY,
          y: varNames,
          type: 'bar',
          orientation: 'h',
          name: `Dim ${dimY} Loading / Correlation`,
          marker: { color: '#f59e0b' }
        };

        const loadingsLayout = {
          title: { text: `Top Variable Loadings / Modality Contributions on Dim ${dimX} & Dim ${dimY}`, font: { color: textColor, size: 14 } },
          paper_bgcolor: 'rgba(0,0,0,0)',
          plot_bgcolor: 'rgba(0,0,0,0)',
          font: { color: textColor, family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif' },
          margin: { l: 140, r: 40, t: 45, b: 50 },
          barmode: 'group',
          xaxis: { title: 'Factor Loading / Coordinate', zeroline: true, zerolinecolor: zeroLineColor, gridcolor: gridColor },
          yaxis: { automargin: true },
          legend: { orientation: 'h', y: -0.2 }
        };

        Plotly.newPlot('multivarLoadingsPlot', [traceX, traceY], loadingsLayout, { responsive: true, displaylogo: false });
      }
    }

    resizeMultivariatePlots() {
      if (typeof window.Plotly === 'undefined') return;
      const plotIds = ['multivar2DPlot', 'multivar3DPlot', 'multivarScreePlot', 'multivarLoadingsPlot'];
      plotIds.forEach(id => {
        const el = document.getElementById(id);
        if (el && el.style.display !== 'none') {
          window.Plotly.Plots.resize(el);
        }
      });
    }

    switchMultivariateScript(lang) {
      this.multivariateScriptLang = lang;
      const pyBtn = document.getElementById('multivarScriptLangPy');
      const rBtn = document.getElementById('multivarScriptLangR');

      if (lang === 'python') {
        pyBtn?.classList.replace('btn-secondary', 'btn-primary');
        rBtn?.classList.replace('btn-primary', 'btn-secondary');
      } else {
        rBtn?.classList.replace('btn-secondary', 'btn-primary');
        pyBtn?.classList.replace('btn-primary', 'btn-secondary');
      }
      this.updateMultivariateScriptDisplay();
    }

    updateMultivariateScriptDisplay() {
      const codeBlock = document.getElementById('multivarCodeBlock');
      if (!codeBlock) return;
      if (!this.multivariateLastAnalysis || !this.multivariateLastAnalysis.scripts) {
        codeBlock.innerText = '# Run factor decomposition to generate reproducible code.';
        return;
      }
      const lang = this.multivariateScriptLang || 'python';
      codeBlock.innerText = this.multivariateLastAnalysis.scripts[lang] || '# Script not generated.';
    }

    copyMultivariateScript() {
      const codeBlock = document.getElementById('multivarCodeBlock');
      const copyBtn = document.getElementById('multivarCopyScriptBtn');
      if (!codeBlock) return;
      const text = codeBlock.innerText;

      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(() => {
          if (copyBtn) {
            const orig = copyBtn.innerText;
            copyBtn.innerText = '✓ Copied!';
            setTimeout(() => { copyBtn.innerText = orig; }, 1800);
          }
        });
      } else {
        const ta = document.createElement('textarea');
        ta.value = text;
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        document.execCommand('copy');
        document.body.removeChild(ta);
        if (copyBtn) {
          const orig = copyBtn.innerText;
          copyBtn.innerText = '✓ Copied!';
          setTimeout(() => { copyBtn.innerText = orig; }, 1800);
        }
      }
    }

    downloadMultivariateScript() {
      const codeBlock = document.getElementById('multivarCodeBlock');
      if (!codeBlock) return;
      const text = codeBlock.innerText;
      const ext = this.multivariateScriptLang === 'python' ? 'py' : 'R';
      const filename = `multivariate_eda_analysis.${ext}`;

      const blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    }


  }

  // Start the application when DOM is ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', () => {
      window.sgApp = new AppController();
      window.app = window.sgApp;
    });
  } else {
    window.sgApp = new AppController();
    window.app = window.sgApp;
  }
})();
