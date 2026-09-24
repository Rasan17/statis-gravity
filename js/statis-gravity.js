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
        testName: 'Mann-Whitney U Test',
        n1, n2,
        rankSumA, rankSumB,
        u1, u2,
        statistic: U,
        zScore: z,
        pValue,
        rankBiserial,
        isSignificant: pValue < 0.05
      };
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

      const processed = groups.map(g => {
        const stats = Descriptive.calculate(g.data);
        return { name: g.name || 'Group', stats, data: stats.values };
      }).filter(g => g.stats.n > 0);

      const k = processed.length;
      if (k < 2) return { error: 'At least 2 groups must have valid data.' };

      const totalN = processed.reduce((acc, g) => acc + g.stats.n, 0);
      const grandSum = processed.reduce((acc, g) => acc + g.stats.sum, 0);
      const grandMean = grandSum / totalN;

      const dfBetween = k - 1;
      const dfWithin = totalN - k;
      if (dfWithin <= 0) return { error: 'Insufficient degrees of freedom.' };

      let ssBetween = 0;
      let ssWithin = 0;
      for (const g of processed) {
        ssBetween += g.stats.n * Math.pow(g.stats.mean - grandMean, 2);
        for (const val of g.data) ssWithin += Math.pow(val - g.stats.mean, 2);
      }
      const ssTotal = ssBetween + ssWithin;

      const msBetween = ssBetween / dfBetween;
      const msWithin = ssWithin / dfWithin;
      const F = msWithin === 0 ? 0 : msBetween / msWithin;
      const pValue = Distributions.fPValue(F, dfBetween, dfWithin);

      const etaSquared = ssTotal === 0 ? 0 : ssBetween / ssTotal;
      const omegaSquared = (ssTotal + msWithin) === 0 ? 0 :
        (ssBetween - dfBetween * msWithin) / (ssTotal + msWithin);

      // Tukey's HSD Post-Hoc Pairwise Contrasts
      const pairwise = [];
      const numPairs = (k * (k - 1)) / 2;
      const pooledSD = Math.sqrt(msWithin);

      for (let i = 0; i < k; i++) {
        for (let j = i + 1; j < k; j++) {
          const gA = processed[i];
          const gB = processed[j];
          const meanDiff = gA.stats.mean - gB.stats.mean;
          const seDiff = Math.sqrt(msWithin * (1 / gA.stats.n + 1 / gB.stats.n));
          const seTukey = Math.sqrt((msWithin / 2) * (1 / gA.stats.n + 1 / gB.stats.n));
          const q = seTukey === 0 ? 0 : Math.abs(meanDiff) / seTukey;
          const tEquiv = q / Math.SQRT2;
          const pRaw = Distributions.tPValue(tEquiv, dfWithin);
          const pAdjusted = Math.min(1.0, pRaw * numPairs);

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
            tStatistic: tEquiv,
            pValue: pAdjusted,
            pValueRaw: pRaw,
            ci95,
            cohensD,
            isSignificant: pAdjusted < 0.05
          });
        }
      }

      return {
        testName: 'One-Way ANOVA',
        k, totalN, grandMean,
        groups: processed,
        dfBetween, dfWithin,
        ssBetween, ssWithin, ssTotal,
        msBetween, msWithin,
        fStatistic: F,
        pValue,
        etaSquared,
        omegaSquared: Math.max(0, omegaSquared),
        pairwise,
        isSignificant: pValue < 0.05
      };
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

    d.addHeading2('Comparative Inferential Test Results');
    d.addTable(
      ['Inferential Parameter', 'Calculated Value', 'Clinical Interpretation / Benchmark'],
      [
        ['Test Statistic', `${data.testName.includes('Mann-Whitney') ? 'U = ' : 't = '}${data.statistic.toFixed(3)}`, 'Standardized difference between cohort locations'],
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
      .addSubTitle('Module: Multi-Cohort Variance & Tukey HSD Post-Hoc Pairwise Analysis')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Cohorts Information & Input Data')
      .addParagraph(`Number of Independent Cohorts (k): ${data.k || data.groups.length}`)
      .addParagraph(`Total Analyzed Sample Size (N): ${data.totalN} patients/specimens`)
      .addParagraph(`Grand Mean across All Cohorts: ${data.grandMean.toFixed(2)}`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
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

    d.addHeading2('One-Way ANOVA Summary Table');
    d.addTable(
      ['Source of Variation', 'Sum of Squares (SS)', 'Degrees of Freedom (df)', 'Mean Square (MS)', 'F-Statistic', 'p-Value', 'Omega-Squared (ω²)'],
      [
        ['Between Groups (Treatment)', `${(data.ssBetween || 0).toFixed(2)}`, `${data.dfBetween || 0}`, `${(data.msBetween || 0).toFixed(2)}`, `F = ${(data.fStatistic || 0).toFixed(2)}`, `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + (data.pValue || 0).toFixed(4)}`, `${(data.omegaSquared || 0).toFixed(3)}`],
        ['Within Groups (Residual/Error)', `${(data.ssWithin || 0).toFixed(2)}`, `${data.dfWithin || 0}`, `${(data.msWithin || 0).toFixed(2)}`, '-', '-', `Eta² (η²) = ${(data.etaSquared || 0).toFixed(3)}`],
        ['Total', `${(data.ssTotal !== undefined ? data.ssTotal : ((data.ssBetween || 0) + (data.ssWithin || 0))).toFixed(2)}`, `${(data.dfBetween || 0) + (data.dfWithin || 0)}`, '-', '-', '-', '-']
      ]
    );

    if (data.pairwise && data.pairwise.length > 0) {
      d.addHeading2('Tukey\'s HSD Post-Hoc Pairwise Contrasts');
      const pairRows = data.pairwise.map(p => [
        p.comparison,
        `${(p.meanDiff >= 0 ? '+' : '')}${p.meanDiff.toFixed(2)}`,
        `${p.seDiff.toFixed(3)}`,
        `q = ${p.qStatistic.toFixed(2)}`,
        `${p.pValue < 0.001 ? 'p < .001' : 'p = ' + p.pValue.toFixed(4)}`,
        `[${p.ci95[0].toFixed(2)}, ${p.ci95[1].toFixed(2)}]`,
        `d = ${p.cohensD.toFixed(2)}`,
        p.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)'
      ]);
      d.addTable(['Pairwise Contrast', 'Mean Diff (ΔM)', 'Std Error', 'Tukey q', 'Adjusted p', '95% CI of Diff', 'Cohen\'s d', 'Significance'], pairRows);
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'ANOVA & Post-Hoc APA Clinical Summary',
      data.reportText || 'ANOVA narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('One-Way Omnibus ANOVA: Selected because testing multiple cohorts with uncorrected pairwise t-tests results in severe Family-Wise Error Rate inflation (FWER). With 3 cohorts, 3 comparisons yield α_FW = 1 - (1 - 0.05)³ = 14.3%; with 5 cohorts (10 comparisons), α_FW exceeds 40%. ANOVA provides a rigorous omnibus test that simultaneously assesses whether any between-cohort variance exceeds within-cohort residual variation.');
    d.addBullet('Tukey\'s Honest Significant Difference (HSD): Chosen as the post-hoc method because it utilizes the Studentized Range distribution (q) to strictly bound the overall Family-Wise Error Rate at α = 0.05 across all possible pairwise comparisons, while preserving substantially greater statistical power than overly conservative Bonferroni adjustments.');
    d.addBullet('Omega-Squared (ω²) Reporting: Included alongside Eta-squared (η²) because Eta-squared represents a sample proportion of variance that is positively biased in small clinical samples. Omega-squared provides an unbiased population effect size estimate.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('ANOVA partitions the total sum of squares into treatment (between) and error (within) components: SS_Total = SS_Between + SS_Within.');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Between-Groups Mean Square: MS_B = SS_B / (k - 1).');
    d.addBullet('Within-Groups Mean Square: MS_W = SS_W / (N - k).');
    d.addBullet('F-Ratio: F = MS_B / MS_W ~ F(k-1, N-k).');
    d.addBullet('Tukey Studentized Range: q = |x̄A - x̄B| / √[ (MS_W / 2) (1/nA + 1/nB) ].');
    d.addParagraph('Key Academic References:');
    d.addBullet('Fisher RA (1925). Statistical Methods for Research Workers. Oliver and Boyd, Edinburgh.');
    d.addBullet('Tukey JW (1949). Comparing individual means in the analysis of variance. Biometrics, 5(2): 99–114.');
    d.addBullet('Hayter AJ (1984). A proof of the conjecture that the Tukey-Kramer multiple comparisons procedure is conservative. Annals of Statistics, 12(1): 61–75.');

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

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Pedagogical Synthesis & Clinical Trial Relevance',
      data.reportText || 'The Central Limit Theorem and Student\'s t convergence demonstrate the mathematical foundations of parametric testing in clinical trials.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Foundation of Inferential Biostatistics: Parametric hypothesis tests (Student t-test, ANOVA, ordinary least squares regression) mathematically assume normally distributed errors or sample means. The Central Limit Theorem provides the mathematical justification for deploying these tests in clinical trials with n ≥ 30 even when raw clinical metrics (e.g. ICU stay, recovery hours) are skewed.');
    d.addBullet('Gosset\'s Student\'s t Adjustment: In small clinical cohorts (n < 30), estimating population variance σ² using sample variance s² introduces substantial stochastic instability into the test statistic denominator. Using Gaussian critical values (z = 1.96) severely inflates the Type I error rate (e.g. to 14.5% at n = 4). Student\'s t distribution compensates for this extra uncertainty by thickening the tails and demanding a higher critical threshold (t = 3.182 at n = 4).');
    d.addBullet('The n ≥ 31 Clinical Threshold: As demonstrated by the simulation, when sample size reaches n ≥ 31 (degrees of freedom ν ≥ 30), the critical t cutoff drops to 2.042 (only 4.2% wider than 1.960), and tail probability converges close to 5.0%. This mathematical threshold explains why sample sizes of 30 or greater historically permit Gaussian approximation in medical trial protocols.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Classical Lindberg-Lévy Central Limit Theorem: Let X₁, X₂, ..., X_n be independent and identically distributed (i.i.d.) random variables with mean μ and finite variance σ². Then as n → ∞: √n (X̄_n - μ) / σ → N(0, 1).');
    d.addBullet('Student\'s t Distribution Density: f(t; ν) = [ Γ((ν+1)/2) / (√(πν) Γ(ν/2)) ] · [ 1 + t²/ν ]^{-(ν+1)/2}. As ν → ∞, [ 1 + t²/ν ]^{-(ν+1)/2} → exp(-t²/2), converging to Standard Normal N(0, 1).');
    d.addBullet('Standard Error of the Mean: SE = σ / √n. Quadrupling patient enrollment cuts the estimation uncertainty in half.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Student [Gosset WS] (1908). The probable error of a mean. Biometrika, 6(1): 1–25.');
    d.addBullet('Laplace PS (1810). Mémoire sur les approximations des formules qui sont fonctions de très grands nombres et sur leur application aux probabilités. Mémoires de l\'Académie Royale des Sciences de Paris.');
    d.addBullet('Gauss CF (1809). Theoria motus corporum coelestium in sectionibus conicis solem ambientium. Hamburg: Perthes et Besser.');
    d.addBullet('Altman DG, Bland JM (1995). Statistics Notes: The normal distribution. BMJ, 310(6975): 298–299.');

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
  // 11. CANVAS CHART ENGINE & PLOTS
  // ==========================================
  class ChartEngine {
    constructor(canvasId, options = {}) {
      this.canvas = typeof canvasId === 'string' ? document.getElementById(canvasId) : canvasId;
      if (!this.canvas) return;
      this.ctx = this.canvas.getContext('2d');
      this.options = Object.assign({
        padding: { top: 35, right: 25, bottom: 45, left: 55 },
        theme: 'dark'
      }, options);

      this.colors = {
        dark: {
          bg: '#090d16',
          grid: 'rgba(255, 255, 255, 0.07)',
          axis: 'rgba(255, 255, 255, 0.22)',
          text: '#94a3b8',
          textBold: '#f8fafc',
          primary: '#00d2ff',
          secondary: '#a855f7',
          accent: '#10b981',
          danger: '#f43f5e'
        },
        light: {
          bg: '#ffffff',
          grid: 'rgba(0, 0, 0, 0.06)',
          axis: 'rgba(0, 0, 0, 0.2)',
          text: '#64748b',
          textBold: '#0f172a',
          primary: '#0284c7',
          secondary: '#9333ea',
          accent: '#059669',
          danger: '#e11d48'
        }
      };

      this.initHiDPI();
    }

    get palette() {
      return this.colors[this.options.theme] || this.colors.dark;
    }

    initHiDPI() {
      if (!this.canvas) return;
      const rect = this.canvas.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      this.width = rect.width || 600;
      this.height = rect.height || 320;

      this.canvas.width = this.width * dpr;
      this.canvas.height = this.height * dpr;
      if (this.ctx.resetTransform) this.ctx.resetTransform();
      this.ctx.scale(dpr, dpr);
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
      this.init();
    }

    init() {
      this.applyTheme(this.theme);
      this.initTabs();
      this.initEngines();
      this.bindEvents();
      this.loadInitialData();
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
              if (canvas && this.engines[canvas.id]?.lastRender) {
                this.engines[canvas.id].initHiDPI();
                this.engines[canvas.id].lastRender();
              }
            });
          }
        });
      });
    }

    initEngines() {
      ['descCanvas', 'descBoxCanvas', 'descViolinCanvas', 'hypoCanvas', 'anovaCanvas', 'corrCanvas', 'rocCanvas', 'teachingDistCanvas', 'teachingCltParentCanvas', 'teachingCltSamplingCanvas', 'teachingTCanvas'].forEach(id => {
        const el = document.getElementById(id);
        if (el) {
          this.engines[id] = new ChartEngine(el, { theme: this.theme });
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
      document.getElementById('hypoSampleBtn')?.addEventListener('click', () => {
        document.getElementById('hypoGroupA').value = DataParser.samples.icpDynamics.groupA.join(', ');
        document.getElementById('hypoGroupB').value = DataParser.samples.icpDynamics.groupB.join(', ');
        this.runHypo();
      });
      document.getElementById('hypoErrorBarMode')?.addEventListener('change', () => {
        this.runHypo();
      });

      // 3. ANOVA
      document.getElementById('anovaComputeBtn')?.addEventListener('click', () => this.runAnova());
      document.getElementById('anovaSampleBtn')?.addEventListener('click', () => {
        const ca = DataParser.samples.cranialAsymmetry;
        document.getElementById('anovaG1').value = ca[0].join(', ');
        document.getElementById('anovaG2').value = ca[1].join(', ');
        document.getElementById('anovaG3').value = ca[2].join(', ');
        this.runAnova();
      });
      document.getElementById('anovaErrorBarMode')?.addEventListener('change', () => {
        this.runAnova();
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
      const test = document.getElementById('hypoTestType')?.value || 'welch';

      const a = DataParser.parseSeries(rawA);
      const b = DataParser.parseSeries(rawB);

      let res;
      if (test === 'student') res = Hypothesis.independentTTest(a, b);
      else if (test === 'welch') res = Hypothesis.welchTTest(a, b);
      else if (test === 'paired') res = Hypothesis.pairedTTest(a, b);
      else res = Hypothesis.mannWhitneyUTest(a, b);

      if (res.error) {
        alert(res.error);
        return;
      }

      res.groupA = Object.assign(res.groupA || Descriptive.calculate(a), { name: nameA });
      res.groupB = Object.assign(res.groupB || Descriptive.calculate(b), { name: nameB });

      document.getElementById('hypoStat').innerText = (res.statistic || res.zScore || 0).toFixed(2);
      document.getElementById('hypoPVal').innerText = Exporter.formatP(res.pValue);
      const badge = document.getElementById('hypoPValBadge');
      badge.className = `badge ${res.isSignificant ? 'badge-sig' : 'badge-ns'}`;
      badge.innerText = res.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)';
      document.getElementById('hypoEffect').innerText = (res.cohensD !== undefined ? res.cohensD.toFixed(2) : (res.rankBiserial || 0).toFixed(2));
      document.getElementById('hypoDiff').innerText = res.meanDiff !== undefined ? res.meanDiff.toFixed(2) : 'N/A';

      const report = `A ${res.testName} demonstrated a ${res.isSignificant ? 'statistically significant' : 'non-significant'} difference between ${nameA} (M = ${res.groupA.mean.toFixed(2)}, SD = ${res.groupA.sd.toFixed(2)}) and ${nameB} (M = ${res.groupB.mean.toFixed(2)}, SD = ${res.groupB.sd.toFixed(2)}), ${Exporter.formatP(res.pValue)}, Cohen's d = ${(res.cohensD || 0).toFixed(2)}.`;
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
      this.results = this.results || {}; this.results.hypothesis = res;
    }

    runAnova() {
      const g1 = DataParser.parseSeries(document.getElementById('anovaG1')?.value || '');
      const g2 = DataParser.parseSeries(document.getElementById('anovaG2')?.value || '');
      const g3 = DataParser.parseSeries(document.getElementById('anovaG3')?.value || '');

      const groups = [
        { name: document.getElementById('anovaName1')?.value || 'Cohort 1', data: g1 },
        { name: document.getElementById('anovaName2')?.value || 'Cohort 2', data: g2 },
        { name: document.getElementById('anovaName3')?.value || 'Cohort 3', data: g3 }
      ];

      const res = Anova.oneWay(groups);
      if (res.error) {
        alert(res.error);
        return;
      }

      document.getElementById('anovaF').innerText = res.fStatistic.toFixed(2);
      document.getElementById('anovaP').innerText = Exporter.formatP(res.pValue);
      document.getElementById('anovaEta').innerText = res.etaSquared.toFixed(3);
      document.getElementById('anovaOmega').innerText = res.omegaSquared.toFixed(3);

      // Render Post-Hoc Pairwise Table
      const tbody = document.getElementById('anovaPostHocBody');
      if (tbody) {
        if (res.pairwise && res.pairwise.length > 0) {
          tbody.innerHTML = res.pairwise.map(p => `
            <tr>
              <td style="font-weight: 600; color: var(--text-main);">${p.comparison}</td>
              <td>${p.meanDiff >= 0 ? '+' : ''}${p.meanDiff.toFixed(2)}</td>
              <td>${p.seDiff.toFixed(2)}</td>
              <td>q = ${p.qStatistic.toFixed(2)} (t = ${p.tStatistic.toFixed(2)})</td>
              <td style="font-weight: 600; color: ${p.isSignificant ? 'var(--cyan-primary)' : 'var(--text-muted)'};">${Exporter.formatP(p.pValue)}</td>
              <td>[${p.ci95[0].toFixed(2)}, ${p.ci95[1].toFixed(2)}]</td>
              <td>${p.cohensD.toFixed(2)}</td>
              <td>
                <span class="badge ${p.isSignificant ? 'badge-sig' : 'badge-ns'}">
                  ${p.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)'}
                </span>
              </td>
            </tr>
          `).join('');
        } else {
          tbody.innerHTML = `<tr><td colspan="8" style="text-align: center; color: var(--text-dim);">No pairwise contrasts available</td></tr>`;
        }
      }

      // Detailed Clinical / Academic Summary
      let report = `A one-way between-subjects ANOVA was conducted across ${res.k} cohorts. `;
      report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} omnibus effect: F(${res.dfBetween}, ${res.dfWithin}) = ${res.fStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, η² = ${res.etaSquared.toFixed(3)}, ω² = ${res.omegaSquared.toFixed(3)}.\n\n`;

      if (res.pairwise && res.pairwise.length > 0) {
        report += `Tukey's HSD post-hoc pairwise contrasts revealed that:\n`;
        res.pairwise.forEach(p => {
          if (p.isSignificant) {
            report += `• ${p.comparison}: Statistically significant difference (ΔM = ${p.meanDiff.toFixed(2)}, 95% CI [${p.ci95[0].toFixed(2)}, ${p.ci95[1].toFixed(2)}], ${Exporter.formatP(p.pValue)}, Cohen's d = ${p.cohensD.toFixed(2)}).\n`;
          } else {
            report += `• ${p.comparison}: No statistically significant difference (ΔM = ${p.meanDiff.toFixed(2)}, 95% CI [${p.ci95[0].toFixed(2)}, ${p.ci95[1].toFixed(2)}], ${Exporter.formatP(p.pValue)}, ns).\n`;
          }
        });
      }
      document.getElementById('anovaReportText').innerText = report;

      // Render Dispersion Plot (95% CI, SEM, SD, or IQR Box & Whiskers)
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
      this.results = this.results || {}; this.results.anova = res;
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

        const res = Diagnostic.evaluate2x2(a, b, c, d);
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

        const res = Categorical.twoByTwo(a, b, c, d);
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
        exportData = {
          k: res.k,
          totalN: res.totalN,
          grandMean: res.grandMean,
          ssBetween: res.ssBetween,
          dfBetween: res.dfBetween,
          msBetween: res.msBetween,
          ssWithin: res.ssWithin,
          dfWithin: res.dfWithin,
          msWithin: res.msWithin,
          ssTotal: res.ssTotal,
          fStatistic: res.fStatistic,
          pValue: res.pValue,
          etaSquared: res.etaSquared,
          omegaSquared: res.omegaSquared,
          reportText: document.getElementById('anovaReportText')?.innerText,
          groups: res.groups,
          pairwise: res.pairwise
        };
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
          reportText: document.getElementById('teachingReportText')?.innerText
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
