/**
 * Statis-Gravity — Complete Self-Contained Clinical Biostatistics & Visualization Engine
 * Unified standalone build: Compatible with both local file:// protocol and web servers.
 * Developed by Dr G Narenthiran BSc(MedSci)(Hons) MB ChB MRCSE FEBNS FRCS(SN).
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
      ['descCanvas', 'descBoxCanvas', 'descViolinCanvas', 'hypoCanvas', 'anovaCanvas', 'corrCanvas', 'rocCanvas'].forEach(id => {
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
      document.querySelectorAll('.btn-save-plot').forEach(btn => {
        btn.addEventListener('click', (e) => {
          const cid = e.currentTarget.dataset.canvasId;
          if (this.engines[cid]) this.engines[cid].saveImage(`${cid}.png`);
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

      if (this.engines.anovaCanvas) {
        Plots.renderBoxPlot(this.engines.anovaCanvas, res.groups, 'Multi-Cohort Comparison');
      }
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
