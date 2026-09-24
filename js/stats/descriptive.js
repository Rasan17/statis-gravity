/**
 * Statis-Gravity - Descriptive Statistics Module
 * Calculates central tendency, dispersion, quartiles, confidence intervals, and normality metrics.
 */

import { Distributions } from './distributions.js';

export const Descriptive = {
  /**
   * Cleans an array of raw values into valid finite numbers.
   */
  cleanData(raw) {
    if (!Array.isArray(raw)) return [];
    return raw
      .map(v => (typeof v === 'number' ? v : parseFloat(v)))
      .filter(v => typeof v === 'number' && !isNaN(v) && isFinite(v));
  },

  /**
   * Computes comprehensive descriptive statistics on a dataset.
   */
  calculate(data) {
    const values = this.cleanData(data);
    const n = values.length;

    if (n === 0) {
      return { n: 0, error: 'Dataset contains no valid numeric values.' };
    }

    const sorted = [...values].sort((a, b) => a - b);
    const sum = sorted.reduce((acc, val) => acc + val, 0);
    const mean = sum / n;

    // Median
    let median;
    if (n % 2 === 0) {
      median = (sorted[n / 2 - 1] + sorted[n / 2]) / 2;
    } else {
      median = sorted[Math.floor(n / 2)];
    }

    // Min, Max, Range
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

    // Quartiles & Percentiles (Linear interpolation / R Type 7)
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

    // Outlier bounds (Tukey's fences)
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

    // Skewness & Kurtosis (Sample adjusted, Fisher-Pearson)
    let skewness = 0;
    let kurtosis = 0;
    if (n >= 3 && sd > 0) {
      skewness = (n * sumCubeDiff) / ((n - 1) * (n - 2) * Math.pow(sd, 3));
    }
    if (n >= 4 && sd > 0) {
      const term1 = (n * (n + 1) * sumQuadDiff) / ((n - 1) * (n - 2) * (n - 3) * Math.pow(sd, 4));
      const term2 = (3 * Math.pow(n - 1, 2)) / ((n - 2) * (n - 3));
      kurtosis = term1 - term2; // Excess kurtosis
    }

    // Skewness interpretation (Bulmer's rules)
    let skewnessInterpretation;
    if (skewness > 1.0) skewnessInterpretation = 'High Right Skew (Positive)';
    else if (skewness > 0.5) skewnessInterpretation = 'Moderate Right Skew';
    else if (skewness < -1.0) skewnessInterpretation = 'High Left Skew (Negative)';
    else if (skewness < -0.5) skewnessInterpretation = 'Moderate Left Skew';
    else skewnessInterpretation = 'Approximately Symmetric';

    // Kurtosis interpretation (Excess Kurtosis relative to normal = 0)
    let kurtosisInterpretation;
    if (kurtosis >= 1.0) kurtosisInterpretation = 'Leptokurtic (Heavy Tails / Peaked)';
    else if (kurtosis <= -1.0) kurtosisInterpretation = 'Platykurtic (Light Tails / Flat)';
    else kurtosisInterpretation = 'Mesokurtic (Normal Tails)';

    // 95% and 99% Confidence Interval for the Mean (Student's t critical)
    let ci95 = [mean, mean];
    let ci99 = [mean, mean];
    if (n > 1) {
      const df = n - 1;
      // t* for alpha=0.05 (two-tailed) via inverse normal or t approximation
      const z95 = Distributions.invNormalCDF(0.975);
      const z99 = Distributions.invNormalCDF(0.995);
      // Correction factor for small sample t-distribution
      const t95Factor = df > 30 ? z95 : z95 * (1 + 1 / (4 * df) + 1 / (32 * df * df));
      const t99Factor = df > 30 ? z99 : z99 * (1 + 1 / (4 * df) + 1 / (32 * df * df));

      const margin95 = t95Factor * sem;
      const margin99 = t99Factor * sem;
      ci95 = [mean - margin95, mean + margin95];
      ci99 = [mean - margin99, mean + margin99];
    }

    // Jarque-Bera Normality Test statistic & p-value
    let jbStat = 0;
    let jbPValue = 1.0;
    if (n >= 10 && sd > 0) {
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
