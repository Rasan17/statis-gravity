/**
 * Statis-Gravity - Sample Size & Power Calculation Module
 * A-priori statistical power and sample size determination for clinical studies.
 */

import { Distributions } from './distributions.js';

export const PowerAnalysis = {
  /**
   * Sample size for comparing two independent continuous means
   * @param {number} m1 - Mean 1
   * @param {number} m2 - Mean 2
   * @param {number} sd - Common standard deviation (or pooled SD)
   * @param {number} alpha - Type I error rate (e.g. 0.05)
   * @param {number} power - Statistical power (e.g. 0.80 or 0.90)
   * @param {number} ratio - Allocation ratio n2 / n1 (default 1)
   */
  sampleSizeMeans(m1, m2, sd, alpha = 0.05, power = 0.80, ratio = 1.0) {
    const diff = Math.abs(m1 - m2);
    if (diff === 0 || sd <= 0) {
      return { error: 'Mean difference must be greater than zero and SD positive.' };
    }

    const cohensD = diff / sd;
    const zAlpha = Distributions.invNormalCDF(1 - alpha / 2);
    const zBeta = Distributions.invNormalCDF(power);

    // Standard formula: n1 = (z_alpha + z_beta)^2 * sd^2 * (1 + 1/ratio) / diff^2
    const n1 = Math.ceil(
      Math.pow(zAlpha + zBeta, 2) * Math.pow(sd, 2) * (1 + 1 / ratio) / Math.pow(diff, 2)
    );
    const n2 = Math.ceil(n1 * ratio);

    return {
      test: 'Two Independent Means Comparison (t-test)',
      m1, m2, diff, sd,
      cohensD,
      alpha,
      power,
      ratio,
      zAlpha,
      zBeta,
      nPerGroup: ratio === 1 ? n1 : null,
      n1,
      n2,
      totalN: n1 + n2
    };
  },

  /**
   * Sample size for comparing two independent proportions
   * @param {number} p1 - Proportion in group 1 (0 < p1 < 1)
   * @param {number} p2 - Proportion in group 2 (0 < p2 < 1)
   * @param {number} alpha - Type I error (default 0.05)
   * @param {number} power - Power (default 0.80)
   */
  sampleSizeProportions(p1, p2, alpha = 0.05, power = 0.80) {
    if (p1 <= 0 || p1 >= 1 || p2 <= 0 || p2 >= 1 || p1 === p2) {
      return { error: 'Proportions must be between 0 and 1, and p1 != p2.' };
    }

    const zAlpha = Distributions.invNormalCDF(1 - alpha / 2);
    const zBeta = Distributions.invNormalCDF(power);

    const pBar = (p1 + p2) / 2;
    const qBar = 1 - pBar;

    const term1 = zAlpha * Math.sqrt(2 * pBar * qBar);
    const term2 = zBeta * Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));
    const n = Math.ceil(Math.pow(term1 + term2, 2) / Math.pow(p1 - p2, 2));

    return {
      test: 'Two Independent Proportions (Chi-Square/Z-test)',
      p1, p2,
      diff: Math.abs(p1 - p2),
      alpha,
      power,
      nPerGroup: n,
      totalN: 2 * n
    };
  }
};
