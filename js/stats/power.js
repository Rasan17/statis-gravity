/**
 * Statis-Gravity - Sample Size & Power Calculation Module
 * A-priori sample size planning and post-hoc statistical power determination
 * Supports:
 *  1. Two Independent Means (Two-Sample t-Test)
 *  2. Paired Mean Study Design (Before vs After / Paired t-Test)
 *  3. Clinical Study Design: 2x2 Contingency Table (Chi-Square & Fisher's Exact Test)
 */

import { Distributions } from './distributions.js';

export const PowerAnalysis = {
  /**
   * 1. Two Independent Means Comparison (Two-Sample t-Test)
   * @param {Object} params
   * @param {number} params.m1 - Mean in Group 1
   * @param {number} params.m2 - Mean in Group 2
   * @param {number} params.sd - Common standard deviation (or pooled SD)
   * @param {number} params.alpha - Type I error rate (default 0.05)
   * @param {number} params.power - Target power (default 0.80)
   * @param {number|null} params.nPerGroup - If provided, calculates achieved statistical power
   * @param {number} params.ratio - Allocation ratio n2 / n1 (default 1.0)
   */
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
      // Post-hoc power calculation
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
      // A-priori sample size calculation
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

  /**
   * 2. Paired Mean Study Design (Before vs After / Matched Pairs t-Test)
   * @param {Object} params
   * @param {number} params.m1 - Baseline / Pre mean
   * @param {number} params.m2 - Post-intervention mean
   * @param {number} params.sdDiff - Standard deviation of paired differences (σ_d)
   * @param {number} params.alpha - Significance level (default 0.05)
   * @param {number} params.power - Target statistical power (default 0.80)
   * @param {number|null} params.nPairs - If provided, calculates achieved statistical power
   */
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

    const dz = diff / sdDiff; // Cohen's d_z for paired differences
    const zAlpha = Distributions.invNormalCDF(1 - alpha / 2);

    if (nPairs !== null && nPairs > 0) {
      // Post-hoc power calculation for paired design
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
        totalN: n, // N pairs = N participants in paired design
        isPostHoc: true
      };
    } else {
      // A-priori sample size calculation for paired design
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

  /**
   * 3. Clinical Study Design: 2x2 Contingency Table (Test vs Control Event Rates)
   * Evaluated by Pearson Chi-Square or Fisher's Exact Test
   * @param {Object} params
   * @param {number} params.p1 - Event rate in Test (Intervention) group
   * @param {number} params.p2 - Event rate in Control group
   * @param {number} params.alpha - Significance level (default 0.05)
   * @param {number} params.power - Desired statistical power (default 0.80)
   * @param {number|null} params.nPerGroup - If provided, calculates achieved statistical power
   * @param {string} params.testType - 'chisq' (Standard asymptotic) or 'fisher' (Continuity-corrected / Fisher's exact)
   */
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

    // Clinical effect size measures
    const arr = diff; // Absolute Risk Reduction
    const rr = p1 / p2; // Relative Risk
    const rrr = p2 > 0 ? (p2 - p1) / p2 : 0; // Relative Risk Reduction
    const or = (p1 / (1 - p1)) / (p2 / (1 - p2)); // Odds Ratio
    const nnt = 1 / arr; // Number Needed to Treat
    const cohenH = 2 * Math.abs(Math.asin(Math.sqrt(p1)) - Math.asin(Math.sqrt(p2))); // Cohen's h for proportions

    const sigma0 = Math.sqrt(2 * pBar * qBar);
    const sigma1 = Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));

    if (nPerGroup !== null && nPerGroup > 0) {
      // Post-hoc power calculation given sample size nPerGroup
      const n = Math.round(nPerGroup);

      // Standard Chi-Square power
      const zBetaChisq = (Math.sqrt(n) * diff - zAlpha * sigma0) / sigma1;
      const powerChisq = Math.min(0.9999, Math.max(0.0001, Distributions.normalCDF(zBetaChisq)));

      // Continuity-Corrected / Fisher's exact power (Casagrande-Pike-Smith adjustment)
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
      // A-priori sample size calculation
      const zBeta = Distributions.invNormalCDF(power);
      const term1 = zAlpha * Math.sqrt(2 * pBar * qBar);
      const term2 = zBeta * Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));
      const uncorrectedN = Math.ceil(Math.pow(term1 + term2, 2) / Math.pow(diff, 2));

      // Fleiss / Casagrande-Pike-Smith continuity correction for Fisher's exact test
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

  /**
   * Backwards compatible aliases
   */
  sampleSizeMeans(m1, m2, sd, alpha = 0.05, power = 0.80) {
    return this.twoIndependentMeans({ m1, m2, sd, alpha, power });
  },

  sampleSizeProportions(p1, p2, alpha = 0.05, power = 0.80) {
    return this.contingency2x2({ p1, p2, alpha, power, testType: 'chisq' });
  }
};
