/**
 * Statis-Gravity - Clinical Diagnostic & ROC Analysis Module
 * 2x2 Diagnostic performance, Youden's J, Likelihood ratios, Empirical ROC Curve, and AUC.
 */

import { Distributions } from './distributions.js';

export const Diagnostic = {
  /**
   * Diagnostic 2x2 performance metrics
   * @param {number} tp - True Positives
   * @param {number} fp - False Positives
   * @param {number} fn - False Negatives
   * @param {number} tn - True Negatives
   */
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

    if (total === 0) {
      return { error: 'Diagnostic table cannot be completely empty.' };
    }

    const sensitivity = diseased > 0 ? tp / diseased : 0;
    const specificity = nonDiseased > 0 ? tn / nonDiseased : 0;
    const ppv = testPositive > 0 ? tp / testPositive : 0;
    const npv = testNegative > 0 ? tn / testNegative : 0;
    const accuracy = (tp + tn) / total;
    const prevalence = diseased / total;

    // Likelihood Ratios
    const plr = (1 - specificity) > 0 ? sensitivity / (1 - specificity) : Infinity;
    const nlr = specificity > 0 ? (1 - sensitivity) / specificity : Infinity;

    // Youden's J Statistic = Sensitivity + Specificity - 1
    const youdenJ = sensitivity + specificity - 1;

    // Wilson Score 95% Confidence Interval for proportions
    const wilsonCI = (p, n) => {
      if (n === 0) return [0, 0];
      const z = Distributions.invNormalCDF(0.975); // 1.95996
      const denom = 1 + (z * z) / n;
      const center = (p + (z * z) / (2 * n)) / denom;
      const margin = (z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))) / denom;
      return [Math.max(0, center - margin), Math.min(1, center + margin)];
    };

    return {
      tp, fp, fn, tn,
      diseased, nonDiseased, testPositive, testNegative, total,
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
      plr,
      nlr,
      youdenJ
    };
  },

  /**
   * Computes empirical ROC curve from continuous biomarker values and binary true status (1=diseased, 0=control)
   * @param {Array<{score: number, status: number}>} samples
   */
  computeROC(samples) {
    const valid = samples.filter(
      s => typeof s.score === 'number' && !isNaN(s.score) && (s.status === 0 || s.status === 1)
    );

    const nPos = valid.filter(s => s.status === 1).length;
    const nNeg = valid.filter(s => s.status === 0).length;

    if (nPos === 0 || nNeg === 0) {
      return { error: 'ROC requires at least one positive and one negative case.' };
    }

    // Sort descending by score
    valid.sort((a, b) => b.score - a.score);

    // Collect distinct thresholds
    const thresholds = [-Infinity];
    for (let i = 0; i < valid.length; i++) {
      thresholds.push(valid[i].score);
    }
    thresholds.push(Infinity);
    const uniqueThresholds = Array.from(new Set(thresholds)).sort((a, b) => b - a);

    const rocPoints = [];
    let bestCutoff = { threshold: uniqueThresholds[0], youdenJ: -1, sens: 0, spec: 0 };

    for (const t of uniqueThresholds) {
      let tp = 0;
      let fp = 0;
      for (const s of valid) {
        if (s.score >= t) {
          if (s.status === 1) tp++;
          else fp++;
        }
      }
      const fn = nPos - tp;
      const tn = nNeg - fp;

      const tpr = tp / nPos; // Sensitivity
      const fpr = fp / nNeg; // 1 - Specificity
      const youden = tpr - fpr;

      if (youden > bestCutoff.youdenJ) {
        bestCutoff = {
          threshold: t === -Infinity ? 'Lowest' : t === Infinity ? 'Highest' : t,
          youdenJ: youden,
          sens: tpr,
          spec: 1 - fpr
        };
      }

      rocPoints.push({ threshold: t, fpr, tpr, tp, fp, fn, tn });
    }

    // Sort points by FPR ascending, then TPR ascending
    rocPoints.sort((a, b) => a.fpr - b.fpr || a.tpr - b.tpr);

    // Calculate AUC via trapezoidal rule
    let auc = 0;
    for (let i = 1; i < rocPoints.length; i++) {
      const deltaX = rocPoints[i].fpr - rocPoints[i - 1].fpr;
      const avgY = (rocPoints[i].tpr + rocPoints[i - 1].tpr) / 2;
      auc += deltaX * avgY;
    }
    auc = Math.min(1.0, Math.max(0.0, auc));

    // Standard Error of AUC (Hanley & McNeil approximation)
    const q1 = auc / (2 - auc);
    const q2 = (2 * auc * auc) / (1 + auc);
    const seNumerator = auc * (1 - auc) + (nPos - 1) * (q1 - auc * auc) + (nNeg - 1) * (q2 - auc * auc);
    const seAUC = Math.sqrt(Math.max(0, seNumerator / (nPos * nNeg)));
    const z95 = Distributions.invNormalCDF(0.975);
    const aucCI95 = [Math.max(0, auc - z95 * seAUC), Math.min(1, auc + z95 * seAUC)];

    return {
      nPos,
      nNeg,
      total: nPos + nNeg,
      auc,
      seAUC,
      aucCI95,
      points: rocPoints,
      optimalCutoff: bestCutoff
    };
  }
};
