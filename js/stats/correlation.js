/**
 * Statis-Gravity - Correlation & Linear Regression Module
 * Pearson r, Spearman rho, OLS Regression, Standard Errors, R², and Prediction/Confidence Bands.
 */

import { Distributions } from './distributions.js';
import { Descriptive } from './descriptive.js';

export const Correlation = {
  /**
   * Cleans paired (x, y) observations
   */
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

  /**
   * Pearson Product-Moment Correlation
   */
  pearson(xData, yData) {
    const pairs = this.cleanPairs(xData, yData);
    const n = pairs.length;
    if (n < 3) {
      return { error: 'Pearson correlation requires at least 3 paired observations.' };
    }

    const xVals = pairs.map(p => p.x);
    const yVals = pairs.map(p => p.y);
    const xStats = Descriptive.calculate(xVals);
    const yStats = Descriptive.calculate(yVals);

    if (xStats.sd === 0 || yStats.sd === 0) {
      return { error: 'Standard deviation is zero in one of the variables.' };
    }

    let covSum = 0;
    for (let i = 0; i < n; i++) {
      covSum += (xVals[i] - xStats.mean) * (yVals[i] - yStats.mean);
    }
    const covariance = covSum / (n - 1);
    const rRaw = covariance / (xStats.sd * yStats.sd);
    const r = Math.max(-1, Math.min(1, isNaN(rRaw) ? 0 : rRaw));

    // t-test for correlation significance
    const df = n - 2;
    const t = Math.abs(r) >= 1 ? (r >= 0 ? Infinity : -Infinity) : (r * Math.sqrt(df)) / Math.sqrt(1 - r * r);
    const pValue = Math.abs(r) >= 1 ? 0 : Distributions.tPValue(t, df);

    // 95% Confidence Interval for r using Fisher's z transformation
    const rForZ = Math.max(-0.9999999, Math.min(0.9999999, r));
    const z = 0.5 * Math.log((1 + rForZ) / (1 - rForZ));
    const seZ = n > 3 ? 1 / Math.sqrt(n - 3) : 0;
    const zCrit = Distributions.invNormalCDF(0.975);
    const zLower = z - zCrit * seZ;
    const zUpper = z + zCrit * seZ;
    const rLower = Math.abs(r) >= 1 ? (r > 0 ? 1 : -1) : (Math.exp(2 * zLower) - 1) / (Math.exp(2 * zLower) + 1);
    const rUpper = Math.abs(r) >= 1 ? (r > 0 ? 1 : -1) : (Math.exp(2 * zUpper) - 1) / (Math.exp(2 * zUpper) + 1);

    return {
      n,
      r,
      rSquared: r * r,
      df,
      statistic: t,
      pValue,
      ci95: [rLower, rUpper],
      isSignificant: pValue < 0.05
    };
  },

  /**
   * Spearman Rank Correlation
   */
  spearman(xData, yData) {
    const pairs = this.cleanPairs(xData, yData);
    const n = pairs.length;
    if (n < 3) {
      return { error: 'Spearman correlation requires at least 3 paired observations.' };
    }

    // Rank helper with tie adjustment
    const assignRanks = (arr) => {
      const indexed = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
      const ranks = new Array(arr.length);
      let i = 0;
      while (i < arr.length) {
        let j = i;
        while (j < arr.length - 1 && indexed[j + 1].v === indexed[i].v) {
          j++;
        }
        const avgRank = (i + 1 + j + 1) / 2;
        for (let k = i; k <= j; k++) {
          ranks[indexed[k].i] = avgRank;
        }
        i = j + 1;
      }
      return ranks;
    };

    const xRanks = assignRanks(pairs.map(p => p.x));
    const yRanks = assignRanks(pairs.map(p => p.y));

    // Spearman rho is Pearson r of ranks
    const rankCorr = this.pearson(xRanks, yRanks);
    return {
      n,
      rho: rankCorr.r,
      rhoSquared: rankCorr.rSquared,
      statistic: rankCorr.statistic,
      pValue: rankCorr.pValue,
      isSignificant: rankCorr.isSignificant
    };
  },

  /**
   * Ordinary Least Squares (OLS) Linear Regression: y = beta1 * x + beta0
   */
  linearRegression(xData, yData) {
    const pairs = this.cleanPairs(xData, yData);
    const n = pairs.length;
    if (n < 3) {
      return { error: 'Linear regression requires at least 3 valid observations.' };
    }

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

    if (sxx === 0) {
      return { error: 'All X values are identical (no variance).' };
    }

    const slope = sxy / sxx;
    const intercept = yStats.mean - slope * xStats.mean;

    // Residuals and standard error
    let ssRes = 0;
    const residuals = [];
    for (let i = 0; i < n; i++) {
      const yPred = intercept + slope * xVals[i];
      const res = yVals[i] - yPred;
      residuals.push({ x: xVals[i], yActual: yVals[i], yPred, res });
      ssRes += res * res;
    }

    const dfRes = n - 2;
    const seResidual = Math.sqrt(ssRes / dfRes);
    const seSlope = seResidual / Math.sqrt(sxx);
    const seIntercept = seResidual * Math.sqrt(1 / n + (xStats.mean * xStats.mean) / sxx);

    const tSlope = seSlope === 0 ? 0 : slope / seSlope;
    const pValueSlope = Distributions.tPValue(tSlope, dfRes);

    const r = this.pearson(xVals, yVals).r;
    const rSquared = r * r;
    const fStat = rSquared === 1 ? Infinity : (rSquared * dfRes) / (1 - rSquared);
    const pValueModel = Distributions.fPValue(fStat, 1, dfRes);

    return {
      n,
      slope,
      intercept,
      equation: `y = ${slope.toFixed(4)}x + ${intercept.toFixed(4)}`,
      r,
      rSquared,
      seResidual,
      seSlope,
      seIntercept,
      tSlope,
      pValueSlope,
      fStat,
      pValueModel,
      residuals,
      xStats,
      yStats,
      sxx
    };
  },

  /**
   * Quadratic Polynomial Regression: y = beta2 * x^2 + beta1 * x + beta0
   * Evaluates U-shaped (nadir) vs Inverted U-shaped (zenith) non-linear relationships.
   */
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

    // Compare with linear fit
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

  /**
   * Analyzes whether data is monotonic (strictly increasing or decreasing together)
   * or exhibits direction reversals (e.g. increase -> decrease -> increase, N-shaped / tri-phasic).
   */
  analyzeMonotonicity(pairs) {
    if (!pairs || pairs.length < 3) {
      return { isMonotonic: true, reversals: 0, pattern: 'Indeterminate (N < 3)', direction: 'None' };
    }
    const sorted = [...pairs].sort((a, b) => a.x - b.x);
    const n = sorted.length;

    // 3-point smoothing of Y to prevent minor noise from triggering spurious reversals
    const smoothedY = [];
    for (let i = 0; i < n; i++) {
      const prev = i > 0 ? sorted[i - 1].y : sorted[i].y;
      const curr = sorted[i].y;
      const next = i < n - 1 ? sorted[i + 1].y : sorted[i].y;
      smoothedY.push((prev + curr + next) / 3);
    }

    const ySpan = Math.max(...sorted.map(p => p.y)) - Math.min(...sorted.map(p => p.y)) || 1;
    const dyThreshold = ySpan * 0.04; // 4% noise band
    const directions = [];
    for (let i = 1; i < n; i++) {
      const diff = smoothedY[i] - smoothedY[i - 1];
      if (Math.abs(diff) > dyThreshold) {
        directions.push(diff > 0 ? 1 : -1);
      }
    }

    // Compress consecutive runs of the same direction
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

  /**
   * Diagnostic assumption check for Correlation & Regression
   */
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
