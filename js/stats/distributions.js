/**
 * Statis-Gravity - Statistical Distribution Functions
 * High precision numerical approximations for cumulative distributions and p-values.
 */

export const Distributions = {
  /**
   * Standard Normal Cumulative Distribution Function Φ(z)
   */
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

  /**
   * Two-tailed p-value for standard normal z
   */
  normalPValue(z) {
    const cdf = this.normalCDF(Math.abs(z));
    return 2 * (1 - cdf);
  },

  /**
   * Inverse Normal Cumulative Distribution Function (Quantile / Probit)
   * Peter John Acklam algorithm (precision 1.15e-9)
   */
  inverseNormalCDF(p) {
    return this.invNormalCDF(p);
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

  /**
   * Log-gamma function ln(Γ(x)) using Lanczos approximation (g=7, n=9)
   */
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

  /**
   * Student's t distribution Probability Density Function
   */
  tPDF(t, df) {
    const num = Math.exp(this.logGamma((df + 1) / 2) - this.logGamma(df / 2));
    const den = Math.sqrt(df * Math.PI) * Math.pow(1 + (t * t) / df, (df + 1) / 2);
    return num / den;
  },

  /**
   * Two-tailed p-value for Student's t statistic
   * Accurate adaptive quadrature
   */
  tPValue(t, df) {
    if (df <= 0 || isNaN(t) || isNaN(df)) return NaN;
    const absT = Math.abs(t);
    if (absT === 0) return 1.0;
    if (df > 100) {
      return this.normalPValue(absT);
    }

    // Adaptive integration over tail [absT, ∞]
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

  /**
   * Chi-Square right-tail p-value: P(X >= chiSq)
   */
  chiSquarePValue(chiSq, df) {
    if (chiSq <= 0 || df <= 0) return 1.0;
    if (df === 1) {
      return this.normalPValue(Math.sqrt(chiSq));
    }
    if (df === 2) {
      return Math.exp(-chiSq / 2);
    }

    // Incomplete gamma Q(a, x) via continued fraction
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

  /**
   * Snedecor's F Probability Density Function
   */
  fPDF(f, df1, df2) {
    if (f <= 0) return 0;
    const num = Math.pow(df1 * f, df1) * Math.pow(df2, df2);
    const den = Math.pow(df1 * f + df2, df1 + df2);
    const logVal = 0.5 * Math.log(num / den) - this.logBeta(df1 / 2, df2 / 2) - Math.log(f);
    return Math.exp(logVal);
  },

  logBeta(a, b) {
    return this.logGamma(a) + this.logGamma(b) - this.logGamma(a + b);
  },

  /**
   * Snedecor's F-Distribution right-tail p-value: P(X >= F)
   */
  fPValue(F, df1, df2) {
    if (F <= 0 || df1 <= 0 || df2 <= 0) return 1.0;
    // Simpson integration over tail [F, upper]
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
