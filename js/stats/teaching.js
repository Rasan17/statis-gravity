/**
 * Statis-Gravity - Teaching & Simulation Suite
 * Computer-generated distributions, probability density functions, and Central Limit Theorem simulation.
 */

import { Distributions } from './distributions.js';
import { Descriptive } from './descriptive.js';

export const Teaching = {
  /**
   * Random Number Generators & Distributions
   */
  generators: {
    /**
     * Standard Normal Variate (Box-Muller Transform)
     */
    standardNormal() {
      let u = 0, v = 0;
      while (u === 0) u = Math.random();
      while (v === 0) v = Math.random();
      return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
    },

    /**
     * Normal (Gaussian) Distribution N(mean, sd)
     */
    normal(n, mean = 0, sd = 1) {
      const data = new Array(n);
      for (let i = 0; i < n; i++) {
        data[i] = mean + sd * this.standardNormal();
      }
      return data;
    },

    /**
     * Student's t-Distribution t(df)
     * t = Z / sqrt(V / df), where Z ~ N(0,1) and V ~ Chi-Square(df)
     */
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

    /**
     * Continuous Uniform Distribution U(min, max)
     */
    uniform(n, min = 0, max = 10) {
      const data = new Array(n);
      const span = max - min;
      for (let i = 0; i < n; i++) {
        data[i] = min + Math.random() * span;
      }
      return data;
    },

    /**
     * Exponential Distribution Exp(rate)
     * Inverse transform: X = -ln(1 - U) / lambda
     */
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

    /**
     * Log-Normal Distribution LogNormal(mu, sigma)
     */
    logNormal(n, mu = 1.0, sigma = 0.6) {
      const data = new Array(n);
      for (let i = 0; i < n; i++) {
        const z = this.standardNormal();
        data[i] = Math.exp(mu + sigma * z);
      }
      return data;
    },

    /**
     * Bimodal Mixture Distribution: p * N(m1, s1) + (1-p) * N(m2, s2)
     */
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

    /**
     * Poisson Distribution Pois(lambda)
     */
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

    /**
     * Chi-Square Distribution Chi2(df)
     */
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

  /**
   * Theoretical Probability Density Functions (PDF)
   */
  pdf: {
    normal(x, mean = 0, sd = 1) {
      if (sd <= 0) return 0;
      const z = (x - mean) / sd;
      return (1.0 / (sd * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
    },

    studentsT(x, df = 5, mean = 0, scale = 1) {
      if (scale <= 0 || df <= 0) return 0;
      const t = (x - mean) / scale;
      // Gamma((df+1)/2) / (sqrt(df*pi) * Gamma(df/2))
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
      // Poisson PMF: (lambda^k * exp(-lambda)) / k!
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

  /**
   * Educational distribution metadata: Clinical examples, theoretical properties & formulas
   */
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

  /**
   * Central Limit Theorem (CLT) Simulation Manager
   */
  clt: {
    // Current configuration
    populationType: 'uniform', // 'uniform' | 'exponential' | 'bimodal' | 'ushaped'
    sampleSize: 30, // n per draw
    sampleMeans: [], // accumulated means
    lastSample: [], // individual values from last draw

    populations: {
      uniform: {
        name: 'Uniform Parent U(0, 10)',
        min: 0,
        max: 10,
        mean: 5.0,
        sd: Math.sqrt(100 / 12), // ~2.8868
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
        mean: 2.0, // 1/0.5
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
        // Var = p*s1^2 + (1-p)*s2^2 + p*(1-p)*(m1-m2)^2 = 0.36 + 0.25 * 36 = 9.36 -> SD ~ 3.059
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
          // Arc-sine / Beta(0.5, 0.5) distribution
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

      // Calculate empirical metrics on sampleMeans
      const stats = Descriptive.calculate(this.sampleMeans);

      return {
        population: pop,
        sampleSize: n,
        samplesDrawn: k,
        theoreticalMean,
        theoreticalSE,
        observedMean: stats.mean,
        observedSE: stats.sd, // SD of sample means is empirical Standard Error!
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
   * Explains how two distributions overlap and separate as SD and SEM (via sample size n)
   * change, how alpha = 0.05 defines the boundary of significance, and what happens when alpha changes.
   */
  significanceOverlap: {
    /**
     * Compute two-tailed critical value for Student's t distribution at arbitrary alpha
     */
    getTCriticalValue(df, alpha = 0.05) {
      if (df <= 0) return NaN;
      const sigAlpha = Math.min(0.20, Math.max(0.0001, alpha));
      const p = 1 - sigAlpha / 2;
      const z = Distributions.invNormalCDF(p);
      if (df >= 500) return z;
      if (df === 1) return Math.tan((p - 0.5) * Math.PI);
      if (df === 2) {
        // Exact formula for df=2
        return Math.sqrt(2 / (4 * (1 - p) * p) - 2);
      }
      // Cornish-Fisher 4th-order asymptotic expansion
      const nu = df;
      const z2 = z * z, z3 = z2 * z, z5 = z3 * z2, z7 = z5 * z2, z9 = z7 * z2;
      const a = (z3 + z) / (4 * nu);
      const b = (5 * z5 + 16 * z3 + 3 * z) / (96 * nu * nu);
      const c = (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * Math.pow(nu, 3));
      const d = (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * Math.pow(nu, 4));
      return Math.max(z, z + a + b + c + d);
    },

    /**
     * Numerical computation of Weitzman's Overlap Coefficient (OVL)
     * between two normal distributions with arbitrary means and standard deviations/errors.
     */
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

    /**
     * Compute full analytical metrics for two-sample overlap and alpha boundary
     * Supports both shared and group-specific SD and SEM/n.
     */
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

      // Group 1 & Group 2 SDs
      const baseSD = parseFloat(sd) || 2.5;
      const sigma1 = Math.max(0.2, parseFloat(sd1 !== undefined && sd1 !== null ? sd1 : baseSD));
      const sigma2 = Math.max(0.2, parseFloat(sd2 !== undefined && sd2 !== null ? sd2 : baseSD));

      // Group 1 sample size & SEM
      let sampleN1, sError1;
      if (sem1 !== null && sem1 !== undefined && !isNaN(parseFloat(sem1))) {
        sError1 = Math.max(0.01, parseFloat(sem1));
        sampleN1 = Math.max(2, Math.min(1000, Math.round(Math.pow(sigma1 / sError1, 2))));
      } else {
        sampleN1 = Math.max(2, Math.round(n1 !== undefined && n1 !== null ? n1 : (n || 16)));
        sError1 = sigma1 / Math.sqrt(sampleN1);
      }

      // Group 2 sample size & SEM
      let sampleN2, sError2;
      if (sem2 !== null && sem2 !== undefined && !isNaN(parseFloat(sem2))) {
        sError2 = Math.max(0.01, parseFloat(sem2));
        sampleN2 = Math.max(2, Math.min(1000, Math.round(Math.pow(sigma2 / sError2, 2))));
      } else {
        sampleN2 = Math.max(2, Math.round(n2 !== undefined && n2 !== null ? n2 : (n || 16)));
        sError2 = sigma2 / Math.sqrt(sampleN2);
      }

      const sigAlpha = Math.min(0.20, Math.max(0.001, parseFloat(alpha) || 0.05));

      // Precision metrics with Welch-Satterthwaite approximation
      const v1 = Math.pow(sError1, 2);
      const v2 = Math.pow(sError2, 2);
      const seDiff = Math.sqrt(v1 + v2);

      const dfNum = Math.pow(v1 + v2, 2);
      const dfDenom = ((sampleN1 > 1) ? Math.pow(v1, 2) / (sampleN1 - 1) : 0) +
                      ((sampleN2 > 1) ? Math.pow(v2, 2) / (sampleN2 - 1) : 0);
      const df = dfDenom > 0 ? Math.max(1, dfNum / dfDenom) : (sampleN1 + sampleN2 - 2);

      // Critical values
      const zCrit = Distributions.invNormalCDF(1 - sigAlpha / 2);
      const tCrit = this.getTCriticalValue(df, sigAlpha);

      // Critical Difference threshold: minimum mean difference needed to reject H0 at alpha
      const deltaCrit = tCrit * seDiff;
      const deltaCritZ = zCrit * seDiff;

      // Test statistics
      const tStat = seDiff > 0 ? dMu / seDiff : 0;
      const zStat = seDiff > 0 ? dMu / seDiff : 0;
      const pValue = Distributions.tPValue(tStat, df);
      const isSignificant = pValue < sigAlpha;

      // Effect Size (Cohen's d using pooled SD)
      const pooledSD = Math.sqrt(
        ((sampleN1 - 1) * Math.pow(sigma1, 2) + (sampleN2 - 1) * Math.pow(sigma2, 2)) /
        Math.max(1, (sampleN1 + sampleN2 - 2))
      );
      const cohensD = pooledSD > 0 ? dMu / pooledSD : 0;

      // Overlap Coefficients (Weitzman's OVL for arbitrary SDs & SEMs)
      const patientOVL = this.computeOverlap(mu1, sigma1, mu2, sigma2);
      const meansOVL = this.computeOverlap(mu1, sError1, mu2, sError2);

      // Confidence Intervals for the means: mu +/- tCrit * SEM
      const moe1 = tCrit * sError1;
      const moe2 = tCrit * sError2;
      const ci1 = [mu1 - moe1, mu1 + moe1];
      const ci2 = [mu2 - moe2, mu2 + moe2];
      const ciOverlapDist = Math.max(0, ci1[1] - ci2[0]);

      // Null hypothesis bounds (H0: difference centered at 0 with SE = seDiff)
      const nullCritLeft = -deltaCrit;
      const nullCritRight = deltaCrit;

      // Pedagogical explanation synthesis
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
  }
};

