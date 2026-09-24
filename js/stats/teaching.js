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
  }
};
