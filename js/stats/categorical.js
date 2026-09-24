/**
 * Statis-Gravity - Categorical & Clinical Risk Metrics Module
 * Chi-Square test, Fisher's Exact Test, Odds Ratio (OR), Relative Risk (RR), ARR, and NNT.
 */

import { Distributions } from './distributions.js';

export const Categorical = {
  /**
   * 2x2 Contingency Table Analysis
   * Table layout:
   *              Exposed / Treat    Control / Unexposed    Total
   * Disease +          a                      b             a+b (R1)
   * Disease -          c                      d             c+d (R2)
   * Total             a+c (C1)               b+d (C2)        N
   */
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
      return { error: 'Contingency table contains empty rows or columns.' };
    }

    // Expected Frequencies
    const eA = (r1 * c1) / n;
    const eB = (r1 * c2) / n;
    const eC = (r2 * c1) / n;
    const eD = (r2 * c2) / n;

    // Pearson's Chi-Square
    const chiSqStandard =
      Math.pow(a - eA, 2) / eA +
      Math.pow(b - eB, 2) / eB +
      Math.pow(c - eC, 2) / eC +
      Math.pow(d - eD, 2) / eD;
    const pValueStandard = Distributions.chiSquarePValue(chiSqStandard, 1);

    // Yates' Continuity Corrected Chi-Square
    const chiSqYates =
      n * Math.pow(Math.max(0, Math.abs(a * d - b * c) - n / 2), 2) / (r1 * r2 * c1 * c2);
    const pValueYates = Distributions.chiSquarePValue(chiSqYates, 1);

    // Fisher's Exact Test (Hypergeometric two-tailed p-value)
    const pFisher = this.fishersExact2x2(a, b, c, d);

    // Odds Ratio (OR) & 95% Woolf CI
    // Zero cell correction if needed
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

    // Relative Risk (RR) / Risk Ratio
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

    // Absolute Risk Difference / Reduction (ARR)
    const arr = riskTreated - riskControl;
    const seARR = Math.sqrt((riskTreated * (1 - riskTreated)) / c1 + (riskControl * (1 - riskControl)) / c2);
    const arrCI95 = [arr - z95 * seARR, arr + z95 * seARR];

    // Number Needed to Treat (NNT)
    const nnt = arr !== 0 ? Math.abs(1 / arr) : Infinity;

    return {
      table: { a, b, c, d, r1, r2, c1, c2, n },
      expected: { eA, eB, eC, eD },
      chiSquare: {
        standard: chiSqStandard,
        pValueStandard,
        yates: chiSqYates,
        pValueYates,
        df: 1
      },
      fishersExact: {
        pValue: pFisher
      },
      riskMetrics: {
        oddsRatio,
        orCI95,
        relativeRisk,
        rrCI95,
        riskTreated,
        riskControl,
        arr,
        arrCI95,
        nnt
      }
    };
  },

  /**
   * Fisher's Exact Test for 2x2 table using log-gamma factorials
   */
  fishersExact2x2(a, b, c, d) {
    const r1 = a + b;
    const r2 = c + d;
    const c1 = a + c;
    const c2 = b + d;
    const n = a + b + c + d;

    // Probability of a given configuration under hypergeometric distribution
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

    const observedProb = hyperProb(a);
    const minX = Math.max(0, c1 - r2);
    const maxX = Math.min(r1, c1);

    let twoTailedP = 0;
    for (let x = minX; x <= maxX; x++) {
      const p = hyperProb(x);
      // Small tolerance for floating point imprecision
      if (p <= observedProb + 1e-12) {
        twoTailedP += p;
      }
    }
    return Math.min(1.0, Math.max(0.0, twoTailedP));
  }
};
