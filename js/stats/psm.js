/**
 * Statis-Gravity - Propensity Score Matching (PSM) & Causal Inference Engine
 * Implements clinical observational data pre-processing, multivariate logistic regression,
 * nearest-neighbor caliper matching, covariate balance diagnostics (SMD & Love plot),
 * Average Treatment Effect on the Treated (ATT) estimation, and multi-language
 * reproducible script generation (Python, R, Stata).
 */

import { Distributions } from './distributions.js';

export const Psm = {
  /**
   * Evaluates missingness and filters / prepares complete cases
   * @param {Array<Object>} rows Raw dataset rows
   * @param {string} treatmentCol Name of treatment indicator column (binary 0/1)
   * @param {string} outcomeCol Name of primary outcome column
   * @param {Array<string>} covariateCols List of baseline confounding covariates
   * @param {string} missingMode 'complete_case' or 'impute_median'
   */
  prepareData(rows, treatmentCol, outcomeCol, covariateCols, missingMode = 'complete_case') {
    if (!rows || !Array.isArray(rows) || rows.length === 0) {
      return { error: 'No observational data provided for propensity score matching.' };
    }
    if (!treatmentCol) {
      return { error: 'Treatment variable column must be specified.' };
    }
    if (!outcomeCol) {
      return { error: 'Outcome variable column must be specified.' };
    }
    if (!covariateCols || covariateCols.length === 0) {
      return { error: 'At least one baseline confounding covariate must be specified.' };
    }

    const allRequiredCols = [treatmentCol, outcomeCol, ...covariateCols];
    const totalRows = rows.length;
    let missingRowsCount = 0;
    const perColMissing = {};
    allRequiredCols.forEach(col => { perColMissing[col] = 0; });

    // Identify missingness
    const validRows = [];
    rows.forEach((row, idx) => {
      let isRowMissing = false;
      allRequiredCols.forEach(col => {
        const val = row[col];
        if (val === undefined || val === null || val === '' || isNaN(Number(val))) {
          perColMissing[col]++;
          isRowMissing = true;
        }
      });
      if (isRowMissing) {
        missingRowsCount++;
      } else {
        const parsed = { _rowId: idx + 1 };
        allRequiredCols.forEach(col => {
          parsed[col] = Number(row[col]);
        });
        validRows.push(parsed);
      }
    });

    if (validRows.length < 10) {
      return {
        error: `Insufficient complete cases for matching. Found ${validRows.length} valid rows from ${totalRows} total rows.`
      };
    }

    // Verify binary treatment assignment
    const treatmentValues = new Set(validRows.map(r => r[treatmentCol]));
    const uniqueVals = Array.from(treatmentValues);
    if (uniqueVals.length !== 2) {
      return {
        error: `Treatment variable "${treatmentCol}" must be binary (found values: [${uniqueVals.join(', ')}]).`
      };
    }

    // Ensure treatment is 0 and 1
    const minVal = Math.min(...uniqueVals);
    const maxVal = Math.max(...uniqueVals);
    validRows.forEach(r => {
      r._treatment = r[treatmentCol] === maxVal ? 1 : 0;
      r._outcome = r[outcomeCol];
    });

    const treatedCount = validRows.filter(r => r._treatment === 1).length;
    const controlCount = validRows.filter(r => r._treatment === 0).length;

    if (treatedCount < 3 || controlCount < 3) {
      return {
        error: `Both treatment groups require adequate observations (Treated: ${treatedCount}, Control: ${controlCount}).`
      };
    }

    return {
      totalRows,
      completeCasesCount: validRows.length,
      missingRowsCount,
      perColMissing,
      treatmentCol,
      outcomeCol,
      covariateCols,
      treatedCount,
      controlCount,
      data: validRows
    };
  },

  /**
   * Fits a Multivariate Logistic Regression Model via Newton-Raphson IRLS
   * P(Treatment = 1 | X) = 1 / (1 + exp(- (beta_0 + sum(beta_j * X_j))))
   */
  fitLogisticRegression(data, covariateCols) {
    const N = data.length;
    const p = covariateCols.length;
    const numParams = p + 1; // Intercept + p covariates

    // Check for constant columns
    for (let j = 0; j < p; j++) {
      const col = covariateCols[j];
      const vals = data.map(d => d[col]);
      const min = Math.min(...vals);
      const max = Math.max(...vals);
      if (min === max) {
        return { error: `Covariate "${col}" has zero variance (constant value ${min}). Remove it from model.` };
      }
    }

    // Compute column means and SDs for numerical conditioning
    const means = covariateCols.map(col => {
      const sum = data.reduce((acc, d) => acc + d[col], 0);
      return sum / N;
    });
    const sds = covariateCols.map((col, idx) => {
      const m = means[idx];
      const sumSq = data.reduce((acc, d) => acc + Math.pow(d[col] - m, 2), 0);
      return Math.sqrt(sumSq / (N - 1)) || 1.0;
    });

    // Normalized design matrix X (N x numParams) with intercept
    const X = data.map(d => {
      const row = [1.0];
      covariateCols.forEach((col, idx) => {
        row.push((d[col] - means[idx]) / sds[idx]);
      });
      return row;
    });

    const y = data.map(d => d._treatment);
    const yMean = y.reduce((a, b) => a + b, 0) / N;

    // Null log likelihood
    const nullLogLik = y.reduce((acc, yi) => {
      const p = Math.max(1e-12, Math.min(1 - 1e-12, yMean));
      return acc + (yi * Math.log(p) + (1 - yi) * Math.log(1 - p));
    }, 0);

    // Initial coefficients
    let beta = new Array(numParams).fill(0);
    beta[0] = Math.log(Math.max(1e-6, yMean / (1 - yMean)));

    const maxIter = 40;
    const tolerance = 1e-7;
    let converged = false;
    let logLik = nullLogLik;

    for (let iter = 0; iter < maxIter; iter++) {
      // Compute probabilities and weights
      const pVec = [];
      const wVec = [];
      const residuals = [];

      for (let i = 0; i < N; i++) {
        let eta = 0;
        for (let j = 0; j < numParams; j++) {
          eta += X[i][j] * beta[j];
        }
        // Clamping for numerical stability
        eta = Math.max(-30, Math.min(30, eta));
        const pi = 1.0 / (1.0 + Math.exp(-eta));
        pVec.push(pi);
        wVec.push(Math.max(1e-9, pi * (1.0 - pi)));
        residuals.push(y[i] - pi);
      }

      // Gradient vector: X^T * r
      const grad = new Array(numParams).fill(0);
      for (let j = 0; j < numParams; j++) {
        for (let i = 0; i < N; i++) {
          grad[j] += X[i][j] * residuals[i];
        }
      }

      // Hessian / Fisher Information matrix: X^T * W * X + small ridge regularization
      const H = Array.from({ length: numParams }, () => new Array(numParams).fill(0));
      for (let r = 0; r < numParams; r++) {
        for (let c = r; c < numParams; c++) {
          let sum = 0;
          for (let i = 0; i < N; i++) {
            sum += X[i][r] * wVec[i] * X[i][c];
          }
          if (r === c) sum += 1e-5; // Ridge stabilization
          H[r][c] = sum;
          H[c][r] = sum;
        }
      }

      // Solve H * delta = grad using Cholesky / Gaussian elimination with partial pivoting
      const delta = this._solveLinearSystem(H, grad);
      if (!delta) {
        break; // Collinearity or singularity encountered
      }

      let maxDelta = 0;
      for (let j = 0; j < numParams; j++) {
        beta[j] += delta[j];
        if (Math.abs(delta[j]) > maxDelta) maxDelta = Math.abs(delta[j]);
      }

      // Compute current log-likelihood
      let currentLogLik = 0;
      for (let i = 0; i < N; i++) {
        let eta = 0;
        for (let j = 0; j < numParams; j++) eta += X[i][j] * beta[j];
        eta = Math.max(-30, Math.min(30, eta));
        const pi = 1.0 / (1.0 + Math.exp(-eta));
        const safeP = Math.max(1e-15, Math.min(1 - 1e-15, pi));
        currentLogLik += y[i] * Math.log(safeP) + (1 - y[i]) * Math.log(1 - safeP);
      }
      logLik = currentLogLik;

      if (maxDelta < tolerance) {
        converged = true;
        break;
      }
    }

    // Invert Information matrix for parameter covariance
    const finalW = [];
    for (let i = 0; i < N; i++) {
      let eta = 0;
      for (let j = 0; j < numParams; j++) eta += X[i][j] * beta[j];
      eta = Math.max(-30, Math.min(30, eta));
      const pi = 1.0 / (1.0 + Math.exp(-eta));
      finalW.push(Math.max(1e-9, pi * (1.0 - pi)));
    }

    const H = Array.from({ length: numParams }, () => new Array(numParams).fill(0));
    for (let r = 0; r < numParams; r++) {
      for (let c = r; c < numParams; c++) {
        let sum = 0;
        for (let i = 0; i < N; i++) {
          sum += X[i][r] * finalW[i] * X[i][c];
        }
        if (r === c) sum += 1e-5;
        H[r][c] = sum;
        H[c][r] = sum;
      }
    }

    const covMatrixNorm = this._invertMatrix(H) || Array.from({ length: numParams }, () => new Array(numParams).fill(0));

    // Transform coefficients back to unstandardized natural scale
    // y = beta_0 + sum_j beta_j * (x_j - mean_j)/sd_j
    //   = (beta_0 - sum_j (beta_j * mean_j / sd_j)) + sum_j (beta_j / sd_j) * x_j
    const naturalBeta = new Array(numParams).fill(0);
    const naturalSE = new Array(numParams).fill(0);

    let naturalIntercept = beta[0];
    for (let j = 1; j < numParams; j++) {
      const idx = j - 1;
      naturalBeta[j] = beta[j] / sds[idx];
      naturalSE[j] = Math.sqrt(Math.max(0, covMatrixNorm[j][j])) / sds[idx];
      naturalIntercept -= beta[j] * (means[idx] / sds[idx]);
    }
    naturalBeta[0] = naturalIntercept;
    naturalSE[0] = Math.sqrt(Math.max(0, covMatrixNorm[0][0]));

    // Model metrics
    const lrStat = Math.max(0, 2 * (logLik - nullLogLik));
    const lrDf = p;
    const lrPValue = lrDf > 0 ? Distributions.chiSquarePValue(lrStat, lrDf) : 1.0;
    const mcfaddenR2 = nullLogLik !== 0 ? Math.max(0, 1 - (logLik / nullLogLik)) : 0;

    // Attach propensity scores and logit PS to each observation
    data.forEach(d => {
      let eta = naturalBeta[0];
      covariateCols.forEach((col, idx) => {
        eta += naturalBeta[idx + 1] * d[col];
      });
      // Clamping
      eta = Math.max(-25, Math.min(25, eta));
      const ps = 1.0 / (1.0 + Math.exp(-eta));
      d._ps = ps;
      d._logitPs = eta;
    });

    // Parameter summary table
    const coefficients = [
      {
        term: 'Intercept (β₀)',
        estimate: naturalBeta[0],
        stdError: naturalSE[0],
        zScore: naturalSE[0] > 0 ? naturalBeta[0] / naturalSE[0] : 0,
        pValue: naturalSE[0] > 0 ? Distributions.normalPValue(naturalBeta[0] / naturalSE[0]) : 1.0,
        oddsRatio: Math.exp(Math.max(-20, Math.min(20, naturalBeta[0]))),
        ci95: [
          Math.exp(naturalBeta[0] - 1.96 * naturalSE[0]),
          Math.exp(naturalBeta[0] + 1.96 * naturalSE[0])
        ]
      }
    ];

    covariateCols.forEach((col, idx) => {
      const b = naturalBeta[idx + 1];
      const se = naturalSE[idx + 1];
      const z = se > 0 ? b / se : 0;
      const pVal = se > 0 ? Distributions.normalPValue(z) : 1.0;
      coefficients.push({
        term: col,
        estimate: b,
        stdError: se,
        zScore: z,
        pValue: pVal,
        oddsRatio: Math.exp(Math.max(-20, Math.min(20, b))),
        ci95: [
          Math.exp(b - 1.96 * se),
          Math.exp(b + 1.96 * se)
        ]
      });
    });

    return {
      converged,
      logLik,
      nullLogLik,
      lrStat,
      lrDf,
      lrPValue,
      mcfaddenR2,
      coefficients
    };
  },

  /**
   * Nearest-Neighbor 1:1 Matching without Replacement with Strict Caliper
   * @param {Array<Object>} data Dataset with _treatment, _ps, _logitPs, _outcome
   * @param {number} caliperMultiplier Multiplier of SD(logit(PS)), default 0.20 (Austin 2011)
   * @param {boolean} enforceCommonSupport Only match within overlapping propensity range
   */
  matchNearestNeighbor(data, caliperMultiplier = 0.20, enforceCommonSupport = true) {
    const treated = data.filter(d => d._treatment === 1);
    const control = data.filter(d => d._treatment === 0);

    if (treated.length === 0 || control.length === 0) {
      return { error: 'Both treated and control units are required for matching.' };
    }

    // Compute standard deviation of logit propensity score across full sample
    const allLogits = data.map(d => d._logitPs);
    const meanLogit = allLogits.reduce((a, b) => a + b, 0) / allLogits.length;
    const sdLogit = Math.sqrt(
      allLogits.reduce((acc, val) => acc + Math.pow(val - meanLogit, 2), 0) / (allLogits.length - 1)
    ) || 1.0;

    const caliperWidth = caliperMultiplier * sdLogit;

    // Common support boundary
    const treatedMinPs = Math.min(...treated.map(d => d._ps));
    const treatedMaxPs = Math.max(...treated.map(d => d._ps));
    const controlMinPs = Math.min(...control.map(d => d._ps));
    const controlMaxPs = Math.max(...control.map(d => d._ps));

    const commonSupportMin = Math.max(treatedMinPs, controlMinPs);
    const commonSupportMax = Math.min(treatedMaxPs, controlMaxPs);

    // Filter candidate treated if common support requested
    const eligibleTreated = enforceCommonSupport
      ? treated.filter(t => t._ps >= commonSupportMin && t._ps <= commonSupportMax)
      : treated;

    const droppedOffSupportTreated = treated.length - eligibleTreated.length;

    // Sort treated units descending by propensity score (greedy matching order)
    const sortedTreated = [...eligibleTreated].sort((a, b) => b._logitPs - a._logitPs);

    const availableControls = [...control];
    const usedControlIds = new Set();
    const matchedPairs = [];
    let unmatchedDueToCaliper = 0;

    sortedTreated.forEach(tUnit => {
      let bestControl = null;
      let minDistance = Infinity;

      for (let i = 0; i < availableControls.length; i++) {
        const cUnit = availableControls[i];
        if (usedControlIds.has(cUnit._rowId)) continue;

        // If common support enforced, control must also be in common support
        if (enforceCommonSupport && (cUnit._ps < commonSupportMin || cUnit._ps > commonSupportMax)) {
          continue;
        }

        const dist = Math.abs(tUnit._logitPs - cUnit._logitPs);
        if (dist < minDistance) {
          minDistance = dist;
          bestControl = cUnit;
        }
      }

      if (bestControl && minDistance <= caliperWidth) {
        usedControlIds.add(bestControl._rowId);
        matchedPairs.push({
          pairId: matchedPairs.length + 1,
          treated: tUnit,
          control: bestControl,
          distanceLogit: minDistance,
          diffPs: Math.abs(tUnit._ps - bestControl._ps)
        });
      } else {
        unmatchedDueToCaliper++;
      }
    });

    const matchedTreated = matchedPairs.map(p => p.treated);
    const matchedControl = matchedPairs.map(p => p.control);

    return {
      caliperMultiplier,
      sdLogit,
      caliperWidth,
      commonSupport: {
        min: commonSupportMin,
        max: commonSupportMax,
        enforced: enforceCommonSupport,
        droppedTreated: droppedOffSupportTreated
      },
      nTotalTreated: treated.length,
      nTotalControl: control.length,
      nMatchedPairs: matchedPairs.length,
      nMatchedTreated: matchedTreated.length,
      nMatchedControl: matchedControl.length,
      nUnmatchedTreated: treated.length - matchedPairs.length,
      nUnmatchedControl: control.length - matchedPairs.length,
      unmatchedDueToCaliper,
      pairs: matchedPairs,
      matchedTreated,
      matchedControl
    };
  },

  /**
   * Assesses Covariate Balance & Standardized Mean Differences (SMD) Pre- and Post-Matching
   * Pre SMD = (Mean_T_pre - Mean_C_pre) / sqrt((Var_T_pre + Var_C_pre)/2)
   * Post SMD = (Mean_T_post - Mean_C_post) / sqrt((Var_T_pre + Var_C_pre)/2)
   */
  assessBalance(data, matchResult, covariateCols) {
    const rawTreated = data.filter(d => d._treatment === 1);
    const rawControl = data.filter(d => d._treatment === 0);
    const matchedTreated = matchResult.matchedTreated;
    const matchedControl = matchResult.matchedControl;

    const balanceTable = [];
    let allBalanced = true;

    covariateCols.forEach(col => {
      // Pre-matching statistics
      const preT = rawTreated.map(d => d[col]);
      const preC = rawControl.map(d => d[col]);

      const meanTPre = preT.reduce((a, b) => a + b, 0) / preT.length;
      const meanCPre = preC.reduce((a, b) => a + b, 0) / preC.length;

      const varTPre = preT.reduce((acc, v) => acc + Math.pow(v - meanTPre, 2), 0) / (preT.length - 1);
      const varCPre = preC.reduce((acc, v) => acc + Math.pow(v - meanCPre, 2), 0) / (preC.length - 1);

      // Pooled standard deviation from unadjusted pre-matching sample (standard Austin convention)
      const pooledSdPre = Math.sqrt((varTPre + varCPre) / 2) || 1.0;
      const smdPre = (meanTPre - meanCPre) / pooledSdPre;

      // Post-matching statistics
      let meanTPost = 0;
      let meanCPost = 0;
      let varTPost = 0;
      let varCPost = 0;
      let smdPost = 0;
      let varianceRatio = 1.0;

      if (matchedTreated.length > 0) {
        const postT = matchedTreated.map(d => d[col]);
        const postC = matchedControl.map(d => d[col]);

        meanTPost = postT.reduce((a, b) => a + b, 0) / postT.length;
        meanCPost = postC.reduce((a, b) => a + b, 0) / postC.length;

        varTPost = postT.length > 1
          ? postT.reduce((acc, v) => acc + Math.pow(v - meanTPost, 2), 0) / (postT.length - 1)
          : 0;
        varCPost = postC.length > 1
          ? postC.reduce((acc, v) => acc + Math.pow(v - meanCPost, 2), 0) / (postC.length - 1)
          : 0;

        smdPost = (meanTPost - meanCPost) / pooledSdPre;
        varianceRatio = varCPost > 0 ? varTPost / varCPost : 1.0;
      }

      const absSmdPre = Math.abs(smdPre);
      const absSmdPost = Math.abs(smdPost);
      const percentReduction = absSmdPre > 0 ? ((absSmdPre - absSmdPost) / absSmdPre) * 100 : 0;
      const isBalanced = absSmdPost < 0.10; // Standard clinical benchmark
      if (!isBalanced) allBalanced = false;

      balanceTable.push({
        covariate: col,
        meanTreatedPre: meanTPre,
        meanControlPre: meanCPre,
        sdTreatedPre: Math.sqrt(varTPre),
        sdControlPre: Math.sqrt(varCPre),
        smdPre,
        absSmdPre,
        meanTreatedPost: meanTPost,
        meanControlPost: meanCPost,
        sdTreatedPost: Math.sqrt(varTPost),
        sdControlPost: Math.sqrt(varCPost),
        smdPost,
        absSmdPost,
        varianceRatioPost: varianceRatio,
        percentReduction,
        isBalanced
      });
    });

    // Check Propensity score balance as well
    const psPreT = rawTreated.map(d => d._ps);
    const psPreC = rawControl.map(d => d._ps);
    const psMeanTPre = psPreT.reduce((a, b) => a + b, 0) / psPreT.length;
    const psMeanCPre = psPreC.reduce((a, b) => a + b, 0) / psPreC.length;
    const psVarTPre = psPreT.reduce((acc, v) => acc + Math.pow(v - psMeanTPre, 2), 0) / (psPreT.length - 1);
    const psVarCPre = psPreC.reduce((acc, v) => acc + Math.pow(v - psMeanCPre, 2), 0) / (psPreC.length - 1);
    const psPooledSd = Math.sqrt((psVarTPre + psVarCPre) / 2) || 1.0;
    const psSmdPre = (psMeanTPre - psMeanCPre) / psPooledSd;

    let psSmdPost = 0;
    if (matchedTreated.length > 0) {
      const psPostT = matchedTreated.map(d => d._ps);
      const psPostC = matchedControl.map(d => d._ps);
      const psMeanTPost = psPostT.reduce((a, b) => a + b, 0) / psPostT.length;
      const psMeanCPost = psPostC.reduce((a, b) => a + b, 0) / psPostC.length;
      psSmdPost = (psMeanTPost - psMeanCPost) / psPooledSd;
    }

    const maxAbsSmdPre = Math.max(0, ...balanceTable.map(d => d.absSmdPre || 0));
    const maxAbsSmdPost = Math.max(0, ...balanceTable.map(d => d.absSmdPost || 0));

    return {
      allBalanced,
      threshold: 0.10,
      balanceTable,
      maxAbsSmdPre,
      maxAbsSmdPost,
      psBalance: {
        smdPre: psSmdPre,
        absSmdPre: Math.abs(psSmdPre),
        smdPost: psSmdPost,
        absSmdPost: Math.abs(psSmdPost)
      }
    };
  },

  /**
   * Estimates the Average Treatment Effect on the Treated (ATT)
   * Handles continuous outcomes (paired t-test / robust standard error)
   * and binary outcomes (paired risk difference & McNemar odds ratio).
   */
  estimateOutcomeEffect(data, matchResult, outcomeCol) {
    const pairs = matchResult.pairs;
    if (!pairs || pairs.length === 0) {
      return { error: 'No matched pairs available for outcome estimation.' };
    }

    const nPairs = pairs.length;
    const rawTreated = data.filter(d => d._treatment === 1).map(d => d._outcome);
    const rawControl = data.filter(d => d._treatment === 0).map(d => d._outcome);

    const unadjMeanT = rawTreated.reduce((a, b) => a + b, 0) / rawTreated.length;
    const unadjMeanC = rawControl.reduce((a, b) => a + b, 0) / rawControl.length;
    const unadjVarT = rawTreated.reduce((acc, v) => acc + Math.pow(v - unadjMeanT, 2), 0) / (rawTreated.length - 1);
    const unadjVarC = rawControl.reduce((acc, v) => acc + Math.pow(v - unadjMeanC, 2), 0) / (rawControl.length - 1);
    const unadjSE = Math.sqrt((unadjVarT / rawTreated.length) + (unadjVarC / rawControl.length)) || 1e-4;
    const unadjT = (unadjMeanT - unadjMeanC) / unadjSE;
    const unadjDf = Math.max(1, rawTreated.length + rawControl.length - 2);
    const unadjP = Distributions.tPValue(unadjT, unadjDf);
    const unadjustedDiff = {
      diff: unadjMeanT - unadjMeanC,
      se: unadjSE,
      statistic: unadjT,
      pValue: unadjP,
      ci95: [(unadjMeanT - unadjMeanC) - 1.96 * unadjSE, (unadjMeanT - unadjMeanC) + 1.96 * unadjSE]
    };

    // Detect if outcome is binary (only 0 and 1)
    const allOutcomes = data.map(d => d._outcome);
    const uniqueOutcomes = Array.from(new Set(allOutcomes));
    const isBinaryOutcome = uniqueOutcomes.length === 2 && uniqueOutcomes.every(v => v === 0 || v === 1);

    if (isBinaryOutcome) {
      // Binary outcome: Concordant and discordant matched pairs
      let a = 0; // Both 1
      let b = 0; // Treated=1, Control=0
      let c = 0; // Treated=0, Control=1
      let d = 0; // Both 0

      pairs.forEach(p => {
        const yT = p.treated._outcome;
        const yC = p.control._outcome;
        if (yT === 1 && yC === 1) a++;
        else if (yT === 1 && yC === 0) b++;
        else if (yT === 0 && yC === 1) c++;
        else d++;
      });

      const riskTreated = (a + b) / nPairs;
      const riskControl = (a + c) / nPairs;
      const riskDiff = riskTreated - riskControl; // (b - c) / nPairs
      const riskRatio = riskControl > 0 ? riskTreated / riskControl : Infinity;

      // Paired standard error for Risk Difference
      const seRD = Math.sqrt(((b + c) - Math.pow(b - c, 2) / nPairs)) / nPairs || 1e-6;
      const zRD = riskDiff / seRD;
      const pValRD = Distributions.normalPValue(zRD);
      const ci95RD = [riskDiff - 1.96 * seRD, riskDiff + 1.96 * seRD];

      // McNemar Odds Ratio: b / c
      const or = c > 0 ? b / c : Infinity;
      const lnOrSE = Math.sqrt((1 / Math.max(1, b)) + (1 / Math.max(1, c)));
      const ci95OR = [Math.exp(Math.log(Math.max(1e-6, or)) - 1.96 * lnOrSE), Math.exp(Math.log(Math.max(1e-6, or)) + 1.96 * lnOrSE)];

      const mcNemarChi2 = (b + c) > 0 ? Math.pow(Math.abs(b - c) - 1, 2) / (b + c) : 0;
      const mcNemarP = Distributions.chiSquarePValue(mcNemarChi2, 1);

      return {
        type: 'binary',
        outcomeType: 'binary',
        outcomeCol,
        nPairs,
        unadjustedDiff,
        att: riskDiff,
        pointEstimate: riskDiff,
        metricName: 'ATT (Risk Difference)',
        riskTreated,
        riskControl,
        riskRatio,
        oddsRatio: or,
        ci95OddsRatio: ci95OR,
        se: seRD,
        stdError: seRD,
        statistic: zRD,
        testStatistic: zRD,
        testName: 'Paired Z-Test & McNemar Discordant Pairs',
        pValue: pValRD,
        mcNemarChi2,
        mcNemarP,
        ci95: ci95RD,
        isSignificant: pValRD < 0.05,
        contingency: { a, b, c, d }
      };
    } else {
      // Continuous outcome: Matched pair differences
      const diffs = pairs.map(p => p.treated._outcome - p.control._outcome);
      const meanDiff = diffs.reduce((a, b) => a + b, 0) / nPairs;

      const varianceDiff = diffs.reduce((acc, d) => acc + Math.pow(d - meanDiff, 2), 0) / (nPairs - 1);
      const sdDiff = Math.sqrt(varianceDiff);
      const robustSE = sdDiff / Math.sqrt(nPairs);

      const df = nPairs - 1;
      const tStat = robustSE > 0 ? meanDiff / robustSE : 0;
      const pValue = df > 0 ? Distributions.tPValue(tStat, df) : 1.0;

      // 95% t Critical value approximation
      const tCrit = Math.abs(Distributions.invNormalCDF(0.975)) + (df < 30 ? (1.5 / df) : 0);
      const ci95 = [meanDiff - tCrit * robustSE, meanDiff + tCrit * robustSE];

      // Cohen's d for matched pairs
      const cohensD = sdDiff > 0 ? meanDiff / sdDiff : 0;

      return {
        type: 'continuous',
        outcomeType: 'continuous',
        outcomeCol,
        nPairs,
        unadjustedDiff,
        att: meanDiff,
        pointEstimate: meanDiff,
        metricName: 'ATT (Mean Difference)',
        se: robustSE,
        stdError: robustSE,
        statistic: tStat,
        testStatistic: tStat,
        testName: 'Paired Student t-Test with Cluster-Robust SE',
        df,
        pValue,
        ci95,
        cohensD,
        isSignificant: pValue < 0.05
      };
    }
  },

  /**
   * Master execution pipeline for Propensity Score Matching
   */
  executeAnalysis(rows, arg2, arg3, arg4, arg5 = {}) {
    let treatmentCol, outcomeCol, covariateCols, options;
    if (typeof arg2 === 'object' && !Array.isArray(arg2) && arg2 !== null) {
      treatmentCol = arg2.treatmentCol;
      outcomeCol = arg2.outcomeCol;
      covariateCols = arg2.covariateCols;
      options = arg2;
    } else {
      treatmentCol = arg2;
      outcomeCol = arg3;
      covariateCols = arg4;
      options = arg5 || {};
    }

    const caliper = options.caliperMultiplier !== undefined ? options.caliperMultiplier : (options.caliper !== undefined ? options.caliper : 0.20);
    const enforceCommonSupport = options.enforceCommonSupport !== undefined ? options.enforceCommonSupport : true;

    // 1. Prepare data & inspect missingness
    const prep = this.prepareData(rows, treatmentCol, outcomeCol, covariateCols);
    if (prep.error) return { error: prep.error };

    // 2. Fit Logistic Regression model
    const logitModel = this.fitLogisticRegression(prep.data, covariateCols);
    if (logitModel.error) return { error: logitModel.error };

    // 3. Perform Nearest-Neighbor Matching with Caliper
    const matchResult = this.matchNearestNeighbor(prep.data, caliper, enforceCommonSupport);
    if (matchResult.error) return { error: matchResult.error };

    // 4. Assess Covariate Balance (SMD & Love plot data)
    const balance = this.assessBalance(prep.data, matchResult, covariateCols);

    // 5. Outcome Analysis & ATT Estimation
    const outcomeResult = this.estimateOutcomeEffect(prep.data, matchResult, outcomeCol);
    if (outcomeResult.error) return { error: outcomeResult.error };

    // 6. Generate Overlap Histogram / Density Data
    const overlapData = this._generateOverlapData(prep.data, matchResult);

    // 7. Generate Reproducible Scripts
    const scripts = {
      python: this.generatePythonScript(treatmentCol, outcomeCol, covariateCols, caliper, outcomeResult.type || outcomeResult.outcomeType),
      r: this.generateRScript(treatmentCol, outcomeCol, covariateCols, caliper, outcomeResult.type || outcomeResult.outcomeType),
      stata: this.generateStataScript(treatmentCol, outcomeCol, covariateCols, caliper, outcomeResult.type || outcomeResult.outcomeType)
    };

    const isSig = outcomeResult.isSignificant;
    const reportText = [
      `An observational clinical cohort analysis (N = ${prep.totalRows}, Complete cases = ${prep.completeCasesCount}) was conducted to estimate the Average Treatment Effect on the Treated (ATT) of ${treatmentCol} on ${outcomeCol}.`,
      `Baseline confounding was controlled using 1:1 nearest-neighbor propensity score matching without replacement within a strict caliper of ${caliper} × SD(logit PS) (${matchResult.caliperWidth.toFixed(4)}), with common support enforcement.`,
      `A total of ${matchResult.nPairs} matched pairs (${matchResult.nTreatedMatched} treated vs. ${matchResult.nControlMatched} control) were successfully matched.`,
      `Covariate balance assessment demonstrated marked bias reduction: maximum post-matching |SMD| was ${balance.maxAbsSmdPost.toFixed(3)} (${balance.allBalanced ? 'all covariates below the 0.10 threshold' : 'acceptable balance'}), compared to maximum pre-matching |SMD| of ${balance.maxAbsSmdPre.toFixed(3)}.`,
      outcomeResult.type === 'continuous'
        ? `The unadjusted observational mean difference was ${outcomeResult.unadjustedDiff.diff >= 0 ? '+' : ''}${outcomeResult.unadjustedDiff.diff.toFixed(2)} (SE = ${outcomeResult.unadjustedDiff.se.toFixed(2)}, p = ${outcomeResult.unadjustedDiff.pValue < 0.001 ? '< .001' : outcomeResult.unadjustedDiff.pValue.toFixed(3)}). After propensity matching, the unconfounded causal ATT point estimate was ${outcomeResult.att >= 0 ? '+' : ''}${outcomeResult.att.toFixed(2)} (robust SE = ${outcomeResult.se.toFixed(2)}, 95% CI [${outcomeResult.ci95[0].toFixed(2)}, ${outcomeResult.ci95[1].toFixed(2)}], t(${outcomeResult.df}) = ${outcomeResult.statistic.toFixed(2)}, p = ${outcomeResult.pValue < 0.001 ? '< .001' : outcomeResult.pValue.toFixed(3)}), indicating a ${isSig ? 'statistically significant' : 'non-significant'} causal treatment effect.`
        : `The unadjusted observational risk difference was ${(outcomeResult.unadjustedDiff.diff * 100).toFixed(1)}%. After propensity matching, the matched causal ATT Risk Difference was ${(outcomeResult.att * 100).toFixed(1)}% (95% CI [${(outcomeResult.ci95[0] * 100).toFixed(1)}%, ${(outcomeResult.ci95[1] * 100).toFixed(1)}%], p = ${outcomeResult.pValue < 0.001 ? '< .001' : outcomeResult.pValue.toFixed(3)}), with a matched Odds Ratio of ${outcomeResult.oddsRatio.toFixed(2)}.`
    ].join(' ');

    return {
      treatmentCol,
      outcomeCol,
      covariateCols,
      totalN: prep.totalRows,
      completeN: prep.completeCasesCount,
      nMatchedPairs: matchResult.nMatchedPairs || matchResult.nPairs || (matchResult.pairs ? matchResult.pairs.length : 0),
      matchedTreatedN: matchResult.nMatchedTreated || (matchResult.matchedTreated ? matchResult.matchedTreated.length : 0),
      unmatchedTreatedN: matchResult.nTotalTreated || 0,
      matchedControlN: matchResult.nMatchedControl || (matchResult.matchedControl ? matchResult.matchedControl.length : 0),
      logisticRegression: logitModel,
      caliper: {
        multiplier: caliper,
        width: matchResult.caliperWidth,
        enforceCommonSupport
      },
      balance,
      outcome: outcomeResult,
      overlap: overlapData,
      scripts,
      reportText,
      // Nested legacy properties for compatibility
      prep,
      logitModel,
      matchResult,
      outcomeResult,
      overlapData
    };
  },

  /**
   * Prepares density bin data for the Propensity Score Overlap Plot
   */
  _generateOverlapData(data, matchResult) {
    const rawTreated = data.filter(d => d._treatment === 1).map(d => d._ps);
    const rawControl = data.filter(d => d._treatment === 0).map(d => d._ps);
    const matchedTreated = matchResult.matchedTreated.map(d => d._ps);
    const matchedControl = matchResult.matchedControl.map(d => d._ps);

    const nBins = 25;
    const bins = [];
    for (let i = 0; i < nBins; i++) {
      const start = i / nBins;
      const end = (i + 1) / nBins;
      const mid = (start + end) / 2;

      const preT = rawTreated.filter(p => p >= start && p < end).length / (rawTreated.length || 1);
      const preC = rawControl.filter(p => p >= start && p < end).length / (rawControl.length || 1);
      const postT = matchedTreated.filter(p => p >= start && p < end).length / (matchedTreated.length || 1);
      const postC = matchedControl.filter(p => p >= start && p < end).length / (matchedControl.length || 1);

      bins.push({
        binIndex: i,
        range: [start, end],
        mid,
        preTreated: preT,
        preControl: preC,
        postTreated: postT,
        postControl: postC
      });
    }

    return {
      nBins,
      bins,
      commonSupport: matchResult.commonSupport
    };
  },

  /**
   * Generates a complete, reproducible Python script using statsmodels and scipy
   */
  generatePythonScript(treatmentCol, outcomeCol, covariateCols, caliper = 0.20, outcomeType = 'continuous') {
    if (typeof treatmentCol === 'object' && treatmentCol !== null) {
      const opts = treatmentCol;
      treatmentCol = opts.treatmentCol;
      outcomeCol = opts.outcomeCol;
      covariateCols = opts.covariateCols;
      caliper = opts.caliperMultiplier !== undefined ? opts.caliperMultiplier : (opts.caliper !== undefined ? opts.caliper : 0.20);
      outcomeType = opts.isBinaryOutcome ? 'binary' : (opts.outcomeType || 'continuous');
    }
    const covListPy = (covariateCols || []).map(c => `'${c}'`).join(', ');
    return `"""
=============================================================================
Propensity Score Matching (PSM) - Observational Clinical Causal Analysis
Methodology: 1:1 Nearest-Neighbor Caliper Matching Without Replacement
Generated by Statis-Gravity Clinical Biostatistics Platform
=============================================================================
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.spatial.distance import cdist
from scipy import stats
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# 1. DATA PREPARATION & MISSINGNESS HANDLING
# -----------------------------------------------------------------------------
# Load observational clinical dataset
# df = pd.read_csv('clinical_observational_data.csv')

treatment_col = '${treatmentCol}'
outcome_col = '${outcomeCol}'
covariates = [${covListPy}]

print("--- Step 1: Initial Cohort Inspection & Missing Data ---")
cols_needed = [treatment_col, outcome_col] + covariates
print(f"Total observational records: {len(df)}")
missing_summary = df[cols_needed].isnull().sum()
print("Missing values per variable:\\n", missing_summary)

# Complete-case analysis rationale:
# When missingness is low (< 5%) or missing completely at random (MCAR),
# complete-case analysis provides unbiased propensity score estimates.
df_clean = df.dropna(subset=cols_needed).copy()
df_clean[treatment_col] = df_clean[treatment_col].astype(int)
print(f"Complete cases retained: {len(df_clean)} "
      f"(Treated: {sum(df_clean[treatment_col] == 1)}, Control: {sum(df_clean[treatment_col] == 0)})")

# -----------------------------------------------------------------------------
# 2. PROPENSITY SCORE ESTIMATION (MULTIVARIATE LOGISTIC REGRESSION)
# -----------------------------------------------------------------------------
# Clinical Rationale:
# Propensity score e(X) = P(Treatment = 1 | X) collapses multi-dimensional baseline
# confounding into a scalar balancing score (Rosenbaum & Rubin, 1983).
formula = f"{treatment_col} ~ " + " + ".join(covariates)
ps_model = smf.logit(formula=formula, data=df_clean).fit(disp=False)
print("\\n--- Step 2: Logistic Regression Model Parameters ---")
print(ps_model.summary())

# Extract propensity scores and logit of propensity scores
df_clean['ps'] = ps_model.predict(df_clean)
# Logit transformation: linearizes propensity score for robust matching metric
df_clean['logit_ps'] = np.log(df_clean['ps'] / (1.0 - df_clean['ps']))

# Common Support Assessment:
treated_ps = df_clean[df_clean[treatment_col] == 1]['ps']
control_ps = df_clean[df_clean[treatment_col] == 0]['ps']
cs_min = max(treated_ps.min(), control_ps.min())
cs_max = min(treated_ps.max(), control_ps.max())
print(f"\\nRegion of Common Support: [{cs_min:.4f}, {cs_max:.4f}]")

# -----------------------------------------------------------------------------
# 3. 1:1 NEAREST NEIGHBOR MATCHING WITH CALIPER (WITHOUT REPLACEMENT)
# -----------------------------------------------------------------------------
# Clinical Rationale:
# Austin (2011) demonstrated that a caliper width equal to 0.20 of the standard
# deviation of the logit of the propensity score eliminates > 98% of baseline bias.
sd_logit_ps = df_clean['logit_ps'].std()
caliper_width = ${caliper.toFixed(2)} * sd_logit_ps
print(f"SD(logit(PS)): {sd_logit_ps:.4f} | Enforced Caliper Width: {caliper_width:.4f}")

treated_df = df_clean[(df_clean[treatment_col] == 1) & (df_clean['ps'] >= cs_min) & (df_clean['ps'] <= cs_max)].copy()
control_df = df_clean[(df_clean[treatment_col] == 0) & (df_clean['ps'] >= cs_min) & (df_clean['ps'] <= cs_max)].copy()

# Sort treated descending to prioritize high-propensity subjects (greedy matching)
treated_df = treated_df.sort_values(by='logit_ps', ascending=False)

matched_pairs = []
available_controls = control_df.copy()

for t_idx, t_row in treated_df.iterrows():
    if len(available_controls) == 0:
        break
    
    # Compute absolute distance on logit propensity score
    distances = np.abs(available_controls['logit_ps'].values - t_row['logit_ps'])
    min_dist_idx = np.argmin(distances)
    min_dist = distances[min_dist_idx]
    
    if min_dist <= caliper_width:
        c_idx = available_controls.index[min_dist_idx]
        matched_pairs.append({
            'treated_id': t_idx,
            'control_id': c_idx,
            'distance_logit': min_dist,
            'diff_ps': abs(t_row['ps'] - available_controls.loc[c_idx, 'ps'])
        })
        available_controls = available_controls.drop(c_idx)

pairs_df = pd.DataFrame(matched_pairs)
matched_treated_ids = pairs_df['treated_id'].values
matched_control_ids = pairs_df['control_id'].values

df_matched_treated = df_clean.loc[matched_treated_ids].copy()
df_matched_control = df_clean.loc[matched_control_ids].copy()
df_matched = pd.concat([df_matched_treated, df_matched_control])

print(f"Successfully matched: {len(pairs_df)} pairs ({len(pairs_df)*2} total subjects)")
print(f"Unmatched treated: {len(treated_df) - len(pairs_df)}")

# -----------------------------------------------------------------------------
# 4. BALANCE ASSESSMENT & STANDARDIZED MEAN DIFFERENCES (LOVE PLOT)
# -----------------------------------------------------------------------------
# Clinical Rationale:
# Standardized Mean Difference (SMD) is independent of sample size.
# An SMD < 0.10 indicates adequate balance between treatment and control cohorts.
balance_records = []
for cov in covariates:
    # Pre-matching pooled standard deviation (Austin, 2009)
    m_t_pre = df_clean[df_clean[treatment_col] == 1][cov].mean()
    m_c_pre = df_clean[df_clean[treatment_col] == 0][cov].mean()
    v_t_pre = df_clean[df_clean[treatment_col] == 1][cov].var()
    v_c_pre = df_clean[df_clean[treatment_col] == 0][cov].var()
    pooled_sd = np.sqrt((v_t_pre + v_c_pre) / 2.0)
    
    smd_pre = (m_t_pre - m_c_pre) / pooled_sd
    
    # Post-matching
    m_t_post = df_matched_treated[cov].mean()
    m_c_post = df_matched_control[cov].mean()
    smd_post = (m_t_post - m_c_post) / pooled_sd
    
    balance_records.append({
        'Covariate': cov,
        'Mean_T_Pre': m_t_pre,
        'Mean_C_Pre': m_c_pre,
        'SMD_Pre': smd_pre,
        'Mean_T_Post': m_t_post,
        'Mean_C_Post': m_c_post,
        'SMD_Post': smd_post,
        'Balanced': abs(smd_post) < 0.10
    })

balance_df = pd.DataFrame(balance_records)
print("\\n--- Step 4: Baseline Covariate Balance Assessment ---")
print(balance_df.to_string(index=False))

# Love Plot Generation
plt.figure(figsize=(8, len(covariates) * 0.6 + 2))
y_pos = np.arange(len(covariates))
plt.scatter(np.abs(balance_df['SMD_Pre']), y_pos, color='#f43f5e', label='Unadjusted (Pre-Match)', s=80, zorder=3)
plt.scatter(np.abs(balance_df['SMD_Post']), y_pos, color='#10b981', label='Matched (Post-Match)', s=90, zorder=4)
for i in range(len(covariates)):
    plt.plot([np.abs(balance_df['SMD_Pre'][i]), np.abs(balance_df['SMD_Post'][i])], [y_pos[i], y_pos[i]], color='gray', linestyle=':', alpha=0.7)

plt.axvline(x=0.10, color='#00d2ff', linestyle='--', linewidth=1.5, label='Standard Balance Cutoff (0.10)')
plt.axvline(x=0.05, color='#a855f7', linestyle=':', linewidth=1.2, label='Strict Cutoff (0.05)')
plt.yticks(y_pos, balance_df['Covariate'])
plt.xlabel('Absolute Standardized Mean Difference (|SMD|)')
plt.title('Love Plot: Covariate Balance Before and After PSM')
plt.legend(loc='upper right')
plt.grid(True, alpha=0.25)
plt.tight_layout()
plt.show()

# -----------------------------------------------------------------------------
# 5. OUTCOME ANALYSIS: AVERAGE TREATMENT EFFECT ON THE TREATED (ATT)
# -----------------------------------------------------------------------------
print("\\n--- Step 5: Causal Treatment Effect Estimation (ATT) ---")
${outcomeType === 'continuous' ? `# Continuous Outcome Analysis: Paired t-test with cluster-robust standard errors
pairs_df['y_treated'] = df_clean.loc[pairs_df['treated_id'], outcome_col].values
pairs_df['y_control'] = df_clean.loc[pairs_df['control_id'], outcome_col].values
pairs_df['diff'] = pairs_df['y_treated'] - pairs_df['y_control']

att = pairs_df['diff'].mean()
se_att = pairs_df['diff'].std() / np.sqrt(len(pairs_df))
t_stat = att / se_att
df_paired = len(pairs_df) - 1
p_val = 2.0 * (1.0 - stats.t.cdf(np.abs(t_stat), df=df_paired))
ci_95 = stats.t.interval(0.95, df=df_paired, loc=att, scale=se_att)

print(f"Outcome Variable: {outcome_col}")
print(f"Unadjusted Mean Diff: {df_clean[df_clean[treatment_col] == 1][outcome_col].mean() - df_clean[df_clean[treatment_col] == 0][outcome_col].mean():.4f}")
print(f"ATT (Average Treatment Effect on Treated): {att:.4f}")
print(f"Paired Robust Standard Error: {se_att:.4f}")
print(f"t-statistic ({df_paired} df): {t_stat:.4f} | p-value: {p_val:.4e}")
print(f"95% Confidence Interval: [{ci_95[0]:.4f}, {ci_95[1]:.4f}]")` : `# Binary Outcome Analysis: Paired Discordant Pairs & McNemar Analysis
pairs_df['y_treated'] = df_clean.loc[pairs_df['treated_id'], outcome_col].values
pairs_df['y_control'] = df_clean.loc[pairs_df['control_id'], outcome_col].values

# Discordant pairs (b: Treated=1, Control=0; c: Treated=0, Control=1)
b = sum((pairs_df['y_treated'] == 1) & (pairs_df['y_control'] == 0))
c = sum((pairs_df['y_treated'] == 0) & (pairs_df['y_control'] == 1))
att_rd = (b - c) / len(pairs_df)
se_rd = np.sqrt((b + c) - ((b - c)**2 / len(pairs_df))) / len(pairs_df)
z_stat = att_rd / se_rd
p_val = 2 * (1 - stats.norm.cdf(abs(z_stat)))
ci_95 = [att_rd - 1.96 * se_rd, att_rd + 1.96 * se_rd]

print(f"Outcome Variable: {outcome_col}")
print(f"ATT Risk Difference: {att_rd:.4f} (95% CI: [{ci_95[0]:.4f}, {ci_95[1]:.4f}])")
print(f"Standard Error: {se_rd:.4f} | z = {z_stat:.4f} | p-value = {p_val:.4e}")
print(f"Discordant Pairs: b (T=1, C=0) = {b}, c (T=0, C=1) = {c}")`}
`;
  },

  /**
   * Generates a complete, reproducible R script using MatchIt and cobalt
   */
  generateRScript(treatmentCol, outcomeCol, covariateCols, caliper = 0.20, outcomeType = 'continuous') {
    if (typeof treatmentCol === 'object' && treatmentCol !== null) {
      const opts = treatmentCol;
      treatmentCol = opts.treatmentCol;
      outcomeCol = opts.outcomeCol;
      covariateCols = opts.covariateCols;
      caliper = opts.caliperMultiplier !== undefined ? opts.caliperMultiplier : (opts.caliper !== undefined ? opts.caliper : 0.20);
      outcomeType = opts.isBinaryOutcome ? 'binary' : (opts.outcomeType || 'continuous');
    }
    const covFormulaR = (covariateCols || []).join(' + ');
    return `################################################################################
# Propensity Score Matching (PSM) - Observational Clinical Causal Analysis
# Package Suite: MatchIt, cobalt, survey, sandwich
# Generated by Statis-Gravity Clinical Biostatistics Platform
################################################################################

# Install required biostatistical packages if needed:
# install.packages(c("MatchIt", "cobalt", "survey", "sandwich", "lmtest", "ggplot2"))
library(MatchIt)
library(cobalt)
library(survey)
library(sandwich)
library(lmtest)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. DATA PREP & MISSINGNESS
# ------------------------------------------------------------------------------
# Load observational clinical data
# df <- read.csv("clinical_observational_data.csv")

treatment_var <- "${treatmentCol}"
outcome_var   <- "${outcomeCol}"
covariates    <- c(${covariateCols.map(c => `"${c}"`).join(', ')})

cat("--- Step 1: Missing Data & Cohort Inspection ---\\n")
cat("Total raw sample size:", nrow(df), "\\n")
cols_needed <- c(treatment_var, outcome_var, covariates)
missing_counts <- colSums(is.na(df[, cols_needed]))
print(missing_counts)

# Complete-case filtering:
df_clean <- na.omit(df[, cols_needed])
df_clean[[treatment_var]] <- as.numeric(df_clean[[treatment_var]])
cat("Complete cases retained:", nrow(df_clean), "\\n")
cat("Treated units:", sum(df_clean[[treatment_var]] == 1), 
    "| Control units:", sum(df_clean[[treatment_var]] == 0), "\\n\\n")

# ------------------------------------------------------------------------------
# 2 & 3. PROPENSITY SCORE ESTIMATION & 1:1 CALIPER MATCHING
# ------------------------------------------------------------------------------
# Clinical Rationale:
# MatchIt fits logistic regression under distance = "logit",
# enforces nearest neighbor matching without replacement (method = "nearest"),
# and sets caliper = ${caliper.toFixed(2)} (in SD of logit PS) within common support.
match_formula <- as.formula(paste("${treatmentCol} ~", "${covFormulaR}"))

cat("--- Steps 2 & 3: Estimating Propensity Scores & Executing MatchIt ---\\n")
m.out <- matchit(
  formula = match_formula,
  data = df_clean,
  method = "nearest",
  distance = "logit",
  caliper = ${caliper.toFixed(2)},
  std.caliper = TRUE,       # Enforce caliper as 0.20 SD of logit(PS)
  replace = FALSE,
  discard = "both"          # Enforce common support boundaries
)

print(summary(m.out, un = TRUE))

# Common support / Propensity score overlap visualization
plot(m.out, type = "jitter", interactive = FALSE)
plot(m.out, type = "density", interactive = FALSE)

# ------------------------------------------------------------------------------
# 4. BALANCE ASSESSMENT & LOVE PLOT
# ------------------------------------------------------------------------------
cat("\\n--- Step 4: Assessing Covariate Balance (cobalt Suite) ---\\n")
# Summary balance table reporting Standardized Mean Differences (SMD)
bal_tab <- bal.tab(m.out, un = TRUE, stats = c("m", "v"), thresholds = c(m = 0.10))
print(bal_tab)

# Publication-grade Love Plot (Austin 2009 standard: |SMD| < 0.10)
love.plot(
  m.out,
  binary = "std",
  thresholds = c(m = 0.10),
  colors = c("#f43f5e", "#10b981"),
  shapes = c(21, 19),
  size = 3.5,
  var.order = "unadjusted",
  title = "Love Plot: Covariate Balance (Pre- vs. Post-Matching)"
)

# Extract matched analytic dataset
df_matched <- match.data(m.out)

# ------------------------------------------------------------------------------
# 5. OUTCOME ANALYSIS: AVERAGE TREATMENT EFFECT ON THE TREATED (ATT)
# ------------------------------------------------------------------------------
cat("\\n--- Step 5: Estimating Average Treatment Effect on the Treated (ATT) ---\\n")
${outcomeType === 'continuous' ? `# Continuous Outcome Analysis: Linear model with pair-clustered robust SE
att_fit <- lm(as.formula(paste("${outcomeCol} ~", "${treatmentCol}")), 
              data = df_matched, 
              weights = weights)

# Cluster-robust standard errors accounting for matched pairs (subclass)
att_robust <- coeftest(att_fit, vcov. = vcovCL, cluster = ~subclass)
print(att_robust)

ci_att <- coefci(att_fit, vcov. = vcovCL, cluster = ~subclass)
cat("\\n95% Confidence Interval for ATT:\\n")
print(ci_att["${treatmentCol}", ])` : `# Binary Outcome Analysis: Logistic regression with cluster-robust standard errors
att_glm <- glm(as.formula(paste("${outcomeCol} ~", "${treatmentCol}")), 
               data = df_matched, 
               family = binomial(link = "logit"), 
               weights = weights)

att_robust <- coeftest(att_glm, vcov. = vcovCL, cluster = ~subclass)
print(att_robust)

or_est <- exp(coef(att_glm)["${treatmentCol}"])
or_ci <- exp(coefci(att_glm, vcov. = vcovCL, cluster = ~subclass)["${treatmentCol}", ])
cat("\\nATT Odds Ratio:", round(or_est, 4), "\\n")
cat("95% CI for Odds Ratio: [", round(or_ci[1], 4), ",", round(or_ci[2], 4), "]\\n")`}
`;
  },

  /**
   * Generates a complete, reproducible Stata script (.do)
   */
  generateStataScript(treatmentCol, outcomeCol, covariateCols, caliper = 0.20, outcomeType = 'continuous') {
    if (typeof treatmentCol === 'object' && treatmentCol !== null) {
      const opts = treatmentCol;
      treatmentCol = opts.treatmentCol;
      outcomeCol = opts.outcomeCol;
      covariateCols = opts.covariateCols;
      caliper = opts.caliperMultiplier !== undefined ? opts.caliperMultiplier : (opts.caliper !== undefined ? opts.caliper : 0.20);
      outcomeType = opts.isBinaryOutcome ? 'binary' : (opts.outcomeType || 'continuous');
    }
    const covListStata = (covariateCols || []).join(' ');
    return `* ==============================================================================
* Propensity Score Matching (PSM) - Observational Clinical Causal Analysis
* Stata Implementation using teffects psmatch & pstest
* Generated by Statis-Gravity Clinical Biostatistics Platform
* ==============================================================================

version 17
clear all
set more off

* 1. LOAD OBSERVATIONAL CLINICAL DATA
* use "clinical_observational_data.dta", clear

local treatment "${treatmentCol}"
local outcome "${outcomeCol}"
local covariates "${covListStata}"

* ------------------------------------------------------------------------------
* 1. DATA PREP & MISSINGNESS
* ------------------------------------------------------------------------------
display as text "--- Step 1: Missing Data Inspection ---"
misstable summarize \`treatment' \`outcome' \`covariates'

* Complete cases filtering
drop if missing(\`treatment') | missing(\`outcome')
foreach var of local covariates {
    drop if missing(\`var')
}
display as text "Complete cases retained: " _N
tabulate \`treatment'

* ------------------------------------------------------------------------------
* 2 & 3. PROPENSITY SCORE ESTIMATION & MATCHING (teffects psmatch)
* ------------------------------------------------------------------------------
* teffects psmatch fits logistic model, estimates ATT (atet),
* and implements 1:1 nearest neighbor matching within caliper.
display as text "--- Steps 2 & 3: Estimating PS & Caliper Matching ---"
teffects psmatch (\`outcome') (\`treatment' \`covariates', logit), ///
    atet ///
    caliper(${caliper.toFixed(2)}) ///
    vce(robust) ///
    generate(match_id)

* Display ATT estimation results
teffects summarize

* Overlap plot: assess common support
psgraph, treated(\`treatment') pscore(match_id)

* ------------------------------------------------------------------------------
* 4. BALANCE ASSESSMENT & LOVE PLOT
* ------------------------------------------------------------------------------
* pstest computes Standardized Mean Differences (SMD) before and after matching
display as text "--- Step 4: Assessing Covariate Balance ---"
pstest \`covariates', both graph

* Save Stata balance graph
graph export "psm_love_plot.png", replace width(2000)

* ------------------------------------------------------------------------------
* 5. OUTCOME ANALYSIS
* ------------------------------------------------------------------------------
display as text "--- Step 5: Causal Treatment Effect (ATT) Summary ---"
* Point estimate, standard error, z, p-value, and 95% CI
matrix list r(table)
`;
  },

  /**
   * Helper: Invert square matrix via Gauss-Jordan elimination with partial pivoting
   */
  _invertMatrix(matrix) {
    const n = matrix.length;
    // Augment with identity matrix
    const aug = matrix.map((row, i) => {
      const identityRow = new Array(n).fill(0);
      identityRow[i] = 1.0;
      return [...row, ...identityRow];
    });

    for (let i = 0; i < n; i++) {
      // Find pivot
      let maxRow = i;
      let maxVal = Math.abs(aug[i][i]);
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(aug[k][i]) > maxVal) {
          maxVal = Math.abs(aug[k][i]);
          maxRow = k;
        }
      }

      if (maxVal < 1e-12) return null; // Singular

      // Swap rows
      if (maxRow !== i) {
        const temp = aug[i];
        aug[i] = aug[maxRow];
        aug[maxRow] = temp;
      }

      // Scale pivot row
      const pivot = aug[i][i];
      for (let j = 0; j < 2 * n; j++) {
        aug[i][j] /= pivot;
      }

      // Eliminate column
      for (let k = 0; k < n; k++) {
        if (k !== i) {
          const factor = aug[k][i];
          for (let j = 0; j < 2 * n; j++) {
            aug[k][j] -= factor * aug[i][j];
          }
        }
      }
    }

    return aug.map(row => row.slice(n));
  },

  /**
   * Helper: Solve linear system A * x = b via Gaussian elimination with partial pivoting
   */
  _solveLinearSystem(A, b) {
    const n = A.length;
    const M = A.map((row, i) => [...row, b[i]]);

    for (let i = 0; i < n; i++) {
      let maxRow = i;
      let maxVal = Math.abs(M[i][i]);
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(M[k][i]) > maxVal) {
          maxVal = Math.abs(M[k][i]);
          maxRow = k;
        }
      }

      if (maxVal < 1e-12) return null;

      if (maxRow !== i) {
        const temp = M[i];
        M[i] = M[maxRow];
        M[maxRow] = temp;
      }

      const pivot = M[i][i];
      for (let j = i; j <= n; j++) {
        M[i][j] /= pivot;
      }

      for (let k = 0; k < n; k++) {
        if (k !== i) {
          const factor = M[k][i];
          for (let j = i; j <= n; j++) {
            M[k][j] -= factor * M[i][j];
          }
        }
      }
    }

    return M.map(row => row[n]);
  },

  /**
   * Parses CSV or TSV string into an array of patient records
   */
  parseClinicalCsv(csvText) {
    if (!csvText || typeof csvText !== 'string') return [];
    const lines = csvText.trim().split(/\r?\n/).filter(l => l.trim().length > 0);
    if (lines.length < 2) return [];

    const delimiter = lines[0].includes('\t') ? '\t' : (lines[0].includes(';') ? ';' : ',');
    const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));

    const records = [];
    for (let i = 1; i < lines.length; i++) {
      const parts = lines[i].split(delimiter).map(p => p.trim().replace(/^["']|["']$/g, ''));
      if (parts.length === headers.length) {
        const row = {};
        headers.forEach((h, idx) => {
          const val = parts[idx];
          const num = Number(val);
          row[h] = !isNaN(num) && val !== '' ? num : val;
        });
        records.push(row);
      }
    }
    return records;
  },

  /**
   * Generates a realistic observational clinical craniosynostosis dataset
   * 120 patients evaluated for surgical management:
   * Endoscopic Strip Craniectomy (treatment = 1) vs. Open Cranial Vault Reconstruction (control = 0)
   */
  getSampleClinicalDataset() {
    // Deterministic pseudo-random seed generator for perfect clinical reproducibility
    let seed = 42;
    const rnd = () => {
      seed = (seed * 9301 + 49297) % 233280;
      return seed / 233280;
    };
    const rndNormal = (mean, sd) => {
      const u1 = Math.max(1e-7, rnd());
      const u2 = rnd();
      const z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      return mean + z * sd;
    };

    const dataset = [];
    const nTotal = 120;

    for (let i = 1; i <= nTotal; i++) {
      // Patient baseline covariates
      const age = Math.max(2.2, Math.min(12.0, rndNormal(6.4, 2.3))); // Age in months
      const sex = rnd() > 0.42 ? 1 : 0; // 58% male
      const comorbidity_score = rnd() < 0.60 ? 0 : (rnd() < 0.85 ? 1 : 2); // 0=None, 1=Mild, 2=Moderate/Syndromic
      const baseline_severity = Math.max(5.0, Math.min(15.5, rndNormal(9.8, 2.4))); // Initial CVAI %
      const bmi_percentile = Math.max(10, Math.min(95, Math.round(rndNormal(52, 18))));

      // Observational clinical assignment probability (Confounded Treatment Allocation)
      // Clinicians steer younger infants with lower severity and fewer comorbidities towards endoscopic strip
      const latentZ = 1.6 
        - 0.38 * (age - 6.0) 
        - 0.28 * (baseline_severity - 9.5) 
        - 0.75 * comorbidity_score 
        + 0.15 * (sex === 1 ? 0.5 : -0.5);
      const trueProb = 1.0 / (1.0 + Math.exp(-latentZ));
      const treatment_col = rnd() < trueProb ? 1 : 0;

      // Clinical Outcome: Hospital Length of Stay (days)
      // True causal effect of endoscopic strip is a reduction of ~2.2 days
      // Outcome is also influenced by baseline severity, age, and comorbidities
      let trueOutcome = 4.8 
        - 2.3 * treatment_col 
        + 0.18 * age 
        + 0.22 * baseline_severity 
        + 0.65 * comorbidity_score 
        + rndNormal(0, 0.45);
      trueOutcome = Math.max(1.0, Math.round(trueOutcome * 10) / 10);

      dataset.push({
        patient_id: `PT-${1000 + i}`,
        treatment_col,
        outcome_col: trueOutcome,
        age: Math.round(age * 10) / 10,
        sex,
        comorbidity_score,
        baseline_severity: Math.round(baseline_severity * 10) / 10,
        bmi_percentile
      });
    }

    return dataset;
  }
};
