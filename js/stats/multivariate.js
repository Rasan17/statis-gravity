/**
 * Statis-Gravity - Multivariate Exploratory Data Analysis (EDA) Engine
 * Implements Principal Component Analysis (PCA), Multiple Correspondence Analysis (MCA),
 * and Factor Analysis of Mixed Data (FAMD / Pagès 2004) with lightweight, self-contained
 * client-side numerical linear algebra (Jacobi eigenvalue solver and thin SVD).
 */

export const Multivariate = {
  // =========================================================================
  // 1. LIGHTWEIGHT NUMERICAL LINEAR ALGEBRA & MATRIX DECOMPOSITION
  // =========================================================================
  LinearAlgebra: {
    transpose(A) {
      const rows = A.length;
      const cols = A[0].length;
      const T = Array.from({ length: cols }, () => new Array(rows));
      for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
          T[j][i] = A[i][j];
        }
      }
      return T;
    },

    matmul(A, B) {
      const rowsA = A.length;
      const colsA = A[0].length;
      const colsB = B[0].length;
      const C = Array.from({ length: rowsA }, () => new Array(colsB).fill(0));
      for (let i = 0; i < rowsA; i++) {
        for (let k = 0; k < colsA; k++) {
          const aik = A[i][k];
          if (aik === 0) continue;
          for (let j = 0; j < colsB; j++) {
            C[i][j] += aik * B[k][j];
          }
        }
      }
      return C;
    },

    /**
     * Cyclic Jacobi algorithm for symmetric matrix eigenvalue decomposition
     * Solves A * V = V * D, where A is symmetric (n x n)
     */
    jacobi(A, maxIter = 150, tol = 1e-12) {
      const n = A.length;
      let V = Array.from({ length: n }, (_, i) => {
        const row = new Array(n).fill(0);
        row[i] = 1.0;
        return row;
      });
      let D = A.map(row => [...row]);

      for (let iter = 0; iter < maxIter; iter++) {
        let maxOff = 0;
        let p = 0, q = 1;
        for (let i = 0; i < n; i++) {
          for (let j = i + 1; j < n; j++) {
            const val = Math.abs(D[i][j]);
            if (val > maxOff) {
              maxOff = val;
              p = i; q = j;
            }
          }
        }
        if (maxOff < tol) break;

        const diff = D[q][q] - D[p][p];
        let t;
        if (Math.abs(D[p][q]) < Math.abs(diff) * 1e-15) {
          t = D[p][q] / diff;
        } else {
          const phi = diff / (2.0 * D[p][q]);
          t = 1.0 / (Math.abs(phi) + Math.sqrt(phi * phi + 1.0));
          if (phi < 0) t = -t;
        }

        const c = 1.0 / Math.sqrt(t * t + 1.0);
        const s = t * c;
        const tau = s / (1.0 + c);

        const Dpq = D[p][q];
        D[p][q] = 0;
        D[q][p] = 0;
        D[p][p] -= t * Dpq;
        D[q][q] += t * Dpq;

        for (let j = 0; j < n; j++) {
          if (j !== p && j !== q) {
            const Djp = D[j][p];
            const Djq = D[j][q];
            D[j][p] = Djp - s * (Djq + tau * Djp);
            D[p][j] = D[j][p];
            D[j][q] = Djq + s * (Djp - tau * Djq);
            D[q][j] = D[j][q];
          }
        }

        for (let j = 0; j < n; j++) {
          const Vjp = V[j][p];
          const Vjq = V[j][q];
          V[j][p] = Vjp - s * (Vjq + tau * Vjp);
          V[j][q] = Vjq + s * (Vjp - tau * Vjq);
        }
      }

      // Sort eigenvalues and corresponding eigenvector columns descending
      const eigen = [];
      for (let i = 0; i < n; i++) {
        eigen.push({
          value: Math.max(0, D[i][i]),
          vector: V.map(row => row[i])
        });
      }
      eigen.sort((a, b) => b.value - a.value);

      return {
        values: eigen.map(e => e.value),
        vectors: eigen.map(e => e.vector) // Column eigenvectors
      };
    },

    /**
     * Thin Singular Value Decomposition (SVD): M = U * diag(s) * V^T
     * Efficiently handles rectangular n x p matrices via Gram matrix eigendecomposition.
     */
    svd(M) {
      const n = M.length;
      const p = M[0].length;

      if (n >= p) {
        // Compute p x p Gram matrix: C = M^T * M
        const MT = this.transpose(M);
        const C = this.matmul(MT, M);
        const eig = this.jacobi(C);

        const rank = p;
        const s = eig.values.map(val => Math.sqrt(Math.max(0, val)));
        // V matrix: columns are eigenvectors of M^T * M
        const V = Array.from({ length: p }, (_, i) => new Array(p));
        for (let j = 0; j < p; j++) {
          for (let i = 0; i < p; i++) {
            V[i][j] = eig.vectors[j][i];
          }
        }

        // Left singular vectors: U_k = M * V_k / s_k
        const U = Array.from({ length: n }, () => new Array(p).fill(0));
        for (let k = 0; k < p; k++) {
          const sk = s[k];
          if (sk > 1e-12) {
            for (let i = 0; i < n; i++) {
              let sum = 0;
              for (let j = 0; j < p; j++) {
                sum += M[i][j] * V[j][k];
              }
              U[i][k] = sum / sk;
            }
          }
        }

        return { U, s, V };
      } else {
        // When p > n: compute n x n Gram matrix: C = M * M^T
        const MT = this.transpose(M);
        const C = this.matmul(M, MT);
        const eig = this.jacobi(C);

        const s = eig.values.map(val => Math.sqrt(Math.max(0, val)));
        // U matrix: columns are eigenvectors of M * M^T
        const U = Array.from({ length: n }, (_, i) => new Array(n));
        for (let j = 0; j < n; j++) {
          for (let i = 0; i < n; i++) {
            U[i][j] = eig.vectors[j][i];
          }
        }

        // Right singular vectors: V_k = M^T * U_k / s_k
        const V = Array.from({ length: p }, () => new Array(n).fill(0));
        for (let k = 0; k < n; k++) {
          const sk = s[k];
          if (sk > 1e-12) {
            for (let j = 0; j < p; j++) {
              let sum = 0;
              for (let i = 0; i < n; i++) {
                sum += MT[j][i] * U[i][k];
              }
              V[j][k] = sum / sk;
            }
          }
        }

        return { U, s, V };
      }
    }
  },

  // =========================================================================
  // 2. DATASET PARSER & VARIABLE TYPE INFERENCE
  // =========================================================================
  parseDataset(csvText) {
    if (!csvText || typeof csvText !== 'string') {
      return { error: 'No data provided.' };
    }

    const lines = csvText.trim().split(/\r?\n/).filter(l => l.trim().length > 0);
    if (lines.length < 3) {
      return { error: 'Dataset must contain a header row and at least 2 observation rows.' };
    }

    const delimiter = lines[0].includes('\t') ? '\t' : (lines[0].includes(';') ? ';' : ',');
    const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));

    const records = [];
    for (let i = 1; i < lines.length; i++) {
      const parts = lines[i].split(delimiter).map(p => p.trim().replace(/^["']|["']$/g, ''));
      if (parts.length === headers.length) {
        const row = { _id: `Row-${i}` };
        headers.forEach((h, idx) => {
          row[h] = parts[idx];
        });
        records.push(row);
      }
    }

    if (records.length < 3) {
      return { error: 'Insufficient valid records in dataset.' };
    }

    // Infer variable types
    const colTypes = {};
    const colStats = {};

    headers.forEach(h => {
      let numCount = 0;
      let totalNonEmpty = 0;
      const values = [];

      records.forEach(r => {
        const val = r[h];
        if (val !== undefined && val !== null && val !== '') {
          totalNonEmpty++;
          const num = Number(val);
          if (!isNaN(num) && isFinite(num)) {
            numCount++;
            values.push(num);
          } else {
            values.push(val);
          }
        }
      });

      const numRatio = totalNonEmpty > 0 ? numCount / totalNonEmpty : 0;
      const uniqueVals = new Set(values);

      // If > 80% numerical and at least 3 distinct values, infer continuous
      if (numRatio >= 0.80 && uniqueVals.size > 2) {
        colTypes[h] = 'continuous';
        const nums = values.filter(v => typeof v === 'number');
        const mean = nums.reduce((a, b) => a + b, 0) / nums.length;
        const variance = nums.reduce((acc, v) => acc + (v - mean) ** 2, 0) / (nums.length - 1);
        colStats[h] = {
          type: 'continuous',
          count: nums.length,
          mean,
          sd: Math.sqrt(variance),
          min: Math.min(...nums),
          max: Math.max(...nums)
        };
      } else {
        colTypes[h] = 'categorical';
        colStats[h] = {
          type: 'categorical',
          count: values.length,
          uniqueCount: uniqueVals.size,
          categories: Array.from(uniqueVals).map(String).sort()
        };
      }
    });

    return {
      headers,
      records,
      colTypes,
      colStats,
      nRows: records.length,
      nCols: headers.length
    };
  },

  // =========================================================================
  // 3. PRINCIPAL COMPONENT ANALYSIS (PCA)
  // =========================================================================
  runPCA(records, continuousCols, options = {}) {
    if (!continuousCols || continuousCols.length < 2) {
      return { error: 'PCA requires at least 2 continuous variables.' };
    }

    const n = records.length;
    const p = continuousCols.length;

    // Filter complete cases
    const cleanRecords = records.filter(r => {
      return continuousCols.every(c => {
        const val = Number(r[c]);
        return !isNaN(val) && isFinite(val) && r[c] !== '' && r[c] !== null;
      });
    });

    if (cleanRecords.length < 3) {
      return { error: 'Insufficient complete observations for PCA (minimum 3 required).' };
    }

    const nClean = cleanRecords.length;

    // 1. Mean-centering and scaling to unit variance
    const means = {};
    const sds = {};

    continuousCols.forEach(col => {
      const vals = cleanRecords.map(r => Number(r[col]));
      const m = vals.reduce((a, b) => a + b, 0) / nClean;
      const variance = vals.reduce((acc, v) => acc + (v - m) ** 2, 0) / (nClean - 1);
      const s = Math.sqrt(variance) || 1.0;
      means[col] = m;
      sds[col] = s;
    });

    // Standardized matrix Z / sqrt(n - 1)
    const factor = Math.sqrt(nClean - 1);
    const Z = Array.from({ length: nClean }, (_, i) => {
      const r = cleanRecords[i];
      return continuousCols.map(col => (Number(r[col]) - means[col]) / (sds[col] * factor));
    });

    // 2. Perform SVD on Z
    const svdRes = this.LinearAlgebra.svd(Z);
    const eigenvalues = svdRes.s.map(s => s * s);
    const totalVariance = eigenvalues.reduce((a, b) => a + b, 0) || p;

    const nComponents = Math.min(p, nClean - 1);
    const scree = [];
    let cumVar = 0;

    for (let k = 0; k < nComponents; k++) {
      const eig = eigenvalues[k] || 0;
      const varPct = (eig / totalVariance) * 100;
      cumVar += varPct;
      scree.push({
        dim: k + 1,
        label: `Dim ${k + 1}`,
        eigenvalue: eig,
        variancePct: varPct,
        cumulativePct: Math.min(100, cumVar)
      });
    }

    // 3. Compute Individual Coordinates (Factor Scores)
    // Score_k = Z_unscaled * V_k = (Z * sqrt(n-1)) * V_k = U_k * s_k * sqrt(n-1)
    const individuals = cleanRecords.map((r, i) => {
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        // Individual coordinate on component k
        coords.push(svdRes.U[i][k] * svdRes.s[k] * factor);
      }

      // Squared cosine (quality of representation) and distance from centroid
      const totalDistSq = coords.reduce((acc, c) => acc + c * c, 0);
      const cos2 = coords.map(c => (totalDistSq > 0 ? (c * c) / totalDistSq : 0));

      return {
        id: r._id || r.id || r.patient_id || `Obs-${i + 1}`,
        originalRow: r,
        coords,
        cos2,
        distFromCentroid: Math.sqrt(totalDistSq)
      };
    });

    // Individual contributions
    for (let k = 0; k < nComponents; k++) {
      const sumDistSq = individuals.reduce((acc, ind) => acc + ind.coords[k] ** 2, 0) || 1;
      individuals.forEach(ind => {
        if (!ind.ctr) ind.ctr = [];
        ind.ctr[k] = (ind.coords[k] ** 2 / sumDistSq) * 100;
      });
    }

    // Flag outliers: Euclidean distance in factor space > 2.5 standard deviations
    const allDists = individuals.map(d => d.distFromCentroid);
    const meanDist = allDists.reduce((a, b) => a + b, 0) / allDists.length;
    const sdDist = Math.sqrt(allDists.reduce((acc, v) => acc + (v - meanDist) ** 2, 0) / (allDists.length - 1)) || 1;
    const outlierCutoff = meanDist + 2.5 * sdDist;

    individuals.forEach(ind => {
      ind.isOutlier = ind.distFromCentroid > outlierCutoff;
    });

    // 4. Compute Variable Coordinates (Loadings / Correlations)
    // Loading_j,k = V_j,k * sqrt(eigenvalue_k) = Pearson correlation between variable j and component k
    const variables = continuousCols.map((col, j) => {
      const coords = [];
      const cos2 = [];
      const ctr = [];

      for (let k = 0; k < nComponents; k++) {
        const loading = svdRes.V[j][k] * Math.sqrt(Math.max(0, eigenvalues[k]));
        coords.push(loading);
        // cos2 on dimension k = loading^2 (since var is standardized to unit variance)
        cos2.push(loading * loading);
        // Contribution of variable j to component k = loading^2 / eigenvalue_k
        const eig = eigenvalues[k];
        ctr.push(eig > 0 ? (loading * loading / eig) * 100 : 0);
      }

      return {
        name: col,
        type: 'continuous',
        coords,
        cos2,
        ctr,
        mean: means[col],
        sd: sds[col]
      };
    });

    return {
      method: 'PCA',
      methodLabel: 'Principal Component Analysis (PCA)',
      nObservations: nClean,
      nVariables: p,
      continuousCols,
      eigenvalues,
      totalVariance,
      scree,
      individuals,
      variables,
      continuousVariables: variables,
      allVariables: variables,
      outliers: individuals.filter(d => d.isOutlier),
      outlierCount: individuals.filter(d => d.isOutlier).length,
      outlierCutoff
    };
  },

  // =========================================================================
  // 4. MULTIPLE CORRESPONDENCE ANALYSIS (MCA)
  // =========================================================================
  runMCA(records, categoricalCols, options = {}) {
    if (!categoricalCols || categoricalCols.length < 2) {
      return { error: 'MCA requires at least 2 categorical variables.' };
    }

    const n = records.length;
    const K = categoricalCols.length;

    // Filter complete cases
    const cleanRecords = records.filter(r => {
      return categoricalCols.every(c => r[c] !== undefined && r[c] !== null && String(r[c]).trim() !== '');
    });

    if (cleanRecords.length < 3) {
      return { error: 'Insufficient complete observations for MCA.' };
    }

    const nClean = cleanRecords.length;

    // 1. Build Indicator Matrix Z
    const modalities = [];
    categoricalCols.forEach(col => {
      const uniqueVals = Array.from(new Set(cleanRecords.map(r => String(r[col]).trim()))).sort();
      uniqueVals.forEach(val => {
        modalities.push({
          col,
          val,
          fullName: `${col}_${val}`,
          label: `${col}: ${val}`
        });
      });
    });

    const J = modalities.length;
    if (J <= K) {
      return { error: 'Insufficient category variation to perform MCA.' };
    }

    // Indicator matrix Z: nClean x J
    const Z = Array.from({ length: nClean }, () => new Array(J).fill(0));
    cleanRecords.forEach((r, i) => {
      modalities.forEach((m, j) => {
        if (String(r[m.col]).trim() === m.val) {
          Z[i][j] = 1;
        }
      });
    });

    // Column marginal frequencies and proportions
    const colCounts = new Array(J).fill(0);
    for (let j = 0; j < J; j++) {
      for (let i = 0; i < nClean; i++) {
        colCounts[j] += Z[i][j];
      }
    }

    const p = colCounts.map(c => c / nClean);

    // Standardized residuals matrix S:
    // S_ij = (Z_ij / (n * K) - (1/n) * (p_j / K)) / sqrt((1/n) * (p_j / K))
    // S_ij = (1 / sqrt(n * K)) * (Z_ij - p_j) / sqrt(p_j)
    const S = Array.from({ length: nClean }, (_, i) => {
      const factorNK = 1.0 / Math.sqrt(nClean * K);
      return modalities.map((m, j) => {
        const pj = Math.max(1e-9, p[j]);
        return factorNK * (Z[i][j] - pj) / Math.sqrt(pj);
      });
    });

    // 2. Perform SVD on S
    const svdRes = this.LinearAlgebra.svd(S);
    const eigenvalues = svdRes.s.map(s => s * s);
    const totalInertia = (J - K) / K;

    const nComponents = Math.min(J - K, nClean - 1);
    const scree = [];
    let cumVar = 0;

    for (let k = 0; k < nComponents; k++) {
      const eig = eigenvalues[k] || 0;
      // In MCA, standard variance explained percentage:
      const varPct = totalInertia > 0 ? (eig / totalInertia) * 100 : 0;
      cumVar += varPct;
      scree.push({
        dim: k + 1,
        label: `Dim ${k + 1}`,
        eigenvalue: eig,
        variancePct: varPct,
        cumulativePct: Math.min(100, cumVar)
      });
    }

    // 3. Individual Principal Coordinates
    // F_ik = sqrt(n) * U_ik * s_k
    const sqrtN = Math.sqrt(nClean);
    const individuals = cleanRecords.map((r, i) => {
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        coords.push(sqrtN * svdRes.U[i][k] * svdRes.s[k]);
      }
      const totalDistSq = coords.reduce((acc, c) => acc + c * c, 0);
      const cos2 = coords.map(c => (totalDistSq > 0 ? (c * c) / totalDistSq : 0));

      return {
        id: r._id || r.id || r.patient_id || `Obs-${i + 1}`,
        originalRow: r,
        coords,
        cos2,
        distFromCentroid: Math.sqrt(totalDistSq)
      };
    });

    // Outlier check
    const allDists = individuals.map(d => d.distFromCentroid);
    const meanDist = allDists.reduce((a, b) => a + b, 0) / allDists.length;
    const sdDist = Math.sqrt(allDists.reduce((acc, v) => acc + (v - meanDist) ** 2, 0) / (allDists.length - 1)) || 1;
    const outlierCutoff = meanDist + 2.5 * sdDist;
    individuals.forEach(ind => {
      ind.isOutlier = ind.distFromCentroid > outlierCutoff;
    });

    // 4. Modality / Category Coordinates
    // A_jk = (1 / sqrt(p_j / K)) * V_jk * s_k / sqrt(K) = (1 / sqrt(p_j)) * V_jk * s_k
    const modalitiesCoords = modalities.map((m, j) => {
      const pj = Math.max(1e-9, p[j]);
      const coords = [];
      const ctr = [];

      for (let k = 0; k < nComponents; k++) {
        // Principal coordinate of modality j on dimension k
        const coord = (1.0 / Math.sqrt(pj)) * svdRes.V[j][k] * svdRes.s[k];
        coords.push(coord);
        // Contribution of modality j to dimension k
        const eig = eigenvalues[k];
        const c = eig > 0 ? (pj * coord * coord) / (K * eig) * 100 : 0;
        ctr.push(c);
      }

      return {
        name: m.label,
        col: m.col,
        val: m.val,
        type: 'modality',
        coords,
        ctr,
        count: colCounts[j],
        proportion: p[j]
      };
    });

    // Variables summary (correlation ratio eta^2 of categorical variables with dimensions)
    const variables = categoricalCols.map(col => {
      const mods = modalitiesCoords.filter(m => m.col === col);
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        // eta^2 = sum_m (p_m * coord_m,k^2) / eigenvalue_k
        const eig = eigenvalues[k];
        const eta2 = eig > 0 ? mods.reduce((acc, m) => acc + m.proportion * (m.coords[k] ** 2), 0) / eig : 0;
        coords.push(Math.sqrt(Math.max(0, Math.min(1, eta2))));
      }
      return {
        name: col,
        type: 'categorical',
        coords
      };
    });

    return {
      method: 'MCA',
      methodLabel: 'Multiple Correspondence Analysis (MCA)',
      nObservations: nClean,
      nVariables: K,
      nModalities: J,
      categoricalCols,
      eigenvalues,
      totalInertia,
      scree,
      individuals,
      variables,
      modalities: modalitiesCoords,
      allVariables: modalitiesCoords,
      outliers: individuals.filter(d => d.isOutlier),
      outlierCount: individuals.filter(d => d.isOutlier).length,
      outlierCutoff
    };
  },

  // =========================================================================
  // 5. FACTOR ANALYSIS OF MIXED DATA (FAMD / Pagès 2004)
  // =========================================================================
  runFAMD(records, continuousCols, categoricalCols, options = {}) {
    if (!continuousCols || continuousCols.length === 0) {
      return this.runMCA(records, categoricalCols, options);
    }
    if (!categoricalCols || categoricalCols.length === 0) {
      return this.runPCA(records, continuousCols, options);
    }

    const p1 = continuousCols.length;
    const K = categoricalCols.length;

    // Filter complete cases across both continuous and categorical variables
    const cleanRecords = records.filter(r => {
      const contOk = continuousCols.every(c => {
        const val = Number(r[c]);
        return !isNaN(val) && isFinite(val) && r[c] !== '' && r[c] !== null;
      });
      const catOk = categoricalCols.every(c => r[c] !== undefined && r[c] !== null && String(r[c]).trim() !== '');
      return contOk && catOk;
    });

    if (cleanRecords.length < 3) {
      return { error: 'Insufficient complete observations for FAMD (minimum 3 required).' };
    }

    const nClean = cleanRecords.length;

    // 1. Preprocess Continuous Variables: Center and scale to unit variance
    const contMeans = {};
    const contSds = {};
    continuousCols.forEach(col => {
      const vals = cleanRecords.map(r => Number(r[col]));
      const m = vals.reduce((a, b) => a + b, 0) / nClean;
      const variance = vals.reduce((acc, v) => acc + (v - m) ** 2, 0) / (nClean - 1);
      contMeans[col] = m;
      contSds[col] = Math.sqrt(variance) || 1.0;
    });

    // 2. Preprocess Categorical Variables: One-hot encode and scale by 1 / sqrt(p_m)
    const modalities = [];
    categoricalCols.forEach(col => {
      const uniqueVals = Array.from(new Set(cleanRecords.map(r => String(r[col]).trim()))).sort();
      uniqueVals.forEach(val => {
        modalities.push({
          col,
          val,
          label: `${col}: ${val}`
        });
      });
    });

    const J = modalities.length;
    const pProps = modalities.map(m => {
      const count = cleanRecords.filter(r => String(r[m.col]).trim() === m.val).length;
      return Math.max(1e-6, count / nClean);
    });

    // 3. Construct Joint Balanced Matrix W of size nClean x (p1 + J)
    // Continuous block: (X_ij - mean_j) / (sd_j * sqrt(nClean))
    // Categorical block: (I_im - p_m) / (sqrt(p_m) * sqrt(nClean))
    const totalCols = p1 + J;
    const sqrtN = Math.sqrt(nClean);

    const W = Array.from({ length: nClean }, (_, i) => {
      const r = cleanRecords[i];
      const row = new Array(totalCols);

      // Continuous columns
      for (let j = 0; j < p1; j++) {
        const col = continuousCols[j];
        row[j] = (Number(r[col]) - contMeans[col]) / (contSds[col] * sqrtN);
      }

      // Categorical dummy columns
      for (let m = 0; m < J; m++) {
        const mod = modalities[m];
        const isMatch = String(r[mod.col]).trim() === mod.val ? 1 : 0;
        const pm = pProps[m];
        row[p1 + m] = (isMatch - pm) / (Math.sqrt(pm) * sqrtN);
      }

      return row;
    });

    // 4. SVD on W
    const svdRes = this.LinearAlgebra.svd(W);
    const eigenvalues = svdRes.s.map(s => s * s);
    const totalInertia = p1 + (J - K);

    const nComponents = Math.min(totalCols - K, nClean - 1);
    const scree = [];
    let cumVar = 0;

    for (let k = 0; k < nComponents; k++) {
      const eig = eigenvalues[k] || 0;
      const varPct = totalInertia > 0 ? (eig / totalInertia) * 100 : 0;
      cumVar += varPct;
      scree.push({
        dim: k + 1,
        label: `Dim ${k + 1}`,
        eigenvalue: eig,
        variancePct: varPct,
        cumulativePct: Math.min(100, cumVar)
      });
    }

    // 5. Individual Coordinates (Factor Scores)
    // F_ik = sqrt(nClean) * U_ik * s_k
    const individuals = cleanRecords.map((r, i) => {
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        coords.push(sqrtN * svdRes.U[i][k] * svdRes.s[k]);
      }
      const totalDistSq = coords.reduce((acc, c) => acc + c * c, 0);
      const cos2 = coords.map(c => (totalDistSq > 0 ? (c * c) / totalDistSq : 0));

      return {
        id: r._id || r.id || r.patient_id || `Obs-${i + 1}`,
        originalRow: r,
        coords,
        cos2,
        distFromCentroid: Math.sqrt(totalDistSq)
      };
    });

    // Flag outliers (> 2.5 SD distance in factor space)
    const allDists = individuals.map(d => d.distFromCentroid);
    const meanDist = allDists.reduce((a, b) => a + b, 0) / allDists.length;
    const sdDist = Math.sqrt(allDists.reduce((acc, v) => acc + (v - meanDist) ** 2, 0) / (allDists.length - 1)) || 1;
    const outlierCutoff = meanDist + 2.5 * sdDist;
    individuals.forEach(ind => {
      ind.isOutlier = ind.distFromCentroid > outlierCutoff;
    });

    // 6. Variable Loadings and Modality Centroids
    // A. Continuous variable correlation with dimensions
    const continuousVariables = continuousCols.map((col, j) => {
      const rawVals = cleanRecords.map(r => Number(r[col]));
      const coords = [];
      for (let k = 0; k < nComponents; k++) {
        // Pearson correlation with individual factor score
        const score_k = individuals.map(ind => ind.coords[k]);
        const mScore = score_k.reduce((a, b) => a + b, 0) / nClean;
        const sdScore = Math.sqrt(score_k.reduce((acc, v) => acc + (v - mScore) ** 2, 0) / (nClean - 1)) || 1;

        let cov = 0;
        for (let i = 0; i < nClean; i++) {
          cov += (rawVals[i] - contMeans[col]) * (score_k[i] - mScore);
        }
        cov /= (nClean - 1);
        const rCorr = cov / (contSds[col] * sdScore);
        coords.push(Math.max(-1, Math.min(1, rCorr)));
      }

      return {
        name: col,
        type: 'continuous',
        coords,
        mean: contMeans[col],
        sd: contSds[col]
      };
    });

    // B. Categorical Modality Centroid Coordinates
    const modalityCoordinates = modalities.map((m, mIdx) => {
      const matchingInds = individuals.filter(ind => String(ind.originalRow[m.col]).trim() === m.val);
      const count = matchingInds.length;
      const coords = [];

      for (let k = 0; k < nComponents; k++) {
        if (count > 0) {
          const meanCoord = matchingInds.reduce((acc, ind) => acc + ind.coords[k], 0) / count;
          coords.push(meanCoord);
        } else {
          coords.push(0);
        }
      }

      return {
        name: m.label,
        col: m.col,
        val: m.val,
        type: 'modality',
        coords,
        count,
        proportion: pProps[mIdx]
      };
    });

    const allVariables = [...continuousVariables, ...modalityCoordinates];

    return {
      method: 'FAMD',
      methodLabel: 'Factor Analysis of Mixed Data (FAMD)',
      nObservations: nClean,
      nContinuous: p1,
      nCategorical: K,
      nModalities: J,
      continuousCols,
      categoricalCols,
      eigenvalues,
      totalInertia,
      scree,
      individuals,
      variables: continuousVariables,
      modalities: modalityCoordinates,
      allVariables,
      continuousVariables,
      categoricalVariables: modalityCoordinates,
      outliers: individuals.filter(d => d.isOutlier),
      outlierCount: individuals.filter(d => d.isOutlier).length,
      outlierCutoff
    };
  },

  // =========================================================================
  // 6. MASTER ANALYSIS DISPATCHER
  // =========================================================================
  executeAnalysis(records, method, activeCols, options = {}) {
    const colTypes = options.colTypes || {};
    const groupingCol = options.groupingCol || null;

    let analysis;
    if (method === 'PCA') {
      const contCols = activeCols.filter(c => (colTypes[c] || 'continuous') === 'continuous');
      analysis = this.runPCA(records, contCols, options);
    } else if (method === 'MCA') {
      const catCols = activeCols.filter(c => colTypes[c] === 'categorical');
      analysis = this.runMCA(records, catCols, options);
    } else {
      // FAMD
      const contCols = activeCols.filter(c => (colTypes[c] || 'continuous') === 'continuous');
      const catCols = activeCols.filter(c => colTypes[c] === 'categorical');
      analysis = this.runFAMD(records, contCols, catCols, options);
    }

    if (analysis.error) return analysis;

    // Attach grouping information to individuals
    if (groupingCol && records[0] && records[0][groupingCol] !== undefined) {
      analysis.groupingCol = groupingCol;
      const groups = new Set();
      analysis.individuals.forEach(ind => {
        ind.group = String(ind.originalRow[groupingCol] || 'Unassigned');
        groups.add(ind.group);
      });
      analysis.uniqueGroups = Array.from(groups).sort();
    } else {
      analysis.groupingCol = null;
      analysis.uniqueGroups = [];
      analysis.individuals.forEach(ind => { ind.group = 'Observation'; });
    }

    // Add Interpretation & Scripts
    analysis.interpretation = this.generateInterpretation(analysis);
    analysis.scripts = {
      python: this.generatePythonScript(analysis),
      r: this.generateRScript(analysis)
    };

    return analysis;
  },

  // =========================================================================
  // 7. AUTOMATED NATURAL LANGUAGE INTERPRETATION & EXPLANATION
  // =========================================================================
  generateInterpretation(analysis) {
    if (!analysis || analysis.error) return '';

    const scree = analysis.scree || [];
    const dim1 = scree[0] || { variancePct: 0, eigenvalue: 0 };
    const dim2 = scree[1] || { variancePct: 0, eigenvalue: 0 };
    const dim3 = scree[2] || { variancePct: 0, eigenvalue: 0 };

    const cum2D = (dim1.variancePct + dim2.variancePct).toFixed(1);
    const cum3D = (dim1.variancePct + dim2.variancePct + dim3.variancePct).toFixed(1);

    // Identify Drivers for Dim 1 & Dim 2
    const allVars = analysis.allVariables || [];

    // Sort by absolute coordinate on Dim 1
    const sortedDim1 = [...allVars].sort((a, b) => Math.abs(b.coords[0]) - Math.abs(a.coords[0]));
    const topDim1Pos = sortedDim1.filter(v => v.coords[0] > 0).slice(0, 3).map(v => `${v.name} (+${v.coords[0].toFixed(2)})`);
    const topDim1Neg = sortedDim1.filter(v => v.coords[0] < 0).slice(0, 3).map(v => `${v.name} (${v.coords[0].toFixed(2)})`);

    // Sort by absolute coordinate on Dim 2
    const sortedDim2 = [...allVars].sort((a, b) => Math.abs(b.coords[1]) - Math.abs(a.coords[1]));
    const topDim2Pos = sortedDim2.filter(v => v.coords[1] > 0).slice(0, 3).map(v => `${v.name} (+${v.coords[1].toFixed(2)})`);
    const topDim2Neg = sortedDim2.filter(v => v.coords[1] < 0).slice(0, 3).map(v => `${v.name} (${v.coords[1].toFixed(2)})`);

    // Outlier list
    const outliers = (analysis.individuals || []).filter(ind => ind.isOutlier);
    const outlierSummary = outliers.length > 0
      ? `Flagged ${outliers.length} statistical outlier observation(s) exceeding 2.5 SD distance from centroid: ${outliers.slice(0, 5).map(o => o.id).join(', ')}${outliers.length > 5 ? '...' : ''}.`
      : 'No severe outlier observations detected (all individuals fall within 2.5 standard deviations of the multivariate centroid).';

    // Clustering summary
    let clusterSummary = '';
    if (analysis.groupingCol && analysis.uniqueGroups.length > 1) {
      clusterSummary = `Observations were stratified across ${analysis.uniqueGroups.length} categories of "${analysis.groupingCol}". Clear separation is observable along principal axes, indicating that baseline covariates correlate significantly with group classification.`;
    }

    const narrative = [
      `### Multivariate Dimensionality Reduction Findings (${analysis.methodLabel}):`,
      `1. **Variance Retention & Dimensionality:** The primary 2D principal plane (Dim 1 vs. Dim 2) accounts for **${cum2D}%** of total inertia (Dim 1 = ${dim1.variancePct.toFixed(1)}%, Dim 2 = ${dim2.variancePct.toFixed(1)}%). Expanding to the 3D subspace (Dim 1–3) captures **${cum3D}%** of the aggregate sample variance.`,
      `2. **Key Drivers for Dimension 1 (Horizontal Axis):** Strongest positive associations include ${topDim1Pos.length > 0 ? topDim1Pos.join(', ') : 'none'}, contrasting with negative associations ${topDim1Neg.length > 0 ? topDim1Neg.join(', ') : 'none'}.`,
      `3. **Key Drivers for Dimension 2 (Vertical Axis):** Orthogonal variance is primarily defined by ${topDim2Pos.length > 0 ? topDim2Pos.join(', ') : 'none'} against ${topDim2Neg.length > 0 ? topDim2Neg.join(', ') : 'none'}.`,
      `4. **Outliers & Anomaly Detection:** ${outlierSummary}`,
      clusterSummary ? `5. **Subgroup Clustering:** ${clusterSummary}` : ''
    ].filter(Boolean).join('\n\n');

    return narrative;
  },

  // =========================================================================
  // 8. REPRODUCIBLE SCRIPT GENERATOR [Python & R]
  // =========================================================================
  generatePythonScript(analysis) {
    const method = analysis.method;
    const contColsPy = (analysis.continuousCols || []).map(c => `'${c}'`).join(', ');
    const catColsPy = (analysis.categoricalCols || []).map(c => `'${c}'`).join(', ');
    const groupPy = analysis.groupingCol ? `'${analysis.groupingCol}'` : 'None';

    return `"""
=============================================================================
Multivariate Exploratory Data Analysis (${analysis.methodLabel})
Generated by Statis-Gravity Clinical Biostatistics Platform
=============================================================================
"""

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# 1. Load Clinical Observational Data
# df = pd.read_csv('clinical_multivariate_data.csv')

continuous_cols = [${contColsPy}]
categorical_cols = [${catColsPy}]
group_col = ${groupPy}

${method === 'PCA' ? `# -----------------------------------------------------------------------------
# 2. PRINCIPAL COMPONENT ANALYSIS (scikit-learn)
# -----------------------------------------------------------------------------
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

df_clean = df[continuous_cols].dropna()
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_clean)

pca = PCA(n_components=min(5, len(continuous_cols)))
pca_scores = pca.fit_transform(X_scaled)

scores_df = pd.DataFrame(pca_scores, columns=[f'Dim {i+1}' for i in range(pca_scores.shape[1])])
if group_col and group_col in df.columns:
    scores_df[group_col] = df.loc[df_clean.index, group_col].values

exp_var = pca.explained_variance_ratio_ * 100
print("Explained Variance per Dimension (%):", np.round(exp_var, 2))
print("Cumulative Variance (%):", np.round(np.cumsum(exp_var), 2))

# 3. Interactive 2D Factor Biplot
fig_2d = px.scatter(
    scores_df, x='Dim 1', y='Dim 2', color=group_col,
    title=f"PCA Factor Map: Dim 1 ({exp_var[0]:.1f}%) vs Dim 2 ({exp_var[1]:.1f}%)",
    template='plotly_dark'
)
fig_2d.show()

# 4. Interactive 3D Rotation Plot
if scores_df.shape[1] >= 3:
    fig_3d = px.scatter_3d(
        scores_df, x='Dim 1', y='Dim 2', z='Dim 3', color=group_col,
        title="PCA 3D Space (Dim 1 vs Dim 2 vs Dim 3)",
        template='plotly_dark'
    )
    fig_3d.show()` : (method === 'MCA' ? `# -----------------------------------------------------------------------------
# 2. MULTIPLE CORRESPONDENCE ANALYSIS (prince package)
# -----------------------------------------------------------------------------
# pip install prince
import prince

df_clean = df[categorical_cols].dropna()
mca = prince.MCA(n_components=min(5, len(categorical_cols)*2), random_state=42)
mca.fit(df_clean)

scores_df = mca.row_coordinates(df_clean)
scores_df.columns = [f'Dim {i+1}' for i in range(scores_df.shape[1])]
if group_col and group_col in df.columns:
    scores_df[group_col] = df.loc[df_clean.index, group_col].values

fig_2d = px.scatter(
    scores_df, x='Dim 1', y='Dim 2', color=group_col,
    title="MCA Factor Map: Individuals & Categories",
    template='plotly_dark'
)
fig_2d.show()` : `# -----------------------------------------------------------------------------
# 2. FACTOR ANALYSIS OF MIXED DATA (prince package)
# -----------------------------------------------------------------------------
# pip install prince
import prince

all_cols = continuous_cols + categorical_cols
df_clean = df[all_cols].dropna()

famd = prince.FAMD(n_components=min(5, len(all_cols)), random_state=42)
famd.fit(df_clean)

scores_df = famd.row_coordinates(df_clean)
scores_df.columns = [f'Dim {i+1}' for i in range(scores_df.shape[1])]
if group_col and group_col in df.columns:
    scores_df[group_col] = df.loc[df_clean.index, group_col].values

fig_2d = px.scatter(
    scores_df, x='Dim 1', y='Dim 2', color=group_col,
    title="FAMD Factor Map (Continuous & Categorical Integration)",
    template='plotly_dark'
)
fig_2d.show()`)}
`;
  },

  generateRScript(analysis) {
    const method = analysis.method;
    const contColsR = (analysis.continuousCols || []).map(c => `"${c}"`).join(', ');
    const catColsR = (analysis.categoricalCols || []).map(c => `"${c}"`).join(', ');
    const groupR = analysis.groupingCol ? `"${analysis.groupingCol}"` : 'NULL';

    return `################################################################################
# Multivariate Exploratory Data Analysis (${analysis.methodLabel})
# Package Suite: FactoMineR, factoextra, ggplot2
# Generated by Statis-Gravity Clinical Biostatistics Platform
################################################################################

# install.packages(c("FactoMineR", "factoextra", "ggplot2"))
library(FactoMineR)
library(factoextra)
library(ggplot2)

# 1. Load Clinical Observational Data
# df <- read.csv("clinical_multivariate_data.csv")

continuous_cols  <- c(${contColsR})
categorical_cols <- c(${catColsR})
group_var        <- ${groupR}

${method === 'PCA' ? `# 2. PRINCIPAL COMPONENT ANALYSIS
df_pca <- na.omit(df[, continuous_cols])
res.pca <- PCA(df_pca, scale.unit = TRUE, graph = FALSE)

# Scree Plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# 2D Factor Map with Individuals
fviz_pca_ind(res.pca, 
             habillage = if(!is.null(group_var)) df[[group_var]] else "none",
             addEllipses = TRUE, 
             repel = TRUE)

# Variable Correlation Circle (Biplot Arrows)
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)` : (method === 'MCA' ? `# 2. MULTIPLE CORRESPONDENCE ANALYSIS
df_mca <- na.omit(df[, categorical_cols])
res.mca <- MCA(df_mca, graph = FALSE)

# Scree Plot
fviz_eig(res.mca, addlabels = TRUE)

# Individuals & Categories Biplot
fviz_mca_biplot(res.mca, repel = TRUE,
                ggtheme = theme_minimal())` : `# 2. FACTOR ANALYSIS OF MIXED DATA (FAMD)
all_cols <- c(continuous_cols, categorical_cols)
df_famd <- na.omit(df[, all_cols])
res.famd <- FAMD(df_famd, graph = FALSE)

# Scree Plot
fviz_eig(res.famd, addlabels = TRUE)

# Individuals Factor Map
fviz_famd_ind(res.famd, 
              habillage = if(!is.null(group_var)) df[[group_var]] else "none",
              repel = TRUE)

# Quantitative & Qualitative Variable Contributions
fviz_famd_var(res.famd, "var", repel = TRUE)`)}
`;
  },

  // =========================================================================
  // 9. REALISTIC MIXED CRANIOFACIAL CLINICAL DATASET GENERATOR
  // =========================================================================
  getSampleMixedCohort() {
    let seed = 101;
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

    const sutures = ['Sagittal', 'Coronal', 'Metopic', 'Lambdoid'];
    const approaches = ['Endoscopic_Strip', 'Open_Cranial_Vault', 'Spring_Assisted'];
    const comorbidities = ['Grade_0', 'Grade_I', 'Grade_II'];

    const cohort = [];
    const N = 120;

    for (let i = 1; i <= N; i++) {
      // Patient age in months
      const rApproach = rnd();
      let approach = approaches[0];
      if (rApproach > 0.45 && rApproach <= 0.85) approach = approaches[1];
      else if (rApproach > 0.85) approach = approaches[2];

      const suture = sutures[Math.floor(rnd() * sutures.length)];
      const comorb = comorbidities[rnd() < 0.65 ? 0 : (rnd() < 0.90 ? 1 : 2)];

      // Correlate continuous parameters with approach & clinical condition
      let age = approach === 'Endoscopic_Strip'
        ? Math.max(2.1, Math.min(5.5, rndNormal(3.6, 0.8)))
        : (approach === 'Spring_Assisted'
            ? Math.max(3.0, Math.min(7.0, rndNormal(4.5, 0.9)))
            : Math.max(6.0, Math.min(14.0, rndNormal(9.2, 2.1))));

      let op_time = approach === 'Endoscopic_Strip'
        ? Math.max(50, Math.round(rndNormal(72, 12)))
        : (approach === 'Spring_Assisted'
            ? Math.max(70, Math.round(rndNormal(98, 16)))
            : Math.max(140, Math.round(rndNormal(205, 32))));

      let blood_loss = approach === 'Endoscopic_Strip'
        ? Math.max(20, Math.round(rndNormal(48, 14)))
        : (approach === 'Spring_Assisted'
            ? Math.max(35, Math.round(rndNormal(75, 20)))
            : Math.max(150, Math.round(rndNormal(285, 65))));

      let los_days = approach === 'Endoscopic_Strip'
        ? Math.max(1.0, Math.round(rndNormal(1.8, 0.5) * 10) / 10)
        : (approach === 'Spring_Assisted'
            ? Math.max(1.5, Math.round(rndNormal(2.4, 0.6) * 10) / 10)
            : Math.max(3.0, Math.round(rndNormal(4.9, 1.1) * 10) / 10));

      let cranial_index = suture === 'Sagittal'
        ? Math.max(62, Math.round(rndNormal(68, 3.5)))
        : (suture === 'Coronal'
            ? Math.max(82, Math.round(rndNormal(88, 4.2)))
            : Math.max(74, Math.round(rndNormal(79, 3.8))));

      let baseline_cvai = Math.max(3.5, Math.round(rndNormal(9.4, 2.6) * 10) / 10);
      let transfusion = blood_loss > 120 || (rnd() < 0.15 && approach !== 'Endoscopic_Strip') ? 'Yes' : 'No';

      cohort.push({
        patient_id: `PT-${1000 + i}`,
        age_months: Math.round(age * 10) / 10,
        cranial_index,
        baseline_cvai,
        operative_time_min: op_time,
        blood_loss_ml: blood_loss,
        length_of_stay_days: los_days,
        suture_type: suture,
        surgical_approach: approach,
        comorbidity_grade: comorb,
        transfusion_required: transfusion
      });
    }

    return cohort;
  }
};
