/**
 * Statis-Gravity - Data Ingestion & Sample Datasets
 * Parses CSV/TSV, whitespace-separated text, and supplies realistic clinical datasets.
 */

export const DataParser = {
  /**
   * Parses string containing numbers (separated by comma, newline, tab, or space)
   */
  parseSeries(text) {
    if (!text || typeof text !== 'string') return [];
    return text
      .split(/[\r\n,;\t\s]+/)
      .map(s => parseFloat(s.trim()))
      .filter(n => !isNaN(n) && isFinite(n));
  },

  /**
   * Parses 2-column or Multi-column CSV/TSV
   */
  parseTable(text) {
    if (!text || typeof text !== 'string') return { headers: [], columns: [] };
    const lines = text.trim().split(/\r?\n/).filter(line => line.trim().length > 0);
    if (lines.length === 0) return { headers: [], columns: [] };

    // Detect delimiter: comma, tab, or semicolon
    const firstLine = lines[0];
    const delimiter = firstLine.includes('\t') ? '\t' : firstLine.includes(';') ? ';' : ',';

    const rawRows = lines.map(line => line.split(delimiter).map(cell => cell.trim()));
    const firstRowHasStrings = rawRows[0].some(cell => isNaN(parseFloat(cell)) && cell.length > 0);

    let headers = [];
    let dataRows = [];

    if (firstRowHasStrings) {
      headers = rawRows[0];
      dataRows = rawRows.slice(1);
    } else {
      headers = rawRows[0].map((_, idx) => `Var ${idx + 1}`);
      dataRows = rawRows;
    }

    const columns = headers.map((name, colIdx) => ({
      name,
      values: dataRows
        .map(row => parseFloat(row[colIdx]))
        .filter(v => !isNaN(v) && isFinite(v))
    }));

    return { headers, columns };
  },

  /**
   * Pre-packaged realistic clinical datasets
   */
  samples: {
    icpDynamics: {
      name: 'Neurosurgical ICP Dynamics (Marmarou Infusion Study)',
      description: 'Continuous intracranial pressure (mmHg) measured pre- and post-bolus CSF infusion.',
      groupA: {
        name: 'Pre-Infusion Baseline (mmHg)',
        data: [10.2, 11.5, 9.8, 12.1, 10.9, 13.4, 11.0, 9.5, 12.8, 10.4, 11.7, 10.0]
      },
      groupB: {
        name: 'Post-Infusion Plateau (mmHg)',
        data: [18.5, 21.0, 19.2, 23.4, 20.1, 25.6, 22.0, 19.8, 24.5, 20.3, 21.8, 19.5]
      }
    },

    drainOutputSkewed: {
      name: 'Postoperative Subdural Drain Output (mL/24h)',
      description: 'Right-skewed postoperative drain volumes with mild and extreme outlier spikes.',
      data: [12, 14, 15, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 38, 62]
    },

    cranialAsymmetry: {
      name: 'Cranial Vault Remodeling (CVAI % Across Cohorts)',
      description: 'Cranial Vault Asymmetry Index (CVAI %) across 3 clinical management pathways.',
      groups: [
        { name: 'Conservative Repositioning', data: [7.2, 6.8, 7.5, 6.9, 8.1, 7.0, 7.4, 6.5, 7.9, 7.1] },
        { name: 'Orthotic Helmet Therapy', data: [4.1, 3.8, 4.5, 3.9, 4.8, 3.6, 4.2, 3.5, 4.0, 3.7] },
        { name: 'Endoscopic Strip Craniectomy', data: [2.5, 2.8, 2.2, 2.6, 3.1, 2.4, 2.9, 2.1, 2.7, 2.3] }
      ]
    },

    diagnosticBiomarker: {
      name: 'CSF Neurofilament Light (NfL) in iNPH vs Controls',
      description: 'Biomarker concentration (pg/mL) and verified clinical outcome status (1 = iNPH Responder, 0 = Non-responder).',
      data: [
        { score: 1850, status: 1 },
        { score: 1720, status: 1 },
        { score: 1640, status: 1 },
        { score: 1590, status: 1 },
        { score: 1530, status: 1 },
        { score: 1480, status: 1 },
        { score: 1390, status: 1 },
        { score: 1320, status: 1 },
        { score: 1280, status: 1 },
        { score: 1210, status: 1 },
        { score: 1150, status: 1 },
        { score: 1100, status: 0 },
        { score: 1040, status: 1 },
        { score: 980, status: 0 },
        { score: 920, status: 0 },
        { score: 860, status: 0 },
        { score: 810, status: 0 },
        { score: 750, status: 0 },
        { score: 690, status: 0 },
        { score: 620, status: 0 }
      ]
    },

    shuntComplications: {
      name: 'CSF Shunt Valve Type vs Proximal Occlusion',
      description: 'Multicenter trial 2x2 contingency cohort comparing programmable vs fixed-pressure valves.',
      a: 14, // Programmable + Occlusion
      b: 36, // Fixed + Occlusion
      c: 186, // Programmable + Patency maintained
      d: 164 // Fixed + Patency maintained
    },

    shunt2x2: { a: 14, b: 36, c: 186, d: 164 },

    diagnostic2x2: {
      name: 'Rapid Point-of-Care Biomarker vs Gold Standard Reference',
      description: 'Prospective diagnostic cohort evaluating novel serum biomarker test against reference pathology.',
      a: 92, // True Positive (TP)
      b: 8,  // False Positive (FP)
      c: 12, // False Negative (FN)
      d: 188 // True Negative (TN)
    }
  }
};
