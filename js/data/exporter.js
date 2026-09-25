/**
 * Statis-Gravity - Clinical Report Generator & Data Exporter
 * Generates APA/ICMJE compliant statistical summaries and downloadable exports.
 */

import { DocxReports } from '../export/docx-generator.js';

export const Exporter = {
  /**
   * Formats a p-value according to AMA/APA guidelines:
   * e.g., p < .001, or p = .034
   */
  formatP(p) {
    if (isNaN(p) || p === null) return 'N/A';
    if (p < 0.001) return 'p < .001';
    return `p = ${p.toFixed(3).replace(/^0\./, '.')}`;
  },

  /**
   * Generates APA-style statistical statement for two-group t-test
   */
  formatTTestReport(res) {
    if (res.error) return res.error;
    const sign = res.isSignificant ? 'statistically significant' : 'not statistically significant';
    return (
      `An independent-samples t-test indicated that the difference between ${res.groupA.name || 'Group A'} ` +
      `(M = ${res.groupA.mean.toFixed(2)}, SD = ${res.groupA.sd.toFixed(2)}, n = ${res.groupA.n}) and ` +
      `${res.groupB.name || 'Group B'} (M = ${res.groupB.mean.toFixed(2)}, SD = ${res.groupB.sd.toFixed(2)}, n = ${res.groupB.n}) ` +
      `was ${sign}, t(${res.df.toFixed(1)}) = ${res.statistic.toFixed(2)}, ${this.formatP(res.pValue)}, ` +
      `Cohen's d = ${res.cohensD.toFixed(2)}, 95% CI for difference [${res.ci95[0].toFixed(2)}, ${res.ci95[1].toFixed(2)}].`
    );
  },

  /**
   * Generates APA-style statement for ANOVA
   */
  formatAnovaReport(res) {
    if (res.error) return res.error;
    const sign = res.isSignificant ? 'statistically significant' : 'not statistically significant';
    return (
      `A one-way between-subjects ANOVA was conducted to compare the effect of intervention across ${res.k} groups. ` +
      `There was a ${sign} effect, F(${res.dfBetween}, ${res.dfWithin}) = ${res.fStatistic.toFixed(2)}, ${this.formatP(res.pValue)}, ` +
      `η² = ${res.etaSquared.toFixed(3)}, ω² = ${res.omegaSquared.toFixed(3)}.`
    );
  },

  /**
   * Generates APA-style statement for Chi-Square and Risk Metrics
   */
  formatCategoricalReport(res) {
    if (res.error) return res.error;
    const t = res.table;
    const r = res.riskMetrics;
    return (
      `A 2x2 contingency analysis (N = ${t.n}) revealed a Pearson Chi-Square of χ²(1) = ${res.chiSquare.standard.toFixed(2)}, ` +
      `${this.formatP(res.chiSquare.pValueStandard)} (Fisher's exact test ${this.formatP(res.fishersExact.pValue)}). ` +
      `The Odds Ratio was ${r.oddsRatio.toFixed(2)} (95% CI [${r.orCI95[0].toFixed(2)}, ${r.orCI95[1].toFixed(2)}]), ` +
      `Relative Risk = ${r.relativeRisk.toFixed(2)} (95% CI [${r.rrCI95[0].toFixed(2)}, ${r.rrCI95[1].toFixed(2)}]), ` +
      `with Absolute Risk Reduction (ARR) = ${(r.arr * 100).toFixed(1)}% (NNT = ${r.nnt.toFixed(1)}).`
    );
  },

  /**
   * Generates APA/ICMJE-style descriptive summary
   */
  formatDescriptiveReport(res) {
    if (res.error) return res.error;
    let modeStr = 'no distinct mode (all values unique)';
    if (res.modes && res.modes.length > 0) {
      modeStr = `mode = ${res.modes.map(m => m.toFixed(2)).join(', ')} (freq = ${res.maxFreq || 2})`;
    }
    const normStr = res.normality.isNormal
      ? `normally distributed via Jarque-Bera omnibus test (JB = ${res.normality.statistic.toFixed(2)}, df = 2, ${this.formatP(res.normality.pValue)}; skewness = ${res.skewness.toFixed(2)} [${res.skewnessInterpretation || 'symmetric'}], excess kurtosis = ${res.kurtosis.toFixed(2)} [${res.kurtosisInterpretation || 'mesokurtic'}])`
      : `significantly non-normal via Jarque-Bera omnibus test (JB = ${res.normality.statistic.toFixed(2)}, df = 2, ${this.formatP(res.normality.pValue)}; skewness = ${res.skewness.toFixed(2)} [${res.skewnessInterpretation || 'skewed'}], excess kurtosis = ${res.kurtosis.toFixed(2)} [${res.kurtosisInterpretation || 'leptokurtic'}])`;

    let outlierStr = 'No outliers detected beyond Tukey fences (1.5 × IQR)';
    if (res.outliers && res.outliers.length > 0) {
      outlierStr = `${res.outliers.length} outlier(s) detected beyond Tukey fences [${res.lowerFence.toFixed(1)}, ${res.upperFence.toFixed(1)}]: ` +
        res.outliers.map(o => `${o.value.toFixed(1)} (Z = ${o.zScore > 0 ? '+' : ''}${o.zScore.toFixed(2)}, ${o.type})`).join(', ');
    }

    return (
      `Continuous variable summary (N = ${res.n}): ` +
      `Mean = ${res.mean.toFixed(2)} (SD = ${res.sd.toFixed(2)}, 95% CI [${res.ci95[0].toFixed(2)}, ${res.ci95[1].toFixed(2)}]); ` +
      `Median = ${res.median.toFixed(2)} (IQR = ${res.iqr.toFixed(2)}, Q1 = ${res.q1.toFixed(2)}, Q3 = ${res.q3.toFixed(2)}); ` +
      `${modeStr}; Range = ${res.min.toFixed(2)} to ${res.max.toFixed(2)}. ` +
      `Distribution: ${normStr}. ` +
      `Outliers: ${outlierStr}. ` +
      (res.outliers && res.outliers.length > 0 || !res.normality.isNormal ? 'Reporting Median and IQR is recommended due to skewness/outliers.' : 'Reporting Mean and SD is appropriate.')
    );
  },

  /**
   * Generates APA-style statement for Diagnostic / ROC performance
   */
  formatDiagnosticReport(eval2x2, roc) {
    let text = `Diagnostic Evaluation (N = ${eval2x2.total}): ` +
      `Sensitivity = ${(eval2x2.sensitivity * 100).toFixed(1)}% (95% CI [${(eval2x2.sensitivityCI95[0]*100).toFixed(1)}%, ${(eval2x2.sensitivityCI95[1]*100).toFixed(1)}%]), ` +
      `Specificity = ${(eval2x2.specificity * 100).toFixed(1)}% (95% CI [${(eval2x2.specificityCI95[0]*100).toFixed(1)}%, ${(eval2x2.specificityCI95[1]*100).toFixed(1)}%]), ` +
      `PPV = ${(eval2x2.ppv * 100).toFixed(1)}%, NPV = ${(eval2x2.npv * 100).toFixed(1)}%, ` +
      `Positive LR = ${eval2x2.plr.toFixed(2)}, Negative LR = ${eval2x2.nlr.toFixed(2)}, ` +
      `Youden's Index J = ${eval2x2.youdenJ.toFixed(3)}.`;

    if (roc && roc.auc) {
      text += ` ROC Area Under Curve (AUC) = ${roc.auc.toFixed(3)} (95% CI [${roc.aucCI95[0].toFixed(3)}, ${roc.aucCI95[1].toFixed(3)}]).`;
    }
    return text;
  },

  /**
   * Copies formatted text to system clipboard
   */
  async copyToClipboard(text) {
    try {
      await navigator.clipboard.writeText(text);
      return true;
    } catch {
      // Fallback
      const ta = document.createElement('textarea');
      ta.value = text;
      document.body.appendChild(ta);
      ta.select();
      document.execCommand('copy');
      document.body.removeChild(ta);
      return true;
    }
  },

  /**
   * Downloads data as CSV file
   */
  downloadCSV(filename, content) {
    const blob = new Blob([content], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    link.click();
  },

  /**
   * Downloads a Blob object with automatic URL cleanup
   */
  downloadBlob(filename, blob) {
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  },

  /**
   * Exports specified tab analysis data to .docx
   */
  exportTabToDocx(tabId, data) {
    let builder;
    switch (tabId) {
      case 'descriptive':
        builder = DocxReports.createDescriptiveDocx(data);
        break;
      case 'hypothesis':
        builder = DocxReports.createHypothesisDocx(data);
        break;
      case 'anova':
        builder = DocxReports.createAnovaDocx(data);
        break;
      case 'categorical':
        builder = DocxReports.createCategoricalDocx(data);
        break;
      case 'correlation':
        builder = DocxReports.createCorrelationDocx(data);
        break;
      case 'diagnostic':
        builder = DocxReports.createDiagnosticDocx(data);
        break;
      case 'power':
        builder = DocxReports.createPowerDocx(data);
        break;
      case 'teaching':
        builder = DocxReports.createTeachingDocx(data);
        break;
      case 'teaching-bayesian':
        builder = DocxReports.createTeachingBayesianDocx ? DocxReports.createTeachingBayesianDocx(data) : DocxReports.createTeachingDocx(data);
        break;
      case 'propensity':
        builder = DocxReports.createPsmDocx(data);
        break;
      case 'multivariate':
        builder = DocxReports.createMultivariateDocx(data);
        break;
      default:
        console.error('Unknown tab for DOCX export:', tabId);
        return;
    }
    const blob = builder.generateBlob();
    this.downloadBlob(`statis-gravity-${tabId}-report.docx`, blob);
  }
};
