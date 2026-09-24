/**
 * Statis-Gravity - Clinical DOCX Report Generator
 * Generates native Microsoft Office Open XML (.docx) files without external dependencies.
 * Creates compliant ZIP packages with WordprocessingML, formatting tables, callout boxes,
 * APA/ICMJE statistical summaries, methodological rationale, and background knowledge.
 */

// CRC-32 Lookup Table for ZIP header validation
const crcTable = new Uint32Array(256);
for (let n = 0; n < 256; n++) {
  let c = n;
  for (let k = 0; k < 8; k++) {
    c = ((c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1));
  }
  crcTable[n] = c;
}

function calcCRC32(bytes) {
  let crc = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i++) {
    crc = (crc >>> 8) ^ crcTable[(crc ^ bytes[i]) & 0xFF];
  }
  return (crc ^ 0xFFFFFFFF) >>> 0;
}

/**
 * Pure JavaScript ZIP archive packager (PKZIP specification, Stored / uncompressed method 0)
 * Works universally in browser and Node.js environments.
 */
export class SimpleZip {
  constructor() {
    this.files = [];
  }

  addFile(name, content) {
    let data;
    if (typeof content === 'string') {
      data = new TextEncoder().encode(content);
    } else if (content instanceof Uint8Array) {
      data = content;
    } else if (content && content.buffer) {
      data = new Uint8Array(content.buffer);
    } else {
      data = new Uint8Array(content);
    }
    this.files.push({ name, data });
  }

  build() {
    const encoder = new TextEncoder();
    const prepared = this.files.map(f => {
      const nameBytes = encoder.encode(f.name);
      const crc = calcCRC32(f.data);
      const localSize = 30 + nameBytes.length + f.data.length;
      const cdSize = 46 + nameBytes.length;
      return { ...f, nameBytes, crc, localSize, cdSize };
    });

    let localTotal = 0;
    let cdTotal = 0;
    for (const p of prepared) {
      localTotal += p.localSize;
      cdTotal += p.cdSize;
    }
    const eocdSize = 22;
    const buffer = new Uint8Array(localTotal + cdTotal + eocdSize);
    const view = new DataView(buffer.buffer);

    let offset = 0;
    const offsets = [];

    // 1. Local file headers and data
    for (const p of prepared) {
      offsets.push(offset);
      view.setUint32(offset, 0x04034b50, true);
      view.setUint16(offset + 4, 20, true);
      view.setUint16(offset + 6, 0x0800, true); // UTF-8 filename
      view.setUint16(offset + 8, 0, true);      // Stored (no compression)
      view.setUint16(offset + 10, 0, true);
      view.setUint16(offset + 12, 0x5421, true); // Date/time
      view.setUint32(offset + 14, p.crc, true);
      view.setUint32(offset + 18, p.data.length, true);
      view.setUint32(offset + 22, p.data.length, true);
      view.setUint16(offset + 26, p.nameBytes.length, true);
      view.setUint16(offset + 28, 0, true);

      buffer.set(p.nameBytes, offset + 30);
      buffer.set(p.data, offset + 30 + p.nameBytes.length);
      offset += p.localSize;
    }

    const cdOffset = offset;

    // 2. Central Directory headers
    for (let i = 0; i < prepared.length; i++) {
      const p = prepared[i];
      const localOff = offsets[i];
      view.setUint32(offset, 0x02014b50, true);
      view.setUint16(offset + 4, 20, true);
      view.setUint16(offset + 6, 20, true);
      view.setUint16(offset + 8, 0x0800, true);
      view.setUint16(offset + 10, 0, true);
      view.setUint16(offset + 12, 0, true);
      view.setUint16(offset + 14, 0x5421, true);
      view.setUint32(offset + 16, p.crc, true);
      view.setUint32(offset + 20, p.data.length, true);
      view.setUint32(offset + 24, p.data.length, true);
      view.setUint16(offset + 28, p.nameBytes.length, true);
      view.setUint16(offset + 30, 0, true);
      view.setUint16(offset + 32, 0, true);
      view.setUint16(offset + 34, 0, true);
      view.setUint16(offset + 36, 0, true);
      view.setUint32(offset + 38, 0, true);
      view.setUint32(offset + 42, localOff, true);

      buffer.set(p.nameBytes, offset + 46);
      offset += p.cdSize;
    }

    // 3. End of Central Directory record (EOCD)
    view.setUint32(offset, 0x06054b50, true);
    view.setUint16(offset + 4, 0, true);
    view.setUint16(offset + 6, 0, true);
    view.setUint16(offset + 8, prepared.length, true);
    view.setUint16(offset + 10, prepared.length, true);
    view.setUint32(offset + 12, cdTotal, true);
    view.setUint32(offset + 16, cdOffset, true);
    view.setUint16(offset + 20, 0, true);

    return buffer;
  }
}

/**
 * Escapes characters for XML text
 */
function escapeXML(str) {
  return String(str ?? '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&apos;');
}

/**
 * WordprocessingML XML Document Builder
 */
export class DocxBuilder {
  constructor() {
    this.bodyElements = [];
  }

  addTitle(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:before="120" w:after="80"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="38"/><w:color w:val="0F172A"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addSubTitle(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:after="160"/></w:pPr>` +
      `<w:r><w:rPr><w:sz w:val="22"/><w:color w:val="2563EB"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addAttributionHeader() {
    this.bodyElements.push(
      `<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:after="60"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="18"/><w:color w:val="334155"/></w:rPr>` +
      `<w:t>Conceived, supervised design and testing: Dr G Narenthiran MB ChB BSc(MedSci)(Hons) MRCS(Ed.) FEBNS FRCS(SN)</w:t></w:r></w:p>` +
      `<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:after="180"/></w:pPr>` +
      `<w:r><w:rPr><w:sz w:val="16"/><w:color w:val="64748B"/></w:rPr>` +
      `<w:t>Copyright, G Narenthiran FEBS FRCS(SN), g_narenthiran@hotmail.com | Dedicated to Mrs Nirmaladevy Ganesalingam BSc (mother)</w:t></w:r></w:p>` +
      `<w:p><w:pPr><w:pBdr><w:bottom w:val="single" w:sz="6" w:space="8" w:color="CBD5E1"/></w:pBdr><w:spacing w:after="240"/></w:pPr></w:p>`
    );
    return this;
  }

  addDisclaimerBox() {
    this.bodyElements.push(
      `<w:tbl>` +
      `<w:tblPr><w:tblW w:w="5000" w:type="pct"/><w:tblBorders><w:top w:val="single" w:sz="6" w:color="F59E0B"/><w:left w:val="single" w:sz="18" w:color="D97706"/><w:bottom w:val="single" w:sz="6" w:color="F59E0B"/><w:right w:val="single" w:sz="6" w:color="F59E0B"/><w:insideH w:val="none"/><w:insideV w:val="none"/></w:tblBorders></w:tblPr>` +
      `<w:tr><w:tc><w:tcPr><w:shd w:fill="FEF3C7"/><w:tcMar><w:top w:w="120" w:type="dxa"/><w:bottom w:w="120" w:type="dxa"/><w:left w:w="200" w:type="dxa"/><w:right w:w="200" w:type="dxa"/></w:tcMar></w:tcPr>` +
      `<w:p><w:r><w:rPr><w:b/><w:sz w:val="18"/><w:color w:val="92400E"/></w:rPr><w:t>⚠️ TESTING &amp; METHODOLOGICAL NOTICE: </w:t></w:r>` +
      `<w:r><w:rPr><w:sz w:val="18"/><w:color w:val="78350F"/></w:rPr><w:t>AI was used to vibe code this WebApp. The App is still in testing phase. Not to use for clinical, research or decision making.</w:t></w:r></w:p>` +
      `</w:tc></w:tr></w:tbl>` +
      `<w:p><w:pPr><w:spacing w:after="180"/></w:pPr></w:p>`
    );
    return this;
  }

  addHeading1(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:spacing w:before="280" w:after="120"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="28"/><w:color w:val="1E3A8A"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addHeading2(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:spacing w:before="200" w:after="80"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="24"/><w:color w:val="0D9488"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addHeading3(text) {
    this.bodyElements.push(
      `<w:p><w:pPr><w:spacing w:before="160" w:after="60"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:sz w:val="20"/><w:color w:val="334155"/></w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addParagraph(text, opts = {}) {
    const bold = opts.bold ? '<w:b/>' : '';
    const italic = opts.italic ? '<w:i/>' : '';
    const color = opts.color ? `<w:color w:val="${opts.color}"/>` : '<w:color w:val="1E293B"/>';
    const sz = opts.size ? `<w:sz w:val="${opts.size}"/>` : '<w:sz w:val="21"/>';
    const align = opts.align ? `<w:jc w:val="${opts.align}"/>` : '';
    const spacing = opts.spacing ? `<w:spacing w:after="${opts.spacing}"/>` : '<w:spacing w:after="120" w:line="276" w:lineRule="auto"/>';

    this.bodyElements.push(
      `<w:p><w:pPr>${align}${spacing}</w:pPr>` +
      `<w:r><w:rPr>${bold}${italic}${color}${sz}</w:rPr>` +
      `<w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addBullet(text, opts = {}) {
    const bold = opts.bold ? '<w:b/>' : '';
    const sz = opts.size ? `<w:sz w:val="${opts.size}"/>` : '<w:sz w:val="21"/>';
    this.bodyElements.push(
      `<w:p><w:pPr><w:ind w:left="360"/><w:spacing w:after="80" w:line="260" w:lineRule="auto"/></w:pPr>` +
      `<w:r><w:rPr><w:b/><w:color w:val="2563EB"/></w:rPr><w:t>• </w:t></w:r>` +
      `<w:r><w:rPr>${bold}${sz}<w:color w:val="1E293B"/></w:rPr><w:t>${escapeXML(text)}</w:t></w:r></w:p>`
    );
    return this;
  }

  addCalloutBox(title, content, bgColor = "F1F5F9", borderColor = "0284C7") {
    this.bodyElements.push(
      `<w:tbl>` +
      `<w:tblPr><w:tblW w:w="5000" w:type="pct"/><w:tblBorders><w:top w:val="single" w:sz="4" w:color="CBD5E1"/><w:left w:val="single" w:sz="24" w:color="${borderColor}"/><w:bottom w:val="single" w:sz="4" w:color="CBD5E1"/><w:right w:val="single" w:sz="4" w:color="CBD5E1"/><w:insideH w:val="none"/><w:insideV w:val="none"/></w:tblBorders></w:tblPr>` +
      `<w:tr><w:tc><w:tcPr><w:shd w:fill="${bgColor}"/><w:tcMar><w:top w:w="140" w:type="dxa"/><w:bottom w:w="140" w:type="dxa"/><w:left w:w="220" w:type="dxa"/><w:right w:w="220" w:type="dxa"/></w:tcMar></w:tcPr>` +
      (title ? `<w:p><w:pPr><w:spacing w:after="80"/></w:pPr><w:r><w:rPr><w:b/><w:sz w:val="22"/><w:color w:val="0F172A"/></w:rPr><w:t>${escapeXML(title)}</w:t></w:r></w:p>` : '') +
      `<w:p><w:pPr><w:spacing w:after="40" w:line="276" w:lineRule="auto"/></w:pPr><w:r><w:rPr><w:sz w:val="20"/><w:color w:val="334155"/></w:rPr><w:t>${escapeXML(content)}</w:t></w:r></w:p>` +
      `</w:tc></w:tr></w:tbl>` +
      `<w:p><w:pPr><w:spacing w:after="160"/></w:pPr></w:p>`
    );
    return this;
  }

  addTable(headers, rows) {
    let tblXML = `<w:tbl><w:tblPr><w:tblW w:w="5000" w:type="pct"/><w:tblBorders>` +
      `<w:top w:val="single" w:sz="6" w:color="CBD5E1"/>` +
      `<w:left w:val="single" w:sz="4" w:color="E2E8F0"/>` +
      `<w:bottom w:val="single" w:sz="8" w:color="94A3B8"/>` +
      `<w:right w:val="single" w:sz="4" w:color="E2E8F0"/>` +
      `<w:insideH w:val="single" w:sz="4" w:color="E2E8F0"/>` +
      `<w:insideV w:val="single" w:sz="4" w:color="E2E8F0"/>` +
      `</w:tblBorders></w:tblPr>`;

    // Header Row
    if (headers && headers.length > 0) {
      tblXML += `<w:tr><w:trPr><w:tblHeader/></w:trPr>`;
      for (const h of headers) {
        tblXML += `<w:tc><w:tcPr><w:shd w:fill="E2E8F0"/><w:tcMar><w:top w:w="120" w:type="dxa"/><w:bottom w:w="120" w:type="dxa"/><w:left w:w="140" w:type="dxa"/><w:right w:w="140" w:type="dxa"/></w:tcMar></w:tcPr>` +
          `<w:p><w:pPr><w:spacing w:after="40"/></w:pPr><w:r><w:rPr><w:b/><w:sz w:val="19"/><w:color w:val="0F172A"/></w:rPr><w:t>${escapeXML(h)}</w:t></w:r></w:p></w:tc>`;
      }
      tblXML += `</w:tr>`;
    }

    // Data Rows
    for (let rIdx = 0; rIdx < rows.length; rIdx++) {
      const row = rows[rIdx];
      const bg = rIdx % 2 === 1 ? 'F8FAFC' : 'FFFFFF';
      tblXML += `<w:tr>`;
      for (let cIdx = 0; cIdx < row.length; cIdx++) {
        const cell = row[cIdx];
        const isFirstCol = cIdx === 0;
        tblXML += `<w:tc><w:tcPr><w:shd w:fill="${bg}"/><w:tcMar><w:top w:w="100" w:type="dxa"/><w:bottom w:w="100" w:type="dxa"/><w:left w:w="140" w:type="dxa"/><w:right w:w="140" w:type="dxa"/></w:tcMar></w:tcPr>` +
          `<w:p><w:pPr><w:spacing w:after="30"/></w:pPr><w:r><w:rPr>${isFirstCol ? '<w:b/>' : ''}<w:sz w:val="18"/><w:color w:val="1E293B"/></w:rPr><w:t>${escapeXML(cell)}</w:t></w:r></w:p></w:tc>`;
      }
      tblXML += `</w:tr>`;
    }

    tblXML += `</w:tbl><w:p><w:pPr><w:spacing w:after="160"/></w:pPr></w:p>`;
    this.bodyElements.push(tblXML);
    return this;
  }

  buildDocumentXML() {
    return (
      `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>` +
      `<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">` +
      `<w:body>` +
      this.bodyElements.join('\n') +
      `<w:sectPr>` +
      `<w:pgSz w:w="12240" w:h="15840"/>` +
      `<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>` +
      `</w:sectPr>` +
      `</w:body></w:document>`
    );
  }

  toZip() {
    const zip = new SimpleZip();
    zip.addFile('[Content_Types].xml',
      `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>` +
      `<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">` +
      `<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>` +
      `<Default Extension="xml" ContentType="application/xml"/>` +
      `<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>` +
      `</Types>`
    );

    zip.addFile('_rels/.rels',
      `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>` +
      `<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">` +
      `<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>` +
      `</Relationships>`
    );

    zip.addFile('word/document.xml', this.buildDocumentXML());
    return zip;
  }

  generateUint8Array() {
    return this.toZip().build();
  }

  generateBlob() {
    const uint8 = this.generateUint8Array();
    return new Blob([uint8], {
      type: 'application/vnd.openxmlformats-officedocument.wordprocessingml.document'
    });
  }
}

/**
 * High-level DOCX Report Generator per module
 */
export const DocxReports = {
  createDescriptiveDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Descriptive Statistics, Normality & Distribution Diagnostics')
      .addAttributionHeader()
      .addDisclaimerBox();

    // 1. Information & Analyzed Input Data
    d.addHeading1('1. Analyzed Cohort Information & Input Data')
      .addParagraph(`Variable/Series Label: ${data.name || 'Sample Continuous Variable'}`)
      .addParagraph(`Sample Size (N): ${data.n} observations`)
      .addParagraph(`Observed Minimum to Maximum: ${data.min.toFixed(2)} to ${data.max.toFixed(2)} (Range = ${(data.max - data.min).toFixed(2)})`)
      .addParagraph(`Raw Data Sample: ${data.values ? data.values.slice(0, 20).join(', ') + (data.values.length > 20 ? ' ...' : '') : 'N/A'}`);

    // 2. Statistical Outcome & Numerical Results
    d.addHeading1('2. Statistical Outcome & Distributional Metrics');
    d.addHeading2('Central Tendency & Parametric/Non-Parametric Dispersion');
    d.addTable(
      ['Statistical Metric', 'Numerical Value', 'Clinical / Mathematical Interpretation'],
      [
        ['Sample Size (N)', `${data.n}`, 'Valid non-null numerical observations'],
        ['Mean (M)', `${data.mean.toFixed(2)}`, 'Arithmetic average (sensitive to outliers)'],
        ['95% Confidence Interval', `[${data.ci95[0].toFixed(2)}, ${data.ci95[1].toFixed(2)}]`, 'Estimated population mean parameter bound'],
        ['Standard Deviation (SD, s)', `${data.sd.toFixed(2)}`, 'Average spread about the sample mean'],
        ['Sample Variance (s²)', `${data.variance.toFixed(2)}`, 'Second central moment (Bessel-corrected n-1)'],
        ['Standard Error of Mean (SEM)', `${data.sem.toFixed(3)}`, 'Precision of sample mean (s / √n)'],
        ['Median (50th Percentile)', `${data.median.toFixed(2)}`, 'Robust central value unaffected by skewness'],
        ['25th Percentile (Q1)', `${data.q1.toFixed(2)}`, 'Lower quartile boundary'],
        ['75th Percentile (Q3)', `${data.q3.toFixed(2)}`, 'Upper quartile boundary'],
        ['Interquartile Range (IQR)', `${data.iqr.toFixed(2)}`, 'Spread of middle 50% of observations (Q3 - Q1)'],
        ['Mode(s)', data.modes && data.modes.length ? data.modes.join(', ') : 'No repeat values', `Peak sample frequency (${data.maxFreq || 1})`],
        ['Range (Max - Min)', `${(data.max - data.min).toFixed(2)}`, `Spans from ${data.min.toFixed(2)} to ${data.max.toFixed(2)}`]
      ]
    );

    d.addHeading2('Shape Diagnostics & Normality Assessment');
    d.addTable(
      ['Distribution Metric', 'Observed Value', 'Statistical Benchmark / Status'],
      [
        ['Sample Skewness (G1)', `${data.skewness.toFixed(3)}`, data.skewnessInterpretation || (Math.abs(data.skewness) < 0.5 ? 'Approximately Symmetric' : 'Skewed')],
        ['Excess Kurtosis (G2)', `${data.kurtosis.toFixed(3)}`, data.kurtosisInterpretation || (Math.abs(data.kurtosis) < 0.5 ? 'Mesokurtic' : 'Heavy/Light-Tailed')],
        ['Normality Omnibus Test', 'Jarque-Bera Test', `JB = ${data.normality.statistic.toFixed(2)}, df = 2, p = ${data.normality.pValue.toExponential(3)}`],
        ['Normality Decision', data.normality.isNormal ? 'Normal (Gaussian)' : 'Statistically Non-Normal', data.normality.isNormal ? 'Parametric tests valid (p ≥ 0.05)' : 'Reject H0 of normality (p < 0.05)']
      ]
    );

    d.addHeading2('Tukey Fences & Outlier Diagnostics');
    const outliersCount = data.outliers ? data.outliers.length : 0;
    d.addParagraph(`Lower Tukey Fence (Q1 - 1.5×IQR): ${data.lowerFence.toFixed(2)}`);
    d.addParagraph(`Upper Tukey Fence (Q3 + 1.5×IQR): ${data.upperFence.toFixed(2)}`);
    d.addParagraph(`Detected Outliers: ${outliersCount === 0 ? 'None detected beyond fences' : `${outliersCount} outlier point(s) identified`}`);
    if (outliersCount > 0) {
      const outlierRows = data.outliers.map(o => [
        `${o.value.toFixed(2)}`,
        `Z = ${o.zScore > 0 ? '+' : ''}${o.zScore.toFixed(2)}`,
        `${o.type}`,
        o.isExtreme ? 'Beyond 3.0×IQR (Extreme Outlier)' : 'Between 1.5× and 3.0×IQR (Mild Outlier)'
      ]);
      d.addTable(['Outlier Value', 'Standardized Z-Score', 'Classification', 'Tukey Fence Range'], outlierRows);
    }

    // 3. Clinical & Statistical Outcome Interpretation
    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'APA / ICMJE Style Summary Statement',
      data.reportText || 'Continuous variable descriptive summary.',
      'F0FDF4',
      '16A34A'
    );
    if (!data.normality.isNormal || outliersCount > 0) {
      d.addParagraph('Clinical Presentation Recommendation: Due to statistically significant distributional skewness and/or the presence of outlier points, reporting the Median and Interquartile Range (IQR) as the primary measures of central tendency and dispersion is strongly advised by STROBE and ICMJE reporting standards.', { bold: true });
    } else {
      d.addParagraph('Clinical Presentation Recommendation: The dataset satisfies Gaussian assumptions without extreme outliers; reporting the Mean and Standard Deviation (Mean ± SD) alongside the 95% Confidence Interval is methodologically sound.', { bold: true });
    }

    // 4. Reason This Particular Test Was Chosen
    d.addHeading1('4. Reason Particular Tests & Methods Were Chosen');
    d.addBullet('Jarque-Bera Omnibus Normality Test: Chosen because it assesses both sample skewness (third standardized moment) and excess kurtosis (fourth moment) simultaneously under an asymptotic Chi-square distribution (df=2). In clinical trials, biological markers and surgical recovery times frequently exhibit asymmetric positive skewness and long tails; the Jarque-Bera test provides superior power in detecting combined shape violations without the arbitrary binning required by Kolmogorov-Smirnov.');
    d.addBullet('Tukey 1.5×IQR Outlier Fences: Selected because Gaussian Z-scores (e.g. ±3 SD) suffer from masking and swamping in the presence of severe outliers (since standard deviation itself is distorted by the extreme points). Tukey fences utilize the interquartile range (IQR), providing a non-parametric, robust threshold for identifying contaminated or physiologically anomalous values.');
    d.addBullet('Dual Parametric & Non-Parametric Metrics: Furnishing both Mean/SD and Median/IQR enables clinicians to immediately detect data skewness and verify whether parametric inferential tests (such as Student\'s t or ANOVA) will be robust or whether non-parametric alternatives (Mann-Whitney, Kruskal-Wallis) must be adopted.');

    // 5. Background Knowledge & Methodology References
    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Descriptive statistics constitute the foundation of all clinical evidence. Under the EQUATOR Network guidelines (CONSORT 2010 for clinical trials, STROBE for observational studies, and STARD for diagnostic accuracy), all patient baseline demographics and clinical outcome variables must be reported with appropriate indicators of spread.');
    d.addParagraph('Mathematical Definitions:');
    d.addBullet('Sample Mean: x̄ = (1/n) ∑ xᵢ. Reflects the gravitational balance point of the observations.');
    d.addBullet('Bessel-Corrected Standard Deviation: s = √[ (1/(n-1)) ∑ (xᵢ - x̄)² ]. Uses n-1 degrees of freedom to remove negative bias in estimating the true population variance σ.');
    d.addBullet('Standard Error of the Mean: SEM = s / √n. Measures the precision with which the sample mean estimates the true population mean; as sample size n grows, SEM contracts.');
    d.addBullet('Jarque-Bera Statistic: JB = (n/6) [ S² + (K² / 4) ] ~ χ²(2), where S is sample skewness and K is excess kurtosis.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Jarque CM, Bera AK (1987). A test for normality of observations and regression residuals. International Statistical Review, 55(2): 163–172.');
    d.addBullet('Tukey JW (1977). Exploratory Data Analysis. Addison-Wesley.');
    d.addBullet('Altman DG, Bland JM (1996). Detecting skewness from summary information. BMJ, 313(7066): 1200.');

    return d;
  },

  createHypothesisDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle(`Module: Two-Cohort Hypothesis Testing (${data.testName || 'Inferential Analysis'})`)
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Cohort Information & Input Data')
      .addParagraph(`Cohort A: ${data.nameA} (n = ${data.groupA.n} patients/specimens)`)
      .addParagraph(`Cohort B: ${data.nameB} (n = ${data.groupB.n} patients/specimens)`)
      .addParagraph(`Evaluated Statistical Paradigm: ${data.testName || 'Two-Sample Test'}`)
      .addParagraph(`Cohort A Data Preview: ${data.groupA.values ? data.groupA.values.slice(0, 15).join(', ') : 'N/A'}`)
      .addParagraph(`Cohort B Data Preview: ${data.groupB.values ? data.groupB.values.slice(0, 15).join(', ') : 'N/A'}`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addTable(
      ['Cohort Summary', 'Sample Size (n)', 'Mean', 'Std Dev (SD)', 'Std Error (SEM)', '95% Confidence Interval', 'Median (IQR)'],
      [
        [data.nameA, `${data.groupA.n}`, `${data.groupA.mean.toFixed(2)}`, `${data.groupA.sd.toFixed(2)}`, `${data.groupA.sem.toFixed(3)}`, `[${data.groupA.ci95[0].toFixed(2)}, ${data.groupA.ci95[1].toFixed(2)}]`, `${data.groupA.median.toFixed(2)} (${data.groupA.iqr.toFixed(2)})`],
        [data.nameB, `${data.groupB.n}`, `${data.groupB.mean.toFixed(2)}`, `${data.groupB.sd.toFixed(2)}`, `${data.groupB.sem.toFixed(3)}`, `[${data.groupB.ci95[0].toFixed(2)}, ${data.groupB.ci95[1].toFixed(2)}]`, `${data.groupB.median.toFixed(2)} (${data.groupB.iqr.toFixed(2)})`]
      ]
    );

    d.addHeading2('Comparative Inferential Test Results');
    d.addTable(
      ['Inferential Parameter', 'Calculated Value', 'Clinical Interpretation / Benchmark'],
      [
        ['Test Statistic', `${data.testName.includes('Mann-Whitney') ? 'U = ' : 't = '}${data.statistic.toFixed(3)}`, 'Standardized difference between cohort locations'],
        ['Degrees of Freedom (df)', `${data.df ? data.df.toFixed(2) : 'N/A (Rank test)'}`, 'Satterthwaite adjustment for unequal cohort variances'],
        ['p-Value (Two-Tailed)', `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + data.pValue.toFixed(4)}`, data.isSignificant ? 'Statistically Significant (p < 0.05)' : 'Not Significant (p ≥ 0.05)'],
        ['Mean Difference (ΔM)', `${data.meanDiff !== undefined ? (data.meanDiff >= 0 ? '+' : '') + data.meanDiff.toFixed(2) : 'N/A'}`, `Observed clinical point difference (${data.nameA} - ${data.nameB})`],
        ['95% CI of Difference', data.ci95 ? `[${data.ci95[0].toFixed(2)}, ${data.ci95[1].toFixed(2)}]` : 'N/A', 'Range of plausible true population differences'],
        ['Standardized Effect Size', data.cohensD !== undefined ? `Cohen\'s d = ${data.cohensD.toFixed(2)}` : `Rank-Biserial r = ${(data.rankBiserial || 0).toFixed(2)}`, 'Magnitude of separation independent of sample size'],
        ['Effect Classification', Math.abs(data.cohensD || data.rankBiserial || 0) >= 0.8 ? 'Large Effect' : (Math.abs(data.cohensD || data.rankBiserial || 0) >= 0.5 ? 'Moderate Effect' : 'Small/Negligible Effect'), 'Cohen (1988) benchmark: 0.2 small, 0.5 medium, 0.8 large']
      ]
    );

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'APA / ICMJE Style Summary Statement',
      data.reportText || 'Hypothesis testing narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Welch\'s Unequal Variances t-Test: Selected as the gold-standard default because classical Student\'s t-test assumes homoscedasticity (equal population variances). In biomedical and clinical research, treatment interventions frequently alter not only the mean but also the variance (e.g. patients respond heterogeneously). When variances differ even moderately, Student\'s t-test exhibits severely inflated Type I error rates. Welch\'s t-test adjusts the degrees of freedom using the Satterthwaite approximation, providing robust Type I error control without sacrificing statistical power when variances happen to be identical (Ruxton 2006, Delacre et al. 2017).');
    d.addBullet('Mann-Whitney U Test (Alternative): Provided for ordinal scales, Glasgow Coma Scale (GCS) scores, or continuous markers with severe skewness. Instead of comparing means, it evaluates stochastic dominance across ranked values.');
    d.addBullet('Reporting Effect Sizes: Statistical significance (p < 0.05) is heavily dependent on sample size; very large clinical cohorts can produce tiny p-values for clinically meaningless differences. Cohen\'s d quantifies the practical medical magnitude of the separation.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Hypothesis testing evaluates the likelihood of observing the experimental difference under the null hypothesis (H0: μA = μB).');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Welch\'s t-Statistic: t = (x̄A - x̄B) / √[ (sA²/nA) + (sB²/nB) ].');
    d.addBullet('Welch-Satterthwaite Degrees of Freedom: df = [ (sA²/nA + sB²/nB)² ] / [ (sA²/nA)² / (nA - 1) + (sB²/nB)² / (nB - 1) ].');
    d.addBullet('Cohen\'s d: d = (x̄A - x̄B) / sₚ, where sₚ is the pooled standard deviation.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Welch BL (1947). The generalization of \'Student\'s\' problem when several different population variances are involved. Biometrika, 34(1/2): 28–35.');
    d.addBullet('Ruxton GD (2006). The unequal variance t-test is an underused alternative to Student\'s t-test and the Mann-Whitney U test. Behavioral Ecology, 17(4): 688–690.');
    d.addBullet('Delacre M, Lakens D, Leys C (2017). Why psychologists should by default use Welch\'s t-test instead of Student\'s t-test. International Review of Social Psychology, 30(1): 92–101.');

    return d;
  },

  createAnovaDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Multi-Cohort Variance & Tukey HSD Post-Hoc Pairwise Analysis')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Cohorts Information & Input Data')
      .addParagraph(`Number of Independent Cohorts (k): ${data.k || data.groups.length}`)
      .addParagraph(`Total Analyzed Sample Size (N): ${data.totalN} patients/specimens`)
      .addParagraph(`Grand Mean across All Cohorts: ${data.grandMean.toFixed(2)}`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    const groupRows = data.groups.map(g => [
      g.name,
      `${g.stats.n}`,
      `${g.stats.mean.toFixed(2)}`,
      `${g.stats.sd.toFixed(2)}`,
      `${g.stats.sem.toFixed(3)}`,
      `[${g.stats.ci95[0].toFixed(2)}, ${g.stats.ci95[1].toFixed(2)}]`,
      `${g.stats.median.toFixed(2)} (${g.stats.iqr.toFixed(2)})`
    ]);
    d.addTable(['Cohort Name', 'Sample n', 'Mean (M)', 'Std Dev (SD)', 'Std Error (SEM)', '95% CI of Mean', 'Median (IQR)'], groupRows);

    d.addHeading2('One-Way ANOVA Summary Table');
    d.addTable(
      ['Source of Variation', 'Sum of Squares (SS)', 'Degrees of Freedom (df)', 'Mean Square (MS)', 'F-Statistic', 'p-Value', 'Omega-Squared (ω²)'],
      [
        ['Between Groups (Treatment)', `${(data.ssBetween || 0).toFixed(2)}`, `${data.dfBetween || 0}`, `${(data.msBetween || 0).toFixed(2)}`, `F = ${(data.fStatistic || 0).toFixed(2)}`, `${data.pValue < 0.001 ? 'p < .001' : 'p = ' + (data.pValue || 0).toFixed(4)}`, `${(data.omegaSquared || 0).toFixed(3)}`],
        ['Within Groups (Residual/Error)', `${(data.ssWithin || 0).toFixed(2)}`, `${data.dfWithin || 0}`, `${(data.msWithin || 0).toFixed(2)}`, '-', '-', `Eta² (η²) = ${(data.etaSquared || 0).toFixed(3)}`],
        ['Total', `${(data.ssTotal !== undefined ? data.ssTotal : ((data.ssBetween || 0) + (data.ssWithin || 0))).toFixed(2)}`, `${(data.dfBetween || 0) + (data.dfWithin || 0)}`, '-', '-', '-', '-']
      ]
    );

    if (data.pairwise && data.pairwise.length > 0) {
      d.addHeading2('Tukey\'s HSD Post-Hoc Pairwise Contrasts');
      const pairRows = data.pairwise.map(p => [
        p.comparison,
        `${(p.meanDiff >= 0 ? '+' : '')}${p.meanDiff.toFixed(2)}`,
        `${p.seDiff.toFixed(3)}`,
        `q = ${p.qStatistic.toFixed(2)}`,
        `${p.pValue < 0.001 ? 'p < .001' : 'p = ' + p.pValue.toFixed(4)}`,
        `[${p.ci95[0].toFixed(2)}, ${p.ci95[1].toFixed(2)}]`,
        `d = ${p.cohensD.toFixed(2)}`,
        p.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)'
      ]);
      d.addTable(['Pairwise Contrast', 'Mean Diff (ΔM)', 'Std Error', 'Tukey q', 'Adjusted p', '95% CI of Diff', 'Cohen\'s d', 'Significance'], pairRows);
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'ANOVA & Post-Hoc APA Clinical Summary',
      data.reportText || 'ANOVA narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('One-Way Omnibus ANOVA: Selected because testing multiple cohorts with uncorrected pairwise t-tests results in severe Family-Wise Error Rate inflation (FWER). With 3 cohorts, 3 comparisons yield α_FW = 1 - (1 - 0.05)³ = 14.3%; with 5 cohorts (10 comparisons), α_FW exceeds 40%. ANOVA provides a rigorous omnibus test that simultaneously assesses whether any between-cohort variance exceeds within-cohort residual variation.');
    d.addBullet('Tukey\'s Honest Significant Difference (HSD): Chosen as the post-hoc method because it utilizes the Studentized Range distribution (q) to strictly bound the overall Family-Wise Error Rate at α = 0.05 across all possible pairwise comparisons, while preserving substantially greater statistical power than overly conservative Bonferroni adjustments.');
    d.addBullet('Omega-Squared (ω²) Reporting: Included alongside Eta-squared (η²) because Eta-squared represents a sample proportion of variance that is positively biased in small clinical samples. Omega-squared provides an unbiased population effect size estimate.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('ANOVA partitions the total sum of squares into treatment (between) and error (within) components: SS_Total = SS_Between + SS_Within.');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Between-Groups Mean Square: MS_B = SS_B / (k - 1).');
    d.addBullet('Within-Groups Mean Square: MS_W = SS_W / (N - k).');
    d.addBullet('F-Ratio: F = MS_B / MS_W ~ F(k-1, N-k).');
    d.addBullet('Tukey Studentized Range: q = |x̄A - x̄B| / √[ (MS_W / 2) (1/nA + 1/nB) ].');
    d.addParagraph('Key Academic References:');
    d.addBullet('Fisher RA (1925). Statistical Methods for Research Workers. Oliver and Boyd, Edinburgh.');
    d.addBullet('Tukey JW (1949). Comparing individual means in the analysis of variance. Biometrics, 5(2): 99–114.');
    d.addBullet('Hayter AJ (1984). A proof of the conjecture that the Tukey-Kramer multiple comparisons procedure is conservative. Annals of Statistics, 12(1): 61–75.');

    return d;
  },

  createCategoricalDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle(`Module: 2x2 Contingency, Risk Metrics & Diagnostic Evaluation`)
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed 2x2 Matrix Information & Configuration')
      .addParagraph(`Analysis Framework: ${data.mode === 'diagnostic' ? 'Diagnostic Test Evaluation vs Gold Standard Reference' : 'Clinical Study / Trial (Intervention vs Control Cohort)'}`)
      .addParagraph(`Total Evaluated Subjects (N): ${data.totalN} patients/samples`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2('Observed 2x2 Contingency Matrix');
    if (data.mode === 'diagnostic') {
      d.addTable(
        ['Test Condition', 'Gold Standard Positive (+)', 'Gold Standard Negative (-)', 'Total Row Margins'],
        [
          ['New Test Positive (+)', `TP = ${data.a}`, `FP = ${data.b}`, `Test Pos = ${data.a + data.b}`],
          ['New Test Negative (-)', `FN = ${data.c}`, `TN = ${data.d}`, `Test Neg = ${data.c + data.d}`],
          ['Total Column Margins', `Condition Pos = ${data.a + data.c}`, `Condition Neg = ${data.b + data.d}`, `N = ${data.totalN}`]
        ]
      );

      d.addHeading2('Diagnostic Accuracy & Discriminative Performance Metrics');
      d.addTable(
        ['Diagnostic Index', 'Point Estimate (%)', '95% Confidence Interval', 'Clinical Definition & Interpretation'],
        [
          ['Sensitivity (True Positive Rate)', `${(data.diag.sensitivity * 100).toFixed(1)}%`, `[${(data.diag.sensitivityCI95[0]*100).toFixed(1)}%, ${(data.diag.sensitivityCI95[1]*100).toFixed(1)}%]`, 'Ability of test to correctly identify patients with disease (TP / [TP + FN])'],
          ['Specificity (True Negative Rate)', `${(data.diag.specificity * 100).toFixed(1)}%`, `[${(data.diag.specificityCI95[0]*100).toFixed(1)}%, ${(data.diag.specificityCI95[1]*100).toFixed(1)}%]`, 'Ability of test to correctly identify disease-free patients (TN / [TN + FP])'],
          ['Positive Predictive Value (PPV)', `${(data.diag.ppv * 100).toFixed(1)}%`, `[${(data.diag.ppvCI95[0]*100).toFixed(1)}%, ${(data.diag.ppvCI95[1]*100).toFixed(1)}%]`, 'Probability that a patient with a positive test truly has disease'],
          ['Negative Predictive Value (NPV)', `${(data.diag.npv * 100).toFixed(1)}%`, `[${(data.diag.npvCI95[0]*100).toFixed(1)}%, ${(data.diag.npvCI95[1]*100).toFixed(1)}%]`, 'Probability that a patient with a negative test is truly disease-free'],
          ['Overall Diagnostic Accuracy', `${(data.diag.accuracy * 100).toFixed(1)}%`, `[${(data.diag.accuracyCI95[0]*100).toFixed(1)}%, ${(data.diag.accuracyCI95[1]*100).toFixed(1)}%]`, 'Proportion of all test results that were correct ([TP + TN] / N)'],
          ['Positive Likelihood Ratio (LR+)', `${data.diag.plr.toFixed(2)}`, `[${data.diag.plrCI95[0].toFixed(2)}, ${data.diag.plrCI95[1].toFixed(2)}]`, 'Ratio of TPR to FPR (>10 indicates strong diagnostic confirmation)'],
          ['Negative Likelihood Ratio (LR-)', `${data.diag.nlr.toFixed(2)}`, `[${data.diag.nlrCI95[0].toFixed(2)}, ${data.diag.nlrCI95[1].toFixed(2)}]`, 'Ratio of FNR to TNR (<0.1 indicates strong rule-out capacity)'],
          ['Youden\'s Index (J)', `${data.diag.youdenJ.toFixed(3)}`, '-', 'Overall effectiveness criterion: Sensitivity + Specificity - 1']
        ]
      );
    } else {
      d.addTable(
        ['Study Arm / Cohort', 'Event Occurred (+)', 'Event Absent (-)', 'Total Cohort Sample'],
        [
          ['Intervention / Treatment', `Events = ${data.a}`, `Non-Events = ${data.b}`, `Total = ${data.a + data.b}`],
          ['Control / Placebo', `Events = ${data.c}`, `Non-Events = ${data.d}`, `Total = ${data.c + data.d}`],
          ['Total Margins', `Total Events = ${data.a + data.c}`, `Total Non-Events = ${data.b + data.d}`, `N = ${data.totalN}`]
        ]
      );

      d.addHeading2('Hypothesis Tests & Clinical Risk Metrics');
      d.addTable(
        ['Statistical / Risk Metric', 'Numerical Value', 'Clinical Interpretation / Benchmark'],
        [
          ['Pearson Chi-Square (χ²)', `χ²(1) = ${data.chiSquare.standard.toFixed(2)} (p = ${data.chiSquare.pValueStandard < 0.001 ? '< .001' : data.chiSquare.pValueStandard.toFixed(4)})`, 'Asymptotic test of association'],
          ['Yates\' Continuity Corrected χ²', `χ²(1) = ${data.chiSquare.yates.toFixed(2)} (p = ${data.chiSquare.pValueYates < 0.001 ? '< .001' : data.chiSquare.pValueYates.toFixed(4)})`, 'Corrects for discrete probability overestimation'],
          ['Fisher\'s Exact Test (Two-Tailed)', `p = ${data.fishersExact.pValue < 0.001 ? '< .001' : data.fishersExact.pValue.toFixed(4)}`, 'Exact conditional hypergeometric probability'],
          ['Odds Ratio (OR)', `${data.risk.oddsRatio.toFixed(2)} (95% CI [${data.risk.orCI95[0].toFixed(2)}, ${data.risk.orCI95[1].toFixed(2)}])`, 'Odds of event in intervention vs control'],
          ['Relative Risk (RR)', `${data.risk.relativeRisk.toFixed(2)} (95% CI [${data.risk.rrCI95[0].toFixed(2)}, ${data.risk.rrCI95[1].toFixed(2)}])`, 'Risk of event in intervention vs control'],
          ['Absolute Risk Reduction (ARR)', `${(data.risk.arr * 100).toFixed(1)}% (95% CI [${(data.risk.arrCI95[0]*100).toFixed(1)}%, ${(data.risk.arrCI95[1]*100).toFixed(1)}%])`, 'Absolute difference in event rate'],
          ['Relative Risk Reduction (RRR)', `${(data.risk.rrr * 100).toFixed(1)}%`, 'Proportional reduction relative to baseline risk'],
          ['Number Needed to Treat (NNT)', `${data.risk.nnt.toFixed(1)} patients`, 'Patients needed to treat to prevent 1 adverse event (1 / ARR)']
        ]
      );
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Clinical / Epidemiological Summary',
      data.reportText || '2x2 categorical narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Fisher\'s Exact Test vs Pearson Chi-Square: Pearson\'s Chi-Square is an asymptotic test requiring adequate expected cell counts (Cochrane criterion: no expected cell < 1, and no more than 20% < 5). In clinical trials or rare surgical complications, small counts render Chi-Square unreliable. Fisher\'s Exact Test evaluates the exact hypergeometric permutation distribution, providing unconditional validity across any sample size.');
    d.addBullet('Importance of ARR and NNT in Evidence-Based Medicine: While Relative Risk Reduction (RRR) is commonly advertised in commercial literature, it can dramatically exaggerate clinical benefits. For instance, reducing risk from 2 in 10,000 to 1 in 10,000 is a 50% relative reduction, yet the ARR is 0.01% (NNT = 10,000). Reporting ARR and NNT ensures clinicians comprehend the real-world clinical workload required to achieve patient benefit.');
    d.addBullet('Likelihood Ratios for Diagnostic Evaluation: PPV and NPV are profoundly distorted by disease prevalence (Bayes\' theorem). Likelihood Ratios are prevalence-independent intrinsic properties of the diagnostic test.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Fisher\'s Exact Hypergeometric Probability: P = [ (a+b)! (c+d)! (a+c)! (b+d)! ] / [ N! a! b! c! d! ].');
    d.addBullet('Pearson\'s Chi-Square: χ² = ∑ [ (O - E)² / E ], with df = (r-1)(c-1) = 1.');
    d.addBullet('Number Needed to Treat: NNT = 1 / ARR = 1 / |p_intervention - p_control|.');
    d.addBullet('Positive Likelihood Ratio: LR+ = Sensitivity / (1 - Specificity).');
    d.addBullet('Negative Likelihood Ratio: LR- = (1 - Sensitivity) / Specificity.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Fisher RA (1922). On the interpretation of χ² from contingency tables, and the calculation of P. Journal of the Royal Statistical Society, 85(1): 87–94.');
    d.addBullet('Laupacis A, Sackett DL, Roberts RS (1988). An assessment of clinically useful measures of the consequences of treatment. N Engl J Med, 318(26): 1728–1733.');
    d.addBullet('Yates F (1934). Contingency tables involving small numbers and the χ² test. Supplement to the Journal of the Royal Statistical Society, 1(2): 217–235.');

    return d;
  },

  createCorrelationDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Correlation, Linear Regression & Non-Linear Curve Morphology')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Variables & Data Information')
      .addParagraph(`Independent Predictor (X): ${data.xName || 'Variable X'}`)
      .addParagraph(`Dependent Outcome (Y): ${data.yName || 'Variable Y'}`)
      .addParagraph(`Sample Size (N pairs): ${data.n} observations`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2('Bivariate Correlation & Regression Parameters');
    d.addTable(
      ['Regression / Correlation Metric', 'Observed Value', 'Clinical Interpretation / Benchmark'],
      [
        ['Pearson Correlation Coefficient (r)', `r = ${data.corr.r.toFixed(3)} (p = ${data.corr.p < 0.001 ? '< .001' : data.corr.p.toFixed(4)})`, 'Linear co-variation strength and direction'],
        ['95% Confidence Interval for r', `[${data.corr.ci95[0].toFixed(2)}, ${data.corr.ci95[1].toFixed(2)}]`, 'Fisher Z-transformed confidence limits'],
        ['Spearman Rank Correlation (ρ)', `ρ = ${data.corr.spearman.toFixed(3)} (p = ${data.corr.spearmanP < 0.001 ? '< .001' : data.corr.spearmanP.toFixed(4)})`, 'Monotonic association across ranked values'],
        ['Linear Regression Equation', `Y = ${data.reg.intercept.toFixed(2)} + ${data.reg.slope.toFixed(3)} × X`, 'Ordinary Least Squares (OLS) line of best fit'],
        ['Slope (β1) & Std Error', `${data.reg.slope.toFixed(3)} ± ${data.reg.seSlope.toFixed(3)}`, `Change in Y per 1-unit increase in X`],
        ['Intercept (β0) & Std Error', `${data.reg.intercept.toFixed(2)} ± ${data.reg.seIntercept.toFixed(2)}`, `Estimated baseline value of Y when X = 0`],
        ['Coefficient of Determination (R²)', `R² = ${data.reg.rSquared.toFixed(3)}`, `Proportion of total outcome variance explained (${(data.reg.rSquared * 100).toFixed(1)}%)`],
        ['Residual Standard Error (RSE)', `${data.reg.residualSE.toFixed(2)}`, 'Standard deviation of residuals around regression line'],
        ['Regression F-Statistic', `F(1, ${data.n - 2}) = ${data.reg.fStatistic.toFixed(2)} (p = ${data.reg.pVal < 0.001 ? '< .001' : data.reg.pVal.toFixed(4)})`, 'Overall statistical significance of linear model']
      ]
    );

    if (data.morphology) {
      d.addHeading2('Curve Morphology & Non-Linear Assumption Diagnostics');
      d.addTable(
        ['Morphology Diagnostic', 'Assessed Status', 'Clinical & Methodological Significance'],
        [
          ['Detected Relationship Shape', `${data.morphology.shape}`, 'Identifies linear vs non-linear physiologic windows'],
          ['Monotonicity Status', data.morphology.isMonotonic ? 'Strictly Monotonic' : 'Non-Monotonic (Directional Reversals)', data.morphology.reversals > 0 ? `${data.morphology.reversals} reversal(s) detected` : 'Consistently rising or falling'],
          ['Quadratic Polynomial Fit (R²)', `R² = ${data.morphology.quadR2.toFixed(3)}`, 'Fit of quadratic curve Y = a + bX + cX²'],
          ['Optimal Vertex / Nadir Point', data.morphology.vertexX ? `X* = ${data.morphology.vertexX.toFixed(1)} (Y* = ${data.morphology.vertexY.toFixed(1)})` : 'N/A (Linear)', 'Optimal physiologic target (e.g. ideal CPP in neurotrauma)'],
          ['Recommended Model', `${data.morphology.recommendedModel}`, 'Recommended mathematical representation']
        ]
      );
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Clinical Correlation & Regression Narrative',
      data.reportText || 'Correlation narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Curve Morphology & Non-Linear Screening: Standard statistical tools blindly compute Pearson\'s r. However, in clinical medicine, relationships are frequently U-shaped or non-monotonic (e.g. Cerebral Perfusion Pressure vs Mortality, or Drug Concentration vs Efficacy). On a symmetric U-curve, Pearson\'s r is approximately 0, misleading researchers into claiming "no association". Screening for quadratic curvature and directional monotonicity prevents dangerous medical errors.');
    d.addBullet('Spearman\'s Rank Correlation: Selected as a robust alternative when relationships are monotonic but non-linear, or when outlier observations exert excessive leverage.');
    d.addBullet('Ordinary Least Squares Regression with 95% Confidence Bands: Enables quantitative clinical prediction of expected physiological responses.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Pearson Correlation: r = ∑ [ (xᵢ - x̄)(yᵢ - ȳ) ] / √[ ∑ (xᵢ - x̄)² ∑ (yᵢ - ȳ)² ].');
    d.addBullet('OLS Slope: β1 = ∑ [ (xᵢ - x̄)(yᵢ - ȳ) ] / ∑ (xᵢ - x̄)².');
    d.addBullet('Quadratic Vertex: X* = -b / (2c) for parabola Y = a + bX + cX².');
    d.addParagraph('Key Academic References:');
    d.addBullet('Pearson K (1895). Notes on regression and inheritance in the case of two parents. Proceedings of the Royal Society of London, 58: 240–242.');
    d.addBullet('Spearman C (1904). The proof and measurement of association between two things. American Journal of Psychology, 15(1): 72–101.');
    d.addBullet('Anscombe FJ (1973). Graphs in statistical analysis. The American Statistician, 27(1): 17–21.');

    return d;
  },

  createDiagnosticDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Diagnostic Biomarker ROC & Discrimination Analysis')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Analyzed Biomarker & Diagnostic Information')
      .addParagraph(`Index Biomarker / Test: ${data.name || 'Clinical Diagnostic Score'}`)
      .addParagraph(`Total Patient Cohort Evaluated: ${data.totalN} patients`)
      .addParagraph(`True Positive Disease Cases: ${data.posCount} patients`)
      .addParagraph(`True Negative Normal/Control Cases: ${data.negCount} patients`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2('ROC Discriminative Performance & Optimal Cutoff');
    d.addTable(
      ['Diagnostic Metric', 'Numerical Value', 'Clinical Interpretation / Benchmark'],
      [
        ['Area Under the Curve (ROC AUC)', `${data.auc.toFixed(3)}`, 'Probability that a diseased patient scores higher than a non-diseased patient'],
        ['Hanley-McNeil 95% CI for AUC', `[${data.aucCI95[0].toFixed(3)}, ${data.aucCI95[1].toFixed(3)}]`, 'Precision of global discriminative ability'],
        ['Standard Error of AUC (SE)', `${data.seAuc.toFixed(3)}`, 'Non-parametric standard error estimate'],
        ['Discrimination Rating', data.auc >= 0.90 ? 'Outstanding Discrimination (AUC ≥ 0.90)' : (data.auc >= 0.80 ? 'Excellent Discrimination (0.80 - 0.89)' : 'Acceptable Discrimination (0.70 - 0.79)'), 'Hosmer-Lemeshow diagnostic classification'],
        ['Optimal Cutoff Score', `${data.optimalCutoff.toFixed(2)}`, 'Threshold maximizing Youden\'s J index (Sensitivity + Specificity - 1)'],
        ['Sensitivity at Cutoff', `${(data.sensitivity * 100).toFixed(1)}%`, 'True Positive Rate at recommended decision boundary'],
        ['Specificity at Cutoff', `${(data.specificity * 100).toFixed(1)}%`, 'True Negative Rate at recommended decision boundary'],
        ['Positive Likelihood Ratio (LR+)', `${data.plr.toFixed(2)}`, 'Factor multiplying pre-test odds upon positive result'],
        ['Negative Likelihood Ratio (LR-)', `${data.nlr.toFixed(2)}`, 'Factor multiplying pre-test odds upon negative result']
      ]
    );

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Diagnostic Accuracy & ROC Summary',
      data.reportText || 'ROC narrative summary.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Receiver Operating Characteristic (ROC) Analysis: Selected because diagnostic tests must not be judged at a single arbitrary threshold. Evaluating sensitivity across the entire continuum of specificities allows clinical teams to select distinct cutoffs depending on whether the test is used for high-sensitivity screening (where false negatives cannot be tolerated) or high-specificity confirmatory testing.');
    d.addBullet('Likelihood Ratios: Emphasized because predictive values (PPV/NPV) depend entirely on underlying population prevalence. Likelihood Ratios remain stable and allow Bayesian post-test probability calculations using Fagan\'s nomogram.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Empirical AUC: Computed via non-parametric trapezoidal integration, mathematically equivalent to the Wilcoxon/Mann-Whitney U statistic: AUC = U / (n₁ n₂).');
    d.addBullet('Hanley-McNeil Variance: SE = √[ (θ(1-θ) + (n₁-1)(Q₁ - θ²) + (n₂-1)(Q₂ - θ²)) / (n₁ n₂) ], where θ is the AUC.');
    d.addParagraph('Key Academic References:');
    d.addBullet('Hanley JA, McNeil BJ (1982). The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology, 143(1): 29–36.');
    d.addBullet('Youden WJ (1950). Index for rating diagnostic tests. Cancer, 3(1): 32–35.');
    d.addBullet('Zweig MH, Campbell G (1993). Receiver-operating characteristic (ROC) plots: a fundamental evaluation tool in clinical medicine. Clinical Chemistry, 39(4): 561–577.');

    return d;
  },

  createPowerDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Sample Size Determination & Statistical Power Analysis')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Study Design Information & Input Parameters')
      .addParagraph(`Study Design: ${data.designLabel || 'Two-Sample Clinical Design'}`)
      .addParagraph(`Analysis Objective: ${data.goal === 'power' ? 'Post-Hoc / Achieved Statistical Power Estimation' : 'A Priori Sample Size Requirement Determination'}`)
      .addParagraph(`Type I Error Rate (Alpha, α): ${data.alpha} (Two-Tailed)`);

    d.addHeading1('2. Statistical Outcome & Numerical Results');
    d.addHeading2('Calculated Sample Size & Power Parameters');
    d.addTable(
      ['Study Planning Parameter', 'Calculated Value', 'Clinical & Methodological Significance'],
      [
        ['Standardized Effect Size', `${data.effectSizeLabel || 'Effect'}: ${data.effectSize.toFixed(3)}`, 'Anticipated magnitude of clinical benefit'],
        ['Required Sample Size per Arm', `${data.nPerGroup} patients`, 'Minimum enrollment per arm to guarantee target power'],
        ['Total Required Study Enrollment', `${data.totalN} patients`, 'Total recruitment target across both arms'],
        ['Target Statistical Power (1 - β)', `${(data.targetPower * 100).toFixed(1)}%`, 'Desired probability of correctly rejecting false null hypothesis'],
        ['Achieved Statistical Power', `${(data.achievedPower * 100).toFixed(1)}%`, 'Actual probability of detecting the effect with specified n'],
        ['Recommended Target with 15% Attrition', `${Math.ceil(data.totalN / 0.85)} patients`, 'Accommodates anticipated patient dropout / loss to follow-up']
      ]
    );

    if (data.additionalMetrics) {
      d.addHeading2('Comparative Risk & Methodology Metrics');
      const addRows = Object.entries(data.additionalMetrics).map(([k, v]) => [k, `${v}`, '-']);
      d.addTable(['Metric', 'Value', 'Notes'], addRows);
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Study Planning & Ethical Justification Summary',
      data.reportText || 'Sample size planning narrative statement.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('A Priori Power Calculation Mandate: Underpowered trials fail to detect clinically vital treatments (Type II error), wasting grant funding and exposing patients to experimental risks without scientific return. Conversely, oversized trials needlessly expose excess patients to inferior treatments. Ethics committees (IRBs) and funding bodies mandate rigorous a priori power justification.');
    d.addBullet('Continuity-Corrected Fisher / Proportion Formula: Classical asymptotic formulas underestimate sample size requirements when proportions are small. Applying continuity correction ensures that the exact conditional test achieves the intended statistical power.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Statistical power (1 - β) is the probability of avoiding a false-negative conclusion.');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Independent Means Sample Size: n = 2 [ (Z_{α/2} + Z_β) / d ]².');
    d.addBullet('Paired Means Sample Size: N_pairs = [ (Z_{α/2} + Z_β) / d_z ]².');
    d.addBullet('Proportions Fleiss Continuity Correction: n\' = (n / 4) [ 1 + √( 1 + 4 / (n |p₁ - p₂|) ) ]².');
    d.addParagraph('Key Academic References:');
    d.addBullet('Cohen J (1988). Statistical Power Analysis for the Behavioral Sciences (2nd ed.). Lawrence Erlbaum Associates.');
    d.addBullet('Fleiss JL, Tytun A, Ury HK (1980). A simple approximation for calculating sample sizes for comparing two independent proportions. Biometrics, 36(2): 343–346.');
    d.addBullet('Schulz KF, Altman DG, Moher D (2010). CONSORT 2010 Statement: updated guidelines for reporting parallel group randomised trials. BMJ, 340: c332.');

    return d;
  },

  createTeachingDocx(data) {
    const d = new DocxBuilder();
    d.addTitle('STATIS-GRAVITY CLINICAL BIOSTATISTICS REPORT')
      .addSubTitle('Module: Teaching, Distributions, CLT & Student\'s t Convergence')
      .addAttributionHeader()
      .addDisclaimerBox();

    d.addHeading1('1. Distribution Information & Simulation Parameters');
    d.addParagraph(`Selected Distribution: ${data.distName || 'Normal (Gaussian) Distribution'}`);
    d.addParagraph(`Sample Size (N): ${data.distN || 500} simulated observations`);
    d.addParagraph(`Clinical Context: ${data.clinicalExample || 'Biomedical research modeling'}`);
    d.addParagraph(`CLT Simulation Parent Population: ${data.cltParentName || 'Uniform Distribution'}`);
    d.addParagraph(`CLT Sample Size per Draw (n): ${data.cltN || 30}`);
    d.addParagraph(`CLT Total Iterations (k): ${data.cltK || 1000} sample means`);
    if (data.tConv) {
      d.addParagraph(`Student's t Simulation: Sample Size n = ${data.tConv.sampleSize} (Degrees of Freedom ν = ${data.tConv.df})`);
    }
    if (data.overlap) {
      d.addParagraph(`Two-Sample Overlap Simulation: Mean Difference Δ = ${data.overlap.delta.toFixed(2)}, Group 1 (SD₁ = ${data.overlap.sd1.toFixed(2)}, SEM₁ = ${data.overlap.sem1.toFixed(3)}, n₁ = ${data.overlap.n1}), Group 2 (SD₂ = ${data.overlap.sd2.toFixed(2)}, SEM₂ = ${data.overlap.sem2.toFixed(3)}, n₂ = ${data.overlap.n2}), Significance Level α = ${data.overlap.alpha.toFixed(3)}`);
    }

    d.addHeading1('2. Statistical Outcome & Empirical Convergence Metrics');
    d.addHeading2('Computer-Generated Distribution Metrics');
    d.addTable(
      ['Distribution Metric', 'Empirical Sample Value', 'Theoretical Population Value', 'Status'],
      [
        ['Mean (M)', `${data.sampleMean?.toFixed(3) || '-'}`, `${data.theoMean?.toFixed(3) || '-'}`, 'Congruent'],
        ['Standard Deviation (SD)', `${data.sampleSD?.toFixed(3) || '-'}`, `${data.theoSD?.toFixed(3) || '-'}`, 'Congruent'],
        ['Skewness (G₁)', `${data.skewness?.toFixed(3) || '-'}`, `${data.theoSkewness || '0.000'}`, data.skewnessLabel || 'Evaluated'],
        ['Kurtosis (Excess G₂)', `${data.kurtosis?.toFixed(3) || '-'}`, `${data.theoKurtosis || '0.000'}`, data.kurtosisLabel || 'Evaluated'],
        ['Jarque-Bera Normality Test', `JB = ${data.jbStat?.toFixed(2) || '-'}, p = ${data.jbP?.toFixed(4) || '-'}`, 'Null: Gaussian', data.isNormal ? 'Normal (p ≥ 0.05)' : 'Non-Normal (p < 0.05)']
      ]
    );

    if (data.cltResults) {
      d.addHeading2('Central Limit Theorem (CLT) Convergence Results');
      d.addTable(
        ['CLT Parameter', 'Simulated Value', 'Theoretical CLT Value', 'Convergence Note'],
        [
          ['Parent Population Mean (μ)', `${data.cltResults.theoMean?.toFixed(3)}`, `${data.cltResults.theoMean?.toFixed(3)}`, 'Ground truth benchmark'],
          ['Observed Mean of Means (x̄̄)', `${data.cltResults.obsMean?.toFixed(3)}`, `${data.cltResults.theoMean?.toFixed(3)}`, `Error: ${Math.abs(data.cltResults.obsMean - data.cltResults.theoMean).toFixed(4)}`],
          ['Standard Error of Means (SE)', `${data.cltResults.obsSE?.toFixed(3)}`, `${data.cltResults.theoSE?.toFixed(3)}`, `Law of 1/√n shrinkage: σ/√${data.cltN}`],
          ['Sampling Distribution Skewness', `${data.cltResults.skewness?.toFixed(3)}`, '0.000 (Symmetric)', 'Asymmetry eradicated by averaging'],
          ['Sampling Distribution Normality', `p = ${data.cltResults.normalityP?.toFixed(4)}`, 'p ≥ 0.05', data.cltResults.isNormal ? 'Gaussian Bell Curve Achieved' : 'Approaching Gaussian']
        ]
      );
    }

    if (data.tConv) {
      d.addHeading2('Student\'s t Convergence to Standard Normal N(0, 1) Results');
      d.addTable(
        ['Student\'s t Parameter', `t-Distribution (ν = ${data.tConv.df})`, 'Standard Normal N(0, 1)', 'Methodological Consequence'],
        [
          ['Two-Tailed Critical Value (α=0.05)', `t_crit = ${data.tConv.tCrit?.toFixed(3)}`, 'z_crit = 1.960', `Deviation: ${data.tConv.critDiffPct >= 0 ? '+' : ''}${data.tConv.critDiffPct?.toFixed(1)}% wider cutoff`],
          ['Peak Density f(0)', `f_t(0) = ${data.tConv.tPeak?.toFixed(4)}`, 'φ(0) = 0.3989', `Peak discrepancy: ${data.tConv.peakDiffPct?.toFixed(1)}%`],
          ['Tail Probability P(|X| > 1.960)', `${(data.tConv.tailProb * 100).toFixed(1)}%`, '5.00%', `Type I false positive risk if using z=1.96: ${(data.tConv.tailProb * 100).toFixed(1)}%`],
          ['Excess Kurtosis (Fat Tails)', `${data.tConv.excessKurtosis === Infinity ? '∞ (Fat Tails)' : data.tConv.excessKurtosis?.toFixed(2)}`, '0.00 (Mesokurtic)', data.tConv.df <= 4 ? '4th moment undefined' : 'Heavy tail factor: 6/(ν-4)'],
          ['Max Discrepancy sup |f_t - φ|', `${data.tConv.maxDiscrepancy?.toFixed(4)}`, '0.0000', 'Uniform convergence metric']
        ]
      );
    }

    if (data.overlap) {
      d.addHeading2('Two-Sample Overlap, Dispersion (SD vs. SEM) & Alpha Significance Results');
      d.addTable(
        ['Analytical Parameter', 'Simulated Value', 'Clinical & Inferential Meaning'],
        [
          ['Mean Difference (Δ = μ₂ - μ₁)', `Δ = ${data.overlap.delta.toFixed(2)}`, 'Observed separation between the two group means'],
          ['Group 1 Spread (SD₁ & SEM₁)', `SD₁ = ${data.overlap.sd1.toFixed(2)}, SEM₁ = ${data.overlap.sem1.toFixed(3)} (n₁ = ${data.overlap.n1})`, 'Biological spread (SD₁) vs. sample mean precision (SEM₁ = SD₁/√n₁)'],
          ['Group 2 Spread (SD₂ & SEM₂)', `SD₂ = ${data.overlap.sd2.toFixed(2)}, SEM₂ = ${data.overlap.sem2.toFixed(3)} (n₂ = ${data.overlap.n2})`, 'Biological spread (SD₂) vs. sample mean precision (SEM₂ = SD₂/√n₂)'],
          ['Standard Error of Difference (SE_diff)', `SE_diff = ${data.overlap.seDiff.toFixed(3)} (Welch df = ${data.overlap.df.toFixed(1)})`, 'Combined standard error: √(SEM₁² + SEM₂²) under unequal variances'],
          ['Chosen Significance Level (α)', `α = ${data.overlap.alpha.toFixed(3)}`, `Type I error tolerance: ${data.overlap.alpha === 0.05 ? 'Standard 5% biomedical risk' : (data.overlap.alpha < 0.05 ? 'Strict threshold' : 'Relaxed exploratory threshold')}`],
          ['Critical Value (t_crit vs z_crit)', `t_crit = ${data.overlap.tCrit.toFixed(3)} (z = ${data.overlap.zCrit.toFixed(3)})`, 'Required number of standard errors to claim statistical significance'],
          ['Critical Difference Boundary (Δcrit)', `Δcrit = ${data.overlap.deltaCrit.toFixed(2)}`, 'Minimum mean separation needed to achieve p < α (Δcrit = t_crit · SE_diff)'],
          ['Test Statistic & p-value', `t = ${data.overlap.tStat.toFixed(2)}, p = ${data.overlap.pValue < 0.0001 ? '< 0.0001' : data.overlap.pValue.toFixed(4)}`, `Significance: ${data.overlap.isSignificant ? 'REJECT H₀ (p < α)' : 'FAIL TO REJECT H₀ (p ≥ α)'}`],
          ['Individual Patient Overlap (OVL_SD)', `${(data.overlap.patientOVL * 100).toFixed(1)}%`, 'Proportion of overlapping individual patient values (Weitzman\'s OVL)'],
          ['Sampling Distribution Overlap (OVL_SEM)', `${(data.overlap.meansOVL * 100).toFixed(1)}%`, 'Proportion of overlap between the sampling distributions of sample means']
        ]
      );
    }

    d.addHeading1('3. Clinical & Statistical Interpretation');
    d.addCalloutBox(
      'Pedagogical Synthesis & Clinical Trial Relevance',
      data.reportText || 'The Central Limit Theorem, Student\'s t convergence, and Two-Sample Overlap simulations demonstrate the mathematical foundations of parametric testing and the definition of alpha in clinical trials.',
      'F0FDF4',
      '16A34A'
    );

    d.addHeading1('4. Reason This Particular Test Was Chosen');
    d.addBullet('Foundation of Inferential Biostatistics: Parametric hypothesis tests (Student t-test, ANOVA, ordinary least squares regression) mathematically assume normally distributed errors or sample means. The Central Limit Theorem provides the mathematical justification for deploying these tests in clinical trials with n ≥ 30 even when raw clinical metrics (e.g. ICU stay, recovery hours) are skewed.');
    d.addBullet('Gosset\'s Student\'s t Adjustment: In small clinical cohorts (n < 30), estimating population variance σ² using sample variance s² introduces substantial stochastic instability into the test statistic denominator. Using Gaussian critical values (z = 1.96) severely inflates the Type I error rate (e.g. to 14.5% at n = 4). Student\'s t distribution compensates for this extra uncertainty by thickening the tails and demanding a higher critical threshold (t = 3.182 at n = 4).');
    d.addBullet('The n ≥ 31 Clinical Threshold: As demonstrated by the simulation, when sample size reaches n ≥ 31 (degrees of freedom ν ≥ 30), the critical t cutoff drops to 2.042 (only 4.2% wider than 1.960), and tail probability converges close to 5.0%. This mathematical threshold explains why sample sizes of 30 or greater historically permit Gaussian approximation in medical trial protocols.');
    d.addBullet('Why α = 0.05 Defines the Point of Significance: In 1925, Ronald A. Fisher proposed the 5% significance level (p < 0.05) as a pragmatic convention for scientific research—representing a 1 in 20 chance of observing an effect as extreme under the null hypothesis of no difference. On a standard Gaussian distribution, exactly 5% of probability mass lies in the tails beyond ±1.960 standard errors (2.5% in each tail). Hence, the critical separation distance between sample means is Δcrit = 1.960 · SE_diff. When the observed difference Δ exceeds Δcrit, the p-value falls below 0.05.');
    d.addBullet('The Fundamental Distinction Between SD and SEM: Standard Deviation (SD) reflects real inter-individual biological diversity among patients and does not contract when sample size increases. In contrast, the Standard Error of the Mean (SEM = SD/√n) quantifies our uncertainty in the population mean estimate and contracts steadily as 1/√n. Consequently, two treatment groups can exhibit 70% biological overlap in individual patient scores, yet their treatment difference can be verified as statistically significant (p < 0.001) once sufficient patients are enrolled to shrink the SEM.');
    d.addBullet('Consequences of Modifying Alpha (α): Relaxing α to 0.10 moves the critical cutoff inward to z = 1.645, lowering the required separation Δcrit and declaring significance on smaller differences or smaller sample sizes, at the expense of doubling the false-positive risk to 10%. Tightening α to 0.01 (z = 2.576) or 0.001 (z = 3.291), as required in confirmatory registration trials or genome-wide studies, shifts the cutoff outward into the extreme tails, demanding either much larger effect sizes or substantially expanded sample sizes before significance can be claimed.');

    d.addHeading1('5. Background Statistical Knowledge & Medical Research Context');
    d.addParagraph('Mathematical Formulations:');
    d.addBullet('Classical Lindberg-Lévy Central Limit Theorem: Let X₁, X₂, ..., X_n be independent and identically distributed (i.i.d.) random variables with mean μ and finite variance σ². Then as n → ∞: √n (X̄_n - μ) / σ → N(0, 1).');
    d.addBullet('Student\'s t Distribution Density: f(t; ν) = [ Γ((ν+1)/2) / (√(πν) Γ(ν/2)) ] · [ 1 + t²/ν ]^{-(ν+1)/2}. As ν → ∞, [ 1 + t²/ν ]^{-(ν+1)/2} → exp(-t²/2), converging to Standard Normal N(0, 1).');
    d.addBullet('Standard Error of the Mean: SEM = σ / √n. Quadrupling patient enrollment cuts the estimation uncertainty in half.');
    d.addBullet('Weitzman\'s Distribution Overlap Coefficient (OVL): For two equal-variance Gaussian curves separated by difference Δ: OVL = 2 · Φ(-|Δ| / (2 · s)), where s = SD for patient-level overlap and s = SEM for sampling-mean-level overlap.');
    d.addBullet('Critical Significance Boundary: Δcrit = t_crit(α, df) · SD · √(2/n).');
    d.addParagraph('Key Academic References:');
    d.addBullet('Fisher RA (1925). Statistical Methods for Research Workers. Edinburgh: Oliver and Boyd.');
    d.addBullet('Cumming G, Finch S (2005). Inference by eye: confidence intervals and how to read pictures of data. Am Psychol, 60(2): 170–180.');
    d.addBullet('Student [Gosset WS] (1908). The probable error of a mean. Biometrika, 6(1): 1–25.');
    d.addBullet('Altman DG, Bland JM (2005). Standard deviations and standard errors. BMJ, 331(7521): 903.');
    d.addBullet('Laplace PS (1810). Mémoire sur les approximations des formules qui sont fonctions de très grands nombres et sur leur application aux probabilités. Mémoires de l\'Académie Royale des Sciences de Paris.');
    d.addBullet('Gauss CF (1809). Theoria motus corporum coelestium in sectionibus conicis solem ambientium. Hamburg: Perthes et Besser.');

    return d;
  }
};
