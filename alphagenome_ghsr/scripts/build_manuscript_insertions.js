const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell,
  AlignmentType, HeadingLevel, WidthType, ShadingType, BorderStyle,
  PageBreak, LevelFormat,
} = require("docx");

const FA = "Tahoma";
const EN = "Times New Roman";

// ---------- helpers ----------
const fa = (text, o = {}) =>
  new Paragraph({
    bidirectional: true,
    alignment: AlignmentType.RIGHT,
    spacing: { before: o.before ?? 60, after: o.after ?? 60, line: 300 },
    children: [new TextRun({ text, font: FA, size: o.size ?? 21, bold: o.bold, rightToLeft: true, color: o.color })],
  });

// Persian heading for each section of the guide
const faHead = (text) =>
  new Paragraph({
    bidirectional: true,
    alignment: AlignmentType.RIGHT,
    spacing: { before: 360, after: 120 },
    keepNext: true,
    border: { bottom: { style: BorderStyle.SINGLE, size: 8, color: "2A78D6", space: 2 } },
    children: [new TextRun({ text, font: FA, size: 28, bold: true, rightToLeft: true, color: "1C5CAB" })],
  });

// Location box: Persian instruction + optional English anchor sentence
const LOC_FILL = "EEF3FA";
function loc(numberLabel, instructionFa, anchorEn, anchorLabelFa = "بعد از این جمله:") {
  const paras = [
    new Paragraph({
      bidirectional: true,
      alignment: AlignmentType.RIGHT,
      shading: { type: ShadingType.CLEAR, fill: LOC_FILL, color: "auto" },
      spacing: { before: 200, after: 0, line: 300 },
      keepNext: true,
      border: { right: { style: BorderStyle.SINGLE, size: 24, color: "2A78D6", space: 6 } },
      children: [
        new TextRun({ text: `محل درج ${numberLabel}: `, font: FA, size: 21, bold: true, rightToLeft: true, color: "1C5CAB" }),
        new TextRun({ text: instructionFa, font: FA, size: 21, rightToLeft: true }),
      ],
    }),
  ];
  if (anchorEn) {
    paras.push(
      new Paragraph({
        bidirectional: true,
        alignment: AlignmentType.RIGHT,
        shading: { type: ShadingType.CLEAR, fill: LOC_FILL, color: "auto" },
        spacing: { before: 0, after: 0 },
        keepNext: true,
        border: { right: { style: BorderStyle.SINGLE, size: 24, color: "2A78D6", space: 6 } },
        children: [new TextRun({ text: anchorLabelFa, font: FA, size: 19, rightToLeft: true, color: "52514E" })],
      }),
      new Paragraph({
        alignment: AlignmentType.LEFT,
        shading: { type: ShadingType.CLEAR, fill: LOC_FILL, color: "auto" },
        spacing: { before: 0, after: 120 },
        keepNext: true,
        border: { right: { style: BorderStyle.SINGLE, size: 24, color: "2A78D6", space: 6 } },
        children: [new TextRun({ text: `“${anchorEn}”`, font: EN, size: 20, italics: true, color: "52514E" })],
      })
    );
  }
  return paras;
}

// Label introducing the English text to insert
const insertLabel = (text = "متن جدید برای درج (انگلیسی):") =>
  new Paragraph({
    bidirectional: true,
    alignment: AlignmentType.RIGHT,
    spacing: { before: 120, after: 60 },
    keepNext: true,
    children: [new TextRun({ text, font: FA, size: 19, bold: true, rightToLeft: true, color: "008300" })],
  });

// English manuscript text (what actually goes into the paper)
const en = (text, o = {}) =>
  new Paragraph({
    alignment: o.align ?? AlignmentType.JUSTIFIED,
    spacing: { before: o.before ?? 0, after: o.after ?? 140, line: 360 },
    indent: o.indent,
    children: (Array.isArray(text) ? text : [text]).map((t) =>
      typeof t === "string" ? new TextRun({ text: t, font: EN, size: 24 }) : new TextRun({ font: EN, size: 24, ...t })
    ),
  });

// English sub-heading to be inserted into the manuscript
const enHead = (text) =>
  new Paragraph({
    spacing: { before: 120, after: 120 },
    keepNext: true,
    children: [new TextRun({ text, font: EN, size: 24, bold: true })],
  });

const rule = () =>
  new Paragraph({
    spacing: { before: 120, after: 120 },
    border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: "C3C2B7", space: 1 } },
    children: [],
  });

// ---------- Table 17 ----------
const cellBorder = { style: BorderStyle.SINGLE, size: 4, color: "898781" };
const borders = { top: cellBorder, bottom: cellBorder, left: cellBorder, right: cellBorder };
const colW = [1500, 1500, 1350, 1250, 1950, 2016]; // sum 9566 (A4 text width ~9026-9638)
const TW = colW.reduce((a, b) => a + b, 0);
function cell(text, w, head = false) {
  const parts = Array.isArray(text) ? text : [text];
  return new TableCell({
    width: { size: w, type: WidthType.DXA },
    borders,
    margins: { top: 60, bottom: 60, left: 80, right: 80 },
    shading: head ? { type: ShadingType.CLEAR, fill: "E1E0D9", color: "auto" } : undefined,
    children: parts.map(
      (p) =>
        new Paragraph({
          spacing: { before: 0, after: 0, line: 240 },
          children: [new TextRun({ text: p, font: EN, size: 17, bold: head })],
        })
    ),
  });
}
const T17 = [
  ["Variant", "Predicted maximum change in GHSR expression (|log2FC|)*", "Median percentile within the amplicon ISM distribution†", "Direction across 16 conditions (ALT higher / lower)", "Genotypes detected (HW / Het / HM)", "Statistical association (present study)"],
  [["rs2922126", "chr3:172,449,471 T>A", "1,015 bp upstream"], "0.027 (≈1.9%)", "43.7", "15 / 1", ["Cases: 31 / 16 / 3", "Controls: 28 / 18 / 4"], "No association with cognitive screening status (adjusted P = 0.862); nominal biochemical associations not significant after FDR correction"],
  [["rs115459563", "chr3:172,449,444 G>A", "988 bp upstream"], "0.016 (≈1.1%)", "18.9", "5 / 11", ["Cases: 50 / 0 / 0", "Controls: 43 / 7 / 0"], "Genotype distribution differed nominally (P = 0.012); regression models did not converge"],
  [["rs34136054", "chr3:172,449,372 T>C‡", "916 bp upstream"], "0.025 (≈1.7%)", "42.1", "15 / 1", ["Cases: 50 / 0 / 0", "Controls: 49 / 1 / 0"], "Not evaluable (one carrier)"],
];
const table17 = new Table({
  width: { size: TW, type: WidthType.DXA },
  columnWidths: colW,
  rows: T17.map((r, i) => new TableRow({ tableHeader: i === 0, children: r.map((c, j) => cell(c, colW[j], i === 0)) })),
});

// ---------- Figure ----------
const figPath = require("path").join(__dirname, "..", "results", "figures", "figure2_predicted_effects.png");
const figure = new Paragraph({
  alignment: AlignmentType.CENTER,
  spacing: { before: 120, after: 120 },
  children: [new ImageRun({ type: "png", data: fs.readFileSync(figPath), transformation: { width: 600, height: 285 } })],
});

// ---------- content ----------
const children = [];

// Title + intro (Persian)
children.push(
  new Paragraph({
    bidirectional: true, alignment: AlignmentType.CENTER, spacing: { after: 120 },
    children: [new TextRun({ text: "متن‌های پیشنهادی برای افزودن تحلیل AlphaGenome 2 به مقاله", font: FA, size: 32, bold: true, rightToLeft: true })],
  }),
  new Paragraph({
    alignment: AlignmentType.CENTER, spacing: { after: 240 },
    children: [new TextRun({ text: "Ghrelin Receptor Gene Polymorphisms, Cognitive Screening Status, and Osteoporosis-Related Parameters", font: EN, size: 20, italics: true, color: "52514E" })],
  }),
  fa("این فایل فقط متن‌هایی را دارد که باید به مقاله اضافه یا در آن اصلاح شوند. برای هر مورد، ابتدا در کادر آبی گفته شده که دقیقاً کجای مقاله باید قرار بگیرد و جملهٔ مرجع از خود مقاله (ARTICLE__revised-v2) نقل شده است تا جای آن را راحت پیدا کنید. سپس متن انگلیسی آمادهٔ درج آمده است."),
  fa("شماره‌گذاری: مقاله در حال حاضر ۵ شکل، ۱۶ جدول و ۸۳ رفرنس دارد؛ بنابراین مورد جدید «Figure 6»، «Table 17» و رفرنس‌های ۸۴ تا ۸۸ نام‌گذاری شده‌اند. اگر شماره‌ها در نسخهٔ نهایی تغییر کرد، آن‌ها را هماهنگ کنید."),
  fa("نکتهٔ مهم در نگارش: این تحلیل بعد از ژنوتیپینگ انجام شده است. در همهٔ متن‌ها آن را «گذشته‌نگر و مکمل» معرفی کرده‌ام و ننوشته‌ام که نتایج مقاله را «اثبات» می‌کند؛ نوشته‌ام که با نتایج «هم‌سو» است و توضیح زیستی ممکنی برای آن‌ها ارائه می‌دهد. این نگارش در برابر داوران قابل دفاع‌تر است.", { color: "52514E" }),
);

// ===================== 1. METHODS =====================
children.push(faHead("۱. بخش روش‌ها (Methods)"));

// 1a correction
children.push(
  ...loc("۱-الف (اصلاح یک جمله)", "بخش Methods، زیربخش «Genetic Analysis». این جمله را حذف و متن زیر را جایگزین آن کنید. دلیل: در ژنوم مرجع، آلل اصلی rs34136054 حرف T است نه C، و جهت رشتهٔ DNA برای سه SNP یکسان نوشته نشده بود.",
    "The three target variants located on chromosome 3 (GRCh38/hg38), rs2922126 (A>T), rs115459563 (G>A), and rs34136054 (C>T), were captured within the amplified 399-bp promoter/regulatory fragment of GHSR, selected based on its role in transcriptional regulation and prior metabolic associations.",
    "جملهٔ فعلی که باید جایگزین شود:"),
  insertLabel("متن جایگزین (انگلیسی):"),
  en("The three target variants located on chromosome 3 (GRCh38/hg38), rs2922126 (chr3:172,449,471 T>A), rs115459563 (chr3:172,449,444 G>A) and rs34136054 (chr3:172,449,372 T>C/G), were captured within the amplified 399-bp upstream regulatory fragment of GHSR (chr3:172,449,261–172,449,659), located 916–1,015 bp upstream of the GHSR transcription start site. This fragment was selected based on its role in transcriptional regulation and prior metabolic associations. Alleles are reported on the GRCh38 forward strand; because GHSR is transcribed from the reverse strand, chromatograms read in the direction of transcription show the complementary bases."),
);

// 1b new subsection
children.push(
  ...loc("۱-ب (زیربخش جدید)", "بخش Methods، بعد از پایان زیربخش «Genetic Analysis» و درست قبل از عنوان «Statistical Analysis». یک زیربخش کامل جدید با عنوان زیر اضافه کنید.",
    "Genotype distributions in the control group conformed to Hardy–Weinberg equilibrium for rs2922126 (p = 0.62), rs115459563 (p = 0.59), and rs34136054 (p = 0.94)."),
  insertLabel(),
  enHead("In silico prediction of regulatory variant effects"),
  en("To explore whether the genotyped variants could plausibly alter GHSR transcriptional regulation, a retrospective, complementary in silico analysis was performed after the genotyping and statistical analyses had been completed; it did not inform variant selection, study design or the statistical analysis plan. Only public variant coordinates and public ontology identifiers were used; no participant-level data were used in, or transmitted during, this analysis."),
  en("Variant coordinates and alleles were retrieved independently from NCBI dbSNP (build 157) (86) and Ensembl (GRCh38) (87), and reference alleles were verified against the GRCh38 reference sequence. Variant consequences relative to GHSR transcripts were obtained with the Ensembl Variant Effect Predictor (87), overlap with annotated regulatory elements was assessed using the Ensembl Regulatory Build, and population allele frequencies were obtained from gnomAD (88). The published primer pair (Table 1) was mapped onto GRCh38 to confirm the location of the sequenced fragment. Because dbSNP lists rs34136054 as tri-allelic (T>C and T>G), both alternate alleles were analysed."),
  en("Regulatory effects were predicted with AlphaGenome 2 (Google DeepMind) (84, 85), a deep-learning model that predicts gene expression, chromatin accessibility, transcription-factor binding, histone modifications and other regulatory read-outs directly from DNA sequence in a specified tissue or cell type. The model was accessed through the official Python package (alphagenome version 2.0.0) on 14 November 2026. For each variant, predictions for the reference and alternate alleles were compared using all 22 officially recommended variant-scoring methods within a 1-Mb sequence window centred on the variant. Effects are expressed as log2 fold changes (alternate versus reference) of predicted signal."),
  en("Sixteen tissues and cell types were selected to reflect the phenotypes measured in this study: seven brain or neuronal conditions (hypothalamus, hippocampus, hippocampal astrocyte, whole brain, temporal lobe, brain neuron, glutamatergic neuron), pituitary gland, four bone-related conditions (osteoblast, osteocyte, mesenchymal stem cell, and blood monocyte as a surrogate for the osteoclast lineage, for which no model condition was available) and four reference tissues (stomach, testis, thyroid gland, ovary). For each prediction, the model indicates whether the underlying data were experimentally observed or imputed, and all summaries were repeated using observed data only."),
  en("Because the model's scores have no intrinsic significance threshold, two internal reference analyses were performed. First, in silico mutagenesis was used to score all 1,197 possible single-nucleotide substitutions across the 399-bp sequenced fragment, providing a local reference distribution against which the three variants could be compared. Second, all predictions were repeated at four sequence-window lengths (16 kb to 1 Mb) to assess whether the direction of predicted effects was stable."),
);

// 1c stats sentence
children.push(
  ...loc("۱-ج (یک جملهٔ کوتاه)", "بخش Methods، انتهای زیربخش «Statistical Analysis». این جمله را بعد از آخرین جملهٔ همین زیربخش اضافه کنید.",
    "The statistical significance level was set at P<0.05."),
  insertLabel(),
  en("No inferential statistics, P-values or significance thresholds were applied to the in silico predictions, which are model outputs rather than measurements; they are reported descriptively."),
);

// ===================== 2. RESULTS =====================
children.push(faHead("۲. بخش نتایج (Results)"));
children.push(
  ...loc("۲ (زیربخش جدید)", "بخش Results، انتهای نتایج؛ بعد از آخرین جملهٔ مربوط به رگرسیون لجستیک و درست قبل از عنوان «Discussion» (این جمله دو بار در مقاله آمده؛ منظور مورد دوم است که آخرین جملهٔ Results است). جدول ۱۷ و شکل ۶ (در ادامهٔ همین فایل) هم در این زیربخش ارجاع داده شده‌اند.",
    "Since SNP3 was observed in only one individual, it could not be evaluated, and logistic regression analysis was not performed for this variant."),
  insertLabel(),
  enHead("In silico prediction of regulatory effects"),
  en("Coordinates from dbSNP and Ensembl were in complete agreement, and all reference alleles matched the GRCh38 sequence. Mapping the primer pair onto GRCh38 yielded a 399-bp product (chr3:172,449,261–172,449,659), reproducing the observed product size and containing all three variants. Because GHSR is transcribed from the reverse strand, the variants lie upstream of the gene, 1,015 bp (rs2922126), 988 bp (rs115459563) and 916 bp (rs34136054) from the transcription start site, and were classified as upstream gene variants for all GHSR transcripts. None overlapped an annotated regulatory element; the nearest were an enhancer 90 bp and a CTCF-binding site 349 bp from rs2922126. The nearest other protein-coding gene (TNFSF10) lies 74 kb away."),
  en("Predicted GHSR expression was highest in the pituitary gland, followed by the hypothalamus (41% of the pituitary value) and other brain regions, and was lowest in osteoblasts, osteocytes and mesenchymal stem cells (≤1% of the maximum), consistent with the known tissue distribution of GHSR."),
  en("All three variants showed only very small predicted effects (Table 17, Figure 6). The largest predicted change in GHSR expression in any of the 16 conditions was a log2 fold change of 0.027 for rs2922126 (≈1.9%), 0.016 for rs115459563 and 0.025 for rs34136054; the largest effect on any regulatory read-out was 0.051 (≈3.6%). In contrast, other substitutions within the same 399-bp fragment produced predicted effects 4.7- to 95-fold larger, showing that the model was responsive at this locus. None of the three variants exceeded the 90th percentile of this local distribution in any of 48 tissue–read-out combinations (median percentiles 43.7, 18.9 and 42.1, respectively)."),
  en("No variant showed a consistent direction of effect across tissues (Table 17), and the magnitude of the predicted effect followed the predicted GHSR expression level of each tissue (Spearman ρ = 0.63–0.83). The direction of the predicted expression effect was preserved across all four sequence-window lengths in only 13 of 48 variant–tissue combinations. Predicted effects in bone-related cells were the smallest observed (|log2 fold change| < 0.0004). Restricting the analysis to experimentally observed data identified two small but directionally consistent signals: rs2922126 was predicted to reduce CTCF binding signal in all 12 available observed tracks across eight tissues, and rs34136054 was predicted to slightly increase active-promoter histone marks (H3K4me3, H3K9ac, H3K27ac) in temporal lobe tissue."),
);

// Table 17
children.push(
  ...loc("جدول ۱۷", "در بخش Results، داخل همین زیربخش جدید و بعد از پاراگراف‌های بالا (یا طبق قالب مجله، در انتهای مقاله کنار سایر جدول‌ها)."),
  insertLabel("جدول ۱۷ (انگلیسی):"),
  en([{ text: "Table 17 ", bold: true }, "In silico predicted regulatory effects of the three GHSR variants alongside the experimental genotyping and association findings of the present study."], { align: AlignmentType.LEFT }),
  table17,
  en("HW, homozygous wild-type; Het, heterozygous; HM, homozygous mutant; ISM, in silico mutagenesis; FDR, false discovery rate. *Largest absolute value across 16 tissues and cell types; percentages are the corresponding change in predicted signal. †Percentile of the variant's absolute predicted effect among all 1,197 possible single-nucleotide substitutions in the sequenced 399-bp fragment (median over four tissues and four read-outs). ‡The alternative T>G allele gave a similarly small maximum effect (0.023). Computational predictions are model outputs, not measurements; no P-values are attached to them.", { before: 80 }),
);
// fix the footnote font size
children[children.length - 1] = new Paragraph({
  spacing: { before: 80, after: 200, line: 260 },
  children: [new TextRun({ text: "HW, homozygous wild-type; Het, heterozygous; HM, homozygous mutant; ISM, in silico mutagenesis; FDR, false discovery rate. *Largest absolute value across 16 tissues and cell types; percentages are the corresponding change in predicted signal. †Percentile of the variant's absolute predicted effect among all 1,197 possible single-nucleotide substitutions in the sequenced 399-bp fragment (median over four tissues and four read-outs). ‡The alternative T>G allele gave a similarly small maximum effect (0.023). Computational predictions are model outputs, not measurements; no P-values are attached to them.", font: EN, size: 18 })],
});

// Figure 6
children.push(
  ...loc("شکل ۶", "در بخش Results، داخل همین زیربخش جدید، بعد از جدول ۱۷. فایل با کیفیت بالای شکل (PDF / SVG / PNG ۶۰۰ dpi) در پوشهٔ results/figures با نام figure2_predicted_effects موجود است."),
  insertLabel("شکل ۶ و زیرنویس آن (انگلیسی):"),
  figure,
  new Paragraph({
    alignment: AlignmentType.JUSTIFIED,
    spacing: { before: 60, after: 200, line: 300 },
    children: [
      new TextRun({ text: "Figure 6 In silico predicted regulatory effects of the three GHSR variants. ", font: EN, size: 22, bold: true }),
      new TextRun({ text: "(A) Largest absolute predicted effect (log2 fold change, alternate versus reference allele) of each variant on each regulatory read-out across 16 tissues and cell types (logarithmic axis). The dashed line marks a log2 fold change of 0.1 (≈7% change), which no variant reached. (B) The three variants (circles) compared with all 1,197 possible single-nucleotide substitutions across the sequenced 399-bp fragment (violins) in pituitary gland, scored in the same analysis. Other substitutions in the fragment produced much larger predicted effects, whereas the three variants fell within the central part of the distribution. Predictions were generated with AlphaGenome 2 (84, 85).", font: EN, size: 22 }),
    ],
  }),
);

// ===================== 3. DISCUSSION =====================
children.push(faHead("۳. بخش بحث (Discussion)"));

children.push(
  ...loc("۳-الف (یک جمله در Principal Findings)", "بخش Discussion، زیربخش «Principal Findings»، انتهای پاراگراف سوم. این جمله را بعد از آن اضافه کنید.",
    "Instead, they highlight several preliminary observations that require confirmation in larger studies with adequate statistical power and clearly defined clinical outcomes."),
  insertLabel(),
  en("A complementary in silico analysis was consistent with these findings: none of the three variants was predicted to alter GHSR regulation to an appreciable extent in brain, pituitary or bone-related cells."),
);

children.push(
  ...loc("۳-ب (پاراگراف جدید)", "بخش Discussion، زیربخش «GHSR Polymorphisms and Cognitive-Group Classification». بعد از پاراگرافی که با «However, the present study evaluated genetic variants at the DNA level...» شروع می‌شود و با جملهٔ زیر تمام می‌شود.",
    "The lack of an association between rs2922126 and cognitive-group classification neither confirms nor excludes a role for the ghrelin–GHSR pathway in neurodegeneration."),
  insertLabel(),
  en("To address the absence of functional data at least computationally, we retrospectively evaluated the three variants with AlphaGenome 2, a model that predicts regulatory activity directly from DNA sequence (84, 85). Across 16 tissues and cell types, including hypothalamus, hippocampus, temporal lobe and pituitary gland, none of the variants was predicted to change GHSR expression by more than about 2%. This near-null result is informative rather than a failure of the method: other substitutions within the same 399-bp fragment were predicted to have effects up to 95-fold larger, and none of the three variants ranked among the top 10% of this local distribution. Moreover, no variant showed a consistent direction of effect across tissues, and the small predicted differences largely followed the baseline expression of GHSR in each tissue rather than suggesting tissue-specific regulation."),
  en("These predictions are therefore consistent with the absence of an association between the investigated variants and cognitive-group classification, and they offer a plausible molecular explanation for it: if the variants do not appreciably alter GHSR transcription, a detectable effect on cognitive screening status would not be expected. However, the two analyses address different questions. The in silico analysis concerns predicted cis-regulatory activity, whereas the association analysis concerns phenotype in 100 participants, and neither can confirm the other. Given the low achieved power of the association analysis, the computational result should be regarded as supporting rather than proving the absence of a meaningful effect."),
);

children.push(
  ...loc("۳-ج (پاراگراف جدید دربارهٔ استخوان)", "بخش Discussion، زیربخش «GHSR Polymorphisms, Bone Mineral Density, and Thyroid Function». بعد از پاراگرافی که با جملهٔ زیر تمام می‌شود.",
    "Larger studies are needed to determine whether specific GHSR variants are associated with skeletal outcomes."),
  insertLabel(),
  en("The in silico analysis provides additional context for these skeletal findings. AlphaGenome 2 predicted GHSR expression to be very low in osteoblasts, osteocytes and mesenchymal stem cells (≤1% of the level predicted for the pituitary gland), and the predicted allelic effects in these cells were the smallest of all tissues examined. If GHSR expression in bone cells is indeed low, a direct, cell-autonomous effect of these upstream variants on bone cells would be unlikely, and any influence of the ghrelin–GHSR system on bone would more plausibly act indirectly, for example through pituitary growth hormone secretion or hypothalamic regulation of energy balance. This is consistent with the absence of significant associations between the variants and BMD in the present study. These are predictions rather than measurements, however, and GHSR expression has been reported experimentally in human osteoblasts (56)."),
);

children.push(
  ...loc("۳-د (محدودیت هفتم)", "بخش Discussion، زیربخش «Strengths and Limitations». بعد از پاراگراف «Sixth, ...» که با جملهٔ زیر تمام می‌شود، و قبل از پاراگراف «Taken together, these limitations indicate...».",
    "GHSR expression, circulating ghrelin concentrations, and receptor function were not measured, and the findings have not been replicated in an independent population."),
  insertLabel(),
  en("Seventh, the in silico analysis is based on model predictions rather than experimental measurements, was performed retrospectively, and cannot exclude regulatory mechanisms that the model does not capture, such as distal regulatory elements, combined effects of several variants, or effects in cell types not represented in the model. No osteoclast-lineage condition was available, and some predictions, particularly for transcription-factor binding, were based largely on imputed rather than observed data. The absence of a predicted effect therefore does not demonstrate that the variants are biologically inert, and these predictions are not clinical or diagnostic evidence."),
);

children.push(
  ...loc("۳-هـ (ادامهٔ Future Research Directions)", "بخش «Future Research Directions»، انتهای همان پاراگراف. این جمله‌ها را بعد از آخرین جملهٔ آن اضافه کنید.",
    "...are warranted to clarify whether the ghrelin signaling pathway plays a causal role in the intersection of neurodegeneration and age-related bone loss."),
  insertLabel(),
  en("The in silico predictions also suggest which functional experiments would be most informative. Because the predicted effects are very small, allele-specific expression analysis of GHSR in pituitary, hypothalamic or cortical tissue from heterozygous individuals would be more sensitive than conventional reporter assays. The small but consistent predicted reduction in CTCF binding associated with rs2922126 could be tested directly by allele-specific chromatin immunoprecipitation. Any reporter assay should be performed in GHSR-expressing cells, such as pituitary or hypothalamic neuronal cell lines, rather than in bone cells."),
);

// ===================== 4. CONCLUSION =====================
children.push(faHead("۴. بخش نتیجه‌گیری (Conclusion)"));
children.push(
  ...loc("۴", "بخش Conclusion، پاراگراف اول. بعد از جملهٔ زیر و قبل از پاراگراف «Overall, the present results do not establish...».",
    "Higher educational attainment was independently associated with lower odds of case-group membership, although differences in education and the use of cognitive screening thresholds should be considered when interpreting this finding."),
  insertLabel(),
  en("Consistent with these findings, a complementary in silico analysis using AlphaGenome 2 predicted no appreciable effect of the three variants on GHSR regulation in brain, pituitary or bone-related cells."),
);

// ===================== 5. ABSTRACT (optional) =====================
children.push(faHead("۵. چکیده (Abstract) — اختیاری"));
children.push(
  fa("فقط اگر محدودیت تعداد کلمات مجله اجازه می‌دهد. جمله‌های مربوط به توان آماری و محدودیت‌ها را برای اضافه کردن این جمله حذف نکنید.", { color: "52514E" }),
  ...loc("۵", "چکیده، قبل از جمله‌ای که با «In conclusion, this exploratory study...» شروع می‌شود.",
    "...each additional year of education was associated with lower odds of membership in the high 6-CIT score group (OR = 0.68, p < 0.001) and osteoporosis versus normal BMD remained associated with group classification (OR = 7.25, p = 0.048)."),
  insertLabel(),
  en("A retrospective in silico analysis using AlphaGenome 2 predicted no appreciable effect of any of the three variants on GHSR regulation across 16 brain, pituitary, bone-related and reference tissues."),
);

// ===================== 6. DATA AVAILABILITY / DISCLOSURE =====================
children.push(faHead("۶. بخش‌های اداری مقاله"));
children.push(
  ...loc("۶-الف (Availability of data and materials)", "انتهای بخش «Availability of data and materials». نشانی مخزن کد را جایگزین [URL] کنید.",
    "...are fully available from the corresponding author upon reasonable request."),
  insertLabel(),
  en("The code used for the in silico analysis, the complete model outputs and the supplementary tables are available at [URL]. Only public variant information was used for this analysis."),
  ...loc("۶-ب (Acknowledgments) — پیشنهادی", "بخش Acknowledgments، پاراگرافی که دربارهٔ Grammarly و ابزارهای مدل زبانی است. چون این تحلیل با کمک ابزار هوش مصنوعی انجام شده، پیشنهاد می‌کنم آن پاراگراف را به این صورت گسترش دهید (طبق سیاست مجله تنظیم شود).",
    "Language editing and grammatical proofreading of the manuscript were assisted by Grammarly and large language model tools; the authors reviewed, verified, and take full responsibility for the content of the publication."),
  insertLabel("متن جایگزین پیشنهادی (انگلیسی):"),
  en("Language editing and grammatical proofreading of the manuscript were assisted by Grammarly and large language model tools. An AI coding assistant was used to help implement and document the in silico analysis pipeline. The authors reviewed and verified all analyses and text and take full responsibility for the content of the publication."),
);

// ===================== 7. REFERENCES =====================
children.push(faHead("۷. رفرنس‌های جدید (References)"));
children.push(
  ...loc("۷", "انتهای فهرست References، بعد از رفرنس شمارهٔ ۸۳ (Roe CM, et al. Neurology 2007)."),
  insertLabel(),
);
const refs = [
  "84.\tAvsec Ž, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature. 2026;649:1206–1218. doi:10.1038/s41586-025-10014-0.",
  "85.\tGoogle DeepMind. AlphaGenome 2: deciphering how cells read DNA to better understand disease [Internet]. 2026 Nov 11. Available from: https://deepmind.google/blog/alphagenome-2-deciphering-how-cells-read-dna-to-better-understand-disease/. AlphaGenome Python SDK version 2.0.0: https://github.com/google-deepmind/alphagenome (accessed 14 November 2026).",
  "86.\tSherry ST, Ward MH, Kholodov M, et al. dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001;29(1):308–311.",
  "87.\tMcLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016;17:122.",
  "88.\tChen S, Francioli LC, Goodrich JK, et al. A genomic mutational constraint map using variation in 76,156 human genomes. Nature. 2024;625:92–100.",
];
refs.forEach((r) =>
  children.push(new Paragraph({
    spacing: { before: 0, after: 100, line: 300 },
    indent: { left: 540, hanging: 540 },
    children: [new TextRun({ text: r, font: EN, size: 22 })],
  }))
);
children.push(
  fa("دربارهٔ رفرنس ۸۵: AlphaGenome 2 در ۱۱ نوامبر ۲۰۲۶ معرفی شده و من مقالهٔ داوری‌شده یا پیش‌چاپی برای آن پیدا نکردم. اگر تا زمان ارسال مقاله، گزارش فنی یا مقالهٔ رسمی آن منتشر شد، آن را جایگزین رفرنس ۸۵ کنید. رفرنس ۸۴ (مقالهٔ Nature) مرجع اصلی خانوادهٔ مدل AlphaGenome است.", { before: 160, color: "52514E" }),
  fa("رفرنس ۸۸ به نسخهٔ gnomAD اشاره دارد. اگر مجله ترجیح می‌دهد، می‌توانید به‌جای آن مقالهٔ مربوط به gnomAD v4 را بیاورید.", { color: "52514E" }),
);

// ===================== 8. Supplementary =====================
children.push(faHead("۸. جدول‌های تکمیلی (Supplementary) — اختیاری"));
children.push(
  fa("اگر مجله فایل تکمیلی می‌پذیرد، این سه جدول آماده هستند (پوشهٔ results/tables):"),
  fa("• Supplementary Table S1 — مشخصات و تأیید مختصات سه SNP (فایل table_S1_variant_annotation.csv)"),
  fa("• Supplementary Table S2 — ۱۶ بافت انتخاب‌شده و کیفیت داده‌ها (فایل table_S3_condition_selection_and_qc.csv)"),
  fa("• Supplementary Table S3 — مقایسه با همهٔ ۱۱۹۷ جایگزینی ممکن (فایل table_S5_ism_reference_distribution.csv)"),
  fa("در این صورت، در انتهای پاراگراف اول زیربخش جدید Results بنویسید: «(Supplementary Table S1)» و بعد از جملهٔ مربوط به ۹۰٪ بنویسید: «(Supplementary Table S3)»."),
);

const doc = new Document({
  styles: { default: { document: { run: { font: EN, size: 24 } } } },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1134, bottom: 1134, left: 1134, right: 1134 } } },
    children,
  }],
});

Packer.toBuffer(doc).then((b) => {
  fs.writeFileSync(require("path").join(__dirname, "..", "docs", "Manuscript_insertions_AlphaGenome2.docx"), b);
  console.log("written");
});
