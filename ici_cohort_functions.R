# =============================================================================
# ici_cohort_functions.R   (v2, 2026-11-20)
# Functions that build the harmonized multi-cohort ICI dataset
# (PRJEB23709 / PHS000452 / GSE91061 / SU2C-MARK) in the layout the
# Cohort Explorer app (appv20260823b.R) loads from multiple_cohorts_data.RData.
#
# What the app expects in df_tf (checked against Test_Table_Expr_Clinical.csv):
#   * 10 clinical columns, in this order:
#       SampleID, PatientID, BiopsyMatch, BiopsyTime, Response, OS_Status,
#       OS_Time, Prior_Tx, Regimen, StudyID
#   * then one column per gene, LINEAR TPM (the app applies log2(x+1) itself;
#     a log-scale df_tf would be logged twice)
#   * then the 22 CIBERSORT fractions (rescaled to 100) + P.value, Correlation,
#     RMSE, Absolute.score..sig.score.
#   * StudyID values: PRJEB23709, PHS000452, GSE91061, PHS000452_GSE91061
#     (pooled melanoma cohort, SampleID suffixed "_C2n5"), SU2CMARK.
#   * Any other column is treated by the app as a gene, so extra clinical
#     fields (Response_Group, Cancer, PDL1_TPS, ...) are kept in clinical_all
#     and NOT written into df_tf.
#
# Other conventions
#   * Expression matrices are genes x samples, rownames = HGNC symbol.
#   * OS_Time in DAYS, OS_Status 1 = death / 0 = censored.
#   * Response keeps the trimmed RECIST/BOR label (CR/PR/SD/PD/NE/MR), as in the
#     old table; Response_Group (PRCR/SDPD) is added in clinical_all.
#   * Every dplyr verb is namespaced: AnnotationDbi::select masks dplyr::select.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

# -----------------------------------------------------------------------------
# 0. Schema and helpers
# -----------------------------------------------------------------------------

#' The 10 clinical columns of df_tf, in the order the old table used.
CLIN_COLS <- c("SampleID", "PatientID", "BiopsyMatch", "BiopsyTime", "Response",
               "OS_Status", "OS_Time", "Prior_Tx", "Regimen", "StudyID")

#' Extra clinical fields kept in clinical_all only.
CLIN_EXTRA <- c("Cancer", "Response_Group")

CIBERSORT_CELLS <- c("B.cells.naive", "B.cells.memory", "Plasma.cells", "T.cells.CD8",
  "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
  "T.cells.follicular.helper", "T.cells.regulatory..Tregs.", "T.cells.gamma.delta",
  "NK.cells.resting", "NK.cells.activated", "Monocytes", "Macrophages.M0",
  "Macrophages.M1", "Macrophages.M2", "Dendritic.cells.resting",
  "Dendritic.cells.activated", "Mast.cells.resting", "Mast.cells.activated",
  "Eosinophils", "Neutrophils")
CIBERSORT_STATS <- c("P.value", "Correlation", "RMSE", "Absolute.score..sig.score.")

POOLED_ID     <- "PHS000452_GSE91061"
POOLED_SUFFIX <- "_C2n5"

log_step <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

assert_cols <- function(df, cols, where = deparse(substitute(df))) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss)) stop(where, ": missing column(s): ", paste(miss, collapse = ", "),
                         "\n  available: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Find one column by regex (case-insensitive); errors if 0 or >1 match.
#' Spreadsheet headers with line breaks get different names from readxl vs
#' xlsx::read.xlsx ("Dead.Alive...Dead...True." vs "Dead.Alive..Dead...True."),
#' so cohort code looks columns up by pattern instead of exact name.
find_col <- function(df, pattern, where = "table") {
  hit <- grep(pattern, colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(hit) != 1) stop(where, ": expected one column matching '", pattern, "', found ",
                             length(hit), " (", paste(hit, collapse = ", "), ")", call. = FALSE)
  hit
}

#' Trimmed, upper-case RECIST label (fixes "PD " / "SD " in PRJEB23709).
clean_response <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y %in% c("", "NA")] <- NA_character_
  y
}

#' RECIST label -> PRCR / SDPD. NE -> NA unless ne_as = "SDPD" (old behaviour).
response_group <- function(x, ne_as = NA_character_) {
  x <- clean_response(x)
  dplyr::case_when(
    x %in% c("CR", "PR")       ~ "PRCR",
    x %in% c("SD", "PD", "MR") ~ "SDPD",
    x == "NE"                  ~ ne_as,
    TRUE                       ~ NA_character_)
}

#' Death indicator (logical / "TRUE" / "Dead" / 0-1) -> 0/1.
as_event01 <- function(x) {
  if (is.logical(x) || is.numeric(x)) return(as.integer(x))
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(y))
  out[y %in% c("TRUE", "1", "DEAD", "DECEASED", "YES", "Y")] <- 1L
  out[y %in% c("FALSE", "0", "ALIVE", "NO", "N")]            <- 0L
  out[grepl("^DEAD", y)] <- 1L
  out
}

#' TRUE when the patient has both a Pre and an On sample in the table given.
#' Call it on the table AFTER aligning to expression, so a patient whose second
#' biopsy has no RNA-seq is not flagged as matched.
add_biopsy_match <- function(clin) {
  clin %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(BiopsyMatch = all(c("Pre", "On") %in% BiopsyTime)) %>%
    dplyr::ungroup()
}

#' Order/type the clinical table: CLIN_COLS, CLIN_EXTRA, then anything else.
finalize_clinical <- function(clin) {
  for (col in setdiff(c(CLIN_COLS, CLIN_EXTRA), colnames(clin))) clin[[col]] <- NA
  clin %>%
    dplyr::mutate(SampleID = as.character(SampleID), PatientID = as.character(PatientID),
                  OS_Time = as.numeric(OS_Time), OS_Status = as.integer(OS_Status),
                  BiopsyMatch = as.logical(BiopsyMatch)) %>%
    dplyr::select(dplyr::all_of(c(CLIN_COLS, CLIN_EXTRA)), dplyr::everything()) %>%
    as.data.frame()
}

# -----------------------------------------------------------------------------
# 1. Gene annotation
# -----------------------------------------------------------------------------

#' Protein-coding gene table (ENSEMBL, ENTREZID, SYMBOL) from org.Hs.eg.db.
get_protein_coding_genes <- function() {
  stopifnot(requireNamespace("org.Hs.eg.db", quietly = TRUE),
            requireNamespace("AnnotationDbi", quietly = TRUE))
  db  <- org.Hs.eg.db::org.Hs.eg.db
  ens <- AnnotationDbi::keys(db, keytype = "ENSEMBL")
  info <- suppressMessages(AnnotationDbi::select(
    db, keys = ens, columns = c("ENTREZID", "SYMBOL", "GENETYPE"), keytype = "ENSEMBL"))
  info[info$GENETYPE == "protein-coding" & !is.na(info$SYMBOL), ]
}

#' Gene length = width of the union of exons (featureCounts style), GRCh38.
get_gene_lengths_exon_union <- function() {
  for (p in c("TxDb.Hsapiens.UCSC.hg38.knownGene", "GenomicFeatures", "org.Hs.eg.db", "AnnotationDbi"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Install Bioconductor package: ", p)
  txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  exons <- GenomicFeatures::exonsBy(txdb, by = "gene")
  len   <- sum(BiocGenerics::width(GenomicRanges::reduce(exons)))
  df <- data.frame(EntrezID = names(len), Length_bp = as.numeric(len))
  df$Gene_Symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = df$EntrezID,
                                          column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  df %>%
    dplyr::filter(!is.na(Gene_Symbol)) %>%
    dplyr::arrange(Gene_Symbol, dplyr::desc(Length_bp)) %>%
    dplyr::distinct(Gene_Symbol, .keep_all = TRUE) %>%
    dplyr::mutate(Length_kb = Length_bp / 1000) %>%
    dplyr::select(Gene_Symbol, Length_bp, Length_kb)
}

# -----------------------------------------------------------------------------
# 2. Expression-matrix utilities
# -----------------------------------------------------------------------------

looks_logged <- function(mat) {
  v <- as.numeric(mat[is.finite(mat)])
  mean(v < 0) > 0.01 || stats::quantile(v, 0.999, na.rm = TRUE) < 30
}

check_scale <- function(mat, expect = c("linear", "log"), name = "matrix") {
  expect <- match.arg(expect)
  lg <- looks_logged(mat)
  if (expect == "linear" && lg)
    stop(name, " looks log-transformed (>1% negative values or 99.9th pct < 30). ",
         "Back-transform (2^x - 1) before treating it as TPM.", call. = FALSE)
  if (expect == "log" && !lg) warning(name, " does not look log-transformed.", call. = FALSE)
  invisible(lg)
}

#' Collapse duplicated gene symbols.
#'  "mean" = what the old notebooks did (group_by(Symbol) %>% summarise_all(mean));
#'           use it to reproduce the old table exactly.
#'  "sum"  = keeps each sample's TPM total at 1e6.
#'  "max_mean" = keep the most highly expressed row.
collapse_duplicate_genes <- function(mat, genes, fun = c("mean", "sum", "max_mean")) {
  fun <- match.arg(fun)
  keep <- !is.na(genes) & genes != ""
  mat <- as.matrix(mat[keep, , drop = FALSE]); genes <- genes[keep]
  storage.mode(mat) <- "double"
  if (!anyDuplicated(genes)) { rownames(mat) <- genes; return(mat) }
  n_dup <- sum(duplicated(genes))
  if (fun == "max_mean") {
    o <- order(-rowMeans(mat)); mat <- mat[o, , drop = FALSE]; genes <- genes[o]
    out <- mat[!duplicated(genes), , drop = FALSE]; rownames(out) <- genes[!duplicated(genes)]
  } else {
    out <- rowsum(mat, genes, reorder = TRUE)
    if (fun == "mean") out <- out / as.vector(table(genes)[rownames(out)])
  }
  log_step("collapsed ", n_dup, " duplicated gene rows (", fun, ")")
  out
}

#' Counts / CPM / size-factor-normalised counts -> TPM (linear input only).
counts_to_tpm <- function(mat, lengths_kb) {
  mat <- as.matrix(mat)
  stopifnot(nrow(mat) == length(lengths_kb), all(lengths_kb > 0, na.rm = TRUE))
  check_scale(mat, "linear", "counts/CPM input")
  rpk <- mat / lengths_kb
  sweep(rpk, 2, colSums(rpk, na.rm = TRUE), "/") * 1e6
}

log2tpm <- function(mat) log2(mat + 1)

summarize_matrix <- function(mat, name) {
  cs <- colSums(mat, na.rm = TRUE)
  data.frame(Matrix = name, Genes = nrow(mat), Samples = ncol(mat),
             Logged = looks_logged(mat), Median = signif(stats::median(mat, na.rm = TRUE), 4),
             ColSum_min = signif(min(cs), 4), ColSum_max = signif(max(cs), 4),
             Pct_NA = round(100 * mean(is.na(mat)), 2))
}

# -----------------------------------------------------------------------------
# 3. Expression loaders (genes x samples, linear TPM, protein-coding symbols)
# -----------------------------------------------------------------------------

#' PRJEB23709 (Gide et al.): normalizedCounts by gene symbol -> TPM.
#' The old table's PRJEB23709 values are ~4x the other cohorts (median 34 vs
#' 6-9, IRF2 46 vs 7-23); check_scale() stops here if the file is log-scale.
load_expr_PRJEB23709 <- function(path, gene_lengths, coding_symbols, dup_fun = "mean") {
  raw <- readr::read_delim(path, show_col_types = FALSE)
  assert_cols(raw, "Gene", "PRJEB23709 expression")
  raw <- raw %>% dplyr::select(-dplyr::any_of(c("ID", "Description")))
  mat <- collapse_duplicate_genes(as.matrix(raw[, -1, drop = FALSE]), raw$Gene, dup_fun)
  has_len <- rownames(mat) %in% gene_lengths$Gene_Symbol
  log_step("PRJEB23709: ", sum(!has_len), " genes dropped (no exon length)")
  mat <- mat[has_len, , drop = FALSE]
  len <- gene_lengths$Length_kb[match(rownames(mat), gene_lengths$Gene_Symbol)]
  tpm <- counts_to_tpm(mat, len)
  tpm[rownames(tpm) %in% coding_symbols, , drop = FALSE]
}

#' PHS000452 (Liu et al.): samples x genes, first column = sample ID.
load_expr_PHS000452 <- function(path, coding_symbols, input_is_log = NULL, dup_fun = "mean") {
  raw <- utils::read.delim(path, check.names = FALSE)          # keeps "HLA-A"
  ids <- as.character(raw[[1]])
  mat <- t(as.matrix(raw[, -1, drop = FALSE])); colnames(mat) <- ids
  if (is.null(input_is_log)) input_is_log <- looks_logged(mat)
  if (input_is_log) { log_step("PHS000452: input looks logged -> 2^x - 1"); mat <- 2^mat - 1 }
  mat <- collapse_duplicate_genes(mat, rownames(mat), dup_fun)
  mat[rownames(mat) %in% coding_symbols, , drop = FALSE]
}

#' GSE91061 (Riaz et al.): NCBI-generated TPM keyed by Entrez GeneID.
load_expr_GSE91061 <- function(expr_path, annot_path, dup_fun = "mean") {
  raw  <- utils::read.delim(expr_path, check.names = FALSE)
  anno <- utils::read.delim(annot_path, check.names = FALSE)
  assert_cols(raw, "GeneID"); assert_cols(anno, c("GeneID", "Symbol", "GeneType"))
  anno <- anno[anno$GeneType == "protein-coding", ]
  sym  <- anno$Symbol[match(raw$GeneID, anno$GeneID)]
  mat  <- as.matrix(raw[, setdiff(colnames(raw), "GeneID"), drop = FALSE])
  collapse_duplicate_genes(mat, sym, dup_fun)
}

#' SU2C-MARK (Ravi et al.): RNA-SeQC TPM keyed by versioned ENSEMBL ID.
#' Accepts the GCT (SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct) or the
#' supplementary xlsx (sheet Table_S13_RNA_TPM). The two hold the same values:
#' the old table's SU2CMARK genes equal the GCT exactly.
load_expr_SU2CMARK <- function(path, coding_genes, dup_fun = "mean",
                               sheet = "Table_S13_RNA_TPM", skip = 2) {
  raw <- if (grepl("\\.gct$", path, ignore.case = TRUE)) {
    utils::read.delim(path, skip = 2, check.names = FALSE)
  } else {
    as.data.frame(readxl::read_excel(path, sheet = sheet, skip = skip))
  }
  assert_cols(raw, "Name", "SU2C-MARK expression")
  # "ENSG...12_PAR_Y" rows (Y copies of pseudo-autosomal genes, all ~0) keep
  # their suffix so they do not match and are dropped; stripping it would
  # average CD99, CSF2RA, IL3RA ... with a zero row and halve them.
  ens <- sub("\\.\\d+$", "", raw$Name)
  sym <- coding_genes$SYMBOL[match(ens, coding_genes$ENSEMBL)]
  mat <- as.matrix(raw[, setdiff(colnames(raw), c("Name", "Description")), drop = FALSE])
  collapse_duplicate_genes(mat, sym, dup_fun)
}

# -----------------------------------------------------------------------------
# 4. Clinical harmonizers (one row per RNA-seq sample)
# -----------------------------------------------------------------------------

#' GSE91061: clinical xlsx + GEO mapping (PatientID like "Pt1_Pre_AD101148-6").
clin_GSE91061 <- function(clin_xlsx, mapping_xlsx, ne_as = NA_character_) {
  clin <- as.data.frame(readxl::read_excel(clin_xlsx))
  map  <- as.data.frame(readxl::read_excel(mapping_xlsx))
  assert_cols(map, c("PatientID", "geo_accession"), "GSE91061 mapping")
  dead_col <- find_col(clin, "^Dead", "GSE91061 clinical")
  time_col <- find_col(clin, "^Time.to.Death", "GSE91061 clinical")
  assert_cols(clin, c("Patient", "Cohort", "Response"), "GSE91061 clinical")
  parts <- do.call(rbind, strsplit(map$PatientID, "_"))
  map$Patient <- parts[, 1]; map$BiopsyTime <- parts[, 2]
  dt <- map %>% dplyr::select(geo_accession, Patient, BiopsyTime) %>%
    dplyr::inner_join(clin, by = "Patient")
  dt %>%
    dplyr::transmute(
      SampleID       = geo_accession,
      PatientID      = Patient,
      BiopsyTime     = dplyr::case_when(toupper(BiopsyTime) == "PRE" ~ "Pre",
                                        toupper(BiopsyTime) == "ON"  ~ "On"),
      Response       = clean_response(Response),
      OS_Status      = as_event01(.data[[dead_col]]),
      OS_Time        = round(as.numeric(.data[[time_col]]) * 7, 6),   # weeks -> days
      Prior_Tx       = dplyr::case_when(Cohort == "NIV3-NAIVE" ~ "Ipi_Naive",
                                        Cohort == "NIV3-PROG"  ~ "Ipi_Exp"),
      Regimen        = "PD1",
      StudyID        = "GSE91061",
      Cancer         = "Melanoma",
      Response_Group = response_group(Response, ne_as)) %>%
    add_biopsy_match() %>% finalize_clinical()
}

#' PHS000452: Liu et al. supplementary table.
#' prior_def  "ctla4"      Ipi_Exp if priorCTLA4 == 1 (Liu et al.)
#'            "ctla4_mapk" Ipi_Exp if prior CTLA4 OR prior MAPK therapy (old notebook)
#'            "any"        Ipi_Exp if numPriorTherapies > 0
#' All PHS000452 patients received anti-PD-1 monotherapy (Tx = Nivolumab /
#' Pembrolizumab). postCombinedCTLA_PD1 is a SUBSEQUENT therapy, so Regimen is
#' "PD1" for everyone (the old table was right; v1 of this file was not).
clin_PHS000452 <- function(clin_xlsx, prior_def = c("ctla4", "ctla4_mapk", "any"),
                           ne_as = NA_character_) {
  prior_def <- match.arg(prior_def)
  clin <- as.data.frame(readxl::read_excel(clin_xlsx))
  ctx_col <- find_col(clin, "^biopsyContext", "PHS000452 clinical")
  assert_cols(clin, c("SampleID", "OS", "dead", "BR", "priorCTLA4", "priorMAPKTx",
                      "numPriorTherapies"), "PHS000452 clinical")
  clin %>%
    dplyr::filter(!is.na(OS)) %>%
    dplyr::transmute(
      SampleID       = SampleID,
      PatientID      = SampleID,
      # biopsyContext 1 = Pre-Ipi, 2 = On-Ipi, 3 = Pre-PD1, 4 = On-PD1
      BiopsyTime     = ifelse(as.numeric(.data[[ctx_col]]) < 4, "Pre", "On"),
      Response       = clean_response(BR),
      OS_Status      = as_event01(dead),
      OS_Time        = as.numeric(OS),
      Prior_Tx       = switch(prior_def,
                         ctla4      = ifelse(priorCTLA4 == 1, "Ipi_Exp", "Ipi_Naive"),
                         ctla4_mapk = ifelse(priorCTLA4 == 0 & priorMAPKTx == 0, "Ipi_Naive", "Ipi_Exp"),
                         any        = ifelse(numPriorTherapies > 0, "Ipi_Exp", "Ipi_Naive")),
      Regimen        = "PD1",
      StudyID        = "PHS000452",
      Cancer         = "Melanoma",
      Response_Group = response_group(BR, ne_as)) %>%
    add_biopsy_match() %>% finalize_clinical()
}

#' PRJEB23709: Gide et al. clinical table (Table_ClinicalData.tsv).
#' The tsv lists only PRE samples, but the expression file also holds 18
#' on-treatment "*_EDT" samples. Pass expr_samples = colnames(expression) so:
#'  * include_edt = FALSE (default, = old table): EDT samples stay out, and
#'    BiopsyMatch = TRUE for PRE samples whose patient has an EDT sample
#'    (reproduces the old table's 16 flags exactly);
#'  * include_edt = TRUE: EDT samples are added as BiopsyTime "On" rows that
#'    carry their patient's outcome, and BiopsyMatch flags both rows.
#' PatientID keeps the arm ("PD1_26" vs "ipiPD1_26"): the bare number used in
#' the old table merges two different patients from the two arms.
#' prior_label: the old table used "Tx_Naive"; the patients are ipi-naive.
clin_PRJEB23709 <- function(clin_tsv, prior_label = "Tx_Naive", expr_samples = NULL,
                            include_edt = FALSE, ne_as = NA_character_) {
  clin <- as.data.frame(readr::read_tsv(clin_tsv, show_col_types = FALSE))
  alive_col <- find_col(clin, "^Alive", "PRJEB23709 clinical")
  assert_cols(clin, c("Sample_id", "Timing", "RECIST", "Overall_survival", "ICI_target"),
              "PRJEB23709 clinical")
  out <- clin %>%
    dplyr::transmute(
      SampleID       = Sample_id,
      PatientID      = sub("_(PRE|EDT|ON)$", "", Sample_id, ignore.case = TRUE),
      BiopsyTime     = ifelse(toupper(Timing) == "PRE", "Pre", "On"),
      Response       = clean_response(RECIST),
      OS_Status      = 1L - as.integer(.data[[alive_col]]),   # Alive = 1 -> event 0
      OS_Time        = as.numeric(Overall_survival),
      Prior_Tx       = prior_label,
      Regimen        = ifelse(ICI_target == "PD1", "PD1", "Combo"),
      StudyID        = "PRJEB23709",
      Cancer         = "Melanoma",
      Response_Group = response_group(RECIST, ne_as)) %>%
    add_biopsy_match()
  if (!is.null(expr_samples)) {
    edt    <- grep("_EDT$", expr_samples, value = TRUE, ignore.case = TRUE)
    pre_of <- sub("_EDT$", "_PRE", edt, ignore.case = TRUE)
    keep   <- pre_of %in% out$SampleID
    if (any(!keep)) log_step("PRJEB23709: ", sum(!keep), " EDT samples without a PRE clinical row ignored")
    edt <- edt[keep]; pre_of <- pre_of[keep]
    if (include_edt) {
      on <- out[match(pre_of, out$SampleID), ]
      on$SampleID <- edt; on$BiopsyTime <- "On"
      out <- add_biopsy_match(dplyr::bind_rows(out, on))
    } else {
      out$BiopsyMatch <- out$SampleID %in% pre_of
    }
  }
  finalize_clinical(out)
}

#' SU2C-MARK. Accepts either
#'  * the processed table Processed_Table_ClinicalData_SUCMARK.csv
#'    (SampleID, QC, TNM, PDL1_TPS, Line_of_Therapy, ..., PD1_Category, BORcb,
#'    OS_Status, OS_Time, Tx) -- 124 RNA-seq patients, used for the old table; or
#'  * the raw Table_S1_Clinical_Annotations.xlsx (same 124 after RNA_All == 1).
#' qc_keep_only = TRUE drops samples whose Pre-treatment RNA QC is "Flag"
#' (12 of 124); the old table kept them.
clin_SU2CMARK <- function(path, qc_keep_only = FALSE, ne_as = NA_character_, skip = 2) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    d <- utils::read.csv(path, check.names = FALSE)
  } else {
    raw <- as.data.frame(readxl::read_excel(path, sheet = 1, skip = skip))
    assert_cols(raw, c("Harmonized_SU2C_RNA_Tumor_Sample_ID_v2", "RNA_All",
                       "Pre-treatment_RNA_Sample_QC", "Initial_Stage", "PDL1_TPS",
                       "Line_of_Therapy", "Prior_Platinum", "Prior_TKI", "Agent_PD1_Category",
                       "Harmonized_Confirmed_BOR", "Harmonized_OS_Event", "Harmonized_OS_Days"),
                "SU2C-MARK Table S1")
    d <- raw %>% dplyr::filter(RNA_All == 1) %>%
      dplyr::transmute(SampleID = Harmonized_SU2C_RNA_Tumor_Sample_ID_v2,
                       QC = `Pre-treatment_RNA_Sample_QC`, TNM = Initial_Stage,
                       PDL1_TPS, Line_of_Therapy, Prior_Platinum, Prior_TKI,
                       PD1_Category = Agent_PD1_Category, BORcb = Harmonized_Confirmed_BOR,
                       OS_Status = Harmonized_OS_Event, OS_Time = Harmonized_OS_Days,
                       Tx = ifelse(Line_of_Therapy > 1, "Tx_Exp", "Tx_Naive"))
  }
  assert_cols(d, c("SampleID", "QC", "PD1_Category", "BORcb", "OS_Status", "OS_Time", "Tx"),
              "SU2C-MARK clinical")
  if (qc_keep_only) d <- d[d$QC == "Keep", , drop = FALSE]
  d %>%
    dplyr::transmute(
      SampleID       = SampleID,
      PatientID      = SampleID,
      BiopsyMatch    = FALSE,
      BiopsyTime     = "Pre",
      Response       = clean_response(BORcb),
      OS_Status      = as_event01(OS_Status),
      OS_Time        = as.numeric(OS_Time),
      Prior_Tx       = Tx,
      Regimen        = ifelse(PD1_Category == "PD(L)1", "PD1", "Combo"),
      StudyID        = "SU2CMARK",
      Cancer         = "NSCLC",
      Response_Group = response_group(BORcb, ne_as),
      PD1_Category, QC, PDL1_TPS = suppressWarnings(as.numeric(PDL1_TPS)),
      Line_of_Therapy) %>%
    finalize_clinical()
}

# -----------------------------------------------------------------------------
# 5. Validation, alignment, patient-level selection
# -----------------------------------------------------------------------------

validate_clinical <- function(clin, name = unique(clin$StudyID)[1]) {
  assert_cols(clin, CLIN_COLS, name)
  if (anyDuplicated(clin$SampleID)) stop(name, ": duplicated SampleID(s): ",
    paste(head(unique(clin$SampleID[duplicated(clin$SampleID)])), collapse = ", "), call. = FALSE)
  issues <- character()
  if (any(!clin$OS_Status %in% c(0L, 1L, NA))) issues <- c(issues, "OS_Status not 0/1")
  if (any(clin$OS_Time < 0, na.rm = TRUE))     issues <- c(issues, "negative OS_Time")
  if (any(is.na(clin$BiopsyTime)))             issues <- c(issues, "BiopsyTime NA")
  if (any(is.na(clin$Regimen)))                issues <- c(issues, "Regimen NA")
  med <- stats::median(clin$OS_Time, na.rm = TRUE)
  if (is.finite(med) && (med < 60 || med > 3650))
    issues <- c(issues, sprintf("median OS_Time = %.0f - check units are days", med))
  if (length(issues)) warning(name, ": ", paste(issues, collapse = "; "), call. = FALSE)
  invisible(issues)
}

#' Keep samples present in both tables; recompute BiopsyMatch on what is left.
align_cohort <- function(expr, clin) {
  name <- unique(clin$StudyID)[1]
  common <- intersect(clin$SampleID, colnames(expr))
  log_step(sprintf("%s: %d clinical, %d expression, %d matched (dropped %d clin-only, %d expr-only)",
                   name, nrow(clin), ncol(expr), length(common),
                   nrow(clin) - length(common), ncol(expr) - length(common)))
  clin <- clin[match(common, clin$SampleID), , drop = FALSE]
  if (!all(clin$BiopsyTime == "Pre")) clin <- finalize_clinical(add_biopsy_match(clin))
  list(expr = expr[, common, drop = FALSE], clin = clin)
}

#' One sample per patient (Pre preferred). Patients with only On biopsies are
#' dropped unless keep_on_only = TRUE.
select_one_per_patient <- function(clin, prefer = "Pre", keep_on_only = FALSE) {
  out <- clin %>%
    dplyr::mutate(.pref = BiopsyTime == prefer) %>%
    dplyr::arrange(StudyID, PatientID, dplyr::desc(.pref), SampleID) %>%
    dplyr::group_by(StudyID, PatientID) %>% dplyr::slice(1) %>% dplyr::ungroup()
  if (!keep_on_only) out <- out %>% dplyr::filter(.pref)
  as.data.frame(dplyr::select(out, -.pref))
}

# -----------------------------------------------------------------------------
# 6. Batch correction and the pooled PHS000452 + GSE91061 cohort
# -----------------------------------------------------------------------------

combine_matrices <- function(mat_list) {
  genes <- Reduce(intersect, lapply(mat_list, rownames))
  log_step("combine: ", length(genes), " shared genes across ", length(mat_list), " matrices")
  ids <- unlist(lapply(mat_list, colnames))
  if (anyDuplicated(ids)) stop("SampleID collision across cohorts: ",
                               paste(head(ids[duplicated(ids)]), collapse = ", "))
  do.call(cbind, lapply(mat_list, function(m) m[genes, , drop = FALSE]))
}

#' Rows that ComBat can handle (non-constant in every batch).
combat_ok_rows <- function(log_mat, batch, min_sd = 1e-3) {
  Reduce(`&`, lapply(split(seq_along(batch), batch), function(ix)
    apply(log_mat[, ix, drop = FALSE], 1, stats::sd) > min_sd))
}

batch_correct_combat <- function(log_mat, batch, covars = NULL, par_prior = TRUE) {
  stopifnot(requireNamespace("sva", quietly = TRUE), length(batch) == ncol(log_mat))
  check_scale(log_mat, "log", "ComBat input")
  mod <- if (is.null(covars)) NULL else stats::model.matrix(~ ., data = covars)
  sva::ComBat(dat = log_mat, batch = batch, mod = mod, par.prior = par_prior)
}

#' PC variance explained by batch (R^2); near 0 after correction = batch removed.
batch_variance_explained <- function(log_mat, batch, n_pc = 5, n_top = 2000) {
  v <- apply(log_mat, 1, stats::var)
  top <- log_mat[order(-v)[seq_len(min(n_top, length(v)))], , drop = FALSE]
  pc <- stats::prcomp(t(top), center = TRUE, scale. = FALSE)
  n_pc <- min(n_pc, ncol(pc$x))
  ve <- (pc$sdev^2 / sum(pc$sdev^2))[seq_len(n_pc)]
  r2 <- sapply(seq_len(n_pc), function(i) summary(stats::lm(pc$x[, i] ~ factor(batch)))$r.squared)
  data.frame(PC = paste0("PC", seq_len(n_pc)), VarExplained = round(ve, 3), R2_batch = round(r2, 3))
}

#' Build the pooled "PHS000452_GSE91061" cohort the app shows as "Cohort 1".
#' method = "combat_seq" (default, = the old table): sva::ComBat_seq() on the
#'            combined PHS000452 + GSE91061 TPM matrix (genes shared by the two),
#'            batch = StudyID, no group. Re-running this on the old table's own
#'            uncorrected rows gives per-gene correlation 0.9998 with its pooled
#'            rows and the same residual cohort offset (0.27 log2); exact
#'            equality needs the same gene set, because ComBat_seq's dispersion
#'            estimate is shared across genes. Output is integer-valued.
#'            Caveat: ComBat_seq models counts; on TPM it rounds low values
#'            (TPM < 0.5 -> 0) and only halves the cohort offset.
#'          "combat": ComBat on log2(TPM+1), back-transformed to TPM.
#'          "none"  : plain stacking.
#' Returns list(expr = linear genes x samples, clin, diagnostics).
make_pooled_cohort <- function(expr, clin, studies = c("PHS000452", "GSE91061"),
                               method = c("combat_seq", "combat", "none"),
                               covars = NULL, group = NULL) {
  method <- match.arg(method)
  lin <- combine_matrices(expr[studies])
  cl  <- dplyr::bind_rows(clin[studies])
  cl  <- cl[match(colnames(lin), cl$SampleID), , drop = FALSE]
  batch <- cl$StudyID
  off <- function(m) round(stats::median(abs(apply(log2tpm(m), 1, function(x)
    mean(x[batch == studies[1]]) - mean(x[batch == studies[2]])))), 3)
  diag <- list(offset_before = off(lin))
  if (method == "combat_seq") {
    stopifnot(requireNamespace("sva", quietly = TRUE))
    grp <- if (is.null(group)) NULL else cl[[group]]
    lin <- suppressMessages(sva::ComBat_seq(counts = lin, batch = batch, group = grp))
  } else if (method == "combat") {
    lg  <- log2tpm(lin)
    ok  <- combat_ok_rows(lg, batch)
    log_step("pooled ComBat: correcting ", sum(ok), " / ", length(ok), " genes")
    cov_df <- if (is.null(covars)) NULL else cl[, covars, drop = FALSE]
    lg[ok, ] <- batch_correct_combat(lg[ok, , drop = FALSE], batch, cov_df)
    lin <- pmax(2^lg - 1, 0)
  }
  diag$offset_after <- off(lin)
  log_step("pooled ", method, ": per-gene cohort offset (log2) ", diag$offset_before, " -> ", diag$offset_after)
  colnames(lin) <- paste0(colnames(lin), POOLED_SUFFIX)
  cl$SampleID <- paste0(cl$SampleID, POOLED_SUFFIX)
  cl$StudyID  <- POOLED_ID
  list(expr = lin, clin = cl, diagnostics = diag)
}

# -----------------------------------------------------------------------------
# 7. CIBERSORT in / out
# -----------------------------------------------------------------------------

#' Write the CIBERSORT mixture file (Symbol + one column per sample, linear TPM,
#' LM22 genes shared by all cohorts). Reproduces the old Mixture_ImmuneRegions.txt
#' values exactly. Use QN = FALSE for RNA-seq TPM: QN = TRUE (meant for arrays)
#' quantile-normalizes across ALL mixture samples, so a sample's fractions change
#' whenever the sample set changes, and it makes the absolute score a constant.
write_cibersort_mixture <- function(expr, lm22_symbols, path) {
  genes <- Reduce(intersect, c(lapply(expr, rownames), list(lm22_symbols)))
  mix <- do.call(cbind, lapply(expr, function(m) m[genes, , drop = FALSE]))
  out <- data.frame(Symbol = genes, mix, check.names = FALSE)
  utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
  log_step("CIBERSORT mixture: ", length(genes), " genes x ", ncol(mix), " samples -> ", path)
  invisible(genes)
}

#' Read CIBERSORT results (data.frame returned by CIBERSORT(), or the
#' CIBERSORT-results.txt it writes), rescale the 22 fractions to sum to 100
#' (as normal100() did) and return SampleID + 22 cells + 4 stats.
read_cibersort <- function(x) {
  res <- if (is.character(x)) utils::read.delim(x, check.names = FALSE) else as.data.frame(x)
  if (!"SampleID" %in% colnames(res)) {
    if (colnames(res)[1] %in% c("Mixture", "Input Sample", "")) colnames(res)[1] <- "SampleID"
    else res <- tibble::rownames_to_column(res, "SampleID")
  }
  colnames(res) <- make.names(colnames(res))      # "T cells CD8" -> "T.cells.CD8", "P-value" -> "P.value"
  colnames(res)[colnames(res) == "SampleID"] <- "SampleID"
  assert_cols(res, c("SampleID", CIBERSORT_CELLS), "CIBERSORT results")
  for (s in setdiff(CIBERSORT_STATS, colnames(res))) res[[s]] <- NA_real_
  frac <- as.matrix(res[, CIBERSORT_CELLS])
  tot  <- rowSums(frac)
  # absolute = TRUE with QN = FALSE gives all-zero rows when a sample's median
  # over the signature genes is 0 (7 of 473 in the 2026 mixture) -> NA, not NaN
  if (any(tot == 0)) warning(sum(tot == 0), " CIBERSORT rows sum to 0 (set to NA): ",
                             paste(head(res$SampleID[tot == 0]), collapse = ", "), call. = FALSE)
  tot[tot == 0] <- NA
  res[, CIBERSORT_CELLS] <- frac / tot * 100
  res[, c("SampleID", CIBERSORT_CELLS, CIBERSORT_STATS)]
}

# -----------------------------------------------------------------------------
# 8. The app table
# -----------------------------------------------------------------------------

#' Build df_tf: one row per sample, CLIN_COLS + genes (linear TPM) + CIBERSORT.
#' expr: named list of linear genes x samples matrices (pooled cohort included)
#' clin: data.frame with CLIN_COLS for every sample in expr
#' cibersort: output of read_cibersort() or NULL (columns then left out; the app
#'            falls back to the last 22 columns, so leave them out only for testing)
#' genes: restrict to these; default = genes shared by every matrix
make_app_table <- function(expr, clin, cibersort = NULL, genes = NULL) {
  if (is.null(genes)) genes <- Reduce(intersect, lapply(expr, rownames))
  log_step("df_tf: ", length(genes), " genes")
  ex <- do.call(cbind, lapply(expr, function(m) m[genes, , drop = FALSE]))
  if (anyDuplicated(colnames(ex))) stop("duplicated SampleID across expression matrices")
  ex <- as.data.frame(t(ex), check.names = FALSE) %>% tibble::rownames_to_column("SampleID")
  out <- clin[, CLIN_COLS] %>% dplyr::inner_join(ex, by = "SampleID")
  if (!is.null(cibersort)) {
    # exact SampleID first; pooled rows then fall back to their source sample
    # (old table: CIBERSORT ran on per-cohort TPM) and single-cohort rows to
    # their pooled copy (CIBERSORT run on the pooled matrix instead)
    idx <- match(out$SampleID, cibersort$SampleID)
    alt <- ifelse(grepl(paste0(POOLED_SUFFIX, "$"), out$SampleID),
                  sub(paste0(POOLED_SUFFIX, "$"), "", out$SampleID), paste0(out$SampleID, POOLED_SUFFIX))
    idx[is.na(idx)] <- match(alt[is.na(idx)], cibersort$SampleID)
    cs  <- cibersort[idx, c(CIBERSORT_CELLS, CIBERSORT_STATS)]
    if (any(is.na(cs$T.cells.CD8))) warning(sum(is.na(cs$T.cells.CD8)), " samples have no CIBERSORT result")
    out <- cbind(out, cs)
  }
  out$BiopsyMatch <- as.logical(out$BiopsyMatch)
  out
}

#' Cross-tab of the clinical table.
cohort_summary <- function(clin) {
  clin %>%
    dplyr::group_by(StudyID, BiopsyTime, Prior_Tx, Regimen) %>%
    dplyr::summarise(Samples = dplyr::n(), Patients = dplyr::n_distinct(PatientID),
                     Deaths = sum(OS_Status == 1, na.rm = TRUE),
                     PRCR = sum(Response %in% c("CR", "PR")),
                     Median_OS_days = stats::median(OS_Time, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
}

# -----------------------------------------------------------------------------
# 9. Audit / compare against an existing table
# -----------------------------------------------------------------------------

#' Integrity report for a df_tf-style table.
audit_wide_table <- function(df, n_clin = 10) {
  cl <- df[, seq_len(n_clin)]
  num <- setdiff(colnames(df)[-seq_len(n_clin)], c(CIBERSORT_CELLS, CIBERSORT_STATS))
  g  <- as.matrix(df[, num])
  by_study <- split(seq_len(nrow(df)), cl$StudyID)
  overview <- cl %>%
    dplyr::group_by(StudyID) %>%
    dplyr::summarise(Rows = dplyr::n(), Unique_samples = dplyr::n_distinct(SampleID),
                     Patients = dplyr::n_distinct(PatientID),
                     Patients_multi_sample = sum(table(PatientID) > 1),
                     On_biopsies = sum(BiopsyTime == "On", na.rm = TRUE),
                     Median_OS_days = stats::median(OS_Time, na.rm = TRUE),
                     Response_whitespace = sum(grepl("^\\s|\\s$", Response)), .groups = "drop")
  scale_tbl <- do.call(rbind, lapply(names(by_study), function(s) {
    v <- g[by_study[[s]], , drop = FALSE]
    data.frame(StudyID = s, Median = stats::median(v, na.rm = TRUE), Max = max(v, na.rm = TRUE),
               IRF2_median = if ("IRF2" %in% colnames(v)) stats::median(v[, "IRF2"]) else NA,
               Looks_logged = looks_logged(v))
  }))
  cat("== Overview\n"); print(as.data.frame(overview))
  cat("\n== Expression scale (genes only)\n"); print(scale_tbl)
  cat("\n== Non-gene columns the app would treat as genes:",
      paste(setdiff(colnames(df)[-seq_len(n_clin)][!vapply(df[-seq_len(n_clin)], is.numeric, TRUE)], ""),
            collapse = ", "), "\n")
  invisible(list(overview = overview, scale = scale_tbl))
}

#' Compare a new df_tf to the old one, sample by sample.
#' Returns per-column agreement for the clinical fields and, per study, the
#' share of gene values that match within `tol` (relative).
compare_to_reference <- function(new, ref, tol = 1e-3, genes = NULL) {
  ref <- ref[!duplicated(paste(ref$StudyID, ref$SampleID)), ]
  key_n <- paste(new$StudyID, new$SampleID); key_r <- paste(ref$StudyID, ref$SampleID)
  common <- intersect(key_r, key_n)
  n <- new[match(common, key_n), ]; r <- ref[match(common, key_r), ]
  clin_cmp <- do.call(rbind, lapply(setdiff(CLIN_COLS, c("SampleID", "StudyID")), function(cc) {
    a <- as.character(n[[cc]]); b <- trimws(as.character(r[[cc]]))
    if (cc == "OS_Time") { eq <- abs(as.numeric(a) - as.numeric(b)) < 0.5 } else eq <- a == b
    data.frame(Column = cc, StudyID = n$StudyID, Equal = eq)
  })) %>%
    dplyr::group_by(StudyID, Column) %>%
    dplyr::summarise(Agree = paste0(sum(Equal, na.rm = TRUE), "/", dplyr::n()), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Column, values_from = Agree) %>% as.data.frame()
  # The old table was round-tripped through read.csv(), which turned every "-"
  # in a gene name into "." (HLA-A -> HLA.A). Match on make.names() of the new names.
  new_names <- stats::setNames(colnames(new), make.names(colnames(new)))
  if (is.null(genes)) genes <- intersect(names(new_names), colnames(ref))
  genes <- setdiff(genes, c(CLIN_COLS, CIBERSORT_CELLS, CIBERSORT_STATS))
  if (!length(genes)) return(list(matched_samples = length(common), clinical = clin_cmp))
  expr_cmp <- do.call(rbind, lapply(split(seq_len(nrow(n)), n$StudyID), function(ix) {
    a <- as.matrix(n[ix, new_names[genes]]); b <- as.matrix(r[ix, genes])
    rel <- abs(a - b) / pmax(abs(b), 1)
    data.frame(StudyID = n$StudyID[ix[1]], Samples = length(ix), Genes = length(genes),
               Pct_equal = round(100 * mean(rel < tol, na.rm = TRUE), 2),
               Median_ratio_new_over_ref = round(stats::median((a + 1) / (b + 1), na.rm = TRUE), 3))
  }))
  list(matched_samples = length(common), clinical = clin_cmp, expression = expr_cmp,
       genes_only_in_ref = setdiff(colnames(ref), c(names(new_names), CIBERSORT_CELLS, CIBERSORT_STATS)),
       genes_only_in_new = unname(new_names[setdiff(names(new_names), colnames(ref))]))
}
