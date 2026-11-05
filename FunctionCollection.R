# ==============================================================================
# Standalone Kaplan-Meier Plotting Script (for R Notebook use)
# Adapted from: appv20260709.R (Melanoma and NSCLC Cohorts Explorer)
#
# This script reproduces the app's core single-factor K-M plotting logic
# (High vs Low split on a continuous variable, e.g. gene expression or cell
# fraction) as a standalone function, outside of Shiny.
#
# EXPECTED INPUT TABLE FORMAT (data.frame / tibble), one row per sample:
#   SampleID   : character/factor, unique sample identifier
#   OS_Time    : numeric, time-to-event in DAYS
#   OS_Status  : numeric/integer, 0 = censored, 1 = event
#   <value_col>: numeric, the variable used to split samples into High/Low
#                (e.g. a gene's expression, or a cell type fraction 0-1)
#
# Optional:
#   A pre-existing categorical grouping column (e.g. "Response", "Tx") can
#   be supplied instead of/alongside a continuous split — see
#   plot_km_by_group() below.
# ==============================================================================

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. GENERIC HELPERS  (unchanged conventions from appv20260709.R)
# ------------------------------------------------------------------------------

# Normalize a continuous variable prior to splitting.
# IMPORTANT (per app conventions): never apply this to cell-type fraction
# variables that are already on a 0-1 scale.
normalize_vector <- function(raw_vec, norm_type = c("none", "log2", "zscore", "robust")) {
  norm_type <- match.arg(norm_type)
  switch(norm_type,
         log2   = log2(raw_vec + 1),
         zscore = as.numeric(scale(raw_vec, center = TRUE, scale = TRUE)),
         robust = (raw_vec - median(raw_vec, na.rm = TRUE)) / IQR(raw_vec, na.rm = TRUE),
         none   = as.numeric(raw_vec))
}

# Administrative censoring at a max follow-up (in years).
apply_admin_censoring <- function(df, max_years) {
  max_time_days <- max_years * 365.25
  df %>% dplyr::mutate(
    OS_Status = ifelse(OS_Time > max_time_days, 0, OS_Status),
    OS_Time   = ifelse(OS_Time > max_time_days, max_time_days, OS_Time))
}

# X-axis breaks/labels in whole years, for plots plotted on a day scale.
year_axis_breaks <- function(max_years) {
  max_time_days <- max_years * 365.25
  breaks_days   <- seq(0, max_time_days, by = 365.25)
  list(days = breaks_days, labels = round(breaks_days / 365.25), max_days = max_time_days)
}

# Split a cohort into High / Low using a symmetric tail-cut percentage.
# split_pct = 50 -> median split. split_pct = 25 -> top 25% vs bottom 25%
# (middle 50% dropped).
assign_global_group <- function(df, val_col, split_pct) {
  q_high <- quantile(df[[val_col]], 1 - (split_pct / 100), na.rm = TRUE)
  q_low  <- quantile(df[[val_col]], split_pct / 100,       na.rm = TRUE)
  df_temp <- df %>% dplyr::mutate(
    Group = dplyr::case_when(
      .data[[val_col]] > q_high ~ "High",
      .data[[val_col]] <= q_low  ~ "Low",
      TRUE ~ "Middle"))
  df_out       <- df_temp %>% dplyr::filter(Group != "Middle")
  df_out$Group <- factor(df_out$Group, levels = c("Low", "High"))
  return(df_out)
}

# ------------------------------------------------------------------------------
# 2. CORE K-M PLOT FUNCTION (continuous variable -> High/Low split)
#    Mirrors the Feature 3 "Factor-Based Profiling" single-factor K-M plot,
#    including: tryCatch safety, ordered-unnamed-palette convention (avoids
#    color misalignment vs. named vectors + legend.labs), and a caption with
#    median OS / landmark OS rate.
# ------------------------------------------------------------------------------

#' Draw a High vs Low Kaplan-Meier plot for a continuous variable.
#'
#' @param df            data.frame with SampleID, OS_Time (days), OS_Status, value_col
#' @param value_col     name of the continuous column to split on (string)
#' @param split_pct     tail-cut percentage (5-50; 50 = median split)
#' @param norm_type     "none" | "log2" | "zscore" | "robust" — applied to value_col
#'                      before splitting. Set to "none" for 0-1 scale variables
#'                      such as immune cell fractions.
#' @param max_years     administrative censoring / x-axis cutoff, in years
#' @param show_metrics  "both" | "median" | "landmark" | "none"
#' @param landmark_years landmark year(s) used for the landmark OS rate caption
#' @param title         plot title (defaults to value_col name)
#' @param plot_title_suffix optional extra text appended to the title (e.g. cohort name)
#'
#' @details Risk table is disabled. Each subgroup's median survival time is
#'   marked with a vertical dashed reference line dropping to the x-axis
#'   (via ggsurvplot's \code{surv.median.line = "v"}), colored to match its
#'   stratum. The log-rank p-value is shown with the total N (after
#'   High/Low cutoff filtering) directly beneath it.
#'
#' @return a survminer ggsurvplot object (invisibly prints it); NULL if the
#'         plot could not be generated (e.g. insufficient samples/events).
plot_km_high_low <- function(df,
                             value_col,
                             split_pct       = 50,
                             norm_type       = c("none", "log2", "zscore", "robust"),
                             max_years       = NULL,
                             show_metrics    = c("both", "median", "landmark", "none"),
                             landmark_years  = 1,
                             title           = NULL,
                             plot_title_suffix = NULL) {
  
  norm_type    <- match.arg(norm_type)
  show_metrics <- match.arg(show_metrics)
  
  req_cols <- c("SampleID", "OS_Time", "OS_Status", value_col)
  missing_cols <- setdiff(req_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input table is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  df <- df %>% dplyr::filter(!is.na(.data[[value_col]]), !is.na(OS_Time), !is.na(OS_Status))
  if (nrow(df) < 5) {
    message("Insufficient samples (< 5) for K-M plot.")
    return(NULL)
  }
  
  # Default max_years to the observed max follow-up if not supplied
  if (is.null(max_years)) {
    max_years <- ceiling(max(df$OS_Time, na.rm = TRUE) / 365.25)
  }
  df <- apply_admin_censoring(df, max_years)
  
  # Normalize (skip for pre-scaled fractions by using norm_type = "none")
  df$Expr_Val <- normalize_vector(df[[value_col]], norm_type)
  
  # High/Low split
  df_split <- tryCatch(assign_global_group(df, "Expr_Val", split_pct), error = function(e) NULL)
  if (is.null(df_split) || nrow(df_split) < 5 || length(unique(df_split$Group)) < 2 ||
      length(unique(df_split$OS_Status)) < 2) {
    message("Insufficient group separation or event variation for K-M plot.")
    return(NULL)
  }
  
  plot_title <- if (is.null(title)) paste("Survival by", value_col) else title
  if (!is.null(plot_title_suffix)) plot_title <- paste0(plot_title, "\n", plot_title_suffix)
  
  result <- tryCatch({
    km_fit       <- survminer::surv_fit(Surv(OS_Time, OS_Status) ~ Group, data = df_split)
    # Ordered UNNAMED palette derived from the fitted strata order — avoids
    # color/legend misalignment that can occur with named vectors + legend.labs.
    strata_order <- gsub("Group=", "", names(km_fit$strata))
    km_pal       <- unname(ifelse(strata_order == "High", "#E41A1C", "#377EB8"))
    
    caption_lines <- c()
    if (show_metrics %in% c("median", "both")) {
      ms <- survminer::surv_median(km_fit)
      caption_lines <- c(caption_lines, paste(
        "Median OS \u2192",
        paste(sapply(seq_len(nrow(ms)), function(i) {
          lbl <- gsub("Group=", "", ms$strata[i])
          val <- ifelse(is.na(ms$median[i]), "NE", paste(round(ms$median[i] / 365.25, 2), "yrs"))
          paste0(lbl, ": ", val)
        }), collapse = " | ")))
    }
    if (show_metrics %in% c("landmark", "both")) {
      ls <- summary(km_fit, times = landmark_years * 365.25)
      if (length(ls$strata) > 0) {
        caption_lines <- c(caption_lines, paste0(
          landmark_years, "-Year OS \u2192 ",
          paste(sapply(seq_along(ls$strata), function(i) {
            lbl <- gsub("Group=", "", ls$strata[i])
            paste0(lbl, ": ", round(ls$surv[i] * 100, 1), "%")
          }), collapse = " | ")))
      }
    }
    caption_text <- if (length(caption_lines) > 0) paste(caption_lines, collapse = "\n") else NULL
    
    # Log-rank p-value + N (after High/Low cutoff filtering), stacked as a
    # single annotation so N appears directly under the p-value.
    sdiff    <- tryCatch(survdiff(Surv(OS_Time, OS_Status) ~ Group, data = df_split), error = function(e) NULL)
    pval_num <- if (!is.null(sdiff)) pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE) else NA
    pval_lbl <- paste0("p = ", ifelse(is.na(pval_num), "NA", format.pval(pval_num, digits = 3)),
                       "\nN = ", nrow(df_split))
    
    p <- survminer::ggsurvplot(
      km_fit, data = df_split, pval = pval_lbl, risk.table = FALSE,
      surv.median.line = "v",
      title = plot_title, legend.title = value_col,
      palette = km_pal, legend.labs = strata_order, ggtheme = theme_bw(),
      xscale = 'd_y', break.time.by = 365.25, xlim = c(0, max_years * 365.25))
    
    ax <- year_axis_breaks(max_years)
    p$plot  <- p$plot  + scale_x_continuous(breaks = ax$days, labels = ax$labels) +
      labs(x = "Time (Years)", caption = caption_text) +
      theme(plot.caption = element_text(face = "bold", size = 11, hjust = 1, color = "grey30"))
    if (!is.null(p$table)) {
      p$table <- p$table + scale_x_continuous(breaks = ax$days, labels = ax$labels) +
        labs(x = "Time (Years)")
    }
    p
  }, error = function(e) {
    message("K-M plot error: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(result)) {
    plot.new()
    title(main = paste0("Plot unavailable\n(", value_col, ")"), col.main = "grey50", cex.main = 1)
    return(invisible(NULL))
  }
  
  print(result)
  invisible(result)
}

# ------------------------------------------------------------------------------
# 3. K-M PLOT FROM AN EXISTING CATEGORICAL GROUPING COLUMN
#    Use this when the table already has a group/subgroup column (e.g.
#    "Response", "Tx", "StudyID") rather than a continuous variable to split.
# ------------------------------------------------------------------------------

#' Draw a Kaplan-Meier plot for an existing categorical grouping column.
#'
#' @param df           data.frame with SampleID, OS_Time (days), OS_Status, group_col
#' @param group_col    name of the categorical column to stratify by (string)
#' @param max_years    administrative censoring / x-axis cutoff, in years
#' @param show_metrics "both" | "median" | "landmark" | "none"
#' @param landmark_years landmark year(s) used for the landmark OS rate caption
#' @param title        plot title
#' @param palette      color palette name or vector passed to ggsurvplot (default "Dark2")
plot_km_by_group <- function(df,
                             group_col,
                             max_years       = NULL,
                             show_metrics    = c("both", "median", "landmark", "none"),
                             landmark_years  = 1,
                             title           = NULL,
                             palette         = "Dark2") {
  
  show_metrics <- match.arg(show_metrics)
  req_cols <- c("SampleID", "OS_Time", "OS_Status", group_col)
  missing_cols <- setdiff(req_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input table is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  df <- df %>% dplyr::filter(!is.na(.data[[group_col]]), !is.na(OS_Time), !is.na(OS_Status))
  df[[group_col]] <- droplevels(factor(df[[group_col]]))
  if (nrow(df) < 5 || nlevels(df[[group_col]]) < 2) {
    message("Insufficient samples or groups (< 2) for K-M plot.")
    return(invisible(NULL))
  }
  
  if (is.null(max_years)) max_years <- ceiling(max(df$OS_Time, na.rm = TRUE) / 365.25)
  df <- apply_admin_censoring(df, max_years)
  
  plot_title <- if (is.null(title)) paste("Survival by", group_col) else title
  
  result <- tryCatch({
    frm <- as.formula(paste0("Surv(OS_Time, OS_Status) ~ ", group_col))
    km_fit <- survminer::surv_fit(frm, data = df)
    
    sub_text <- "Log-Rank p-value: NA (need >= 2 groups)"
    if (nlevels(df[[group_col]]) >= 2) {
      sdiff <- tryCatch(survdiff(frm, data = df), error = function(e) NULL)
      if (!is.null(sdiff)) {
        pv <- pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)
        sub_text <- paste0("Log-Rank p-value: ", format.pval(pv, digits = 3))
      }
    }
    
    caption_lines <- c()
    if (show_metrics %in% c("median", "both")) {
      ms <- survminer::surv_median(km_fit)
      caption_lines <- c(caption_lines, paste(
        "Median OS \u2192",
        paste(sapply(seq_len(nrow(ms)), function(i) {
          lbl <- gsub(paste0(group_col, "="), "", ms$strata[i])
          val <- ifelse(is.na(ms$median[i]), "NE", paste(round(ms$median[i] / 365.25, 2), "yrs"))
          paste0(lbl, ": ", val)
        }), collapse = " | ")))
    }
    if (show_metrics %in% c("landmark", "both")) {
      ls <- summary(km_fit, times = landmark_years * 365.25)
      if (length(ls$strata) > 0) {
        caption_lines <- c(caption_lines, paste0(
          landmark_years, "-Year OS \u2192 ",
          paste(sapply(seq_along(ls$strata), function(i) {
            lbl <- gsub(paste0(group_col, "="), "", ls$strata[i])
            paste0(lbl, ": ", round(ls$surv[i] * 100, 1), "%")
          }), collapse = " | ")))
      }
    }
    caption_text <- if (length(caption_lines) > 0) paste(caption_lines, collapse = "\n") else NULL
    
    p <- survminer::ggsurvplot(
      km_fit, data = df, pval = FALSE, risk.table = TRUE, risk.table.height = 0.3,
      legend = "bottom", legend.title = group_col,
      title = plot_title, subtitle = sub_text, palette = palette,
      xscale = 'd_y', break.time.by = 365.25, xlim = c(0, max_years * 365.25))
    
    ax <- year_axis_breaks(max_years)
    p$plot  <- p$plot  + scale_x_continuous(breaks = ax$days, labels = ax$labels) +
      labs(x = "Time (Years)", caption = caption_text) +
      theme(plot.subtitle = element_text(size = 13, face = "italic", color = "grey30"),
            plot.caption  = element_text(face = "bold", size = 11, hjust = 1,
                                         color = "grey30", lineheight = 1.2))
    p$table <- p$table + scale_x_continuous(breaks = ax$days, labels = ax$labels) +
      labs(x = "Time (Years)")
    p
  }, error = function(e) {
    message("K-M plot error: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(result)) {
    plot.new()
    title(main = paste0("Plot unavailable\n(", group_col, ")"), col.main = "grey50", cex.main = 1)
    return(invisible(NULL))
  }
  
  print(result)
  invisible(result)
}

# ------------------------------------------------------------------------------
# 5. MULTI-CURVE TAIL-CUT K-M PLOT (e.g. Top 50% vs Top 25% vs Top 10%)
#    Mirrors the app's Feature 3 "K-M Survival Curves (Multi Tail-Cuts)"
#    subtab: build_tailcut_groups() + plot_tailcut_km(). Draws one curve per
#    cut percentage (all drawn from the same "side" of the distribution, so
#    e.g. Top 10% is a strict subset of Top 25%, which is a subset of Top
#    50%) and reports pairwise log-rank p-values for every pair of cuts.
# ------------------------------------------------------------------------------

# Build the stacked long-format data for a set of tail-cut percentages on one
# side ("High" = top X%, "Low" = bottom X%) of a (normalized) score column.
build_tailcut_groups <- function(df, score_col, cut_pcts, side = c("High", "Low")) {
  side     <- match.arg(side)
  cut_pcts <- sort(unique(cut_pcts), decreasing = TRUE)
  df_list  <- list()
  for (pct in cut_pcts) {
    df_g <- tryCatch(assign_global_group(df, score_col, pct), error = function(e) NULL)
    if (is.null(df_g)) next
    df_side <- df_g %>% dplyr::filter(Group == side)
    if (nrow(df_side) == 0) next
    lbl              <- paste0(ifelse(side == "High", "Top ", "Bottom "), pct, "%")
    df_side$CutLabel <- lbl
    df_list[[lbl]]   <- df_side
  }
  if (length(df_list) == 0) return(NULL)
  lbl_levels <- names(df_list)
  out <- dplyr::bind_rows(df_list)
  out$CutLabel <- factor(out$CutLabel, levels = lbl_levels)
  out
}

#' Draw a multi-curve K-M plot comparing several tail-cut percentages of a
#' gene (or any continuous variable), e.g. Top 50% vs Top 25% vs Top 10%.
#'
#' @param df            data.frame with SampleID, OS_Time (days), OS_Status, value_col
#' @param value_col     name of the continuous column to cut on, e.g. a gene (string)
#' @param cut_pcts      numeric vector of tail-cut percentages to overlay,
#'                      e.g. c(50, 25, 10) for Top 50% / Top 25% / Top 10%.
#'                      Fully customizable — pass any set of percentages (1-50).
#' @param side          "High" (top X%, default) or "Low" (bottom X%)
#' @param norm_type     "none" | "log2" | "zscore" | "robust" — applied to value_col
#'                      before cutting. Use "none" for 0-1 scale variables
#'                      such as immune cell fractions.
#' @param max_years     administrative censoring / x-axis cutoff, in years
#' @param title         plot title (defaults to value_col name)
#'
#' @details Risk table is disabled. The plot caption shows the pairwise
#'   log-rank p-value for every pair of cut groups shown, all on a single
#'   line (e.g. for cut_pcts = c(50, 25, 10): 50% vs 25%, 50% vs 10%, and
#'   25% vs 10%), with significance stars (*p<0.05, **p<0.01, ***p<0.001).
#'   Colors are an ordered, unnamed palette so the tightest cut (e.g. Top
#'   10%) is always the darkest shade.
#'
#' @return a survminer ggsurvplot object (invisibly prints it); NULL if the
#'         plot could not be generated (e.g. insufficient samples/events).
plot_km_tailcut <- function(df,
                            value_col,
                            cut_pcts        = c(50, 25, 10),
                            side            = c("High", "Low"),
                            norm_type       = c("none", "log2", "zscore", "robust"),
                            max_years       = NULL,
                            title           = NULL) {
  
  side         <- match.arg(side)
  norm_type    <- match.arg(norm_type)
  
  req_cols <- c("SampleID", "OS_Time", "OS_Status", value_col)
  missing_cols <- setdiff(req_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input table is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  if (length(cut_pcts) < 1) stop("cut_pcts must contain at least one percentage.")
  if (any(cut_pcts <= 0 | cut_pcts > 50)) stop("cut_pcts must be in (0, 50].")
  
  df <- df %>% dplyr::filter(!is.na(.data[[value_col]]), !is.na(OS_Time), !is.na(OS_Status))
  if (nrow(df) < 5) {
    message("Insufficient samples (< 5) for K-M plot.")
    return(invisible(NULL))
  }
  
  if (is.null(max_years)) max_years <- ceiling(max(df$OS_Time, na.rm = TRUE) / 365.25)
  df <- apply_admin_censoring(df, max_years)
  
  df$Score_Val <- normalize_vector(df[[value_col]], norm_type)
  
  df_grp <- build_tailcut_groups(df, "Score_Val", cut_pcts, side)
  if (is.null(df_grp) || nlevels(droplevels(df_grp$CutLabel)) < 1) {
    message("Insufficient samples in one or more cut groups.")
    return(invisible(NULL))
  }
  df_grp$CutLabel <- droplevels(df_grp$CutLabel)
  cut_labels <- levels(df_grp$CutLabel)
  
  plot_title <- if (is.null(title)) paste0(value_col, " \u2014 ",
                                           ifelse(side == "High", "Top-% Cuts", "Bottom-% Cuts")) else title
  
  result <- tryCatch({
    fit <- survfit(Surv(OS_Time, OS_Status) ~ CutLabel, data = df_grp)
    
    # Pairwise log-rank p-values across every pair of displayed cut groups,
    # stacked into a single caption line with significance stars
    # (*p<0.05, **p<0.01, ***p<0.001).
    caption_text <- "Pairwise Log-Rank: NA (need >= 2 groups)"
    if (length(cut_labels) >= 2) {
      pairs <- combn(cut_labels, 2, simplify = FALSE)
      pval_parts <- sapply(pairs, function(pr) {
        df_pair <- df_grp %>% dplyr::filter(CutLabel %in% pr)
        sdiff <- tryCatch(survdiff(Surv(OS_Time, OS_Status) ~ CutLabel, data = df_pair), error = function(e) NULL)
        pv <- if (!is.null(sdiff)) pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE) else NA
        signif_stars <- dplyr::case_when(
          is.na(pv)   ~ "",
          pv < 0.001  ~ "***",
          pv < 0.01   ~ "**",
          pv < 0.05   ~ "*",
          TRUE        ~ "")
        paste0(pr[1], " vs ", pr[2], ": ", ifelse(is.na(pv), "NA", format.pval(pv, digits = 3)), signif_stars)
      })
      caption_text <- paste0("\u2014", paste(pval_parts, collapse = " | "))
    }
    
    # N per cut group, appended after the title so each curve's sample size is visible
    n_by_cut  <- table(df_grp$CutLabel)
    n_lbl     <- paste(paste0(names(n_by_cut), ": N=", as.integer(n_by_cut)), collapse = " | ")
    plot_title_full <- paste0(plot_title, "\n", n_lbl)
    
    # Ordered, unnamed palette: tighter cuts get darker shades of the same hue
    base_pal <- if (side == "High") c("#a50f15", "#fb6a4a", "#fcbba1", "#fee0d2", "#fff5f0") else
      c("#08519c", "#6baed6", "#c6dbef", "#deebf7", "#f7fbff")
    pal <- unname(base_pal[seq_along(cut_labels)])
    
    p <- survminer::ggsurvplot(
      fit, data = df_grp, pval = FALSE, risk.table = FALSE,
      surv.median.line = NULL,
      palette = pal, legend.labs = cut_labels, legend = "top",
      title    = plot_title_full,
      ggtheme  = theme_bw(),
      xscale   = 'd_y', break.time.by = 365.25, xlim = c(0, max_years * 365.25))
    
    ax <- year_axis_breaks(max_years)
    p$plot <- p$plot + scale_x_continuous(breaks = ax$days, labels = ax$labels) +
      labs(x = "Time (Years)", caption = caption_text) +
      theme(plot.caption = element_text(face = "bold", size = 11, hjust = 1,
                                        color = "grey30", lineheight = 1.2))
    if (!is.null(p$table)) {
      p$table <- p$table + scale_x_continuous(breaks = ax$days, labels = ax$labels) +
        labs(x = "Time (Years)")
    }
    p
  }, error = function(e) {
    message("K-M tail-cut plot error: ", conditionMessage(e))
    NULL
  })
  
  if (is.null(result)) {
    plot.new()
    title(main = paste0("Plot unavailable\n(", value_col, ")"), col.main = "grey50", cex.main = 1)
    return(invisible(NULL))
  }
  
  print(result)
  invisible(result)
}

# ==============================================================================
# 6. EXAMPLE USAGE (run in an R notebook chunk)
# ==============================================================================
# --- Example A: continuous variable (e.g. gene expression), median split ---
# set.seed(1)
# n <- 200
# example_df <- data.frame(
#   SampleID  = paste0("PT-", 1:n),
#   OS_Time   = sample(100:1500, n, replace = TRUE),
#   OS_Status = sample(0:1, n, replace = TRUE),
#   IRF1      = rnorm(n, 5, 2)
# )
# plot_km_high_low(example_df, value_col = "IRF1", split_pct = 50,
#                   norm_type = "log2", show_metrics = "both", landmark_years = 1)
#
# --- Example B: immune cell fraction (0-1 scale), top/bottom 25% ---
# example_df$T.cells.CD8 <- runif(n, 0, 0.4)
# plot_km_high_low(example_df, value_col = "T.cells.CD8", split_pct = 25,
#                   norm_type = "none")
#
# --- Example C: existing categorical grouping column ---
# example_df$Response <- sample(c("Responder", "Non-Responder"), n, replace = TRUE)
# plot_km_by_group(example_df, group_col = "Response", show_metrics = "both")
#
# --- Example D: three-curve tail-cut comparison (Top 50% / 25% / 10%) ---
# plot_km_tailcut(example_df, value_col = "IRF1", cut_pcts = c(50, 25, 10),
#                  side = "High", norm_type = "log2")
#
# --- Example E: fully customizable cut percentages, e.g. Top 40%/20%/5% ---
# plot_km_tailcut(example_df, value_col = "IRF1", cut_pcts = c(40, 20, 5),
#                  side = "High", norm_type = "log2")
# ==============================================================================

# ==============================================================================
# Standalone helper functions
#   1) plot_km_nested_subset()  — K-M plot of Gene B (Top vs Bottom split),
#                                  restricted to patients in the Top/Bottom
#                                  split of Gene A. Mirrors the app's own
#                                  nested-stratification logic (Feature 3).
#   2) plot_gene_expr_barplot() — Cross-cohort gene expression barplot with
#                                  pairwise Wilcoxon p-value brackets
#
# DEPENDENCIES (already defined in appv20260709.R — NOT redefined here):
#   - assign_global_group(df, val_col, split_pct)
#   - assign_group_subgroup(df, val_col, split_pct, subgroup_col = NULL)
#   - year_axis_breaks(max_years)
#   - normalize_vector(raw_vec, norm_type)
# These functions must be sourced from appv20260709.R (or pasted above this
# file) before calling either function below.
#
# Written to match the conventions already used in appv20260709.R:
#   - tryCatch-wrapped survfit()/ggsurvplot() (KM rendering safety)
#   - unnamed, strata-ordered palette vectors (avoids ggsurvplot color misalignment)
#   - `[` instead of `[[` for named-vector lookups (safer NA on mismatch)
#   - normalize_vector() is applied per-gene, never to cell-fraction (0-1) columns
#   - year_axis_breaks() for consistent x-axis year labeling
#   - caption built from surv_median()/summary(fit, times=) the same way as
#     plot_subset_km() / plot_tailcut_km()
# ==============================================================================

# ------------------------------------------------------------------------------
# 1) plot_km_nested_subset()
#
# Nested two-gene stratification, same idea as prep_nested_stratification() +
# plot_subset_km() in appv20260709.R, returned as a single self-contained
# plotting function with THREE curves:
#
#   Stage 1 (Gene A): split the WHOLE cohort into Top splitA% / Bottom splitA%
#                      by geneA, using the app's assign_group_subgroup().
#   Subset:            keep only the requested Gene-A group (default: "High",
#                      i.e. patients in the TOP splitA% of Gene A).
#   Stage 2 (Gene B): within that Gene-A subset, split by geneB into
#                      Top splitB% / Bottom splitB%, again via
#                      assign_group_subgroup().
#   Plot:              3 curves --
#                        - Gene-B Top splitB%    (red,  solid)
#                        - Gene-B Bottom splitB% (blue, solid)
#                        - Entire Gene-A subset  (black, dashed) -- i.e. ALL
#                          patients in the top/bottom splitA% of Gene A,
#                          regardless of their Gene-B status. These patients
#                          overlap with the Top-B/Bottom-B curves by design,
#                          same overlay pattern as the "Overall (Reference)"
#                          curve in plot_subset_km().
#                      Log-rank p-value is computed ONLY between Gene-B Top
#                      vs. Bottom (the black/reference curve is excluded from
#                      the test).
#
# Args:
#   df             - data frame with OS_Time, OS_Status, geneA, geneB columns
#                     (already admin-censored via apply_admin_censoring() by
#                     the caller, consistent with plot_subset_km/plot_tailcut_km)
#   geneA, geneB   - column names of the two genes
#   splitA, splitB - split percentages (e.g. 50 for a 50/50 median split)
#   subsetA        - which Gene-A group to restrict to: "High" (top splitA%,
#                     the default — matches "patients in the top 50% of gene A")
#                     or "Low" (bottom splitA%)
#   subgroup_col   - optional column for per-cohort/subgroup-wise splitting,
#                     passed straight through to assign_group_subgroup();
#                     NULL = split across the whole supplied cohort (global)
#   norm_type      - "none" | "log2" | "zscore" | "robust", applied to geneA
#                     and geneB before splitting (via normalize_vector()).
#                     Normalization is SKIPPED automatically for whichever of
#                     geneA/geneB is listed in `cell_types` (already a 0-1
#                     scale fraction), same rule used in
#                     plot_cell_fraction_boxplot()/get_cox_zscore_pval() --
#                     this lets you pass a cell-fraction column as geneB (or
#                     geneA) without it being log2/z-score transformed.
#   cell_types     - character vector of cell-fraction column names, used
#                     only to decide which of geneA/geneB to skip
#                     normalization for; default none (normalize whichever
#                     of geneA/geneB norm_type applies to)
#   max_years, show_metrics, landmark_years - same as the other KM helpers
# ------------------------------------------------------------------------------
plot_km_nested_subset <- function(df, geneA, geneB, splitA = 50, splitB = 50,
                                  subsetA = c("High", "Low"), subgroup_col = NULL,
                                  norm_type = c("none", "log2", "zscore", "robust"),
                                  cell_types = character(0),
                                  max_years, show_metrics = "both", landmark_years = 1) {
  subsetA   <- match.arg(subsetA)
  norm_type <- match.arg(norm_type)
  
  if (is.null(df) || !all(c(geneA, geneB) %in% colnames(df))) return(NULL)
  df <- df %>%
    dplyr::filter(!is.na(.data[[geneA]]), !is.na(.data[[geneB]]),
                  !is.na(OS_Time), !is.na(OS_Status))
  if (nrow(df) < 6) return(NULL)
  
  # normalize geneA/geneB unless each is itself a cell-fraction column
  # (already on a 0-1 scale) -- same rule plot_cell_fraction_boxplot() applies
  if (!(geneA %in% cell_types)) df[[geneA]] <- normalize_vector(df[[geneA]], norm_type)
  if (!(geneB %in% cell_types)) df[[geneB]] <- normalize_vector(df[[geneB]], norm_type)
  
  # --- Stage 1: split the whole cohort by Gene A, keep only the requested arm ---
  df_A   <- assign_group_subgroup(df, geneA, splitA, subgroup_col) %>%
    dplyr::rename(GroupA = Group)
  df_sub <- df_A %>% dplyr::filter(GroupA == subsetA)
  if (nrow(df_sub) < 6) return(NULL)
  
  # --- Stage 2: within the Gene-A subset, split by Gene B ---
  df_stage2 <- assign_group_subgroup(df_sub, geneB, splitB, subgroup_col) %>%
    dplyr::rename(GroupB = Group)
  if (length(unique(df_stage2$GroupB)) < 2) return(NULL)
  
  a_lbl     <- paste0(if (subsetA == "High") "Top " else "Bottom ", splitA, "% of ", geneA)
  top_b_lbl <- paste0("Top ",    splitB, "% ", geneB)
  bot_b_lbl <- paste0("Bottom ", splitB, "% ", geneB)
  ref_lbl   <- paste0("All (", a_lbl, ")")
  
  # Build the three display groups:
  #   - Top-B / Bottom-B: the Gene-B split within the Gene-A subset
  #   - ref_lbl: the ENTIRE Gene-A subset (regardless of Gene-B status) --
  #     patients here legitimately overlap with Top-B/Bottom-B, same overlay
  #     pattern as the "Overall (Reference)" curve in plot_subset_km().
  df_stage2 <- df_stage2 %>%
    dplyr::mutate(PlotGroup = dplyr::if_else(GroupB == "High", top_b_lbl, bot_b_lbl))
  df_ref  <- df_sub %>% dplyr::mutate(PlotGroup = ref_lbl)
  df_plot <- dplyr::bind_rows(df_stage2, df_ref)
  
  group_levels     <- c(top_b_lbl, ref_lbl, bot_b_lbl)
  df_plot$PlotGroup <- factor(df_plot$PlotGroup,
                              levels = intersect(group_levels, unique(df_plot$PlotGroup)))
  
  max_time_days <- max_years * 365.25
  
  # --- log-rank test: Gene-B High vs Low only (ref/black curve excluded from the test) ---
  df_tb <- df_stage2 %>% dplyr::filter(PlotGroup %in% c(top_b_lbl, bot_b_lbl))
  sdiff <- tryCatch(survdiff(Surv(OS_Time, OS_Status) ~ PlotGroup, data = df_tb),
                    error = function(e) NULL)
  pval  <- if (!is.null(sdiff)) pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE) else NA
  
  fit <- tryCatch(survfit(Surv(OS_Time, OS_Status) ~ PlotGroup, data = df_plot),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  # --- p-value + N annotation, same style/placement as plot_km_high_low():
  # "p = <value>\nN = <value>", using ggsurvplot's default pval position
  # (no pval.coord override). N = total sample size actually used in the
  # log-rank comparison (Top-B + Bottom-B; the overlapping reference/black
  # curve is excluded, same as it's excluded from the log-rank test itself). ---
  pval_lbl <- paste0("p = ", ifelse(is.na(pval), "NA", format.pval(pval, digits = 3)),
                     "\nN = ", nrow(df_tb))
  
  # --- caption: median OS / landmark OS rate only ---
  caption_lines <- c()
  if (show_metrics %in% c("median", "both")) {
    ms <- survminer::surv_median(fit)
    caption_lines <- c(caption_lines, paste("Median OS \u2192", paste(
      sapply(seq_len(nrow(ms)), function(i) {
        sg  <- gsub("PlotGroup=", "", ms$strata[i])
        val <- ms$median[i]
        paste0(sg, ": ", ifelse(is.na(val), "NE", paste(round(val / 365.25, 2), "yrs")))
      }), collapse = " | ")))
  }
  if (show_metrics %in% c("landmark", "both")) {
    land_summ <- summary(fit, times = landmark_years * 365.25)
    if (length(land_summ$strata) > 0) {
      caption_lines <- c(caption_lines, paste0(landmark_years, "-Year OS Rate \u2192 ", paste(
        sapply(seq_along(land_summ$strata), function(i) {
          sg <- gsub("PlotGroup=", "", land_summ$strata[i])
          paste0(sg, ": ", round(land_summ$surv[i] * 100, 1), "%")
        }), collapse = " | ")))
    } else {
      caption_lines <- c(caption_lines, paste0(landmark_years, "-Year OS Rate \u2192 NE"))
    }
  }
  caption_text <- if (length(caption_lines) > 0) paste(caption_lines, collapse = "\n") else NULL
  
  # --- palette/linetype: strata-order-safe lookup via `[`, not `[[` ---
  strata_order <- gsub("PlotGroup=", "", names(fit$strata))
  col_map <- c(setNames("#E41A1C", top_b_lbl), setNames("black", ref_lbl), setNames("#377EB8", bot_b_lbl))
  lty_map <- c(setNames("solid",   top_b_lbl), setNames("dashed", ref_lbl), setNames("solid",   bot_b_lbl))
  pal_ordered <- unname(col_map[strata_order])
  lty_ordered <- unname(lty_map[strata_order])
  legend_labs <- strata_order  # plain group names -- N is shown in-plot instead
  
  p <- tryCatch(
    ggsurvplot(fit, data = df_plot, risk.table = TRUE, risk.table.height = 0.3,
               pval = pval_lbl,
               palette = pal_ordered, linetype = lty_ordered,
               legend.labs = legend_labs, legend = "top",
               title    = paste0(geneB, " Survival within ", a_lbl),
               xscale = 'd_y', break.time.by = 365.25, xlim = c(0, max_time_days)),
    error = function(e) NULL)
  if (is.null(p)) return(NULL)
  
  ax <- year_axis_breaks(max_years)
  p$plot  <- p$plot  + scale_x_continuous(breaks = ax$days, labels = ax$labels) +
    labs(x = "Time (Years)", caption = caption_text) +
    theme(plot.caption = element_text(face = "bold", size = 11, hjust = 1,
                                      color = "grey30", lineheight = 1.2))
  p$table <- p$table + scale_x_continuous(breaks = ax$days, labels = ax$labels) +
    labs(x = "Time (Years)")
  
  return(p)
}


# ------------------------------------------------------------------------------
# 2) plot_gene_expr_barplot()
#
# Bar graph of gene expression (chosen normalization) on the y-axis, faceted
# by gene, with selected cohorts on the x-axis. Pairwise Wilcoxon p-values
# between selected cohorts are shown as brackets above the bars, using
# ggpubr::stat_compare_means() the same way the existing cross-cohort
# "Gene Expression Barplot" subtab (f3bar_plot) does — the only differences
# here are (a) normalization is configurable rather than fixed to log2, and
# (b) p-values are printed as exact values ("p.format") rather than
# significance stars, per the request.
#
# Args:
#   df               - filtered data frame (e.g. df_filtered())
#   genes             - character vector of gene column names to plot
#   norm_type         - "none" | "log2" | "zscore" | "robust" (uses normalize_vector())
#   group_col         - column to use for the x-axis cohort grouping (e.g. "StudyID",
#                        or a pre-computed "Hierarchical_Subgroup")
#   selected_cohorts  - optional character vector to subset/restrict which
#                        cohort levels are shown & compared; NULL = use all
#                        levels present in df
# ------------------------------------------------------------------------------
plot_gene_expr_barplot <- function(df, genes, norm_type = c("none", "log2", "zscore", "robust"),
                                   group_col = "StudyID", selected_cohorts = NULL) {
  norm_type <- match.arg(norm_type)
  
  genes <- intersect(genes, colnames(df))
  if (length(genes) == 0 || !(group_col %in% colnames(df))) return(NULL)
  
  req_cols <- unique(c("SampleID", group_col, genes))
  df <- df %>% dplyr::select(dplyr::all_of(req_cols))
  
  if (!is.null(selected_cohorts)) {
    df <- df %>% dplyr::filter(.data[[group_col]] %in% selected_cohorts)
  }
  if (nrow(df) == 0) return(NULL)
  
  # Normalize per-gene (never pool across genes before normalizing, since
  # different genes/cells live on different scales — same rule as elsewhere
  # in the app for zscore/robust normalization).
  df_long <- df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(genes), names_to = "Gene", values_to = "RawExpr") %>%
    dplyr::filter(!is.na(RawExpr)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(NormExpr = normalize_vector(RawExpr, norm_type)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Gene = factor(Gene, levels = genes))
  
  df_long[[group_col]] <- droplevels(factor(df_long[[group_col]]))
  cohort_lvls <- levels(df_long[[group_col]])
  if (length(cohort_lvls) == 0) return(NULL)
  comp_list <- if (length(cohort_lvls) >= 2) combn(cohort_lvls, 2, simplify = FALSE) else NULL
  
  y_lab <- switch(norm_type,
                  none   = "Expression (raw)",
                  log2   = "Log2(Expression + 1)",
                  zscore = "Expression (Z-score)",
                  robust = "Expression (Robust: median/IQR)")
  
  p <- ggpubr::ggbarplot(df_long, x = group_col, y = "NormExpr",
                         add = c("mean_se", "jitter"), fill = group_col,
                         facet.by = "Gene", palette = "Dark2",
                         xlab = "Cohort", ylab = y_lab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  if (!is.null(comp_list)) {
    p <- p + ggpubr::stat_compare_means(comparisons = comp_list, method = "wilcox.test",
                                        label = "p.format")
  }
  
  return(p)
}

# ==============================================================================
#   3) get_cox_zscore_pval()      — Cox PH Wald Z-score + log-rank p-value for
#                                    ONE gene within ONE cohort/subset
#   4) compute_survival_heatmap_data() — assembles Gene x Cohort Z-score/p-value
#                                    grid (calls get_cox_zscore_pval per cell)
#   5) plot_survival_heatmap()    — renders the Gene x Cohort heatmap
#
# DEPENDENCY / CORRECTION NOTE:
#   These reuse assign_global_group() and normalize_vector() from
#   appv20260709.R rather than re-implementing the High/Low split. This
#   matters because assign_global_group() contains the factor-level fix
#   `factor(df_out$Group, levels = c("Low","High"))` -- forcing "Low" as the
#   Cox PH reference level so the Wald Z-score sign consistently reads as
#   "High vs. Low" (positive Z = High group has worse survival). This is the
#   same convention already used by the app's own get_km_zscore_pval() /
#   f2_heatmap_df(); re-deriving the split without it would silently risk
#   flipping the sign of every Z-score in the heatmap.
#
#   NOTE ALSO: the app applies NO multiple-testing correction anywhere
#   (no p.adjust() calls exist in appv20260709.R) -- the significance
#   asterisks below use raw p-values, matching that existing convention.
#   Flag if you'd like a Benjamini-Hochberg (FDR) correction added across
#   the Gene x Cohort grid; it isn't currently part of the app's pattern.
# ==============================================================================

# NOTE: kept as a SEPARATE, uniquely-named pair from get_km_zscore_pval() /
# compute_survival_heatmap() / plot_survival_heatmap() further below (which use
# a different Target/Subgroup_Label column schema and add HR/CI + two_factor
# support). An earlier revision here tried to consolidate these into a single
# plot_survival_heatmap() name and broke calls to compute_survival_heatmap_data()
# (this function) -- two functions sharing one name silently collide when
# sourced (last definition wins), so compute_survival_heatmap_data()'s own
# plotting function is now named plot_survival_heatmap_data() instead, to
# guarantee it can never collide with plot_survival_heatmap() again regardless
# of source order.
#
# ------------------------------------------------------------------------------
# 3) get_cox_zscore_pval()
#
# Same logic as the app's own get_km_zscore_pval() in "single" mode, just
# renamed/simplified for a plain Gene-vs-Cohort survival heatmap (no
# anchor/two-factor branch needed here).
#
# Args:
#   df_sub     - data for ONE cohort (already filtered to that cohort)
#   gene_col   - gene column name to split High/Low on
#   split_pct  - split percentage (50 = median split)
#   norm_type  - "none" | "log2" | "zscore" | "robust", applied to gene_col
#                before splitting (via normalize_vector()); "none" is the
#                right choice for already-normalized expression data
#
# Returns: list(Z, Pval, N, Events) -- Z is the Cox PH Wald Z-score for
#   Group (High vs Low), Pval is the log-rank p-value for the same split.
# ------------------------------------------------------------------------------
get_cox_zscore_pval <- function(df_sub, gene_col, split_pct = 50,
                                norm_type = c("none", "log2", "zscore", "robust")) {
  norm_type <- match.arg(norm_type)
  
  df_eval <- df_sub %>% dplyr::filter(!is.na(.data[[gene_col]]), !is.na(OS_Time), !is.na(OS_Status))
  if (nrow(df_eval) < 10) return(list(Z = NA, Pval = NA, N = NA, Events = NA))
  
  df_eval[[gene_col]] <- normalize_vector(df_eval[[gene_col]], norm_type)
  
  df_eval <- assign_global_group(df_eval, gene_col, split_pct)
  if (nrow(df_eval) < 5 || length(unique(df_eval$Group)) < 2 ||
      length(unique(df_eval$OS_Status)) < 2) {
    return(list(Z = NA, Pval = NA, N = NA, Events = NA))
  }
  
  sdiff <- tryCatch(survdiff(Surv(OS_Time, OS_Status) ~ Group, data = df_eval), error = function(e) NULL)
  pval  <- if (!is.null(sdiff)) pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE) else NA
  
  cox_fit <- tryCatch(coxph(Surv(OS_Time, OS_Status) ~ Group, data = df_eval), error = function(e) NULL)
  if (!is.null(cox_fit)) {
    scox    <- summary(cox_fit)
    z_score <- scox$coefficients[1, "z"]
  } else {
    z_score <- NA
  }
  
  events <- sum(df_eval$OS_Status == 1, na.rm = TRUE)
  list(Z = z_score, Pval = pval, N = nrow(df_eval), Events = events)
}


# ------------------------------------------------------------------------------
# 4) compute_survival_heatmap_data()
#
# Builds the Gene x Cohort grid of Cox PH Z-scores / p-values, one call to
# get_cox_zscore_pval() per (cohort, gene) cell -- mirrors the app's own
# f2_heatmap_df() loop structure and its Z-capping / asterisk / label logic.
#
# Args:
#   df          - full (pre-filtered) patient-level data frame
#   genes       - character vector of gene column names (x-axis, in this order)
#   cohort_col  - column identifying the cohort/study grouping (y-axis)
#   cohorts     - optional character vector to restrict/order which cohort
#                 levels are included; NULL = all levels present in df,
#                 alphabetically sorted
#   split_pct   - High/Low split percentage passed to get_cox_zscore_pval()
#   norm_type   - normalization passed to get_cox_zscore_pval()
#
# Returns a data frame with one row per (Cohort, Gene) cell, including:
#   Z, Pval, N, Events, Z_capped (clipped to [-5, 5]), Signif ("" / "*" / "**" /
#   "***"), Label (heatmap cell text: Z-score + newline + asterisks),
#   plus Gene/Cohort as ordered factors ready for plot_survival_heatmap_data().
# ------------------------------------------------------------------------------
compute_survival_heatmap_data <- function(df, genes, cohort_col = "StudyID", cohorts = NULL,
                                          split_pct = 50,
                                          norm_type = c("none", "log2", "zscore", "robust")) {
  norm_type <- match.arg(norm_type)
  
  genes <- intersect(genes, colnames(df))
  if (length(genes) == 0 || !(cohort_col %in% colnames(df))) return(NULL)
  
  if (!is.null(cohorts)) df <- df %>% dplyr::filter(.data[[cohort_col]] %in% cohorts)
  cohort_lvls <- if (!is.null(cohorts)) cohorts else sort(unique(as.character(df[[cohort_col]])))
  if (length(cohort_lvls) == 0) return(NULL)
  
  res_list <- list()
  for (co in cohort_lvls) {
    df_sub <- df %>% dplyr::filter(.data[[cohort_col]] == co)
    for (g in genes) {
      stats <- get_cox_zscore_pval(df_sub, g, split_pct, norm_type)
      res_list[[length(res_list) + 1]] <- data.frame(
        Cohort = co, Gene = g,
        N = stats$N, Events = stats$Events, Z = stats$Z, Pval = stats$Pval)
    }
  }
  
  dplyr::bind_rows(res_list) %>%
    dplyr::mutate(
      Z_capped = ifelse(is.na(Z), NA, pmin(pmax(Z, -5), 5)),
      Signif   = dplyr::case_when(is.na(Pval) ~ "", Pval < 0.001 ~ "***",
                                  Pval < 0.01 ~ "**", Pval < 0.05 ~ "*", TRUE ~ ""),
      Label    = ifelse(is.na(Z_capped), "NA", paste0(round(Z_capped, 2),
                                                      ifelse(Signif == "", "", paste0("\n", Signif)))),
      Gene     = factor(Gene,   levels = genes),
      Cohort   = factor(Cohort, levels = cohort_lvls))
}


# ------------------------------------------------------------------------------
# 5) plot_survival_heatmap_data()
#
# Renders the Gene x Cohort heatmap from compute_survival_heatmap_data():
# genes across the top (x-axis), cohorts down the left side (y-axis), cell
# fill = capped Z-score (-5 to 5, diverging blue-white-red), cell text =
# Z-score with significance asterisks (blank if not significant). Same
# visual convention as the app's f2_heatmap_ui (geom_tile + geom_text +
# scale_fill_gradient2), returned as a plain ggplot object here -- wrap in
# plotly::ggplotly() yourself if you want the click-to-drill-down behavior
# the app's version has.
#
# NAMED plot_survival_heatmap_data() (not plot_survival_heatmap()) specifically
# so it can never collide with the Target/Subgroup_Label plot_survival_heatmap()
# defined later in this file -- see the note at the top of this section.
#
# Args:
#   hdf         - output of compute_survival_heatmap_data() (Gene/Cohort schema)
#   title       - subtitle text shown under the plot title area
#   orientation - "landscape" (default; Gene on x-axis, Cohort on y-axis) or
#                 "portrait" (axes swapped: Cohort on x-axis, Gene on y-axis)
# ------------------------------------------------------------------------------
plot_survival_heatmap_data <- function(hdf, title = "Survival Heatmap (High vs Low)",
                                       orientation = c("landscape", "portrait")) {
  if (is.null(hdf) || nrow(hdf) == 0) return(NULL)
  orientation <- match.arg(orientation)
  x_var <- if (orientation == "landscape") "Gene" else "Cohort"
  y_var <- if (orientation == "landscape") "Cohort" else "Gene"

  ggplot(hdf, aes(x = .data[[x_var]], y = .data[[y_var]], fill = Z_capped)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Label), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                         midpoint = 0, limits = c(-5, 5), oob = scales::squish,
                         na.value = "grey80", name = "Z-score") +
    scale_y_discrete(limits = rev(levels(hdf[[y_var]]))) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
          axis.text.y = element_text(face = "bold", size = 12), panel.grid = element_blank()) +
    labs(subtitle = title, x = NULL, y = NULL)
}

# ------------------------------------------------------------------------------
# 6) plot_cell_fraction_boxplot()
#
# Template: the app's own f3_cell_fraction_ui / observe block
# (ggboxplot(x="CellType", color="Plot_Group") + stat_compare_means(label=
# "p.signif")). Reuses assign_group_subgroup() (same as plot_km_nested_subset())
# for the gene High/Low split -- including its "Low" reference-level factor
# convention -- and normalize_vector() for the cut-off gene, exactly like the
# app does at df$ScoreA/target_B normalization. No new shared helper needed.
#
# Draws immune cell fraction (y) by cell type (x), with two side-by-side boxes
# per cell type (High vs Low of the cut-off gene), and a Wilcoxon significance
# label ("ns" / "*" / "**" / "***" / "****") above each pair.
#
# Args:
#   df           - patient-level data frame (SampleID, gene_col, cell_types, ...)
#   gene_col     - gene used as the High/Low cut-off (e.g. "IRF2")
#   cell_types   - character vector of immune-cell-fraction column names
#                  (x-axis order); pass the app's own `cell_types` vector
#   cohort_col   - column identifying cohort (e.g. "StudyID")
#   cohort       - cohort value(s) to restrict to; NULL = use df as-is
#                  (i.e. already filtered/selected by the caller)
#   split_pct    - High/Low cut-off percentage (50 = median split)
#   subgroup_col - optional column for per-subgroup-wise splitting, passed
#                  straight through to assign_group_subgroup(); NULL = one
#                  global split across the (cohort-filtered) data
#   norm_type    - normalization applied to gene_col before splitting (skipped
#                  automatically if gene_col is itself a cell-fraction column,
#                  same rule the app applies to ScoreA/target_B)
#
# Returns a ggplot (ggpubr::ggboxplot) object, or NULL if there isn't enough
# data / group separation to plot.
# ------------------------------------------------------------------------------
plot_cell_fraction_boxplot <- function(df, gene_col, cell_types,
                                       cohort_col = "StudyID", cohort = NULL,
                                       split_pct = 50, subgroup_col = NULL,
                                       norm_type = c("none", "log2", "zscore", "robust")) {
  norm_type <- match.arg(norm_type)
  
  cell_types <- intersect(cell_types, colnames(df))
  if (length(cell_types) == 0 || !(gene_col %in% colnames(df))) return(NULL)
  
  if (!is.null(cohort) && cohort_col %in% colnames(df)) {
    df <- df %>% dplyr::filter(.data[[cohort_col]] %in% cohort)
  }
  df <- df %>% dplyr::filter(!is.na(.data[[gene_col]]))
  if (nrow(df) < 6) return(NULL)
  
  # normalize the cut-off gene unless it's itself a cell-fraction column
  # (already on a 0-1 scale) -- same rule the app applies at line ~1360
  if (!(gene_col %in% cell_types)) {
    df[[gene_col]] <- normalize_vector(df[[gene_col]], norm_type)
  }
  
  df_grp <- assign_group_subgroup(df, gene_col, split_pct, subgroup_col) %>%
    dplyr::rename(Plot_Group = Group)
  if (length(unique(df_grp$Plot_Group)) < 2) return(NULL)
  
  df_box <- df_grp %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cell_types), names_to = "CellType", values_to = "Fraction") %>%
    dplyr::filter(!is.na(Fraction)) %>%
    dplyr::mutate(CellType = factor(CellType, levels = cell_types))
  
  # qualitative palette for CellType fill, same construction method the app
  # uses for its `cell_colors` vector (RColorBrewer "qual" palettes concatenated)
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
  col_vector    <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  cell_pal      <- unname(col_vector[seq_along(cell_types)])
  grp_pal       <- c("Low" = "#377EB8", "High" = "#E41A1C")
  
  cohort_lbl <- if (!is.null(cohort)) paste(cohort, collapse = ", ") else "Selected Cohort"
  
  ggpubr::ggboxplot(df_box, x = "CellType", y = "Fraction", fill = "CellType", color = "Plot_Group",
                    bxp.errorbar = TRUE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(values = cell_pal) +
    ggplot2::scale_color_manual(values = grp_pal) +
    ggpubr::stat_compare_means(ggplot2::aes(group = Plot_Group), method = "wilcox.test", label = "p.signif") +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(subtitle = paste0(cohort_lbl, " \u2014 Split by ", gene_col,
                                    " (Top/Bottom ", split_pct, "%)"),
                  x = "Immune Cell Type", y = "Cell Fraction")
}

# ==============================================================================
#   7) compute_gene_cell_correlation() — Spearman rho / p-value per
#                                          (Cohort, CellType), same
#                                          Kendall-fallback convention as
#                                          the app's f1_heatmap_corr()
#   8) build_gene_cell_correlation_table() — pivots that into a DT table with
#                                          a Cohort column spanning two
#                                          sub-columns (rho, p-value) per
#                                          cohort, Cell Type as the row key
#
# DEPENDENCY NOTE: reuses normalize_vector() (for the gene column) -- no new
# shared helper is introduced beyond these two functions, since nothing in
# the earlier tasks already computes a correlation grid or a grouped-header
# table.
# ==============================================================================

# ------------------------------------------------------------------------------
# 7) compute_gene_cell_correlation()
#
# Template: the app's own f1_heatmap_corr() (lines ~903-937). Same logic,
# just keyed on Cohort x CellType instead of Subgroup x Target:
#   - Spearman's rho by default
#   - automatically falls back to Kendall's tau for a given (cohort, cell
#     type) pair if that cell type's fraction is >zero_frac_cutoff (default
#     30%) zero within that cohort -- identical rule to the app's Correlation
#     Plots tab, so results stay consistent with what's shown there.
#
# Args:
#   df                - patient-level data frame (SampleID, gene_col,
#                       cell_types, cohort_col, ...)
#   gene_col          - gene column to correlate against each cell type
#   cell_types        - character vector of cell-fraction column names
#   cohort_col        - column identifying cohort (e.g. "StudyID")
#   cohorts           - optional character vector to restrict/order cohorts;
#                       NULL = all cohorts present in df, alphabetical
#   norm_type         - normalization applied to gene_col before correlating
#                       (via normalize_vector()); skipped automatically if
#                       gene_col is itself a cell-fraction column
#   zero_frac_cutoff  - zero-value fraction threshold that triggers the
#                       Spearman -> Kendall fallback (default 0.30, matching
#                       the app)
#
# Returns a long-format data frame: Cohort, CellType, N, Method ("spearman" /
# "kendall" / "Insufficient" / "Error"), Rho, Pval -- one row per
# (Cohort, CellType) pair.
# ------------------------------------------------------------------------------
compute_gene_cell_correlation <- function(df, gene_col, cell_types, cohort_col = "StudyID",
                                          cohorts = NULL,
                                          norm_type = c("none", "log2", "zscore", "robust"),
                                          zero_frac_cutoff = 0.30) {
  norm_type <- match.arg(norm_type)
  
  cell_types <- intersect(cell_types, colnames(df))
  if (length(cell_types) == 0 || !(gene_col %in% colnames(df)) || !(cohort_col %in% colnames(df))) {
    return(NULL)
  }
  
  if (!is.null(cohorts)) df <- df %>% dplyr::filter(.data[[cohort_col]] %in% cohorts)
  cohort_lvls <- if (!is.null(cohorts)) cohorts else sort(unique(as.character(df[[cohort_col]])))
  if (length(cohort_lvls) == 0) return(NULL)
  
  if (!(gene_col %in% cell_types)) {
    df[[gene_col]] <- normalize_vector(df[[gene_col]], norm_type)
  }
  
  res_list <- list()
  for (co in cohort_lvls) {
    df_sub <- df %>% dplyr::filter(.data[[cohort_col]] == co)
    x_all  <- df_sub[[gene_col]]
    for (ct in cell_types) {
      y_all <- df_sub[[ct]]
      valid <- !is.na(x_all) & !is.na(y_all)
      x_val <- x_all[valid]; y_val <- y_all[valid]
      
      if (length(x_val) < 5) {
        res_list[[length(res_list) + 1]] <- data.frame(
          Cohort = co, CellType = ct, N = length(x_val), Method = "Insufficient", Rho = NA, Pval = NA)
        next
      }
      
      prop_zero  <- sum(y_val == 0) / length(y_val)
      cor_method <- ifelse(prop_zero > zero_frac_cutoff, "kendall", "spearman")
      test <- tryCatch(cor.test(x_val, y_val, method = cor_method, exact = FALSE), error = function(e) NULL)
      
      if (is.null(test)) {
        res_list[[length(res_list) + 1]] <- data.frame(
          Cohort = co, CellType = ct, N = length(x_val), Method = "Error", Rho = NA, Pval = NA)
      } else {
        res_list[[length(res_list) + 1]] <- data.frame(
          Cohort = co, CellType = ct, N = length(x_val), Method = cor_method,
          Rho = unname(test$estimate), Pval = test$p.value)
      }
    }
  }
  
  dplyr::bind_rows(res_list) %>%
    dplyr::mutate(CellType = factor(CellType, levels = cell_types),
                  Cohort   = factor(Cohort,   levels = cohort_lvls))
}

# ------------------------------------------------------------------------------
# Internal helper: pivots compute_gene_cell_correlation()'s long-format output
# into wide format (one row per CellType, two columns "<cohort>__rho" /
# "<cohort>__pval" per cohort). Shared by build_gene_cell_correlation_table()
# and the export functions below, so the pivot logic exists in exactly one
# place.
# ------------------------------------------------------------------------------
.pivot_corr_wide <- function(corr_df, cell_types, cohorts = NULL, digits_rho = 3, digits_pval = 3) {
  if (is.null(corr_df) || nrow(corr_df) == 0) return(NULL)
  
  cohort_lvls <- if (!is.null(cohorts)) cohorts else levels(corr_df$Cohort)
  
  wide_list <- lapply(cohort_lvls, function(co) {
    sub <- corr_df %>% dplyr::filter(as.character(Cohort) == co) %>%
      dplyr::select(CellType, Rho, Pval)
    names(sub)[2:3] <- paste0(co, c("__rho", "__pval"))
    sub
  })
  wide_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "CellType"), wide_list)
  
  cell_lvls <- intersect(cell_types, unique(as.character(wide_df$CellType)))
  wide_df <- wide_df %>%
    dplyr::mutate(CellType = factor(CellType, levels = cell_lvls)) %>%
    dplyr::arrange(CellType)
  
  rho_cols  <- grep("__rho$",  names(wide_df), value = TRUE)
  pval_cols <- grep("__pval$", names(wide_df), value = TRUE)
  wide_df[rho_cols]  <- lapply(wide_df[rho_cols],  function(x) round(x, digits_rho))
  wide_df[pval_cols] <- lapply(wide_df[pval_cols], function(x) signif(x, digits_pval))
  
  list(wide_df = wide_df, cohort_lvls = cohort_lvls)
}

# 8) build_gene_cell_correlation_table()
#
# Pivots compute_gene_cell_correlation()'s long-format output into a DT table
# with Cell Type as the row key and each Cohort spanning two sub-columns
# (rho, p-value), using a custom DT header "sketch" (the app already depends
# on DT for all its other tables, so no new package is introduced).
#
# Args:
#   corr_df      - output of compute_gene_cell_correlation()
#   cell_types   - character vector giving the row order (cell types)
#   cohorts      - optional character vector giving the cohort/column order;
#                  NULL = order found in corr_df$Cohort's factor levels
#   digits_rho   - decimal places for rho (default 3)
#   digits_pval  - significant figures for p-value (default 3)
#
# Returns a DT::datatable htmlwidget.
# ------------------------------------------------------------------------------
build_gene_cell_correlation_table <- function(corr_df, cell_types, cohorts = NULL,
                                              digits_rho = 3, digits_pval = 3) {
  pv <- .pivot_corr_wide(corr_df, cell_types, cohorts, digits_rho, digits_pval)
  if (is.null(pv)) return(NULL)
  wide_df     <- pv$wide_df %>% dplyr::rename(`Cell Type` = CellType)
  cohort_lvls <- pv$cohort_lvls
  
  # two-row header: row 1 = cohort names spanning 2 columns each, row 2 = rho / p-value
  sketch <- htmltools::withTags(table(
    class = "display",
    thead(
      tr(
        th(rowspan = 2, "Cell Type"),
        lapply(cohort_lvls, function(co) th(colspan = 2, co))
      ),
      tr(
        unlist(lapply(cohort_lvls, function(co) list(th("rho"), th("p-value"))), recursive = FALSE)
      )
    )
  ))
  
  DT::datatable(wide_df, container = sketch, rownames = FALSE,
                options = list(dom = "Bfrtip", scrollX = TRUE, scrollY = "500px",
                               scrollCollapse = TRUE, paging = FALSE))
}

# ==============================================================================
#   9) export_corr_to_csv()    — save compute_gene_cell_correlation() output
#                                 to .csv (long or wide format)
#  10) export_corr_to_excel()  — save it to .xlsx with a real two-row merged
#                                 header (Cohort spanning "rho"/"p-value"),
#                                 matching build_gene_cell_correlation_table()'s
#                                 layout
#  11) export_gene_cell_correlation() — convenience wrapper that picks CSV vs.
#                                 Excel from the file extension
#
# DEPENDENCY NOTE: reuses the .pivot_corr_wide() helper defined above (same
# pivot logic as build_gene_cell_correlation_table()) -- no duplicated
# pivoting code. Excel export requires the 'openxlsx' package (not currently
# a dependency of appv20260709.R); install with install.packages("openxlsx")
# if you don't have it. CSV export uses only base R.
# ==============================================================================

# ------------------------------------------------------------------------------
# 9) export_corr_to_csv()
#
# Args:
#   corr_df      - output of compute_gene_cell_correlation()
#   file_path    - destination path, e.g. "gene_cell_correlation.csv"
#   cell_types   - row order (only needed for format = "wide")
#   cohorts      - optional column order (only needed for format = "wide")
#   digits_rho   - decimal places for rho
#   digits_pval  - significant figures for p-value
#   format       - "wide" (Cell Type rows x Cohort rho/p-value column pairs,
#                  flattened to single-row headers like "CohortA (rho)",
#                  "CohortA (p-value)" since CSV has no merged cells) or
#                  "long" (one row per Cohort x CellType, as returned by
#                  compute_gene_cell_correlation() directly)
#
# Returns the exported data frame invisibly.
# ------------------------------------------------------------------------------
export_corr_to_csv <- function(corr_df, file_path, cell_types = NULL, cohorts = NULL,
                               digits_rho = 3, digits_pval = 3,
                               format = c("wide", "long")) {
  format <- match.arg(format)
  if (is.null(corr_df) || nrow(corr_df) == 0) {
    message("No correlation data to export.")
    return(invisible(NULL))
  }
  
  if (format == "long") {
    out_df <- corr_df %>%
      dplyr::mutate(Rho = round(Rho, digits_rho), Pval = signif(Pval, digits_pval)) %>%
      dplyr::select(Cohort, CellType, N, Method, Rho, Pval)
  } else {
    if (is.null(cell_types)) stop("cell_types is required when format = 'wide'.")
    pv <- .pivot_corr_wide(corr_df, cell_types, cohorts, digits_rho, digits_pval)
    if (is.null(pv)) { message("No correlation data to export."); return(invisible(NULL)) }
    out_df <- pv$wide_df %>% dplyr::rename(`Cell Type` = CellType)
    names(out_df) <- gsub("__rho$",  " (rho)",    names(out_df))
    names(out_df) <- gsub("__pval$", " (p-value)", names(out_df))
  }
  
  utils::write.csv(out_df, file_path, row.names = FALSE, na = "NA")
  message("Saved correlation table to: ", file_path)
  invisible(out_df)
}


# ------------------------------------------------------------------------------
# 10) export_corr_to_excel()
#
# Writes a real merged two-row header to .xlsx (row 1: cohort names spanning
# 2 columns each; row 2: "rho" / "p-value") using openxlsx, so the exported
# file visually matches build_gene_cell_correlation_table()'s on-screen layout.
#
# Args: same as export_corr_to_csv(), minus `format` (always wide here) plus:
#   sheet_name - worksheet name (default "Gene-Cell Correlation")
#
# Returns the wide data frame invisibly.
# ------------------------------------------------------------------------------
export_corr_to_excel <- function(corr_df, file_path, cell_types, cohorts = NULL,
                                 digits_rho = 3, digits_pval = 3,
                                 sheet_name = "Gene-Cell Correlation") {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for Excel export. Install it with install.packages('openxlsx').")
  }
  
  pv <- .pivot_corr_wide(corr_df, cell_types, cohorts, digits_rho, digits_pval)
  if (is.null(pv)) { message("No correlation data to export."); return(invisible(NULL)) }
  wide_df     <- pv$wide_df
  cohort_lvls <- pv$cohort_lvls
  n_cohorts   <- length(cohort_lvls)
  n_cols      <- 1 + 2 * n_cohorts   # Cell Type + 2 cols per cohort
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet_name)
  
  # Row 1-2: "Cell Type" merged vertically; each cohort name merged across its 2 cols
  openxlsx::writeData(wb, sheet_name, "Cell Type", startRow = 1, startCol = 1, colNames = FALSE)
  openxlsx::mergeCells(wb, sheet_name, cols = 1, rows = 1:2)
  
  col_idx <- 2
  for (co in cohort_lvls) {
    openxlsx::writeData(wb, sheet_name, co, startRow = 1, startCol = col_idx, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet_name, cols = col_idx:(col_idx + 1), rows = 1)
    openxlsx::writeData(wb, sheet_name, "rho",     startRow = 2, startCol = col_idx,     colNames = FALSE)
    openxlsx::writeData(wb, sheet_name, "p-value", startRow = 2, startCol = col_idx + 1, colNames = FALSE)
    col_idx <- col_idx + 2
  }
  
  header_style <- openxlsx::createStyle(
    fontName = "Arial", textDecoration = "bold", halign = "center", valign = "center",
    border = "TopBottomLeftRight", fgFill = "#D9E1F2")
  openxlsx::addStyle(wb, sheet_name, header_style, rows = 1:2, cols = 1:n_cols,
                     gridExpand = TRUE, stack = TRUE)
  
  body_style <- openxlsx::createStyle(fontName = "Arial", border = "TopBottomLeftRight")
  openxlsx::addStyle(wb, sheet_name, body_style, rows = 3:(2 + nrow(wide_df)), cols = 1:n_cols,
                     gridExpand = TRUE, stack = TRUE)
  
  # Data body starts at row 3, directly under the two header rows
  openxlsx::writeData(wb, sheet_name, wide_df, startRow = 3, startCol = 1,
                      colNames = FALSE, rowNames = FALSE)
  
  openxlsx::setColWidths(wb, sheet_name, cols = 1:n_cols, widths = c(24, rep(12, 2 * n_cohorts)))
  openxlsx::freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 2)
  
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  message("Saved correlation table to: ", file_path)
  invisible(wide_df)
}


# ------------------------------------------------------------------------------
# 11) export_gene_cell_correlation()
#
# Convenience wrapper: picks CSV vs. Excel export based on file_path's
# extension (.csv -> export_corr_to_csv(); .xlsx/.xls -> export_corr_to_excel()).
# ------------------------------------------------------------------------------
export_gene_cell_correlation <- function(corr_df, file_path, cell_types = NULL, cohorts = NULL,
                                         digits_rho = 3, digits_pval = 3,
                                         csv_format = c("wide", "long")) {
  csv_format <- match.arg(csv_format)
  ext <- tolower(tools::file_ext(file_path))
  
  if (ext %in% c("xlsx", "xls")) {
    export_corr_to_excel(corr_df, file_path, cell_types, cohorts, digits_rho, digits_pval)
  } else if (ext == "csv") {
    export_corr_to_csv(corr_df, file_path, cell_types, cohorts, digits_rho, digits_pval, format = csv_format)
  } else {
    stop("Unsupported file extension: '.", ext, "'. Use '.csv', '.xlsx', or '.xls'.")
  }
}

# ==============================================================================
#  12) plot_gene_gene_scatter() — Gene B (y) vs. Gene A (x) dot plot with a
#                                  linear regression line + 95% CI band,
#                                  faceted by cohort, Spearman's rho/p-value
#                                  annotated per cohort panel.
#
# Template: the app's own single-target scatter block inside the observe()
# loop (~lines 866-901) -- geom_point() + smooth line + correlation
# annotation, per Subgroup. This standalone version generalizes that to any
# two genes across any set of selected cohorts, using
# ggpubr::stat_cor(method = "spearman") (ggpubr is already a loaded app
# dependency) to compute and display rho/p per facet automatically, rather
# than hand-rolling per-cohort cor.test() + annotate() calls.
#
# DEPENDENCY NOTE: reuses normalize_vector() and the cell_types
# normalization-skip rule from plot_km_nested_subset() -- no other new
# helpers needed.
#
# Args:
#   df          - patient-level data frame (SampleID, geneA, geneB, cohort_col, ...)
#   geneA       - gene plotted on the x-axis
#   geneB       - gene plotted on the y-axis
#   cohort_col  - column identifying cohort (e.g. "StudyID")
#   cohorts     - character vector of cohorts to include (facet panels, in
#                 this order); NULL = all cohorts present in df, alphabetical
#   norm_type   - "none" | "log2" | "zscore" | "robust", applied to geneA/geneB
#                 before plotting (via normalize_vector()); skipped
#                 automatically for whichever of geneA/geneB is listed in
#                 `cell_types` (same rule as plot_km_nested_subset())
#   cell_types  - character vector of cell-fraction column names, used only
#                 to decide which of geneA/geneB to skip normalization for
#   free_scales - if TRUE (default), each cohort panel gets its own x/y
#                 range (facet_wrap(scales = "free")); FALSE locks all
#                 panels to the same axes
#
# Returns a ggplot object (geom_point + lm smooth + ggpubr::stat_cor,
# faceted by cohort), or NULL if there isn't enough data to plot.
# ------------------------------------------------------------------------------
plot_gene_gene_scatter <- function(df, geneA, geneB, cohort_col = "StudyID", cohorts = NULL,
                                   norm_type = c("none", "log2", "zscore", "robust"),
                                   cell_types = character(0), free_scales = TRUE) {
  norm_type <- match.arg(norm_type)
  
  if (is.null(df) || !all(c(geneA, geneB) %in% colnames(df)) || !(cohort_col %in% colnames(df))) {
    return(NULL)
  }
  
  if (!is.null(cohorts)) df <- df %>% dplyr::filter(.data[[cohort_col]] %in% cohorts)
  cohort_lvls <- if (!is.null(cohorts)) cohorts else sort(unique(as.character(df[[cohort_col]])))
  if (length(cohort_lvls) == 0) return(NULL)
  
  df <- df %>% dplyr::filter(!is.na(.data[[geneA]]), !is.na(.data[[geneB]]))
  if (nrow(df) < 6) return(NULL)
  
  # normalize geneA/geneB unless each is itself a cell-fraction column
  # (already on a 0-1 scale) -- same rule plot_km_nested_subset() applies
  if (!(geneA %in% cell_types)) df[[geneA]] <- normalize_vector(df[[geneA]], norm_type)
  if (!(geneB %in% cell_types)) df[[geneB]] <- normalize_vector(df[[geneB]], norm_type)
  
  df[[cohort_col]] <- factor(as.character(df[[cohort_col]]), levels = cohort_lvls)
  
  ax_suffix <- if (norm_type != "none") paste0(" (", norm_type, ")") else ""
  
  ggplot(df, ggplot2::aes(x = .data[[geneA]], y = .data[[geneB]])) +
    geom_point(alpha = 0.6, size = 2.2, color = "#2c3e50") +
    geom_smooth(method = "loess", se = TRUE, level = 0.95, color = "red", fill = "grey70", linewidth = 1) +
    ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", label.sep = "\n",
                     size = 4.2, color = "black") +
    facet_wrap(stats::as.formula(paste0("~`", cohort_col, "`")),
               scales = if (free_scales) "free" else "fixed") +
    labs(title = paste0(geneB, " vs. ", geneA),
         x = paste0(geneA, ax_suffix), y = paste0(geneB, ax_suffix)) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))
}

# ==============================================================================
# app_functions_for_notebook.R
#
# Standalone, notebook-ready versions of the analysis logic that lives inside
# appv20260823b.R's server() function (Features 1-3). Everything here is a
# plain function that takes a data frame + parameters and returns a
# plot/table/data.frame -- no `input$`, `reactive()`, `output$`, or Shiny
# dependency anywhere, so it can be sourced into any R notebook (including
# Figure9.Rmd) and called directly, then tweaked per-figure.
#
# Written to match the conventions already established in Figure9.Rmd:
#   - roxygen-style @param/@return docs
#   - tryCatch-wrapped model/plot calls, NULL return + message() on failure
#   - unnamed, strata-ordered palettes for ggsurvplot (avoids legend/color
#     misalignment)
#   - normalize_vector() applied per-column, never to cell-fraction (0-1)
#     columns
#   - year_axis_breaks() for consistent x-axis year labeling
#
# ------------------------------------------------------------------------------
# REDUNDANCY MAP -- what NOT to redefine
# ------------------------------------------------------------------------------
# Figure9.Rmd already defines the following, VERBATIM or functionally
# identical to appv20260823b.R. Source Figure9.Rmd (or copy just these
# chunks) and do not duplicate them again here:
#
#   normalize_vector()          - identical in both files
#   apply_admin_censoring()     - identical in both files
#   year_axis_breaks()          - identical in both files
#   assign_global_group()       - identical in both files
#   build_tailcut_groups()      - identical in both files
#   plot_km_tailcut()           - identical in both files
#   plot_gene_expr_barplot()    - identical in both files (only the app's
#                                 default xlab differs by one word; harmless)
#   plot_km_high_low()          - Figure9-only; this *is* the notebook-ready
#                                 version of the app's single-factor K-M plot
#                                 (app never split this out as its own
#                                 function -- it's inlined 3x in server(),
#                                 once each for F2 single-mode, F3
#                                 single-mode, and the F2 drill-down modal).
#                                 Use plot_km_high_low() for all three.
#   plot_km_by_group()          - Figure9-only; use for any pre-existing
#                                 categorical column (Response, Tx, StudyID...)
#
# Figure9.Rmd ALSO defines two functions that are *partial reimplementations*
# of app logic and should be UPGRADED/REPLACED by the versions in this file
# rather than kept as-is, because they silently drop functionality that
# exists in the app:
#
#   get_cox_zscore_pval()        -> replaced by get_km_zscore_pval() below.
#                                  Figure9's version only returns Z/Pval/N/
#                                  Events and only supports a single-factor
#                                  split. The app's version (used by every
#                                  cell of the Feature 2 heatmap) also
#                                  returns HR/CI_lo/CI_hi and supports the
#                                  two-factor nested/global anchor split.
#                                  compute_survival_heatmap() below is built
#                                  on the full version so the notebook
#                                  heatmap can carry an HR column too.
#
#   plot_km_nested_subset()      -> covers similar ground to
#                                  prep_nested_stratification() +
#                                  plot_subset_km() below, but only ever
#                                  shows 2 curves (Gene-B High/Low) + 1
#                                  reference curve, and only computes ONE
#                                  log-rank p-value (High vs Low). The app's
#                                  plot_subset_km() (reproduced below) is
#                                  what actually powers Feature 2's two-factor
#                                  drill-down and Feature 3's two-factor K-M
#                                  panels, and reports 3 p-values
#                                  (High-vs-Ref, Low-vs-Ref, High-vs-Low).
#                                  Keep plot_km_nested_subset() only if you
#                                  specifically want its simpler 3-curve look;
#                                  otherwise prefer prep_nested_stratification()
#                                  + plot_subset_km() for parity with the app.
#
# Everything below this line is NEW -- i.e. logic that currently only exists
# inlined inside appv20260823b.R's server() (reactives / renderPlot /
# renderDT blocks) and has no notebook-callable equivalent anywhere yet.
# ==============================================================================

library(tidyverse)
library(survival)
library(survminer)
library(ggpubr)
library(RColorBrewer)
library(DT)

# ==============================================================================
# SECTION 1 -- COHORT / SUBGROUP LABELING HELPERS
# (app Section 1 "GENERIC HELPERS"; not present in Figure9.Rmd at all)
# ==============================================================================

#' Map a raw StudyID to its human-readable cohort alias.
#' @param study_id character vector of StudyID values
#' @return character vector of aliases (unmatched IDs pass through unchanged)
rename_cohort <- function(study_id) {
  dplyr::case_when(
    study_id == "PRJEB23709" ~ "Cohort 2",
    study_id %in% c("PHS000452", "phs000452") ~ "Cohort 1.1",
    study_id == "GSE91061" ~ "Cohort 1.2",
    study_id == "PHS000452_GSE91061" ~ "Cohort 1",
    study_id %in% c("SU2CMARK", "SU2C-MARK") ~ "Cohort NSCLC N124",
    TRUE ~ as.character(study_id)
  )
}

#' Order a set of StudyID values into the app's canonical cohort display order.
#' @param studies character vector of StudyID values present in the data
#' @return `studies`, reordered (known IDs first in canonical order, any
#'   unrecognized IDs appended afterward in their original order)
order_studies <- function(studies) {
  target_study_order <- c("PRJEB23709", "PHS000452", "phs000452", "GSE91061",
                          "PHS000452_GSE91061", "SU2CMARK", "SU2C-MARK", "SU2CMARK_N124")
  matched   <- intersect(target_study_order, studies)
  unmatched <- setdiff(studies, target_study_order)
  c(matched, unmatched)
}

#' Collapse the distinct Biopsy Timing values within a group into one
#' readable string (e.g. c("On","Pre") -> "Pre & On"), so a paired-biopsy
#' group isn't misread as containing only one timepoint.
#' @param vals character/factor vector of BiopsyTime values within one group
#' @return character vector, "Pre..." labels first, then "On...", then others
order_biopsy_timing <- function(vals) {
  vals   <- unique(as.character(vals))
  vals   <- vals[!is.na(vals)]
  is_pre <- grepl("pre", vals, ignore.case = TRUE)
  is_on  <- grepl("on",  vals, ignore.case = TRUE) & !is_pre
  c(sort(vals[is_pre]), sort(vals[is_on]), sort(vals[!is_pre & !is_on]))
}

#' assign_global_group(), but optionally split independently within each
#' level of a subgroup column instead of once globally.
#' @param df            data.frame with OS_Time/OS_Status + val_col
#' @param val_col       continuous column to split on
#' @param split_pct     tail-cut percentage (5-50; 50 = median split)
#' @param subgroup_col  optional column name; when supplied (and not "none"),
#'   the High/Low cutoffs are computed separately within each level, and
#'   results are row-bound back together. NULL/"none" = one global split.
#' @return data.frame with a new `Group` factor column ("Low","High")
assign_group_subgroup <- function(df, val_col, split_pct, subgroup_col = NULL) {
  if (is.null(subgroup_col) || subgroup_col == "none") {
    return(assign_global_group(df, val_col, split_pct))
  }
  df_list <- df %>%
    dplyr::filter(!is.na(.data[[val_col]]) & !is.na(OS_Time) & !is.na(OS_Status)) %>%
    split(.[[subgroup_col]])
  purrr::map_dfr(df_list, ~assign_global_group(.x, val_col, split_pct))
}

#' Build a single "Hierarchical_Subgroup" factor column that mirrors Feature
#' 3's subgroup structure: "<Cohort Alias> (<StudyID>)" when multiple
#' cohorts are present, optionally suffixed with " - <subgroup value>" when
#' a subgroup_col is supplied. With one cohort and no subgroup_col, every row
#' gets the same single-level label instead of "All Cohorts".
#' @param df       data.frame with a StudyID column
#' @param sub_col  optional column name to further split each cohort by
#'   (e.g. "Prior_Tx", "BiopsyTime"); NULL = cohort-level grouping only
#' @return `df` with an added ordered factor column `Hierarchical_Subgroup`
add_hierarchical_subgroup <- function(df, sub_col = NULL) {
  selected_studies <- unique(as.character(df$StudyID))
  multi_cohort      <- length(selected_studies) > 1
  ordered_st        <- order_studies(selected_studies)

  if (is.null(sub_col) || sub_col == "none") {
    df$Hierarchical_Subgroup <- if (multi_cohort) "All Cohorts" else
      paste0(rename_cohort(selected_studies[1]), " (", selected_studies[1], ")")
    level_order <- unique(df$Hierarchical_Subgroup)
  } else if (sub_col == "StudyID") {
    df$Hierarchical_Subgroup <- paste0(rename_cohort(df$StudyID), " (", df$StudyID, ")")
    level_order <- purrr::map_chr(ordered_st, ~paste0(rename_cohort(.x), " (", .x, ")"))
    level_order <- intersect(level_order, unique(df$Hierarchical_Subgroup))
  } else if (multi_cohort) {
    df$Hierarchical_Subgroup <- paste0(rename_cohort(df$StudyID), " (", df$StudyID, ") - ", df[[sub_col]])
    level_order <- c()
    for (st in ordered_st) {
      st_label <- paste0(rename_cohort(st), " (", st, ")")
      sub_vals <- sort(unique(as.character(df[[sub_col]][df$StudyID == st])))
      if (length(sub_vals) > 0) level_order <- c(level_order, paste0(st_label, " - ", sub_vals))
    }
  } else {
    st_label <- paste0(rename_cohort(selected_studies[1]), " (", selected_studies[1], ")")
    df$Hierarchical_Subgroup <- paste0(st_label, " - ", df[[sub_col]])
    sub_vals <- sort(unique(as.character(df[[sub_col]])))
    level_order <- paste0(st_label, " - ", sub_vals)
  }
  df$Hierarchical_Subgroup <- factor(df$Hierarchical_Subgroup, levels = level_order)
  df
}

#' Expand a data frame into a stacked "Subgroup" factor column: one copy of
#' the whole cohort per StudyID (in canonical order), optionally with an
#' "All Cohorts" copy prepended and/or per-subgroup-value copies nested
#' under each StudyID. This is the shared shape both Feature 1 (f1_data) and
#' Feature 2 (f2_data_proc) build by hand inline in the app -- consolidated
#' here into one function so notebook figures don't have to re-derive it.
#' @param df          patient-level data.frame with a StudyID column
#' @param subgroup_col optional column name to nest under each StudyID
#'   ("none"/NULL = cohort-level rows only)
#' @param show_all    prepend one row-set labeled "All Cohorts" (or
#'   "All Cohorts (N = ...)" if `label_n = TRUE`) covering every row
#' @param label_n     append " (N = <n>)" to every Subgroup label (Feature 1
#'   style); FALSE reproduces Feature 2's plain labels
#' @return `df` row-bound across all requested subsets, with an added
#'   ordered factor column `Subgroup`
build_cohort_subgroups <- function(df, subgroup_col = NULL, show_all = FALSE, label_n = TRUE) {
  sub_col <- if (is.null(subgroup_col) || subgroup_col == "none") NULL else subgroup_col
  mk_lbl  <- function(base, n) if (label_n) paste0(base, " (N = ", n, ")") else base

  df_list <- list()
  if (show_all) {
    df_all <- df
    df_all$Subgroup <- mk_lbl(if (label_n) "All Cohorts" else "All Cohorts", nrow(df_all))
    df_list[[length(df_list) + 1]] <- df_all
  }
  for (st in order_studies(unique(as.character(df$StudyID)))) {
    df_st <- df %>% dplyr::filter(StudyID == st)
    if (nrow(df_st) == 0) next
    mapped_st <- rename_cohort(st)
    df_st$Subgroup <- mk_lbl(mapped_st, nrow(df_st))
    df_list[[length(df_list) + 1]] <- df_st
    if (!is.null(sub_col)) {
      for (sv in sort(unique(as.character(df_st[[sub_col]])))) {
        if (is.na(sv)) next
        df_sv <- df_st %>% dplyr::filter(.data[[sub_col]] == sv)
        if (nrow(df_sv) == 0) next
        df_sv$Subgroup <- mk_lbl(paste0(mapped_st, " | ", sv), nrow(df_sv))
        df_list[[length(df_list) + 1]] <- df_sv
      }
    }
  }
  if (length(df_list) == 0) return(NULL)
  dplyr::bind_rows(df_list) %>% dplyr::mutate(Subgroup = factor(Subgroup, levels = unique(Subgroup)))
}

# ==============================================================================
# SECTION 2 -- FEATURE 1: GENE-CELL / GENE-GENE CORRELATION
# (app Section: server "FEATURE 1"; Figure9.Rmd only has the *scatter* half
#  of this, via plot_gene_gene_scatter() -- the heatmap/matrix half below is
#  entirely new)
# ==============================================================================

#' Compute Spearman/Kendall correlation of one anchor column against several
#' target columns, across every level of a Subgroup factor. Kendall's tau is
#' used automatically in place of Spearman's rho whenever a target is a cell
#' fraction with >30% zero values (matches the app's Correlation Plots tab).
#' @param df            output of build_cohort_subgroups() (needs a
#'   `Subgroup` factor column) with an `Expr_Val` (pre-normalized anchor)
#'   column already present, OR pass `anchor_col` + `norm_type` and this
#'   function will normalize it for you
#' @param targets       character vector of target column names (genes or
#'   cell types)
#' @param anchor_col    name of the raw anchor column (used only if
#'   `Expr_Val` is not already present in `df`)
#' @param norm_type     normalization applied to anchor_col (and to any
#'   target not in `cell_types`) if not already normalized
#' @param cell_types    character vector of cell-fraction column names
#'   (used for the zero-fraction Kendall fallback + normalization skip)
#' @return long data.frame: Subgroup, Target, Rho, Pval, N, Method, plus
#'   Signif/Label/Subgroup_Label/Target ready for plot_corr_heatmap()
compute_correlation_by_subgroup <- function(df, targets, anchor_col = NULL,
                                            norm_type = c("none", "log2", "zscore", "robust"),
                                            cell_types = character(0)) {
  norm_type <- match.arg(norm_type)
  targets   <- intersect(targets, colnames(df))
  if (length(targets) == 0) return(NULL)

  if (!("Expr_Val" %in% colnames(df))) {
    if (is.null(anchor_col)) stop("Supply anchor_col, or pre-compute an Expr_Val column.")
    df$Expr_Val <- normalize_vector(df[[anchor_col]], norm_type)
  }
  for (tgt in targets) if (!(tgt %in% cell_types)) df[[tgt]] <- normalize_vector(df[[tgt]], norm_type)

  grps <- levels(droplevels(factor(df$Subgroup)))
  res_list <- list()
  for (g in grps) {
    df_sub <- df %>% dplyr::filter(as.character(Subgroup) == g)
    for (tgt in targets) {
      x_val <- df_sub$Expr_Val; y_val <- df_sub[[tgt]]
      valid <- !is.na(x_val) & !is.na(y_val); x_val <- x_val[valid]; y_val <- y_val[valid]
      if (length(x_val) < 5) {
        res_list[[length(res_list) + 1]] <- data.frame(Subgroup = g, Target = tgt, Rho = NA, Pval = NA, N = NA, Method = "Insufficient")
        next
      }
      prop_zero  <- if (tgt %in% cell_types) sum(y_val == 0) / length(y_val) else 0
      cor_method <- ifelse(prop_zero > 0.30, "kendall", "spearman")
      test <- tryCatch(cor.test(x_val, y_val, method = cor_method, exact = FALSE), error = function(e) NULL)
      res_list[[length(res_list) + 1]] <- if (is.null(test))
        data.frame(Subgroup = g, Target = tgt, Rho = NA, Pval = NA, N = length(x_val), Method = "Error")
      else
        data.frame(Subgroup = g, Target = tgt, Rho = unname(test$estimate), Pval = test$p.value,
                   N = length(x_val), Method = cor_method)
    }
  }
  dplyr::bind_rows(res_list) %>%
    dplyr::mutate(Signif = dplyr::case_when(is.na(Pval) ~ "", Pval < 0.001 ~ "***",
                                            Pval < 0.01 ~ "**", Pval < 0.05 ~ "*", TRUE ~ ""),
                 Label          = ifelse(is.na(Rho), "NA", paste0(round(Rho, 2), "\n", Signif)),
                 Target         = factor(Target,  levels = targets),
                 Subgroup_Label = factor(Subgroup, levels = grps))
}

#' Draw the Target x Subgroup correlation heatmap from
#' compute_correlation_by_subgroup()'s output. Plain ggplot (no plotly) --
#' wrap in plotly::ggplotly() yourself if you want interactive drill-down.
#' @param corr_df output of compute_correlation_by_subgroup()
#' @param title   subtitle text
#' @return a ggplot object, or NULL
plot_corr_heatmap <- function(corr_df, title = "Spearman \u03C1 or Kendall \u03C4 (*p<0.05, **p<0.01, ***p<0.001)") {
  if (is.null(corr_df) || nrow(corr_df) == 0) return(NULL)
  ggplot(corr_df, aes(x = Target, y = Subgroup_Label, fill = Rho)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Label), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey80", name = "Correlation") +
    scale_y_discrete(limits = rev(levels(corr_df$Subgroup_Label))) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
          axis.text.y = element_text(face = "bold", size = 12), panel.grid = element_blank()) +
    labs(subtitle = title, x = NULL, y = NULL)
}

#' Single anchor-vs-target scatter plot for one Subgroup, with a loess fit
#' and a rho/p-value annotation (auto Spearman/Kendall, same 30%-zero rule
#' as compute_correlation_by_subgroup()). Reproduces the app's single-target
#' Correlation Plots panel (and its modal drill-down scatter) as one
#' reusable function.
#' @param df_sub    single-Subgroup slice with `Expr_Val` + `target` columns
#' @param target    target column name (already normalized if desired)
#' @param anchor_label,target_label,norm_label axis/title labels
#' @param cell_types cell-fraction column names, for the Kendall fallback rule
#' @return a ggplot object
plot_corr_scatter <- function(df_sub, target, anchor_label = "Expr_Val", target_label = target,
                              norm_label = "", cell_types = character(0), subtitle = NULL) {
  prop_zero  <- if (target %in% cell_types) sum(df_sub[[target]] == 0, na.rm = TRUE) / nrow(df_sub) else 0
  cor_method <- ifelse(prop_zero > 0.30, "kendall", "spearman")
  stat_label <- ifelse(prop_zero > 0.30, "Kendall \u03C4", "Spearman \u03C1")
  note_text  <- if (prop_zero > 0.30) sprintf("Note: Kendall's Tau used (%.1f%% zeroes)", prop_zero * 100) else "(Spearman's Rank)"
  cor_res    <- tryCatch(cor.test(df_sub$Expr_Val, df_sub[[target]], method = cor_method, exact = FALSE),
                         error = function(e) list(estimate = NA, p.value = NA))
  annot_text <- sprintf("%s = %.3f\nP-value = %.2e\n%s", stat_label, cor_res$estimate, cor_res$p.value, note_text)

  ggplot(df_sub, aes(x = Expr_Val, y = .data[[target]])) +
    geom_point(alpha = 0.6, size = 2.5, color = "#2c3e50") +
    geom_smooth(method = "loess", color = "#E41A1C", fill = "grey70", linewidth = 1) +
    labs(title = paste(anchor_label, "vs", target_label), subtitle = subtitle,
        x = paste0(anchor_label, " (", norm_label, ")"), y = paste0(target_label, " (", norm_label, ")")) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
             label = annot_text, size = 5, fontface = "bold", color = "black") +
    theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold"))
}

#' Reshape a LONG correlation table (Subgroup x Target, with Rho/Pval) into
#' WIDE format: one row per Target, two columns (rho, p-value) per Subgroup.
#' Generic on column names -- reuse for gene-cell OR gene-gene results.
#' @return list(wide_df, group_levels), or NULL
pivot_corr_wide <- function(corr_df, group_col, item_col, rho_col = "Rho", pval_col = "Pval",
                            item_levels, group_levels, digits_rho = 3, digits_pval = 3) {
  if (is.null(corr_df) || nrow(corr_df) == 0) return(NULL)
  wide_list <- lapply(group_levels, function(g) {
    sub <- corr_df %>% dplyr::filter(as.character(.data[[group_col]]) == g) %>%
      dplyr::select(dplyr::all_of(c(item_col, rho_col, pval_col)))
    names(sub) <- c(item_col, paste0(g, "__rho"), paste0(g, "__pval"))
    sub
  })
  wide_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = item_col), wide_list)
  wide_df[[item_col]] <- factor(as.character(wide_df[[item_col]]), levels = item_levels)
  wide_df <- wide_df %>% dplyr::arrange(.data[[item_col]])
  rho_cols  <- grep("__rho$",  names(wide_df), value = TRUE)
  pval_cols <- grep("__pval$", names(wide_df), value = TRUE)
  wide_df[rho_cols]  <- lapply(wide_df[rho_cols],  function(x) round(x, digits_rho))
  wide_df[pval_cols] <- lapply(wide_df[pval_cols], function(x) signif(x, digits_pval))
  list(wide_df = wide_df, group_levels = group_levels)
}

#' Write a pivot_corr_wide() result to .xlsx with a real merged two-row
#' header (row 1: group names spanning 2 cols; row 2: "rho"/"p-value").
#' Requires 'openxlsx' (not loaded by default here).
write_wide_corr_excel <- function(pv, file, row_label = "Item") {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Install 'openxlsx' for this export.")
  wide_df <- pv$wide_df; group_lvls <- pv$group_levels
  n_cols  <- 1 + 2 * length(group_lvls)
  wb <- openxlsx::createWorkbook(); sheet_name <- "Correlation Matrix"
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, row_label, startRow = 1, startCol = 1, colNames = FALSE)
  openxlsx::mergeCells(wb, sheet_name, cols = 1, rows = 1:2)
  col_idx <- 2
  for (g in group_lvls) {
    openxlsx::writeData(wb, sheet_name, g, startRow = 1, startCol = col_idx, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet_name, cols = col_idx:(col_idx + 1), rows = 1)
    openxlsx::writeData(wb, sheet_name, "rho",     startRow = 2, startCol = col_idx,     colNames = FALSE)
    openxlsx::writeData(wb, sheet_name, "p-value", startRow = 2, startCol = col_idx + 1, colNames = FALSE)
    col_idx <- col_idx + 2
  }
  header_style <- openxlsx::createStyle(fontName = "Arial", textDecoration = "bold", halign = "center",
                                        valign = "center", border = "TopBottomLeftRight", fgFill = "#D9E1F2")
  openxlsx::addStyle(wb, sheet_name, header_style, rows = 1:2, cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
  openxlsx::writeData(wb, sheet_name, wide_df, startRow = 3, startCol = 1, colNames = FALSE, rowNames = FALSE)
  openxlsx::setColWidths(wb, sheet_name, cols = 1:n_cols, widths = c(24, rep(12, 2 * length(group_lvls))))
  openxlsx::freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 2)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

# ==============================================================================
# SECTION 3 -- FEATURE 2 / SURVIVAL HEATMAP CORE STATS
# get_km_zscore_pval() and compute_median_os_diff() are the two pieces of
# real analytic logic in the app that Figure9.Rmd does not have a full
# equivalent for. get_km_zscore_pval() supersedes Figure9's
# get_cox_zscore_pval() (see redundancy note at top of file).
# ==============================================================================

#' Cox PH Wald Z-score + HR/CI + log-rank p-value for one target column,
#' optionally restricted to one arm of an "anchor" two-factor split. This is
#' the exact per-cell computation behind every cell of the Feature 2
#' Survival Heatmap.
#' @param df_sub     data for one cohort/subgroup slice
#' @param anchor     Factor-A column name (two_factor mode only)
#' @param target     column to split High/Low and test
#' @param split_pct  Factor-A (or single-mode) tail-cut percentage
#' @param mode       "single" or "two_factor"
#' @param view       which Factor-A arm to restrict to in two_factor mode
#'   ("High"/"Low")
#' @param tf_method  "nested" (split B within the A arm) or "global" (split
#'   B across everyone, then intersect with the A arm)
#' @param split_pct_B Factor-B tail-cut percentage (two_factor mode)
#' @return list(Z, Pval, N, Events, HR, CI_lo, CI_hi)
get_km_zscore_pval <- function(df_sub, anchor, target, split_pct,
                               mode = c("single", "two_factor"), view = "High",
                               tf_method = c("nested", "global"), split_pct_B = NULL) {
  mode <- match.arg(mode); tf_method <- match.arg(tf_method)
  na_out <- list(Z = NA, Pval = NA, N = NA, Events = NA, HR = NA, CI_lo = NA, CI_hi = NA)

  df_eval <- df_sub %>% dplyr::filter(!is.na(.data[[target]]))
  if (mode == "two_factor") df_eval <- df_eval %>% dplyr::filter(!is.na(.data[[anchor]]))
  if (nrow(df_eval) < 10) return(na_out)

  if (mode == "single") {
    df_eval <- assign_global_group(df_eval, target, split_pct)
  } else if (tf_method == "nested") {
    df_A     <- assign_global_group(df_eval, anchor, split_pct) %>% dplyr::rename(GroupA = Group)
    df_A_str <- df_A %>% dplyr::filter(GroupA == view)
    df_eval  <- assign_global_group(df_A_str, target, split_pct_B)
  } else {
    df_A <- assign_global_group(df_eval, anchor, split_pct) %>% dplyr::rename(GroupA = Group)
    df_B <- assign_global_group(df_eval, target, split_pct_B)
    common  <- intersect(df_A$SampleID, df_B$SampleID)
    df_eval <- df_eval %>% dplyr::filter(SampleID %in% common) %>%
      dplyr::left_join(df_A %>% dplyr::select(SampleID, GroupA), by = "SampleID") %>%
      dplyr::left_join(df_B %>% dplyr::select(SampleID, Group),  by = "SampleID") %>%
      dplyr::filter(GroupA == view)
  }

  if (nrow(df_eval) < 5 || length(unique(df_eval$Group)) < 2 ||
      length(unique(df_eval$OS_Status)) < 2) return(na_out)

  sdiff   <- tryCatch(survdiff(Surv(OS_Time, OS_Status) ~ Group, data = df_eval), error = function(e) NULL)
  pval    <- if (!is.null(sdiff)) pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE) else NA
  cox_fit <- tryCatch(coxph(Surv(OS_Time, OS_Status) ~ Group, data = df_eval), error = function(e) NULL)
  if (!is.null(cox_fit)) {
    scox <- summary(cox_fit)
    list(Z = scox$coefficients[1, "z"], Pval = pval, N = nrow(df_eval),
         Events = sum(df_eval$OS_Status == 1, na.rm = TRUE),
         HR = scox$conf.int[1, "exp(coef)"], CI_lo = scox$conf.int[1, "lower .95"],
         CI_hi = scox$conf.int[1, "upper .95"])
  } else {
    list(Z = NA, Pval = pval, N = nrow(df_eval), Events = sum(df_eval$OS_Status == 1, na.rm = TRUE),
         HR = NA, CI_lo = NA, CI_hi = NA)
  }
}

#' Assemble the Target x Subgroup survival heatmap grid (single- or
#' two-factor), calling get_km_zscore_pval() once per cell. Supersedes
#' Figure9's compute_survival_heatmap_data(): adds two_factor support and
#' the HR/CI_lo/CI_hi columns.
#' @param df          data.frame with a Subgroup factor column, e.g. from
#'   build_cohort_subgroups()
#' @param targets     character vector of target columns (x-axis)
#' @param mode        "single" or "two_factor"
#' @param anchor      Factor-A column name (two_factor mode only)
#' @param split_pct   Factor-A (or single-mode) tail-cut percentage
#' @param split_pct_B Factor-B tail-cut percentage (two_factor mode)
#' @param view        Factor-A arm to restrict to (two_factor mode)
#' @param tf_method   "nested" or "global" (two_factor mode)
#' @return long data.frame ready for plot_survival_heatmap(), with
#'   Z_capped/Signif/Label/Target/Subgroup_Label plus HR/CI_lo/CI_hi
compute_survival_heatmap <- function(df, targets, mode = c("single", "two_factor"),
                                     anchor = NULL, split_pct = 50, split_pct_B = NULL,
                                     view = "High", tf_method = c("nested", "global")) {
  mode <- match.arg(mode); tf_method <- match.arg(tf_method)
  targets <- intersect(targets, colnames(df))
  if (length(targets) == 0 || !("Subgroup" %in% colnames(df))) return(NULL)
  grps <- levels(droplevels(factor(df$Subgroup)))

  anchor_strata_lbl <- if (mode == "single") "All" else paste0("A-", view)
  comparison_lbl    <- if (mode == "single") "High vs Low" else "B-High vs B-Low"
  reference_lbl     <- if (mode == "single") "Low" else "B-Low"

  res_list <- list()
  for (g in grps) {
    df_sub <- df %>% dplyr::filter(as.character(Subgroup) == g)
    for (tgt in targets) {
      stats <- get_km_zscore_pval(df_sub, anchor, tgt, split_pct, mode, view, tf_method, split_pct_B)
      res_list[[length(res_list) + 1]] <- data.frame(
        Subgroup = g, Target = tgt, Anchor_Strata = anchor_strata_lbl,
        Comparison = comparison_lbl, Reference = reference_lbl,
        N = stats$N, Events = stats$Events, HR = stats$HR, CI_lo = stats$CI_lo, CI_hi = stats$CI_hi,
        Z = stats$Z, Pval = stats$Pval)
    }
  }
  dplyr::bind_rows(res_list) %>%
    dplyr::mutate(Z_capped = ifelse(is.na(Z), NA, pmin(pmax(Z, -5), 5)),
                 Signif   = dplyr::case_when(is.na(Pval) ~ "", Pval < 0.001 ~ "***",
                                             Pval < 0.01 ~ "**", Pval < 0.05 ~ "*", TRUE ~ ""),
                 Label    = ifelse(is.na(Z_capped), "NA", paste0(round(Z_capped, 2), "\n", Signif)),
                 Target         = factor(Target,  levels = targets),
                 Subgroup_Label = factor(Subgroup, levels = grps))
}

#' Render a Target x Subgroup survival heatmap from compute_survival_heatmap()
#' output. Same visual convention as plot_corr_heatmap()/the app's f2_heatmap_ui.
#' @param hdf         output of compute_survival_heatmap()
#' @param title       subtitle text
#' @param orientation "landscape" (default; Target on x-axis, Subgroup on
#'   y-axis -- matches the app) or "portrait" (axes swapped: Subgroup on
#'   x-axis, Target on y-axis -- useful when you have many targets and few
#'   subgroups, or need a tall/narrow figure)
plot_survival_heatmap <- function(hdf, title = "Survival Heatmap (High vs Low)",
                                  orientation = c("landscape", "portrait")) {
  if (is.null(hdf) || nrow(hdf) == 0) return(NULL)
  orientation <- match.arg(orientation)
  x_var <- if (orientation == "landscape") "Target" else "Subgroup_Label"
  y_var <- if (orientation == "landscape") "Subgroup_Label" else "Target"

  ggplot(hdf, aes(x = .data[[x_var]], y = .data[[y_var]], fill = Z_capped)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = Label), size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                         midpoint = 0, limits = c(-5, 5), oob = scales::squish,
                         na.value = "grey80", name = "Z-score") +
    scale_y_discrete(limits = rev(levels(hdf[[y_var]]))) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
          axis.text.y = element_text(face = "bold", size = 12), panel.grid = element_blank()) +
    labs(subtitle = title, x = NULL, y = NULL)
}

#' High-vs-Low difference in median OS, with not-reached (NR) handling.
#' If either arm's median is NR (surv_median() returns NA), substitutes that
#' ARM'S OWN last observed follow-up time (its true ">=" bound) and flags the
#' result as a bounded estimate rather than an exact one, so callers never
#' silently treat a median-of-NA as a real number, and never silently treat
#' a bounded estimate as if it were exact.
#'
#' BUG FIX (this revision): the previous version substituted the pooled
#' admin-censoring end-of-study date (`max_years`, the x-axis ceiling shared
#' by BOTH arms) for an NR arm's median, instead of that specific arm's own
#' last observed time. Because apply_admin_censoring() only ever truncates
#' follow-up, never extends it, a small/short-follow-up arm (e.g. a tight
#' 15% or 35% expression tail-cut) can have its longest-followed patient end
#' well before the pooled cohort's max -- especially when the single
#' longest-followed patient in the whole cohort happens to fall in the OTHER
#' arm. Using the pooled end-of-study date for the NR arm silently inflated
#' its imputed "median" (and hence Diff_yrs) whenever that happened. Each
#' arm's own last observed follow-up time is now used instead, giving the
#' tightest defensible ">=" bound for that arm specifically. THIS FUNCTION IS
#' NOT PRESENT ANYWHERE IN Figure9.Rmd -- it is new in the 2026-08-23 app update.
#'
#' TWO WAYS TO CALL THIS:
#'   1) Pre-grouped: `df` already has a categorical grouping column (e.g.
#'      you built it yourself with assign_global_group(), or it's the
#'      Plot_Group/GroupB column from compute_factor_profiling_data() /
#'      prep_nested_stratification()). Leave `value_col = NULL` (default)
#'      and set `group_col` to that column's name -- this is the original
#'      calling convention and still works exactly as before.
#'   2) Raw, same input as plot_km_high_low(): pass the SAME `df`,
#'      `value_col`, `split_pct`, `norm_type`, and `max_years` you passed
#'      (or would pass) to plot_km_high_low(), and leave `group_col` at its
#'      default. This function then reproduces plot_km_high_low()'s exact
#'      preprocessing (apply_admin_censoring() -> normalize_vector() ->
#'      assign_global_group()) before computing the median difference, so
#'      the numbers are guaranteed to match the plot -- there is no need to
#'      call plot_km_high_low() first or capture its return value.
#'
#' @param df         data.frame. Either already grouped (see mode 1 above),
#'   or raw patient-level data with `value_col` present (mode 2).
#' @param value_col  NULL (default, mode 1) or the name of the continuous
#'   column to split on (mode 2) -- same argument as plot_km_high_low().
#'   Supplying this switches the function into mode 2 and takes precedence
#'   over `group_col`/pre-existing grouping.
#' @param split_pct  tail-cut percentage (mode 2 only; default 50 = median
#'   split, same default as plot_km_high_low())
#' @param norm_type  "none" | "log2" | "zscore" | "robust" (mode 2 only;
#'   applied to value_col before splitting, same as plot_km_high_low())
#' @param group_col  categorical column with (at least) group_high/group_low
#'   (mode 1 only; ignored in mode 2, where the derived split is always
#'   named "Group", matching assign_global_group()'s output)
#' @param time_col,status_col OS time (days) / event status columns
#' @param max_years  admin-censoring cutoff, in years, used only to truncate
#'   follow-up before splitting (mode 2) / as a required sanity bound (mode
#'   1). No longer used as the NR substitute value itself -- see bug-fix note
#'   above. Mode 2 defaults this the same way plot_km_high_low() does
#'   (ceiling of the observed max follow-up) when left NULL/missing; mode 1
#'   requires it explicitly, as before.
#' @param group_high,group_low the two group labels to compare (defaults
#'   match both assign_global_group()'s and plot_km_high_low()'s "High"/"Low")
#' @return one-row data.frame: Median_High_days/yrs, Median_Low_days/yrs,
#'   Diff_days/yrs, Is_Bounded_Estimate, Bound_Note -- or NULL
compute_median_os_diff <- function(df, group_col = "Group", time_col = "OS_Time",
                                   status_col = "OS_Status", max_years = NULL,
                                   group_high = "High", group_low = "Low",
                                   value_col = NULL, split_pct = 50,
                                   norm_type = c("none", "log2", "zscore", "robust")) {
  norm_type <- match.arg(norm_type)
  
  if (!is.null(value_col)) {
    # ---- Mode 2: raw df + value_col, mirrors plot_km_high_low()'s pipeline ----
    # apply_admin_censoring() (shared with plot_km_high_low()) hardcodes the
    # column names OS_Time/OS_Status, same as everywhere else in the app --
    # so mode 2 only supports that pair, not arbitrary time_col/status_col.
    if (!identical(time_col, "OS_Time") || !identical(status_col, "OS_Status")) {
      stop("Mode 2 (value_col supplied) requires time_col = 'OS_Time' and status_col = 'OS_Status', ",
           "since apply_admin_censoring() -- shared with plot_km_high_low() -- assumes those names.")
    }
    req_cols <- c(time_col, status_col, value_col)
    missing_cols <- setdiff(req_cols, colnames(df))
    if (length(missing_cols) > 0) stop("Input table is missing required column(s): ", paste(missing_cols, collapse = ", "))
    
    df <- df[!is.na(df[[value_col]]) & !is.na(df[[time_col]]) & !is.na(df[[status_col]]), ]
    if (nrow(df) < 5) return(NULL)
    
    if (is.null(max_years)) max_years <- ceiling(max(df[[time_col]], na.rm = TRUE) / 365.25)
    df <- apply_admin_censoring(df, max_years)  # assumes OS_Time/OS_Status; fine since time_col/status_col default to those
    df$Expr_Val <- normalize_vector(df[[value_col]], norm_type)
    
    df <- tryCatch(assign_global_group(df, "Expr_Val", split_pct), error = function(e) NULL)
    if (is.null(df) || nrow(df) < 2 || length(unique(df$Group)) < 2) return(NULL)
    group_col <- "Group"
  }
  
  if (is.null(max_years)) stop("max_years is required in mode 1 (pre-grouped df).")
  df <- df[!is.na(df[[time_col]]) & !is.na(df[[status_col]]) & !is.na(df[[group_col]]), ]
  if (nrow(df) < 2 || length(unique(df[[group_col]])) < 2) return(NULL)
  fmla <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", group_col))
  fit  <- tryCatch(survminer::surv_fit(fmla, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  med_df <- survminer::surv_median(fit)
  med_df$Group <- gsub(paste0("^", group_col, "="), "", med_df$strata)

  # For an NR arm, the tightest defensible ">=" bound is THAT ARM'S OWN last
  # observed follow-up time (whether the last row was an event or a
  # censoring) -- NOT the pooled/global admin-censoring cutoff (max_years),
  # which can badly overstate a small or short-follow-up arm's true
  # follow-up whenever the single longest-followed patient in the cohort
  # happens to sit in the OTHER arm. See bug-fix note above.
  get_med <- function(g) {
    val <- med_df$median[med_df$Group == g]
    group_days     <- df[[time_col]][df[[group_col]] == g]
    group_max_days <- suppressWarnings(max(group_days, na.rm = TRUE))
    if (!is.finite(group_max_days)) group_max_days <- NA_real_  # shouldn't happen post-filtering; guards div-by-empty
    if (length(val) == 0 || is.na(val)) list(days = group_max_days, nr = TRUE) else list(days = val, nr = FALSE)
  }
  high <- get_med(group_high); low <- get_med(group_low)
  diff_days <- high$days - low$days

  bound_note <- dplyr::case_when(
    high$nr && low$nr ~ "Both arms NR \u2014 difference not meaningful",
    high$nr            ~ paste0("Lower-bound only: High median NR (last follow-up in High arm used), true difference is \u2265 ", round(diff_days / 365.25, 2), " yrs"),
    low$nr              ~ paste0("Upper-bound only: Low median NR (last follow-up in Low arm used), true difference is \u2264 ", round(diff_days / 365.25, 2), " yrs"),
    TRUE                ~ "Exact")
  
  data.frame(Median_High_days = high$days, Median_Low_days = low$days,
             Median_High_yrs  = round(high$days / 365.25, 2), Median_Low_yrs = round(low$days / 365.25, 2),
             Diff_days = diff_days, Diff_yrs = round(diff_days / 365.25, 2),
             Is_Bounded_Estimate = high$nr || low$nr, Bound_Note = bound_note, stringsAsFactors = FALSE)
}

# ==============================================================================
# SECTION 4 -- TWO-FACTOR STRATIFICATION + DRILL-DOWN K-M
# prep_nested_stratification/prep_global_stratification/plot_subset_km are
# what actually power every two-factor K-M plot in the app (Feature 2 modal
# AND Feature 3). NOT present in Figure9.Rmd. Prefer these over Figure9's
# plot_km_nested_subset() for exact parity with the app (see redundancy note).
# ==============================================================================

#' Nested two-factor split: split the WHOLE cohort by Factor A, then split
#' Factor B independently WITHIN each Factor-A arm.
#' @param df           data.frame (already admin-censored/NA-filtered)
#' @param genes_A      character vector of Factor-A columns (summed if >1)
#' @param target_B     Factor-B column name
#' @param splitA,splitB tail-cut percentages
#' @param subgroup_col optional per-subgroup split column
#' @return data.frame with GroupA, GroupB, Stratum ("A-High / B-High", ...)
prep_nested_stratification <- function(df, genes_A, target_B, splitA, splitB, subgroup_col = NULL) {
  df_proc <- df %>% dplyr::mutate(FactorA_SumScore = rowSums(dplyr::select(., dplyr::all_of(genes_A)), na.rm = TRUE))
  df_proc <- assign_group_subgroup(df_proc, "FactorA_SumScore", splitA, subgroup_col) %>% dplyr::rename(GroupA = Group)
  df_list <- df_proc %>% split(.$GroupA)
  if (!is.null(df_list$High) && nrow(df_list$High) > 0)
    df_list$High <- assign_group_subgroup(df_list$High, target_B, splitB, subgroup_col) %>% dplyr::rename(GroupB = Group)
  if (!is.null(df_list$Low) && nrow(df_list$Low) > 0)
    df_list$Low  <- assign_group_subgroup(df_list$Low,  target_B, splitB, subgroup_col) %>% dplyr::rename(GroupB = Group)
  df_stratified <- dplyr::bind_rows(df_list) %>%
    dplyr::mutate(GroupB  = factor(GroupB, levels = c("Low", "High")),
                 Stratum = factor(paste0("A-", GroupA, " / B-", GroupB),
                                  levels = c("A-High / B-High", "A-High / B-Low", "A-Low / B-High", "A-Low / B-Low")))
  cols_to_keep <- unique(c("SampleID", "OS_Time", "OS_Status", "Stratum", "GroupA", "GroupB",
                           "FactorA_SumScore", genes_A, target_B, subgroup_col))
  df_stratified %>% dplyr::select(dplyr::any_of(cols_to_keep))
}

#' Global two-factor split: split Factor A and Factor B each independently
#' across the WHOLE cohort, then intersect (patient is High-A AND High-B,
#' etc.) -- answers a different question than nested (co-occurrence of two
#' independent biomarker extremes, vs. "does B matter within A's extreme?").
prep_global_stratification <- function(df, genes_A, target_B, splitA, splitB, subgroup_col = NULL) {
  df_proc <- df %>% dplyr::mutate(FactorA_SumScore = rowSums(dplyr::select(., dplyr::all_of(genes_A)), na.rm = TRUE))
  df_A <- assign_group_subgroup(df_proc, "FactorA_SumScore", splitA, subgroup_col) %>%
    dplyr::rename(GroupA = Group) %>% dplyr::select(SampleID, GroupA)
  df_B <- assign_group_subgroup(df_proc, target_B, splitB, subgroup_col) %>%
    dplyr::rename(GroupB = Group) %>% dplyr::select(SampleID, GroupB)
  df_proc %>%
    dplyr::inner_join(df_A, by = "SampleID") %>% dplyr::inner_join(df_B, by = "SampleID") %>%
    dplyr::filter(!is.na(GroupA), !is.na(GroupB)) %>%
    dplyr::mutate(GroupA  = factor(GroupA, levels = c("Low", "High")),
                 GroupB  = factor(GroupB, levels = c("Low", "High")),
                 Stratum = factor(paste0("A-", GroupA, " / B-", GroupB),
                                  levels = c("A-High / B-High", "A-High / B-Low", "A-Low / B-High", "A-Low / B-Low")))
}

#' Draw the two-factor drill-down K-M plot for ONE Factor-A arm: Factor-B
#' High vs Low vs an "Overall (Reference)" curve (everyone in that A arm,
#' regardless of B). Reports 3 log-rank p-values (B-High vs Ref, B-Low vs
#' Ref, B-High vs B-Low). This is the exact function behind the Feature 2
#' modal drill-down and Feature 3's two-factor K-M panels.
#' @param df_strat  output of prep_nested_stratification() or
#'   prep_global_stratification()
#' @param subset_A  which Factor-A arm to plot ("High" or "Low")
#' @param labelA,labelB display names for Factor A / Factor B
#' @param max_years admin-censoring / x-axis cutoff, in years
#' @param show_metrics "both" | "median" | "landmark" | "none"
#' @param landmark_years landmark year(s) for the landmark OS rate caption
#' @return a survminer ggsurvplot object, or NULL
plot_subset_km <- function(df_strat, subset_A, labelA, labelB, max_years,
                           show_metrics = "both", landmark_years = 1) {
  df_sub <- df_strat %>% dplyr::filter(GroupA == subset_A)
  if (nrow(df_sub) < 2 || length(unique(df_sub$GroupB)) < 2) return(NULL)
  max_time_days    <- max_years * 365.25
  year_breaks_days <- seq(0, max_time_days, by = 365.25)
  year_labels      <- round(year_breaks_days / 365.25)
  df_ref  <- df_sub %>% dplyr::mutate(GroupB = "Overall (Reference)")
  df_plot <- dplyr::bind_rows(df_sub, df_ref)
  valid_b <- intersect(c("High", "Low", "Overall (Reference)"), unique(df_plot$GroupB))
  df_plot$GroupB <- factor(df_plot$GroupB, levels = valid_b)

  get_pval <- function(df_pair) {
    if (length(unique(df_pair$GroupB)) < 2) return(NA)
    sdiff <- tryCatch(survdiff(Surv(OS_Time, OS_Status) ~ GroupB, data = df_pair), error = function(e) NULL)
    if (is.null(sdiff)) return(NA)
    pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)
  }
  p_high <- get_pval(df_plot %>% dplyr::filter(GroupB %in% c("High", "Overall (Reference)")))
  p_low  <- get_pval(df_plot %>% dplyr::filter(GroupB %in% c("Low",  "Overall (Reference)")))
  p_hl   <- get_pval(df_plot %>% dplyr::filter(GroupB %in% c("High", "Low")))
  sub_text <- paste0("Log-Rank p-values:\nB-High vs B-Low: ", ifelse(is.na(p_hl), "NA", format.pval(p_hl, digits = 3)),
                     "\nB-High vs Ref: ", ifelse(is.na(p_high), "NA", format.pval(p_high, digits = 3)),
                     " | B-Low vs Ref: ", ifelse(is.na(p_low), "NA", format.pval(p_low, digits = 3)))

  fit <- tryCatch(survfit(Surv(OS_Time, OS_Status) ~ GroupB, data = df_plot), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  caption_lines <- c()
  if (show_metrics %in% c("median", "both")) {
    med_summary <- survminer::surv_median(fit)
    med_text <- sapply(seq_len(nrow(med_summary)), function(i) {
      val <- med_summary$median[i]; sg <- gsub("GroupB=", "", med_summary$strata[i])
      sg  <- if (sg == "Overall (Reference)") "Ref" else paste0("B-", sg)
      paste0(sg, ": ", ifelse(is.na(val), "NE", paste(round(val / 365.25, 2), "yrs")))
    })
    caption_lines <- c(caption_lines, paste("Median OS \u2192", paste(med_text, collapse = " | ")))
  }
  if (show_metrics %in% c("landmark", "both")) {
    land_summ <- summary(fit, times = landmark_years * 365.25)
    if (length(land_summ$strata) > 0) {
      land_text <- sapply(seq_along(land_summ$strata), function(i) {
        sg <- gsub("GroupB=", "", land_summ$strata[i]); sg <- if (sg == "Overall (Reference)") "Ref" else paste0("B-", sg)
        paste0(sg, ": ", round(land_summ$surv[i] * 100, 1), "%")
      })
      caption_lines <- c(caption_lines, paste0(landmark_years, "-Year OS Rate \u2192 ", paste(land_text, collapse = " | ")))
    } else caption_lines <- c(caption_lines, paste0(landmark_years, "-Year OS Rate \u2192 NE"))
  }
  caption_text <- if (length(caption_lines) > 0) paste(caption_lines, collapse = "\n") else NULL

  cols <- c("High" = "#E41A1C", "Low" = "#377EB8", "Overall (Reference)" = "black")
  lty  <- c("High" = "solid",   "Low" = "solid",   "Overall (Reference)" = "dashed")
  lbls <- sapply(valid_b, function(x) if (x == "Overall (Reference)") paste0("Reference (A-", subset_A, ")") else paste0("B-", x))

  p <- tryCatch(
    ggsurvplot(fit, data = df_plot, pval = FALSE, risk.table = TRUE, risk.table.height = 0.35,
               palette = unname(cols[valid_b]), linetype = unname(lty[valid_b]),
               legend.labs = unname(lbls), legend = "bottom",
               title = paste0("A-", subset_A, " Cohort\n(", labelB, " Impact)"), subtitle = sub_text,
               xscale = 'd_y', break.time.by = 365.25, xlim = c(0, max_time_days)),
    error = function(e) NULL)
  if (is.null(p)) return(NULL)
  p$plot <- p$plot + scale_x_continuous(breaks = year_breaks_days, labels = year_labels) +
    labs(x = "Time (Years)", caption = caption_text) +
    theme(plot.subtitle = element_text(size = 10, face = "italic", color = "grey30"),
          plot.caption  = element_text(face = "bold", size = 11, hjust = 1, color = "grey30", lineheight = 1.2))
  p$table <- p$table + scale_x_continuous(breaks = year_breaks_days, labels = year_labels) + labs(x = "Time (Years)")
  p
}

# ==============================================================================
# SECTION 5 -- FEATURE 3: COMPREHENSIVE PROFILING ORCHESTRATION
# compute_factor_profiling_data() extracts f3_data_proc(); the Cox stats
# table, expr boxplot, cell-fraction boxplot, and cell-fraction stacked bar
# each extract one Feature-3 render block. NONE of this exists in
# Figure9.Rmd. Note: for the resulting single-factor K-M plot itself, reuse
# Figure9's plot_km_by_group(df_out, group_col = "Plot_Group") rather than
# writing a new one -- it is functionally identical to the app's Feature 3
# single-mode K-M render.
# ==============================================================================

#' Build Feature 3's stratified analysis table: a Factor-A sum-score (single
#' gene, or summed if genes_A has >1 element), grouped into a
#' Hierarchical_Subgroup, then split single- or two-factor.
#' @param df           patient-level data.frame (already admin-censored)
#' @param genes_A      character vector of Factor-A columns (summed if >1;
#'   never normalized if a column is itself a cell-fraction)
#' @param mode         "single" or "two_factor"
#' @param target_B     Factor-B column name (two_factor mode only)
#' @param split_pct    Factor-A tail-cut percentage
#' @param split_pct_B  Factor-B tail-cut percentage (two_factor mode)
#' @param tf_method    "nested" or "global" (two_factor mode)
#' @param subgroup_col optional column for the Hierarchical_Subgroup nesting
#' @param norm_type    normalization applied to genes_A/target_B (skipped
#'   for any column in `cell_types`)
#' @param cell_types   character vector of cell-fraction column names
#' @return list(df, mode, labelA, labelB, grps) -- `df` carries Plot_Group
#'   (single mode) or GroupA/GroupB/Plot_Group (two_factor mode)
compute_factor_profiling_data <- function(df, genes_A, mode = c("single", "two_factor"),
                                          target_B = NULL, split_pct = 50, split_pct_B = NULL,
                                          tf_method = c("nested", "global"), subgroup_col = NULL,
                                          norm_type = c("none", "log2", "zscore", "robust"),
                                          cell_types = character(0)) {
  mode <- match.arg(mode); tf_method <- match.arg(tf_method); norm_type <- match.arg(norm_type)

  df <- add_hierarchical_subgroup(df, subgroup_col)
  df$ScoreA <- 0
  for (g in genes_A) if (g %in% colnames(df)) {
    val <- if (g %in% cell_types) df[[g]] else normalize_vector(df[[g]], norm_type)
    df$ScoreA <- df$ScoreA + val
  }
  labelA <- if (length(genes_A) > 1) "Factor A (Sum Score)" else genes_A[1]

  if (mode == "single") {
    df_out <- assign_group_subgroup(df, "ScoreA", split_pct, "Hierarchical_Subgroup") %>% dplyr::rename(Plot_Group = Group)
    return(list(df = df_out, mode = "single", labelA = labelA, grps = levels(df$Hierarchical_Subgroup)))
  }

  stopifnot(!is.null(target_B))
  if (!(target_B %in% cell_types)) df[[target_B]] <- normalize_vector(df[[target_B]], norm_type)

  if (tf_method == "nested") {
    df_A <- assign_group_subgroup(df, "ScoreA", split_pct, "Hierarchical_Subgroup") %>% dplyr::rename(GroupA = Group)
    dl <- df_A %>% split(.$GroupA)
    if (!is.null(dl$High) && nrow(dl$High) > 0)
      dl$High <- assign_group_subgroup(dl$High, target_B, split_pct_B, "Hierarchical_Subgroup") %>% dplyr::rename(GroupB = Group)
    if (!is.null(dl$Low) && nrow(dl$Low) > 0)
      dl$Low  <- assign_group_subgroup(dl$Low,  target_B, split_pct_B, "Hierarchical_Subgroup") %>% dplyr::rename(GroupB = Group)
    df_out <- dplyr::bind_rows(dl)
  } else {
    dA <- assign_group_subgroup(df, "ScoreA", split_pct, "Hierarchical_Subgroup") %>%
      dplyr::rename(GroupA = Group) %>% dplyr::select(SampleID, GroupA)
    dB <- assign_group_subgroup(df, target_B, split_pct_B, "Hierarchical_Subgroup") %>%
      dplyr::rename(GroupB = Group) %>% dplyr::select(SampleID, GroupB)
    df_out <- df %>% dplyr::inner_join(dA, by = "SampleID") %>% dplyr::inner_join(dB, by = "SampleID") %>%
      dplyr::filter(!is.na(GroupA), !is.na(GroupB))
  }
  df_out <- df_out %>%
    dplyr::mutate(GroupB     = factor(GroupB, levels = c("Low", "High")),
                 Plot_Group = factor(paste0("A-", GroupA, " / B-", GroupB),
                                     levels = c("A-High / B-High", "A-High / B-Low", "A-Low / B-High", "A-Low / B-Low"))) %>%
    dplyr::filter(!is.na(Plot_Group))
  list(df = df_out, mode = "two_factor", labelA = labelA, labelB = target_B, grps = levels(df$Hierarchical_Subgroup))
}

#' Cox PH stats table for compute_factor_profiling_data() output: one row
#' per Hierarchical_Subgroup (single mode), or up to 3 rows per (Subgroup,
#' Factor-A arm) in two_factor mode (B-High-vs-B-Low, plus each B arm vs. an
#' "A-Median (Ref)" pseudo-baseline built from the whole A arm).
#' @param res  output of compute_factor_profiling_data()
#' @return data.frame: Subgroup, Anchor_Strata, Comparison, Reference, N,
#'   Events, HR, CI_lo, CI_hi, P_Val
compute_cox_stats_table <- function(res) {
  df <- res$df
  grps <- res$grps[res$grps %in% df$Hierarchical_Subgroup]

  run_cox_row <- function(df_model, group_col, ref_level, subgroup_lbl, anchor_strata_lbl, comparison_lbl, ref_lbl) {
    if (nrow(df_model) < 5 || length(unique(df_model$OS_Status)) < 2 ||
        length(unique(df_model[[group_col]])) < 2) return(NULL)
    df_model[[group_col]] <- relevel(factor(df_model[[group_col]]), ref = ref_level)
    cox_fit <- tryCatch(coxph(as.formula(paste0("Surv(OS_Time, OS_Status) ~ ", group_col)), data = df_model), error = function(e) NULL)
    if (is.null(cox_fit)) return(NULL)
    sm <- summary(cox_fit)
    data.frame(Subgroup = subgroup_lbl, Anchor_Strata = anchor_strata_lbl, Comparison = comparison_lbl,
              Reference = ref_lbl, N = sm$n, Events = sm$nevent, HR = round(sm$conf.int[1, 1], 3),
              CI_lo = round(sm$conf.int[1, 3], 3), CI_hi = round(sm$conf.int[1, 4], 3),
              P_Val = format.pval(sm$coefficients[1, 5], digits = 4))
  }

  purrr::map_dfr(grps, function(g) {
    df_sub <- df %>% dplyr::filter(Hierarchical_Subgroup == g)
    if (res$mode == "single") {
      run_cox_row(df_sub, "Plot_Group", "Low", g, "All", "High vs Low", "Low")
    } else {
      purrr::map_dfr(c("High", "Low"), function(a_strata) {
        df_A <- df_sub %>% dplyr::filter(GroupA == a_strata); if (nrow(df_A) < 5) return(NULL)
        rows <- list(run_cox_row(df_A, "GroupB", "Low", g, paste0("A-", a_strata), "B-High vs B-Low", "B-Low"))
        for (b_level in c("High", "Low")) {
          df_test <- df_A %>% dplyr::filter(GroupB == b_level); if (nrow(df_test) < 2) next
          df_paired <- dplyr::bind_rows(df_A %>% dplyr::mutate(CoxGroup = "A-Median (Ref)"),
                                        df_test %>% dplyr::mutate(CoxGroup = paste0("B-", b_level)))
          rows[[length(rows) + 1]] <- run_cox_row(df_paired, "CoxGroup", "A-Median (Ref)", g, paste0("A-", a_strata),
                                                  paste0("B-", b_level, " vs A-", a_strata, " Median"),
                                                  paste0("A-", a_strata, " Median"))
        }
        dplyr::bind_rows(rows)
      })
    }
  })
}

#' Boxplot of the Factor-A score (ScoreA) by Plot_Group, for ONE
#' Hierarchical_Subgroup slice. Single mode adds a Wilcoxon High-vs-Low
#' bracket automatically.
#' @param df_sub  single-subgroup slice of compute_factor_profiling_data()$df
#' @param mode    "single" or "two_factor"
#' @param labelA  Factor-A display name (used in the title/y-axis)
#' @param norm_label normalization label shown on the y-axis
#' @param subtitle plot subtitle (e.g. the subgroup name)
plot_factor_score_boxplot <- function(df_sub, mode = c("single", "two_factor"), labelA, norm_label = "", subtitle = NULL) {
  mode <- match.arg(mode)
  fill_pal <- if (mode == "single") c(Low = "#377EB8", High = "#E41A1C") else
    c("A-High / B-High" = "#a50f15", "A-High / B-Low" = "#fb6a4a", "A-Low / B-High" = "#6baed6", "A-Low / B-Low" = "#08519c")
  p <- ggplot(df_sub, aes(x = Plot_Group, y = ScoreA, fill = Plot_Group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "grey30") +
    theme_bw(base_size = 14) + scale_fill_manual(values = fill_pal) +
    labs(title = paste(labelA, "Expression Distribution"), subtitle = subtitle,
        x = "Stratification Group", y = paste0(labelA, " (", norm_label, ")")) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  if (mode == "single") p <- p + ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = list(c("High", "Low")))
  p
}

#' Immune cell fraction boxplot (x = CellType) colored by an arbitrary
#' Plot_Group factor (2 levels for single-factor mode, 4 for two-factor).
#' Generalizes Figure9's plot_cell_fraction_boxplot(), which only supports
#' a 2-level High/Low grouping -- pass a 4-level `color_pal` here to also
#' cover Feature 3's two-factor Plot_Group.
#' @param df_sub      single-subgroup slice with cell_types + a grouping col
#' @param cell_types  character vector of cell-fraction column names
#' @param group_col   grouping column name (e.g. "Plot_Group")
#' @param color_pal   named color vector, one entry per level of group_col
#' @param cell_colors named/unnamed color vector for the CellType fill
#' @param subtitle    plot subtitle
plot_cell_fraction_by_group <- function(df_sub, cell_types, group_col = "Plot_Group",
                                        color_pal, cell_colors, subtitle = NULL) {
  df_box <- df_sub %>% tidyr::pivot_longer(cols = dplyr::all_of(cell_types), names_to = "CellType", values_to = "fraction")
  ggpubr::ggboxplot(df_box, x = "CellType", y = "fraction", fill = "CellType", color = group_col, bxp.errorbar = TRUE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::scale_fill_manual(values = cell_colors) + ggplot2::scale_color_manual(values = color_pal) +
    ggpubr::stat_compare_means(ggplot2::aes(group = .data[[group_col]]), label = "p.signif") +
    ggplot2::guides(fill = "none") + ggplot2::labs(subtitle = subtitle)
}

#' Stacked bar of per-sample cell-type composition, faceted by an arbitrary
#' Plot_Group factor. NOT present anywhere in Figure9.Rmd.
#' @param df_sub      single-subgroup slice with cell_types + group_col +
#'   a sample-id column
#' @param cell_types  character vector of cell-fraction column names
#' @param group_col   column to facet by (e.g. "Plot_Group")
#' @param sample_col  per-sample x-axis column (e.g. "SampleID")
#' @param cell_colors named/unnamed color vector for the CellType fill
#' @param subtitle    plot subtitle
plot_cell_fraction_stackedbar <- function(df_sub, cell_types, group_col = "Plot_Group",
                                          sample_col = "SampleID", cell_colors, subtitle = NULL) {
  df_bar <- df_sub %>% tidyr::pivot_longer(cols = dplyr::all_of(cell_types), names_to = "CellType", values_to = "fraction")
  ggplot(df_bar, aes(x = .data[[sample_col]], y = fraction, fill = CellType)) +
    geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values = cell_colors) +
    facet_grid(stats::reformulate(group_col, "."), scales = "free", space = "free") +
    theme_minimal() + theme(axis.text.x = element_blank()) + labs(subtitle = subtitle)
}

# ==============================================================================
# SECTION 6 -- COHORT DISTRIBUTION SUMMARY TABLE
# Extracted from output$cohort_summary_table. NOT present in Figure9.Rmd.
# ==============================================================================

#' Rebuild the app's "Active Cohort Distribution Summary" table: one row per
#' (StudyID, Prior_Tx, BiopsyMatch, Regimen) group, with Biopsy Timing
#' values collapsed via order_biopsy_timing() so a matched-pair group shows
#' "Pre & On" instead of being split across two rows.
#' @param df filtered patient-level data.frame with StudyID/Prior_Tx/
#'   BiopsyTime/BiopsyMatch/Regimen columns
#' @return data.frame ready for DT::datatable() or plain printing
build_cohort_summary_table <- function(df) {
  summary_df <- df %>%
    dplyr::mutate(Cohort_Alias = rename_cohort(StudyID)) %>%
    dplyr::group_by(StudyID, Cohort_Alias, Prior_Tx, BiopsyMatch, Regimen) %>%
    dplyr::summarise(Biopsy_Timing = paste(order_biopsy_timing(BiopsyTime), collapse = " & "),
                     Patients = dplyr::n(), .groups = "drop")
  ordered_studies <- order_studies(unique(as.character(summary_df$StudyID)))
  summary_df %>%
    dplyr::mutate(StudyID = factor(StudyID, levels = ordered_studies)) %>%
    dplyr::arrange(StudyID, Prior_Tx) %>%
    dplyr::select(`Cohort Alias` = Cohort_Alias, `Original StudyID` = StudyID, `Prior Treatment` = Prior_Tx,
                 `Biopsy Timing` = Biopsy_Timing, `Biopsy Match` = BiopsyMatch, Regimen = Regimen, N = Patients)
}

# ==============================================================================
# SECTION 6B -- ONE BIOPSY PER PATIENT (survival-analysis deduplication)
# Extracted from df_filtered()'s dedup branch + output$dedup_patient_warning.
# NOT present anywhere in Figure9.Rmd. This is what the app's "One Biopsy per
# Patient (for Survival Analysis)" toggle actually does: for any PatientID
# with more than one row (e.g. matched Pre/Post biopsies, which share the
# same OS_Time/OS_Status), keep exactly one row -- the `prefer`-timepoint row
# if that patient has one, otherwise whatever other row they DO have. This is
# NOT the same as simply filtering to `BiopsyTime == prefer`: a plain filter
# also DROPS every patient who has no row at the preferred timepoint (e.g. an
# On-only, unmatched patient), silently shrinking N, whereas the app's actual
# behavior only removes the *extra* rows from patients who have more than one.
# ==============================================================================

#' Collapse a patient-level table to exactly one row per patient, preferring
#' a given biopsy timepoint when a patient has more than one row. Patients
#' with only one row (at any timepoint) are always kept, unlike a plain
#' `filter(BiopsyTime == prefer)`, which would drop them.
#' @param df           patient-level data.frame with `patient_col` and
#'   `biopsy_col` columns (e.g. PatientID, BiopsyTime)
#' @param prefer       biopsy_col value to keep when a patient has >1 row
#'   (default "Pre", matching the app's default `dedup_bx_preference`)
#' @param patient_col  column identifying the patient (default "PatientID")
#' @param biopsy_col   column identifying the biopsy timepoint (default
#'   "BiopsyTime")
#' @return `df` with exactly one row per distinct `patient_col` value
dedup_one_biopsy_per_patient <- function(df, prefer = "Pre",
                                         patient_col = "PatientID", biopsy_col = "BiopsyTime") {
  req_cols <- c(patient_col, biopsy_col)
  missing_cols <- setdiff(req_cols, colnames(df))
  if (length(missing_cols) > 0) stop("Input table is missing required column(s): ", paste(missing_cols, collapse = ", "))

  df %>%
    dplyr::mutate(.dedup_rank = ifelse(.data[[biopsy_col]] == prefer, 0, 1)) %>%
    dplyr::arrange(.data[[patient_col]], .dedup_rank) %>%
    dplyr::distinct(.data[[patient_col]], .keep_all = TRUE) %>%
    dplyr::select(-.dedup_rank)
}

#' Count (and optionally list) patients with more than one row in `df` --
#' the same check behind the app's "Warning: N patient(s) currently have
#' more than one biopsy row..." banner. Run this BEFORE and AFTER
#' dedup_one_biopsy_per_patient() (or before/after a manual subset() filter)
#' to confirm deduplication worked and that no unintended patients were
#' dropped rather than just de-duplicated.
#' @param df          data.frame with `patient_col`
#' @param patient_col column identifying the patient (default "PatientID")
#' @return list(n_duplicated = <int>, duplicated_ids = <character vector>)
count_duplicate_patients <- function(df, patient_col = "PatientID") {
  dup_tbl <- df %>% dplyr::count(.data[[patient_col]], name = "n_rows") %>% dplyr::filter(n_rows > 1)
  list(n_duplicated = nrow(dup_tbl), duplicated_ids = as.character(dup_tbl[[patient_col]]))
}

# --- THE SCALED BOXPLOT FUNCTION ---
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)

plot_comparable_boxplot <- function(expr_matrix, # Symbols in rownames, samples in columns
                                    metadata_df,
                                    target_genes,
                                    title,
                                    ref_group = "Normal",
                                    transform = c("log2", "zscore", "none")) {
  
  # 0. Initial Setup
  transform <- match.arg(transform)
  colnames(expr_matrix) <- gsub("^X", "", colnames(expr_matrix))
  
  # 1. Sample Matching Logic
  id_cols <- c("geo_accession", "title", "barcode")
  common_samples <- character(0)
  
  for (col in id_cols) {
    if (col %in% colnames(metadata_df)) {
      test_samples <- intersect(colnames(expr_matrix), metadata_df[[col]])
      if (length(test_samples) > 0) {
        common_samples <- test_samples
        # Important: set rownames based on the match to ensure alignment
        temp_meta <- metadata_df
        rownames(temp_meta) <- temp_meta[[col]]
        metadata_df <- temp_meta
        message("Successfully matched using column: ", col)
        break
      }
    }
  }
  
  if (length(common_samples) == 0) {
    stop("Error: No overlapping samples found.")
  }
  
  # 2. Subset and Relevel
  expr_sub <- expr_matrix[, common_samples, drop = FALSE]
  meta_sub <- metadata_df[common_samples, , drop = FALSE]
  meta_sub$group <- as.factor(meta_sub$group)
  if (ref_group %in% levels(meta_sub$group)) {
    meta_sub$group <- relevel(meta_sub$group, ref = ref_group)
  }
  
  # 3. Data Transformation
  present_genes <- intersect(target_genes, rownames(expr_sub))
  subset_expr <- expr_sub[present_genes, , drop = FALSE]
  
  if (transform == "log2") {
    # Calculate log2(TPM + 1)
    processed_expr <- log2(subset_expr + 1)
    y_label <- "Expression [log2(TPM+1)]"
  } else if (transform == "zscore") {
    # Log transform first, then Z-score (standard practice)
    log_expr <- log2(subset_expr + 1)
    processed_expr <- t(apply(log_expr, 1, scale))
    colnames(processed_expr) <- colnames(subset_expr)
    y_label <- "Expression (Z-score)"
  } else {
    processed_expr <- subset_expr
    y_label <- "Expression (Raw)"
  }
  
  # 4. Pivot Data for ggplot
  data_long <- as.data.frame(t(processed_expr)) %>%
    mutate(Group = meta_sub$group) %>%
    pivot_longer(cols = all_of(present_genes), names_to = "Gene", values_to = "Val") %>%
    mutate(Gene = factor(Gene, levels = present_genes))  # keep facet order stable
  
  # 5. Define Comparison Pairs (ref_group vs each other level)
  group_names <- levels(meta_sub$group)
  my_comparisons <- list()
  if (length(group_names) >= 2) {
    for (i in 2:length(group_names)) {
      my_comparisons[[i - 1]] <- c(group_names[1], group_names[i])
    }
  }
  
  # 6. Pre-compute stats + bracket positions PER FACET (per Gene)
  # This is the key fix: stat_compare_means()'s automatic label.y placement
  # is NOT reliably facet-aware when facet_wrap(scales = "free_y") is used,
  # so brackets/labels can be calculated using the wrong y-range and end up
  # drawn outside the visible panel (i.e. "missing"). We compute an explicit
  # y.position for every comparison ourselves, separately within each Gene's
  # own data range, then stat_pvalue_manual() draws exactly at those
  # coordinates. add_significance() guarantees every pair gets a label
  # ("ns" through "****"), never blank.
  # (Note: deliberately NOT using rstatix::add_xy_position() here since it
  # requires rstatix >= 0.6.0 and isn't available in all environments; the
  # manual y.position calculation below achieves the same result.)
  stat.test <- data_long %>%
    group_by(Gene) %>%
    wilcox_test(Val ~ Group, comparisons = my_comparisons) %>%
    ungroup() %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  gene_ymax <- data_long %>%
    group_by(Gene) %>%
    summarise(y_base = max(Val, na.rm = TRUE), .groups = "drop")
  
  stat.test <- stat.test %>%
    left_join(gene_ymax, by = "Gene") %>%
    group_by(Gene) %>%
    mutate(y.position = y_base + y_base * 0.08 + (row_number() - 1) * (y_base * 0.12)) %>%
    ungroup() %>%
    dplyr::select(-y_base)
  
  # 7. Generate Plot
  p <- ggboxplot(data_long,
                 x = "Group",
                 y = "Val",
                 color = "Group",
                 palette = "npj",
                 add = "jitter",                    # show individual sample points
                 add.params = list(size = 1.4, alpha = 0.5, jitter = 0.2),
                 outlier.shape = NA,                 # avoid double-plotting outliers; jitter already shows all points
                 short.panel.labs = TRUE,
                 panel.labs.font = list(face = "bold")) +
    stat_pvalue_manual(stat.test,
                       label = "p.adj.signif",      # always draws a label: ns, *, **, ***, ****
                       hide.ns = FALSE,
                       tip.length = 0.01,
                       bracket.size = 0.4) +
    labs(title = title,
         y = y_label,
         x = NULL) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 10),
                       expand = expansion(mult = c(0, 0.08))) +  # small extra headroom above computed brackets
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          strip.background = element_rect(fill = "gray95", color = "gray80"),
          strip.text = element_text(size = 11),
          legend.position = "none") +
    facet_wrap(~Gene)
  
  return(p)
}


# --- CV (COEFFICIENT OF VARIATION) CALCULATION & COMPARISON HELPER ---
# Companion to plot_comparable_boxplot(): uses the same sample-matching
# logic (geo_accession / title / barcode) and group releveling so it can
# be dropped straight into km_and_barplot_helpers.R.
#
# Requires: install.packages("cvequality")  (Feltz & Miller / Krishnamoorthy & Lee tests)

library(dplyr)
library(tidyr)
library(purrr)
library(cvequality)

calc_cv_stats <- function(expr_matrix,       # Symbols in rownames, samples in columns (TPM)
                          metadata_df,
                          target_genes,
                          ref_group = "Normal",
                          scale = c("raw", "log2"),
                          test_method = c("asymptotic", "mslr"),
                          n_resamples = 1000) {
  
  scale <- match.arg(scale)
  test_method <- match.arg(test_method)
  colnames(expr_matrix) <- gsub("^X", "", colnames(expr_matrix))
  
  # 1. Sample matching (identical logic to plot_comparable_boxplot)
  id_cols <- c("geo_accession", "title", "barcode")
  common_samples <- character(0)
  for (col in id_cols) {
    if (col %in% colnames(metadata_df)) {
      test_samples <- intersect(colnames(expr_matrix), metadata_df[[col]])
      if (length(test_samples) > 0) {
        common_samples <- test_samples
        temp_meta <- metadata_df
        rownames(temp_meta) <- temp_meta[[col]]
        metadata_df <- temp_meta
        message("Successfully matched using column: ", col)
        break
      }
    }
  }
  if (length(common_samples) == 0) stop("Error: No overlapping samples found.")
  
  # 2. Subset and relevel
  expr_sub <- expr_matrix[, common_samples, drop = FALSE]
  meta_sub <- metadata_df[common_samples, , drop = FALSE]
  meta_sub$group <- as.factor(meta_sub$group)
  if (ref_group %in% levels(meta_sub$group)) {
    meta_sub$group <- relevel(meta_sub$group, ref = ref_group)
  }
  
  present_genes <- intersect(target_genes, rownames(expr_sub))
  subset_expr <- expr_sub[present_genes, , drop = FALSE]
  
  # NOTE ON SCALE: CV = SD / mean is only meaningful on a true ratio scale
  # with a meaningful zero (raw TPM). log2(TPM+1) shifts/compresses the
  # scale and z-scoring forces mean = 0 (CV undefined), so both distort CV.
  # Default is "raw"; "log2" is offered only if you have a specific reason
  # to compare dispersion on the log scale (interpret with caution).
  if (scale == "log2") {
    subset_expr <- log2(subset_expr + 1)
  }
  
  data_long <- as.data.frame(t(subset_expr)) %>%
    mutate(Group = meta_sub$group) %>%
    pivot_longer(cols = all_of(present_genes), names_to = "Gene", values_to = "Val") %>%
    filter(!is.na(Val))
  
  # 3. Descriptive CV per Gene x Group
  cv_summary <- data_long %>%
    group_by(Gene, Group) %>%
    summarise(
      n      = dplyr::n(),
      mean   = mean(Val, na.rm = TRUE),
      sd     = sd(Val, na.rm = TRUE),
      CV     = sd / mean,      # unitless ratio; undefined/unstable if mean ~ 0
      CV_pct = CV * 100,
      .groups = "drop"
    ) %>%
    arrange(Gene, Group)
  
  # 4. Omnibus test for equality of CVs across ALL groups, per gene
  #    - "asymptotic": Feltz & Miller (1996) chi-squared test, fast, best for
  #      larger n per group.
  #    - "mslr": Krishnamoorthy & Lee (2014) modified signed log-likelihood
  #      ratio test with resampled critical values; more reliable for small
  #      or unequal group sizes.
  run_cv_test <- function(x, y) {
    if (test_method == "asymptotic") {
      res <- cvequality::asymptotic_test(x = x, y = y)
      tibble(statistic = res$D_AD, df = length(unique(y)) - 1,
             p.value = res$p_value, method = "Feltz-Miller asymptotic")
    } else {
      res <- cvequality::mslr_test(nr = n_resamples, x = x, y = y)
      tibble(statistic = res$MSLRT, df = NA_real_,
             p.value = res$p_value, method = "Krishnamoorthy-Lee MSLRT")
    }
  }
  
  cv_test_overall <- data_long %>%
    group_by(Gene) %>%
    group_modify(~ run_cv_test(.x$Val, as.character(.x$Group))) %>%
    ungroup() %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"))  # across genes
  
  # 5. Pairwise CV comparisons: ref_group vs each other group, per gene
  other_groups <- setdiff(levels(meta_sub$group), ref_group)
  
  cv_test_pairwise <- purrr::map_dfr(present_genes, function(g) {
    purrr::map_dfr(other_groups, function(grp) {
      sub <- data_long %>% filter(Gene == g, Group %in% c(ref_group, grp))
      out <- run_cv_test(sub$Val, as.character(sub$Group))
      out %>% mutate(Gene = g, Comparison = paste(ref_group, "vs", grp), .before = 1)
    })
  }) %>%
    group_by(Gene) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%  # BH within each gene's pairwise set
    ungroup()
  
  list(
    cv_summary       = cv_summary,       # n, mean, sd, CV, CV_pct per Gene x Group
    cv_test_overall  = cv_test_overall,  # k-group omnibus test per Gene
    cv_test_pairwise = cv_test_pairwise  # ref_group vs each other group, per Gene
  )
}

# ==============================================================================
# SECTION 7 -- EXAMPLE USAGE
# ==============================================================================
# cell_types <- c("B.cells.naive", "B.cells.memory", "T.cells.CD8", ...)   # as in the app
#
# --- Feature 1: gene-cell correlation heatmap, per-cohort ---
# df1 <- df_tf %>% dplyr::filter(!is.na(IRF1)) %>% dplyr::mutate(Expr_Val = normalize_vector(IRF1, "log2"))
# df1 <- build_cohort_subgroups(df1, subgroup_col = "Prior_Tx", show_all = TRUE, label_n = TRUE)
# corr_df <- compute_correlation_by_subgroup(df1, targets = cell_types, cell_types = cell_types)
# plot_corr_heatmap(corr_df)
#
# --- Feature 2: single-factor survival heatmap ---
# df2 <- df_tf %>% apply_admin_censoring(10) %>% dplyr::filter(!is.na(OS_Time), !is.na(OS_Status))
# df2 <- build_cohort_subgroups(df2, subgroup_col = "StudyID", show_all = TRUE, label_n = FALSE)
# hdf <- compute_survival_heatmap(df2, targets = c("IRF1","STAT1"), mode = "single", split_pct = 50)
# plot_survival_heatmap(hdf)
# compute_median_os_diff(assign_global_group(df2, "IRF1", 50) %>% dplyr::rename(Expr_Group = Group),
#                        group_col = "Expr_Group", max_years = 10)
#
# --- Feature 2/3: two-factor drill-down K-M ---
# strat <- prep_nested_stratification(df2, genes_A = "IRF1", target_B = "STAT1", splitA = 50, splitB = 50)
# plot_subset_km(strat, subset_A = "High", labelA = "IRF1", labelB = "STAT1", max_years = 10)
#
# --- Feature 3: comprehensive profiling, single factor ---
# res <- compute_factor_profiling_data(df2, genes_A = "IRF1", mode = "single", split_pct = 50)
# plot_km_by_group(res$df, group_col = "Plot_Group", title = paste("Survival by", res$labelA))  # from Figure9.Rmd
# compute_cox_stats_table(res)
# plot_factor_score_boxplot(res$df %>% dplyr::filter(Hierarchical_Subgroup == res$grps[1]), "single", res$labelA)
#
# --- One Biopsy per Patient, e.g. for the Cohort 1 (PHS000452_GSE91061) survival subset ---
# df_C1_raw <- df_pre %>% subset(StudyID == datasets[1])                  # no BiopsyTime filter yet
# count_duplicate_patients(df_C1_raw)                                    # how many patients have >1 row?
# df_C1     <- dedup_one_biopsy_per_patient(df_C1_raw, prefer = "Pre")     # keep 1 row/patient, prefer Pre
# count_duplicate_patients(df_C1)                                        # should now report n_duplicated = 0
# nrow(df_C1) == dplyr::n_distinct(df_C1_raw$PatientID)                  # sanity check: 1 row per patient
# ==============================================================================
