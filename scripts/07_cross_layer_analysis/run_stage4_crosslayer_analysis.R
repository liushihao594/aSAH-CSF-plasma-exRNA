options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("revision_stage1/R_library", .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

root <- "revision_stage4"
dir.create(root, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "scripts"), showWarnings = FALSE)
dir.create(file.path(root, "intermediate"), showWarnings = FALSE)
dir.create(file.path(root, "logs"), showWarnings = FALSE)
dir.create(file.path(root, "plots"), showWarnings = FALSE)

seed <- 20261121L
B <- 2000L
sd_floor <- 0.5
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

obj <- readRDS("revision_stage1/intermediate/prepared_stage1_inputs.rds")
inputs <- obj$inputs
strict <- read.csv("revision_stage1/stage1b_full_refit/STAGE1B_STRICT_CORE_ROBUST.csv", check.names = FALSE)
strict_m <- strict$feature[strict$layer == "miRNA"]
strict_g <- strict$feature[strict$layer == "GENCODE"]
stopifnot(length(strict_m) == 16L, length(strict_g) == 33L)

hist_dir <- file.path(
  "revision_stage0d/source_extract/TRACKB_REVISION_TRANSFER_PACKAGE",
  "11_ANALYSIS_RESULTS_AND_AUDITS/Primary_Results/Phase_1B_20260928_205553_Tables"
)
hist_scores <- read.csv(file.path(hist_dir, "TRACKB_CROSSLAYER_DIVERGENCE_SCORES_FINAL.csv"), check.names = FALSE)
hist_slopes <- read.csv(file.path(hist_dir, "TRACKB_CROSSLAYER_SUBJECT_SLOPES.csv"), check.names = FALSE)
hist_final <- read.csv(file.path(hist_dir, "TRACKB_CROSSLAYER_FINAL.csv"), check.names = FALSE)
hist_influence <- read.csv(file.path(hist_dir, "TRACKB_CROSSLAYER_INFLUENCE_SENSITIVITY.csv"), check.names = FALSE)

# Recover the original unfiltered feature universes. The historical LOSO procedure
# applied its prevalence filter inside each training fold, so restricting input to the
# 507/1,079 full-cohort retained sets is not numerically equivalent near the threshold.
miRNA_all <- as.matrix(read.delim(gzfile(file.path(
  "revision_stage0/source_extract/TRACKB_SOURCE_DATA_PACKAGE/01_GEO_processed_expression",
  "GSE121868_KJENS1-aSAH_Project-2016-10-13_exceRpt_miRNA_ReadsPerMillion.txt.gz"
)), row.names = 1, check.names = FALSE))
GENCODE_all <- as.matrix(read.delim(gzfile(file.path(
  "revision_stage0/source_extract/TRACKB_SOURCE_DATA_PACKAGE/01_GEO_processed_expression",
  "GSE121868_KJENS1-aSAH_Project-2016-10-13_exceRpt_gencode_ReadsPerMillion.txt.gz"
)), row.names = 1, check.names = FALSE))
stopifnot(ncol(miRNA_all) == 523L, ncol(GENCODE_all) == 523L,
          setequal(colnames(miRNA_all), colnames(GENCODE_all)))
input_cache <- new.env(parent = emptyenv())

prepare_unfiltered_input <- function(layer, dataset, omitted = character()) {
  key <- paste(layer, dataset, paste(sort(omitted), collapse = ";"), sep = "|")
  if (exists(key, envir = input_cache, inherits = FALSE)) return(get(key, envir = input_cache))
  frozen <- inputs[[paste0(layer, "_", dataset)]]
  pp <- frozen$pairs[!frozen$pairs$subject %in% omitted, , drop = FALSE]
  full <- if (layer == "miRNA") miRNA_all else GENCODE_all
  base_eligible <- if (layer == "miRNA") rep(TRUE, nrow(full)) else sub("^.*:", "", rownames(full)) != "miRNA"
  csf_cols <- match(pp$CSF_expression_id, colnames(full))
  plasma_cols <- match(pp$plasma_expression_id, colnames(full))
  stopifnot(!anyNA(csf_cols), !anyNA(plasma_cols))
  all_cols <- c(csf_cols, plasma_cols)
  total_positive <- rowSums(full[, all_cols, drop = FALSE] > 1)
  candidate <- rep(FALSE, nrow(full))
  for (id in unique(pp$subject)) {
    held <- pp$subject == id
    held_cols <- c(csf_cols[held], plasma_cols[held])
    training_n <- length(all_cols) - length(held_cols)
    train_positive <- total_positive - rowSums(full[, held_cols, drop = FALSE] > 1)
    candidate <- candidate | (base_eligible & train_positive / training_n >= 0.2)
  }
  rpm <- full[candidate, , drop = FALSE]
  csf <- rpm[, csf_cols, drop = FALSE]
  plasma <- rpm[, plasma_cols, drop = FALSE]
  raw <- log2(csf + 1) - log2(plasma + 1)
  reference_center <- colMeans(frozen$raw[, match(pp$pair_id, frozen$pairs$pair_id), drop = FALSE])
  clr <- sweep(raw, 2L, reference_center, "-")
  out <- list(layer = layer, dataset = dataset, pairs = pp, raw = raw, clr = clr, csf = csf, plasma = plasma)
  assign(key, out, envir = input_cache)
  out
}

crossfit_layer <- function(inp, gradient_name, eligible_features, layer, scheme, omitted = character()) {
  pp <- inp$pairs[!inp$pairs$subject %in% omitted, , drop = FALSE]
  keep_col <- !inp$pairs$subject %in% omitted
  grad <- inp[[gradient_name]][, keep_col, drop = FALSE]
  csf <- inp$csf[, keep_col, drop = FALSE]
  plasma <- inp$plasma[, keep_col, drop = FALSE]
  stopifnot(nrow(pp) == ncol(grad), identical(rownames(grad), rownames(csf)), identical(rownames(grad), rownames(plasma)))
  eligible <- rownames(grad) %in% eligible_features
  subjects <- sort(unique(pp$subject))
  values <- rep(NA_real_, nrow(pp))
  weight_rows <- vector("list", length(subjects))
  for (ii in seq_along(subjects)) {
    id <- subjects[ii]
    tr <- pp$subject != id
    ho <- !tr
    training_rpm <- cbind(csf[, tr, drop = FALSE], plasma[, tr, drop = FALSE])
    use <- eligible & rowMeans(training_rpm > 1) >= 0.2
    yy <- grad[use, , drop = FALSE]
    train_subjects <- subjects[subjects != id]
    td <- vapply(train_subjects, function(d) rowMeans(yy[, pp$subject == d, drop = FALSE]), numeric(sum(use)))
    if (is.null(dim(td))) td <- matrix(td, ncol = 1L)
    direction <- sign(rowMeans(td))
    scale <- pmax(apply(yy[, tr, drop = FALSE], 1L, sd), sd_floor)
    values[ho] <- colMeans(sweep(sweep(yy[, ho, drop = FALSE], 1L, scale, "/"), 1L, direction, "*"))
    weight_rows[[ii]] <- data.frame(
      scheme = scheme, layer = layer, heldout_subject = id,
      feature_id = rownames(yy), direction = direction, scale = scale,
      training_subject_N = length(train_subjects), training_pair_N = sum(tr),
      stringsAsFactors = FALSE
    )
  }
  list(score = values, weights = do.call(rbind, weight_rows), pairs = pp)
}

crossfit_two_layers <- function(dataset, scale, scheme = "ORIGINAL_FULL_FEATURE", omitted = character(), strict_only = FALSE) {
  im <- prepare_unfiltered_input("miRNA", dataset, omitted)
  ig <- prepare_unfiltered_input("GENCODE", dataset, omitted)
  stopifnot(identical(im$pairs$pair_id, ig$pairs$pair_id))
  grad_name <- tolower(scale)
  elig_m <- if (strict_only) strict_m else rownames(im[[grad_name]])
  elig_g <- if (strict_only) strict_g else rownames(ig[[grad_name]])[sub("^.*:", "", rownames(ig[[grad_name]])) != "miRNA"]
  ms <- crossfit_layer(im, grad_name, elig_m, "miRNA", scheme, character())
  gs <- crossfit_layer(ig, grad_name, elig_g, "GENCODE", scheme, character())
  stopifnot(identical(ms$pairs$pair_id, gs$pairs$pair_id))
  pp <- ms$pairs
  scores <- data.frame(
    subject = pp$subject, collection_label = pp$collection_label,
    day = pp$centered_label + 4, pair_id = pp$pair_id,
    miRNA_score = ms$score, GENCODE_score = gs$score,
    dataset = dataset, scale = scale, scoring_scheme = scheme,
    stringsAsFactors = FALSE
  )
  list(scores = scores, weights = rbind(ms$weights, gs$weights))
}

patient_slopes <- function(sc, id_col = "subject") {
  ids <- unique(sc[[id_col]])
  do.call(rbind, lapply(ids, function(id) {
    z <- sc[sc[[id_col]] == id, , drop = FALSE]
    dx <- z$day - mean(z$day)
    data.frame(
      id = id,
      original_subject = z$subject[1],
      n_labels = nrow(z),
      miRNA_slope = sum(dx * z$miRNA_score) / sum(dx^2),
      GENCODE_slope = sum(dx * z$GENCODE_score) / sum(dx^2),
      stringsAsFactors = FALSE
    )
  }))
}

score_metrics <- function(sc, id_col = "subject") {
  ids <- sc[[id_col]]
  mc <- sc$miRNA_score - ave(sc$miRNA_score, ids, FUN = mean)
  gc <- sc$GENCODE_score - ave(sc$GENCODE_score, ids, FUN = mean)
  d <- data.frame(miRNA_score = sc$miRNA_score, GENCODE_score = sc$GENCODE_score,
                  id = factor(ids), day = factor(sc$day))
  mr <- residuals(lm(miRNA_score ~ id + day, d))
  gr <- residuals(lm(GENCODE_score ~ id + day, d))
  sl <- patient_slopes(sc, id_col)
  same <- sum(sign(sl$miRNA_slope) == sign(sl$GENCODE_slope))
  c(
    pair_Pearson = safe_cor(sc$miRNA_score, sc$GENCODE_score),
    pair_Spearman = safe_cor(sc$miRNA_score, sc$GENCODE_score, "spearman"),
    within_subject_Pearson = safe_cor(mc, gc),
    within_subject_Spearman = safe_cor(mc, gc, "spearman"),
    subject_label_adjusted_Pearson = safe_cor(mr, gr),
    subject_label_adjusted_Spearman = safe_cor(mr, gr, "spearman"),
    slope_Pearson = safe_cor(sl$miRNA_slope, sl$GENCODE_slope),
    slope_Spearman = safe_cor(sl$miRNA_slope, sl$GENCODE_slope, "spearman"),
    slope_same_direction_N = same,
    slope_same_direction_fraction = same / nrow(sl)
  )
}

bootstrap_scores <- function(sc, branch, B = 2000L, seed_offset = 0L) {
  subjects <- sort(unique(sc$subject))
  set.seed(seed + seed_offset)
  out <- vector("list", B)
  for (b in seq_len(B)) {
    draw <- sample(subjects, length(subjects), replace = TRUE)
    pieces <- lapply(seq_along(draw), function(k) {
      z <- sc[sc$subject == draw[k], , drop = FALSE]
      z$boot_subject <- paste0(draw[k], "__copy", k)
      z
    })
    bs <- do.call(rbind, pieces)
    met <- score_metrics(bs, "boot_subject")
    out[[b]] <- data.frame(branch = branch, replicate = b, metric = names(met), estimate = as.numeric(met))
  }
  do.call(rbind, out)
}

summarize_branch <- function(sc, branch, boot) {
  point <- score_metrics(sc)
  slopes <- patient_slopes(sc)
  do.call(rbind, lapply(names(point), function(metric) {
    vals <- boot$estimate[boot$branch == branch & boot$metric == metric]
    vals <- vals[is.finite(vals)]
    data.frame(
      branch = branch, dataset = sc$dataset[1], scale = sc$scale[1], scoring_scheme = sc$scoring_scheme[1],
      N_subjects = length(unique(sc$subject)), N_pairs = nrow(sc), metric = metric,
      estimate = unname(point[metric]), CI_lower = if (length(vals)) unname(quantile(vals, .025)) else NA_real_,
      CI_upper = if (length(vals)) unname(quantile(vals, .975)) else NA_real_,
      bootstrap_requested = B, bootstrap_valid = length(vals), bootstrap_invalid = B - length(vals),
      exact_binomial_P = if (metric == "slope_same_direction_N") binom.test(unname(point[metric]), nrow(slopes), .5)$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

# Primary cross-fitted scores.
primary_fits <- list()
primary_weights <- list()
for (dataset in c("FULL", "QC_PASS_ONLY")) for (scale in c("RAW", "CLR")) {
  key <- paste(dataset, scale, sep = "_")
  fit <- crossfit_two_layers(dataset, scale, "ORIGINAL_FULL_FEATURE")
  primary_fits[[key]] <- fit$scores
  primary_weights[[key]] <- fit$weights
}
primary_scores <- do.call(rbind, primary_fits)
rownames(primary_scores) <- NULL

# Exact original reproduction gate.
raw_full <- primary_fits[["FULL_RAW"]]
hist_key <- paste(hist_scores$subject, hist_scores$Dxx, sep = "|")
new_key <- paste(raw_full$subject, raw_full$collection_label, sep = "|")
ord <- match(hist_key, new_key)
stopifnot(!anyNA(ord))
score_diff_m <- max(abs(hist_scores$miRNA_score - raw_full$miRNA_score[ord]))
score_diff_g <- max(abs(hist_scores$GENCODE_score - raw_full$GENCODE_score[ord]))
new_metrics <- score_metrics(raw_full)
new_slopes <- patient_slopes(raw_full)
hs_ord <- match(hist_slopes$subject, new_slopes$original_subject)
slope_diff_m <- max(abs(hist_slopes$miRNA_slope - new_slopes$miRNA_slope[hs_ord]))
slope_diff_g <- max(abs(hist_slopes$GENCODE_slope - new_slopes$GENCODE_slope[hs_ord]))
metric_map <- c(pair_Pearson = "pair_Pearson", pair_Spearman = "pair_Spearman",
                slope_Pearson = "slope_Pearson", slope_Spearman = "slope_Spearman",
                slope_same_direction_N = "slope_same_direction_N",
                within_subject_Pearson = "within_subject_Pearson", within_subject_Spearman = "within_subject_Spearman",
                subject_label_adjusted_Pearson = "subject_label_adjusted_Pearson",
                subject_label_adjusted_Spearman = "subject_label_adjusted_Spearman")
repro <- data.frame(
  item = c("miRNA_score_max_abs_diff", "GENCODE_score_max_abs_diff", "miRNA_subject_slope_max_abs_diff", "GENCODE_subject_slope_max_abs_diff",
           names(metric_map)),
  historical = c(0, 0, 0, 0, as.numeric(hist_final[1, unname(metric_map)])),
  reproduced = c(score_diff_m, score_diff_g, slope_diff_m, slope_diff_g, as.numeric(new_metrics[names(metric_map)])),
  stringsAsFactors = FALSE
)
repro$absolute_difference <- abs(repro$historical - repro$reproduced)
repro$tolerance <- 1e-12
repro$status <- ifelse(repro$absolute_difference <= repro$tolerance, "PASS", "FAIL")
write_csv(repro, file.path(root, "STAGE4_ORIGINAL_SCORE_REPRODUCTION.csv"))
if (any(repro$status != "PASS")) stop("Original cross-layer reproduction gate failed")

# Primary cluster bootstrap.
boot_list <- list()
summary_list <- list()
offset <- 0L
for (key in names(primary_fits)) {
  branch <- paste0("PRIMARY_", key)
  bt <- bootstrap_scores(primary_fits[[key]], branch, B, offset)
  boot_list[[branch]] <- bt
  summary_list[[branch]] <- summarize_branch(primary_fits[[key]], branch, bt)
  offset <- offset + 1000L
}

# Subject12 fixed omission and complete refolding for FULL RAW/CLR.
subject12 <- "EXR-KJENS1Sample12-DO"
s12_rows <- list()
s12_scores <- list()
for (scale in c("RAW", "CLR")) {
  base <- primary_fits[[paste0("FULL_", scale)]]
  fixed <- base[base$subject != subject12, , drop = FALSE]
  fixed$scoring_scheme <- "ORIGINAL_FULL_FEATURE_FIXED_SCORE_OMIT_SUBJECT12"
  refold_fit <- crossfit_two_layers("FULL", scale, "ORIGINAL_FULL_FEATURE_REFOLD_SUBJECT12", omitted = subject12)
  refold <- refold_fit$scores
  for (type in c("FULL", "FIXED_SCORE_OMISSION", "COMPLETE_REFOLDING")) {
    sc <- switch(type, FULL = base, FIXED_SCORE_OMISSION = fixed, COMPLETE_REFOLDING = refold)
    branch <- paste("SUBJECT12", scale, type, sep = "_")
    bt <- bootstrap_scores(sc, branch, B, offset)
    boot_list[[branch]] <- bt
    sm <- summarize_branch(sc, branch, bt)
    summary_list[[branch]] <- sm
    s12_rows[[branch]] <- sm
    s12_scores[[branch]] <- sc
    offset <- offset + 1000L
  }
}
s12_summary <- do.call(rbind, s12_rows)
s12_full_ref <- merge(
  s12_summary[s12_summary$branch %in% c("SUBJECT12_RAW_FULL", "SUBJECT12_CLR_FULL"), c("scale", "metric", "estimate")],
  s12_summary[!s12_summary$branch %in% c("SUBJECT12_RAW_FULL", "SUBJECT12_CLR_FULL"), ],
  by = c("scale", "metric"), suffixes = c("_FULL", "")
)
s12_full_ref$change_from_FULL <- s12_full_ref$estimate - s12_full_ref$estimate_FULL
write_csv(s12_full_ref, file.path(root, "STAGE4_SUBJECT12_SENSITIVITY.csv"))

# All-subject fixed-score LOO; identify largest non-subject12 slope-Pearson influence per scale.
loo_rows <- list()
for (scale in c("RAW", "CLR")) {
  base <- primary_fits[[paste0("FULL_", scale)]]
  base_met <- score_metrics(base)
  for (id in sort(unique(base$subject))) {
    z <- base[base$subject != id, , drop = FALSE]
    met <- score_metrics(z)
    loo_rows[[length(loo_rows) + 1L]] <- data.frame(
      scale = scale, omitted_subject = id, analysis_type = "FIXED_SCORE_OMISSION",
      metric = names(met), estimate = as.numeric(met), full_estimate = as.numeric(base_met[names(met)]),
      change = as.numeric(met - base_met[names(met)]), direction_flip = sign(met) != sign(base_met[names(met)]),
      stringsAsFactors = FALSE
    )
  }
}
loo <- do.call(rbind, loo_rows)
top_non12 <- do.call(rbind, lapply(c("RAW", "CLR"), function(scale) {
  z <- loo[loo$scale == scale & loo$metric == "slope_Pearson" & loo$omitted_subject != subject12, ]
  z[order(-abs(z$change), z$omitted_subject)[1], , drop = FALSE]
}))
for (i in seq_len(nrow(top_non12))) {
  scale <- top_non12$scale[i]
  id <- top_non12$omitted_subject[i]
  refold <- crossfit_two_layers("FULL", scale, paste0("ORIGINAL_FULL_FEATURE_REFOLD_", id), omitted = id)$scores
  met <- score_metrics(refold)
  base_met <- score_metrics(primary_fits[[paste0("FULL_", scale)]])
  loo <- rbind(loo, data.frame(
    scale = scale, omitted_subject = id, analysis_type = "COMPLETE_REFOLDING_TOP_NON12",
    metric = names(met), estimate = as.numeric(met), full_estimate = as.numeric(base_met[names(met)]),
    change = as.numeric(met - base_met[names(met)]), direction_flip = sign(met) != sign(base_met[names(met)]),
    stringsAsFactors = FALSE
  ))
}
# Add subject12 complete-refolding metrics explicitly.
for (scale in c("RAW", "CLR")) {
  refold <- s12_scores[[paste0("SUBJECT12_", scale, "_COMPLETE_REFOLDING")]]
  met <- score_metrics(refold); base_met <- score_metrics(primary_fits[[paste0("FULL_", scale)]])
  loo <- rbind(loo, data.frame(
    scale = scale, omitted_subject = subject12, analysis_type = "COMPLETE_REFOLDING_SUBJECT12",
    metric = names(met), estimate = as.numeric(met), full_estimate = as.numeric(base_met[names(met)]),
    change = as.numeric(met - base_met[names(met)]), direction_flip = sign(met) != sign(base_met[names(met)]),
    stringsAsFactors = FALSE
  ))
}
write_csv(loo, file.path(root, "STAGE4_LOO_INFLUENCE.csv"))

# Strict-core post-selection exploratory scores and bootstrap.
strict_fits <- list(); strict_weights <- list()
for (dataset in c("FULL", "QC_PASS_ONLY")) for (scale in c("RAW", "CLR")) {
  key <- paste(dataset, scale, sep = "_")
  fit <- crossfit_two_layers(dataset, scale, "STRICT_CORE_POST_SELECTION_EXPLORATORY", strict_only = TRUE)
  strict_fits[[key]] <- fit$scores; strict_weights[[key]] <- fit$weights
  branch <- paste0("STRICT_CORE_", key)
  bt <- bootstrap_scores(fit$scores, branch, B, offset)
  boot_list[[branch]] <- bt
  summary_list[[branch]] <- summarize_branch(fit$scores, branch, bt)
  offset <- offset + 1000L
}

all_boot <- do.call(rbind, boot_list)
all_summary <- do.call(rbind, summary_list)
write_csv(primary_scores, file.path(root, "STAGE4_RAW_CLR_CROSSFITTED_SCORES.csv"))
write_csv(all_boot, file.path(root, "STAGE4_PATIENT_BOOTSTRAP.csv"))
write_csv(all_summary, file.path(root, "STAGE4_CORRELATION_SUMMARY.csv"))
within <- all_summary[all_summary$metric %in% c("within_subject_Pearson", "within_subject_Spearman",
                                               "subject_label_adjusted_Pearson", "subject_label_adjusted_Spearman"), ]
write_csv(within, file.path(root, "STAGE4_WITHIN_SUBJECT_ANALYSIS.csv"))
strict_summary <- all_summary[grepl("^STRICT_CORE_", all_summary$branch), ]
strict_summary$interpretation <- "POST_SELECTION_EXPLORATORY; feature set selected in same 30-subject cohort"
write_csv(strict_summary, file.path(root, "STAGE4_STRICT_CORE_EXPLORATORY.csv"))

# Store scores and weights for full auditability.
write_csv(do.call(rbind, lapply(names(strict_fits), function(k) strict_fits[[k]])),
          file.path(root, "intermediate", "STAGE4_STRICT_CORE_SCORES.csv"))
saveRDS(list(primary_scores = primary_fits, primary_weights = primary_weights,
             strict_scores = strict_fits, strict_weights = strict_weights,
             subject12_scores = s12_scores, summary = all_summary, bootstrap = all_boot, loo = loo,
             seed = seed, B = B, sd_floor = sd_floor),
        file.path(root, "intermediate", "STAGE4_RESULTS_OBJECT.rds"))

# Exploratory figures.
save_both <- function(name, p, width, height) {
  ggsave(file.path(root, "plots", paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(root, "plots", paste0(name, ".svg")), p, width = width, height = height)
}

p1 <- ggplot(primary_scores, aes(miRNA_score, GENCODE_score, color = subject)) +
  geom_point(alpha = .65, size = 1.4) + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = .7) +
  facet_grid(dataset ~ scale, scales = "free") + theme_bw() + theme(legend.position = "none") +
  labs(title = "Cross-fitted paired-label scores", x = "miRNA held-out score", y = "GENCODE non-miRNA held-out score",
       caption = "Repeated observations: FULL n=30 subjects/153 pairs; QC PASS-only n=30/141; points are not independent patients")
save_both("STAGE4_FIGURE_A_PAIRED_SCORE_SCATTER", p1, 10, 7)

slope_plot <- do.call(rbind, lapply(c("RAW", "CLR"), function(scale) {
  z <- patient_slopes(primary_fits[[paste0("FULL_", scale)]])
  z$scale <- scale; z
}))
p2 <- ggplot(slope_plot, aes(miRNA_slope, GENCODE_slope)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + geom_point(color = "#0072B2", size = 2) +
  geom_text(data = slope_plot[slope_plot$original_subject == subject12, ], aes(label = "subject12"),
            hjust = 0, vjust = -1, nudge_x = .01, color = "#D55E00") +
  geom_smooth(method = "lm", se = FALSE, color = "grey30") + facet_wrap(~scale, scales = "free") + theme_bw() +
  labs(title = "Subject-specific cross-layer score slopes", x = "miRNA slope", y = "GENCODE slope",
       caption = "30 subjects; each slope uses the subject's 3-7 observed collection labels")
save_both("STAGE4_FIGURE_B_PATIENT_SLOPES", p2, 9, 5)

boot_plot <- all_boot[all_boot$branch %in% paste0("PRIMARY_FULL_", c("RAW", "CLR")) &
                        all_boot$metric %in% c("pair_Pearson", "pair_Spearman", "slope_Pearson", "slope_Spearman"), ]
p3 <- ggplot(boot_plot, aes(estimate, fill = branch)) + geom_histogram(bins = 45, alpha = .65, position = "identity") +
  geom_vline(xintercept = 0, linetype = 2) + facet_wrap(~metric, scales = "free_y") + theme_bw() +
  labs(title = "Patient-level cluster bootstrap distributions", x = "Correlation", y = "Bootstrap replicates",
       caption = paste0("2,000 patient-resampling replicates per branch; seed ", seed))
save_both("STAGE4_FIGURE_C_BOOTSTRAP_DISTRIBUTIONS", p3, 11, 7)

loo_plot <- loo[loo$analysis_type == "FIXED_SCORE_OMISSION" & loo$metric %in% c("pair_Pearson", "slope_Pearson"), ]
loo_plot$subject_short <- sub("EXR-KJENS1Sample", "S", sub("-DO", "", loo_plot$omitted_subject))
p4 <- ggplot(loo_plot, aes(subject_short, estimate, color = metric, group = metric)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + geom_point() + facet_wrap(~scale, ncol = 1) + theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Fixed-score leave-one-subject-out influence", x = "Omitted subject", y = "Correlation",
       caption = "Fixed held-out scores; complete refolding separately retained for subject12 and the strongest non-subject12 influence")
save_both("STAGE4_FIGURE_D_LOO_INFLUENCE", p4, 12, 8)

forest <- all_summary[grepl("^PRIMARY_(FULL|QC_PASS_ONLY)_", all_summary$branch) &
                        all_summary$metric %in% c("pair_Pearson", "pair_Spearman", "slope_Pearson", "slope_Spearman"), ]
p5 <- ggplot(forest, aes(estimate, branch, xmin = CI_lower, xmax = CI_upper, color = metric)) +
  geom_vline(xintercept = 0, linetype = 2) + geom_errorbarh(height = .18, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) + facet_wrap(~metric, scales = "free_y") + theme_bw() +
  labs(title = "FULL and QC PASS-only cross-layer correlations", x = "Correlation (patient-cluster bootstrap 95% CI)", y = NULL)
save_both("STAGE4_FIGURE_E_FULL_QC_FOREST", p5, 12, 8)

within_plot <- within[grepl("^PRIMARY_", within$branch), ]
p6 <- ggplot(within_plot, aes(estimate, branch, xmin = CI_lower, xmax = CI_upper, color = metric)) +
  geom_vline(xintercept = 0, linetype = 2) + geom_errorbarh(height = .15, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) + facet_wrap(~metric, scales = "free_y") + theme_bw() +
  labs(title = "Within-subject and subject/label-adjusted correlations", x = "Correlation (patient-cluster bootstrap 95% CI)", y = NULL)
save_both("STAGE4_FIGURE_F_WITHIN_ADJUSTED", p6, 12, 8)

compare <- all_summary[(grepl("^PRIMARY_FULL_", all_summary$branch) | grepl("^STRICT_CORE_FULL_", all_summary$branch)) &
                         all_summary$metric %in% c("pair_Pearson", "slope_Pearson"), ]
p7 <- ggplot(compare, aes(estimate, branch, xmin = CI_lower, xmax = CI_upper, color = scoring_scheme)) +
  geom_vline(xintercept = 0, linetype = 2) + geom_errorbarh(height = .18) + geom_point(size = 2) +
  facet_wrap(~metric, scales = "free_y") + theme_bw() +
  labs(title = "All-feature and strict-core exploratory scores", x = "Correlation (patient-cluster bootstrap 95% CI)", y = NULL,
       caption = "Strict-core features were selected in the same cohort and are post-selection exploratory")
save_both("STAGE4_FIGURE_G_STRICT_CORE_COMPARISON", p7, 11, 6)

writeLines(capture.output(sessionInfo()), file.path(root, "logs", "STAGE4_sessionInfo.txt"))
writeLines(c(
  paste("completed", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste("bootstrap_seed", seed), paste("bootstrap_replicates_per_branch", B), paste("sd_floor", sd_floor),
  paste("original_reproduction_max_abs_difference", max(repro$absolute_difference)),
  paste("primary_score_rows", nrow(primary_scores)), paste("bootstrap_rows", nrow(all_boot)),
  paste("subject12", subject12), paste("top_non12_RAW", top_non12$omitted_subject[top_non12$scale == "RAW"]),
  paste("top_non12_CLR", top_non12$omitted_subject[top_non12$scale == "CLR"])
), file.path(root, "logs", "STAGE4_execution.log"))

cat(sprintf("STAGE4_ANALYSIS_COMPLETE scores=%d bootstrap_rows=%d repro_max=%.3g\n",
            nrow(primary_scores), nrow(all_boot), max(repro$absolute_difference)))
