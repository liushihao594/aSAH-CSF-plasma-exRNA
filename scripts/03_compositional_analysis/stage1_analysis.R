options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[[1]] else "pilot"
allowed_phases <- c("pilot", "full", "finalize")
if (!phase %in% allowed_phases) stop("phase must be one of: ", paste(allowed_phases, collapse = ", "))

# Keep paths relative because this Windows R build cannot reliably convert the
# Unicode project root to the active C locale. The launcher fixes the cwd.
project_root <- "."
out_dir <- file.path(project_root, "revision_stage1")
dir.create(file.path(out_dir, "intermediate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "diagnostic_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "provenance"), recursive = TRUE, showWarnings = FALSE)

stage_seed <- 20261120L
set.seed(stage_seed)

suppressPackageStartupMessages({
  library(limma)
  library(lme4)
  library(lmerTest)
})

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

direction <- function(x, tol = .Machine$double.eps^0.5) {
  ifelse(is.na(x), NA_character_, ifelse(x > tol, "positive", ifelse(x < -tol, "negative", "zero")))
}

bh_fixed <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  out[ok] <- p.adjust(p[ok], method = "BH", n = length(p))
  out
}

read_arm <- function(layer, arm) {
  readRDS(file.path(project_root, "revision_stage0c", "data", paste0(layer, "_", arm, "_paired_data.rds")))
}

make_clr_gradient <- function(csf, plasma, reference_features, pseudocount = 1) {
  stopifnot(identical(rownames(csf), rownames(plasma)))
  idx <- match(reference_features, rownames(csf))
  if (anyNA(idx)) stop("Reference features not found: ", paste(reference_features[is.na(idx)], collapse = ";"))
  lcsf <- log2(csf + pseudocount)
  lplasma <- log2(plasma + pseudocount)
  sweep(lcsf, 2L, colMeans(lcsf[idx, , drop = FALSE]), "-") -
    sweep(lplasma, 2L, colMeans(lplasma[idx, , drop = FALSE]), "-")
}

make_reference_centered_from_raw <- function(raw_gradient, reference_features) {
  idx <- match(reference_features, rownames(raw_gradient))
  if (anyNA(idx)) stop("Reference feature mismatch")
  sweep(raw_gradient, 2L, colMeans(raw_gradient[idx, , drop = FALSE]), "-")
}

fit_limma_linear <- function(mat, pairs) {
  stopifnot(ncol(mat) == nrow(pairs))
  design <- model.matrix(~ centered_label, data = pairs)
  dup <- duplicateCorrelation(mat, design = design, block = pairs$subject)
  fit <- lmFit(mat, design = design, block = pairs$subject, correlation = dup$consensus.correlation)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  coef_name <- "centered_label"
  j <- match(coef_name, colnames(fit$coefficients))
  if (is.na(j)) stop("centered_label coefficient missing")
  se <- fit$stdev.unscaled[, j] * fit$sigma
  p <- fit$p.value[, j]
  data.frame(
    feature = rownames(mat),
    slope = fit$coefficients[, j],
    SE = se,
    CI_lower = fit$coefficients[, j] - 1.96 * se,
    CI_upper = fit$coefficients[, j] + 1.96 * se,
    t = fit$t[, j],
    df = fit$df.total,
    P = p,
    FDR = p.adjust(p, method = "BH"),
    duplicate_correlation = dup$consensus.correlation,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

subject_slope_sd <- function(mat, pairs) {
  subjects <- unique(pairs$subject)
  slopes <- vapply(subjects, function(s) {
    idx <- which(pairs$subject == s)
    x <- pairs$centered_label[idx]
    w <- (x - mean(x)) / sum((x - mean(x))^2)
    as.numeric(mat[, idx, drop = FALSE] %*% w)
  }, numeric(nrow(mat)))
  apply(slopes, 1L, sd, na.rm = TRUE)
}

extract_conv <- function(fit) {
  msgs <- fit@optinfo$conv$lme4$messages
  opt <- fit@optinfo$conv$opt
  list(
    opt_code = if (length(opt)) paste(opt, collapse = ";") else "0",
    messages = if (length(msgs)) paste(msgs, collapse = " | ") else "",
    bad = (length(opt) && any(as.numeric(opt) != 0)) || length(msgs) > 0L
  )
}

fit_lmm_once <- function(dat, correlated = TRUE, optimizer = "nloptwrap") {
  f <- if (correlated) {
    y ~ label_centered + (1 + label_centered | subject)
  } else {
    y ~ label_centered + (1 | subject) + (0 + label_centered | subject)
  }
  warnings <- character()
  optimizer_control <- if (optimizer == "nloptwrap") list(maxeval = 200000) else list(maxfun = 200000)
  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        f, data = dat, REML = TRUE,
        control = lmerControl(
          optimizer = optimizer,
          optCtrl = optimizer_control,
          check.conv.singular = lme4:::.makeCC(action = "ignore", tol = 1e-4)
        )
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(fit = NULL, error = conditionMessage(fit), warnings = paste(unique(warnings), collapse = " | ")))
  }
  list(fit = fit, error = "", warnings = paste(unique(warnings), collapse = " | "))
}

fit_feature_lmm <- function(y, pairs) {
  dat <- data.frame(
    y = as.numeric(y),
    label_centered = as.numeric(pairs$centered_label),
    subject = factor(pairs$subject)
  )
  attempts <- list(
    list(correlated = TRUE, optimizer = "nloptwrap", trigger = "initial"),
    list(correlated = TRUE, optimizer = "bobyqa", trigger = "retry_correlated"),
    list(correlated = FALSE, optimizer = "nloptwrap", trigger = "fallback_uncorrelated"),
    list(correlated = FALSE, optimizer = "bobyqa", trigger = "retry_uncorrelated")
  )
  attempt_log <- character()
  final <- NULL
  final_meta <- NULL
  correlated_needs_fallback <- FALSE

  for (k in seq_along(attempts)) {
    a <- attempts[[k]]
    if (k == 2L && length(attempt_log) == 1L) {
      # run retry only because the first related fit was invalid
    }
    if (k == 3L && !correlated_needs_fallback) next
    if (k == 4L && !correlated_needs_fallback) next
    res <- fit_lmm_once(dat, correlated = a$correlated, optimizer = a$optimizer)
    if (is.null(res$fit)) {
      attempt_log <- c(attempt_log, paste(a$trigger, a$optimizer, "ERROR", res$error, res$warnings, sep = ":"))
      if (a$correlated && k == 2L) correlated_needs_fallback <- TRUE
      next
    }
    conv <- extract_conv(res$fit)
    singular <- isSingular(res$fit, tol = 1e-4)
    sm <- tryCatch(summary(res$fit)$coefficients, error = function(e) NULL)
    coef_ok <- !is.null(sm) && "label_centered" %in% rownames(sm) && all(is.finite(sm["label_centered", c("Estimate", "Std. Error")]))
    invalid <- conv$bad || !coef_ok
    attempt_log <- c(attempt_log, paste(a$trigger, a$optimizer,
      paste0("conv_bad=", conv$bad), paste0("singular=", singular),
      paste0("warning=", res$warnings), paste0("message=", conv$messages), sep = ":"))

    if (a$correlated) {
      if (!invalid && !singular) {
        final <- res$fit
        final_meta <- list(a = a, conv = conv, singular = singular, warnings = res$warnings)
        break
      }
      if (k == 2L) correlated_needs_fallback <- TRUE
    } else {
      if (!invalid && !singular) {
        final <- res$fit
        final_meta <- list(a = a, conv = conv, singular = singular, warnings = res$warnings)
        break
      }
    }
  }

  if (is.null(final)) {
    return(list(
      result = data.frame(
        slope = NA_real_, SE = NA_real_, df = NA_real_, t = NA_real_, P = NA_real_,
        CI_lower = NA_real_, CI_upper = NA_real_, model_valid = FALSE,
        final_structure = "MODEL_NOT_EVALUABLE", optimizer = NA_character_,
        singular = NA, convergence_ok = FALSE, n_attempts = length(attempt_log),
        random_intercept_variance = NA_real_, random_slope_variance = NA_real_,
        intercept_slope_correlation = NA_real_, residual_variance = NA_real_,
        max_abs_scaled_gradient = NA_real_, convergence_messages = paste(attempt_log, collapse = " || "),
        stringsAsFactors = FALSE
      ),
      fit = NULL
    ))
  }

  co <- summary(final)$coefficients
  row <- co["label_centered", , drop = FALSE]
  p_col <- grep("Pr\\(", colnames(co), value = TRUE)[1]
  df_col <- if ("df" %in% colnames(co)) "df" else NA_character_
  t_col <- grep("t value", colnames(co), value = TRUE)[1]
  vc <- as.data.frame(VarCorr(final))
  int_var <- vc$vcov[vc$var1 == "(Intercept)" & is.na(vc$var2)][1]
  slope_var <- vc$vcov[vc$var1 == "label_centered" & is.na(vc$var2)][1]
  corr <- vc$sdcor[!is.na(vc$var2) &
    ((vc$var1 == "(Intercept)" & vc$var2 == "label_centered") |
       (vc$var1 == "label_centered" & vc$var2 == "(Intercept)"))][1]
  resid_var <- vc$vcov[vc$grp == "Residual"][1]
  derivs <- final@optinfo$derivs
  max_grad <- NA_real_
  if (!is.null(derivs$gradient) && !is.null(derivs$Hessian)) {
    max_grad <- tryCatch(max(abs(solve(chol(derivs$Hessian), derivs$gradient))), error = function(e) NA_real_)
  }
  est <- unname(row[1, "Estimate"])
  se <- unname(row[1, "Std. Error"])
  list(
    result = data.frame(
      slope = est, SE = se,
      df = if (!is.na(df_col)) unname(row[1, df_col]) else NA_real_,
      t = if (!is.na(t_col)) unname(row[1, t_col]) else NA_real_,
      P = if (!is.na(p_col)) unname(row[1, p_col]) else NA_real_,
      CI_lower = est - 1.96 * se, CI_upper = est + 1.96 * se,
      model_valid = TRUE,
      final_structure = if (final_meta$a$correlated) "CORRELATED_RANDOM_SLOPE" else "UNCORRELATED_RANDOM_SLOPE_FALLBACK",
      optimizer = final_meta$a$optimizer,
      singular = final_meta$singular, convergence_ok = TRUE, n_attempts = length(attempt_log),
      random_intercept_variance = int_var, random_slope_variance = slope_var,
      intercept_slope_correlation = corr, residual_variance = resid_var,
      max_abs_scaled_gradient = max_grad,
      convergence_messages = paste(attempt_log, collapse = " || "),
      stringsAsFactors = FALSE
    ),
    fit = final
  )
}

prepare_inputs <- function() {
  layers <- c("miRNA", "GENCODE")
  all_inputs <- list()
  dataset_summary <- list()
  equivalence <- list()
  limma_results <- list()
  sensitivity_results <- list()
  pilot_rows <- list()
  center_rows <- list()
  top_features <- list()

  old_qc <- read.csv(file.path(project_root, "revision_stage0c", "QC_SENSITIVITY_RESULTS.csv"), check.names = FALSE)

  for (layer in layers) {
    full <- read_arm(layer, "FULL")
    qc <- read_arm(layer, "QC_PASS_ONLY")
    if (!identical(rownames(full$gradient), rownames(qc$gradient))) stop("Feature mismatch across arms: ", layer)
    features <- rownames(full$gradient)
    dominant <- if (layer == "miRNA") {
      match("hsa-miR-486-5p", features)
    } else {
      z <- grep("^RNY4(?=:|$)", features, perl = TRUE)
      if (length(z) != 1L) stop("RNY4 mapping is not unique: ", paste(features[z], collapse = ";"))
      z
    }
    if (is.na(dominant)) stop("Dominant reference feature absent: ", layer)
    abundance <- rowMeans(cbind(full$CSF_RPM, full$plasma_RPM))
    top5 <- features[order(abundance, decreasing = TRUE)[seq_len(5L)]]
    top_features[[layer]] <- data.frame(layer = layer, rank = 1:5, feature = top5, mean_RPM = abundance[top5])
    ref_all <- features
    ref_excl <- features[-dominant]
    ref_top5 <- setdiff(features, top5)

    arms <- list(FULL = full, QC_PASS_ONLY = qc)
    for (arm_name in names(arms)) {
      a <- arms[[arm_name]]
      raw <- log2(a$CSF_RPM + 1) - log2(a$plasma_RPM + 1)
      max_raw_diff <- max(abs(raw - a$gradient))
      clr <- make_clr_gradient(a$CSF_RPM, a$plasma_RPM, ref_all, 1)
      clr_formula <- make_reference_centered_from_raw(raw, ref_all)
      eq_all <- max(abs(clr - clr_formula))
      clr_excl <- make_clr_gradient(a$CSF_RPM, a$plasma_RPM, ref_excl, 1)
      excl_formula <- make_reference_centered_from_raw(raw, ref_excl)
      eq_excl <- max(abs(clr_excl - excl_formula))
      clr_pc05 <- make_clr_gradient(a$CSF_RPM, a$plasma_RPM, ref_all, 0.5)
      clr_top5 <- make_clr_gradient(a$CSF_RPM, a$plasma_RPM, ref_top5, 1)

      equivalence[[paste(layer, arm_name, sep = "_")]] <- data.frame(
        layer = layer, dataset = arm_name,
        raw_vs_saved_max_abs_diff = max_raw_diff,
        clr_identity_max_abs_diff = eq_all,
        excluded_reference_identity_max_abs_diff = eq_excl,
        excluded_reference_feature = features[dominant],
        stringsAsFactors = FALSE
      )

      m1 <- fit_limma_linear(raw, a$pairs)
      m2 <- fit_limma_linear(clr, a$pairs)
      s_excl <- fit_limma_linear(clr_excl, a$pairs)
      s_pc05 <- fit_limma_linear(clr_pc05, a$pairs)
      s_top5 <- fit_limma_linear(clr_top5, a$pairs)
      m1$layer <- m2$layer <- layer
      m1$dataset <- m2$dataset <- arm_name
      m1$model <- "M1"
      m2$model <- "M2"
      limma_results[[paste(layer, arm_name, "M1", sep = "_")]] <- m1
      limma_results[[paste(layer, arm_name, "M2", sep = "_")]] <- m2
      sensitivity_results[[paste(layer, arm_name, sep = "_")]] <- data.frame(
        feature = features, layer = layer, dataset = arm_name,
        excluded_reference_feature = features[dominant],
        excluded_ref_slope = s_excl$slope, excluded_ref_SE = s_excl$SE,
        excluded_ref_P = s_excl$P, excluded_ref_FDR = s_excl$FDR,
        pc0_5_slope = s_pc05$slope, pc0_5_SE = s_pc05$SE,
        pc0_5_P = s_pc05$P, pc0_5_FDR = s_pc05$FDR,
        top5_excluded_slope = s_top5$slope, top5_excluded_SE = s_top5$SE,
        top5_excluded_P = s_top5$P, top5_excluded_FDR = s_top5$FDR,
        stringsAsFactors = FALSE
      )

      center_raw <- colMeans(raw)
      center_excl <- colMeans(raw[match(ref_excl, rownames(raw)), , drop = FALSE])
      center_top5 <- colMeans(raw[match(ref_top5, rownames(raw)), , drop = FALSE])
      center_fit <- function(y) {
        d <- data.frame(y = y, centered_label = a$pairs$centered_label, subject = a$pairs$subject)
        des <- model.matrix(~ centered_label, d)
        dc <- duplicateCorrelation(matrix(y, nrow = 1L), des, block = d$subject)
        ft <- eBayes(lmFit(matrix(y, nrow = 1L), des, block = d$subject, correlation = dc$consensus.correlation), trend = TRUE, robust = TRUE)
        se <- ft$stdev.unscaled[1, 2] * ft$sigma[1]
        c(slope = unname(ft$coefficients[1, 2]), SE = unname(se), P = unname(ft$p.value[1, 2]))
      }
      for (nm in c("all_reference", "excluded_dominant", "excluded_top5")) {
        yy <- switch(nm, all_reference = center_raw, excluded_dominant = center_excl, excluded_top5 = center_top5)
        cf <- center_fit(yy)
        center_rows[[length(center_rows) + 1L]] <- data.frame(layer = layer, dataset = arm_name, reference = nm,
          slope = cf[["slope"]], SE = cf[["SE"]], P = cf[["P"]])
      }

      csf_zero <- rowMeans(a$CSF_RPM == 0)
      plasma_zero <- rowMeans(a$plasma_RPM == 0)
      combined_zero <- rowMeans(cbind(a$CSF_RPM, a$plasma_RPM) == 0)
      combined_low <- rowMeans(cbind(a$CSF_RPM, a$plasma_RPM) <= 1)
      labels <- c("ALL", sort(unique(a$pairs$collection_label)))
      for (lab in labels) {
        idx <- if (lab == "ALL") seq_len(nrow(a$pairs)) else which(a$pairs$collection_label == lab)
        dataset_summary[[length(dataset_summary) + 1L]] <- data.frame(
          dataset = arm_name, layer = layer, collection_label = lab,
          n_subjects = length(unique(a$pairs$subject[idx])), n_pairs = length(idx), n_biosamples = 2L * length(idx),
          n_features = length(features),
          overall_CSF_zero_fraction = mean(a$CSF_RPM[, idx, drop = FALSE] == 0),
          overall_plasma_zero_fraction = mean(a$plasma_RPM[, idx, drop = FALSE] == 0),
          median_feature_zero_fraction_both_fluids = median(combined_zero),
          features_over_50pct_zero = sum(combined_zero > 0.5),
          features_all_zero = sum(combined_zero == 1),
          overall_RPM_le1_fraction = mean(cbind(a$CSF_RPM[, idx, drop = FALSE], a$plasma_RPM[, idx, drop = FALSE]) <= 1),
          stringsAsFactors = FALSE
        )
      }

      all_inputs[[paste(layer, arm_name, sep = "_")]] <- list(
        layer = layer, dataset = arm_name, pairs = a$pairs, raw = raw, clr = clr,
        csf = a$CSF_RPM, plasma = a$plasma_RPM,
        feature_zero_fraction = combined_zero, feature_low_le1_fraction = combined_low,
        reference_all = ref_all, reference_excluded = ref_excl, reference_top5 = ref_top5,
        dominant_feature = features[dominant], top5 = top5
      )
    }

    full_m1 <- limma_results[[paste(layer, "FULL", "M1", sep = "_")]]
    raw_full <- all_inputs[[paste(layer, "FULL", sep = "_")]]$raw
    pairs_full <- all_inputs[[paste(layer, "FULL", sep = "_")]]$pairs
    sds <- subject_slope_sd(raw_full, pairs_full)
    zfrac <- all_inputs[[paste(layer, "FULL", sep = "_")]]$feature_zero_fraction
    add_pick <- function(idx, reason) {
      if (length(idx)) pilot_rows[[length(pilot_rows) + 1L]] <<- data.frame(layer = layer, feature = features[idx], selection_reason = reason)
    }
    sig_neg <- which(full_m1$FDR < .05 & full_m1$slope < 0)
    if (length(sig_neg)) add_pick(sig_neg[which.min(full_m1$P[sig_neg])], "strongest_old_significant_negative")
    sig_pos <- which(full_m1$FDR < .05 & full_m1$slope > 0)
    if (length(sig_pos)) {
      take <- sig_pos[order(full_m1$P[sig_pos])][seq_len(min(length(sig_pos), if (layer == "miRNA") 3L else 2L))]
      add_pick(take, "old_significant_positive")
    }
    sig <- which(full_m1$FDR < .05)
    ns <- which(full_m1$FDR >= .05)
    if (length(sig)) add_pick(sig[which.max(full_m1$FDR[sig])], "old_q_nearest_0.05_significant")
    if (length(ns)) add_pick(ns[which.min(full_m1$FDR[ns])], "old_q_nearest_0.05_nonsignificant")
    add_pick(which.max(zfrac), "highest_zero_fraction")
    add_pick(which.max(sds), "largest_subject_slope_SD")
    reps <- if (layer == "miRNA") {
      c("hsa-miR-338-3p", "hsa-miR-381-3p", "hsa-miR-410-3p", "hsa-miR-375")
    } else {
      c("MT-TL2:Mt_tRNA", "RN7SL1:antisense", "IL8:protein_coding", "HBA1:protein_coding",
        "TMSB4X:protein_coding", "EIF4A1:retained_intron")
    }
    for (r in reps) add_pick(match(r, features), "Figure_3_or_4_representative")
  }

  pilot_df <- do.call(rbind, pilot_rows)
  pilot_df <- aggregate(selection_reason ~ layer + feature, pilot_df, function(x) paste(unique(x), collapse = ";"))
  pilot_df <- pilot_df[order(pilot_df$layer, pilot_df$feature), ]
  names(all_inputs) <- names(all_inputs)
  saveRDS(list(
    inputs = all_inputs,
    limma = limma_results,
    sensitivities = sensitivity_results,
    equivalence = do.call(rbind, equivalence),
    dataset_summary = do.call(rbind, dataset_summary),
    pilot_features = pilot_df,
    center_slopes = do.call(rbind, center_rows),
    top_features = do.call(rbind, top_features),
    old_qc = old_qc
  ), file.path(out_dir, "intermediate", "prepared_stage1_inputs.rds"))
  write_csv(do.call(rbind, dataset_summary), file.path(out_dir, "STAGE1_DATASET_SUMMARY.csv"))
  write_csv(do.call(rbind, equivalence), file.path(out_dir, "intermediate", "CLR_EQUIVALENCE_CHECK.csv"))
  write_csv(pilot_df, file.path(out_dir, "intermediate", "PILOT_FEATURE_SELECTION.csv"))
  write_csv(do.call(rbind, center_rows), file.path(out_dir, "intermediate", "CLR_CENTER_SLOPES.csv"))
  write_csv(do.call(rbind, top_features), file.path(out_dir, "intermediate", "TOP5_HIGH_ABUNDANCE_FEATURES.csv"))
  invisible(TRUE)
}

plot_pilot_fit <- function(fit, y, pairs, feature, model, layer) {
  if (is.null(fit)) {
    plot.new(); title(main = paste(layer, model, feature)); text(.5, .5, "MODEL NOT EVALUABLE")
    return()
  }
  dat <- data.frame(y = y, time = pairs$centered_label, subject = factor(pairs$subject))
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(dat$time, dat$y, type = "n", xlab = "Collection index centered at D04", ylab = "Gradient", main = paste(layer, model, feature))
  for (s in levels(dat$subject)) {
    ii <- dat$subject == s
    lines(dat$time[ii], dat$y[ii], col = grDevices::adjustcolor("grey40", .35))
  }
  points(dat$time, dat$y, pch = 16, cex = .55)
  abline(fixef(fit)[1], fixef(fit)[2], col = "red", lwd = 2)
  plot(fitted(fit), resid(fit), pch = 16, cex = .65, xlab = "Fitted", ylab = "Residual", main = "Residual vs fitted")
  abline(h = 0, lty = 2)
  qqnorm(resid(fit), pch = 16, cex = .65, main = "Residual Q-Q")
  qqline(resid(fit), col = "red")
  re <- ranef(fit)$subject
  slope_col <- grep("label_centered", colnames(re), value = TRUE)[1]
  if (!is.na(slope_col)) {
    hist(re[[slope_col]], breaks = "FD", main = "Random slope BLUP", xlab = "Subject deviation")
  } else {
    plot.new(); text(.5, .5, "Random slope not extracted")
  }
}

run_pilot <- function() {
  prepare_inputs()
  obj <- readRDS(file.path(out_dir, "intermediate", "prepared_stage1_inputs.rds"))
  pilots <- obj$pilot_features
  results <- list()
  fit_store <- list()
  for (layer in c("miRNA", "GENCODE")) {
    input <- obj$inputs[[paste(layer, "FULL", sep = "_")]]
    selected <- pilots[pilots$layer == layer, , drop = FALSE]
    for (model in c("M3", "M4")) {
      mat <- if (model == "M3") input$raw else input$clr
      pdf_path <- file.path(out_dir, "diagnostic_plots", paste0("PILOT_", layer, "_", model, "_diagnostics.pdf"))
      pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
      for (i in seq_len(nrow(selected))) {
        feature <- selected$feature[i]
        idx <- match(feature, rownames(mat))
        fit <- fit_feature_lmm(mat[idx, ], input$pairs)
        rr <- cbind(data.frame(layer = layer, dataset = "FULL", model = model, feature = feature,
          selection_reason = selected$selection_reason[i], stringsAsFactors = FALSE), fit$result)
        results[[length(results) + 1L]] <- rr
        plot_pilot_fit(fit$fit, mat[idx, ], input$pairs, feature, model, layer)
        fit_store[[paste(layer, model, feature, sep = "::")]] <- fit$fit
      }
      dev.off()
    }
  }
  res <- do.call(rbind, results)
  res$FDR_pilot_descriptive_only <- ave(res$P, res$layer, res$model, FUN = bh_fixed)
  write_csv(res, file.path(out_dir, "intermediate", "PILOT_MODEL_DIAGNOSTICS.csv"))
  saveRDS(fit_store, file.path(out_dir, "intermediate", "PILOT_FITS.rds"))
  gate <- aggregate(model_valid ~ layer + model, res, function(x) c(n = length(x), valid = sum(x), invalid = sum(!x), valid_rate = mean(x)))
  gate <- data.frame(layer = gate$layer, model = gate$model,
    n = gate$model_valid[, "n"], valid = gate$model_valid[, "valid"], invalid = gate$model_valid[, "invalid"],
    valid_rate = gate$model_valid[, "valid_rate"])
  error_text <- res$convergence_messages[grepl(":ERROR:", res$convergence_messages, fixed = TRUE)]
  common_errors <- table(error_text[nzchar(error_text)])
  env_failure <- length(common_errors) && max(common_errors) >= ceiling(.25 * nrow(res))
  gate_pass <- all(gate$valid_rate >= .5) && !env_failure
  gate$environment_error_trigger <- env_failure
  gate$gate <- if (gate_pass) "PASS_TO_FULL" else "STOP_SYSTEMATIC_FAILURE"
  write_csv(gate, file.path(out_dir, "intermediate", "PILOT_GATE.csv"))
  cat("PILOT_GATE=", if (gate_pass) "PASS_TO_FULL" else "STOP_SYSTEMATIC_FAILURE", "\n", sep = "")
  print(gate)
  if (!gate_pass) stop("Pilot gate failed under frozen rules; full mixed-model batch not started.")
  invisible(TRUE)
}

fit_matrix_lmm <- function(mat, pairs, layer, dataset, model) {
  out <- vector("list", nrow(mat))
  for (i in seq_len(nrow(mat))) {
    ff <- fit_feature_lmm(mat[i, ], pairs)
    out[[i]] <- cbind(data.frame(feature = rownames(mat)[i], layer = layer, dataset = dataset, model = model,
      stringsAsFactors = FALSE), ff$result)
    if (i %% 50L == 0L) cat(format(Sys.time()), layer, dataset, model, i, "/", nrow(mat), "\n")
  }
  ans <- do.call(rbind, out)
  ans$FDR <- bh_fixed(ans$P)
  ans$direction <- direction(ans$slope)
  ans
}

run_full <- function() {
  gate_file <- file.path(out_dir, "intermediate", "PILOT_GATE.csv")
  if (!file.exists(gate_file)) stop("Pilot gate file missing")
  gate <- read.csv(gate_file, check.names = FALSE)
  if (!all(gate$gate == "PASS_TO_FULL")) stop("Pilot gate did not pass")
  obj <- readRDS(file.path(out_dir, "intermediate", "prepared_stage1_inputs.rds"))
  tasks <- expand.grid(layer = c("miRNA", "GENCODE"), dataset = c("FULL", "QC_PASS_ONLY"), model = c("M3", "M4"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(tasks))) {
    tt <- tasks[i, ]
    key <- paste(tt$layer, tt$dataset, sep = "_")
    input <- obj$inputs[[key]]
    mat <- if (tt$model == "M3") input$raw else input$clr
    ans <- fit_matrix_lmm(mat, input$pairs, tt$layer, tt$dataset, tt$model)
    path <- file.path(out_dir, "intermediate", paste0("LMM_", tt$layer, "_", tt$dataset, "_", tt$model, ".csv"))
    write_csv(ans, path)
    cat("WROTE", path, "\n")
  }
}

prefix_cols <- function(d, prefix, keep = c("feature", "layer")) {
  names(d)[!names(d) %in% keep] <- paste0(prefix, names(d)[!names(d) %in% keep])
  d
}

safe_merge <- function(x, y, by = c("feature", "layer")) merge(x, y, by = by, all = TRUE, sort = FALSE)

run_finalize <- function() {
  obj <- readRDS(file.path(out_dir, "intermediate", "prepared_stage1_inputs.rds"))
  all_lmm <- list()
  for (layer in c("miRNA", "GENCODE")) for (dataset in c("FULL", "QC_PASS_ONLY")) for (model in c("M3", "M4")) {
    p <- file.path(out_dir, "intermediate", paste0("LMM_", layer, "_", dataset, "_", model, ".csv"))
    if (!file.exists(p)) stop("Missing LMM output: ", p)
    all_lmm[[paste(layer, dataset, model, sep = "_")]] <- read.csv(p, check.names = FALSE)
  }

  diagnostics <- do.call(rbind, all_lmm)
  write_csv(diagnostics, file.path(out_dir, "STAGE1_MODEL_DIAGNOSTICS.csv"))

  all_combined <- list()
  qc_rows <- list()
  robust_rows <- list()
  rep_rows <- list()
  direction_rows <- list()
  old_set_rows <- list()

  for (layer in c("miRNA", "GENCODE")) {
    full_m1 <- obj$limma[[paste(layer, "FULL", "M1", sep = "_")]]
    full_m2 <- obj$limma[[paste(layer, "FULL", "M2", sep = "_")]]
    qc_m1 <- obj$limma[[paste(layer, "QC_PASS_ONLY", "M1", sep = "_")]]
    qc_m2 <- obj$limma[[paste(layer, "QC_PASS_ONLY", "M2", sep = "_")]]
    full_m3 <- all_lmm[[paste(layer, "FULL", "M3", sep = "_")]]
    full_m4 <- all_lmm[[paste(layer, "FULL", "M4", sep = "_")]]
    qc_m3 <- all_lmm[[paste(layer, "QC_PASS_ONLY", "M3", sep = "_")]]
    qc_m4 <- all_lmm[[paste(layer, "QC_PASS_ONLY", "M4", sep = "_")]]
    sens_full <- obj$sensitivities[[paste(layer, "FULL", sep = "_")]]
    sens_qc <- obj$sensitivities[[paste(layer, "QC_PASS_ONLY", sep = "_")]]

    base <- data.frame(feature = full_m1$feature, layer = layer, stringsAsFactors = FALSE)
    add_limma <- function(base, d, pre) {
      z <- d[, c("feature", "slope", "SE", "CI_lower", "CI_upper", "P", "FDR", "duplicate_correlation")]
      names(z)[-1] <- paste0(pre, names(z)[-1])
      merge(base, z, by = "feature", all.x = TRUE, sort = FALSE)
    }
    add_lmm <- function(base, d, pre) {
      keep <- setdiff(names(d), c("layer", "dataset", "model"))
      z <- d[, keep]
      names(z)[names(z) != "feature"] <- paste0(pre, names(z)[names(z) != "feature"])
      merge(base, z, by = "feature", all.x = TRUE, sort = FALSE)
    }
    base <- add_limma(base, full_m1, "FULL_M1_")
    base <- add_limma(base, full_m2, "FULL_M2_")
    base <- add_lmm(base, full_m3, "FULL_M3_")
    base <- add_lmm(base, full_m4, "FULL_M4_")
    base <- add_limma(base, qc_m1, "QC_M1_")
    base <- add_limma(base, qc_m2, "QC_M2_")
    base <- add_lmm(base, qc_m3, "QC_M3_")
    base <- add_lmm(base, qc_m4, "QC_M4_")
    sf <- sens_full
    names(sf)[!names(sf) %in% c("feature", "layer", "dataset")] <- paste0("FULL_SENS_", names(sf)[!names(sf) %in% c("feature", "layer", "dataset")])
    sf$dataset <- NULL; sf$layer <- NULL
    sq <- sens_qc
    names(sq)[!names(sq) %in% c("feature", "layer", "dataset")] <- paste0("QC_SENS_", names(sq)[!names(sq) %in% c("feature", "layer", "dataset")])
    sq$dataset <- NULL; sq$layer <- NULL
    base <- merge(base, sf, by = "feature", all.x = TRUE, sort = FALSE)
    base <- merge(base, sq, by = "feature", all.x = TRUE, sort = FALSE)

    old_sig <- base$FULL_M1_FDR < .05
    same_full_direction <- direction(base$FULL_M1_slope) == direction(base$FULL_M3_slope) & direction(base$FULL_M1_slope) == direction(base$FULL_M4_slope)
    qc_direction_ok <- direction(base$QC_M3_slope) == direction(base$FULL_M1_slope) & direction(base$QC_M4_slope) == direction(base$FULL_M1_slope)
    qc_delta_m3 <- abs(base$QC_M3_slope - base$FULL_M3_slope)
    qc_delta_m4 <- abs(base$QC_M4_slope - base$FULL_M4_slope)
    qc_threshold_m3 <- pmax(.5 * abs(base$FULL_M3_slope), .1)
    qc_threshold_m4 <- pmax(.5 * abs(base$FULL_M4_slope), .1)
    qc_dependency <- !qc_direction_ok | qc_delta_m3 > qc_threshold_m3 | qc_delta_m4 > qc_threshold_m4
    model_evaluable <- base$FULL_M3_model_valid & base$FULL_M4_model_valid & base$QC_M3_model_valid & base$QC_M4_model_valid
    core <- old_sig & base$FULL_M3_FDR < .05 & base$FULL_M4_FDR < .05 & same_full_direction & model_evaluable & !qc_dependency
    comp_sensitive <- old_sig & base$FULL_M4_model_valid & (base$FULL_M4_FDR >= .05 | direction(base$FULL_M4_slope) != direction(base$FULL_M1_slope))
    random_sensitive <- old_sig & base$FULL_M3_model_valid & (base$FULL_M3_FDR >= .05 | direction(base$FULL_M3_slope) != direction(base$FULL_M1_slope))
    qc_sensitive <- old_sig & base$FULL_M3_FDR < .05 & base$FULL_M4_FDR < .05 & same_full_direction & model_evaluable & qc_dependency
    category <- ifelse(!model_evaluable, "MODEL_NOT_EVALUABLE",
      ifelse(core, "CORE_ROBUST",
        ifelse(comp_sensitive, "COMPOSITION_SENSITIVE",
          ifelse(random_sensitive, "RANDOM_SLOPE_SENSITIVE",
            ifelse(qc_sensitive, "QC_SENSITIVE", "OTHER_NOT_CORE")))))
    base$original_M1_significant <- old_sig
    base$full_direction_consistent_M1_M3_M4 <- same_full_direction
    base$qc_direction_consistent <- qc_direction_ok
    base$qc_delta_M3 <- qc_delta_m3
    base$qc_delta_M4 <- qc_delta_m4
    base$qc_dependency_threshold_M3 <- qc_threshold_m3
    base$qc_dependency_threshold_M4 <- qc_threshold_m4
    base$qc_dependency_flag <- qc_dependency
    base$all_random_slope_models_evaluable <- model_evaluable
    base$core_robust <- core
    base$composition_sensitive <- comp_sensitive
    base$random_slope_sensitive <- random_sensitive
    base$qc_sensitive <- qc_sensitive
    base$robustness_category <- category
    input <- obj$inputs[[paste(layer, "FULL", sep = "_")]]
    base$full_feature_zero_fraction <- input$feature_zero_fraction[match(base$feature, names(input$feature_zero_fraction))]
    base$full_feature_RPM_le1_fraction <- input$feature_low_le1_fraction[match(base$feature, names(input$feature_low_le1_fraction))]
    all_combined[[layer]] <- base

    qc_rows[[layer]] <- base[, c("feature", "layer", grep("^FULL_M[1-4]_(slope|SE|FDR|model_valid)$|^QC_M[1-4]_(slope|SE|FDR|model_valid)$|^qc_", names(base), value = TRUE)), drop = FALSE]
    robust_rows[[layer]] <- base[, c("feature", "layer", "robustness_category", "core_robust", "composition_sensitive", "random_slope_sensitive", "qc_sensitive", "all_random_slope_models_evaluable",
      "FULL_M1_slope", "FULL_M1_FDR", "FULL_M3_slope", "FULL_M3_FDR", "FULL_M4_slope", "FULL_M4_FDR",
      "QC_M3_slope", "QC_M3_FDR", "QC_M4_slope", "QC_M4_FDR", "qc_dependency_flag")]

    # Direction summaries for all features and significant features.
    for (dataset in c("FULL", "QC")) for (model in paste0("M", 1:4)) {
      sc <- paste0(dataset, "_", model, "_slope")
      fc <- paste0(dataset, "_", model, "_FDR")
      vc <- if (model %in% c("M3", "M4")) paste0(dataset, "_", model, "_model_valid") else NULL
      eval_idx <- is.finite(base[[sc]])
      if (!is.null(vc)) eval_idx <- eval_idx & base[[vc]]
      sig_idx <- eval_idx & base[[fc]] < .05
      for (set_name in c("all_evaluable", "FDR_lt_0.05")) {
        ii <- if (set_name == "all_evaluable") eval_idx else sig_idx
        dd <- direction(base[[sc]][ii])
        direction_rows[[length(direction_rows) + 1L]] <- data.frame(
          layer = layer, dataset = if (dataset == "QC") "QC_PASS_ONLY" else "FULL", model = model, feature_set = set_name,
          n_evaluable = sum(ii), n_negative = sum(dd == "negative", na.rm = TRUE), n_positive = sum(dd == "positive", na.rm = TRUE),
          n_zero = sum(dd == "zero", na.rm = TRUE), slope_median = if (sum(ii)) median(base[[sc]][ii]) else NA_real_,
          slope_Q1 = if (sum(ii)) unname(quantile(base[[sc]][ii], .25)) else NA_real_,
          slope_Q3 = if (sum(ii)) unname(quantile(base[[sc]][ii], .75)) else NA_real_, stringsAsFactors = FALSE)
      }
    }

    old_features <- base$feature[old_sig]
    for (model in paste0("M", 1:4)) {
      sc <- paste0("FULL_", model, "_slope")
      fc <- paste0("FULL_", model, "_FDR")
      ii <- match(old_features, base$feature)
      dd <- direction(base[[sc]][ii])
      old_set_rows[[length(old_set_rows) + 1L]] <- data.frame(layer = layer, model = model,
        original_significant_set_n = length(ii), n_negative = sum(dd == "negative", na.rm = TRUE),
        n_positive = sum(dd == "positive", na.rm = TRUE), n_direction_changed_vs_M1 = sum(dd != direction(base$FULL_M1_slope[ii]), na.rm = TRUE),
        n_still_FDR_lt_0_05 = sum(base[[fc]][ii] < .05, na.rm = TRUE), stringsAsFactors = FALSE)
    }

    rep_targets <- if (layer == "miRNA") c("hsa-miR-338-3p", "hsa-miR-381-3p", "hsa-miR-410-3p", "hsa-miR-375") else {
      unlist(lapply(c("MT-TL2", "RN7SL1", "IL8", "HBA1", "TMSB4X", "EIF4A1"), function(p) grep(paste0("^", p, "(?=:|$)"), base$feature, value = TRUE, perl = TRUE)), use.names = FALSE)
    }
    rr <- base[match(rep_targets, base$feature), ]
    rr$representative_status <- ifelse(rr$all_random_slope_models_evaluable & rr$full_direction_consistent_M1_M3_M4 & !rr$qc_dependency_flag,
      ifelse(rr$FULL_M1_FDR < .05 & rr$FULL_M3_FDR < .05 & rr$FULL_M4_FDR < .05, "RETAIN_STRONG", "RETAIN_WITH_QUALIFICATION"),
      ifelse(!rr$all_random_slope_models_evaluable, "NOT_EVALUABLE_RANDOM_SLOPE", "QC_OR_DIRECTION_SENSITIVE"))
    rep_rows[[layer]] <- rr
  }

  write_csv(all_combined$miRNA, file.path(out_dir, "STAGE1_MODEL_RESULTS_miRNA.csv"))
  write_csv(all_combined$GENCODE, file.path(out_dir, "STAGE1_MODEL_RESULTS_GENCODE.csv"))
  write_csv(do.call(rbind, qc_rows), file.path(out_dir, "STAGE1_QC_SENSITIVITY.csv"))
  write_csv(do.call(rbind, robust_rows), file.path(out_dir, "STAGE1_ROBUST_SIGNATURE.csv"))
  direction_out <- rbind(do.call(rbind, direction_rows), transform(do.call(rbind, old_set_rows), dataset = "FULL", feature_set = "ORIGINAL_M1_SIGNIFICANT", n_evaluable = original_significant_set_n, n_zero = NA, slope_median = NA, slope_Q1 = NA, slope_Q3 = NA)[, names(do.call(rbind, direction_rows))])
  write_csv(direction_out, file.path(out_dir, "STAGE1_SLOPE_DIRECTION_SUMMARY.csv"))

  reps <- do.call(rbind, rep_rows)
  keep_rep <- c("feature", "layer", "representative_status", "robustness_category", "core_robust", "qc_dependency_flag",
    grep("^(FULL|QC)_M[1-4]_(slope|SE|CI_lower|CI_upper|P|FDR|model_valid|final_structure|singular|convergence_ok)$", names(reps), value = TRUE))
  write_csv(reps[, unique(keep_rep), drop = FALSE], file.path(out_dir, "STAGE1_REPRESENTATIVE_FEATURES.csv"))
  write_csv(obj$equivalence, file.path(out_dir, "intermediate", "CLR_EQUIVALENCE_CHECK.csv"))
  write_csv(do.call(rbind, obj$sensitivities), file.path(out_dir, "intermediate", "REFERENCE_SENSITIVITY_RESULTS.csv"))

  # Compact summary tables used to construct the report.
  summary_models <- do.call(rbind, lapply(c("miRNA", "GENCODE"), function(layer) {
    b <- all_combined[[layer]]
    do.call(rbind, lapply(c("FULL", "QC"), function(ds) do.call(rbind, lapply(paste0("M", 1:4), function(m) {
      fdr <- b[[paste0(ds, "_", m, "_FDR")]]
      slope <- b[[paste0(ds, "_", m, "_slope")]]
      valid <- if (m %in% c("M3", "M4")) b[[paste0(ds, "_", m, "_model_valid")]] else rep(TRUE, nrow(b))
      data.frame(layer = layer, dataset = if (ds == "QC") "QC_PASS_ONLY" else "FULL", model = m,
        n_features = nrow(b), n_model_valid = sum(valid, na.rm = TRUE), n_significant = sum(fdr < .05 & valid, na.rm = TRUE),
        significant_negative = sum(fdr < .05 & valid & slope < 0, na.rm = TRUE),
        significant_positive = sum(fdr < .05 & valid & slope > 0, na.rm = TRUE), stringsAsFactors = FALSE)
    }))))
  }))
  write_csv(summary_models, file.path(out_dir, "intermediate", "MODEL_SUMMARY_COUNTS.csv"))
  robust_summary <- do.call(rbind, lapply(all_combined, function(b) as.data.frame(table(layer = b$layer, category = b$robustness_category), stringsAsFactors = FALSE)))
  write_csv(robust_summary, file.path(out_dir, "intermediate", "ROBUST_CATEGORY_COUNTS.csv"))

  # Comparison plots (analysis figures only, not manuscript figures).
  png(file.path(out_dir, "diagnostic_plots", "NEW_VS_OLD_SLOPES.png"), width = 2400, height = 1800, res = 200)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (layer in c("miRNA", "GENCODE")) {
    b <- all_combined[[layer]]
    plot(b$FULL_M1_slope, b$FULL_M3_slope, pch = 16, cex = .45, xlab = "M1 slope", ylab = "M3 slope", main = paste(layer, "raw: old vs random slope"))
    abline(0, 1, col = "red")
    plot(b$FULL_M2_slope, b$FULL_M4_slope, pch = 16, cex = .45, xlab = "M2 slope", ylab = "M4 slope", main = paste(layer, "CLR: limma vs random slope"))
    abline(0, 1, col = "red")
  }
  dev.off()
  png(file.path(out_dir, "diagnostic_plots", "FULL_VS_QC_SLOPES.png"), width = 2400, height = 1800, res = 200)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (layer in c("miRNA", "GENCODE")) for (m in c("M3", "M4")) {
    b <- all_combined[[layer]]
    plot(b[[paste0("FULL_", m, "_slope")]], b[[paste0("QC_", m, "_slope")]], pch = 16, cex = .45,
      xlab = paste("FULL", m), ylab = paste("QC PASS-only", m), main = paste(layer, m, "QC sensitivity"))
    abline(0, 1, col = "red")
  }
  dev.off()

  capture.output(sessionInfo(), file = file.path(out_dir, "logs", "STAGE1_sessionInfo.txt"))
  cat("FINALIZE_COMPLETE\n")
}

if (phase == "pilot") run_pilot()
if (phase == "full") run_full()
if (phase == "finalize") run_finalize()

