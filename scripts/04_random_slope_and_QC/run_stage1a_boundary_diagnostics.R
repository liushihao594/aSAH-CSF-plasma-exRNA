options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(clubSandwich)
})

root <- "."
stage1 <- file.path(root, "revision_stage1")
out <- file.path(stage1, "stage1a_boundary_diagnostics")
dir.create(file.path(out, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "intermediate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "plots"), recursive = TRUE, showWarnings = FALSE)

tol_singular <- 1e-4
seed <- 20261121L
set.seed(seed)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
direction <- function(x, tol = sqrt(.Machine$double.eps)) {
  ifelse(is.na(x), NA_character_, ifelse(x > tol, "positive", ifelse(x < -tol, "negative", "zero")))
}
safe_ratio <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b != 0, a / b, NA_real_)

obj <- readRDS(file.path(stage1, "intermediate", "prepared_stage1_inputs.rds"))
old_selection <- read.csv(file.path(stage1, "intermediate", "PILOT_FEATURE_SELECTION.csv"), check.names = FALSE)
old_diag <- read.csv(file.path(stage1, "STAGE1_MODEL_DIAGNOSTICS.csv"), check.names = FALSE)

# Expanded pilot selection is deterministic from Stage 1 M1/M2 and expression structure.
selection_rows <- list()
for (layer in c("miRNA", "GENCODE")) {
  old <- old_selection[old_selection$layer == layer, c("layer", "feature", "selection_reason")]
  selection_rows[[length(selection_rows) + 1L]] <- old
  m1 <- obj$limma[[paste(layer, "FULL", "M1", sep = "_")]]
  m2 <- obj$limma[[paste(layer, "FULL", "M2", sep = "_")]]
  inp <- obj$inputs[[paste(layer, "FULL", sep = "_")]]
  existing <- old$feature

  discord <- which(m1$FDR < .05 & direction(m1$slope) != direction(m2$slope))
  if (length(discord)) {
    discord <- discord[order(abs(m1$slope[discord]), decreasing = TRUE)]
    add <- setdiff(m1$feature[discord], existing)
    add <- head(add, 3L)
    if (length(add)) {
      selection_rows[[length(selection_rows) + 1L]] <- data.frame(layer = layer, feature = add,
        selection_reason = "expanded_M1_significant_M1_M2_direction_discordant_top_abs_M1", stringsAsFactors = FALSE)
      existing <- c(existing, add)
    }
  }

  z <- inp$feature_zero_fraction
  add <- setdiff(names(sort(z, decreasing = TRUE)), existing)
  add <- head(add, 1L)
  if (length(add)) {
    selection_rows[[length(selection_rows) + 1L]] <- data.frame(layer = layer, feature = add,
      selection_reason = "expanded_next_highest_zero_fraction", stringsAsFactors = FALSE)
    existing <- c(existing, add)
  }

  clr_sd <- apply(inp$clr, 1L, sd)
  add <- setdiff(names(sort(clr_sd, decreasing = FALSE)), existing)
  add <- head(add[is.finite(clr_sd[add])], 1L)
  if (length(add)) {
    selection_rows[[length(selection_rows) + 1L]] <- data.frame(layer = layer, feature = add,
      selection_reason = "expanded_lowest_CLR_gradient_SD", stringsAsFactors = FALSE)
  }
}
selection <- do.call(rbind, selection_rows)
selection <- aggregate(selection_reason ~ layer + feature, selection,
  function(x) paste(unique(x), collapse = ";"))
selection <- selection[order(selection$layer, selection$feature), ]
write_csv(selection, file.path(out, "intermediate", "EXPANDED_PILOT_SELECTION.csv"))

model_formulas <- list(
  A = y ~ centered_label + (1 + centered_label | subject),
  B = y ~ centered_label + (1 + centered_label || subject),
  C = y ~ centered_label + (0 + centered_label | subject),
  D = y ~ centered_label + (1 | subject)
)
optimizers <- c("nloptwrap", "bobyqa", "nlminbwrap")

extract_opt_code <- function(fit) {
  z <- fit@optinfo$conv$opt
  if (is.null(z) || !length(z)) return("0")
  paste(z, collapse = ";")
}

fit_once <- function(dat, structure, optimizer) {
  warnings <- character()
  opt_ctrl <- if (optimizer == "nloptwrap") list(maxeval = 200000) else list(maxfun = 200000)
  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        model_formulas[[structure]], data = dat, REML = TRUE,
        control = lmerControl(
          optimizer = optimizer,
          optCtrl = opt_ctrl,
          check.conv.singular = lme4:::.makeCC(action = "ignore", tol = tol_singular)
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
    return(list(fit = NULL, row = data.frame(
      fit_error = conditionMessage(fit), warnings = paste(unique(warnings), collapse = " | "),
      optimizer_code = NA_character_, lme4_convergence_messages = NA_character_, convergence_ok = FALSE,
      fixed_effect_finite = FALSE, slope = NA_real_, SE = NA_real_, df = NA_real_, t = NA_real_, P = NA_real_,
      CI_lower = NA_real_, CI_upper = NA_real_, random_intercept_variance = NA_real_,
      random_slope_variance = NA_real_, intercept_slope_correlation = NA_real_, residual_variance = NA_real_,
      RI_SD_ratio = NA_real_, RS_SD_ratio_per_label = NA_real_, RS_SD_ratio_observed_span = NA_real_,
      covariance_eigen_min = NA_real_, covariance_eigen_max = NA_real_, covariance_eigen_ratio = NA_real_,
      rePCA_min_SD = NA_real_, rePCA_max_SD = NA_real_, rePCA_SD_ratio = NA_real_, rePCA_rank = NA_integer_,
      hessian_eigen_min = NA_real_, hessian_eigen_max = NA_real_, hessian_positive = FALSE,
      max_abs_scaled_gradient = NA_real_, isSingular_tol_1e_4 = NA, stringsAsFactors = FALSE)))
  }

  conv_msgs <- fit@optinfo$conv$lme4$messages
  conv_msgs <- if (length(conv_msgs)) paste(conv_msgs, collapse = " | ") else ""
  opt_code <- extract_opt_code(fit)
  opt_bad <- suppressWarnings(any(as.numeric(unlist(strsplit(opt_code, ";", fixed = TRUE))) != 0, na.rm = TRUE))
  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  fixed_ok <- !is.null(co) && "centered_label" %in% rownames(co) &&
    all(is.finite(co["centered_label", intersect(c("Estimate", "Std. Error"), colnames(co))]))
  est <- se <- ddf <- tv <- pv <- NA_real_
  if (fixed_ok) {
    est <- unname(co["centered_label", "Estimate"])
    se <- unname(co["centered_label", "Std. Error"])
    if ("df" %in% colnames(co)) ddf <- unname(co["centered_label", "df"])
    tc <- grep("t value", colnames(co), value = TRUE)[1]
    pc <- grep("Pr\\(", colnames(co), value = TRUE)[1]
    if (!is.na(tc)) tv <- unname(co["centered_label", tc])
    if (!is.na(pc)) pv <- unname(co["centered_label", pc])
  }

  vc <- as.data.frame(VarCorr(fit))
  get_v <- function(var) {
    z <- vc$vcov[vc$var1 == var & is.na(vc$var2)]
    if (length(z)) z[1] else NA_real_
  }
  ri_var <- get_v("(Intercept)")
  rs_var <- get_v("centered_label")
  resid_var <- vc$vcov[vc$grp == "Residual"][1]
  corr_idx <- !is.na(vc$var2) &
    ((vc$var1 == "(Intercept)" & vc$var2 == "centered_label") |
       (vc$var1 == "centered_label" & vc$var2 == "(Intercept)"))
  corr <- if (any(corr_idx)) vc$sdcor[corr_idx][1] else NA_real_
  resid_sd <- sqrt(resid_var)
  ri_ratio <- if (is.finite(ri_var)) sqrt(max(ri_var, 0)) / resid_sd else NA_real_
  rs_ratio <- if (is.finite(rs_var)) sqrt(max(rs_var, 0)) / resid_sd else NA_real_

  active_vars <- switch(structure,
    A = c(ri_var, rs_var), B = c(ri_var, rs_var), C = rs_var, D = ri_var)
  G <- NULL
  if (structure == "A" && all(is.finite(c(ri_var, rs_var, corr)))) {
    cv <- corr * sqrt(max(ri_var, 0) * max(rs_var, 0))
    G <- matrix(c(ri_var, cv, cv, rs_var), 2, 2)
  } else if (structure == "B" && all(is.finite(c(ri_var, rs_var)))) {
    G <- diag(c(ri_var, rs_var), 2)
  } else if (structure == "C" && is.finite(rs_var)) {
    G <- matrix(rs_var, 1, 1)
  } else if (structure == "D" && is.finite(ri_var)) {
    G <- matrix(ri_var, 1, 1)
  }
  cov_eig <- if (!is.null(G)) eigen((G + t(G))/2, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  cov_min <- if (all(is.finite(cov_eig))) min(cov_eig) else NA_real_
  cov_max <- if (all(is.finite(cov_eig))) max(cov_eig) else NA_real_

  rp <- tryCatch(lme4::rePCA(fit), error = function(e) NULL)
  rp_sd <- if (!is.null(rp)) unlist(lapply(rp, function(x) x$sd), use.names = FALSE) else NA_real_
  rp_sd <- rp_sd[is.finite(rp_sd)]
  rp_min <- if (length(rp_sd)) min(rp_sd) else NA_real_
  rp_max <- if (length(rp_sd)) max(rp_sd) else NA_real_
  rp_rank <- if (length(rp_sd) && rp_max > 0) sum(rp_sd / rp_max > tol_singular) else 0L

  derivs <- fit@optinfo$derivs
  h_eig <- NA_real_
  scaled_grad <- NA_real_
  if (!is.null(derivs$Hessian)) {
    h_eig <- tryCatch(eigen((derivs$Hessian + t(derivs$Hessian))/2, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_)
    if (!is.null(derivs$gradient)) {
      scaled_grad <- tryCatch(max(abs(solve(chol(derivs$Hessian), derivs$gradient))), error = function(e) NA_real_)
    }
  }
  h_min <- if (all(is.finite(h_eig))) min(h_eig) else NA_real_
  h_max <- if (all(is.finite(h_eig))) max(h_eig) else NA_real_
  h_pos <- if (is.finite(h_min) && is.finite(h_max)) h_min > max(1e-8, 1e-8 * abs(h_max)) else NA
  conv_ok <- !opt_bad && !nzchar(conv_msgs) && fixed_ok && !isFALSE(h_pos)

  list(fit = fit, row = data.frame(
    fit_error = "", warnings = paste(unique(warnings), collapse = " | "),
    optimizer_code = opt_code, lme4_convergence_messages = conv_msgs, convergence_ok = conv_ok,
    fixed_effect_finite = fixed_ok, slope = est, SE = se, df = ddf, t = tv, P = pv,
    CI_lower = est - 1.96 * se, CI_upper = est + 1.96 * se,
    random_intercept_variance = ri_var, random_slope_variance = rs_var,
    intercept_slope_correlation = corr, residual_variance = resid_var,
    RI_SD_ratio = ri_ratio, RS_SD_ratio_per_label = rs_ratio,
    RS_SD_ratio_observed_span = 6 * rs_ratio,
    covariance_eigen_min = cov_min, covariance_eigen_max = cov_max,
    covariance_eigen_ratio = safe_ratio(cov_min, cov_max),
    rePCA_min_SD = rp_min, rePCA_max_SD = rp_max, rePCA_SD_ratio = safe_ratio(rp_min, rp_max),
    rePCA_rank = rp_rank,
    hessian_eigen_min = h_min, hessian_eigen_max = h_max, hessian_positive = h_pos,
    max_abs_scaled_gradient = scaled_grad,
    isSingular_tol_1e_4 = isSingular(fit, tol = tol_singular), stringsAsFactors = FALSE))
}

optimizer_rows <- list()
patient_rows <- list()
fit_counter <- 0L

for (layer in c("miRNA", "GENCODE")) {
  selected <- selection[selection$layer == layer, ]
  for (dataset in c("FULL", "QC_PASS_ONLY")) {
    inp <- obj$inputs[[paste(layer, dataset, sep = "_")]]
    pairs <- inp$pairs
    stopifnot(all(is.finite(pairs$centered_label)), length(unique(pairs$subject)) == 30L)
    for (gradient_model in c("M3", "M4")) {
      mat <- if (gradient_model == "M3") inp$raw else inp$clr
      for (i in seq_len(nrow(selected))) {
        feature <- selected$feature[i]
        y <- as.numeric(mat[match(feature, rownames(mat)), ])
        dat <- data.frame(y = y, centered_label = as.numeric(pairs$centered_label), subject = factor(pairs$subject))
        response_sd <- sd(y)
        near_constant <- response_sd <= sqrt(.Machine$double.eps) * max(1, abs(mean(y)))

        # Patient-level slope feasibility for this feature/arm/gradient.
        by_sub <- split(seq_len(nrow(dat)), dat$subject)
        patient_slopes <- lapply(names(by_sub), function(s) {
          ix <- by_sub[[s]]
          lmfit <- lm(y ~ centered_label, data = dat[ix, , drop = FALSE])
          sm <- summary(lmfit)$coefficients
          data.frame(subject = s, n_timepoints = length(ix), n_unique_labels = length(unique(dat$centered_label[ix])),
            slope = unname(coef(lmfit)["centered_label"]),
            slope_SE = if (nrow(sm) >= 2L) unname(sm["centered_label", "Std. Error"]) else NA_real_, stringsAsFactors = FALSE)
        })
        ps <- do.call(rbind, patient_slopes)
        slope_mean <- mean(ps$slope)
        slope_se_mean <- sd(ps$slope) / sqrt(nrow(ps))
        slope_t <- slope_mean / slope_se_mean
        slope_df <- nrow(ps) - 1L
        slope_p <- 2 * pt(-abs(slope_t), df = slope_df)
        slope_ci <- slope_mean + c(-1, 1) * qt(.975, df = slope_df) * slope_se_mean

        marginal <- lm(y ~ centered_label + subject, data = dat)
        cr2 <- tryCatch(clubSandwich::coef_test(marginal, vcov = "CR2", cluster = dat$subject, test = "Satterthwaite"),
          error = function(e) e)
        cr2_est <- cr2_se <- cr2_df <- cr2_p <- NA_real_
        cr2_error <- ""
        if (inherits(cr2, "error")) {
          cr2_error <- conditionMessage(cr2)
        } else {
          cr2_row <- cr2[rownames(cr2) == "centered_label", , drop = FALSE]
          if (nrow(cr2_row) == 1L) {
            cr2_est <- cr2_row$beta
            cr2_se <- cr2_row$SE
            cr2_df <- cr2_row$df_Satt
            cr2_p <- cr2_row$p_Satt
          }
        }
        patient_rows[[length(patient_rows) + 1L]] <- data.frame(
          layer = layer, dataset = dataset, gradient_model = gradient_model, feature = feature,
          n_subjects = nrow(ps), min_timepoints = min(ps$n_timepoints), max_timepoints = max(ps$n_timepoints),
          all_subject_slopes_estimable = all(is.finite(ps$slope)),
          individual_slope_mean = slope_mean, individual_slope_SE_of_mean = slope_se_mean,
          individual_slope_CI_lower = slope_ci[1], individual_slope_CI_upper = slope_ci[2],
          individual_slope_t = slope_t, individual_slope_df = slope_df, individual_slope_P = slope_p,
          individual_slope_median = median(ps$slope), individual_slope_min = min(ps$slope), individual_slope_max = max(ps$slope),
          individual_slope_negative_n = sum(ps$slope < 0), individual_slope_positive_n = sum(ps$slope > 0),
          individual_slope_SE_min = min(ps$slope_SE, na.rm = TRUE), individual_slope_SE_median = median(ps$slope_SE, na.rm = TRUE),
          individual_slope_SE_max = max(ps$slope_SE, na.rm = TRUE),
          CR2_slope = cr2_est, CR2_SE = cr2_se, CR2_df_Satterthwaite = cr2_df, CR2_P = cr2_p, CR2_error = cr2_error,
          stringsAsFactors = FALSE)

        for (structure in names(model_formulas)) {
          for (optimizer in optimizers) {
            fr <- fit_once(dat, structure, optimizer)
            fit_counter <- fit_counter + 1L
            optimizer_rows[[fit_counter]] <- cbind(data.frame(
              layer = layer, dataset = dataset, gradient_model = gradient_model, feature = feature,
              selection_reason = selected$selection_reason[i], structure = structure, optimizer = optimizer,
              n_subjects = length(unique(dat$subject)), n_pairs = nrow(dat), n_unique_collection_labels = length(unique(dat$centered_label)),
              centered_label_mean = mean(dat$centered_label), centered_label_min = min(dat$centered_label), centered_label_max = max(dat$centered_label),
              subject_factor_levels = nlevels(dat$subject), missing_y = sum(!is.finite(dat$y)), response_mean = mean(dat$y),
              response_SD = response_sd, response_nearly_constant = near_constant,
              stringsAsFactors = FALSE), fr$row)
          }
        }
        if (fit_counter %% 120L == 0L) cat(format(Sys.time()), "fits", fit_counter, "\n")
      }
    }
  }
}

optimizer_diag <- do.call(rbind, optimizer_rows)
write_csv(optimizer_diag, file.path(out, "OPTIMIZER_DIAGNOSTICS.csv"))
patient_diag <- do.call(rbind, patient_rows)
write_csv(patient_diag, file.path(out, "intermediate", "PATIENT_LEVEL_PILOT_RESULTS.csv"))

# Choose the best numerical fit within each structure without using P values.
group_key <- interaction(optimizer_diag$layer, optimizer_diag$dataset, optimizer_diag$gradient_model,
  optimizer_diag$feature, optimizer_diag$structure, drop = TRUE, lex.order = TRUE)
groups <- split(seq_len(nrow(optimizer_diag)), group_key)
comparison_rows <- lapply(groups, function(ix) {
  z <- optimizer_diag[ix, , drop = FALSE]
  finite <- is.finite(z$slope) & is.finite(z$SE)
  stable_direction <- length(unique(direction(z$slope[finite]))) <= 1L
  med_se <- if (any(finite)) median(z$SE[finite]) else NA_real_
  slope_norm_range <- if (sum(finite) >= 2L) diff(range(z$slope[finite])) / med_se else NA_real_
  se_ratio <- if (sum(finite) >= 2L && min(z$SE[finite]) > 0) max(z$SE[finite]) / min(z$SE[finite]) else NA_real_
  optimizer_stable <- sum(finite) >= 2L && stable_direction && is.finite(slope_norm_range) && slope_norm_range <= .10 &&
    is.finite(se_ratio) && se_ratio <= 1.10
  grad_rank <- ifelse(is.finite(z$max_abs_scaled_gradient), z$max_abs_scaled_gradient, Inf)
  hess_rank <- ifelse(is.na(z$hessian_positive), 1L, ifelse(z$hessian_positive, 0L, 2L))
  conv_rank <- ifelse(z$convergence_ok & finite, 0L, 1L)
  priority <- match(z$optimizer, optimizers)
  best <- order(conv_rank, hess_rank, grad_rank, priority)[1]
  b <- z[best, , drop = FALSE]
  b$optimizer_finite_n <- sum(finite)
  b$optimizer_direction_consistent <- stable_direction
  b$optimizer_slope_range_over_median_SE <- slope_norm_range
  b$optimizer_SE_max_min_ratio <- se_ratio
  b$optimizer_stable <- optimizer_stable
  b$selected_optimizer_reason <- "convergence_then_Hessian_then_scaled_gradient_then_fixed_priority"
  b
})
comparison <- do.call(rbind, comparison_rows)
rownames(comparison) <- NULL
write_csv(comparison, file.path(out, "RANDOM_EFFECT_MODEL_COMPARISON.csv"))

# Boundary classification based on selected A-D fits.
class_key <- interaction(comparison$layer, comparison$dataset, comparison$gradient_model, comparison$feature,
  drop = TRUE, lex.order = TRUE)
class_groups <- split(seq_len(nrow(comparison)), class_key)
class_rows <- lapply(class_groups, function(ix) {
  z <- comparison[ix, , drop = FALSE]
  get <- function(s) z[z$structure == s, , drop = FALSE]
  A <- get("A"); B <- get("B"); C <- get("C"); Dm <- get("D")
  if (!all(c(nrow(A), nrow(B), nrow(C), nrow(Dm)) == 1L)) stop("Structure row missing")
  ri_zero <- is.finite(B$RI_SD_ratio) && B$RI_SD_ratio <= tol_singular
  rs_zero <- is.finite(B$RS_SD_ratio_per_label) && B$RS_SD_ratio_per_label <= tol_singular
  corr_boundary <- is.finite(A$intercept_slope_correlation) && abs(A$intercept_slope_correlation) >= .999
  A_rank_boundary <- isTRUE(A$isSingular_tol_1e_4) && is.finite(A$rePCA_SD_ratio) && A$rePCA_SD_ratio <= tol_singular
  numerical_failure <- (!A$convergence_ok && !B$convergence_ok) || (!A$optimizer_stable && !B$optimizer_stable) ||
    (!A$fixed_effect_finite && !B$fixed_effect_finite)
  primary <- if (numerical_failure) {
    "E_NUMERICAL_NONCONVERGENCE_OR_UNSTABLE_FIXED_EFFECT"
  } else if (ri_zero && rs_zero) {
    "C_BOTH_RANDOM_VARIANCES_AT_BOUNDARY"
  } else if (ri_zero && !rs_zero) {
    "A_RANDOM_INTERCEPT_AT_BOUNDARY"
  } else if (!ri_zero && rs_zero) {
    "B_RANDOM_SLOPE_AT_BOUNDARY"
  } else if ((corr_boundary || A_rank_boundary) && B$convergence_ok && B$optimizer_stable && !isTRUE(B$isSingular_tol_1e_4)) {
    "D_CORRELATION_BOUNDARY_RESOLVED_BY_UNCORRELATED"
  } else if (A$convergence_ok && A$optimizer_stable && !isTRUE(A$isSingular_tol_1e_4)) {
    "REGULAR_CORRELATED_RANDOM_SLOPE"
  } else if (B$convergence_ok && B$optimizer_stable && !isTRUE(B$isSingular_tol_1e_4)) {
    "REGULAR_UNCORRELATED_RANDOM_SLOPE"
  } else {
    "E_UNRESOLVED_BOUNDARY_OR_IDENTIFIABILITY"
  }
  C_eligible <- primary == "A_RANDOM_INTERCEPT_AT_BOUNDARY" && C$convergence_ok && C$optimizer_stable &&
    !isTRUE(C$isSingular_tol_1e_4) && is.finite(C$RS_SD_ratio_per_label) && C$RS_SD_ratio_per_label > tol_singular
  random_slope_identifiable <- primary %in% c("REGULAR_CORRELATED_RANDOM_SLOPE", "REGULAR_UNCORRELATED_RANDOM_SLOPE",
    "D_CORRELATION_BOUNDARY_RESOLVED_BY_UNCORRELATED") || C_eligible
  old <- old_diag[old_diag$layer == z$layer[1] & old_diag$dataset == z$dataset[1] &
    old_diag$model == z$gradient_model[1] & old_diag$feature == z$feature[1], , drop = FALSE]
  old_available <- nrow(old) == 1L
  old_valid <- if (old_available) as.logical(old$model_valid[1]) else NA
  old_overconservative <- old_available && !old_valid && random_slope_identifiable
  recommendation <- if (primary == "REGULAR_CORRELATED_RANDOM_SLOPE") "A" else if (
    primary %in% c("D_CORRELATION_BOUNDARY_RESOLVED_BY_UNCORRELATED", "REGULAR_UNCORRELATED_RANDOM_SLOPE")) "B" else if (
      C_eligible) "C_ELIGIBLE_POSTHOC" else if (primary == "B_RANDOM_SLOPE_AT_BOUNDARY") "RANDOM_SLOPE_UNRESOLVED_RS_ZERO" else
        "RANDOM_SLOPE_UNRESOLVED"
  data.frame(
    layer = z$layer[1], dataset = z$dataset[1], gradient_model = z$gradient_model[1], feature = z$feature[1],
    n_subjects = z$n_subjects[1], n_pairs = z$n_pairs[1], n_unique_collection_labels = z$n_unique_collection_labels[1],
    response_SD = z$response_SD[1], response_nearly_constant = z$response_nearly_constant[1],
    A_slope = A$slope, A_SE = A$SE, A_singular = A$isSingular_tol_1e_4, A_correlation = A$intercept_slope_correlation,
    A_rePCA_SD_ratio = A$rePCA_SD_ratio, A_optimizer_stable = A$optimizer_stable, A_convergence_ok = A$convergence_ok,
    B_slope = B$slope, B_SE = B$SE, B_singular = B$isSingular_tol_1e_4,
    B_RI_variance = B$random_intercept_variance, B_RS_variance = B$random_slope_variance,
    B_RI_SD_ratio = B$RI_SD_ratio, B_RS_SD_ratio_per_label = B$RS_SD_ratio_per_label,
    B_RS_SD_ratio_observed_span = B$RS_SD_ratio_observed_span,
    B_covariance_eigen_min = B$covariance_eigen_min, B_covariance_eigen_max = B$covariance_eigen_max,
    B_rePCA_SD_ratio = B$rePCA_SD_ratio, B_optimizer_stable = B$optimizer_stable, B_convergence_ok = B$convergence_ok,
    C_slope = C$slope, C_SE = C$SE, C_singular = C$isSingular_tol_1e_4,
    C_RS_variance = C$random_slope_variance, C_RS_SD_ratio_per_label = C$RS_SD_ratio_per_label,
    C_optimizer_stable = C$optimizer_stable, C_convergence_ok = C$convergence_ok,
    D_slope = Dm$slope, D_SE = Dm$SE, D_singular = Dm$isSingular_tol_1e_4,
    ri_zero_standardized = ri_zero, rs_zero_standardized = rs_zero,
    correlation_boundary = corr_boundary || A_rank_boundary,
    primary_boundary_class = primary, C_structure_eligible = C_eligible,
    random_slope_identifiable = random_slope_identifiable,
    old_stage1_result_available = old_available, old_stage1_model_valid = old_valid,
    old_blanket_singular_rule_overconservative = old_overconservative,
    recommended_structure_for_this_diagnostic = recommendation,
    stringsAsFactors = FALSE)
})
classification <- do.call(rbind, class_rows)
rownames(classification) <- NULL
write_csv(classification, file.path(out, "BOUNDARY_CLASSIFICATION.csv"))

# Compact diagnostic summaries for reporting.
summary_class <- as.data.frame(table(classification$dataset, classification$gradient_model,
  classification$layer, classification$primary_boundary_class), stringsAsFactors = FALSE)
names(summary_class) <- c("dataset", "gradient_model", "layer", "boundary_class", "n")
summary_class <- summary_class[summary_class$n > 0, ]
write_csv(summary_class, file.path(out, "intermediate", "BOUNDARY_CLASS_COUNTS.csv"))

old_fail <- classification[classification$dataset == "FULL" & classification$old_stage1_result_available &
  !classification$old_stage1_model_valid, ]
old_fail_summary <- as.data.frame(table(old_fail$gradient_model, old_fail$layer, old_fail$primary_boundary_class), stringsAsFactors = FALSE)
names(old_fail_summary) <- c("gradient_model", "layer", "boundary_class", "n")
old_fail_summary <- old_fail_summary[old_fail_summary$n > 0, ]
write_csv(old_fail_summary, file.path(out, "intermediate", "OLD_FAILED_MODEL_BOUNDARY_COUNTS.csv"))

structure_summary <- aggregate(cbind(convergence_ok = as.integer(comparison$convergence_ok),
  nonsingular = as.integer(!comparison$isSingular_tol_1e_4),
  optimizer_stable = as.integer(comparison$optimizer_stable)) ~ dataset + gradient_model + layer + structure,
  data = comparison, FUN = function(x) c(n = length(x), sum = sum(x, na.rm = TRUE), rate = mean(x, na.rm = TRUE)))
flatten <- function(df, col) {
  m <- df[[col]]
  df[[paste0(col, "_n")]] <- m[, "n"]
  df[[paste0(col, "_sum")]] <- m[, "sum"]
  df[[paste0(col, "_rate")]] <- m[, "rate"]
  df[[col]] <- NULL
  df
}
structure_summary <- flatten(structure_summary, "convergence_ok")
structure_summary <- flatten(structure_summary, "nonsingular")
structure_summary <- flatten(structure_summary, "optimizer_stable")
write_csv(structure_summary, file.path(out, "intermediate", "STRUCTURE_DIAGNOSTIC_RATES.csv"))

# Basic model-comparison plots for diagnostic use only.
png(file.path(out, "plots", "A_VS_B_FIXED_SLOPES.png"), width = 2200, height = 2100, res = 180)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
for (layer in c("miRNA", "GENCODE")) for (gm in c("M3", "M4")) {
  cc <- classification[classification$dataset == "FULL" & classification$layer == layer & classification$gradient_model == gm, ]
  plot(cc$A_slope, cc$B_slope, pch = 16, cex = .8, xlab = "Model A slope", ylab = "Model B slope",
    main = paste(layer, gm, "FULL"))
  abline(0, 1, col = "red", lwd = 1.5)
}
dev.off()

png(file.path(out, "plots", "STANDARDIZED_RANDOM_EFFECT_COMPONENTS.png"), width = 2200, height = 1100, res = 180)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (gm in c("M3", "M4")) {
  cc <- classification[classification$dataset == "FULL" & classification$gradient_model == gm, ]
  plot(pmax(cc$B_RI_SD_ratio, 1e-10), pmax(cc$B_RS_SD_ratio_per_label, 1e-10), log = "xy", pch = 16,
    col = ifelse(cc$layer == "miRNA", "#2C7FB8", "#D95F0E"),
    xlab = "Random-intercept SD / residual SD", ylab = "Random-slope SD / residual SD",
    main = paste(gm, "FULL Model B"))
  abline(v = tol_singular, h = tol_singular, lty = 2, col = "grey45")
  legend("bottomright", legend = c("miRNA", "GENCODE"), col = c("#2C7FB8", "#D95F0E"), pch = 16, bty = "n")
}
dev.off()

capture.output({
  cat("Stage 1A session captured: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
  cat("Singularity tolerance:", tol_singular, "\n")
  cat("Seed:", seed, "\n")
  cat("Key package versions:\n")
  for (p in c("lme4", "lmerTest", "clubSandwich", "pbkrtest")) {
    cat(p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT_AVAILABLE_IN_THIS_LIBRARY_PATH", "\n")
  }
  cat("\n")
  print(sessionInfo())
}, file = file.path(out, "logs", "STAGE1A_sessionInfo.txt"))

cat("STAGE1A_DIAGNOSTICS_COMPLETE\n")
cat("Selected features:", nrow(selection), "\n")
cat("Optimizer fits:", nrow(optimizer_diag), "\n")
cat("Comparison rows:", nrow(comparison), "\n")
cat("Classification rows:", nrow(classification), "\n")

