options(stringsAsFactors = FALSE)

root <- "revision_stage1/stage1a_boundary_diagnostics"
checks <- list()
add <- function(name, pass, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, status = if (isTRUE(pass)) "PASS" else "FAIL",
    observed = as.character(observed), expected = as.character(expected), stringsAsFactors = FALSE)
}

required <- c("STAGE1_PLAN_DEVIATION_01.md", "BOUNDARY_CLASSIFICATION.csv",
  "RANDOM_EFFECT_MODEL_COMPARISON.csv", "OPTIMIZER_DIAGNOSTICS.csv",
  "PATIENT_LEVEL_SENSITIVITY_FEASIBILITY.md", "STAGE1A_DECISION_RULES.md", "STAGE1A_FINAL_REPORT.md")
for (f in required) add(paste0("required_", f), file.exists(file.path(root, f)), file.exists(file.path(root, f)), TRUE)

d <- read.csv(file.path(root, "BOUNDARY_CLASSIFICATION.csv"), check.names = FALSE)
cmp <- read.csv(file.path(root, "RANDOM_EFFECT_MODEL_COMPARISON.csv"), check.names = FALSE)
opt <- read.csv(file.path(root, "OPTIMIZER_DIAGNOSTICS.csv"), check.names = FALSE)
pat <- read.csv(file.path(root, "intermediate", "PATIENT_LEVEL_PILOT_RESULTS.csv"), check.names = FALSE)
sel <- read.csv(file.path(root, "intermediate", "EXPANDED_PILOT_SELECTION.csv"), check.names = FALSE)

add("selected_features", nrow(sel) == 30L, nrow(sel), 30)
add("optimizer_rows", nrow(opt) == 1440L, nrow(opt), 1440)
add("comparison_rows", nrow(cmp) == 480L, nrow(cmp), 480)
add("classification_rows", nrow(d) == 120L, nrow(d), 120)
add("patient_rows", nrow(pat) == 120L, nrow(pat), 120)
add("classification_unique_units", nrow(unique(d[c("layer","dataset","gradient_model","feature")])) == 120L,
  nrow(unique(d[c("layer","dataset","gradient_model","feature")])), 120)
add("optimizer_unique_units", nrow(unique(opt[c("layer","dataset","gradient_model","feature","structure","optimizer")])) == 1440L,
  nrow(unique(opt[c("layer","dataset","gradient_model","feature","structure","optimizer")])), 1440)
add("all_subject_counts_30", all(d$n_subjects == 30L), paste(unique(d$n_subjects), collapse = ";"), 30)
add("FULL_pairs_153", all(d$n_pairs[d$dataset == "FULL"] == 153L), paste(unique(d$n_pairs[d$dataset == "FULL"]), collapse = ";"), 153)
add("QC_pairs_141", all(d$n_pairs[d$dataset == "QC_PASS_ONLY"] == 141L), paste(unique(d$n_pairs[d$dataset == "QC_PASS_ONLY"]), collapse = ";"), 141)
add("centered_label", all(opt$centered_label_min == -3 & opt$centered_label_max == 3),
  paste(range(opt$centered_label_min, opt$centered_label_max), collapse = ";"), "-3 to 3")
add("no_missing_response", all(opt$missing_y == 0L), max(opt$missing_y), 0)
add("all_four_structures", identical(sort(unique(cmp$structure)), c("A","B","C","D")),
  paste(sort(unique(cmp$structure)), collapse = ";"), "A;B;C;D")
add("all_three_optimizers", identical(sort(unique(opt$optimizer)), sort(c("nloptwrap","bobyqa","nlminbwrap"))),
  paste(sort(unique(opt$optimizer)), collapse = ";"), "bobyqa;nlminbwrap;nloptwrap")

old <- d[d$dataset == "FULL" & d$old_stage1_result_available & !d$old_stage1_model_valid, ]
add("old_failed_n", nrow(old) == 16L, nrow(old), 16)
add("old_failed_RI_zero", sum(old$ri_zero_standardized & !old$rs_zero_standardized) == 3L,
  sum(old$ri_zero_standardized & !old$rs_zero_standardized), 3)
add("old_failed_RS_zero", sum(old$rs_zero_standardized & !old$ri_zero_standardized) == 13L,
  sum(old$rs_zero_standardized & !old$ri_zero_standardized), 13)
add("old_failed_correlation_boundary", sum(old$correlation_boundary) == 16L, sum(old$correlation_boundary), 16)
add("old_failed_numerical", sum(grepl("^E_", old$primary_boundary_class)) == 0L,
  sum(grepl("^E_", old$primary_boundary_class)), 0)
add("C_eligible_FULL", sum(d$dataset == "FULL" & d$C_structure_eligible) == 5L,
  sum(d$dataset == "FULL" & d$C_structure_eligible), 5)
add("FULL_identifiable", sum(d$dataset == "FULL" & d$random_slope_identifiable) == 39L,
  sum(d$dataset == "FULL" & d$random_slope_identifiable), 39)
add("QC_identifiable", sum(d$dataset == "QC_PASS_ONLY" & d$random_slope_identifiable) == 40L,
  sum(d$dataset == "QC_PASS_ONLY" & d$random_slope_identifiable), 40)

bfull <- cmp[cmp$dataset == "FULL" & cmp$structure == "B", ]
bqc <- cmp[cmp$dataset == "QC_PASS_ONLY" & cmp$structure == "B", ]
add("B_FULL_converged", sum(bfull$convergence_ok) == 60L, sum(bfull$convergence_ok), 60)
add("B_FULL_optimizer_stable", sum(bfull$optimizer_stable) == 57L, sum(bfull$optimizer_stable), 57)
add("B_QC_converged", sum(bqc$convergence_ok) == 60L, sum(bqc$convergence_ok), 60)
add("B_QC_optimizer_stable", sum(bqc$optimizer_stable) == 59L, sum(bqc$optimizer_stable), 59)
add("patient_slopes_all_estimable", all(pat$all_subject_slopes_estimable), sum(pat$all_subject_slopes_estimable), 120)
add("CR2_all_finite", all(is.finite(pat$CR2_slope)), sum(is.finite(pat$CR2_slope)), 120)

root_hash <- tools::md5sum("revision_stage1/STAGE1_PLAN_DEVIATION_01.md")
copy_hash <- tools::md5sum(file.path(root, "STAGE1_PLAN_DEVIATION_01.md"))
add("deviation_copy_matches_root", identical(unname(root_hash), unname(copy_hash)), unname(copy_hash), unname(root_hash))

report <- paste(readLines(file.path(root, "STAGE1A_FINAL_REPORT.md"), warn = FALSE), collapse = "\n")
add("decision_GO_TO_STAGE1B", grepl("GO_TO_STAGE1B_FULL_REFIT", report, fixed = TRUE),
  grepl("GO_TO_STAGE1B_FULL_REFIT", report, fixed = TRUE), TRUE)
add("session_info_exists", file.exists(file.path(root, "logs", "STAGE1A_sessionInfo.txt")),
  file.exists(file.path(root, "logs", "STAGE1A_sessionInfo.txt")), TRUE)

ans <- do.call(rbind, checks)
write.csv(ans, file.path(root, "intermediate", "STAGE1A_OUTPUT_VALIDATION.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat(paste(ans$status, ans$check, "observed=", ans$observed, "expected=", ans$expected), sep = "\n",
  file = file.path(root, "logs", "validate_stage1a_outputs.log"))
if (any(ans$status != "PASS")) stop("Stage 1A validation failed: ", paste(ans$check[ans$status != "PASS"], collapse = ", "))
cat("ALL_STAGE1A_CHECKS_PASS\n")

