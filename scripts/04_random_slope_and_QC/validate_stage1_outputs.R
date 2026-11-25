options(stringsAsFactors = FALSE)
root <- "revision_stage1"
checks <- list()
add <- function(name, pass, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, status = if (isTRUE(pass)) "PASS" else "FAIL",
    observed = as.character(observed), expected = as.character(expected), stringsAsFactors = FALSE)
}

required <- c(
  "STAGE1_FROZEN_ANALYSIS_PLAN.md", "STAGE1_DATASET_SUMMARY.csv",
  "STAGE1_MODEL_RESULTS_miRNA.csv", "STAGE1_MODEL_RESULTS_GENCODE.csv",
  "STAGE1_MODEL_DIAGNOSTICS.csv", "STAGE1_QC_SENSITIVITY.csv",
  "STAGE1_ROBUST_SIGNATURE.csv", "STAGE1_SLOPE_DIRECTION_SUMMARY.csv",
  "STAGE1_REPRESENTATIVE_FEATURES.csv", "STAGE1_FINAL_REPORT.md", "STAGE2_RECOMMENDATION.md"
)
for (f in required) add(paste0("required_file_", f), file.exists(file.path(root, f)), file.exists(file.path(root, f)), TRUE)

mi <- read.csv(file.path(root, "STAGE1_MODEL_RESULTS_miRNA.csv"), check.names = FALSE)
ge <- read.csv(file.path(root, "STAGE1_MODEL_RESULTS_GENCODE.csv"), check.names = FALSE)
dg <- read.csv(file.path(root, "STAGE1_MODEL_DIAGNOSTICS.csv"), check.names = FALSE)
qc <- read.csv(file.path(root, "STAGE1_QC_SENSITIVITY.csv"), check.names = FALSE)
rb <- read.csv(file.path(root, "STAGE1_ROBUST_SIGNATURE.csv"), check.names = FALSE)
rp <- read.csv(file.path(root, "STAGE1_REPRESENTATIVE_FEATURES.csv"), check.names = FALSE)
ds <- read.csv(file.path(root, "STAGE1_DATASET_SUMMARY.csv"), check.names = FALSE)
eq <- read.csv(file.path(root, "intermediate", "CLR_EQUIVALENCE_CHECK.csv"), check.names = FALSE)
repro <- read.csv(file.path(root, "intermediate", "M1_REPRODUCTION_NUMERIC_CHECK.csv"), check.names = FALSE)
rec <- read.csv(file.path(root, "intermediate", "RECENTERING_EQUIVALENCE_NUMERIC_CHECK.csv"), check.names = FALSE)
sumc <- read.csv(file.path(root, "intermediate", "MODEL_SUMMARY_COUNTS.csv"), check.names = FALSE)
gate <- read.csv(file.path(root, "intermediate", "PILOT_GATE.csv"), check.names = FALSE)

add("miRNA_result_rows", nrow(mi) == 507L, nrow(mi), 507)
add("GENCODE_result_rows", nrow(ge) == 1079L, nrow(ge), 1079)
add("miRNA_exact_unique_features", length(unique(mi$feature)) == 507L, length(unique(mi$feature)), 507)
add("GENCODE_exact_unique_features", length(unique(ge$feature)) == 1079L, length(unique(ge$feature)), 1079)
add("QC_rows", nrow(qc) == 1586L, nrow(qc), 1586)
add("robust_rows", nrow(rb) == 1586L, nrow(rb), 1586)
add("representative_rows", nrow(rp) == 10L, nrow(rp), 10)
add("diagnostics_pilot_scope", all(dg$analysis_scope == "PILOT_ONLY"), paste(unique(dg$analysis_scope), collapse = ";"), "PILOT_ONLY")
add("FULL_pairs", ds$n_pairs[ds$dataset == "FULL" & ds$layer == "miRNA" & ds$collection_label == "ALL"] == 153L,
  ds$n_pairs[ds$dataset == "FULL" & ds$layer == "miRNA" & ds$collection_label == "ALL"], 153)
add("QC_pairs", ds$n_pairs[ds$dataset == "QC_PASS_ONLY" & ds$layer == "miRNA" & ds$collection_label == "ALL"] == 141L,
  ds$n_pairs[ds$dataset == "QC_PASS_ONLY" & ds$layer == "miRNA" & ds$collection_label == "ALL"], 141)
add("M1_numeric_reproduction", max(unlist(repro[-1])) < 1e-12, max(unlist(repro[-1])), "<1e-12")
add("CLR_algebra_identity", max(eq$clr_identity_max_abs_diff) < 1e-12, max(eq$clr_identity_max_abs_diff), "<1e-12")
add("old_recenter_numeric_equivalence", max(rec$old_recenter_slope_max_abs_diff) < 1e-12,
  max(rec$old_recenter_slope_max_abs_diff), "<1e-12")

expected_counts <- data.frame(layer = c("miRNA","GENCODE","miRNA","GENCODE","miRNA","GENCODE","miRNA","GENCODE"),
  dataset = c("FULL","FULL","QC_PASS_ONLY","QC_PASS_ONLY","FULL","FULL","QC_PASS_ONLY","QC_PASS_ONLY"),
  model = c("M1","M1","M1","M1","M2","M2","M2","M2"),
  n = c(192,200,176,165,128,173,105,139))
for (i in seq_len(nrow(expected_counts))) {
  z <- sumc[sumc$layer == expected_counts$layer[i] & sumc$dataset == expected_counts$dataset[i] & sumc$model == expected_counts$model[i], ]
  add(paste("count", expected_counts$layer[i], expected_counts$dataset[i], expected_counts$model[i], sep = "_"),
    nrow(z) == 1L && z$n_significant == expected_counts$n[i], if (nrow(z)) z$n_significant else "missing", expected_counts$n[i])
}
add("pilot_gate_stopped", all(gate$gate == "STOP_SYSTEMATIC_FAILURE"), paste(unique(gate$gate), collapse = ";"), "STOP_SYSTEMATIC_FAILURE")
z <- gate[gate$layer == "GENCODE" & gate$model == "M4", ]
add("GENCODE_M4_pilot_below_50pct", nrow(z) == 1L && z$valid_rate < .5, if (nrow(z)) z$valid_rate else "missing", "<0.5")
add("full_M3_M4_not_run", all(!mi$FULL_M3_model_valid) && all(!ge$FULL_M4_model_valid),
  paste(unique(c(mi$FULL_M3_status, ge$FULL_M4_status)), collapse = ";"), "NOT_RUN_PILOT_GATE_FAILED")

plots <- list.files(file.path(root, "diagnostic_plots"), full.names = TRUE)
add("diagnostic_plot_count", length(plots) >= 6L, length(plots), ">=6")
add("diagnostic_plots_nonempty", all(file.info(plots)$size > 1000), min(file.info(plots)$size), ">1000 bytes")

out <- do.call(rbind, checks)
write.csv(out, file.path(root, "provenance", "STAGE1_OUTPUT_VALIDATION.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat(paste(out$status, out$check, "observed=", out$observed, "expected=", out$expected), sep = "\n",
  file = file.path(root, "logs", "STAGE1_output_validation.log"))
if (any(out$status != "PASS")) stop("Output validation failed: ", paste(out$check[out$status != "PASS"], collapse = ", "))
cat("ALL_OUTPUT_CHECKS_PASS\n")
