options(stringsAsFactors = FALSE, warn = 1)

out <- file.path("revision_stage1", "stage1b_full_refit")
z <- readRDS(file.path(out, "intermediate", "STAGE1B_COMPACT_RESULTS.rds"))
checks <- data.frame(check = character(), status = character(), evidence = character(), stringsAsFactors = FALSE)
add <- function(check, ok, evidence) {
  checks <<- rbind(checks, data.frame(check = check, status = ifelse(ok, "PASS", "FAIL"),
    evidence = as.character(evidence), stringsAsFactors = FALSE))
}

add("task_count", nrow(z$results_long) == 6344L, nrow(z$results_long))
add("task_id_unique", length(unique(z$results_long$task_id)) == 6344L, length(unique(z$results_long$task_id)))
add("task_errors", sum(nzchar(z$results_long$task_error), na.rm = TRUE) == 0L,
  sum(nzchar(z$results_long$task_error), na.rm = TRUE))
add("miRNA_full_rows", nrow(z$full_tables$miRNA) == 507L, nrow(z$full_tables$miRNA))
add("GENCODE_full_rows", nrow(z$full_tables$GENCODE) == 1079L, nrow(z$full_tables$GENCODE))
add("strict_core_rows", nrow(z$strict_core) == 49L,
  sprintf("%d (miRNA=%d; GENCODE=%d)", nrow(z$strict_core), sum(z$strict_core$layer == "miRNA"), sum(z$strict_core$layer == "GENCODE")))
add("boundary_supported_rows", nrow(z$boundary_supported) == 65L,
  sprintf("%d (miRNA=%d; GENCODE=%d)", nrow(z$boundary_supported), sum(z$boundary_supported$layer == "miRNA"), sum(z$boundary_supported$layer == "GENCODE")))
add("representative_rows", nrow(z$representatives) == 10L, nrow(z$representatives))

for (layer in c("miRNA", "GENCODE")) {
  x <- z$full_tables[[layer]]
  nfam <- if (layer == "miRNA") 507L else 1079L
  add(paste0(layer, "_features_unique"), length(unique(x$feature)) == nfam, length(unique(x$feature)))
  for (m in c("M3", "M4")) {
    p <- x[[paste0(m, "_P")]]
    expected <- p.adjust(p, method = "BH", n = nfam)
    observed <- x[[paste0(m, "_FDR")]]
    delta <- suppressWarnings(max(abs(expected - observed), na.rm = TRUE))
    if (!is.finite(delta)) delta <- 0
    add(paste(layer, m, "BH_family", sep = "_"), delta < 1e-14,
      sprintf("n=%d; max_abs_diff=%.3g", nfam, delta))
    bad <- x[[paste0(m, "_model_status")]] %in% c("NUMERICAL_UNSTABLE", "MODEL_NOT_EVALUABLE")
    add(paste(layer, m, "unstable_P_is_NA", sep = "_"), all(is.na(p[bad])),
      sprintf("unstable=%d; nonNA_P=%d", sum(bad), sum(!is.na(p[bad]))))
  }
}

required <- c(
  "STAGE1B_FULL_MODEL_RESULTS_miRNA.csv", "STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv",
  "STAGE1B_QC_MODEL_RESULTS.csv", "STAGE1B_CR2_RESULTS.csv",
  "STAGE1B_MODEL_BOUNDARY_DIAGNOSTICS.csv", "STAGE1B_MODEL_COMPARISON.csv",
  "STAGE1B_STRICT_CORE_ROBUST.csv", "STAGE1B_BOUNDARY_SUPPORTED.csv",
  "STAGE1B_REPRESENTATIVE_FEATURES.csv", "STAGE1B_FINAL_REPORT.md", "STAGE2_GO_NO_GO.md")
for (f in required) add(paste0("required_file_", f), file.exists(file.path(out, f)), f)

write.csv(checks, file.path(out, "STAGE1B_OUTPUT_VALIDATION.csv"), row.names = FALSE,
  na = "", fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()), file.path(out, "logs", "STAGE1B_validation_sessionInfo.txt"), useBytes = TRUE)
if (any(checks$status != "PASS")) {
  print(checks[checks$status != "PASS", , drop = FALSE])
  stop("Stage 1B output validation failed")
}
cat(sprintf("VALIDATION_PASS: %d checks\n", nrow(checks)))
