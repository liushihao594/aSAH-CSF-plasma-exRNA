options(stringsAsFactors = FALSE)

root <- "revision_stage3"
read_out <- function(x) read.csv(file.path(root, x), check.names = FALSE)
checks <- list()
add_check <- function(name, observed, expected, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, observed = as.character(observed), expected = as.character(expected),
    status = if (isTRUE(pass)) "PASS" else "FAIL", detail = detail,
    stringsAsFactors = FALSE
  )
}

targets <- read.csv(file.path(root, "intermediate", "STAGE3_TARGET_FEATURES.csv"), check.names = FALSE)
slopes <- read_out("STAGE3_CSF_PLASMA_SLOPES.csv")
inter <- read_out("STAGE3_BIOFLUID_INTERACTION.csv")
mir_prof <- read_out("STAGE3_CORE_MIRNA_PROFILES.csv")
gen_prof <- read_out("STAGE3_GENCODE_PROFILES.csv")
decomp <- read_out("STAGE3_GRADIENT_DECOMPOSITION.csv")
endpoint <- read_out("STAGE3_ENDPOINT_SENSITIVITY.csv")
endpoint_subjects <- read.csv(file.path(root, "intermediate", "STAGE3_ENDPOINT_SUBJECTS.csv"))
complete7 <- read.csv(file.path(root, "intermediate", "STAGE3_COMPLETE7_SUBJECTS.csv"))

add_check("target_total", nrow(targets), 54, nrow(targets) == 54)
add_check("target_miRNA", sum(targets$layer == "miRNA"), 19, sum(targets$layer == "miRNA") == 19)
add_check("target_GENCODE", sum(targets$layer == "GENCODE"), 35, sum(targets$layer == "GENCODE") == 35)
add_check("separate_model_rows", nrow(slopes), 432, nrow(slopes) == 432)
add_check("interaction_rows", nrow(inter), 216, nrow(inter) == 216)
add_check("decomposition_rows", nrow(decomp), 108, nrow(decomp) == 108)
add_check("miRNA_profile_rows", nrow(mir_prof), 896, nrow(mir_prof) == 896)
add_check("GENCODE_profile_rows", nrow(gen_prof), 1848, nrow(gen_prof) == 1848)
add_check("endpoint_rows", nrow(endpoint), 4212, nrow(endpoint) == 4212)
add_check("endpoint_subjects", nrow(endpoint_subjects), 12, nrow(endpoint_subjects) == 12)
add_check("complete_seven_label_subjects", nrow(complete7), 4, nrow(complete7) == 4)

key_slopes <- paste(slopes$feature, slopes$layer, slopes$dataset, slopes$scale, slopes$biofluid, sep = "|")
key_inter <- paste(inter$feature, inter$layer, inter$dataset, inter$scale, sep = "|")
add_check("separate_model_unique_keys", sum(duplicated(key_slopes)), 0, !anyDuplicated(key_slopes))
add_check("interaction_unique_keys", sum(duplicated(key_inter)), 0, !anyDuplicated(key_inter))
add_check("all_separate_models_reportable", sum(slopes$model_status %in% c("NUMERICAL_UNSTABLE", "MODEL_NOT_EVALUABLE")), 0,
          !any(slopes$model_status %in% c("NUMERICAL_UNSTABLE", "MODEL_NOT_EVALUABLE")))
add_check("data_identity", max(abs(inter$data_identity_max_abs), na.rm = TRUE), "<1e-12",
          max(abs(inter$data_identity_max_abs), na.rm = TRUE) < 1e-12)
add_check("joint_gradient_CR2_identity", max(abs(inter$joint_minus_gradient_CR2), na.rm = TRUE), "<1e-12",
          max(abs(inter$joint_minus_gradient_CR2), na.rm = TRUE) < 1e-12)

recalc <- inter$interaction_FDR
for (d in unique(inter$dataset)) for (s in unique(inter$scale)) for (l in unique(inter$layer)) {
  idx <- inter$dataset == d & inter$scale == s & inter$layer == l
  recalc[idx] <- p.adjust(inter$interaction_P[idx], method = "BH", n = sum(idx))
}
add_check("interaction_BH_recalculation", max(abs(recalc - inter$interaction_FDR), na.rm = TRUE), "<1e-12",
          max(abs(recalc - inter$interaction_FDR), na.rm = TRUE) < 1e-12)

recalc_s <- slopes$FDR
for (d in unique(slopes$dataset)) for (s in unique(slopes$scale)) for (l in unique(slopes$layer)) for (b in unique(slopes$biofluid)) {
  idx <- slopes$dataset == d & slopes$scale == s & slopes$layer == l & slopes$biofluid == b
  recalc_s[idx] <- p.adjust(slopes$P[idx], method = "BH", n = sum(idx))
}
add_check("separate_model_BH_recalculation", max(abs(recalc_s - slopes$FDR), na.rm = TRUE), "<1e-12",
          max(abs(recalc_s - slopes$FDR), na.rm = TRUE) < 1e-12)

strict_mir <- decomp[decomp$layer == "miRNA" & decomp$strict_core & decomp$dataset == "FULL" & decomp$scale == "RAW", ]
class_tab <- table(strict_mir$contribution_class)
expected_classes <- c(A_CSF_DECREASE_DOMINANT = 7L,
                      B_PLASMA_INCREASE_DOMINANT = 2L,
                      C_CSF_DECREASE_AND_PLASMA_INCREASE = 7L)
for (nm in names(expected_classes)) {
  obs <- if (nm %in% names(class_tab)) unname(class_tab[[nm]]) else 0L
  add_check(paste0("strict_miRNA_class_", nm), obs, expected_classes[[nm]], obs == expected_classes[[nm]])
}

required <- c(
  "STAGE3_FROZEN_ANALYSIS_PLAN.md", "STAGE3_CSF_PLASMA_SLOPES.csv",
  "STAGE3_BIOFLUID_INTERACTION.csv", "STAGE3_CORE_MIRNA_PROFILES.csv",
  "STAGE3_GENCODE_PROFILES.csv", "STAGE3_GRADIENT_DECOMPOSITION.csv",
  "STAGE3_CLR_SENSITIVITY.csv", "STAGE3_QC_SENSITIVITY.csv",
  "STAGE3_ENDPOINT_SENSITIVITY.csv", "STAGE3_REPRESENTATIVE_FEATURES.csv",
  "STAGE3_FINAL_REPORT.md", "STAGE4_RECOMMENDATION.md",
  "intermediate/STAGE3_RESULTS_OBJECT.rds", "logs/STAGE3_sessionInfo.txt"
)
missing <- required[!file.exists(file.path(root, required))]
add_check("required_outputs_present", length(missing), 0, length(missing) == 0,
          paste(missing, collapse = ";"))

plot_files <- list.files(file.path(root, "plots"), pattern = "\\.(pdf|svg)$", full.names = TRUE)
add_check("plot_file_count", length(plot_files), 8, length(plot_files) == 8)
add_check("plot_minimum_size", if (length(plot_files)) min(file.info(plot_files)$size) else 0, ">10000",
          length(plot_files) == 8 && all(file.info(plot_files)$size > 10000))

out <- do.call(rbind, checks)
write.csv(out, file.path(root, "STAGE3_OUTPUT_VALIDATION.csv"), row.names = FALSE, na = "")
if (any(out$status != "PASS")) {
  print(out[out$status != "PASS", ])
  stop("Stage 3 validation failed")
}
cat(sprintf("STAGE3_VALIDATION_PASS checks=%d\n", nrow(out)))
