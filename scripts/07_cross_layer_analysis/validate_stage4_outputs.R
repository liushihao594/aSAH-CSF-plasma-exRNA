options(stringsAsFactors = FALSE)
root <- "revision_stage4"
read_out <- function(x) read.csv(file.path(root, x), check.names = FALSE)
checks <- list()
add <- function(name, observed, expected, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, observed = as.character(observed), expected = as.character(expected),
    status = if (isTRUE(pass)) "PASS" else "FAIL", detail = detail,
    stringsAsFactors = FALSE
  )
}

repro <- read_out("STAGE4_ORIGINAL_SCORE_REPRODUCTION.csv")
scores <- read_out("STAGE4_RAW_CLR_CROSSFITTED_SCORES.csv")
boot <- read_out("STAGE4_PATIENT_BOOTSTRAP.csv")
summary <- read_out("STAGE4_CORRELATION_SUMMARY.csv")
s12 <- read_out("STAGE4_SUBJECT12_SENSITIVITY.csv")
loo <- read_out("STAGE4_LOO_INFLUENCE.csv")
within <- read_out("STAGE4_WITHIN_SUBJECT_ANALYSIS.csv")
strict <- read_out("STAGE4_STRICT_CORE_EXPLORATORY.csv")

add("original_reproduction_all_pass", sum(repro$status != "PASS"), 0, all(repro$status == "PASS"))
add("original_reproduction_max_abs_diff", max(repro$absolute_difference), "<1e-12",
    max(repro$absolute_difference) < 1e-12)
add("primary_score_rows", nrow(scores), 588, nrow(scores) == 588)
expected_counts <- c(FULL_RAW = 153L, FULL_CLR = 153L, QC_PASS_ONLY_RAW = 141L, QC_PASS_ONLY_CLR = 141L)
for (nm in names(expected_counts)) {
  p <- strsplit(nm, "_")[[1]]
  scale <- tail(p, 1); dataset <- sub(paste0("_", scale, "$"), "", nm)
  obs <- sum(scores$dataset == dataset & scores$scale == scale)
  add(paste0("score_rows_", nm), obs, expected_counts[[nm]], obs == expected_counts[[nm]])
}
key <- paste(scores$dataset, scores$scale, scores$pair_id, sep = "|")
add("primary_score_unique_keys", sum(duplicated(key)), 0, !anyDuplicated(key))
add("primary_scores_finite", sum(!is.finite(scores$miRNA_score) | !is.finite(scores$GENCODE_score)), 0,
    all(is.finite(scores$miRNA_score)) && all(is.finite(scores$GENCODE_score)))

add("bootstrap_rows", nrow(boot), 280000, nrow(boot) == 280000)
add("bootstrap_branch_count", length(unique(boot$branch)), 14, length(unique(boot$branch)) == 14)
add("bootstrap_metric_count", length(unique(boot$metric)), 10, length(unique(boot$metric)) == 10)
tab <- table(boot$branch, boot$metric)
add("bootstrap_each_branch_metric_2000", paste(range(tab), collapse = ":"), "2000:2000", all(tab == 2000))
add("bootstrap_all_finite", sum(!is.finite(boot$estimate)), 0, all(is.finite(boot$estimate)))
add("correlation_summary_rows", nrow(summary), 140, nrow(summary) == 140)
add("summary_bootstrap_valid", paste(range(summary$bootstrap_valid), collapse = ":"), "2000:2000",
    all(summary$bootstrap_valid == 2000 & summary$bootstrap_invalid == 0))
add("summary_CI_order", sum(summary$CI_lower > summary$estimate | summary$estimate > summary$CI_upper), 0,
    all(summary$CI_lower <= summary$estimate & summary$estimate <= summary$CI_upper))
add("subject12_rows", nrow(s12), 40, nrow(s12) == 40)
add("LOO_rows", nrow(loo), 640, nrow(loo) == 640)
add("LOO_fixed_subjects_each_scale",
    paste(sapply(c("RAW", "CLR"), function(z) length(unique(loo$omitted_subject[loo$scale == z & loo$analysis_type == "FIXED_SCORE_OMISSION"]))), collapse = ":"),
    "30:30",
    all(sapply(c("RAW", "CLR"), function(z) length(unique(loo$omitted_subject[loo$scale == z & loo$analysis_type == "FIXED_SCORE_OMISSION"]))) == 30))
add("within_subject_rows", nrow(within), 56, nrow(within) == 56)
add("strict_core_summary_rows", nrow(strict), 40, nrow(strict) == 40)

raw <- summary[summary$branch == "PRIMARY_FULL_RAW", ]
get_est <- function(metric) raw$estimate[match(metric, raw$metric)]
add("RAW_pair_Pearson_exact", get_est("pair_Pearson"), "0.369168915288493",
    abs(get_est("pair_Pearson") - 0.369168915288493) < 1e-12)
add("RAW_pair_Spearman_exact", get_est("pair_Spearman"), "0.253598568614049",
    abs(get_est("pair_Spearman") - 0.253598568614049) < 1e-12)
add("RAW_slope_Pearson_exact", get_est("slope_Pearson"), "0.477030937887742",
    abs(get_est("slope_Pearson") - 0.477030937887742) < 1e-12)
add("RAW_direction_count_exact", get_est("slope_same_direction_N"), 18,
    get_est("slope_same_direction_N") == 18)

required <- c(
  "STAGE4_FROZEN_ANALYSIS_PLAN.md", "STAGE4_ORIGINAL_SCORE_REPRODUCTION.csv",
  "STAGE4_RAW_CLR_CROSSFITTED_SCORES.csv", "STAGE4_PATIENT_BOOTSTRAP.csv",
  "STAGE4_CORRELATION_SUMMARY.csv", "STAGE4_SUBJECT12_SENSITIVITY.csv",
  "STAGE4_LOO_INFLUENCE.csv", "STAGE4_WITHIN_SUBJECT_ANALYSIS.csv",
  "STAGE4_STRICT_CORE_EXPLORATORY.csv", "STAGE4_FIGURE5_DECISION.md",
  "STAGE4_FINAL_REPORT.md", "STAGE5_RECOMMENDATION.md",
  "intermediate/STAGE4_RESULTS_OBJECT.rds", "logs/STAGE4_sessionInfo.txt", "logs/STAGE4_execution.log"
)
missing <- required[!file.exists(file.path(root, required))]
add("required_outputs_present", length(missing), 0, length(missing) == 0, paste(missing, collapse = ";"))
plots <- list.files(file.path(root, "plots"), pattern = "\\.(pdf|svg)$", full.names = TRUE)
add("plot_file_count", length(plots), 14, length(plots) == 14)
add("plot_minimum_size", min(file.info(plots)$size), ">10000", length(plots) == 14 && all(file.info(plots)$size > 10000))

out <- do.call(rbind, checks)
write.csv(out, file.path(root, "STAGE4_OUTPUT_VALIDATION.csv"), row.names = FALSE, na = "")
if (any(out$status != "PASS")) {
  print(out[out$status != "PASS", ])
  stop("Stage 4 validation failed")
}
cat(sprintf("STAGE4_VALIDATION_PASS checks=%d\n", nrow(out)))
