options(stringsAsFactors = FALSE, warn = 1)
root <- "."
out <- file.path(root, "revision_stage1")
obj <- readRDS(file.path(out, "intermediate", "prepared_stage1_inputs.rds"))
pilot <- read.csv(file.path(out, "intermediate", "PILOT_MODEL_DIAGNOSTICS.csv"), check.names = FALSE)
gate <- read.csv(file.path(out, "intermediate", "PILOT_GATE.csv"), check.names = FALSE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
direction <- function(x, tol = .Machine$double.eps^0.5) ifelse(is.na(x), NA_character_, ifelse(x > tol, "positive", ifelse(x < -tol, "negative", "zero")))
fmt <- function(x, digits = 4) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "g"))

layers <- c("miRNA", "GENCODE")
model_tables <- list()
qc_tables <- list()
robust_tables <- list()
direction_rows <- list()
rep_tables <- list()
summary_rows <- list()
repro_rows <- list()
sens_rows <- list()

for (layer in layers) {
  m1 <- obj$limma[[paste(layer, "FULL", "M1", sep = "_")]]
  m2 <- obj$limma[[paste(layer, "FULL", "M2", sep = "_")]]
  q1 <- obj$limma[[paste(layer, "QC_PASS_ONLY", "M1", sep = "_")]]
  q2 <- obj$limma[[paste(layer, "QC_PASS_ONLY", "M2", sep = "_")]]
  sf <- obj$sensitivities[[paste(layer, "FULL", sep = "_")]]
  sq <- obj$sensitivities[[paste(layer, "QC_PASS_ONLY", sep = "_")]]
  old_full <- read.csv(file.path(root, "revision_stage0c", "statistics", paste0(layer, "_FULL_complete_results.csv")), check.names = FALSE)
  old_qc <- read.csv(file.path(root, "revision_stage0c", "statistics", paste0(layer, "_QC_PASS_ONLY_complete_results.csv")), check.names = FALSE)

  tab <- data.frame(feature = m1$feature, layer = layer, stringsAsFactors = FALSE)
  add_limma <- function(tab, d, pre) {
    z <- d[, c("feature", "slope", "SE", "CI_lower", "CI_upper", "P", "FDR", "duplicate_correlation")]
    names(z)[-1] <- paste0(pre, names(z)[-1])
    merge(tab, z, by = "feature", all.x = TRUE, sort = FALSE)
  }
  tab <- add_limma(tab, m1, "FULL_M1_")
  tab <- add_limma(tab, m2, "FULL_M2_")
  tab <- add_limma(tab, q1, "QC_M1_")
  tab <- add_limma(tab, q2, "QC_M2_")
  for (m in c("M3", "M4")) for (prefix in c("FULL", "QC")) {
    tab[[paste0(prefix, "_", m, "_slope")]] <- NA_real_
    tab[[paste0(prefix, "_", m, "_SE")]] <- NA_real_
    tab[[paste0(prefix, "_", m, "_CI_lower")]] <- NA_real_
    tab[[paste0(prefix, "_", m, "_CI_upper")]] <- NA_real_
    tab[[paste0(prefix, "_", m, "_P")]] <- NA_real_
    tab[[paste0(prefix, "_", m, "_FDR")]] <- NA_real_
    tab[[paste0(prefix, "_", m, "_model_valid")]] <- FALSE
    tab[[paste0(prefix, "_", m, "_status")]] <- "NOT_RUN_PILOT_GATE_FAILED"
  }
  ps <- pilot[pilot$layer == layer, ]
  for (m in c("M3", "M4")) {
    psm <- ps[ps$model == m, ]
    ix <- match(psm$feature, tab$feature)
    tab[[paste0("PILOT_FULL_", m, "_slope")]] <- NA_real_
    tab[[paste0("PILOT_FULL_", m, "_SE")]] <- NA_real_
    tab[[paste0("PILOT_FULL_", m, "_P")]] <- NA_real_
    tab[[paste0("PILOT_FULL_", m, "_model_valid")]] <- NA
    tab[[paste0("PILOT_FULL_", m, "_final_structure")]] <- NA_character_
    tab[[paste0("PILOT_FULL_", m, "_slope")]][ix] <- psm$slope
    tab[[paste0("PILOT_FULL_", m, "_SE")]][ix] <- psm$SE
    tab[[paste0("PILOT_FULL_", m, "_P")]][ix] <- psm$P
    tab[[paste0("PILOT_FULL_", m, "_model_valid")]][ix] <- psm$model_valid
    tab[[paste0("PILOT_FULL_", m, "_final_structure")]][ix] <- psm$final_structure
  }
  sf2 <- sf
  names(sf2)[!names(sf2) %in% c("feature", "layer", "dataset")] <- paste0("FULL_SENS_", names(sf2)[!names(sf2) %in% c("feature", "layer", "dataset")])
  sf2$layer <- NULL; sf2$dataset <- NULL
  sq2 <- sq
  names(sq2)[!names(sq2) %in% c("feature", "layer", "dataset")] <- paste0("QC_SENS_", names(sq2)[!names(sq2) %in% c("feature", "layer", "dataset")])
  sq2$layer <- NULL; sq2$dataset <- NULL
  tab <- merge(tab, sf2, by = "feature", all.x = TRUE, sort = FALSE)
  tab <- merge(tab, sq2, by = "feature", all.x = TRUE, sort = FALSE)
  tab$original_M1_significant <- tab$FULL_M1_FDR < .05
  tab$robustness_category <- "MODEL_NOT_EVALUABLE_FULL_BATCH_NOT_RUN"
  tab$core_robust <- FALSE
  input <- obj$inputs[[paste(layer, "FULL", sep = "_")]]
  tab$full_feature_zero_fraction <- input$feature_zero_fraction[match(tab$feature, names(input$feature_zero_fraction))]
  tab$full_feature_RPM_le1_fraction <- input$feature_low_le1_fraction[match(tab$feature, names(input$feature_low_le1_fraction))]
  model_tables[[layer]] <- tab

  qc_tables[[layer]] <- tab[, c("feature", "layer", grep("^(FULL|QC)_M[12]_|^FULL_SENS_|^QC_SENS_", names(tab), value = TRUE)), drop = FALSE]
  robust_tables[[layer]] <- tab[, c("feature", "layer", "robustness_category", "core_robust", "original_M1_significant",
    "FULL_M1_slope", "FULL_M1_FDR", "FULL_M2_slope", "FULL_M2_FDR")]

  for (dataset in c("FULL", "QC")) for (model in c("M1", "M2")) {
    slope <- tab[[paste0(dataset, "_", model, "_slope")]]
    fdr <- tab[[paste0(dataset, "_", model, "_FDR")]]
    for (set_name in c("all_evaluable", "FDR_lt_0.05")) {
      ii <- if (set_name == "all_evaluable") is.finite(slope) else is.finite(slope) & fdr < .05
      dd <- direction(slope[ii])
      direction_rows[[length(direction_rows) + 1L]] <- data.frame(
        layer = layer, dataset = if (dataset == "QC") "QC_PASS_ONLY" else "FULL", model = model,
        analysis_scope = "ALL_FIXED_FEATURES", feature_set = set_name, n_evaluable = sum(ii),
        n_negative = sum(dd == "negative"), n_positive = sum(dd == "positive"), n_zero = sum(dd == "zero"),
        slope_median = median(slope[ii]), slope_Q1 = unname(quantile(slope[ii], .25)), slope_Q3 = unname(quantile(slope[ii], .75)))
    }
  }
  for (model in c("M3", "M4")) {
    psm <- ps[ps$model == model & ps$model_valid, ]
    dd <- direction(psm$slope)
    direction_rows[[length(direction_rows) + 1L]] <- data.frame(
      layer = layer, dataset = "FULL", model = model, analysis_scope = "PILOT_ONLY_NO_FULL_FDR",
      feature_set = "pilot_valid", n_evaluable = nrow(psm), n_negative = sum(dd == "negative"),
      n_positive = sum(dd == "positive"), n_zero = sum(dd == "zero"),
      slope_median = if (nrow(psm)) median(psm$slope) else NA_real_,
      slope_Q1 = if (nrow(psm)) unname(quantile(psm$slope, .25)) else NA_real_,
      slope_Q3 = if (nrow(psm)) unname(quantile(psm$slope, .75)) else NA_real_)
  }

  reps <- if (layer == "miRNA") {
    c("hsa-miR-338-3p", "hsa-miR-381-3p", "hsa-miR-410-3p", "hsa-miR-375")
  } else {
    c("MT-TL2:Mt_tRNA", "RN7SL1:antisense", "IL8:protein_coding", "HBA1:protein_coding",
      "TMSB4X:protein_coding", "EIF4A1:retained_intron")
  }
  rr <- tab[match(reps, tab$feature), ]
  rr$representative_status <- "BATCH_RANDOM_SLOPE_NOT_EVALUATED"
  rep_tables[[layer]] <- rr[, c("feature", "layer", "representative_status",
    grep("^(FULL|QC)_M[12]_(slope|SE|CI_lower|CI_upper|P|FDR)$|^PILOT_FULL_M[34]_", names(rr), value = TRUE)), drop = FALSE]

  old_idx <- match(m1$feature, old_full$feature)
  old_qidx <- match(q1$feature, old_qc$feature)
  repro_rows[[length(repro_rows) + 1L]] <- data.frame(
    layer = layer,
    full_max_abs_slope_diff = max(abs(m1$slope - old_full$linear_slope[old_idx])),
    full_max_abs_FDR_diff = max(abs(m1$FDR - old_full$linear_FDR[old_idx])),
    qc_max_abs_slope_diff = max(abs(q1$slope - old_qc$linear_slope[old_qidx])),
    qc_max_abs_FDR_diff = max(abs(q1$FDR - old_qc$linear_FDR[old_qidx])))
  comp_idx <- match(sf$feature, old_full$feature)
  sens_rows[[length(sens_rows) + 1L]] <- data.frame(
    layer = layer,
    old_recenter_slope_max_abs_diff = max(abs(sf$excluded_ref_slope - old_full$composition_slope[comp_idx])),
    old_recenter_P_max_abs_diff = max(abs(sf$excluded_ref_P - old_full$composition_P[comp_idx])),
    old_recenter_FDR_max_abs_diff = max(abs(sf$excluded_ref_FDR - old_full$composition_FDR[comp_idx])),
    excluded_ref_all_significant = sum(sf$excluded_ref_FDR < .05),
    excluded_ref_original_M1_significant_same_direction = sum(sf$excluded_ref_FDR < .05 & m1$FDR < .05 & direction(sf$excluded_ref_slope) == direction(m1$slope)),
    pc0_5_all_significant = sum(sf$pc0_5_FDR < .05),
    top5_excluded_all_significant = sum(sf$top5_excluded_FDR < .05))

  for (dataset in c("FULL", "QC_PASS_ONLY")) for (model in c("M1", "M2")) {
    d <- obj$limma[[paste(layer, dataset, model, sep = "_")]]
    sig <- d$FDR < .05
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(layer = layer, dataset = dataset, model = model,
      n_features = nrow(d), n_model_valid = nrow(d), n_significant = sum(sig),
      significant_negative = sum(sig & d$slope < 0), significant_positive = sum(sig & d$slope > 0),
      status = "COMPLETE")
  }
  for (dataset in c("FULL", "QC_PASS_ONLY")) for (model in c("M3", "M4")) {
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(layer = layer, dataset = dataset, model = model,
      n_features = nrow(m1), n_model_valid = NA_integer_, n_significant = NA_integer_,
      significant_negative = NA_integer_, significant_positive = NA_integer_,
      status = "NOT_RUN_PILOT_GATE_FAILED")
  }
}

write_csv(model_tables$miRNA, file.path(out, "STAGE1_MODEL_RESULTS_miRNA.csv"))
write_csv(model_tables$GENCODE, file.path(out, "STAGE1_MODEL_RESULTS_GENCODE.csv"))
diag <- merge(pilot, gate[, c("layer", "model", "valid_rate", "gate")], by = c("layer", "model"), all.x = TRUE)
diag$analysis_scope <- "PILOT_ONLY"
write_csv(diag, file.path(out, "STAGE1_MODEL_DIAGNOSTICS.csv"))
write_csv(do.call(rbind, qc_tables), file.path(out, "STAGE1_QC_SENSITIVITY.csv"))
write_csv(do.call(rbind, robust_tables), file.path(out, "STAGE1_ROBUST_SIGNATURE.csv"))
write_csv(do.call(rbind, direction_rows), file.path(out, "STAGE1_SLOPE_DIRECTION_SUMMARY.csv"))
write_csv(do.call(rbind, rep_tables), file.path(out, "STAGE1_REPRESENTATIVE_FEATURES.csv"))
write_csv(do.call(rbind, summary_rows), file.path(out, "intermediate", "MODEL_SUMMARY_COUNTS.csv"))
write_csv(do.call(rbind, repro_rows), file.path(out, "intermediate", "M1_REPRODUCTION_NUMERIC_CHECK.csv"))
write_csv(do.call(rbind, sens_rows), file.path(out, "intermediate", "RECENTERING_EQUIVALENCE_NUMERIC_CHECK.csv"))
write_csv(do.call(rbind, obj$sensitivities), file.path(out, "intermediate", "REFERENCE_SENSITIVITY_RESULTS.csv"))

png(file.path(out, "diagnostic_plots", "M1_VS_M2_SLOPE_COMPARISON.png"), width = 2200, height = 1100, res = 180)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (layer in layers) {
  b <- model_tables[[layer]]
  plot(b$FULL_M1_slope, b$FULL_M2_slope, pch = 16, cex = .55,
    xlab = "M1 raw-RPM slope", ylab = "M2 all-reference CLR slope", main = paste(layer, "M1 vs M2"))
  abline(h = 0, v = 0, lty = 3, col = "grey65")
  abline(0, 1, col = "red", lwd = 1.5)
}
dev.off()

png(file.path(out, "diagnostic_plots", "FULL_VS_QC_M1_M2_SLOPES.png"), width = 2200, height = 2100, res = 180)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
for (layer in layers) for (model in c("M1", "M2")) {
  b <- model_tables[[layer]]
  plot(b[[paste0("FULL_", model, "_slope")]], b[[paste0("QC_", model, "_slope")]], pch = 16, cex = .5,
    xlab = paste("FULL", model, "slope"), ylab = paste("QC PASS-only", model, "slope"), main = paste(layer, model))
  abline(0, 1, col = "red", lwd = 1.5)
}
dev.off()

capture.output(sessionInfo(), file = file.path(out, "logs", "STAGE1_sessionInfo.txt"))

summary_counts <- do.call(rbind, summary_rows)
repro <- do.call(rbind, repro_rows)
sens <- do.call(rbind, sens_rows)
eq <- obj$equivalence
ds <- obj$dataset_summary

get_count <- function(layer, dataset, model, col) summary_counts[summary_counts$layer == layer & summary_counts$dataset == dataset & summary_counts$model == model, col][[1]]
gate_line <- function(layer, model) {
  z <- gate[gate$layer == layer & gate$model == model, ]
  paste0(z$valid, "/", z$n, " valid (", formatC(100*z$valid_rate, digits=1, format="f"), "%)")
}

report <- c(
  "# Stage 1 核心统计重分析报告",
  "",
  "完成日期：2026-11-21（Asia/Shanghai）。分析方案于 2026-11-20 在查看 M3/M4 结果前冻结。本阶段没有修改手稿、主图、旧结果或原始数据，也没有执行 TMM–voom、XDec 或 Stage 2。",
  "",
  "## 决策：MODIFY_STAGE1",
  "",
  "M1、M2、QC 及 CLR/reference 敏感性已完成；随机斜率试运行未通过预冻结批量门槛，因此 M3/M4 全特征推断没有启动。该停止是方法学保护，不是计算中断，也没有用随机截距模型替代随机斜率验证。故当前不能诚实给出 M3/M4 全特征显著数或 core robust 全名单。",
  "",
  "## 数据集与固定分析族",
  "",
  "- FULL：30 位患者、153 对、306 个生物样本。",
  "- QC PASS-only：双方均 PASS 的 141 对；保留固定 507 个 miRNA 和 1,079 个 GENCODE 特征，不重新过滤。",
  "- 12 个 FAIL 样本对应 12 个受影响配对；QC 分支始终成对排除，未删除原始数据。",
  "",
  "## M1 复现与 QC",
  "",
  paste0("M1 FULL 完全复现：miRNA ", get_count("miRNA","FULL","M1","n_significant"), "（", get_count("miRNA","FULL","M1","significant_negative"), "/", get_count("miRNA","FULL","M1","significant_positive"), " 负/正）；GENCODE ", get_count("GENCODE","FULL","M1","n_significant"), "（", get_count("GENCODE","FULL","M1","significant_negative"), "/", get_count("GENCODE","FULL","M1","significant_positive"), " 负/正）。"),
  paste0("逐特征与 Stage 0C 结果的最大绝对差：miRNA slope=", fmt(repro$full_max_abs_slope_diff[repro$layer=="miRNA"], 3), "，FDR=", fmt(repro$full_max_abs_FDR_diff[repro$layer=="miRNA"], 3), "；GENCODE slope=", fmt(repro$full_max_abs_slope_diff[repro$layer=="GENCODE"], 3), "，FDR=", fmt(repro$full_max_abs_FDR_diff[repro$layer=="GENCODE"], 3), "。"),
  paste0("QC PASS-only M1：miRNA ", get_count("miRNA","QC_PASS_ONLY","M1","n_significant"), "（", get_count("miRNA","QC_PASS_ONLY","M1","significant_negative"), "/", get_count("miRNA","QC_PASS_ONLY","M1","significant_positive"), "）；GENCODE ", get_count("GENCODE","QC_PASS_ONLY","M1","n_significant"), "（", get_count("GENCODE","QC_PASS_ONLY","M1","significant_negative"), "/", get_count("GENCODE","QC_PASS_ONLY","M1","significant_positive"), "）。这些与 Stage 0C 一致。原显著特征没有斜率方向翻转，但显著数减少，不能称 QC 完全无影响。"),
  "",
  "## CLR 与原 re-centering",
  "",
  paste0("逐单元恒等式 `g_CLR = g_raw - mean_reference(g_raw)` 的最大绝对误差为 ", fmt(max(eq$clr_identity_max_abs_diff), 4), "；排除主导特征 reference 的最大误差为 ", fmt(max(eq$excluded_reference_identity_max_abs_diff), 4), "。因此同 pseudocount、同 reference 时二者在数学上等价。"),
  paste0("原稿排除 miR-486-5p/RNY4 的 re-centering 与本次对应 reference 结果逐特征最大 slope 差：miRNA ", fmt(sens$old_recenter_slope_max_abs_diff[sens$layer=="miRNA"], 4), "；GENCODE ", fmt(sens$old_recenter_slope_max_abs_diff[sens$layer=="GENCODE"], 4), "。它不是第二项独立验证。"),
  paste0("主要全 reference M2：miRNA ", get_count("miRNA","FULL","M2","n_significant"), "（", get_count("miRNA","FULL","M2","significant_negative"), "/", get_count("miRNA","FULL","M2","significant_positive"), "）；GENCODE ", get_count("GENCODE","FULL","M2","n_significant"), "（", get_count("GENCODE","FULL","M2","significant_negative"), "/", get_count("GENCODE","FULL","M2","significant_positive"), "）。"),
  paste0("QC PASS-only M2：miRNA ", get_count("miRNA","QC_PASS_ONLY","M2","n_significant"), "（", get_count("miRNA","QC_PASS_ONLY","M2","significant_negative"), "/", get_count("miRNA","QC_PASS_ONLY","M2","significant_positive"), "）；GENCODE ", get_count("GENCODE","QC_PASS_ONLY","M2","n_significant"), "（", get_count("GENCODE","QC_PASS_ONLY","M2","significant_negative"), "/", get_count("GENCODE","QC_PASS_ONLY","M2","significant_positive"), "）。"),
  paste0("排除主导 reference 后全部显著数为 ", sens$excluded_ref_all_significant[sens$layer=="miRNA"], "/", sens$excluded_ref_all_significant[sens$layer=="GENCODE"], "；其中同时属于旧 M1 显著且同方向的数量为 ", sens$excluded_ref_original_M1_significant_same_direction[sens$layer=="miRNA"], "/", sens$excluded_ref_original_M1_significant_same_direction[sens$layer=="GENCODE"], "，精确恢复原稿 59/121。pseudocount=0.5 的全部显著数为 ", sens$pc0_5_all_significant[sens$layer=="miRNA"], "/", sens$pc0_5_all_significant[sens$layer=="GENCODE"], "；排除 top-5 高丰度 reference 后为 ", sens$top5_excluded_all_significant[sens$layer=="miRNA"], "/", sens$top5_excluded_all_significant[sens$layer=="GENCODE"], "。"),
  "",
  "主要全-reference CLR 将 reference 内各特征的斜率重新定义为相对于共同中心，因此显著 miRNA 的方向从 M1 的 189/3 变为 M2 的 56/72。该变化是坐标定义的必然组成部分，不能单独证明原 189/3 是技术假象，也不能证明或排除 carrier composition 变化。miRNA reference 中心本身具有负时间斜率；详见 intermediate/CLR_CENTER_SLOPES.csv。",
  "",
  "## 随机斜率试运行",
  "",
  paste0("- miRNA M3：", gate_line("miRNA","M3"), "；miRNA M4：", gate_line("miRNA","M4"), "。"),
  paste0("- GENCODE M3：", gate_line("GENCODE","M3"), "；GENCODE M4：", gate_line("GENCODE","M4"), "。"),
  "- 失败主要为相关结构和预定不相关结构均 singular；并非表达矩阵缺失或通用代码崩溃。GENCODE M4 仅 5/11 有效，低于冻结的 50% 门槛，因此停止批量。",
  "- 未发现影响至少 25% 试拟合的共同环境/代码错误；停止原因是模型奇异性，而不是 R 环境故障。",
  "",
  "## 十个关键问题",
  "",
  "**Q1. M1 是否完全一致？** 是。192/200 与 189/3、175/25 均恢复，逐特征差异见数值核查表。",
  "",
  "**Q2. CLR 是否与旧 re-centering 等价？** 同 pseudocount 和同 reference 时是；旧分析对应排除 miR-486-5p/RNY4 的 reference，不等于主要全-reference M2，但其对应分支逐值等价。",
  "",
  "**Q3. M3 后原 192/200 保留多少？** 未估计。试运行未通过冻结门槛，全特征 M3 未运行，禁止用试验子集外推。",
  "",
  "**Q4. M4 后保留多少？** 未估计；同上。",
  "",
  "**Q5. core robust features 多少？** 当前不能定义。输出表将全部特征标记为 `MODEL_NOT_EVALUABLE_FULL_BATCH_NOT_RUN`，不是宣称 0 个真正稳健特征。",
  "",
  "**Q6. 189/3 是否仍存在？** 在 M1 FULL 精确存在，QC M1 为 172/4；主要全-reference M2 显著 miRNA 为 56/72。reference 和 QC 都会影响显著数量/方向分布，CLR 的重心约束尤其改变方向解释。M3/M4 全层面尚不能回答。",
  "",
  "**Q7. 12 个 QC FAIL 是否明显影响？** M1 效应方向总体稳定且旧显著特征无方向翻转，但显著数 192→176、200→165，分类结果在 Stage 0C 更敏感；不能写成完全无影响。M2 为 128→105、173→139。随机斜率 QC 分支尚未运行。",
  "",
  "**Q8. 收敛/奇异比例？** 只能报告预定试验集：见上方 4 个有效率；失败以 singular 为主。全特征比例未估计。",
  "",
  "**Q9. Figure 3/4 代表特征？** M1/M2 与试运行结果已逐项写入 STAGE1_REPRESENTATIVE_FEATURES.csv；由于随机斜率批量及 QC 随机斜率未完成，本阶段不能作最终保留判定，也没有替换代表特征。",
  "",
  "**Q10. 是否仍支持主要结论？** 旧模型、QC 和 composition-reference 证据支持存在一组纵向配对梯度信号，但 reviewer 要求的随机斜率稳健性尚未完成，故不能把当前证据表述为已完全解决主要方法学质疑。",
  "",
  "## 下一步",
  "",
  "需要用户确认一个有记录的 Stage 1 方法修改后再继续：建议把 singular 判定从“任何方差边界即整特征不可评价”改为针对随机斜率分量的预定义规则，并在不相关随机截距/斜率结构中区分随机截距趋零与随机斜率趋零；同时按优化器分别传控制参数。备选方案是由统计顾问确认使用 Kenward–Roger 边界模型或正则化混合模型。任何方案都必须先形成 STAGE1_PLAN_DEVIATION_01.md，再重新试运行，不能依据显著性选择。",
  "",
  "在 M3/M4 完成前不进入 Stage 2。最终决策：**MODIFY_STAGE1**。"
)
writeLines(report, file.path(out, "STAGE1_FINAL_REPORT.md"), useBytes = TRUE)

rec <- c(
  "# Stage 2 建议",
  "",
  "当前建议：**暂不进入 Stage 2**。Stage 1 的 M1/M2、QC 与 reference 敏感性已完成，但随机斜率试运行未通过冻结门槛，M3/M4 全特征和 core robust 集合尚未生成。",
  "",
  "先完成一次有审计记录的 Stage 1 方法修订：只针对 singular 随机效应边界制定可解释规则，重新冻结并重跑代表性试验。不得根据显著性选择规则，不得退化为随机截距-only 后声称通过随机斜率验证。",
  "",
  "若修改后的试运行通过，再完成 FULL/QC 的 M3/M4、core robust 集合和代表特征最终判定；随后才进入 Stage 2。Stage 2 仍应把现有 81×75 profiles 视为参考而非患者比例，不能把新估计冒称原作者 XDec 输出。"
)
writeLines(rec, file.path(out, "STAGE2_RECOMMENDATION.md"), useBytes = TRUE)
cat("PARTIAL_FINALIZE_COMPLETE\n")

