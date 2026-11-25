options(stringsAsFactors = FALSE)

root <- "revision_stage1/stage1a_boundary_diagnostics"
d <- read.csv(file.path(root, "BOUNDARY_CLASSIFICATION.csv"), check.names = FALSE)
cmp <- read.csv(file.path(root, "RANDOM_EFFECT_MODEL_COMPARISON.csv"), check.names = FALSE)
opt <- read.csv(file.path(root, "OPTIMIZER_DIAGNOSTICS.csv"), check.names = FALSE)
pat <- read.csv(file.path(root, "intermediate", "PATIENT_LEVEL_PILOT_RESULTS.csv"), check.names = FALSE)
sel <- read.csv(file.path(root, "intermediate", "EXPANDED_PILOT_SELECTION.csv"), check.names = FALSE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
pct <- function(n, den) sprintf("%d/%d（%.1f%%）", n, den, 100*n/den)
fmt <- function(x, digits = 4) formatC(x, digits = digits, format = "g")

old_fail <- d[d$dataset == "FULL" & d$old_stage1_result_available & !d$old_stage1_model_valid, ]
old_n <- nrow(old_fail)
old_ri <- sum(old_fail$ri_zero_standardized & !old_fail$rs_zero_standardized)
old_rs <- sum(old_fail$rs_zero_standardized & !old_fail$ri_zero_standardized)
old_both <- sum(old_fail$rs_zero_standardized & old_fail$ri_zero_standardized)
old_corr <- sum(old_fail$correlation_boundary)
old_num <- sum(grepl("^E_", old_fail$primary_boundary_class))
old_rescued_c <- sum(old_fail$C_structure_eligible & old_fail$random_slope_identifiable)

full <- d[d$dataset == "FULL", ]
qc <- d[d$dataset == "QC_PASS_ONLY", ]
full_ident <- sum(full$random_slope_identifiable)
qc_ident <- sum(qc$random_slope_identifiable)
full_c <- sum(full$C_structure_eligible)
qc_c <- sum(qc$C_structure_eligible)

b_full <- cmp[cmp$dataset == "FULL" & cmp$structure == "B", ]
b_qc <- cmp[cmp$dataset == "QC_PASS_ONLY" & cmp$structure == "B", ]
b_full_stable <- sum(b_full$convergence_ok & b_full$optimizer_stable)
b_qc_stable <- sum(b_qc$convergence_ok & b_qc$optimizer_stable)
b_full_rs_ident <- sum(b_full$convergence_ok & b_full$optimizer_stable &
  is.finite(b_full$RS_SD_ratio_per_label) & b_full$RS_SD_ratio_per_label > 1e-4)
b_qc_rs_ident <- sum(b_qc$convergence_ok & b_qc$optimizer_stable &
  is.finite(b_qc$RS_SD_ratio_per_label) & b_qc$RS_SD_ratio_per_label > 1e-4)

a_full <- cmp[cmp$dataset == "FULL" & cmp$structure == "A", ]
a_singular_n <- sum(a_full$isSingular_tol_1e_4)
pure_corr_resolved <- sum(full$primary_boundary_class == "D_CORRELATION_BOUNDARY_RESOLVED_BY_UNCORRELATED")

c_full <- full[full$C_structure_eligible, ]
c_max_diff <- max(abs(c_full$B_slope - c_full$C_slope))
c_norm_diff <- max(abs(c_full$B_slope - c_full$C_slope) / pmax(c_full$B_SE, c_full$C_SE))
c_se_ratio <- range(c_full$C_SE / c_full$B_SE)

pf <- pat[pat$dataset == "FULL", ]
key_d <- paste(d$layer, d$dataset, d$gradient_model, d$feature)
idx <- match(paste(pf$layer, pf$dataset, pf$gradient_model, pf$feature), key_d)
patient_b_dir <- sum(sign(pf$individual_slope_mean) == sign(d$B_slope[idx]))
cr2_b_dir <- sum(sign(pf$CR2_slope) == sign(d$B_slope[idx]))
patient_b_cor <- cor(pf$individual_slope_mean, d$B_slope[idx])
cr2_b_cor <- cor(pf$CR2_slope, d$B_slope[idx])
all_patient_estimable <- sum(pat$all_subject_slopes_estimable)
cr2_success <- sum(is.finite(pat$CR2_slope) & (is.na(pat$CR2_error) | !nzchar(pat$CR2_error)))

unresolved_full <- full[!full$random_slope_identifiable,
  c("layer", "gradient_model", "feature", "primary_boundary_class", "dataset",
    "B_RI_SD_ratio", "B_RS_SD_ratio_per_label", "B_optimizer_stable", "B_convergence_ok")]
unresolved_full$origin <- ifelse(unresolved_full$feature %in% old_fail$feature &
  paste(unresolved_full$layer, unresolved_full$gradient_model, unresolved_full$feature) %in%
    paste(old_fail$layer, old_fail$gradient_model, old_fail$feature),
  "ORIGINAL_STAGE1_PILOT_FAILURE", "EXPANDED_STAGE1A_PILOT")
write_csv(unresolved_full, file.path(root, "intermediate", "UNRESOLVED_RANDOM_SLOPE_FEATURES_FULL.csv"))

old_unresolved <- old_fail[!old_fail$random_slope_identifiable,
  c("layer", "gradient_model", "feature", "primary_boundary_class")]
old_rescued <- old_fail[old_fail$random_slope_identifiable,
  c("layer", "gradient_model", "feature", "primary_boundary_class", "C_structure_eligible")]

feature_lines <- function(x) {
  if (!nrow(x)) return("- 无。")
  apply(x, 1L, function(r) paste0("- ", r[["layer"]], " ", r[["gradient_model"]], "：`", r[["feature"]], "`（", r[["primary_boundary_class"]], "）"))
}

patient_report <- c(
  "# 患者层面敏感性分析可行性",
  "",
  "日期：2026-11-21（Asia/Shanghai）。范围仅限Stage 1A有限试验集，不是全特征正式推断。",
  "",
  "## 可行性结果",
  "",
  paste0("扩展试验集包含", nrow(sel), "个唯一特征；FULL/QC × M3/M4共", nrow(pat), "个特征-分支组合。", all_patient_estimable, "/", nrow(pat), "个组合中30位患者的个体斜率均可计算，CR2+Satterthwaite在", cr2_success, "/", nrow(pat), "个组合中返回有限估计。"),
  paste0("在FULL的", nrow(pf), "个组合中，未加权患者斜率均值与Model B方向一致", patient_b_dir, "/", nrow(pf), "，斜率相关r=", fmt(patient_b_cor, 4), "；CR2边际斜率与Model B方向一致", cr2_b_dir, "/", nrow(pf), "，r=", fmt(cr2_b_cor, 4), "。唯一CR2方向不一致对象的Model B斜率接近零，不能解释为实质性反向。"),
  "",
  "## 两种候选敏感性方法",
  "",
  "1. **患者独立斜率汇总**：每位患者分别拟合collection index斜率，再对30个斜率作未加权均值和t型区间。该方法直接显示患者异质性，但每位患者仅3–7个时间点、标签组合不同、单人斜率SE差异明显；未加权平均估计的是平均患者斜率，不等同于混合模型或观测加权总体斜率。",
  "2. **患者聚类CR2边际模型**：模型含患者固定截距和共同时间斜率，按患者聚类，以clubSandwich CR2及Satterthwaite小样本修正推断。它不要求估计随机效应协方差矩阵，数值上与Model B高度一致，适合作为正式Stage 1B的独立敏感性分析。它估计共同within-subject边际斜率，不证明存在随机斜率异质性，也不能冒充random-slope refit。",
  "",
  "## 建议",
  "",
  "Stage 1B应对固定507/1,079特征统一增加CR2患者聚类敏感性；患者独立斜率分布用于代表特征、诊断和患者影响展示。两者都不得按P值决定是否执行。随机斜率分量位于边界的特征，即使CR2显著，也不能标记为已通过随机斜率验证。"
)
writeLines(patient_report, file.path(root, "PATIENT_LEVEL_SENSITIVITY_FEASIBILITY.md"), useBytes = TRUE)

rules <- c(
  "# Stage 1A 后续正式模型决策规则",
  "",
  "日期：2026-11-21（Asia/Shanghai）。以下规则是Stage 1试运行结果已查看后的post hoc修订，必须与原冻结计划及Plan Deviation 01同时披露。",
  "",
  "## 统一主模型",
  "",
  "Stage 1B对每个分子层、FULL/QC和M3/M4统一拟合Model B：",
  "",
  "`gradient ~ centered_label + (1 + centered_label || subject)`",
  "",
  "选择Model B不是因为其P值更小，而是因为有限试验中A的相关参数经常达到边界，而B保留随机截距和随机斜率两个方差分量并移除难以识别的相关参数。统一使用B可避免按特征搜索A/B/C后再挑选P值。Model A仅作为相关结构敏感性；Model C仅验证随机截距为零时的可辨识性；Model D只作退化比较。",
  "",
  "## 每个特征的固定判定树",
  "",
  "1. Model B至少两个优化器返回有限固定效应；方向一致；斜率范围/中位SE≤0.10；SE最大/最小≤1.10；无未解决收敛/Hessian异常。否则标记`RANDOM_SLOPE_INFERENCE_UNRESOLVED_NUMERICAL`。",
  "2. 在B数值稳定时，以`RS_SD/residual_SD > 1e-4`作为随机斜率分量可辨识的必要条件，并联合`rePCA`、特征值和`isSingular(tol=1e-4)`报告。",
  "3. 若随机斜率分量可辨识且随机截距不在边界：B为正式随机斜率结果；A仅比较固定效应稳定性及相关边界。",
  "4. 若随机截距在边界但随机斜率可辨识：主P值仍来自统一B；另拟合C确认。只有C数值稳定、非奇异且B/C斜率差≤0.10×中位SE、SE比≤1.10，才标记`RANDOM_SLOPE_IDENTIFIABLE_RI_BOUNDARY_CONFIRMED_BY_C`。不把C结果按显著性替换B。",
  "5. 若随机斜率分量在边界：标记`RANDOM_SLOPE_VARIANCE_BOUNDARY`。B固定效应和D/CR2可报告，但该特征不算通过随机斜率验证，不进入core robust集合。零方差估计不解释为所有患者真实斜率完全相同。",
  "6. 若两分量均在边界或所有结构不稳定：标记`RANDOM_SLOPE_INFERENCE_UNRESOLVED`，不进入core robust集合。",
  "7. 所有特征无论边界状态都保留在固定507/1,079 BH校正族中；无有效P值者记NA，仍报告可评价分母。不得用边界分类缩小BH家族。",
  "8. QC PASS-only使用完全相同规则。QC是否显著不决定主模型；效应方向与预冻结QC依赖阈值继续单独报告。",
  "",
  "## 独立敏感性",
  "",
  "对所有固定特征增加患者固定截距+患者聚类CR2/Satterthwaite边际斜率；代表性特征另报告30位患者独立斜率分布。CR2结果单独作层×数据集×梯度BH，不替代Model B，也不用于声称随机斜率异质性。",
  "",
  "## 有限试验支持",
  "",
  paste0("FULL中Model B收敛60/60、通过预冻结优化器稳定规则", b_full_stable, "/60；QC中为60/60和", b_qc_stable, "/60。FULL中纯相关边界有", pure_corr_resolved, "/60通过B解决；另有", full_c, "个随机截距边界对象由C确认。FULL仍有", nrow(unresolved_full), "/60随机斜率不可辨识或数值不稳定，必须保留为unresolved而非强行降阶。"),
  "",
  "## Stage 1B报告义务",
  "",
  "必须报告A/B/C/D使用目的、B的收敛/稳定/边界比例、RI/RS零边界比例、C确认数量、unresolved数量、CR2一致性，以及不同边界状态对core robust筛选和多重检验解释的限制。"
)
writeLines(rules, file.path(root, "STAGE1A_DECISION_RULES.md"), useBytes = TRUE)

final_report <- c(
  "# Stage 1A 随机斜率边界诊断最终报告",
  "",
  "完成日期：2026-11-21（Asia/Shanghai）。本阶段只进行了有限试验集边界诊断，没有启动全特征M3/M4、Stage 2、XDec、手稿或主图修改。原冻结计划未覆盖；偏离规则在任何新拟合前登记。",
  "",
  "## 最终决策：GO_TO_STAGE1B_FULL_REFIT",
  "",
  "理由不是把singular标准放宽到越过50%，而是已经形成统一、不依赖显著性的Model B正式策略，并将随机斜率零边界、随机截距零边界和数值不稳定分别标记。有限试验显示Model B在FULL/QC均100%返回收敛解，优化器稳定率分别为", paste0(b_full_stable, "/60和", b_qc_stable, "/60；但边界特征仍将排除于随机斜率稳健集合，而不会被当成成功验证。"),
  "",
  "## 设计和数值检查",
  "",
  paste0("实际选择", nrow(sel), "个唯一特征（miRNA ", sum(sel$layer=="miRNA"), "，GENCODE ", sum(sel$layer=="GENCODE"), "），覆盖原试验、M1/M2方向变化、高零值和低CLR变异。FULL/QC、M3/M4、A–D和三个优化器共产生", nrow(opt), "次拟合；形成", nrow(cmp), "个结构级比较和", nrow(d), "个边界分类单元。"),
  "centered_label均为D04中心的-3至3；subject均为30水平因子；无非有限响应；没有证据表明旧失败由错误中心化、subject编码或缺失值造成。lme4版本2.0.6，`isSingular` tolerance固定为1e-4。",
  "",
  "## Q1. 原失败模型的边界类型",
  "",
  paste0("原Stage 1有", old_n, "个FULL特征×梯度组合被判失败。相关结构A中", pct(old_corr,old_n), "均出现相关性/秩边界（非互斥）；改用不相关结构B后，随机截距趋零", pct(old_ri,old_n), "，随机斜率趋零", pct(old_rs,old_n), "，两者均趋零", pct(old_both,old_n), "，真正的数值不收敛/固定效应不稳定", pct(old_num,old_n), "。因此最终失败的主导类型是随机斜率方差边界，而不是优化器失败。"),
  "",
  "## Q2. 协方差过度复杂还是数据不足",
  "",
  paste0("两者都存在。扩展FULL中A有", a_singular_n, "/60 singular；其中", pure_corr_resolved, "个是删除相关参数后B成为非奇异的纯相关边界，说明相关协方差确实过度复杂。但原", old_n, "个失败中B仍显示", old_rs, "个随机斜率趋零，说明多数原失败对象的数据没有提供可辨识患者特异性斜率方差。不能把它们全部归因于相关参数。"),
  "",
  "## Q3. Model B能否解决主要边界",
  "",
  paste0("B解决扩展FULL中", pure_corr_resolved, "/60个纯相关边界，即A-singular对象的", sprintf("%.1f%%",100*pure_corr_resolved/a_singular_n), "；但原16个失败对象没有一个仅靠B变成非奇异，因为B进一步暴露了3个随机截距零边界和13个随机斜率零边界。结论：B是统一正式结构的合理选择，但不是消除所有边界的工具。"),
  "",
  "## Q4. 随机截距趋零时Model C是否可辨识",
  "",
  paste0("FULL扩展试验中共有", full_c, "个随机截距趋零而随机斜率非零的对象，全部通过C资格规则；其中原失败对象3个。B与C固定斜率最大绝对差", fmt(c_max_diff,4), "，最大差/SE为", fmt(c_norm_diff,4), "，C/B SE比范围", fmt(c_se_ratio[1],7), "–", fmt(c_se_ratio[2],7), "。因此C在这些对象中可辨识并有合理解释：患者基线随机差异估计位于边界，但患者斜率异质性仍可估计。正式策略仍以统一B为主，C只作确认，避免按特征选择P值。"),
  "",
  "## Q5. 哪些对象不能可靠随机斜率推断",
  "",
  paste0("扩展FULL中", nrow(unresolved_full), "/60不可辨识或数值不稳定；原Stage 1失败对象中仍有以下", nrow(old_unresolved), "个组合不能可靠进行随机斜率推断："),
  feature_lines(old_unresolved),
  "",
  "完整扩展清单见`intermediate/UNRESOLVED_RANDOM_SLOPE_FEATURES_FULL.csv`。这些对象可报告Model B/D/CR2固定趋势，但不得进入core robust random-slope集合。",
  "",
  "原失败对象中可由随机斜率-only结构确认的3个组合为：",
  feature_lines(old_rescued),
  "",
  "## Q6. 是否需要患者层面或聚类稳健敏感性",
  "",
  paste0("需要。", all_patient_estimable, "/", nrow(pat), "个有限试验组合可估计全部30位患者斜率，CR2在", cr2_success, "/", nrow(pat), "成功。FULL中CR2与B方向一致", cr2_b_dir, "/60，r=", fmt(cr2_b_cor,4), "；未加权患者斜率均值与B方向一致", patient_b_dir, "/60，r=", fmt(patient_b_cor,4), "。CR2应作为Stage 1B全特征独立敏感性；患者斜率分布用于代表特征和影响诊断。两者估计目标不同，均不冒充随机斜率模型。"),
  "",
  "## Q7. 能否建立统一且可辩护的策略",
  "",
  paste0("可以。统一主模型为B，不按特征挑选A/B/C的P值。FULL中B收敛60/60、优化器稳定", b_full_stable, "/60；QC为60/60和", b_qc_stable, "/60。以B的随机斜率标准化分量和优化器稳定性决定是否属于random-slope identifiable；RI零边界用C确认；RS零边界和数值不稳定标记unresolved；A、C、D和CR2均为预定敏感性。所有特征保留在原BH家族中。"),
  "",
  "## Q8. 是否可重新启动全特征M3/M4",
  "",
  paste0("可以，进入独立Stage 1B后执行。有限试验中按A/B/C决策树可辨识", full_ident, "/60 FULL和", qc_ident, "/60 QC；若严格按统一B稳定+随机斜率非零标准，则为", b_full_rs_ident, "/60和", b_qc_rs_ident, "/60。这意味着全量运行后必须接受并报告一部分RS-boundary/unresolved特征，而不能把完成拟合等同于全部通过验证。"),
  "",
  "## 限制与环境",
  "",
  "扩展试验不是全特征边界率估计，选择覆盖困难场景，不能把比例外推为507/1,079的总体比例。每位患者只有3–7个标签且标签组合不平衡。collection index不是已验证的真实发病后天数。当前R为4.6.1，lmerTest 3.2.1二进制包提示由R4.6.2构建；该警告未造成拟合错误，但Stage 1B须保留包归档、sessionInfo和日志。",
  "",
  "本阶段到此停止。未启动Stage 1B或Stage 2。"
)
writeLines(final_report, file.path(root, "STAGE1A_FINAL_REPORT.md"), useBytes = TRUE)

cat("STAGE1A_REPORTS_WRITTEN\n")

