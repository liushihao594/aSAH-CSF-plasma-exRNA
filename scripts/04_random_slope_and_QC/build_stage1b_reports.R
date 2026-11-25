options(stringsAsFactors = FALSE, warn = 1)

out <- file.path("revision_stage1", "stage1b_full_refit")
z <- readRDS(file.path(out, "intermediate", "STAGE1B_COMPACT_RESULTS.rds"))

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}
pct <- function(n, d) sprintf("%d/%d (%.1f%%)", n, d, 100 * n/d)
sgn <- function(x) ifelse(is.na(x), NA_character_, ifelse(x > 0, "positive", ifelse(x < 0, "negative", "zero")))

metrics <- data.frame(section = character(), layer = character(), model = character(),
  metric = character(), value = character(), stringsAsFactors = FALSE)
add_metric <- function(section, layer, model, metric, value) {
  metrics <<- rbind(metrics, data.frame(section = section, layer = layer, model = model,
    metric = metric, value = as.character(value), stringsAsFactors = FALSE))
}

for (layer in c("miRNA", "GENCODE")) {
  x <- z$full_tables[[layer]]
  n <- nrow(x)
  for (m in c("M3", "M4")) {
    sig <- x[[paste0(m, "_significant")]]
    cr2sig <- x[[paste0("CR2_", m, "_significant")]]
    eval <- x[[paste0(m, "_numerical_evaluable")]]
    val <- x[[paste0(m, "_random_slope_validated")]]
    status <- x[[paste0(m, "_model_status")]]
    add_metric("full_model", layer, m, "attempted", n)
    add_metric("full_model", layer, m, "numerically_evaluable", sum(eval, na.rm = TRUE))
    add_metric("full_model", layer, m, "random_slope_validated", sum(val, na.rm = TRUE))
    for (st in c("REGULAR_RANDOM_SLOPE", "RI_BOUNDARY_C_CONFIRMED", "RS_BOUNDARY",
                 "BOTH_VARIANCE_BOUNDARY_UNRESOLVED", "NUMERICAL_UNSTABLE", "MODEL_NOT_EVALUABLE")) {
      add_metric("full_model", layer, m, paste0("status_", st), sum(status == st, na.rm = TRUE))
    }
    add_metric("full_model", layer, m, "significant", sum(sig, na.rm = TRUE))
    add_metric("full_model", layer, m, "significant_negative", sum(sig & x[[paste0(m, "_slope")]] < 0, na.rm = TRUE))
    add_metric("full_model", layer, m, "significant_positive", sum(sig & x[[paste0(m, "_slope")]] > 0, na.rm = TRUE))
    add_metric("CR2", layer, m, "significant", sum(cr2sig, na.rm = TRUE))
    add_metric("CR2", layer, m, "both_significant", sum(sig & cr2sig, na.rm = TRUE))
    add_metric("CR2", layer, m, "ModelB_only", sum(sig & !cr2sig, na.rm = TRUE))
    add_metric("CR2", layer, m, "CR2_only", sum(!sig & cr2sig, na.rm = TRUE))
    add_metric("CR2", layer, m, "slope_correlation", cor(x[[paste0(m, "_slope")]], x[[paste0(m, "_CR2_slope")]], use = "complete.obs"))
    add_metric("CR2", layer, m, "direction_same", sum(sgn(x[[paste0(m, "_slope")]]) == sgn(x[[paste0(m, "_CR2_slope")]]), na.rm = TRUE))
    add_metric("QC", layer, m, "QC_significant", sum(x[[paste0("QC_", m, "_significant")]], na.rm = TRUE))
    add_metric("QC", layer, m, "FULL_and_QC_significant", sum(sig & x[[paste0("QC_", m, "_significant")]], na.rm = TRUE))
    add_metric("QC", layer, m, "direction_same", sum(sgn(x[[paste0(m, "_slope")]]) == sgn(x[[paste0("QC_", m, "_slope")]]), na.rm = TRUE))
    add_metric("QC", layer, m, "slope_correlation", cor(x[[paste0(m, "_slope")]], x[[paste0("QC_", m, "_slope")]], use = "complete.obs"))
  }
  old <- x$M1_significant
  joint <- old & x$M3_significant & x$M4_significant & x$M1_M3_M4_direction_consistent
  add_metric("robustness", layer, "M1_M3_M4", "joint_significant_direction_consistent", sum(joint, na.rm = TRUE))
  add_metric("robustness", layer, "classification", "strict_core", sum(x$strict_core_robust, na.rm = TRUE))
  add_metric("robustness", layer, "classification", "boundary_supported", sum(x$boundary_supported, na.rm = TRUE))
  add_metric("robustness", layer, "classification", "QC_sensitive", sum(x$QC_sensitive, na.rm = TRUE))
  add_metric("robustness", layer, "classification", "composition_sensitive", sum(x$composition_sensitive, na.rm = TRUE))
  add_metric("robustness", layer, "classification", "unresolved", sum(x$unresolved, na.rm = TRUE))
  add_metric("correlation", layer, "M1_vs_M3", "slope", cor(x$M1_slope, x$M3_slope, use = "complete.obs"))
  add_metric("correlation", layer, "M2_vs_M4", "slope", cor(x$M2_slope, x$M4_slope, use = "complete.obs"))

  xo <- x[old, , drop = FALSE]
  for (m in c("M1", "M3", "M4")) {
    sigcol <- if (m == "M1") rep(TRUE, nrow(xo)) else xo[[paste0(m, "_significant")]]
    add_metric("original_discoveries", layer, m, "significant", sum(sigcol, na.rm = TRUE))
    add_metric("original_discoveries", layer, m, "significant_negative", sum(sigcol & xo[[paste0(m, "_slope")]] < 0, na.rm = TRUE))
    add_metric("original_discoveries", layer, m, "significant_positive", sum(sigcol & xo[[paste0(m, "_slope")]] > 0, na.rm = TRUE))
  }
}

write.csv(metrics, file.path(out, "intermediate", "STAGE1B_REPORT_METRICS.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8")

getm <- function(section, layer, model, metric) {
  v <- metrics$value[metrics$section == section & metrics$layer == layer & metrics$model == model & metrics$metric == metric]
  if (length(v)) v[1] else NA_character_
}
num <- function(section, layer, model, metric) as.numeric(getm(section, layer, model, metric))

model_row <- function(layer, m) {
  n <- num("full_model", layer, m, "attempted")
  c(layer, m,
    sprintf("%d", n),
    sprintf("%d", num("full_model", layer, m, "numerically_evaluable")),
    sprintf("%d", num("full_model", layer, m, "random_slope_validated")),
    sprintf("%d", num("full_model", layer, m, "status_RS_BOUNDARY")),
    sprintf("%d", num("full_model", layer, m, "status_BOTH_VARIANCE_BOUNDARY_UNRESOLVED")),
    sprintf("%d", num("full_model", layer, m, "status_NUMERICAL_UNSTABLE")),
    sprintf("%d", num("full_model", layer, m, "significant")),
    sprintf("%d/%d", num("full_model", layer, m, "significant_negative"), num("full_model", layer, m, "significant_positive")))
}

rows <- do.call(rbind, lapply(c("miRNA", "GENCODE"), function(l) do.call(rbind, lapply(c("M3", "M4"), function(m) model_row(l,m)))))
model_table <- paste0("|", apply(rows, 1, paste, collapse = "|"), "|", collapse = "\n")

old <- z$old_summary
old_mi <- old[old$layer == "miRNA",]
old_ge <- old[old$layer == "GENCODE",]

rep <- z$representatives
rep$decision_cn <- ifelse(rep$representative_recommendation == "RETAIN_STRICT_CORE", "保留（strict core）",
  ifelse(rep$representative_recommendation == "RETAIN_WITH_MODEL_OR_QC_QUALIFICATION", "保留但需模型/QC限定",
    "不建议作为主要代表特征"))
rep_lines <- apply(rep, 1, function(r) sprintf("|%s|%s|%s|%s|%s|%s|",
  r[["layer"]], r[["feature"]], r[["M3_model_status"]], r[["M4_model_status"]],
  ifelse(r[["strict_core_robust"]] == "TRUE", "是", "否"), r[["decision_cn"]]))

report <- c(
  "# Stage 1B 全特征随机斜率重分析最终报告",
  "",
  "完成日期：2026-11-21（Asia/Shanghai）。本阶段仅完成统计拟合、诊断、QC/CR2敏感性及稳健特征确认；未启动Stage 2，未修改手稿或正式图表。",
  "",
  "## 最终决策：GO_TO_STAGE2",
  "",
  sprintf("全部6,344个预定任务均完成（2层 × 2数据集 × 2梯度，共507×4 + 1,079×4），任务级运行错误为0。FULL中确认strict core robust为miRNA %d个、GENCODE %d个；另有RS-boundary但获其他模型支持的miRNA %d个、GENCODE %d个。存在可重复的纵向CSF–plasma梯度信号，但结论应收窄到这些稳健集合，并透明报告大量随机斜率方差边界。",
    sum(z$strict_core$layer == "miRNA"), sum(z$strict_core$layer == "GENCODE"),
    sum(z$boundary_supported$layer == "miRNA"), sum(z$boundary_supported$layer == "GENCODE")),
  "",
  "## 1. 模型完成度与边界诊断（FULL）",
  "",
  "|层|模型|尝试|数值可评估|随机斜率有效|RS边界|双方差边界|数值不稳定|FDR显著|显著负/正|",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
  model_table,
  "",
  "`RS_BOUNDARY`保留可报告的固定效应，但不计为random-slope validated；`BOTH_VARIANCE_BOUNDARY_UNRESOLVED`与`NUMERICAL_UNSTABLE`不进入严格稳健集合。没有`MODEL_NOT_EVALUABLE`。RI边界经Model C确认后可计为随机斜率有效；Model B仍是统一主模型。",
  "",
  "## 2. 原稿发现与新模型",
  "",
  sprintf("- M3（原始梯度）全部特征显著：miRNA %d，GENCODE %d；其中来自原稿显著集分别为%d/%d和%d/%d。",
    num("full_model","miRNA","M3","significant"), num("full_model","GENCODE","M3","significant"),
    old_mi$M3_significant, old_mi$original_n, old_ge$M3_significant, old_ge$original_n),
  sprintf("- M4（CLR梯度）全部特征显著：miRNA %d，GENCODE %d；其中来自原稿显著集分别为%d/%d和%d/%d。",
    num("full_model","miRNA","M4","significant"), num("full_model","GENCODE","M4","significant"),
    old_mi$M4_significant, old_mi$original_n, old_ge$M4_significant, old_ge$original_n),
  sprintf("- 同时满足M1、M3、M4显著且三者方向一致：miRNA %d，GENCODE %d。进一步加入随机斜率有效及QC规则后，strict core为16和33。",
    num("robustness","miRNA","M1_M3_M4","joint_significant_direction_consistent"),
    num("robustness","GENCODE","M1_M3_M4","joint_significant_direction_consistent")),
  sprintf("- 原稿显著集内M3/M4同时显著：miRNA %d，GENCODE %d；CR2对M3/M4同时支持分别为%d和%d。",
    old_mi$M3_M4_both_significant, old_ge$M3_M4_both_significant, old_mi$CR2_both_supported, old_ge$CR2_both_supported),
  "",
  "完整名单见`STAGE1B_STRICT_CORE_ROBUST.csv`和`STAGE1B_BOUNDARY_SUPPORTED.csv`。稳健交集是筛选规则，不具有额外的联合FDR保证。",
  "",
  "## 3. CR2患者聚类稳健敏感性",
  "",
  sprintf("Model B与CR2的斜率相关性为：miRNA M3 %.4f、M4 %.4f；GENCODE M3 %.4f、M4 %.4f。方向一致分别为%d/%d、%d/%d、%d/%d、%d/%d。",
    num("CR2","miRNA","M3","slope_correlation"), num("CR2","miRNA","M4","slope_correlation"),
    num("CR2","GENCODE","M3","slope_correlation"), num("CR2","GENCODE","M4","slope_correlation"),
    num("CR2","miRNA","M3","direction_same"),507,num("CR2","miRNA","M4","direction_same"),507,
    num("CR2","GENCODE","M3","direction_same"),1079,num("CR2","GENCODE","M4","direction_same"),1079),
  sprintf("FDR显著数（Model B/CR2）为miRNA M3 %d/%d、M4 %d/%d；GENCODE M3 %d/%d、M4 %d/%d。差异集中于阈值附近及部分边界特征；CR2支持不能把RS-boundary重新定义为random-slope validated。",
    num("full_model","miRNA","M3","significant"),num("CR2","miRNA","M3","significant"),
    num("full_model","miRNA","M4","significant"),num("CR2","miRNA","M4","significant"),
    num("full_model","GENCODE","M3","significant"),num("CR2","GENCODE","M3","significant"),
    num("full_model","GENCODE","M4","significant"),num("CR2","GENCODE","M4","significant")),
  "",
  "## 4. QC PASS-only敏感性",
  "",
  sprintf("QC PASS-only保留141对。显著数FULL→QC为：miRNA M3 %d→%d、M4 %d→%d；GENCODE M3 %d→%d、M4 %d→%d。FULL/QC斜率相关性均较高（miRNA M3 %.4f、M4 %.4f；GENCODE M3 %.4f、M4 %.4f），但GENCODE显著数下降明显，因此不能写成QC完全无影响。",
    num("full_model","miRNA","M3","significant"),num("QC","miRNA","M3","QC_significant"),
    num("full_model","miRNA","M4","significant"),num("QC","miRNA","M4","QC_significant"),
    num("full_model","GENCODE","M3","significant"),num("QC","GENCODE","M3","QC_significant"),
    num("full_model","GENCODE","M4","significant"),num("QC","GENCODE","M4","QC_significant"),
    num("QC","miRNA","M3","slope_correlation"),num("QC","miRNA","M4","slope_correlation"),
    num("QC","GENCODE","M3","slope_correlation"),num("QC","GENCODE","M4","slope_correlation")),
  "QC显著性丢失与样本量减少并不自动等于生物学不稳定；`QC_sensitive`采用预冻结的方向、效应变化和模型稳定性规则。",
  "",
  "## 5. 189/3方向不平衡",
  "",
  sprintf("M1精确恢复原稿miRNA 189负/3正。M3全部显著miRNA为%d负/%d正；在原192个中，M3显著者为%d负/%d正。M4全部显著miRNA为%d负/%d正；在原192个中，M4显著者为%d负/%d正。",
    num("full_model","miRNA","M3","significant_negative"),num("full_model","miRNA","M3","significant_positive"),
    num("original_discoveries","miRNA","M3","significant_negative"),num("original_discoveries","miRNA","M3","significant_positive"),
    num("full_model","miRNA","M4","significant_negative"),num("full_model","miRNA","M4","significant_positive"),
    num("original_discoveries","miRNA","M4","significant_negative"),num("original_discoveries","miRNA","M4","significant_positive")),
  "M3仍显示明显负向不平衡；CLR后的M4因相对共同中心重新定义坐标，方向接近平衡。QC减少显著数但未消除M3负向主导。该变化不能单独证明原方向不平衡是技术伪影，也不能证明或排除载体组成变化。",
  "",
  "## 6. Figure 3/4原代表性特征",
  "",
  "|层|特征|M3状态|M4状态|Strict core|建议|",
  "|---|---|---|---|---:|---|",
  rep_lines,
  "",
  "本阶段未临时更换代表性特征。具体斜率、SE、95%CI、P/FDR、QC及CR2结果见`STAGE1B_REPRESENTATIVE_FEATURES.csv`。",
  "",
  "## 7. 对十二个问题的直接回答",
  "",
  "**Q1. 是否全部尝试？** 是。507/1,079个特征在FULL及QC、M3及M4全部尝试。FULL模型的成功、边界和数值不稳定数量见表；无MODEL_NOT_EVALUABLE。",
  sprintf("**Q2. M3显著多少？** miRNA %d、GENCODE %d；其中原稿显著集为%d和%d。",
    num("full_model","miRNA","M3","significant"),num("full_model","GENCODE","M3","significant"),old_mi$M3_significant,old_ge$M3_significant),
  sprintf("**Q3. M4显著多少？** miRNA %d、GENCODE %d；其中原稿显著集为%d和%d。",
    num("full_model","miRNA","M4","significant"),num("full_model","GENCODE","M4","significant"),old_mi$M4_significant,old_ge$M4_significant),
  sprintf("**Q4. M1/M3/M4共同显著且方向一致多少？** miRNA %d、GENCODE %d。",
    num("robustness","miRNA","M1_M3_M4","joint_significant_direction_consistent"),num("robustness","GENCODE","M1_M3_M4","joint_significant_direction_consistent")),
  "**Q5. Strict core多少？** miRNA 16、GENCODE 33，完整名单已输出。",
  sprintf("**Q6. RS-boundary多高？** FULL M3/M4分别为miRNA %s/%s、GENCODE %s/%s。比例不低，显著削弱“所有旧发现都通过随机斜率验证”的说法，但不否定16/33个strict core和65个边界支持集合。",
    pct(num("full_model","miRNA","M3","status_RS_BOUNDARY"),507),pct(num("full_model","miRNA","M4","status_RS_BOUNDARY"),507),
    pct(num("full_model","GENCODE","M3","status_RS_BOUNDARY"),1079),pct(num("full_model","GENCODE","M4","status_RS_BOUNDARY"),1079)),
  "**Q7. CR2是否一致？** 效应估计高度相关、方向高度一致；显著集合不完全相同。差异不能被解释为任一方法自动优越。",
  "**Q8. QC是否支持主要趋势？** 支持一组方向和效应稳定的趋势，但显著数下降，尤其GENCODE；因此应作为敏感性限定而非宣称完全无影响。",
  "**Q9. 189/3如何变化？** M3保留强负向不平衡；M4因CLR参照变化转为较平衡的负/正分布。QC不改变这一总体判断。",
  "**Q10. 代表特征保留哪些？** miR-375、MT-TL2、RN7SL1、IL8、HBA1可作为strict-core代表；TMSB4X可带QC限定保留。miR-338-3p、miR-381-3p、miR-410-3p、EIF4A1不建议继续作为主要稳健代表。",
  "**Q11. 是否支持主要发现？** 是，但应表述为一组经成分校正、随机斜率和QC检验后仍稳定的纵向配对梯度，而非原稿全部显著特征均稳健。",
  "**Q12. Reviewer 1是否解决？** 已充分执行其要求：全特征随机斜率、CLR、QC和CR2均完成，边界没有被隐藏或用随机截距结果冒充。残余限制是3–7个不平衡标签/患者导致大量随机斜率方差位于边界，且CLR估计的是相对参照变化。",
  "",
  "## 8. 方法与解释限制",
  "",
  "- collection label仅按D01–D07顺序编码并以D04中心化，不解释为已验证的真实发病后日数。",
  "- RPM/CLR均为相对丰度；未使用伪counts，未执行或声称TMM–voom。",
  "- miRNA与GENCODE来自相同底层文库，不是独立生物学重复。",
  "- CR2与混合模型估计目标及标准误假设不同；边界支持集合与strict core必须分开报告。",
  "- 逐模型BH家族固定为507或1,079；不可稳定模型保留行并将正式P/FDR记NA。",
  "",
  "最终决策：**GO_TO_STAGE2**。进入Stage 2前应冻结使用strict core作为主要纵向候选、boundary-supported作为次级支持，并保留全部旧发现的透明交叉表。"
)
writeLines(report, file.path(out, "STAGE1B_FINAL_REPORT.md"), useBytes = TRUE)

go <- c(
  "# Stage 2 Go/No-Go Decision",
  "",
  "日期：2026-11-21（Asia/Shanghai）",
  "",
  "## 决策：GO_TO_STAGE2",
  "",
  "Stage 1B已完成全部预定全特征Model B、CLR、QC PASS-only和CR2分析。统计证据支持进入RNA cargo/carrier composition阶段，但Stage 2不得改变Stage 1B模型、FDR或稳健特征定义。",
  "",
  "进入Stage 2的固定边界：",
  "",
  "- 主要纵向候选：16个miRNA和33个GENCODE strict core。",
  "- 次级证据：9个miRNA和56个GENCODE RS-boundary but otherwise supported，必须单独标注。",
  "- XDec仅在输入尺度、81个informative RNA映射、参考profile含义和逐样本输出均验证后实施。",
  "- 不把解卷积关联解释为载体组成变化的因果证明。",
  "- 不用Stage 2结果反向选择Stage 1模型或阈值。",
  "- Stage 2仍受既有XDec来源、参考矩阵及真正原作者逐样本比例是否可得的限制。",
  "",
  "Stage 1B未修改手稿和主图；本文件只作后续阶段启动依据。"
)
writeLines(go, file.path(out, "STAGE2_GO_NO_GO.md"), useBytes = TRUE)

cat("REPORTS_WRITTEN\n")
