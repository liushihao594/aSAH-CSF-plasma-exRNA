options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(clubSandwich)
  library(ggplot2)
})

set.seed(20261121L)
out <- "revision_stage2"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "intermediate"), recursive = TRUE, showWarnings = FALSE)

tol <- 1e-4
optimizers <- c("nloptwrap", "bobyqa", "nlminbwrap")
obj <- readRDS("revision_stage1/intermediate/prepared_stage1_inputs.rds")

dir_sign <- function(x) ifelse(is.na(x), NA_character_, ifelse(x > sqrt(.Machine$double.eps), "positive",
  ifelse(x < -sqrt(.Machine$double.eps), "negative", "zero")))

fit_once <- function(dat, structure = "B", optimizer = "nloptwrap") {
  f <- switch(structure,
    B = y ~ centered_label + (1 + centered_label || subject),
    C = y ~ centered_label + (0 + centered_label | subject))
  warnings <- character()
  ctrl <- if (optimizer == "nloptwrap") list(maxeval = 200000) else list(maxfun = 200000)
  fit <- tryCatch(withCallingHandlers(
    lmerTest::lmer(f, data = dat, REML = TRUE,
      control = lmerControl(optimizer = optimizer, optCtrl = ctrl,
        check.conv.singular = lme4:::.makeCC(action = "ignore", tol = tol))),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(fit, "error")) return(list(fit = NULL, ok = FALSE, error = conditionMessage(fit)))
  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(co) || !"centered_label" %in% rownames(co)) return(list(fit = fit, ok = FALSE, error = "fixed effect unavailable"))
  msgs <- fit@optinfo$conv$lme4$messages
  code <- fit@optinfo$conv$opt
  code_bad <- length(code) && any(suppressWarnings(as.numeric(code)) != 0, na.rm = TRUE)
  est <- unname(co["centered_label", "Estimate"])
  se <- unname(co["centered_label", "Std. Error"])
  df <- if ("df" %in% colnames(co)) unname(co["centered_label", "df"]) else NA_real_
  tc <- grep("t value", colnames(co), value = TRUE)[1]
  pc <- grep("Pr\\(", colnames(co), value = TRUE)[1]
  tv <- if (!is.na(tc)) unname(co["centered_label", tc]) else NA_real_
  pv <- if (!is.na(pc)) unname(co["centered_label", pc]) else NA_real_
  vc <- as.data.frame(VarCorr(fit))
  getv <- function(v) { q <- vc$vcov[vc$var1 == v & is.na(vc$var2)]; if (length(q)) q[1] else NA_real_ }
  ri <- getv("(Intercept)"); rs <- getv("centered_label"); rv <- vc$vcov[vc$grp == "Residual"][1]
  ri_ratio <- if (is.finite(ri) && is.finite(rv) && rv > 0) sqrt(max(ri,0)/rv) else NA_real_
  rs_ratio <- if (is.finite(rs) && is.finite(rv) && rv > 0) sqrt(max(rs,0)/rv) else NA_real_
  ok <- is.finite(est) && is.finite(se) && se > 0 && !code_bad && !length(msgs)
  list(fit = fit, ok = ok, error = "", warnings = paste(unique(warnings), collapse = " | "),
    estimate = est, SE = se, df = df, t = tv, P = pv, CI_lower = est - 1.96*se, CI_upper = est + 1.96*se,
    RI_variance = ri, RS_variance = rs, residual_variance = rv, RI_SD_ratio = ri_ratio,
    RS_SD_ratio = rs_ratio, singular = isSingular(fit, tol = tol))
}

fit_modelB <- function(y, meta) {
  dat <- data.frame(y = as.numeric(y), centered_label = as.numeric(meta$centered_label),
    subject = factor(meta$subject))
  good <- is.finite(dat$y) & is.finite(dat$centered_label) & !is.na(dat$subject)
  dat <- dat[good,,drop=FALSE]
  fits <- lapply(optimizers, function(o) fit_once(dat, "B", o))
  valid <- vapply(fits, function(q) isTRUE(q$ok), logical(1))
  first <- if (any(valid)) which(valid)[1] else which(vapply(fits, function(q) !is.null(q$fit), logical(1)))[1]
  if (!length(first) || is.na(first)) return(data.frame(n = nrow(dat), n_subjects = nlevels(dat$subject),
    slope = NA, SE = NA, df = NA, P = NA, CI_lower = NA, CI_upper = NA,
    model_status = "MODEL_NOT_EVALUABLE", selected_optimizer = NA, optimizer_valid_n = 0,
    optimizer_stable = FALSE, RI_SD_ratio = NA, RS_SD_ratio = NA, CR2_slope = NA,
    CR2_SE = NA, CR2_df = NA, CR2_P = NA))
  q <- fits[[first]]
  vals <- fits[valid]
  stable <- FALSE
  if (length(vals) >= 2) {
    es <- vapply(vals, `[[`, numeric(1), "estimate")
    ss <- vapply(vals, `[[`, numeric(1), "SE")
    stable <- length(unique(dir_sign(es))) == 1 && diff(range(es))/median(ss) <= .10 && max(ss)/min(ss) <= 1.10
  }
  ri0 <- is.finite(q$RI_SD_ratio) && q$RI_SD_ratio <= tol
  rs0 <- is.finite(q$RS_SD_ratio) && q$RS_SD_ratio <= tol
  status <- if (length(vals) < 2 || !stable) "NUMERICAL_UNSTABLE" else if (ri0 && rs0) {
    "BOTH_VARIANCE_BOUNDARY_UNRESOLVED"
  } else if (rs0) {
    "RS_BOUNDARY"
  } else if (ri0) {
    cfits <- lapply(optimizers, function(o) fit_once(dat, "C", o))
    cv <- cfits[vapply(cfits, function(v) isTRUE(v$ok), logical(1))]
    cstable <- length(cv) >= 2
    if (cstable) {
      ce <- vapply(cv, `[[`, numeric(1), "estimate"); cs <- vapply(cv, `[[`, numeric(1), "SE")
      cstable <- length(unique(dir_sign(ce))) == 1 && diff(range(ce))/median(cs) <= .10 && max(cs)/min(cs) <= 1.10 &&
        abs(q$estimate - ce[1])/median(c(q$SE,cs[1])) <= .10
    }
    if (cstable) "RI_BOUNDARY_C_CONFIRMED" else "RI_BOUNDARY_C_UNRESOLVED"
  } else "REGULAR_RANDOM_SLOPE"
  if (status == "NUMERICAL_UNSTABLE") q$P <- NA_real_

  cr <- tryCatch({
    lmfit <- lm(y ~ centered_label + subject, data = dat)
    ct <- clubSandwich::coef_test(lmfit, vcov = "CR2", cluster = dat$subject, test = "Satterthwaite")
    zz <- ct[rownames(ct) == "centered_label",,drop=FALSE]
    c(slope=zz$beta, SE=zz$SE, df=zz$df_Satt, P=zz$p_Satt)
  }, error = function(e) c(slope=NA,SE=NA,df=NA,P=NA))
  data.frame(n = nrow(dat), n_subjects = nlevels(dat$subject), slope = q$estimate, SE = q$SE,
    df = q$df, P = q$P, CI_lower = q$CI_lower, CI_upper = q$CI_upper,
    model_status = status, selected_optimizer = optimizers[first], optimizer_valid_n = sum(valid),
    optimizer_stable = stable, RI_SD_ratio = q$RI_SD_ratio, RS_SD_ratio = q$RS_SD_ratio,
    CR2_slope = cr["slope"], CR2_SE = cr["SE"], CR2_df = cr["df"], CR2_P = cr["P"])
}

bind_fill <- function(xs) {
  cols <- unique(unlist(lapply(xs, names)))
  do.call(rbind, lapply(xs, function(x) { miss <- setdiff(cols,names(x)); for (m in miss) x[[m]] <- NA; x[,cols,drop=FALSE] }))
}

## Sample metadata and fixed high-abundance set
pairs <- obj$inputs$miRNA_FULL$pairs
sample_meta <- rbind(
  data.frame(pair_id=pairs$pair_id, subject=pairs$subject, collection_label=pairs$collection_label,
    centered_label=pairs$centered_label, expression_id=pairs$CSF_expression_id, GSM=pairs$CSF_GSM,
    biofluid="CSF", QC=pairs$CSF_QC, both_QC_PASS=as.logical(pairs$both_QC_PASS),
    transcriptome_reads=pairs$CSF_transcriptome_reads, reference_genome_reads=pairs$CSF_reference_genome_reads),
  data.frame(pair_id=pairs$pair_id, subject=pairs$subject, collection_label=pairs$collection_label,
    centered_label=pairs$centered_label, expression_id=pairs$plasma_expression_id, GSM=pairs$plasma_GSM,
    biofluid="plasma", QC=pairs$plasma_QC, both_QC_PASS=as.logical(pairs$both_QC_PASS),
    transcriptome_reads=pairs$plasma_transcriptome_reads, reference_genome_reads=pairs$plasma_reference_genome_reads))
sample_meta$label_numeric <- as.integer(sub("D", "", sample_meta$collection_label))
stopifnot(nrow(sample_meta)==306, length(unique(sample_meta$expression_id))==306)

mi_all <- cbind(obj$inputs$miRNA_FULL$csf, obj$inputs$miRNA_FULL$plasma)
mi_all <- mi_all[,sample_meta$expression_id,drop=FALSE]
top10 <- names(sort(apply(mi_all,1,median,na.rm=TRUE), decreasing=TRUE))[1:10]
high_features <- unique(c("hsa-miR-486-5p", "hsa-miR-451a", top10))
high_map <- data.frame(feature=high_features, layer="miRNA", rank_median_RPM=match(high_features,top10))
high_map <- rbind(high_map, data.frame(feature="RNY4:misc_RNA",layer="GENCODE",rank_median_RPM=NA_integer_))
write.csv(high_map,file.path(out,"intermediate","STAGE2_FIXED_HIGH_ABUNDANCE_FEATURES.csv"),row.names=FALSE,na="")

## High-abundance sample rows and longitudinal models
high_rows <- list(); high_models <- list(); hi <- 0; hm <- 0
for (i in seq_len(nrow(high_map))) {
  layer <- high_map$layer[i]; feature <- high_map$feature[i]
  inp <- obj$inputs[[paste0(layer,"_FULL")]]
  if (!feature %in% rownames(inp$csf)) next
  vals <- c(as.numeric(inp$csf[feature,]), as.numeric(inp$plasma[feature,]))
  d0 <- cbind(sample_meta, feature=feature, layer=layer, RPM=vals, log2_RPM_plus1=log2(vals+1))
  for (dataset in c("FULL","QC_PASS_ONLY")) {
    dd <- if(dataset=="FULL")d0 else d0[d0$both_QC_PASS,,drop=FALSE]
    for (fluid in c("CSF","plasma")) {
      d <- dd[dd$biofluid==fluid,,drop=FALSE]
      fm <- fit_modelB(d$log2_RPM_plus1,d)
      hm <- hm+1; high_models[[hm]] <- cbind(data.frame(feature=feature,layer=layer,dataset=dataset,biofluid=fluid),fm)
      for (lab in sort(unique(d$collection_label))) {
        q <- d[d$collection_label==lab,,drop=FALSE]
        hi <- hi+1; high_rows[[hi]] <- data.frame(feature=feature,layer=layer,dataset=dataset,biofluid=fluid,
          collection_label=lab,label_numeric=unique(q$label_numeric),n=nrow(q),n_subjects=length(unique(q$subject)),
          RPM_median=median(q$RPM),RPM_Q1=unname(quantile(q$RPM,.25)),RPM_Q3=unname(quantile(q$RPM,.75)),
          RPM_mean=mean(q$RPM),RPM_SE=sd(q$RPM)/sqrt(nrow(q)),log2_mean=mean(q$log2_RPM_plus1),
          log2_SE=sd(q$log2_RPM_plus1)/sqrt(nrow(q)))
      }
    }
  }
}
high_models <- do.call(rbind,high_models)
high_models$FDR <- ave(high_models$P,high_models$dataset,high_models$biofluid,
  FUN=function(p)p.adjust(p,"BH",n=length(p)))
high_models$CR2_FDR <- ave(high_models$CR2_P,high_models$dataset,high_models$biofluid,
  FUN=function(p)p.adjust(p,"BH",n=length(p)))
high_results <- merge(do.call(rbind,high_rows),high_models,by=c("feature","layer","dataset","biofluid"),all.x=TRUE)
write.csv(high_results,file.path(out,"STAGE2_HIGH_ABUNDANCE_RNA_RESULTS.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(high_models,file.path(out,"intermediate","STAGE2_HIGH_ABUNDANCE_MODELS.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## Biotype matrices, using the same exceRpt RPM products and exact main sample IDs
paths <- c(
  miRNA="revision_stage0/source_extract/TRACKB_SOURCE_DATA_PACKAGE/01_GEO_processed_expression/GSE121868_KJENS1-aSAH_Project-2016-10-13_exceRpt_miRNA_ReadsPerMillion.txt.gz",
  GENCODE="revision_stage0/source_extract/TRACKB_SOURCE_DATA_PACKAGE/01_GEO_processed_expression/GSE121868_KJENS1-aSAH_Project-2016-10-13_exceRpt_gencode_ReadsPerMillion.txt.gz",
  piRNA="revision_stage0b/02_xdec_cargo/GSE121868_piRNA_http_retry.gz",
  tRNA="revision_stage0b/02_xdec_cargo/GSE121868_tRNA_http_retry.gz")
validation <- list(); totals <- matrix(0,nrow=7,ncol=306,dimnames=list(c("miRNA","piRNA","tRNA","Y_RNA","lincRNA","snoRNA","snRNA"),sample_meta$expression_id))
for (nm in names(paths)) {
  d <- read.delim(gzfile(paths[nm]),row.names=1,check.names=FALSE,comment.char="",quote="")
  validation[[nm]] <- data.frame(source_layer=nm,path=paths[nm],features=nrow(d),samples=ncol(d),
    primary_samples_matched=sum(sample_meta$expression_id %in% colnames(d)),NA_count=sum(!is.finite(as.matrix(d))),
    negative_count=sum(as.matrix(d)<0,na.rm=TRUE),min_value=min(as.matrix(d),na.rm=TRUE),max_value=max(as.matrix(d),na.rm=TRUE))
  stopifnot(all(sample_meta$expression_id %in% colnames(d)))
  d <- d[,sample_meta$expression_id,drop=FALSE]
  if(nm=="miRNA") totals["miRNA",] <- colSums(d)
  if(nm=="piRNA") totals["piRNA",] <- colSums(d)
  if(nm=="tRNA") totals["tRNA",] <- colSums(d)
  if(nm=="GENCODE") {
    ids <- rownames(d)
    idx <- list(
      Y_RNA=grepl("^(RNY|Y_RNA\\.)",ids),
      lincRNA=grepl(":lincRNA$",ids),
      snoRNA=grepl(":snoRNA$",ids),
      snRNA=grepl(":snRNA$",ids))
    for(k in names(idx)) totals[k,] <- colSums(d[idx[[k]],,drop=FALSE])
    validation[["GENCODE_subtypes"]] <- data.frame(source_layer=names(idx),path=paths[nm],
      features=vapply(idx,sum,integer(1)),samples=ncol(d),primary_samples_matched=ncol(d),NA_count=0,
      negative_count=0,min_value=vapply(idx,function(ii)min(as.matrix(d[ii,,drop=FALSE])),numeric(1)),
      max_value=vapply(idx,function(ii)max(as.matrix(d[ii,,drop=FALSE])),numeric(1)))
  }
  rm(d); gc(FALSE)
}
write.csv(do.call(rbind,validation),file.path(out,"intermediate","STAGE2_BIOTYPE_INPUT_VALIDATION.csv"),row.names=FALSE,na="")

denom <- colSums(totals)
prop <- sweep(totals,2,denom,"/")
logtot <- log(totals + .5)
clr <- sweep(logtot,2,colMeans(logtot),"-")
bio_base <- do.call(rbind,lapply(rownames(totals),function(bt)data.frame(sample_meta,biotype=bt,
  biotype_RPM_total=as.numeric(totals[bt,]),selected_7type_RPM_total=denom,
  selected_7type_proportion=as.numeric(prop[bt,]),biotype_CLR=as.numeric(clr[bt,]))))
bio_rows <- list(); bio_models <- list(); bi<-0; bm<-0
for(dataset in c("FULL","QC_PASS_ONLY")) for(fluid in c("CSF","plasma")) for(bt in rownames(totals)) {
  d <- bio_base[bio_base$biofluid==fluid & bio_base$biotype==bt & (dataset=="FULL" | bio_base$both_QC_PASS),,drop=FALSE]
  mt <- fit_modelB(log2(d$biotype_RPM_total+1),d); mc <- fit_modelB(d$biotype_CLR,d)
  bm<-bm+1; bio_models[[bm]]<-cbind(data.frame(dataset=dataset,biofluid=fluid,biotype=bt),
    setNames(mt,paste0("total_",names(mt))),setNames(mc,paste0("CLR_",names(mc))))
  d$dataset<-dataset; bi<-bi+1; bio_rows[[bi]]<-d
}
bio_models<-do.call(rbind,bio_models)
bio_models$total_FDR<-ave(bio_models$total_P,bio_models$dataset,bio_models$biofluid,FUN=function(p)p.adjust(p,"BH",n=7))
bio_models$CLR_FDR<-ave(bio_models$CLR_P,bio_models$dataset,bio_models$biofluid,FUN=function(p)p.adjust(p,"BH",n=7))
bio_out<-merge(do.call(rbind,bio_rows),bio_models,by=c("dataset","biofluid","biotype"),all.x=TRUE)
write.csv(bio_out,file.path(out,"STAGE2_BIOTYPE_COMPOSITION.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(bio_models,file.path(out,"intermediate","STAGE2_BIOTYPE_TEMPORAL_MODELS.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## Global paired center shifts and pre-specified associations
bio_wide <- reshape(bio_base[,c("expression_id","biotype","biotype_CLR")],idvar="expression_id",timevar="biotype",direction="wide")
global_pair_rows<-list(); global_model_rows<-list(); corr_rows<-list(); loo_rows<-list(); gi<-gm<-ci<-li<-0
for(layer in c("miRNA","GENCODE")) for(dataset in c("FULL","QC_PASS_ONLY")) {
  inp<-obj$inputs[[paste0(layer,"_",dataset)]]; pp<-inp$pairs
  csf_center<-colMeans(log2(inp$csf+1)); plasma_center<-colMeans(log2(inp$plasma+1)); shift<-colMeans(inp$raw)
  stopifnot(max(abs(shift-(csf_center-plasma_center)))<1e-10)
  d<-data.frame(pp,layer=layer,dataset=dataset,center_shift=shift,CSF_log_center=csf_center,
    plasma_log_center=plasma_center,identity_error=shift-(csf_center-plasma_center))
  ## observed high-RNA gradients
  mi<-obj$inputs[[paste0("miRNA_",dataset)]]; ge<-obj$inputs[[paste0("GENCODE_",dataset)]]
  d$gradient_miR486<-as.numeric(mi$raw["hsa-miR-486-5p",])
  d$gradient_miR451a<-as.numeric(mi$raw["hsa-miR-451a",])
  d$gradient_RNY4<-as.numeric(ge$raw["RNY4:misc_RNA",])
  d$mapped_read_logratio<-log2(d$CSF_transcriptome_reads+1)-log2(d$plasma_transcriptome_reads+1)
  bc<-bio_wide[match(d$CSF_expression_id,bio_wide$expression_id),]
  bp<-bio_wide[match(d$plasma_expression_id,bio_wide$expression_id),]
  for(bt in rownames(totals)) d[[paste0("biotype_CLR_gradient_",bt)]]<-bc[[paste0("biotype_CLR.",bt)]]-bp[[paste0("biotype_CLR.",bt)]]
  gi<-gi+1; global_pair_rows[[gi]]<-d
  for(outcome in c("center_shift","CSF_log_center","plasma_log_center")) {
    fm<-fit_modelB(d[[outcome]],d); gm<-gm+1; global_model_rows[[gm]]<-cbind(data.frame(layer=layer,dataset=dataset,outcome=outcome),fm)
  }
  assoc<-c("gradient_miR486","gradient_miR451a","gradient_RNY4","mapped_read_logratio",paste0("biotype_CLR_gradient_",rownames(totals)))
  for(v in assoc){
    ct<-suppressWarnings(cor.test(d$center_shift,d[[v]],method="spearman",exact=FALSE))
    ci<-ci+1; corr_rows[[ci]]<-data.frame(layer=layer,dataset=dataset,variable=v,n=sum(complete.cases(d$center_shift,d[[v]])),
      spearman_rho=unname(ct$estimate),P=ct$p.value)
  }
  for(s in unique(d$subject)) {
    q<-d[d$subject!=s,,drop=FALSE]; fm<-fit_modelB(q$center_shift,q); li<-li+1
    loo_rows[[li]]<-cbind(data.frame(layer=layer,dataset=dataset,excluded_subject=s),fm)
  }
}
pair_table<-do.call(rbind,global_pair_rows); model_table<-do.call(rbind,global_model_rows); corr_table<-do.call(rbind,corr_rows); loo_table<-do.call(rbind,loo_rows)
model_table$FDR<-ave(model_table$P,model_table$layer,model_table$dataset,FUN=function(p)p.adjust(p,"BH",n=3))
model_table$CR2_FDR<-ave(model_table$CR2_P,model_table$layer,model_table$dataset,FUN=function(p)p.adjust(p,"BH",n=3))
corr_table$FDR<-ave(corr_table$P,corr_table$layer,corr_table$dataset,FUN=function(p)p.adjust(p,"BH",n=length(p)))
global_out<-bind_fill(list(
  cbind(record_type="PAIR",pair_table),
  cbind(record_type="MODEL",model_table),
  cbind(record_type="ASSOCIATION",corr_table),
  cbind(record_type="LOO",loo_table)))
write.csv(global_out,file.path(out,"STAGE2_GLOBAL_SHIFT_ANALYSIS.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
saveRDS(list(pair=pair_table,models=model_table,associations=corr_table,LOO=loo_table,
  high_models=high_models,biotype_models=bio_models),file.path(out,"intermediate","STAGE2_MODEL_OBJECTS_AND_TABLES.rds"))

## QC sensitivity comparison
compare_qc <- function(tab, keys, slope="slope", fdr="FDR", status="model_status", family) {
  f<-tab[tab$dataset=="FULL",c(keys,slope,fdr,status),drop=FALSE]
  q<-tab[tab$dataset=="QC_PASS_ONLY",c(keys,slope,fdr,status),drop=FALSE]
  names(f)[-(seq_along(keys))]<-paste0("FULL_",c("slope","FDR","status"))
  names(q)[-(seq_along(keys))]<-paste0("QC_",c("slope","FDR","status"))
  x<-merge(f,q,by=keys,all=TRUE); x$family<-family
  x$direction_consistent<-dir_sign(x$FULL_slope)==dir_sign(x$QC_slope)
  x$absolute_slope_change<-abs(x$QC_slope-x$FULL_slope)
  x
}
qc_high<-compare_qc(high_models,c("feature","layer","biofluid"),family="HIGH_ABUNDANCE")
qc_bio_total<-compare_qc(transform(bio_models,slope=total_slope,FDR=total_FDR,model_status=total_model_status),c("biofluid","biotype"),family="BIOTYPE_TOTAL")
qc_bio_clr<-compare_qc(transform(bio_models,slope=CLR_slope,FDR=CLR_FDR,model_status=CLR_model_status),c("biofluid","biotype"),family="BIOTYPE_CLR")
qc_global<-compare_qc(model_table,c("layer","outcome"),family="GLOBAL_SHIFT")
write.csv(bind_fill(list(qc_high,qc_bio_total,qc_bio_clr,qc_global)),file.path(out,"STAGE2_QC_SENSITIVITY.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## Stage 1B strict-core assessment; XDec/cargo adjustment deliberately absent after NO-GO
core<-read.csv("revision_stage1/stage1b_full_refit/STAGE1B_STRICT_CORE_ROBUST.csv",check.names=FALSE)
full_mi<-read.csv("revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_miRNA.csv",check.names=FALSE)
full_ge<-read.csv("revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv",check.names=FALSE)
full<-rbind(full_mi,full_ge)
keep<-c("feature","layer","M1_slope","M1_SE","M1_FDR","M2_slope","M2_SE","M2_FDR","M3_slope","M3_SE","M3_FDR",
  "M3_model_status","M4_slope","M4_SE","M4_FDR","M4_model_status","QC_M3_slope","QC_M3_FDR","QC_M4_slope","QC_M4_FDR",
  "M2_reference_sensitivity_excluded_reference_feature","M2_reference_sensitivity_excluded_ref_slope",
  "M2_reference_sensitivity_excluded_ref_FDR","M2_reference_sensitivity_pc0_5_slope","M2_reference_sensitivity_pc0_5_FDR",
  "M2_reference_sensitivity_top5_excluded_slope","M2_reference_sensitivity_top5_excluded_FDR")
core_out<-merge(core[,c("feature","layer","strict_core_robust","primary_robustness_category")],full[,keep],by=c("feature","layer"),all.x=TRUE)
core_out$M1_direction<-dir_sign(core_out$M1_slope); core_out$M3_direction<-dir_sign(core_out$M3_slope); core_out$M4_direction<-dir_sign(core_out$M4_slope)
core_out$M3_to_M4_slope_change<-core_out$M4_slope-core_out$M3_slope
core_out$M4_to_M3_absolute_ratio<-ifelse(abs(core_out$M3_slope)>sqrt(.Machine$double.eps),abs(core_out$M4_slope/core_out$M3_slope),NA)
core_out$fixed_high_abundance_target<-core_out$feature %in% high_map$feature
core_out$cargo_adjusted_status<-"NOT_EVALUATED_XDEC_NO_GO"
core_out$stage2_support_statement<-"SUPPORTED_BY_STAGE1B_RAW_AND_CLR_RANDOM_SLOPE_AND_QC;_NOT_CARRIER_ADJUSTED"
core_out$recommended_as_existing_primary_representative<-core_out$feature %in% c("hsa-miR-375","MT-TL2:Mt_tRNA","RN7SL1:antisense","IL8:protein_coding","HBA1:protein_coding")
write.csv(core_out,file.path(out,"STAGE2_CORE_FEATURE_ASSESSMENT.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## Exploratory diagnostic figures
save_plot <- function(filename, plot, width, height) {
  target <- file.path(out,"plots",filename)
  Cairo::CairoPNG(filename=target,width=round(width*180),height=round(height*180),res=180)
  print(plot)
  grDevices::dev.off()
}
stage1_mi<-full_mi
direction_plot<-rbind(
  data.frame(model="M1",direction=dir_sign(stage1_mi$M1_slope),significant=stage1_mi$M1_FDR<.05),
  data.frame(model="M3",direction=dir_sign(stage1_mi$M3_slope),significant=stage1_mi$M3_FDR<.05),
  data.frame(model="M4",direction=dir_sign(stage1_mi$M4_slope),significant=stage1_mi$M4_FDR<.05))
dp<-aggregate(significant~model+direction,direction_plot,sum,na.rm=TRUE)
p1<-ggplot(dp,aes(model,significant,fill=direction))+geom_col(position="dodge")+
  labs(x=NULL,y="FDR-significant miRNAs",title="Slope direction among significant miRNAs",
    subtitle="M1 limma raw gradient; M3 random-slope raw gradient; M4 random-slope CLR gradient",
    caption="FULL: 30 subjects, 153 paired observations; fixed 507-miRNA family")+theme_bw()
save_plot("STAGE2_MIRNA_SLOPE_DIRECTION.png",p1,width=8,height=5)

hp<-high_results[high_results$dataset=="FULL" & high_results$feature %in% c("hsa-miR-486-5p","hsa-miR-451a","RNY4:misc_RNA"),]
p2<-ggplot(hp,aes(label_numeric,RPM_median,color=biofluid,group=biofluid))+geom_line()+geom_point()+
  facet_wrap(~feature,scales="free_y")+scale_x_continuous(breaks=1:7,labels=paste0("D0",1:7))+
  labs(x="Collection label",y="Median RPM",title="High-abundance RNA trajectories by biofluid",
    caption="FULL; medians across available subjects at each label; RPM is relative abundance")+theme_bw()
save_plot("STAGE2_HIGH_ABUNDANCE_TRAJECTORIES.png",p2,width=10,height=5)

bs<-bio_out[bio_out$dataset=="FULL",]
bs<-aggregate(selected_7type_proportion~biofluid+biotype+label_numeric,bs,median)
p3<-ggplot(bs,aes(label_numeric,selected_7type_proportion,color=biotype,group=biotype))+geom_line()+geom_point()+
  facet_wrap(~biofluid)+scale_x_continuous(breaks=1:7,labels=paste0("D0",1:7))+
  labs(x="Collection label",y="Median proportion within seven selected RPM layers",title="Operational RNA-biotype composition",
    caption="FULL; selected exceRpt RPM layers only; not physical cargo proportions")+theme_bw()+theme(legend.position="bottom")
save_plot("STAGE2_BIOTYPE_COMPOSITION.png",p3,width=11,height=6)

gp<-pair_table[pair_table$dataset=="FULL",]
p4<-ggplot(gp,aes(as.factor(collection_label),center_shift,fill=layer))+geom_boxplot(outlier.size=.8)+
  facet_wrap(~layer,scales="free_y")+labs(x="Collection label",y="Mean raw paired gradient across fixed reference",
    title="Global composition-reference shift",caption="FULL: 153 paired observations; g_CLR = g_raw - displayed center shift")+theme_bw()+theme(legend.position="none")
save_plot("STAGE2_GLOBAL_SHIFT.png",p4,width=10,height=5)

qcp<-merge(model_table[model_table$dataset=="FULL",c("layer","outcome","slope")],
  model_table[model_table$dataset=="QC_PASS_ONLY",c("layer","outcome","slope")],by=c("layer","outcome"),suffixes=c("_FULL","_QC"))
p5<-ggplot(qcp,aes(slope_FULL,slope_QC,color=layer,label=outcome))+geom_abline(slope=1,intercept=0,linetype=2)+geom_point(size=3)+geom_text(vjust=-.7)+
  labs(x="FULL slope",y="QC PASS-only slope",title="QC sensitivity of global-shift outcomes",caption="Model B slopes; QC branch retains 141 paired observations")+theme_bw()
save_plot("STAGE2_GLOBAL_SHIFT_QC.png",p5,width=7,height=6)

writeLines(capture.output(sessionInfo()),file.path(out,"logs","STAGE2_sessionInfo.txt"),useBytes=TRUE)
cat(sprintf("STAGE2_ANALYSIS_COMPLETE high=%d biotype_rows=%d global_records=%d core=%d\n",
  nrow(high_results),nrow(bio_out),nrow(global_out),nrow(core_out)))
