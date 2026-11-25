options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(clubSandwich)
})

args <- commandArgs(trailingOnly = TRUE)
worker_id <- if (length(args) >= 1L) as.integer(args[1]) else 1L
n_workers <- if (length(args) >= 2L) as.integer(args[2]) else 1L
if (!is.finite(worker_id) || !is.finite(n_workers) || worker_id < 1L || worker_id > n_workers) stop("Invalid worker assignment")

stage1 <- "revision_stage1"
out <- file.path(stage1, "stage1b_full_refit")
checkpoint_dir <- file.path(out, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "intermediate"), recursive = TRUE, showWarnings = FALSE)

tol_singular <- 1e-4
optimizers <- c("nloptwrap", "bobyqa", "nlminbwrap")
chunk_size <- 100L
set.seed(20261121L)

obj <- readRDS(file.path(stage1, "intermediate", "prepared_stage1_inputs.rds"))

representatives <- list(
  miRNA = c("hsa-miR-338-3p", "hsa-miR-381-3p", "hsa-miR-410-3p", "hsa-miR-375"),
  GENCODE = c("MT-TL2:Mt_tRNA", "RN7SL1:antisense", "IL8:protein_coding", "HBA1:protein_coding",
    "TMSB4X:protein_coding", "EIF4A1:retained_intron")
)

tasks <- list()
task_id <- 0L
for (layer in c("miRNA", "GENCODE")) for (dataset in c("FULL", "QC_PASS_ONLY")) for (gradient_model in c("M3", "M4")) {
  inp <- obj$inputs[[paste(layer, dataset, sep = "_")]]
  for (feature in rownames(inp$raw)) {
    task_id <- task_id + 1L
    tasks[[task_id]] <- data.frame(task_id = task_id, layer = layer, dataset = dataset,
      gradient_model = gradient_model, feature = feature, stringsAsFactors = FALSE)
  }
}
tasks <- do.call(rbind, tasks)
tasks$chunk_id <- ceiling(tasks$task_id / chunk_size)
write.csv(tasks, file.path(out, "intermediate", "STAGE1B_TASK_TABLE.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

safe_ratio <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b != 0, a/b, NA_real_)
dir_sign <- function(x) ifelse(is.na(x), NA_character_, ifelse(x > sqrt(.Machine$double.eps), "positive",
  ifelse(x < -sqrt(.Machine$double.eps), "negative", "zero")))

formulas <- list(
  A = y ~ centered_label + (1 + centered_label | subject),
  B = y ~ centered_label + (1 + centered_label || subject),
  C = y ~ centered_label + (0 + centered_label | subject),
  D = y ~ centered_label + (1 | subject)
)

fit_once <- function(dat, structure, optimizer) {
  warnings <- character()
  opt_ctrl <- if (optimizer == "nloptwrap") list(maxeval = 200000) else list(maxfun = 200000)
  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(formulas[[structure]], data = dat, REML = TRUE,
        control = lmerControl(optimizer = optimizer, optCtrl = opt_ctrl,
          check.conv.singular = lme4:::.makeCC(action = "ignore", tol = tol_singular))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(fit = NULL, stats = data.frame(
      structure = structure, optimizer = optimizer, fit_error = conditionMessage(fit),
      warnings = paste(unique(warnings), collapse = " | "), optimizer_code = NA_character_,
      convergence_messages = NA_character_, convergence_ok = FALSE, fixed_effect_finite = FALSE,
      slope = NA_real_, SE = NA_real_, df = NA_real_, t = NA_real_, P_raw = NA_real_,
      CI_lower = NA_real_, CI_upper = NA_real_, random_intercept_variance = NA_real_,
      random_slope_variance = NA_real_, intercept_slope_correlation = NA_real_, residual_variance = NA_real_,
      RI_SD_ratio = NA_real_, RS_SD_ratio_per_label = NA_real_, RS_SD_ratio_observed_span = NA_real_,
      covariance_eigen_min = NA_real_, covariance_eigen_max = NA_real_, covariance_eigen_ratio = NA_real_,
      rePCA_min_SD = NA_real_, rePCA_max_SD = NA_real_, rePCA_SD_ratio = NA_real_, rePCA_rank = NA_integer_,
      hessian_eigen_min = NA_real_, hessian_eigen_max = NA_real_, hessian_positive = FALSE,
      max_abs_scaled_gradient = NA_real_, isSingular_tol_1e_4 = NA, stringsAsFactors = FALSE)))
  }

  msgs <- fit@optinfo$conv$lme4$messages
  msgs <- if (length(msgs)) paste(msgs, collapse = " | ") else ""
  code <- fit@optinfo$conv$opt
  code <- if (is.null(code) || !length(code)) "0" else paste(code, collapse = ";")
  code_num <- suppressWarnings(as.numeric(unlist(strsplit(code, ";", fixed = TRUE))))
  opt_bad <- any(is.finite(code_num) & code_num != 0)
  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  fixed_ok <- !is.null(co) && "centered_label" %in% rownames(co) &&
    all(is.finite(co["centered_label", intersect(c("Estimate", "Std. Error"), colnames(co))]))
  est <- se <- ddf <- tv <- pv <- NA_real_
  if (fixed_ok) {
    est <- unname(co["centered_label", "Estimate"]); se <- unname(co["centered_label", "Std. Error"])
    if ("df" %in% colnames(co)) ddf <- unname(co["centered_label", "df"])
    tc <- grep("t value", colnames(co), value = TRUE)[1]; pc <- grep("Pr\\(", colnames(co), value = TRUE)[1]
    if (!is.na(tc)) tv <- unname(co["centered_label", tc]); if (!is.na(pc)) pv <- unname(co["centered_label", pc])
  }
  vc <- as.data.frame(VarCorr(fit))
  getv <- function(v) { z <- vc$vcov[vc$var1 == v & is.na(vc$var2)]; if (length(z)) z[1] else NA_real_ }
  ri <- getv("(Intercept)"); rs <- getv("centered_label"); rv <- vc$vcov[vc$grp == "Residual"][1]
  corr_idx <- !is.na(vc$var2) & ((vc$var1 == "(Intercept)" & vc$var2 == "centered_label") |
    (vc$var1 == "centered_label" & vc$var2 == "(Intercept)"))
  corr <- if (any(corr_idx)) vc$sdcor[corr_idx][1] else NA_real_
  ri_ratio <- if (is.finite(ri)) sqrt(max(ri,0)/rv) else NA_real_
  rs_ratio <- if (is.finite(rs)) sqrt(max(rs,0)/rv) else NA_real_
  G <- NULL
  if (structure == "A" && all(is.finite(c(ri,rs,corr)))) {
    cv <- corr*sqrt(max(ri,0)*max(rs,0)); G <- matrix(c(ri,cv,cv,rs),2,2)
  } else if (structure == "B" && all(is.finite(c(ri,rs)))) G <- diag(c(ri,rs),2)
  else if (structure == "C" && is.finite(rs)) G <- matrix(rs,1,1)
  else if (structure == "D" && is.finite(ri)) G <- matrix(ri,1,1)
  eig <- if (!is.null(G)) eigen((G+t(G))/2,symmetric=TRUE,only.values=TRUE)$values else NA_real_
  emin <- if (all(is.finite(eig))) min(eig) else NA_real_; emax <- if (all(is.finite(eig))) max(eig) else NA_real_
  rp <- tryCatch(rePCA(fit), error=function(e) NULL)
  rpsd <- if (!is.null(rp)) unlist(lapply(rp, function(x)x$sd),use.names=FALSE) else NA_real_
  rpsd <- rpsd[is.finite(rpsd)]; rpmin <- if(length(rpsd))min(rpsd) else NA_real_; rpmax <- if(length(rpsd))max(rpsd) else NA_real_
  rprank <- if(length(rpsd)&&rpmax>0)sum(rpsd/rpmax>tol_singular) else 0L
  der <- fit@optinfo$derivs; heig <- NA_real_; sgr <- NA_real_
  if (!is.null(der$Hessian)) {
    heig <- tryCatch(eigen((der$Hessian+t(der$Hessian))/2,symmetric=TRUE,only.values=TRUE)$values,error=function(e)NA_real_)
    if (!is.null(der$gradient)) sgr <- tryCatch(max(abs(solve(chol(der$Hessian),der$gradient))),error=function(e)NA_real_)
  }
  hmin <- if(all(is.finite(heig)))min(heig) else NA_real_; hmax <- if(all(is.finite(heig)))max(heig) else NA_real_
  hpos <- if(is.finite(hmin)&&is.finite(hmax)) hmin>max(1e-8,1e-8*abs(hmax)) else NA
  conv_ok <- !opt_bad && !nzchar(msgs) && fixed_ok && !isFALSE(hpos)
  list(fit=fit, stats=data.frame(
    structure=structure,optimizer=optimizer,fit_error="",warnings=paste(unique(warnings),collapse=" | "),
    optimizer_code=code,convergence_messages=msgs,convergence_ok=conv_ok,fixed_effect_finite=fixed_ok,
    slope=est,SE=se,df=ddf,t=tv,P_raw=pv,CI_lower=est-1.96*se,CI_upper=est+1.96*se,
    random_intercept_variance=ri,random_slope_variance=rs,intercept_slope_correlation=corr,residual_variance=rv,
    RI_SD_ratio=ri_ratio,RS_SD_ratio_per_label=rs_ratio,RS_SD_ratio_observed_span=6*rs_ratio,
    covariance_eigen_min=emin,covariance_eigen_max=emax,covariance_eigen_ratio=safe_ratio(emin,emax),
    rePCA_min_SD=rpmin,rePCA_max_SD=rpmax,rePCA_SD_ratio=safe_ratio(rpmin,rpmax),rePCA_rank=rprank,
    hessian_eigen_min=hmin,hessian_eigen_max=hmax,hessian_positive=hpos,max_abs_scaled_gradient=sgr,
    isSingular_tol_1e_4=isSingular(fit,tol=tol_singular),stringsAsFactors=FALSE))
}

select_primary <- function(rows) {
  valid <- rows$convergence_ok & rows$fixed_effect_finite & is.finite(rows$slope) & is.finite(rows$SE)
  idx <- which(valid)[1]
  if (!length(idx)) idx <- which(is.finite(rows$slope) & is.finite(rows$SE))[1]
  if (!length(idx)) idx <- 1L
  finite <- valid
  stable_dir <- sum(finite)>=2L && length(unique(dir_sign(rows$slope[finite])))==1L
  medse <- if(sum(finite))median(rows$SE[finite]) else NA_real_
  sr <- if(sum(finite)>=2L)diff(range(rows$slope[finite]))/medse else NA_real_
  ser <- if(sum(finite)>=2L&&min(rows$SE[finite])>0)max(rows$SE[finite])/min(rows$SE[finite]) else NA_real_
  stable <- sum(finite)>=2L && stable_dir && is.finite(sr)&&sr<=.10 && is.finite(ser)&&ser<=1.10
  list(primary=rows[idx,,drop=FALSE], n_valid=sum(valid), direction_consistent=stable_dir,
    slope_range_over_median_SE=sr, SE_max_min_ratio=ser, stable=stable)
}

fit_cr2 <- function(dat) {
  fit <- lm(y ~ centered_label + subject, data=dat)
  ct <- tryCatch(clubSandwich::coef_test(fit,vcov="CR2",cluster=dat$subject,test="Satterthwaite"),error=function(e)e)
  if(inherits(ct,"error")) return(data.frame(CR2_slope=NA_real_,CR2_SE=NA_real_,CR2_df=NA_real_,CR2_t=NA_real_,
    CR2_P=NA_real_,CR2_CI_lower=NA_real_,CR2_CI_upper=NA_real_,CR2_error=conditionMessage(ct)))
  z <- ct[rownames(ct)=="centered_label",,drop=FALSE]
  est<-z$beta; se<-z$SE; df<-z$df_Satt; crit<-qt(.975,df)
  data.frame(CR2_slope=est,CR2_SE=se,CR2_df=df,CR2_t=z$tstat,CR2_P=z$p_Satt,
    CR2_CI_lower=est-crit*se,CR2_CI_upper=est+crit*se,CR2_error="")
}

process_task <- function(task) {
  layer <- task$layer; dataset <- task$dataset; gm <- task$gradient_model; feature <- task$feature
  inp <- obj$inputs[[paste(layer,dataset,sep="_")]]
  mat <- if(gm=="M3")inp$raw else inp$clr
  y <- as.numeric(mat[match(feature,rownames(mat)),])
  dat <- data.frame(y=y,centered_label=as.numeric(inp$pairs$centered_label),subject=factor(inp$pairs$subject))
  meta <- data.frame(task_id=task$task_id,layer=layer,dataset=dataset,gradient_model=gm,feature=feature,
    n_subjects=nlevels(dat$subject),n_pairs=nrow(dat),n_unique_labels=length(unique(dat$centered_label)),
    response_SD=sd(y),missing_y=sum(!is.finite(y)),stringsAsFactors=FALSE)

  fitsB <- lapply(optimizers,function(o)fit_once(dat,"B",o))
  rowsB <- do.call(rbind,lapply(fitsB,function(x)x$stats)); selB<-select_primary(rowsB); B<-selB$primary
  ri_zero <- is.finite(B$RI_SD_ratio)&&B$RI_SD_ratio<=tol_singular
  rs_zero <- is.finite(B$RS_SD_ratio_per_label)&&B$RS_SD_ratio_per_label<=tol_singular

  fitA <- fit_once(dat,"A","nloptwrap"); A<-fitA$stats
  rowsC <- NULL; selC <- NULL; C <- NULL; c_confirmed <- FALSE; bc_slope_norm <- bc_se_ratio <- NA_real_
  if(selB$stable && ri_zero && !rs_zero) {
    fitsC<-lapply(optimizers,function(o)fit_once(dat,"C",o)); rowsC<-do.call(rbind,lapply(fitsC,function(x)x$stats)); selC<-select_primary(rowsC); C<-selC$primary
    bc_slope_norm<-abs(B$slope-C$slope)/median(c(B$SE,C$SE)); bc_se_ratio<-max(B$SE,C$SE)/min(B$SE,C$SE)
    c_confirmed<-selC$stable && C$convergence_ok && !isTRUE(C$isSingular_tol_1e_4) &&
      is.finite(C$RS_SD_ratio_per_label)&&C$RS_SD_ratio_per_label>tol_singular &&
      is.finite(bc_slope_norm)&&bc_slope_norm<=.10 && is.finite(bc_se_ratio)&&bc_se_ratio<=1.10
  }
  D <- NULL
  if(selB$stable && rs_zero) D<-fit_once(dat,"D","nloptwrap")$stats

  if(!any(is.finite(rowsB$slope))) status<-"MODEL_NOT_EVALUABLE"
  else if(!selB$stable) status<-"NUMERICAL_UNSTABLE"
  else if(ri_zero&&rs_zero) status<-"BOTH_VARIANCE_BOUNDARY_UNRESOLVED"
  else if(rs_zero) status<-"RS_BOUNDARY"
  else if(ri_zero&&c_confirmed) status<-"RI_BOUNDARY_C_CONFIRMED"
  else if(ri_zero&&!c_confirmed) status<-"RI_BOUNDARY_C_UNCONFIRMED"
  else status<-"REGULAR_RANDOM_SLOPE"
  rs_valid <- status %in% c("REGULAR_RANDOM_SLOPE","RI_BOUNDARY_C_CONFIRMED")
  numerical_evaluable <- !status %in% c("MODEL_NOT_EVALUABLE","NUMERICAL_UNSTABLE")
  official_P <- if(numerical_evaluable)B$P_raw else NA_real_

  cr2<-fit_cr2(dat)
  result<-cbind(meta,data.frame(
    model_status=status,random_slope_validated=rs_valid,numerical_evaluable=numerical_evaluable,
    selected_optimizer=B$optimizer,optimizer_valid_n=selB$n_valid,optimizer_direction_consistent=selB$direction_consistent,
    optimizer_slope_range_over_median_SE=selB$slope_range_over_median_SE,optimizer_SE_max_min_ratio=selB$SE_max_min_ratio,
    optimizer_stable=selB$stable,slope=B$slope,SE=B$SE,df=B$df,t=B$t,P_raw=B$P_raw,P=official_P,
    CI_lower=B$CI_lower,CI_upper=B$CI_upper,random_intercept_variance=B$random_intercept_variance,
    random_slope_variance=B$random_slope_variance,residual_variance=B$residual_variance,
    RI_SD_ratio=B$RI_SD_ratio,RS_SD_ratio_per_label=B$RS_SD_ratio_per_label,
    RS_SD_ratio_observed_span=B$RS_SD_ratio_observed_span,isSingular_tol_1e_4=B$isSingular_tol_1e_4,
    rePCA_SD_ratio=B$rePCA_SD_ratio,rePCA_rank=B$rePCA_rank,covariance_eigen_min=B$covariance_eigen_min,
    covariance_eigen_max=B$covariance_eigen_max,max_abs_scaled_gradient=B$max_abs_scaled_gradient,
    A_slope=A$slope,A_SE=A$SE,A_P_raw=A$P_raw,A_singular=A$isSingular_tol_1e_4,
    A_intercept_slope_correlation=A$intercept_slope_correlation,A_convergence_ok=A$convergence_ok,
    C_fitted=!is.null(C),C_slope=if(is.null(C))NA_real_ else C$slope,C_SE=if(is.null(C))NA_real_ else C$SE,
    C_P_raw=if(is.null(C))NA_real_ else C$P_raw,C_singular=if(is.null(C))NA else C$isSingular_tol_1e_4,
    C_optimizer_stable=if(is.null(selC))NA else selC$stable,C_confirmed=c_confirmed,
    B_C_slope_diff_over_median_SE=bc_slope_norm,B_C_SE_ratio=bc_se_ratio,
    D_fitted=!is.null(D),D_slope=if(is.null(D))NA_real_ else D$slope,D_SE=if(is.null(D))NA_real_ else D$SE,
    D_P_raw=if(is.null(D))NA_real_ else D$P_raw,stringsAsFactors=FALSE),cr2)
  optrows<-rbind(cbind(meta,rowsB),cbind(meta,A))
  if(!is.null(rowsC))optrows<-rbind(optrows,cbind(meta,rowsC))
  if(!is.null(D))optrows<-rbind(optrows,cbind(meta,D))
  repfit<-NULL
  if(feature%in%representatives[[layer]]) repfit<-list(B= fitsB[[match(B$optimizer,optimizers)]]$fit,
    A=fitA$fit,C=if(!is.null(C))fitsC[[match(C$optimizer,optimizers)]]$fit else NULL,D=NULL,meta=meta,status=status)
  list(result=result,optimizer=optrows,representative_fit=repfit)
}

all_chunks <- sort(unique(tasks$chunk_id))
my_chunks <- all_chunks[((all_chunks-1L)%%n_workers)+1L==worker_id]
cat("worker",worker_id,"of",n_workers,"chunks",length(my_chunks),"\n")
for(ch in my_chunks) {
  dest<-file.path(checkpoint_dir,sprintf("chunk_%03d.rds",ch))
  if(file.exists(dest)) { cat(format(Sys.time()),"skip",ch,"\n"); next }
  tt<-tasks[tasks$chunk_id==ch,,drop=FALSE]
  ans<-vector("list",nrow(tt))
  for(i in seq_len(nrow(tt))) {
    ans[[i]]<-tryCatch(process_task(tt[i,]),error=function(e)list(error=conditionMessage(e),task=tt[i,]))
  }
  tmp<-paste0(dest,".tmp_worker",worker_id)
  saveRDS(ans,tmp,compress=TRUE)
  if(!file.rename(tmp,dest))stop("Could not finalize checkpoint ",dest)
  cat(format(Sys.time()),"completed chunk",ch,"tasks",nrow(tt),"\n")
}
cat("WORKER_COMPLETE",worker_id,"\n")

