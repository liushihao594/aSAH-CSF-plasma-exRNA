options(stringsAsFactors = FALSE, warn = 1)

stage1 <- "revision_stage1"
out <- file.path(stage1, "stage1b_full_refit")
dir.create(file.path(out,"plots"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"intermediate"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"logs"),recursive=TRUE,showWarnings=FALSE)

write_csv <- function(x,path) write.csv(x,path,row.names=FALSE,na="",fileEncoding="UTF-8")
direction <- function(x,tol=sqrt(.Machine$double.eps))ifelse(is.na(x),NA_character_,ifelse(x>tol,"positive",ifelse(x< -tol,"negative","zero")))
bh_full <- function(p,nfam){z<-rep(NA_real_,length(p));ok<-is.finite(p);z[ok]<-p.adjust(p[ok],method="BH",n=nfam);z}

fs<-sort(list.files(file.path(out,"checkpoints"),pattern="^chunk_[0-9]+\\.rds$",full.names=TRUE))
if(length(fs)!=64L)stop("Expected 64 checkpoints, found ",length(fs))
results<-list(); opts<-list(); repfits<-list(); errors<-list(); nr<-no<-ne<-0L
for(f in fs){
  x<-readRDS(f)
  for(z in x){
    if(!is.null(z$error)){ne<-ne+1L;errors[[ne]]<-data.frame(task_id=z$task$task_id,error=z$error)} else {
      nr<-nr+1L;results[[nr]]<-z$result
      no<-no+1L;opts[[no]]<-z$optimizer
      if(!is.null(z$representative_fit))repfits[[paste(z$result$layer,z$result$dataset,z$result$gradient_model,z$result$feature,sep="::")]]<-z$representative_fit
    }
  }
}
if(ne>0L)write_csv(do.call(rbind,errors),file.path(out,"intermediate","TASK_ERRORS.csv"))
if(ne>0L)stop("Task errors present: ",ne)
res<-do.call(rbind,results); rownames(res)<-NULL
opt<-do.call(rbind,opts); rownames(opt)<-NULL
if(nrow(res)!=6344L)stop("Expected 6344 result rows, got ",nrow(res))

for(layer in c("miRNA","GENCODE"))for(dataset in c("FULL","QC_PASS_ONLY"))for(gm in c("M3","M4")){
  ix<-res$layer==layer&res$dataset==dataset&res$gradient_model==gm
  nfam<-if(layer=="miRNA")507L else 1079L
  res$FDR[ix]<-bh_full(res$P[ix],nfam)
  res$CR2_FDR[ix]<-bh_full(res$CR2_P[ix],nfam)
}
res$slope_direction<-direction(res$slope);res$CR2_direction<-direction(res$CR2_slope)
res$CR2_same_direction_as_ModelB<-res$slope_direction==res$CR2_direction
res$CR2_effect_difference<-res$CR2_slope-res$slope
res$CR2_significant<-res$CR2_FDR<.05
res$ModelB_significant<-res$FDR<.05

write_csv(res,file.path(out,"STAGE1B_MODEL_COMPARISON.csv"))
write_csv(res[,c("task_id","layer","dataset","gradient_model","feature","n_subjects","n_pairs","model_status",
  "random_slope_validated","numerical_evaluable","selected_optimizer","optimizer_valid_n","optimizer_direction_consistent",
  "optimizer_slope_range_over_median_SE","optimizer_SE_max_min_ratio","optimizer_stable","random_intercept_variance",
  "random_slope_variance","residual_variance","RI_SD_ratio","RS_SD_ratio_per_label","RS_SD_ratio_observed_span",
  "isSingular_tol_1e_4","rePCA_SD_ratio","rePCA_rank","covariance_eigen_min","covariance_eigen_max",
  "max_abs_scaled_gradient","A_singular","A_intercept_slope_correlation","A_convergence_ok","C_fitted","C_singular",
  "C_optimizer_stable","C_confirmed","B_C_slope_diff_over_median_SE","B_C_SE_ratio","D_fitted")],
  file.path(out,"STAGE1B_MODEL_BOUNDARY_DIAGNOSTICS.csv"))
write_csv(res[,c("task_id","layer","dataset","gradient_model","feature","CR2_slope","CR2_SE","CR2_df","CR2_t",
  "CR2_P","CR2_FDR","CR2_CI_lower","CR2_CI_upper","CR2_error","CR2_direction","CR2_same_direction_as_ModelB",
  "CR2_effect_difference","CR2_significant")],file.path(out,"STAGE1B_CR2_RESULTS.csv"))
write_csv(opt,file.path(out,"intermediate","STAGE1B_ALL_OPTIMIZER_DIAGNOSTICS.csv"))
saveRDS(repfits,file.path(out,"intermediate","STAGE1B_REPRESENTATIVE_MODEL_OBJECTS.rds"),compress=TRUE)

obj<-readRDS(file.path(stage1,"intermediate","prepared_stage1_inputs.rds"))

pick_long <- function(layer,dataset,gm){
  z<-res[res$layer==layer&res$dataset==dataset&res$gradient_model==gm,,drop=FALSE]
  z[match(obj$limma[[paste(layer,dataset,"M1",sep="_")]]$feature,z$feature),,drop=FALSE]
}
prefix_df <- function(z,pre,keep="feature"){
  names(z)[!names(z)%in%keep]<-paste0(pre,names(z)[!names(z)%in%keep]);z
}
limma_part <- function(layer,dataset,model,pre){
  z<-obj$limma[[paste(layer,dataset,model,sep="_")]][,c("feature","slope","SE","CI_lower","CI_upper","P","FDR","duplicate_correlation")]
  prefix_df(z,pre)
}
mixed_part <- function(layer,dataset,gm,pre){
  z<-pick_long(layer,dataset,gm)
  keep<-c("feature","model_status","random_slope_validated","numerical_evaluable","selected_optimizer","optimizer_valid_n",
    "optimizer_direction_consistent","optimizer_slope_range_over_median_SE","optimizer_SE_max_min_ratio","optimizer_stable",
    "slope","SE","df","t","P_raw","P","FDR","CI_lower","CI_upper","random_intercept_variance","random_slope_variance",
    "residual_variance","RI_SD_ratio","RS_SD_ratio_per_label","RS_SD_ratio_observed_span","isSingular_tol_1e_4",
    "C_fitted","C_confirmed","B_C_slope_diff_over_median_SE","B_C_SE_ratio","CR2_slope","CR2_SE","CR2_df","CR2_t",
    "CR2_P","CR2_FDR","CR2_CI_lower","CR2_CI_upper","CR2_same_direction_as_ModelB")
  prefix_df(z[,keep],pre)
}
join_feature <- function(lst){Reduce(function(x,y)merge(x,y,by="feature",all=TRUE,sort=FALSE),lst)}

full_tables<-list();qc_tables<-list();class_tables<-list()
for(layer in c("miRNA","GENCODE")){
  full<-join_feature(list(limma_part(layer,"FULL","M1","M1_"),limma_part(layer,"FULL","M2","M2_"),
    mixed_part(layer,"FULL","M3","M3_"),mixed_part(layer,"FULL","M4","M4_")))
  qc<-join_feature(list(limma_part(layer,"QC_PASS_ONLY","M1","M1_"),limma_part(layer,"QC_PASS_ONLY","M2","M2_"),
    mixed_part(layer,"QC_PASS_ONLY","M3","M3_"),mixed_part(layer,"QC_PASS_ONLY","M4","M4_")))
  full$layer<-layer;qc$layer<-layer
  sf<-obj$sensitivities[[paste(layer,"FULL",sep="_")]]
  sf<-sf[,setdiff(names(sf),c("layer","dataset"))]
  names(sf)[names(sf)!="feature"]<-paste0("M2_reference_sensitivity_",names(sf)[names(sf)!="feature"])
  full<-merge(full,sf,by="feature",all.x=TRUE,sort=FALSE)
  full<-full[match(obj$limma[[paste(layer,"FULL","M1",sep="_")]]$feature,full$feature),]
  qc<-qc[match(obj$limma[[paste(layer,"QC_PASS_ONLY","M1",sep="_")]]$feature,qc$feature),]

  qcm<-qc[,c("feature","M3_slope","M3_FDR","M3_model_status","M3_random_slope_validated","M3_numerical_evaluable",
    "M3_CR2_slope","M3_CR2_FDR","M4_slope","M4_FDR","M4_model_status","M4_random_slope_validated",
    "M4_numerical_evaluable","M4_CR2_slope","M4_CR2_FDR")]
  names(qcm)[-1]<-paste0("QC_",names(qcm)[-1])
  comb<-merge(full,qcm,by="feature",all.x=TRUE,sort=FALSE)
  comb<-comb[match(full$feature,comb$feature),]

  d1<-direction(comb$M1_slope);d3<-direction(comb$M3_slope);d4<-direction(comb$M4_slope)
  qd3<-direction(comb$QC_M3_slope);qd4<-direction(comb$QC_M4_slope)
  m1sig<-comb$M1_FDR<.05;m2sig<-comb$M2_FDR<.05;m3sig<-comb$M3_FDR<.05;m4sig<-comb$M4_FDR<.05
  cr3sig<-comb$M3_CR2_FDR<.05;cr4sig<-comb$M4_CR2_FDR<.05
  dirfull<-d1==d3&d1==d4
  qc_dir<-qd3==d3&qd4==d4
  qcd3<-abs(comb$QC_M3_slope-comb$M3_slope)>pmax(.5*abs(comb$M3_slope),.1)
  qcd4<-abs(comb$QC_M4_slope-comb$M4_slope)>pmax(.5*abs(comb$M4_slope),.1)
  qc_dir_change<-qd3!=d3|qd4!=d4
  qc_stability_drop<-(comb$M3_random_slope_validated&!comb$QC_M3_random_slope_validated)|
    (comb$M4_random_slope_validated&!comb$QC_M4_random_slope_validated)|
    !comb$QC_M3_numerical_evaluable|!comb$QC_M4_numerical_evaluable
  qc_sensitive<-qc_dir_change|qcd3|qcd4|qc_stability_drop
  qc_fdr_cross_only<-((m3sig!=(comb$QC_M3_FDR<.05))|(m4sig!=(comb$QC_M4_FDR<.05)))&!qc_sensitive
  comp_effect<-abs(comb$M4_slope-comb$M3_slope)>pmax(.5*abs(comb$M3_slope),.1)
  comp_sensitive<-d3!=d4|comp_effect|(m3sig!=m4sig)
  strict<-m1sig&m3sig&m4sig&dirfull&comb$M3_random_slope_validated&comb$M4_random_slope_validated&
    comb$QC_M3_random_slope_validated&comb$QC_M4_random_slope_validated&qc_dir&!qc_sensitive
  any_rs_boundary<-comb$M3_model_status=="RS_BOUNDARY"|comb$M4_model_status=="RS_BOUNDARY"
  cr_dir<-direction(comb$M3_CR2_slope)==d3&direction(comb$M4_CR2_slope)==d4
  boundary_supported<-!strict&m1sig&m3sig&m4sig&dirfull&any_rs_boundary&cr3sig&cr4sig&cr_dir&
    comb$M3_numerical_evaluable&comb$M4_numerical_evaluable&comb$QC_M3_numerical_evaluable&comb$QC_M4_numerical_evaluable&qc_dir&!qc_sensitive
  unresolved_status<-comb$M3_model_status%in%c("MODEL_NOT_EVALUABLE","NUMERICAL_UNSTABLE","BOTH_VARIANCE_BOUNDARY_UNRESOLVED","RI_BOUNDARY_C_UNCONFIRMED")|
    comb$M4_model_status%in%c("MODEL_NOT_EVALUABLE","NUMERICAL_UNSTABLE","BOTH_VARIANCE_BOUNDARY_UNRESOLVED","RI_BOUNDARY_C_UNCONFIRMED")
  unresolved<-unresolved_status|(any_rs_boundary&!boundary_supported)
  category<-ifelse(strict,"STRICT_CORE_ROBUST",ifelse(boundary_supported,"RS_BOUNDARY_BUT_OTHERWISE_SUPPORTED",
    ifelse(unresolved,"UNRESOLVED",ifelse(comp_sensitive,"COMPOSITION_SENSITIVE",ifelse(qc_sensitive,"QC_SENSITIVE","OTHER_SUPPORTED_OR_NULL")))))
  cls<-data.frame(feature=comb$feature,layer=layer,M1_significant=m1sig,M2_significant=m2sig,M3_significant=m3sig,M4_significant=m4sig,
    M3_M4_both_significant=m3sig&m4sig,M1_M3_M4_direction_consistent=dirfull,
    class_M3_random_slope_validated=comb$M3_random_slope_validated,class_M4_random_slope_validated=comb$M4_random_slope_validated,
    QC_M3_significant=comb$QC_M3_FDR<.05,QC_M4_significant=comb$QC_M4_FDR<.05,QC_direction_consistent=qc_dir,
    QC_effect_change_M3=qcd3,QC_effect_change_M4=qcd4,QC_model_stability_decline=qc_stability_drop,
    QC_sensitive=qc_sensitive,QC_FDR_crossing_only=qc_fdr_cross_only,composition_sensitive=comp_sensitive,
    CR2_M3_significant=cr3sig,CR2_M4_significant=cr4sig,CR2_both_supported=cr3sig&cr4sig&cr_dir,
    any_RS_boundary=any_rs_boundary,strict_core_robust=strict,boundary_supported=boundary_supported,
    unresolved=unresolved,primary_robustness_category=category,stringsAsFactors=FALSE)
  comb<-merge(comb,cls,by=c("feature","layer"),all.x=TRUE,sort=FALSE)
  comb<-comb[match(full$feature,comb$feature),]
  full_tables[[layer]]<-comb;qc_tables[[layer]]<-qc;class_tables[[layer]]<-cls
}

write_csv(full_tables$miRNA,file.path(out,"STAGE1B_FULL_MODEL_RESULTS_miRNA.csv"))
write_csv(full_tables$GENCODE,file.path(out,"STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv"))
write_csv(do.call(rbind,qc_tables),file.path(out,"STAGE1B_QC_MODEL_RESULTS.csv"))
classes<-do.call(rbind,class_tables)
strict<-classes[classes$strict_core_robust,]
boundary<-classes[classes$boundary_supported,]
write_csv(strict,file.path(out,"STAGE1B_STRICT_CORE_ROBUST.csv"))
write_csv(boundary,file.path(out,"STAGE1B_BOUNDARY_SUPPORTED.csv"))

original_cross<-classes[classes$M1_significant,]
write_csv(original_cross,file.path(out,"intermediate","STAGE1B_ORIGINAL_DISCOVERY_CROSSWALK.csv"))
write_csv(classes,file.path(out,"intermediate","STAGE1B_ALL_FEATURE_CLASSIFICATION.csv"))

rep_ids<-list(miRNA=c("hsa-miR-338-3p","hsa-miR-381-3p","hsa-miR-410-3p","hsa-miR-375"),
  GENCODE=c("MT-TL2:Mt_tRNA","RN7SL1:antisense","IL8:protein_coding","HBA1:protein_coding","TMSB4X:protein_coding","EIF4A1:retained_intron"))
reps<-do.call(rbind,lapply(names(rep_ids),function(layer){
  z<-full_tables[[layer]][match(rep_ids[[layer]],full_tables[[layer]]$feature),]
  status<-ifelse(z$strict_core_robust,"RETAIN_STRICT_CORE",ifelse(z$boundary_supported,"RETAIN_WITH_RS_BOUNDARY_QUALIFICATION",
    ifelse(z$M3_significant&z$M4_significant&z$M1_M3_M4_direction_consistent&!z$unresolved,"RETAIN_WITH_MODEL_OR_QC_QUALIFICATION",
      "NOT_RECOMMENDED_AS_PRIMARY_REPRESENTATIVE")))
  data.frame(z[,c("feature","layer","M1_slope","M1_SE","M1_P","M1_FDR","M2_slope","M2_SE","M2_P","M2_FDR",
    "M3_slope","M3_SE","M3_CI_lower","M3_CI_upper","M3_P","M3_FDR","M3_model_status","M3_random_slope_validated",
    "M3_CR2_slope","M3_CR2_SE","M3_CR2_P","M3_CR2_FDR","M4_slope","M4_SE","M4_CI_lower","M4_CI_upper",
    "M4_P","M4_FDR","M4_model_status","M4_random_slope_validated","M4_CR2_slope","M4_CR2_SE","M4_CR2_P","M4_CR2_FDR",
    "QC_M3_slope","QC_M3_FDR","QC_M3_model_status","QC_M4_slope","QC_M4_FDR","QC_M4_model_status",
    "strict_core_robust","boundary_supported","composition_sensitive","QC_sensitive","unresolved")],
    representative_recommendation=status,stringsAsFactors=FALSE)
}))
write_csv(reps,file.path(out,"STAGE1B_REPRESENTATIVE_FEATURES.csv"))

# Summaries.
summary_rows<-list();k<-0L
for(layer in c("miRNA","GENCODE"))for(dataset in c("FULL","QC_PASS_ONLY"))for(gm in c("M3","M4")){
  z<-res[res$layer==layer&res$dataset==dataset&res$gradient_model==gm,]
  k<-k+1L;summary_rows[[k]]<-data.frame(layer=layer,dataset=dataset,model=gm,n_features=nrow(z),
    n_numerical_evaluable=sum(z$numerical_evaluable),n_random_slope_validated=sum(z$random_slope_validated),
    n_regular=sum(z$model_status=="REGULAR_RANDOM_SLOPE"),n_RI_boundary_C_confirmed=sum(z$model_status=="RI_BOUNDARY_C_CONFIRMED"),
    n_RS_boundary=sum(z$model_status=="RS_BOUNDARY"),n_both_variance_boundary=sum(z$model_status=="BOTH_VARIANCE_BOUNDARY_UNRESOLVED"),
    n_numerical_unstable=sum(z$model_status=="NUMERICAL_UNSTABLE"),n_not_evaluable=sum(z$model_status=="MODEL_NOT_EVALUABLE"),
    n_significant=sum(z$FDR<.05,na.rm=TRUE),significant_negative=sum(z$FDR<.05&z$slope<0,na.rm=TRUE),
    significant_positive=sum(z$FDR<.05&z$slope>0,na.rm=TRUE),n_CR2_significant=sum(z$CR2_FDR<.05,na.rm=TRUE),
    ModelB_CR2_slope_cor=cor(z$slope,z$CR2_slope,use="complete.obs"),ModelB_CR2_direction_same=sum(z$CR2_same_direction_as_ModelB,na.rm=TRUE),
    stringsAsFactors=FALSE)
}
summary_counts<-do.call(rbind,summary_rows);write_csv(summary_counts,file.path(out,"intermediate","STAGE1B_MODEL_SUMMARY_COUNTS.csv"))
cat_counts<-as.data.frame(table(classes$layer,classes$primary_robustness_category),stringsAsFactors=FALSE)
names(cat_counts)<-c("layer","category","n");cat_counts<-cat_counts[cat_counts$n>0,]
write_csv(cat_counts,file.path(out,"intermediate","STAGE1B_ROBUST_CATEGORY_COUNTS.csv"))

old_summary<-do.call(rbind,lapply(c("miRNA","GENCODE"),function(layer){
  z<-full_tables[[layer]];o<-z$M1_FDR<.05
  data.frame(layer=layer,original_n=sum(o),M2_significant=sum(o&z$M2_FDR<.05),M3_significant=sum(o&z$M3_FDR<.05,na.rm=TRUE),
    M4_significant=sum(o&z$M4_FDR<.05,na.rm=TRUE),M3_M4_both_significant=sum(o&z$M3_FDR<.05&z$M4_FDR<.05,na.rm=TRUE),
    M1_M3_M4_direction_consistent=sum(o&z$M1_M3_M4_direction_consistent,na.rm=TRUE),
    strict_core=sum(o&z$strict_core_robust),boundary_supported=sum(o&z$boundary_supported),
    RS_boundary=sum(o&z$any_RS_boundary),QC_sensitive=sum(o&z$QC_sensitive),CR2_both_supported=sum(o&z$CR2_both_supported),
    unresolved=sum(o&z$unresolved),stringsAsFactors=FALSE)
}))
write_csv(old_summary,file.path(out,"intermediate","STAGE1B_ORIGINAL_DISCOVERY_SUMMARY.csv"))

# Diagnostic plots, not manuscript figures.
png(file.path(out,"plots","MODEL_B_VS_CR2_SLOPES.png"),width=2200,height=2100,res=180)
par(mfrow=c(2,2),mar=c(4.5,4.5,3,1))
for(layer in c("miRNA","GENCODE"))for(gm in c("M3","M4")){
 z<-res[res$layer==layer&res$dataset=="FULL"&res$gradient_model==gm,]
 plot(z$slope,z$CR2_slope,pch=16,cex=.45,xlab="Model B slope",ylab="CR2 slope",main=paste(layer,gm));abline(0,1,col="red")
}
dev.off()
png(file.path(out,"plots","FULL_VS_QC_MODEL_B_SLOPES.png"),width=2200,height=2100,res=180)
par(mfrow=c(2,2),mar=c(4.5,4.5,3,1))
for(layer in c("miRNA","GENCODE"))for(gm in c("M3","M4")){
 a<-res[res$layer==layer&res$dataset=="FULL"&res$gradient_model==gm,];b<-res[res$layer==layer&res$dataset=="QC_PASS_ONLY"&res$gradient_model==gm,];b<-b[match(a$feature,b$feature),]
 plot(a$slope,b$slope,pch=16,cex=.45,xlab="FULL slope",ylab="QC PASS-only slope",main=paste(layer,gm));abline(0,1,col="red")
}
dev.off()
png(file.path(out,"plots","MODEL_STATUS_COUNTS.png"),width=2400,height=1400,res=180)
tt<-xtabs(~interaction(layer,dataset,gradient_model,sep=" ")+model_status,res)
barplot(t(tt),beside=FALSE,las=2,cex.names=.75,ylab="Features",legend.text=TRUE,args.legend=list(cex=.7,x="topright"))
dev.off()

saveRDS(list(results_long=res,full_tables=full_tables,qc_tables=qc_tables,classifications=classes,
  strict_core=strict,boundary_supported=boundary,representatives=reps,summary_counts=summary_counts,
  old_summary=old_summary),file.path(out,"intermediate","STAGE1B_COMPACT_RESULTS.rds"),compress=TRUE)

capture.output({
 cat("Stage 1B finalized: ",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),"\n",sep="")
 cat("Tasks: 6344; checkpoints: 64; task errors: 0\n")
 cat("isSingular tolerance: 1e-4\n")
 for(p in c("lme4","lmerTest","clubSandwich","limma","statmod"))cat(p,if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else"NOT_AVAILABLE","\n")
 cat("\n");print(sessionInfo())
},file=file.path(out,"logs","STAGE1B_sessionInfo.txt"))

cat("STAGE1B_FINALIZE_COMPLETE\n")
cat("strict",nrow(strict),"boundary_supported",nrow(boundary),"\n")

