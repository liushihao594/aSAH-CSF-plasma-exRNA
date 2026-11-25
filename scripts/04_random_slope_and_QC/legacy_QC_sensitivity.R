# Stage 0C: authorized legacy-model QC sensitivity only; no new random slope/CLR/XDec inference.
options(stringsAsFactors=FALSE)
out <- 'revision_stage0c'
stopifnot(requireNamespace('limma',quietly=TRUE),requireNamespace('statmod',quietly=TRUE))
readCSV <- function(p) read.csv(p,check.names=FALSE,fileEncoding='UTF-8-BOM')
pairs <- readCSV(file.path(out,'data/FULL_153_pairs.csv'))
pass <- tolower(as.character(pairs$both_QC_PASS))=='true'
stopifnot(nrow(pairs)==153,sum(pass)==141,length(unique(pairs$subject))==30)
manifest <- data.frame(item=c('transformation','features','linear_time','correlation','eBayes','FDR','QC','missingness'),setting=c('log2(RPM+1), CSF minus plasma','Frozen original 507 miRNA /1079 GENCODE IDs, no refiltering','label 1-7 centered at D04','duplicateCorrelation re-estimated separately for each cohort/model/refit','trend=TRUE robust=TRUE','BH within each layer and test family, q<0.05','Both members must have native Atlas PASS for sensitivity branch','No imputation, no alteration to original records'))
write.csv(manifest,file.path(out,'statistics/ANALYSIS_SETTINGS.csv'),row.names=FALSE)
fitBase <- function(y,design,subject){
 rho<-limma::duplicateCorrelation(y,design,block=factor(subject))$consensus.correlation
 fit<-limma::lmFit(y,design,block=factor(subject),correlation=rho)
 list(fit=fit,rho=rho)
}
fitLinear <- function(y,p){
 x<-as.integer(sub('D','',p$collection_label))-4
 f<-fitBase(y,model.matrix(~x),p$subject)
 e<-limma::eBayes(f$fit,trend=TRUE,robust=TRUE)
 se<-e$stdev.unscaled[,2]*sqrt(e$s2.post)
 ci<-qt(.975,df=e$df.total)*se
 list(table=data.frame(feature=rownames(y),linear_slope=e$coefficients[,2],linear_SE=se,linear_CI_lower=e$coefficients[,2]-ci,linear_CI_upper=e$coefficients[,2]+ci,linear_t=e$t[,2],linear_df=e$df.total,linear_P=e$p.value[,2],linear_FDR=p.adjust(e$p.value[,2],'BH')),rho=f$rho)
}
fitAll <- function(y,p,layer){
 linear<-fitLinear(y,p)
 label<-factor(as.integer(sub('D','',p$collection_label)),levels=1:7)
 catfit<-fitBase(y,model.matrix(~label),p$subject)
 ca<-limma::eBayes(catfit$fit,trend=TRUE,robust=TRUE)
 categorical<-limma::topTable(ca,coef=2:7,number=Inf,sort.by='none')
 stopifnot(identical(rownames(categorical),rownames(y)))
 avg<-limma::eBayes(limma::contrasts.fit(catfit$fit,matrix(c(1,rep(1/7,6)),ncol=1)),trend=TRUE,robust=TRUE)
 se<-avg$stdev.unscaled[,1]*sqrt(avg$s2.post);ci<-qt(.975,df=avg$df.total)*se
 ref<-setdiff(rownames(y),if(layer=='miRNA')'hsa-miR-486-5p' else 'RNY4:misc_RNA')
 centered<-sweep(y,2,colMeans(y[ref,,drop=FALSE]),'-')
 composition<-fitLinear(centered,p)
 z<-linear$table
 z$mean_gradient<-avg$coefficients[,1];z$mean_SE<-se;z$mean_CI_lower<-z$mean_gradient-ci;z$mean_CI_upper<-z$mean_gradient+ci
 z$mean_P<-avg$p.value[,1];z$mean_FDR<-p.adjust(z$mean_P,'BH')
 z$categorical_F<-categorical$F;z$categorical_P<-categorical$P.Value;z$categorical_FDR<-p.adjust(z$categorical_P,'BH')
 z$composition_slope<-composition$table$linear_slope;z$composition_P<-composition$table$linear_P;z$composition_FDR<-composition$table$linear_FDR
 list(table=z,rho=data.frame(model=c('linear','categorical','legacy_composition_linear'),rho=c(linear$rho,catfit$rho,composition$rho)))
}
allcomparisons<-list();allsummary<-list();allrepro<-list();allreps<-list();allrho<-list();allloo<-list();allsubj<-list();allpatterns<-list()
for(layer in c('miRNA','GENCODE')){
 cat('\nSTART',layer,'\n');flush.console()
 input<-list.files('revision_stage0/source_extract/TRACKB_SOURCE_DATA_PACKAGE/01_GEO_processed_expression',pattern=if(layer=='miRNA')'_miRNA_' else '_gencode_',full.names=TRUE)
 con<-gzfile(input,'rt');lines<-readLines(con);close(con)
 ids<-readLines(paste0('revision_stage0/stage0_',layer,'_retained_features.txt'),warn=FALSE)
 selected<-c(TRUE,sub('\t.*$','',lines[-1])%in%ids)
 rpm<-read.delim(text=paste(lines[selected],collapse='\n'),row.names=1,check.names=FALSE)
 rm(lines);invisible(gc())
 stopifnot(setequal(rownames(rpm),ids))
 y<-log2(as.matrix(rpm[,pairs$CSF_expression_id])+1)-log2(as.matrix(rpm[,pairs$plasma_expression_id])+1)
 colnames(y)<-pairs$pair_id
 for(arm in c('FULL','QC_PASS_ONLY')){
  idx<-if(arm=='FULL')seq_len(nrow(pairs)) else which(pass)
  ys<-y[,idx,drop=FALSE];ps<-pairs[idx,,drop=FALSE]
  saveRDS(list(layer=layer,arm=arm,feature_set='Frozen original features',pairs=ps,gradient=ys,CSF_RPM=as.matrix(rpm[,ps$CSF_expression_id]),plasma_RPM=as.matrix(rpm[,ps$plasma_expression_id])),file.path(out,'data',paste0(layer,'_',arm,'_paired_data.rds')))
  write.csv(data.frame(feature=rownames(ys),ys,check.names=FALSE),file.path(out,'data',paste0(layer,'_',arm,'_gradients.csv')),row.names=FALSE)
  r<-fitAll(ys,ps,layer);r$table$layer<-layer;r$table$arm<-arm;r$table$n_pairs<-nrow(ps);r$table$n_subjects<-length(unique(ps$subject))
  write.csv(r$table,file.path(out,'statistics',paste0(layer,'_',arm,'_complete_results.csv')),row.names=FALSE)
  allrho[[paste(layer,arm)]]<-data.frame(layer=layer,arm=arm,r$rho)
  if(arm=='FULL')full<-r$table else qc<-r$table
  count<-function(q,d=NULL){a<-q<.05;if(is.null(d))sum(a) else sum(a&d)}
  allsummary[[paste(layer,arm)]]<-data.frame(layer=layer,arm=arm,n_features=nrow(ys),n_pairs=ncol(ys),n_subjects=length(unique(ps$subject)),mean_significant=count(r$table$mean_FDR),linear_significant=count(r$table$linear_FDR),linear_negative=count(r$table$linear_FDR,r$table$linear_slope<0),linear_positive=count(r$table$linear_FDR,r$table$linear_slope>0),categorical_significant=count(r$table$categorical_FDR),composition_linear_same_direction_significant=sum(r$table$linear_FDR<.05&r$table$composition_FDR<.05&r$table$linear_slope*r$table$composition_slope>0))
 }
 old<-readCSV(file.path('revision_stage0',paste0('stage0_legacy_refit_',layer,'.csv')))
 old<-old[match(full$feature,old$feature),]
 for(k in c('linear_slope','linear_P','linear_FDR','categorical_P','categorical_FDR'))allrepro[[paste(layer,k)]]<-data.frame(layer=layer,field=k,max_abs_difference=max(abs(full[[k]]-old[[k]])))
 comp<-full;names(comp)[!names(comp)%in%c('feature','layer')]<-paste0('full_',names(comp)[!names(comp)%in%c('feature','layer')])
 q<-qc;names(q)[!names(q)%in%c('feature','layer')]<-paste0('qc_',names(q)[!names(q)%in%c('feature','layer')])
 comp<-merge(comp,q,by=c('feature','layer'),sort=FALSE)
 comp$delta_slope<-comp$qc_linear_slope-comp$full_linear_slope
 comp$relative_slope_change<-abs(comp$delta_slope)/pmax(abs(comp$full_linear_slope),1e-8)
 comp$same_slope_direction<-sign(comp$full_linear_slope)==sign(comp$qc_linear_slope)
 comp$temporal_status<-ifelse(comp$full_linear_FDR<.05,ifelse(comp$qc_linear_FDR<.05,ifelse(comp$same_slope_direction,'SIGNIFICANT_BOTH_SAME_DIRECTION','SIGNIFICANT_BOTH_OPPOSITE_DIRECTION'),'FULL_ONLY_SIGNIFICANT'),ifelse(comp$qc_linear_FDR<.05,'QC_ONLY_SIGNIFICANT','NEITHER_SIGNIFICANT'))
 allcomparisons[[layer]]<-comp
 original<-readCSV(file.path('\u6295\u7a3f\u6750\u6599 - BMJ',if(layer=='miRNA')'Supplementary Data 1 miRNA results.csv' else 'Supplementary Data 2 GENCODE results.csv'))
 repids<-original$feature[tolower(original$representative_flag)=='true'|(!is.na(original$figure_panel)&nzchar(original$figure_panel))]
 reportids<-union(repids,if(layer=='miRNA')c('hsa-miR-486-5p','hsa-miR-451a') else 'RNY4:misc_RNA')
 reps<-comp[comp$feature%in%reportids,];reps$original_representative<-reps$feature%in%repids
 reps$figure_panel<-original$figure_panel[match(reps$feature,original$feature)]
 allreps[[layer]]<-reps
 for(arm in c('FULL','QC_PASS_ONLY')){
  idx<-if(arm=='FULL')seq_len(nrow(pairs)) else which(pass)
  ps<-pairs[idx,,drop=FALSE];ys<-y[,idx,drop=FALSE];base<-if(arm=='FULL')full else qc
  cat('LOO',layer,arm,'\n');flush.console()
  for(subject in sort(unique(ps$subject))){
   keep<-ps$subject!=subject
   lo<-fitLinear(ys[,keep,drop=FALSE],ps[keep,,drop=FALSE])
   rr<-lo$table;rr$layer<-layer;rr$arm<-arm;rr$omitted_subject<-subject;rr$n_pairs<-sum(keep);rr$rho<-lo$rho
   rr$baseline_slope<-base$linear_slope;rr$baseline_FDR<-base$linear_FDR
   rr$direction_same_as_baseline<-sign(rr$linear_slope)==sign(rr$baseline_slope)
   rr$relative_change<-abs(rr$linear_slope-rr$baseline_slope)/pmax(abs(rr$baseline_slope),1e-8)
   allloo[[paste(layer,arm,subject)]]<-rr
   sig<-base$linear_FDR<.05
   allsubj[[paste(layer,arm,subject)]]<-data.frame(layer=layer,arm=arm,omitted_subject=subject,removed_pairs=sum(!keep),slope_Pearson=cor(rr$linear_slope,rr$baseline_slope),baseline_significant=sum(sig),significant_direction_flips=sum(sig&!rr$direction_same_as_baseline),significant_retaining_FDR=sum(sig&rr$linear_FDR<.05),median_relative_change_baseline_significant=median(rr$relative_change[sig]),max_relative_change_baseline_significant=max(rr$relative_change[sig]))
  }
 }
 # Diagnostic only: original temporal class letters are not silently regenerated with new selection rules.
 allpatterns[[layer]]<-data.frame(layer=layer,feature=original$feature,original_temporal_pattern_class=original$temporal_pattern_class,note='Original class recorded for comparison; formal class reselection not performed in Stage 0C')
 cat('DONE',layer,'\n');flush.console()
}
write.csv(do.call(rbind,allcomparisons),file.path(out,'QC_SENSITIVITY_RESULTS.csv'),row.names=FALSE)
write.csv(do.call(rbind,allsummary),file.path(out,'statistics/QC_summary_counts.csv'),row.names=FALSE)
write.csv(do.call(rbind,allrepro),file.path(out,'statistics/FULL_reproduction_check.csv'),row.names=FALSE)
write.csv(do.call(rbind,allreps),file.path(out,'statistics/representative_feature_comparison.csv'),row.names=FALSE)
write.csv(do.call(rbind,allrho),file.path(out,'statistics/model_correlations.csv'),row.names=FALSE)
write.csv(do.call(rbind,allloo),file.path(out,'statistics/LOO_all_feature_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,allsubj),file.path(out,'statistics/LOO_subject_impact_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,allpatterns),file.path(out,'statistics/original_pattern_class_reference.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'logs/QC_sessionInfo.txt'))
print(do.call(rbind,allsummary))
cat('COMPLETE Stage 0C legacy QC sensitivity\n')
