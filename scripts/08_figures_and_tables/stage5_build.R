options(stringsAsFactors = FALSE, warn = 1)
root <- '.'
out <- file.path(root, 'revision_stage5')
dirs <- c('01_results_freeze','02_main_figures','03_supplementary_figures','04_tables',
          '05_graphical_abstract','06_source_code','07_validation','logs')
for (d in dirs) dir.create(file.path(out,d), recursive=TRUE, showWarnings=FALSE)

readc <- function(p) read.csv(file.path(root,p), check.names=FALSE, fileEncoding='UTF-8-BOM')
writec <- function(x,p) write.csv(x,file.path(out,p),row.names=FALSE,na='')
bool <- function(x) isTRUE(x) || identical(toupper(as.character(x)),'TRUE')
safe_num <- function(x) suppressWarnings(as.numeric(x))
`%||%` <- function(a,b) if (length(a) && !all(is.na(a))) a else b

mi <- readc('revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_miRNA.csv')
ge <- readc('revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv')
core <- readc('revision_stage1/stage1b_full_refit/STAGE1B_STRICT_CORE_ROBUST.csv')
boundary <- readc('revision_stage1/stage1b_full_refit/STAGE1B_BOUNDARY_SUPPORTED.csv')
stage2core <- readc('revision_stage2/STAGE2_CORE_FEATURE_ASSESSMENT.csv')
shift <- readc('revision_stage2/STAGE2_GLOBAL_SHIFT_ANALYSIS.csv')
high <- readc('revision_stage2/STAGE2_HIGH_ABUNDANCE_RNA_RESULTS.csv')
biotype <- readc('revision_stage2/STAGE2_BIOTYPE_COMPOSITION.csv')
slopes <- readc('revision_stage3/STAGE3_CSF_PLASMA_SLOPES.csv')
inter <- readc('revision_stage3/STAGE3_BIOFLUID_INTERACTION.csv')
decomp <- readc('revision_stage3/STAGE3_GRADIENT_DECOMPOSITION.csv')
mirprof <- readc('revision_stage3/STAGE3_CORE_MIRNA_PROFILES.csv')
geprof <- readc('revision_stage3/STAGE3_GENCODE_PROFILES.csv')
endpoint <- readc('revision_stage3/STAGE3_ENDPOINT_SENSITIVITY.csv')
repfeat <- readc('revision_stage3/STAGE3_REPRESENTATIVE_FEATURES.csv')
corr <- readc('revision_stage4/STAGE4_CORRELATION_SUMMARY.csv')
scores <- readc('revision_stage4/STAGE4_RAW_CLR_CROSSFITTED_SCORES.csv')
boot <- readc('revision_stage4/STAGE4_PATIENT_BOOTSTRAP.csv')
loo <- readc('revision_stage4/STAGE4_LOO_INFLUENCE.csv')
s12 <- readc('revision_stage4/STAGE4_SUBJECT12_SENSITIVITY.csv')
pairs <- readc('revision_stage0c/data/FULL_153_pairs.csv')
qcpairs <- readc('revision_stage0c/data/QC_AFFECTED_PAIRS.csv')
qcsum <- readc('revision_stage0c/statistics/QC_summary_counts.csv')
mi_pair <- readRDS(file.path(root,'revision_stage0c/data/miRNA_FULL_paired_data.rds'))
ge_pair <- readRDS(file.path(root,'revision_stage0c/data/GENCODE_FULL_paired_data.rds'))

mean_mi_path <- 'revision_stage0d/source_extract/TRACKB_REVISION_TRANSFER_PACKAGE/11_ANALYSIS_RESULTS_AND_AUDITS/Primary_Results/Phase_1B_20260928_205553_Tables/TRACKB_MIRNA_MEAN_GRADIENT_STATS.csv'
mean_ge_path <- 'revision_stage0d/source_extract/TRACKB_REVISION_TRANSFER_PACKAGE/11_ANALYSIS_RESULTS_AND_AUDITS/Primary_Results/Phase_1B_20260928_205553_Tables/TRACKB_GENCODE_MEAN_GRADIENT_STATS.csv'
mean_mi <- readc(mean_mi_path); mean_ge <- readc(mean_ge_path)
mean_mi_comp <- readc('revision_stage0d/source_extract/TRACKB_REVISION_TRANSFER_PACKAGE/11_ANALYSIS_RESULTS_AND_AUDITS/Primary_Results/Phase_1B_20260928_205553_Tables/TRACKB_MIRNA_SENS_composition_MEAN_GRADIENT_STATS.csv')
mean_ge_comp <- readc('revision_stage0d/source_extract/TRACKB_REVISION_TRANSFER_PACKAGE/11_ANALYSIS_RESULTS_AND_AUDITS/Primary_Results/Phase_1B_20260928_205553_Tables/TRACKB_GENCODE_SENS_composition_MEAN_GRADIENT_STATS.csv')

metric <- function(section, name, value, unit, model, n_subjects, n_obs, family, source, script, figure, status='FROZEN', note='')
  data.frame(section,name,value=as.character(value),unit,model,n_subjects,n_observations=n_obs,
             feature_family=family,source_file=source,reproducible_script=script,figure_or_table=figure,status,note)
master <- do.call(rbind, list(
  metric('cohort','source_records',523,'records','provenance reconstruction',30,523,'all','revision_stage0/STAGE0_AUDIT_REPORT.md','revision_stage0/scripts/stage0_audit.py','Figure 1'),
  metric('cohort','unresolved_duplicate_records_excluded',44,'records','provenance exclusion',30,44,'all','revision_stage0/STAGE0_AUDIT_REPORT.md','revision_stage0/scripts/stage0_audit.py','Figure 1'),
  metric('cohort','primary_subjects',30,'subjects','D01-D07 paired cohort',30,153,'all','revision_stage0c/data/FULL_153_pairs.csv','revision_stage0c/scripts/prepare_qc.ps1','Figure 1'),
  metric('cohort','paired_observations_full',153,'pairs','FULL',30,153,'all','revision_stage0c/data/FULL_153_pairs.csv','revision_stage0c/scripts/prepare_qc.ps1','Figure 1'),
  metric('cohort','biological_samples_full',306,'samples','FULL',30,306,'all','revision_stage0c/data/FULL_153_pairs.csv','revision_stage0c/scripts/prepare_qc.ps1','Figure 1'),
  metric('cohort','paired_observations_qc_pass',141,'pairs','QC PASS-only',30,141,'all','revision_stage0c/data/QC_PASS_ONLY_141_pairs.csv','revision_stage0c/scripts/prepare_qc.ps1','Figure 5'),
  metric('legacy_linear','miRNA_M1_significant',sum(mi$M1_significant,na.rm=TRUE),'features','M1 limma duplicateCorrelation',30,153,'507 miRNA','revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_miRNA.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figure 3/Table 1'),
  metric('legacy_linear','miRNA_M1_direction','189 negative / 3 positive','features','M1 limma duplicateCorrelation',30,153,'507 miRNA','revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_miRNA.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figure 3/Table 1'),
  metric('legacy_linear','GENCODE_M1_significant',sum(ge$M1_significant,na.rm=TRUE),'features','M1 limma duplicateCorrelation',30,153,'1079 GENCODE','revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figure 4/Table 1'),
  metric('random_slope','miRNA_M3_significant',sum(mi$M3_significant,na.rm=TRUE),'features','M3 RAW Model B',30,153,'507 miRNA','revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_miRNA.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figure 3/Table 2'),
  metric('random_slope','miRNA_M4_significant',sum(mi$M4_significant,na.rm=TRUE),'features','M4 CLR Model B',30,153,'507 miRNA','revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_miRNA.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figure 3/Table 2'),
  metric('random_slope','GENCODE_M3_significant',sum(ge$M3_significant,na.rm=TRUE),'features','M3 RAW Model B',30,153,'1079 GENCODE','revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figure 4/Table 2'),
  metric('random_slope','GENCODE_M4_significant',sum(ge$M4_significant,na.rm=TRUE),'features','M4 CLR Model B',30,153,'1079 GENCODE','revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figure 4/Table 2'),
  metric('strict_core','miRNA_strict_core',sum(core$layer=='miRNA'),'features','M1+M3+M4+QC+RS validity',30,153,'507 miRNA','revision_stage1/stage1b_full_refit/STAGE1B_STRICT_CORE_ROBUST.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figures 3/5/Table 2'),
  metric('strict_core','GENCODE_strict_core',sum(core$layer=='GENCODE'),'features','M1+M3+M4+QC+RS validity',30,153,'1079 GENCODE','revision_stage1/stage1b_full_refit/STAGE1B_STRICT_CORE_ROBUST.csv','revision_stage1/stage1b_full_refit/scripts/stage1b_full_refit.R','Figures 4/5/Table 2'),
  metric('mean_gradient','miRNA_mean_significant',sum(mean_mi$BH_FDR<.05,na.rm=TRUE),'features','legacy average gradient',30,153,'507 miRNA',mean_mi_path,'historical Phase_1B_locked_execution_v1.R','Figure 2/Table 1','FROZEN_LEGACY','Not random-slope tested; average estimand has composition sensitivity file.'),
  metric('mean_gradient','GENCODE_mean_significant',sum(mean_ge$BH_FDR<.05,na.rm=TRUE),'features','legacy average gradient',30,153,'1079 GENCODE',mean_ge_path,'historical Phase_1B_locked_execution_v1.R','Figure 4/Table 1','FROZEN_LEGACY','Not random-slope tested; average estimand has composition sensitivity file.'),
  metric('mean_gradient_sensitivity','miRNA_CLR_recentered_mean_significant',sum(mean_mi_comp$BH_FDR<.05,na.rm=TRUE),'features','composition-sensitive recentered average gradient',30,153,'507 miRNA','revision_stage0d/.../TRACKB_MIRNA_SENS_composition_MEAN_GRADIENT_STATS.csv','historical Phase_1B_locked_execution_v1.R','Supplementary Table S2','FROZEN_LEGACY','Same data in a different composition coordinate; not an independent validation.'),
  metric('mean_gradient_sensitivity','GENCODE_CLR_recentered_mean_significant',sum(mean_ge_comp$BH_FDR<.05,na.rm=TRUE),'features','composition-sensitive recentered average gradient',30,153,'1079 GENCODE','revision_stage0d/.../TRACKB_GENCODE_SENS_composition_MEAN_GRADIENT_STATS.csv','historical Phase_1B_locked_execution_v1.R','Supplementary Table S2','FROZEN_LEGACY','Same data in a different composition coordinate; not an independent validation.'),
  metric('categorical_temporal','miRNA_legacy_categorical_significant',139,'features','legacy categorical limma model',30,153,'507 miRNA','revision_stage0c/statistics/miRNA_FULL_complete_results.csv','revision_stage0c/scripts/qc_sensitivity.R','Table 1','FROZEN_LEGACY','Categorical omnibus result; not random-slope validated.'),
  metric('categorical_temporal','GENCODE_legacy_categorical_significant',115,'features','legacy categorical limma model',30,153,'1079 GENCODE','revision_stage0c/statistics/GENCODE_FULL_complete_results.csv','revision_stage0c/scripts/qc_sensitivity.R','Table 1','FROZEN_LEGACY','Categorical omnibus result; not random-slope validated.'),
  metric('categorical_temporal_sensitivity','miRNA_recentered_categorical_significant',60,'features','composition-sensitive recentered categorical model',30,153,'507 miRNA','revision_stage0d/.../TRACKB_MIRNA_SENS_composition_CATEGORICAL_TIME_STATS.csv','historical Phase_1B_locked_execution_v1.R','Supplementary Table S5','FROZEN_LEGACY','QC PASS-only legacy categorical count is 89; not random-slope validated.'),
  metric('categorical_temporal_sensitivity','GENCODE_recentered_categorical_significant',92,'features','composition-sensitive recentered categorical model',30,153,'1079 GENCODE','revision_stage0d/.../TRACKB_GENCODE_SENS_composition_CATEGORICAL_TIME_STATS.csv','historical Phase_1B_locked_execution_v1.R','Supplementary Table S5','FROZEN_LEGACY','QC PASS-only legacy categorical count is 92; not random-slope validated.'),
  metric('cross_layer','RAW_pair_Pearson','0.3692 [0.1065, 0.5697]','r [95% CI]','LOSO score + subject bootstrap',30,153,'shared libraries','revision_stage4/STAGE4_CORRELATION_SUMMARY.csv','revision_stage4/scripts/stage4_crosslayer_analysis.R','Supplementary Figure S7'),
  metric('cross_layer','CLR_pair_Pearson','0.8607 [0.7727, 0.9155]','r [95% CI]','LOSO CLR score + subject bootstrap',30,153,'shared libraries','revision_stage4/STAGE4_CORRELATION_SUMMARY.csv','revision_stage4/scripts/stage4_crosslayer_analysis.R','Supplementary Figure S7','FROZEN','Composition-coordinate dependent; technical consistency only.'),
  metric('cross_layer','RAW_subject_slope_Pearson','0.4770 [-0.4141, 0.8889]','r [95% CI]','subject slope + subject bootstrap',30,30,'shared libraries','revision_stage4/STAGE4_CORRELATION_SUMMARY.csv','revision_stage4/scripts/stage4_crosslayer_analysis.R','Supplementary Figure S7','FROZEN','Influenced by subject12; not robust.'),
  metric('xdec','XDec_formal_analysis','NO-GO','status','technical validation',30,306,'81 informative RNAs','revision_stage2/XDEC_GO_NO_GO.md','revision_stage2/scripts/stage2_analysis.R','Supplementary Table S7','LIMITED','No validated sample-level cargo proportions.')
))
writec(master,'STAGE5_RESULTS_MASTER.csv')
writec(master,'01_results_freeze/STAGE5_RESULTS_MASTER.csv')

# Source tables
table1 <- data.frame(
  Layer=c('miRNA','GENCODE-annotated exRNA'), Retained=c(507,1079),
  Mean_gradient_significant=c(313,711), Mean_CSF_high=c(107,573), Mean_plasma_high=c(206,138),
  Linear_M1_significant=c(192,200), Linear_M1_negative=c(189,175), Linear_M1_positive=c(3,25),
  Categorical_M1_significant=c(139,115),
  Method_note=c('Legacy limma+duplicateCorrelation; average/categorical are not random-slope estimands',
                'Legacy limma+duplicateCorrelation; average/categorical are not random-slope estimands'))
table2 <- data.frame(
  Layer=c('miRNA','GENCODE-annotated exRNA'), M1_linear=c(192,200), M3_RAW_random_slope=c(146,145),
  M4_CLR_random_slope=c(70,143), M1_M3_M4_same_direction_significant=c(34,109), Strict_core=c(16,33),
  Boundary_supported=c(sum(boundary$layer=='miRNA'),sum(boundary$layer=='GENCODE')),
  QC_M3_significant=c(sum(mi$QC_M3_significant,na.rm=T),sum(ge$QC_M3_significant,na.rm=T)),
  QC_M4_significant=c(sum(mi$QC_M4_significant,na.rm=T),sum(ge$QC_M4_significant,na.rm=T)),
  Full_QC_slope_correlation_M3=c(0.9722,0.9746), Full_QC_slope_correlation_M4=c(0.9735,0.9757),
  BH_family=c(507,1079))
writec(table1,'04_tables/Table_1_revised.csv'); writec(table2,'04_tables/Table_2_revised.csv')

labels <- sprintf('D%02d',1:7)
subjects <- sort(unique(pairs$subject))
pair_mat <- outer(subjects,labels,function(s,l) as.integer(paste(s,l,sep='::') %in% pairs$pair_id))
colnames(pair_mat)<-labels; rownames(pair_mat)<-subjects
s1 <- data.frame(subject=subjects,pair_mat,check.names=FALSE)
writec(s1,'04_tables/Supplementary_Table_S1_sampling_architecture.csv')

selected_mi <- c('hsa-miR-375','hsa-miR-323a-3p','hsa-miR-148a-5p','hsa-miR-338-3p','hsa-miR-381-3p','hsa-miR-410-3p')
selected_ge <- c('MT-TL2:Mt_tRNA','RN7SL1:antisense','IL8:protein_coding','HBA1:protein_coding','TMSB4X:protein_coding','EIF4A1:retained_intron')
s2 <- repfeat[repfeat$feature %in% c(selected_mi,selected_ge),]
writec(s2,'04_tables/Supplementary_Table_S2_representative_features.csv')

diag_summary <- rbind(
 data.frame(layer='miRNA',model=c('M3','M4'),attempted=507,evaluable=c(sum(mi$M3_numerical_evaluable),sum(mi$M4_numerical_evaluable)),RS_valid=c(sum(mi$class_M3_random_slope_validated),sum(mi$class_M4_random_slope_validated)),RS_boundary=c(sum(mi$M3_model_status=='RS_BOUNDARY'),sum(mi$M4_model_status=='RS_BOUNDARY')),unstable=c(sum(mi$M3_model_status=='NUMERICAL_UNSTABLE'),sum(mi$M4_model_status=='NUMERICAL_UNSTABLE'))),
 data.frame(layer='GENCODE',model=c('M3','M4'),attempted=1079,evaluable=c(sum(ge$M3_numerical_evaluable),sum(ge$M4_numerical_evaluable)),RS_valid=c(sum(ge$class_M3_random_slope_validated),sum(ge$class_M4_random_slope_validated)),RS_boundary=c(sum(ge$M3_model_status=='RS_BOUNDARY'),sum(ge$M4_model_status=='RS_BOUNDARY')),unstable=c(sum(ge$M3_model_status=='NUMERICAL_UNSTABLE'),sum(ge$M4_model_status=='NUMERICAL_UNSTABLE'))))
writec(diag_summary,'04_tables/Supplementary_Table_S3_model_diagnostics.csv')
writec(corr,'04_tables/Supplementary_Table_S4_cross_layer_bootstrap.csv')
s5 <- data.frame(Item=c('FULL cohort','QC PASS-only','RAW gradient','CLR gradient','M3/M4 random effect','BH families','Collection labels','TMM-voom','XDec'),
 Value=c('30 subjects; 153 complete CSF-plasma pairs','30 subjects; 141 pairs; both samples PASS','log2(CSF RPM+1)-log2(plasma RPM+1)','fixed 507/1079 feature reference; pseudocount 1','(1 + centered_label || subject)','507 miRNA or 1079 GENCODE per model','D01-D07 ordered labels centered at D04; not validated exact days','Not performed; true counts unavailable','NO-GO; no validated sample-level cargo proportions'))
writec(s5,'04_tables/Supplementary_Table_S5_analysis_specifications.csv')
core_full <- rbind(mi[mi$strict_core_robust,],ge[ge$strict_core_robust,])
writec(core_full,'04_tables/Supplementary_Table_S6_strict_core_results.csv')
xdec_tab <- data.frame(Item=c('Reference profiles','Informative RNA mapping','Background transform matrix','Original sample proportions','Formal XDec','Permitted interpretation'),
 Status=c('81 x 75 recovered','81/81 mapped','Missing version-matched Atlas_Update.txt','Not recovered','NO-GO','High-abundance/biotype/center-shift analyses only'),
 Evidence=c('revision_stage2/XDEC_TECHNICAL_VALIDATION.md','revision_stage2/XDEC_TECHNICAL_VALIDATION.md','revision_stage2/XDEC_TECHNICAL_VALIDATION.md','revision_stage2/XDEC_GO_NO_GO.md','revision_stage2/XDEC_GO_NO_GO.md','revision_stage2/STAGE2_FINAL_REPORT.md'))
writec(xdec_tab,'04_tables/Supplementary_Table_S7_XDec_feasibility.csv')
writec(mi,'04_tables/Supplementary_Data_1_miRNA_revised.csv')
writec(ge,'04_tables/Supplementary_Data_2_GENCODE_revised.csv')
file.copy(file.path(root,'revision_stage0d/source_extract/TRACKB_REVISION_TRANSFER_PACKAGE/03_SUPPLEMENTARY_MATERIALS/Supplementary_Data/Supplementary_Data_3_provenance_crosswalk.csv'),
          file.path(out,'04_tables/Supplementary_Data_3_provenance_crosswalk_UNCHANGED.csv'),overwrite=TRUE)

# Figure source data
writec(core_full,'02_main_figures/Figure3_4_strict_core_source_data.csv')
writec(decomp[decomp$dataset=='FULL' & decomp$scale=='RAW',],'02_main_figures/Figure3_4_decomposition_source_data.csv')
writec(mean_mi,'02_main_figures/Figure2_source_data.csv')
writec(corr,'03_supplementary_figures/Supplementary_Figure_S7_source_data.csv')

# Unified visual system
COL <- list(CSF='#0072B2', plasma='#D55E00', RAW='#4D4D4D', CLR='#009E73', neg='#3B6FB6', pos='#C94C4C', gray='#B8B8B8', dark='#263238', pale='#EEF3F6')
theme_base <- function(mar=c(4,4,2.2,1),cex=0.78){par(mar=mar,mgp=c(2.3,.7,0),tcl=-.25,las=1,cex=cex,family='sans',bty='l')}
mark_panel <- function(lab, adj=-0.12){mtext(lab,side=3,line=.35,adj=adj,font=2,cex=1.35,xpd=NA)}
save_multi <- function(name, folder, w, h, fun){
  base <- file.path(out,folder,name)
  svg(paste0(base,'.svg'),width=w,height=h,pointsize=9); fun(); dev.off()
  pdf(paste0(base,'.pdf'),width=w,height=h,pointsize=9,useDingbats=FALSE); fun(); dev.off()
  png(paste0(base,'.png'),width=w,height=h,units='in',res=300,pointsize=9,type='cairo'); fun(); dev.off()
  tiff(paste0(base,'.tiff'),width=w,height=h,units='in',res=300,compression='lzw',pointsize=9,type='cairo'); fun(); dev.off()
}

bar_counts <- function(vals,nms,cols,ylim=NULL,main=''){if(is.null(ylim))ylim=c(0,max(vals)*1.18); bp<-barplot(vals,names.arg=nms,col=cols,border=NA,ylim=ylim,main=main,ylab='Significant features'); text(bp,vals,labels=vals,pos=3,cex=.8)}
forest_two <- function(d, feature_col='feature', m3='M3_slope',m3lo='M3_CI_lower',m3hi='M3_CI_upper',m4='M4_slope',m4lo='M4_CI_lower',m4hi='M4_CI_upper',main=''){
  o<-order(d[[m4]]); d<-d[o,]; y<-seq_len(nrow(d)); rng<-range(c(d[[m3lo]],d[[m3hi]],d[[m4lo]],d[[m4hi]]),na.rm=T)
  plot(NA,xlim=rng,ylim=c(.5,nrow(d)+.5),yaxt='n',xlab='Slope per collection label',ylab='',main=main); abline(v=0,col='#999999',lty=2)
  segments(d[[m3lo]],y-.12,d[[m3hi]],y-.12,col=COL$RAW,lwd=1.5); points(d[[m3]],y-.12,pch=16,col=COL$RAW,cex=.65)
  segments(d[[m4lo]],y+.12,d[[m4hi]],y+.12,col=COL$CLR,lwd=1.5); points(d[[m4]],y+.12,pch=17,col=COL$CLR,cex=.65)
  axis(2,at=y,labels=sub('^hsa-','',d[[feature_col]]),las=1,cex.axis=.56); legend('bottomright',c('M3 RAW','M4 CLR'),col=c(COL$RAW,COL$CLR),pch=c(16,17),bty='n',cex=.65)
}
profile_panel <- function(prof, feature, biofluid, scale='RAW', main='', col=COL$CSF){
  d<-prof[prof$feature==feature & prof$dataset=='FULL' & prof$scale==scale & prof$biofluid==biofluid,]; d<-d[order(d$label_numeric),]
  if(!nrow(d) && scale=='RAW'){
    obj<-if(feature %in% rownames(mi_pair$CSF_RPM)) mi_pair else ge_pair
    mat<-if(biofluid=='CSF') obj$CSF_RPM else obj$plasma_RPM
    if(feature %in% rownames(mat)){
      vals<-log2(as.numeric(mat[feature,])+1); labs<-obj$pairs$collection_label
      d<-data.frame(label_numeric=1:7,mean=sapply(labels,function(z)mean(vals[labs==z],na.rm=T)),SE.x=sapply(labels,function(z)sd(vals[labs==z],na.rm=T)/sqrt(sum(labs==z))))
    }
  }
  if(!nrow(d)){plot.new();title(main);text(.5,.5,'Profile unavailable');return()}
  yl<-range(c(d$mean-d$SE.x,d$mean+d$SE.x),na.rm=T); plot(d$label_numeric,d$mean,type='b',pch=16,col=col,lwd=1.7,xaxt='n',xlab='Collection label',ylab=if(scale=='RAW') 'Mean log2(RPM+1)' else 'Mean CLR',main=main,ylim=yl)
  arrows(d$label_numeric,d$mean-d$SE.x,d$label_numeric,d$mean+d$SE.x,angle=90,code=3,length=.035,col=col); axis(1,1:7,labels)
}
gradient_panel <- function(matfile, feature, main=''){
  d<-readc(matfile); row<-d[d$feature==feature,-1,drop=FALSE]; vals<-as.numeric(row[1,]); ids<-names(row); labs<-sub('.*::','',ids); m<-sapply(labels,function(z)mean(vals[labs==z],na.rm=T)); se<-sapply(labels,function(z)sd(vals[labs==z],na.rm=T)/sqrt(sum(labs==z)))
  yl<-range(m-se,m+se,na.rm=T); plot(1:7,m,type='b',pch=16,col=COL$dark,lwd=1.7,xaxt='n',xlab='Collection label',ylab='Mean CSF-plasma gradient',main=main,ylim=yl); arrows(1:7,m-se,1:7,m+se,angle=90,code=3,length=.035,col=COL$dark); abline(h=0,lty=2,col='#999999');axis(1,1:7,labels)
}

# Figure 1
fig1 <- function(){
  layout(matrix(c(1,2,3),1,3),widths=c(1.05,1.55,1.2));
  theme_base(c(1,1,2.4,1)); plot.new(); plot.window(c(0,1),c(0,1)); title('Provenance reconstruction',font.main=2); mark_panel('A',-.05)
  ys<-c(.86,.65,.44,.23); txt<-c('523 GEO/library records','44 unresolved duplicates excluded','30 subjects / 153 paired labels','306 biological samples');
  for(i in 1:4){rect(.08,ys[i]-.065,.92,ys[i]+.065,col=if(i==2)'#FDE2E2' else COL$pale,border=COL$dark);text(.5,ys[i],txt[i],cex=.82);if(i<4)arrows(.5,ys[i]-.075,.5,ys[i+1]+.075,length=.08)}
  theme_base(c(5.3,5,2.4,1)); image(1:7,seq_along(subjects),t(pair_mat[nrow(pair_mat):1,]),col=c('white',COL$CSF),axes=FALSE,xlab='',ylab='Subjects',main='Paired sampling availability'); mark_panel('B'); axis(1,1:7,paste0(labels,'\nn=',colSums(pair_mat)),line=.25,cex.axis=.57); mtext('Collection label',side=1,line=3.6,cex=.8); axis(2,at=seq_along(subjects),labels=rev(sub('EXR-KJENS1Sample','S',sub('-DO','',subjects))),las=1,cex.axis=.45); box()
  theme_base(c(1,1,2.4,1)); plot.new(); plot.window(c(0,1),c(0,1)); title('Frozen analysis workflow',font.main=2); mark_panel('C',-.05); nodes<-list(c(.5,.86,'Paired CSF-plasma\nRPM'),c(.25,.64,'507 miRNA'),c(.75,.64,'1,079 GENCODE'),c(.25,.40,'RAW + CLR\ngradients'),c(.75,.40,'Random-slope\n+ QC'),c(.5,.16,'16 miRNA +\n33 GENCODE\nstrict core'));
  for(n in nodes) {nx<-as.numeric(n[1]); ny<-as.numeric(n[2]); hh<-if(ny==.16).095 else .078; rect(nx-.20,ny-hh,nx+.20,ny+hh,col=COL$pale,border=COL$dark);text(nx,ny,n[3],cex=.62)}; arrows(.5,.78,.25,.72,length=.07);arrows(.5,.78,.75,.72,length=.07);arrows(.25,.56,.25,.48,length=.07);arrows(.75,.56,.75,.48,length=.07);arrows(.25,.32,.43,.25,length=.07);arrows(.75,.32,.57,.25,length=.07); text(.5,.025,'D01-D07: ordered collection labels',cex=.55,col='#555555')
}
save_multi('Figure_1_revised','02_main_figures',7.0,4.3,fig1)

# Figure 2
stable <- c('hsa-miR-204-5p','hsa-miR-1298-5p','hsa-miR-6852-5p','hsa-miR-122-5p')
fig2 <- function(){layout(matrix(c(1,2,3,4),2,2,byrow=T)); theme_base(); hist(mean_mi$effect,breaks=35,col='#A8C7DC',border='white',main='Average relative gradients',xlab='Equal-label mean CSF-plasma gradient');mark_panel('A');abline(v=0,lty=2)
 theme_base(); cls<-c('CSF-high'=sum(mean_mi$BH_FDR<.05 & mean_mi$effect>0),'Plasma-high'=sum(mean_mi$BH_FDR<.05 & mean_mi$effect<0),'Not FDR-significant'=sum(mean_mi$BH_FDR>=.05)); bar_counts(cls,names(cls),c(COL$CSF,COL$plasma,COL$gray),main='Legacy average-gradient classes');mark_panel('B')
 theme_base(); q<-pmax(mean_mi$BH_FDR,1e-300); col<-ifelse(q<.05,ifelse(mean_mi$effect>0,COL$CSF,COL$plasma),COL$gray); plot(mean_mi$effect,-log10(q),pch=16,cex=.55,col=adjustcolor(col,.7),xlab='Mean gradient',ylab='-log10(BH FDR)',main='Effect and FDR');mark_panel('C');abline(v=0,lty=2);abline(h=-log10(.05),lty=3)
 theme_base(c(4,7,2.2,1)); d<-mean_mi[match(stable,mean_mi$feature_id),]; y<-4:1; plot(d$effect,y,xlim=range(d$CI_lower,d$CI_upper),ylim=c(.5,4.5),yaxt='n',pch=18,cex=1.1,col=ifelse(d$effect>0,COL$CSF,COL$plasma),xlab='Mean gradient (95% CI)',ylab='',main='Stable-gradient representatives');mark_panel('D');segments(d$CI_lower,y,d$CI_upper,y,col=ifelse(d$effect>0,COL$CSF,COL$plasma),lwd=2);axis(2,y,sub('hsa-','',stable),las=1,cex.axis=.7);abline(v=0,lty=2)
}
save_multi('Figure_2_revised','02_main_figures',7.0,6.3,fig2)

# Figure 3
core_mi <- mi[mi$strict_core_robust,]
show_mi <- c('hsa-miR-375','hsa-miR-323a-3p','hsa-miR-148a-5p')
fig3 <- function(){
 layout(matrix(c(1,2,2,3,4,5,6,7,8,9,10,11),4,3,byrow=T),heights=c(1.15,1,1,1));
 theme_base(); bar_counts(c(192,146,70,16),c('M1\nRAW','M3\nRAW+RS','M4\nCLR+RS','Strict\ncore'),c('#8C8C8C',COL$RAW,COL$CLR,'#6A3D9A'),main='miRNA evidence contraction');mark_panel('A')
 theme_base(c(4,7,2.2,1),.7); forest_two(core_mi,main='16 strict-core miRNA: M3 vs M4');mark_panel('B')
 for(i in seq_along(show_mi)){f<-show_mi[i];theme_base();profile_panel(mirprof,f,'CSF','RAW',paste0(sub('hsa-','',f),' | CSF'),COL$CSF);mark_panel(c('C','D','E')[i]);theme_base();profile_panel(mirprof,f,'plasma','RAW',paste0(sub('hsa-','',f),' | plasma'),COL$plasma);theme_base();gradient_panel('revision_stage0c/data/miRNA_FULL_gradients.csv',f,paste0(sub('hsa-','',f),' | paired gradient'))}
}
save_multi('Figure_3_revised','02_main_figures',7.0,8.8,fig3)

# Figure 4
core_ge <- ge[ge$strict_core_robust,]
show_ge <- c('MT-TL2:Mt_tRNA','RN7SL1:antisense','IL8:protein_coding','HBA1:protein_coding')
fig4 <- function(){layout(matrix(c(1,2,2,3,4,5,6,7,8,9,10,11,12,13,14),5,3,byrow=T),heights=c(1,1,1,1,1));
 theme_base(); bar_counts(c(200,145,143,33),c('M1\nRAW','M3\nRAW+RS','M4\nCLR+RS','Strict\ncore'),c('#8C8C8C',COL$RAW,COL$CLR,'#6A3D9A'),main='GENCODE evidence contraction');mark_panel('A')
 theme_base(); plot(core_ge$M3_slope,core_ge$M4_slope,pch=16,col=adjustcolor('#6A3D9A',.55),xlab='M3 RAW slope',ylab='M4 CLR slope',main='All 33 strict-core transcript-derived features');mark_panel('B');abline(0,1,lty=2,col='#777777'); d<-core_ge[match(show_ge,core_ge$feature),]; hcols<-c('#0072B2','#D55E00','#009E73','#CC79A7'); points(d$M3_slope,d$M4_slope,pch=21,bg=hcols,col='white',cex=1.15);legend('topleft',sub(':.*','',show_ge),pch=21,pt.bg=hcols,pt.cex=.85,bty='n',cex=.52,ncol=2)
 for(i in seq_along(show_ge)){f<-show_ge[i];theme_base();profile_panel(geprof,f,'CSF','RAW',paste0(sub(':.*','',f),' | CSF'),COL$CSF);mark_panel(c('C','D','E','F')[i]);theme_base();profile_panel(geprof,f,'plasma','RAW',paste0(sub(':.*','',f),' | plasma'),COL$plasma);theme_base();gradient_panel('revision_stage0c/data/GENCODE_FULL_gradients.csv',f,paste0(sub(':.*','',f),' | paired gradient'))}
}
save_multi('Figure_4_revised','02_main_figures',7.0,8.8,fig4)

# Figure 5 robustness synthesis
fig5 <- function(){layout(matrix(1:4,2,2,byrow=T));
 theme_base(); m<-rbind(c(192,146,70,16),c(200,145,143,33)); bp<-barplot(t(m),beside=T,names.arg=c('miRNA','GENCODE'),col=c('#8C8C8C',COL$RAW,COL$CLR,'#6A3D9A'),ylim=c(0,230),main='Evidence retained across frameworks',ylab='Features');mark_panel('A');legend('topright',c('M1','M3','M4','Strict core'),fill=c('#8C8C8C',COL$RAW,COL$CLR,'#6A3D9A'),bty='n',cex=.7)
 theme_base(); cats<-c('miRNA strict'=16,'miRNA boundary'=sum(boundary$layer=='miRNA'),'GENCODE strict'=33,'GENCODE boundary'=sum(boundary$layer=='GENCODE')); bar_counts(cats,names(cats),c('#6A3D9A','#C8B3D8','#6A3D9A','#C8B3D8'),main='Strict vs boundary-supported');mark_panel('B')
 theme_base(); plot(mi$M3_slope,mi$QC_M3_slope,pch=16,cex=.45,col=adjustcolor(COL$RAW,.45),xlab='FULL M3 slope',ylab='QC PASS-only M3 slope',main='QC sensitivity');mark_panel('C');abline(0,1,col=COL$CSF,lwd=2);legend('topleft',sprintf('miRNA r = %.3f',cor(mi$M3_slope,mi$QC_M3_slope,use='complete')),bty='n',cex=.75)
 theme_base(); rates<-c('mi M3'=mean(mi$M3_model_status=='RS_BOUNDARY'),'mi M4'=mean(mi$M4_model_status=='RS_BOUNDARY'),'GEN M3'=mean(ge$M3_model_status=='RS_BOUNDARY'),'GEN M4'=mean(ge$M4_model_status=='RS_BOUNDARY')); vv<-round(100*rates,1);bp<-barplot(vv,names.arg=names(rates),col=c(COL$RAW,COL$CLR,COL$RAW,COL$CLR),border=NA,ylim=c(0,50),main='Random-slope boundary frequency',ylab='Percent of full feature family');mark_panel('D');text(bp,vv,vv,pos=3,cex=.8)
}
save_multi('Figure_5_revised','02_main_figures',7.0,6.2,fig5)

# Supplementary figures
sup1 <- function(){layout(matrix(1:4,2,2,byrow=T));theme_base();v<-c(153,141);bp<-barplot(v,names.arg=c('FULL','QC PASS'),col=c(COL$gray,COL$CSF),border=NA,ylim=c(0,170),main='Complete pairs',ylab='Complete CSF-plasma pairs');mark_panel('A');text(bp,v,v,pos=3);theme_base();bar_counts(c(146,119,70,46),c('M3 F','M3 QC','M4 F','M4 QC'),c(COL$RAW,COL$RAW,COL$CLR,COL$CLR),main='miRNA significant');mark_panel('B');theme_base();bar_counts(c(145,59,143,86),c('M3 F','M3 QC','M4 F','M4 QC'),c(COL$RAW,COL$RAW,COL$CLR,COL$CLR),main='GENCODE significant');mark_panel('C');theme_base();plot(mi$M4_slope,mi$QC_M4_slope,pch=16,cex=.45,col=adjustcolor(COL$CLR,.45),xlab='FULL M4 slope',ylab='QC M4 slope',main='miRNA M4 effect stability');mark_panel('D');abline(0,1,col=COL$dark,lwd=2)}
save_multi('Supplementary_Figure_S1_QC_sensitivity','03_supplementary_figures',7,6.2,sup1)

sup2 <- function(){layout(matrix(1:4,2,2,byrow=T)); for(k in 1:2){d<-if(k==1)mi else ge; theme_base(); plot(d$M2_slope,d$M2_reference_sensitivity_pc0_5_slope,pch=16,cex=.4,col=adjustcolor(COL$CLR,.4),xlab='CLR slope, pseudocount 1',ylab='CLR slope, pseudocount 0.5',main=c('miRNA','GENCODE')[k]);mark_panel(c('A','B')[k]);abline(0,1)};theme_base();plot(mi$M2_slope,mi$M2_reference_sensitivity_excluded_ref_slope,pch=16,cex=.4,col=adjustcolor(COL$CSF,.4),xlab='Full reference',ylab='Exclude miR-486-5p',main='miRNA reference sensitivity');mark_panel('C');abline(0,1);theme_base();plot(ge$M2_slope,ge$M2_reference_sensitivity_excluded_ref_slope,pch=16,cex=.4,col=adjustcolor(COL$plasma,.4),xlab='Full reference',ylab='Exclude RNY4',main='GENCODE reference sensitivity');mark_panel('D');abline(0,1)}
save_multi('Supplementary_Figure_S2_CLR_sensitivity','03_supplementary_figures',7,6.2,sup2)

sup3 <- function(){layout(matrix(1:4,2,2,byrow=T)); theme_base();barplot(t(as.matrix(diag_summary[,c('RS_valid','RS_boundary','unstable')])),beside=T,names.arg=paste(diag_summary$layer,diag_summary$model),col=c(COL$CLR,'#E69F00','#CC79A7'),main='Model status',ylab='Features');mark_panel('A');legend('topright',c('RS valid','RS boundary','Unstable'),fill=c(COL$CLR,'#E69F00','#CC79A7'),bty='n',cex=.65); for(i in 1:2){d<-if(i==1)mi else ge;theme_base();plot(d$M3_slope,d$M3_CR2_slope,pch=16,cex=.4,col=adjustcolor(COL$RAW,.4),xlab='Model B M3 slope',ylab='CR2 slope',main=c('miRNA','GENCODE')[i]);mark_panel(c('B','C')[i]);abline(0,1)};theme_base();bar_counts(c(146,133,70,52,145,122,143,137),c('mi M3 B','mi M3 C','mi M4 B','mi M4 C','ge M3 B','ge M3 C','ge M4 B','ge M4 C'),rep(c(COL$RAW,COL$CLR),4),main='Model B vs CR2 FDR hits');mark_panel('D')}
save_multi('Supplementary_Figure_S3_random_slope_CR2','03_supplementary_figures',7,6.2,sup3)

sup4 <- function(){layout(matrix(c(1,2,3,3),2,2,byrow=T),heights=c(1,1.05)); theme_base(); sm<-shift[shift$record_type=='MODEL' & shift$dataset=='FULL' & shift$layer=='miRNA' & shift$outcome=='center_shift',]; if(nrow(sm)){plot(1,sm$slope[1],ylim=range(c(sm$CI_lower[1],sm$CI_upper[1],0)),xlim=c(.5,1.5),xaxt='n',pch=18,cex=1.2,col=COL$CLR,xlab='',ylab='Slope',main='miRNA global center shift');mark_panel('A');arrows(1,sm$CI_lower[1],1,sm$CI_upper[1],angle=90,code=3,length=.06);axis(1,1,'FULL');abline(h=0,lty=2)} else {plot.new();title('miRNA global center shift');mark_panel('A')}; theme_base(); hh<-unique(high[high$dataset=='FULL' & high$biofluid=='CSF' & high$feature %in% c('hsa-miR-486-5p','hsa-miR-451a','RNY4:misc_RNA'),c('feature','slope','CI_lower','CI_upper')]); y<-seq_len(nrow(hh));plot(hh$slope,y,yaxt='n',pch=16,xlim=range(hh$CI_lower,hh$CI_upper),xlab='CSF slope',ylab='',main='High-abundance RNA');mark_panel('B');segments(hh$CI_lower,y,hh$CI_upper,y);axis(2,y,sub('hsa-','',hh$feature),las=1,cex.axis=.65);abline(v=0,lty=2); theme_base(c(5.8,4.4,2.2,1)); bt<-unique(biotype[biotype$dataset=='FULL' & biotype$biofluid=='CSF',c('biotype','CLR_slope')]); if(nrow(bt)){bp<-barplot(bt$CLR_slope,names.arg=bt$biotype,col=ifelse(bt$CLR_slope>0,COL$CSF,COL$plasma),las=2,cex.names=.72,main='CSF biotype CLR slopes',ylab='Slope');mark_panel('C',adj=-.03);abline(h=0)} else {plot.new();title('CSF biotype CLR slopes');mark_panel('C',adj=-.03)}}
save_multi('Supplementary_Figure_S4_composition_analysis_recomposed','03_supplementary_figures',7,6.0,sup4)

old_mi <- c('hsa-miR-338-3p','hsa-miR-381-3p','hsa-miR-410-3p')
sup5 <- function(){layout(matrix(1:9,3,3,byrow=T));for(i in seq_along(old_mi)){f<-old_mi[i];theme_base();profile_panel(mirprof,f,'CSF','RAW',paste0(sub('hsa-','',f),' | CSF'),COL$CSF);mark_panel(c('A','B','C')[i]);theme_base();profile_panel(mirprof,f,'plasma','RAW',paste0(sub('hsa-','',f),' | plasma'),COL$plasma);theme_base();gradient_panel('revision_stage0c/data/miRNA_FULL_gradients.csv',f,'Paired gradient')}}
save_multi('Supplementary_Figure_S5_original_miRNA_representatives','03_supplementary_figures',7,7.2,sup5)

sup6 <- function(){theme_base(c(5.2,7.5,2.5,1.3));ep<-endpoint[endpoint$layer=='miRNA' & endpoint$scale=='RAW' & endpoint$outcome=='gradient' & endpoint$record_type=='SUMMARY' & endpoint$feature %in% core_mi$feature,]; if(nrow(ep)){o<-order(ep$mean_change);y<-seq_len(nrow(ep));plot(ep$mean_change[o],y,yaxt='n',pch=16,cex=1.0,col=COL$RAW,xlim=range(c(ep$CI_lower,ep$CI_upper,0)),xlab='Mean D07-D01 paired-gradient change (95% CI)',ylab='',main='Strict-core miRNA endpoint sensitivity (12 subjects)');segments(ep$CI_lower[o],y,ep$CI_upper[o],y,lwd=1.4);axis(2,y,sub('hsa-','',ep$feature[o]),las=1,cex.axis=.72);abline(v=0,lty=2,col='#777777')}else{plot.new();title('Endpoint sensitivity')}}
save_multi('Supplementary_Figure_S6_endpoint_sensitivity_recomposed','03_supplementary_figures',7,5.2,sup6)

sup7 <- function(){layout(matrix(1:4,2,2,byrow=T));for(sc in c('RAW','CLR')){d<-scores[scores$dataset=='FULL' & scores$scale==sc & scores$scoring_scheme=='ORIGINAL_FULL_FEATURE',];theme_base();lab<-if(sc=='RAW')'A' else 'B';plot(d$miRNA_score,d$GENCODE_score,pch=16,cex=.55,col=adjustcolor(if(sc=='RAW')COL$RAW else COL$CLR,.5),xlab='miRNA score',ylab='GENCODE score',main=paste(sc,'paired-label scores'));mark_panel(lab);abline(lm(GENCODE_score~miRNA_score,d),col=if(sc=='RAW')COL$RAW else COL$CLR,lwd=2)};theme_base();b<-boot[boot$branch=='PRIMARY_FULL_RAW' & boot$metric=='slope_Pearson',];hist(b$estimate,breaks=35,col='#BDBDBD',border='white',main='RAW slope bootstrap',xlab='Patient-level Pearson r');mark_panel('C');abline(v=0,lty=2);abline(v=.477,col=COL$plasma,lwd=2);theme_base();lr<-loo[loo$scale=='RAW' & loo$analysis_type=='FIXED_SCORE_OMISSION' & loo$metric=='slope_Pearson',];o<-order(abs(lr$change),decreasing=T);lr<-lr[o[1:min(12,nrow(lr))],];barplot(lr$change,names.arg=sub('EXR-KJENS1Sample','S',sub('-DO','',lr$omitted_subject)),las=2,col=ifelse(lr$direction_flip,COL$plasma,COL$gray),main='Patient influence on RAW slope r',ylab='Change after omission');mark_panel('D')}
save_multi('Supplementary_Figure_S7_cross_layer_technical_consistency','03_supplementary_figures',7,6.2,sup7)

# Graphical abstract
ga <- function(){theme_base(c(.2,.2,.2,.2));plot.new();plot.window(c(0,1),c(0,1));rect(.01,.01,.99,.99,col='white',border=NA);text(.04,.93,'Longitudinal relative exRNA gradients after aSAH',adj=0,font=2,cex=1.25,col=COL$dark)
 rect(.05,.58,.27,.82,col='#DCEEF8',border=COL$CSF,lwd=2);text(.16,.76,'30 patients',font=2);text(.16,.68,'153 paired labels');text(.16,.61,'CSF + plasma');
 arrows(.29,.70,.39,.70,length=.08,lwd=2,col=COL$dark);rect(.40,.55,.62,.84,col='#F3F3F3',border=COL$dark,lwd=2);text(.51,.78,'Relative gradients',font=2);text(.51,.69,'RAW + CLR');text(.51,.62,'Random slopes + QC');
 arrows(.64,.70,.73,.70,length=.08,lwd=2,col=COL$dark);rect(.74,.55,.96,.84,col='#EEE4F3',border='#6A3D9A',lwd=2);text(.85,.77,'Robust features',font=2);text(.85,.68,'16 miRNA');text(.85,.61,'33 transcript-derived');
 rect(.08,.16,.44,.43,col='#E7F3FA',border=NA);text(.26,.37,'Distinct biofluid contributions',font=2);segments(.16,.23,.16,.33,col=COL$CSF,lwd=5);segments(.25,.23,.25,.30,col=COL$plasma,lwd=5);text(.16,.19,'CSF');text(.25,.19,'plasma');text(.36,.28,'paired gradient',cex=.85);arrows(.31,.25,.40,.25,length=.07)
 rect(.55,.16,.92,.43,col='#F7EFE7',border=NA);text(.735,.37,'Interpretation boundary',font=2);text(.735,.29,'Relative abundance, not concentration');text(.735,.23,'Composition-aware; not carrier proof');text(.735,.17,'Ordered labels, not exact post-bleed days',cex=.78)
}
save_multi('Graphical_Abstract_revised','05_graphical_abstract',7.0,4.3,ga)

# Reports and legends
legends <- c(
'# Stage 5 figure legends',
'',
'**Figure 1. Study design, sample provenance, and paired longitudinal sampling.** (A) Reconstruction from 523 GEO/library records; 44 unresolved duplicate records were excluded. (B) Availability of complete CSF-plasma pairs across D01-D07. Numbers below labels are contributing subjects. (C) Frozen analysis flow. Collection labels are ordered identifiers and are not interpreted as verified exact post-bleed days.',
'',
'**Figure 2. Persistent relative CSF-plasma miRNA gradients.** Legacy equal-label mean-gradient estimates are shown for 507 retained miRNAs. These estimates reproduce the original limma/duplicateCorrelation analysis and have composition-sensitive companion analyses, but they are not random-slope temporal estimands. Relative abundance is expressed from RPM-derived data.',
'',
'**Figure 3. Robust longitudinal miRNA gradients and biofluid-specific trajectories.** (A) Contraction of the significant set across M1, M3, M4 and strict-core criteria. (B) M3 RAW and M4 CLR slopes with 95% confidence intervals for all 16 strict-core miRNAs. (C-E) Relative CSF abundance, plasma abundance and paired gradients for miR-375, miR-323a-3p and miR-148a-5p, selected post hoc from the frozen strict-core set to represent bilateral, CSF-dominant and plasma-dominant patterns. Points are collection-label means and error bars are standard errors.',
'',
'**Figure 4. GENCODE-annotated extracellular RNA gradients and robust longitudinal features.** The figure shows model-set contraction and relative CSF, plasma and paired-gradient trajectories for MT-TL2, RN7SL1, IL8 and HBA1. These are transcript-derived features detected by small-RNA sequencing and are not measurements of a complete cellular transcriptome.',
'',
'**Figure 5. Robustness of longitudinal exRNA gradients across analytical frameworks.** Counts are not additive across models. Strict-core and random-slope-boundary-supported features are separated. FULL versus QC PASS-only effects and random-slope boundary frequencies summarize sensitivity without treating FDR changes as changes in biological effect.',
'',
'**Supplementary Figure S7. Cross-layer technical consistency and sensitivity analyses.** RAW and CLR cross-fitted scores are derived from the same source libraries. Bootstrap resampling uses patients as clusters. Patient-specific RAW temporal concordance is imprecise and sensitive to subject12; CLR correlations are composition-coordinate dependent. These analyses are technical consistency checks, not independent biological validation.',
'',
'**Supplementary Figure S1. Source-QC sensitivity.** Complete paired observations, FDR-significant M3/M4 counts and FULL versus QC PASS-only slopes are shown. Both members of an affected pair were excluded only in this sensitivity branch.',
'',
'**Supplementary Figure S2. CLR reference and pseudocount sensitivity.** Pseudocount and high-abundance reference exclusions are alternative coordinates of the same data, not independent experiments.',
'',
'**Supplementary Figure S3. Random-slope boundary and CR2 sensitivity.** Model status, Model B versus CR2 slopes and FDR-significant counts are shown. CR2 support does not reclassify an RS-boundary feature.',
'',
'**Supplementary Figure S4. Composition-reference sensitivity analyses.** (A) The time slope of the global miRNA CLR centering term. (B) CSF time slopes and 95% confidence intervals for selected high-abundance RNA features. (C) CSF CLR time slopes across seven recoverable RNA biotypes. These relative-composition analyses do not identify physical carriers. Formal XDec analysis was not performed because the recovered 81 x 75 profiles were reference profiles rather than sample proportions and the version-matched Atlas background transformation was unavailable.',
'',
'**Supplementary Figure S5. Original non-strict-core miRNA representatives.** Separate CSF, plasma and paired-gradient profiles are shown for miR-338-3p, miR-381-3p and miR-410-3p as continuity displays.',
'',
'**Supplementary Figure S6. D01/D07 endpoint-complete sensitivity.** Mean D07-D01 paired-gradient changes and 95% confidence intervals are shown for the 16 strict-core miRNAs among the 12 subjects with both endpoint pairs. All 16 RAW and CLR strict-core gradient directions were concordant with the primary longitudinal analysis. Four subjects had all seven collection labels. The reduced endpoint set has lower precision and does not establish the missing-data mechanism.',
'',
'**Graphical Abstract.** Paired sampling, composition/random-slope/QC sensitivity, 16/33 robust features and distinct biofluid contributions are summarized without implying carrier proof.'
)
writeLines(legends,file.path(out,'STAGE5_FIGURE_LEGENDS.md'),useBytes=TRUE)

session <- capture.output(sessionInfo()); writeLines(session,file.path(out,'logs','STAGE5_sessionInfo.txt'),useBytes=TRUE)
writeLines(c('Stage 5 build completed',paste('Timestamp:',format(Sys.time(),tz='Asia/Shanghai',usetz=TRUE)),paste('Root:',root)),file.path(out,'logs','stage5_build.log'),useBytes=TRUE)

