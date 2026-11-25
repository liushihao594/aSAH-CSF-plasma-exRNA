options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(lme4);library(lmerTest);library(clubSandwich);library(ggplot2)})
set.seed(20261121L)
out<-"revision_stage3";dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(out,"scripts"),showWarnings=FALSE);dir.create(file.path(out,"logs"),showWarnings=FALSE);dir.create(file.path(out,"plots"),showWarnings=FALSE);dir.create(file.path(out,"intermediate"),showWarnings=FALSE)
tol<-1e-4; optimizers<-c("nloptwrap","bobyqa","nlminbwrap")
obj<-readRDS("revision_stage1/intermediate/prepared_stage1_inputs.rds")
strict<-read.csv("revision_stage1/stage1b_full_refit/STAGE1B_STRICT_CORE_ROBUST.csv",check.names=FALSE)
stage1_mi<-read.csv("revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_miRNA.csv",check.names=FALSE)
stage1_ge<-read.csv("revision_stage1/stage1b_full_refit/STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv",check.names=FALSE)
stage1_all<-rbind(stage1_mi,stage1_ge)

rep_mi<-c("hsa-miR-338-3p","hsa-miR-381-3p","hsa-miR-410-3p","hsa-miR-375")
rep_ge<-c("MT-TL2:Mt_tRNA","RN7SL1:antisense","IL8:protein_coding","HBA1:protein_coding","TMSB4X:protein_coding","EIF4A1:retained_intron")
targets<-unique(rbind(
  data.frame(feature=strict$feature,layer=strict$layer,strict_core=TRUE,original_representative=strict$feature%in%c(rep_mi,rep_ge)),
  data.frame(feature=c(rep_mi,rep_ge),layer=c(rep("miRNA",length(rep_mi)),rep("GENCODE",length(rep_ge))),strict_core=c(rep_mi,rep_ge)%in%strict$feature,original_representative=TRUE)))
targets<-targets[order(targets$layer,targets$feature),]
stopifnot(sum(targets$layer=="miRNA")==19,sum(targets$layer=="GENCODE")==35)
write.csv(targets,file.path(out,"intermediate","STAGE3_TARGET_FEATURES.csv"),row.names=FALSE,na="")

dir_sign<-function(x)ifelse(is.na(x),NA_character_,ifelse(x>sqrt(.Machine$double.eps),"positive",ifelse(x< -sqrt(.Machine$double.eps),"negative","zero")))
bind_fill<-function(xs){cols<-unique(unlist(lapply(xs,names)));do.call(rbind,lapply(xs,function(x){for(m in setdiff(cols,names(x)))x[[m]]<-NA;x[,cols,drop=FALSE]}))}

fit_once<-function(dat,structure="B",optimizer="nloptwrap"){
  f<-if(structure=="B")y~centered_label+(1+centered_label||subject) else y~centered_label+(0+centered_label|subject)
  ctrl<-if(optimizer=="nloptwrap")list(maxeval=200000)else list(maxfun=200000); warns<-character()
  fit<-tryCatch(withCallingHandlers(lmerTest::lmer(f,data=dat,REML=TRUE,control=lmerControl(optimizer=optimizer,optCtrl=ctrl,check.conv.singular=lme4:::.makeCC(action="ignore",tol=tol))),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
  if(inherits(fit,"error"))return(list(fit=NULL,ok=FALSE,error=conditionMessage(fit)))
  co<-tryCatch(summary(fit)$coefficients,error=function(e)NULL);if(is.null(co)||!"centered_label"%in%rownames(co))return(list(fit=fit,ok=FALSE,error="coefficient missing"))
  msgs<-fit@optinfo$conv$lme4$messages;code<-fit@optinfo$conv$opt;code_bad<-length(code)&&any(suppressWarnings(as.numeric(code))!=0,na.rm=TRUE)
  est<-unname(co["centered_label","Estimate"]);se<-unname(co["centered_label","Std. Error"]);df<-if("df"%in%colnames(co))unname(co["centered_label","df"])else NA_real_;pc<-grep("Pr\\(",colnames(co),value=TRUE)[1];pv<-if(!is.na(pc))unname(co["centered_label",pc])else NA_real_
  vc<-as.data.frame(VarCorr(fit));getv<-function(v){q<-vc$vcov[vc$var1==v&is.na(vc$var2)];if(length(q))q[1]else NA_real_};ri<-getv("(Intercept)");rs<-getv("centered_label");rv<-vc$vcov[vc$grp=="Residual"][1]
  rir<-if(is.finite(ri)&&rv>0)sqrt(max(ri,0)/rv)else NA_real_;rsr<-if(is.finite(rs)&&rv>0)sqrt(max(rs,0)/rv)else NA_real_
  list(fit=fit,ok=is.finite(est)&&is.finite(se)&&se>0&&!code_bad&&!length(msgs),error="",estimate=est,SE=se,df=df,P=pv,CI_lower=est-1.96*se,CI_upper=est+1.96*se,RI_SD_ratio=rir,RS_SD_ratio=rsr,singular=isSingular(fit,tol=tol),warnings=paste(unique(warns),collapse=" | "))
}

fit_modelB<-function(y,meta){
  dat<-data.frame(y=as.numeric(y),centered_label=as.numeric(meta$centered_label),subject=factor(meta$subject));dat<-dat[complete.cases(dat),]
  fits<-lapply(optimizers,function(o)fit_once(dat,"B",o));valid<-vapply(fits,function(q)isTRUE(q$ok),logical(1));idx<-if(any(valid))which(valid)[1]else which(vapply(fits,function(q)!is.null(q$fit),logical(1)))[1]
  if(!length(idx)||is.na(idx))return(data.frame(n=nrow(dat),n_subjects=nlevels(dat$subject),slope=NA,SE=NA,df=NA,P=NA,CI_lower=NA,CI_upper=NA,model_status="MODEL_NOT_EVALUABLE",random_slope_validated=FALSE,selected_optimizer=NA,optimizer_valid_n=0,optimizer_stable=FALSE,RI_SD_ratio=NA,RS_SD_ratio=NA,CR2_slope=NA,CR2_SE=NA,CR2_df=NA,CR2_P=NA,CR2_CI_lower=NA,CR2_CI_upper=NA))
  q<-fits[[idx]];vals<-fits[valid];stable<-FALSE
  if(length(vals)>=2){es<-vapply(vals,`[[`,numeric(1),"estimate");ss<-vapply(vals,`[[`,numeric(1),"SE");stable<-length(unique(dir_sign(es)))==1&&diff(range(es))/median(ss)<=.10&&max(ss)/min(ss)<=1.10}
  ri0<-is.finite(q$RI_SD_ratio)&&q$RI_SD_ratio<=tol;rs0<-is.finite(q$RS_SD_ratio)&&q$RS_SD_ratio<=tol
  status<-if(length(vals)<2||!stable)"NUMERICAL_UNSTABLE"else if(ri0&&rs0)"BOTH_VARIANCE_BOUNDARY_UNRESOLVED"else if(rs0)"RS_BOUNDARY"else if(ri0){cf<-lapply(optimizers,function(o)fit_once(dat,"C",o));cv<-cf[vapply(cf,function(x)isTRUE(x$ok),logical(1))];ok<-length(cv)>=2;if(ok){ce<-vapply(cv,`[[`,numeric(1),"estimate");cs<-vapply(cv,`[[`,numeric(1),"SE");ok<-length(unique(dir_sign(ce)))==1&&diff(range(ce))/median(cs)<=.10&&max(cs)/min(cs)<=1.10&&abs(q$estimate-ce[1])/median(c(q$SE,cs[1]))<=.10};if(ok)"RI_BOUNDARY_C_CONFIRMED"else"RI_BOUNDARY_C_UNCONFIRMED"}else"REGULAR_RANDOM_SLOPE"
  officialP<-if(status%in%c("NUMERICAL_UNSTABLE","MODEL_NOT_EVALUABLE"))NA_real_ else q$P
  cr<-tryCatch({lmf<-lm(y~centered_label+subject,data=dat);ct<-coef_test(lmf,vcov="CR2",cluster=dat$subject,test="Satterthwaite");z<-ct[rownames(ct)=="centered_label",,drop=FALSE];crit<-qt(.975,z$df_Satt);c(slope=z$beta,SE=z$SE,df=z$df_Satt,P=z$p_Satt,lo=z$beta-crit*z$SE,hi=z$beta+crit*z$SE)},error=function(e)c(slope=NA,SE=NA,df=NA,P=NA,lo=NA,hi=NA))
  data.frame(n=nrow(dat),n_subjects=nlevels(dat$subject),slope=q$estimate,SE=q$SE,df=q$df,P=officialP,CI_lower=q$CI_lower,CI_upper=q$CI_upper,model_status=status,random_slope_validated=status%in%c("REGULAR_RANDOM_SLOPE","RI_BOUNDARY_C_CONFIRMED"),selected_optimizer=optimizers[idx],optimizer_valid_n=sum(valid),optimizer_stable=stable,RI_SD_ratio=q$RI_SD_ratio,RS_SD_ratio=q$RS_SD_ratio,CR2_slope=cr["slope"],CR2_SE=cr["SE"],CR2_df=cr["df"],CR2_P=cr["P"],CR2_CI_lower=cr["lo"],CR2_CI_upper=cr["hi"])
}

fit_joint<-function(ycsf,yplasma,meta){
  d<-rbind(data.frame(y=as.numeric(ycsf),centered_label=meta$centered_label,subject=factor(meta$subject),biofluid="CSF"),data.frame(y=as.numeric(yplasma),centered_label=meta$centered_label,subject=factor(meta$subject),biofluid="plasma"));d$biofluid<-factor(d$biofluid,levels=c("plasma","CSF"));d$subject_biofluid<-interaction(d$subject,d$biofluid,drop=TRUE);d$csf_time<-d$centered_label*(d$biofluid=="CSF")
  fit<-lm(y~0+subject_biofluid+centered_label+csf_time,data=d);ct<-coef_test(fit,vcov="CR2",cluster=d$subject,test="Satterthwaite");get<-function(n)ct[rownames(ct)==n,,drop=FALSE];p<-get("centered_label");di<-get("csf_time");crit<-qt(.975,di$df_Satt)
  gd<-data.frame(y=as.numeric(ycsf-yplasma),centered_label=meta$centered_label,subject=factor(meta$subject));gf<-lm(y~centered_label+subject,gd);gc<-coef_test(gf,vcov="CR2",cluster=gd$subject,test="Satterthwaite");gz<-gc[rownames(gc)=="centered_label",,drop=FALSE]
  data.frame(n_pairs=nrow(meta),n_subjects=length(unique(meta$subject)),joint_plasma_slope=p$beta,joint_plasma_SE=p$SE,joint_plasma_df=p$df_Satt,joint_plasma_P=p$p_Satt,interaction_CSF_minus_plasma=di$beta,interaction_SE=di$SE,interaction_df=di$df_Satt,interaction_P=di$p_Satt,interaction_CI_lower=di$beta-crit*di$SE,interaction_CI_upper=di$beta+crit*di$SE,joint_CSF_slope=p$beta+di$beta,gradient_CR2_slope=gz$beta,gradient_CR2_SE=gz$SE,gradient_CR2_P=gz$p_Satt,joint_minus_gradient_CR2=di$beta-gz$beta,data_identity_max_abs=max(abs((ycsf-yplasma)-gd$y)))
}

make_scale<-function(inp,scale){
  rc<-log2(inp$csf+1);rp<-log2(inp$plasma+1)
  if(scale=="RAW")return(list(csf=rc,plasma=rp,gradient=rc-rp))
  list(csf=sweep(rc,2,colMeans(rc),"-"),plasma=sweep(rp,2,colMeans(rp),"-"),gradient=inp$clr)
}

## all separate and joint models
slope_rows<-list();joint_rows<-list();si<-ji<-0
for(layer in c("miRNA","GENCODE"))for(dataset in c("FULL","QC_PASS_ONLY"))for(scale in c("RAW","CLR")){
  inp<-obj$inputs[[paste0(layer,"_",dataset)]];m<-make_scale(inp,scale);tf<-targets$feature[targets$layer==layer]
  if(scale=="CLR")stopifnot(max(abs((m$csf-m$plasma)-m$gradient))<1e-10)
  for(f in tf){
    for(fl in c("CSF","plasma")){fm<-fit_modelB(if(fl=="CSF")m$csf[f,]else m$plasma[f,],inp$pairs);si<-si+1;slope_rows[[si]]<-cbind(targets[targets$feature==f&targets$layer==layer,],data.frame(dataset=dataset,scale=scale,biofluid=fl),fm)}
    fj<-fit_joint(m$csf[f,],m$plasma[f,],inp$pairs);s1<-stage1_all[stage1_all$feature==f&stage1_all$layer==layer,];prefix<-if(scale=="RAW")"M3"else"M4";ji<-ji+1;joint_rows[[ji]]<-cbind(targets[targets$feature==f&targets$layer==layer,],data.frame(dataset=dataset,scale=scale),fj,stage1B_gradient_slope=s1[[paste0(prefix,"_slope")]],stage1B_gradient_SE=s1[[paste0(prefix,"_SE")]],stage1B_gradient_FDR=s1[[paste0(prefix,"_FDR")]],stage1B_gradient_CR2_slope=s1[[paste0(prefix,"_CR2_slope")]],stage1B_gradient_CR2_FDR=s1[[paste0(prefix,"_CR2_FDR")]])
  }
}
slopes<-do.call(rbind,slope_rows);joint<-do.call(rbind,joint_rows)
for(layer in c("miRNA","GENCODE"))for(dataset in c("FULL","QC_PASS_ONLY"))for(scale in c("RAW","CLR"))for(fl in c("CSF","plasma")){
  q<-slopes$layer==layer&slopes$dataset==dataset&slopes$scale==scale&slopes$biofluid==fl;n<-sum(q);slopes$FDR[q]<-p.adjust(slopes$P[q],"BH",n=n);slopes$CR2_FDR[q]<-p.adjust(slopes$CR2_P[q],"BH",n=n)
}
for(layer in c("miRNA","GENCODE"))for(dataset in c("FULL","QC_PASS_ONLY"))for(scale in c("RAW","CLR")){
  q<-joint$layer==layer&joint$dataset==dataset&joint$scale==scale;n<-sum(q);joint$interaction_FDR[q]<-p.adjust(joint$interaction_P[q],"BH",n=n)
}
write.csv(slopes,file.path(out,"STAGE3_CSF_PLASMA_SLOPES.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(joint,file.path(out,"STAGE3_BIOFLUID_INTERACTION.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(slopes[,c("feature","layer","dataset","scale","biofluid","model_status","random_slope_validated","optimizer_valid_n","optimizer_stable","RI_SD_ratio","RS_SD_ratio")],file.path(out,"intermediate","STAGE3_MODEL_DIAGNOSTICS.csv"),row.names=FALSE,na="")

## label summaries for strict core profiles
make_profiles<-function(layer){
  fs<-strict$feature[strict$layer==layer];rows<-list();k<-0
  for(dataset in c("FULL","QC_PASS_ONLY"))for(scale in c("RAW","CLR")){
    inp<-obj$inputs[[paste0(layer,"_",dataset)]];m<-make_scale(inp,scale)
    for(f in fs)for(fl in c("CSF","plasma"))for(lab in sort(unique(inp$pairs$collection_label))){idx<-inp$pairs$collection_label==lab;y<-if(fl=="CSF")m$csf[f,idx]else m$plasma[f,idx];k<-k+1;rows[[k]]<-data.frame(feature=f,layer=layer,dataset=dataset,scale=scale,biofluid=fl,collection_label=lab,label_numeric=as.integer(sub("D","",lab)),n=sum(idx),n_subjects=length(unique(inp$pairs$subject[idx])),median=median(y),Q1=unname(quantile(y,.25)),Q3=unname(quantile(y,.75)),mean=mean(y),SE=sd(y)/sqrt(length(y)))}
  }
  x<-do.call(rbind,rows);merge(x,slopes,by=c("feature","layer","dataset","scale","biofluid"),all.x=TRUE)
}
core_mi_profiles<-make_profiles("miRNA");core_ge_profiles<-make_profiles("GENCODE")
write.csv(core_mi_profiles,file.path(out,"STAGE3_CORE_MIRNA_PROFILES.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(core_ge_profiles,file.path(out,"STAGE3_GENCODE_PROFILES.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## decomposition classification
classify<-function(csf,plasma,diff,lo,hi){
  if(any(!is.finite(c(csf,plasma,diff,lo,hi)))||lo<=0&&hi>=0)return("E_UNCLEAR")
  if(diff<0){
    if(csf<0&&plasma>0){r<-(-csf)/(plasma);if(r>=2)return("A_CSF_DECREASE_DOMINANT");if(r<=.5)return("B_PLASMA_INCREASE_DOMINANT");return("C_CSF_DECREASE_AND_PLASMA_INCREASE")}
    if(csf<0&&plasma<=0){if(abs(csf)>=2*abs(plasma))return("A_CSF_DECREASE_DOMINANT");return("D_SAME_DIRECTION_DIFFERENT_MAGNITUDE")}
    if(csf>=0&&plasma>0){if(abs(plasma)>=2*abs(csf))return("B_PLASMA_INCREASE_DOMINANT");return("D_SAME_DIRECTION_DIFFERENT_MAGNITUDE")}
    return("E_UNCLEAR")
  }
  if(csf>0&&plasma<0){r<-csf/(-plasma);if(r>=2)return("A_CSF_CHANGE_DOMINANT");if(r<=.5)return("B_PLASMA_CHANGE_DOMINANT");return("C_OPPOSITE_DIRECTIONS_REINFORCING")}
  if(sign(csf)==sign(plasma))return("D_SAME_DIRECTION_DIFFERENT_MAGNITUDE")
  if(abs(csf)>=2*abs(plasma))"A_CSF_CHANGE_DOMINANT"else if(abs(plasma)>=2*abs(csf))"B_PLASMA_CHANGE_DOMINANT"else"E_UNCLEAR"
}
decomp_rows<-list();di<-0
for(scale in c("RAW","CLR")){
  j<-joint[joint$dataset=="FULL"&joint$scale==scale,];s<-slopes[slopes$dataset=="FULL"&slopes$scale==scale,]
  sw<-reshape(s[,c("feature","layer","biofluid","slope","SE","CI_lower","CI_upper","P","FDR","CR2_slope","CR2_SE","CR2_P","CR2_FDR","model_status")],idvar=c("feature","layer"),timevar="biofluid",direction="wide")
  x<-merge(j,sw,by=c("feature","layer"));x$contribution_class<-mapply(classify,x$slope.CSF,x$slope.plasma,x$interaction_CSF_minus_plasma,x$interaction_CI_lower,x$interaction_CI_upper);x$algebraic_slope_difference<-x$slope.CSF-x$slope.plasma;x$modelB_difference_minus_stage1B<-x$algebraic_slope_difference-x$stage1B_gradient_slope;di<-di+1;decomp_rows[[di]]<-x
}
decomp<-do.call(rbind,decomp_rows)
write.csv(decomp,file.path(out,"STAGE3_GRADIENT_DECOMPOSITION.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## RAW vs CLR sensitivity
raws<-slopes[slopes$scale=="RAW",];clrs<-slopes[slopes$scale=="CLR",];clr_sens<-merge(raws,clrs,by=c("feature","layer","strict_core","original_representative","dataset","biofluid"),suffixes=c("_RAW","_CLR"),all=TRUE);clr_sens$direction_same<-dir_sign(clr_sens$slope_RAW)==dir_sign(clr_sens$slope_CLR);clr_sens$slope_change_CLR_minus_RAW<-clr_sens$slope_CLR-clr_sens$slope_RAW
jr<-joint[joint$scale=="RAW",];jc<-joint[joint$scale=="CLR",];jcomp<-merge(jr,jc,by=c("feature","layer","strict_core","original_representative","dataset"),suffixes=c("_RAW","_CLR"));jcomp$record_type<-"INTERACTION";clr_sens$record_type<-"BIOFLUID"
write.csv(bind_fill(list(clr_sens,jcomp)),file.path(out,"STAGE3_CLR_SENSITIVITY.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## leave-one-subject-out fixed-subject slope influence for strict miRNAs and all representatives
loo_targets<-targets[targets$strict_core&targets$layer=="miRNA"|targets$original_representative,];loo_rows<-list();li<-0
for(i in seq_len(nrow(loo_targets))){layer<-loo_targets$layer[i];f<-loo_targets$feature[i]
  for(scale in c("RAW","CLR")){inp<-obj$inputs[[paste0(layer,"_FULL")]];m<-make_scale(inp,scale)
    for(fl in c("CSF","plasma","gradient")){y<-switch(fl,CSF=m$csf[f,],plasma=m$plasma[f,],gradient=m$gradient[f,]);base<-coef(lm(y~centered_label+factor(subject),data=data.frame(y=y,centered_label=inp$pairs$centered_label,subject=inp$pairs$subject)))["centered_label"]
      vals<-sapply(unique(inp$pairs$subject),function(sub){q<-inp$pairs$subject!=sub;coef(lm(y[q]~inp$pairs$centered_label[q]+factor(inp$pairs$subject[q])))[2]})
      li<-li+1;loo_rows[[li]]<-data.frame(feature=f,layer=layer,scale=scale,biofluid=fl,base_fixed_subject_slope=base,LOO_min=min(vals),LOO_max=max(vals),LOO_max_abs_change=max(abs(vals-base)),LOO_direction_flip=any(dir_sign(vals)!=dir_sign(base)),n_exclusions=length(vals))
    }
  }
}
loo<-do.call(rbind,loo_rows);write.csv(loo,file.path(out,"intermediate","STAGE3_LOO_PATIENT_INFLUENCE.csv"),row.names=FALSE,na="")

## FULL vs QC comparisons
sf<-slopes[slopes$dataset=="FULL",];sq<-slopes[slopes$dataset=="QC_PASS_ONLY",];qc_s<-merge(sf,sq,by=c("feature","layer","strict_core","original_representative","scale","biofluid"),suffixes=c("_FULL","_QC"),all=TRUE);qc_s$direction_same<-dir_sign(qc_s$slope_FULL)==dir_sign(qc_s$slope_QC);qc_s$absolute_slope_change<-abs(qc_s$slope_QC-qc_s$slope_FULL);qc_s$record_type<-"BIOFLUID"
jf<-joint[joint$dataset=="FULL",];jq<-joint[joint$dataset=="QC_PASS_ONLY",];qc_j<-merge(jf,jq,by=c("feature","layer","strict_core","original_representative","scale"),suffixes=c("_FULL","_QC"),all=TRUE);qc_j$direction_same<-dir_sign(qc_j$interaction_CSF_minus_plasma_FULL)==dir_sign(qc_j$interaction_CSF_minus_plasma_QC);qc_j$absolute_slope_change<-abs(qc_j$interaction_CSF_minus_plasma_QC-qc_j$interaction_CSF_minus_plasma_FULL);qc_j$record_type<-"INTERACTION"
qc_out<-bind_fill(list(qc_s,qc_j));qc_out<-merge(qc_out,loo[loo$biofluid!="gradient",],by=c("feature","layer","scale","biofluid"),all.x=TRUE)
loo_gradient<-loo[loo$biofluid=="gradient",];loo_gradient$record_type<-"LOO_GRADIENT"
qc_out<-bind_fill(list(qc_out,loo_gradient))
write.csv(qc_out,file.path(out,"STAGE3_QC_SENSITIVITY.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## D01/D07 endpoint complete subjects
pp<-obj$inputs$miRNA_FULL$pairs;tab<-table(pp$subject,pp$collection_label);endpoint_subjects<-rownames(tab)[tab[,"D01"]>0&tab[,"D07"]>0];complete7<-rownames(tab)[rowSums(tab[,paste0("D0",1:7),drop=FALSE]>0)==7];stopifnot(length(endpoint_subjects)==12,length(complete7)==4)
endpoint_rows<-list();ei<-0
for(layer in c("miRNA","GENCODE"))for(scale in c("RAW","CLR")){
  inp<-obj$inputs[[paste0(layer,"_FULL")]];m<-make_scale(inp,scale);tf<-targets$feature[targets$layer==layer]
  i1<-match(paste(endpoint_subjects,"D01"),paste(inp$pairs$subject,inp$pairs$collection_label));i7<-match(paste(endpoint_subjects,"D07"),paste(inp$pairs$subject,inp$pairs$collection_label));stopifnot(!anyNA(i1),!anyNA(i7))
  for(f in tf){dc<-m$csf[f,i7]-m$csf[f,i1];dp<-m$plasma[f,i7]-m$plasma[f,i1];dg<-dc-dp
    for(outcome in c("CSF","plasma","gradient")){v<-switch(outcome,CSF=dc,plasma=dp,gradient=dg);for(k in seq_along(endpoint_subjects)){ei<-ei+1;endpoint_rows[[ei]]<-data.frame(record_type="SUBJECT_CHANGE",feature=f,layer=layer,scale=scale,outcome=outcome,subject=endpoint_subjects[k],change_D07_minus_D01=v[k],n=12)};crit<-qt(.975,11);ei<-ei+1;endpoint_rows[[ei]]<-data.frame(record_type="SUMMARY",feature=f,layer=layer,scale=scale,outcome=outcome,subject=NA,change_D07_minus_D01=NA,n=12,mean_change=mean(v),SE=sd(v)/sqrt(12),CI_lower=mean(v)-crit*sd(v)/sqrt(12),CI_upper=mean(v)+crit*sd(v)/sqrt(12),median_change=median(v),proportion_positive=mean(v>0),P=t.test(v)$p.value)}
  }
}
endpoint<-bind_fill(endpoint_rows)
for(layer in c("miRNA","GENCODE"))for(scale in c("RAW","CLR"))for(outcome in c("CSF","plasma","gradient")){q<-endpoint$record_type=="SUMMARY"&endpoint$layer==layer&endpoint$scale==scale&endpoint$outcome==outcome;endpoint$FDR[q]<-p.adjust(endpoint$P[q],"BH",n=sum(q))}
ep_sum<-endpoint[endpoint$record_type=="SUMMARY",];ep_sum$full_model_slope<-NA_real_;ep_sum$direction_agreement_FULL<-NA
for(i in seq_len(nrow(ep_sum))){if(ep_sum$outcome[i]=="gradient"){q<-joint$feature==ep_sum$feature[i]&joint$layer==ep_sum$layer[i]&joint$dataset=="FULL"&joint$scale==ep_sum$scale[i];sl<-joint$stage1B_gradient_slope[q]}else{q<-slopes$feature==ep_sum$feature[i]&slopes$layer==ep_sum$layer[i]&slopes$dataset=="FULL"&slopes$scale==ep_sum$scale[i]&slopes$biofluid==ep_sum$outcome[i];sl<-slopes$slope[q]};ep_sum$full_model_slope[i]<-sl;ep_sum$direction_agreement_FULL[i]<-dir_sign(ep_sum$mean_change[i])==dir_sign(sl)}
endpoint<-merge(endpoint,ep_sum[,c("feature","layer","scale","outcome","full_model_slope","direction_agreement_FULL")],by=c("feature","layer","scale","outcome"),all.x=TRUE,suffixes=c("",".summary"));write.csv(endpoint,file.path(out,"STAGE3_ENDPOINT_SENSITIVITY.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(data.frame(endpoint_complete_subjects=endpoint_subjects),file.path(out,"intermediate","STAGE3_ENDPOINT_SUBJECTS.csv"),row.names=FALSE)
write.csv(data.frame(complete_7label_subjects=complete7),file.path(out,"intermediate","STAGE3_COMPLETE7_SUBJECTS.csv"),row.names=FALSE)

## representative table
reps<-targets[targets$original_representative,];rep_out<-merge(reps,decomp,by=c("feature","layer","strict_core","original_representative"),all.x=TRUE);rec<-read.csv("revision_stage1/stage1b_full_refit/STAGE1B_REPRESENTATIVE_FEATURES.csv",check.names=FALSE);rep_out<-merge(rep_out,rec[,c("feature","representative_recommendation")],by="feature",all.x=TRUE);write.csv(rep_out,file.path(out,"STAGE3_REPRESENTATIVE_FEATURES.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

## figures from true observations
build_plot_data<-function(layer,features){inp<-obj$inputs[[paste0(layer,"_FULL")]];m<-make_scale(inp,"RAW");do.call(rbind,lapply(features,function(f)rbind(data.frame(feature=f,subject=inp$pairs$subject,label=as.integer(sub("D","",inp$pairs$collection_label)),outcome="CSF",value=as.numeric(m$csf[f,])),data.frame(feature=f,subject=inp$pairs$subject,label=as.integer(sub("D","",inp$pairs$collection_label)),outcome="plasma",value=as.numeric(m$plasma[f,])),data.frame(feature=f,subject=inp$pairs$subject,label=as.integer(sub("D","",inp$pairs$collection_label)),outcome="gradient",value=as.numeric(m$gradient[f,])))))}
save_both<-function(name,p,width,height){svglite::svglite(file.path(out,"plots",paste0(name,".svg")),width=width,height=height);print(p);dev.off();Cairo::CairoPDF(file.path(out,"plots",paste0(name,".pdf")),width=width,height=height);print(p);dev.off()}
plot_trajectories<-function(layer,features,title,name){d<-build_plot_data(layer,features);avg<-aggregate(value~feature+outcome+label,d,mean);p<-ggplot(d,aes(label,value,group=subject))+geom_line(alpha=.12,color="grey40")+geom_point(alpha=.12,size=.5,color="grey40")+geom_line(data=avg,aes(x=label,y=value,group=1),inherit.aes=FALSE,color="#0072B2",linewidth=.9)+geom_point(data=avg,aes(x=label,y=value),inherit.aes=FALSE,color="#0072B2",size=1.4)+facet_grid(feature~outcome,scales="free_y")+scale_x_continuous(breaks=1:7,labels=paste0("D0",1:7))+labs(x="Collection label",y="Relative abundance: log2(RPM+1) or paired difference",title=title,caption="FULL: 30 subjects, 153 paired observations; thin lines are patient trajectories; blue line is label mean")+theme_bw()+theme(strip.text.y=element_text(angle=0));save_both(name,p,13,2+1.55*length(features))}
plot_trajectories("miRNA",rep_mi,"Original Figure 3 miRNA: separate biofluid trajectories","STAGE3_FIGURE_A_MIRNA_TRIPLE")
plot_trajectories("GENCODE",rep_ge,"Original Figure 4 features: separate biofluid trajectories","STAGE3_FIGURE_C_GENCODE_TRIPLE")

dm<-decomp[decomp$layer=="miRNA"&decomp$strict_core&decomp$scale=="RAW",];forest<-rbind(data.frame(feature=dm$feature,outcome="CSF",slope=dm$slope.CSF,lo=dm$CI_lower.CSF,hi=dm$CI_upper.CSF,class=dm$contribution_class),data.frame(feature=dm$feature,outcome="plasma",slope=dm$slope.plasma,lo=dm$CI_lower.plasma,hi=dm$CI_upper.plasma,class=dm$contribution_class),data.frame(feature=dm$feature,outcome="gradient",slope=dm$stage1B_gradient_slope,lo=dm$stage1B_gradient_slope-1.96*dm$stage1B_gradient_SE,hi=dm$stage1B_gradient_slope+1.96*dm$stage1B_gradient_SE,class=dm$contribution_class));forest$feature<-factor(forest$feature,levels=rev(dm$feature[order(dm$stage1B_gradient_slope)]));pB<-ggplot(forest,aes(slope,feature,color=outcome))+geom_vline(xintercept=0,linetype=2)+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.18)+geom_point()+facet_wrap(~outcome,nrow=1)+labs(x="Slope per collection label (95% CI)",y=NULL,title="Strict-core miRNA biofluid and paired-gradient slopes",caption="FULL RAW; selected features are conditional on Stage 1B gradient screening")+theme_bw()+theme(legend.position="none");save_both("STAGE3_FIGURE_B_CORE_MIRNA_FOREST",pB,12,8)

sens<-rbind(data.frame(feature=slopes$feature[slopes$layer=="miRNA"&slopes$strict_core],condition=paste(slopes$dataset[slopes$layer=="miRNA"&slopes$strict_core],slopes$scale[slopes$layer=="miRNA"&slopes$strict_core],sep="_"),outcome=slopes$biofluid[slopes$layer=="miRNA"&slopes$strict_core],slope=slopes$slope[slopes$layer=="miRNA"&slopes$strict_core]),data.frame(feature=joint$feature[joint$layer=="miRNA"&joint$strict_core],condition=paste(joint$dataset[joint$layer=="miRNA"&joint$strict_core],joint$scale[joint$layer=="miRNA"&joint$strict_core],sep="_"),outcome="gradient",slope=joint$stage1B_gradient_slope[joint$layer=="miRNA"&joint$strict_core]));pD<-ggplot(sens,aes(slope,feature,color=condition,shape=outcome))+geom_vline(xintercept=0,linetype=2)+geom_point(position=position_dodge(width=.6),size=1.8)+labs(x="Slope per collection label",y=NULL,title="Strict-core miRNA RAW/CLR and QC sensitivity",caption="FULL versus paired QC PASS-only; CLR is a re-expression of the same RPM data")+theme_bw()+theme(legend.position="bottom");save_both("STAGE3_FIGURE_D_QC_CLR_SENSITIVITY",pD,12,8)

saveRDS(list(targets=targets,slopes=slopes,interactions=joint,decomposition=decomp,CLR_sensitivity=clr_sens,QC_sensitivity=qc_out,endpoint=endpoint,LOO=loo),file.path(out,"intermediate","STAGE3_RESULTS_OBJECT.rds"))
capture.output(sessionInfo(),file=file.path(out,"logs","STAGE3_sessionInfo.txt"))
cat(sprintf("STAGE3_ANALYSIS_COMPLETE targets=%d slopes=%d interactions=%d endpoint_subjects=%d complete7=%d\n",nrow(targets),nrow(slopes),nrow(joint),length(endpoint_subjects),length(complete7)))
