# Technical checks only. Never run XDec main script, EDec fitting, or estimate patient proportions.
options(stringsAsFactors=FALSE)
out<-'revision_stage0c'
code<-'revision_stage0b/04_code_environment/XDec.R'
ast<-parse(code)
allowed<-c('OriginalSet.RNA','OptiPrep.RNA','quantile_normalisation','logistic','Predictions.CTSubtype')
env<-new.env(parent=baseenv())
for(expr in ast){if(is.call(expr)&&as.character(expr[[1]])%in%c('=','<-')&&is.symbol(expr[[2]])&&as.character(expr[[2]])%in%allowed)eval(expr,env)}
ref<-read.delim('revision_stage0b/02_xdec_cargo/Atlas_DeconvolutionResults.txt',row.names=1,check.names=FALSE)
target<-read.csv(file.path(out,'xdec/XDEC_81RNA_306samples_RPM.csv'),row.names=1,check.names=FALSE,fileEncoding='UTF-8-BOM')
stopifnot(identical(rownames(ref),rownames(target)),setequal(rownames(ref),env$OriginalSet.RNA),length(env$Predictions.CTSubtype)==ncol(ref))
singular<-svd(as.matrix(ref))$d
cr<-cor(ref);diag(cr)<-NA
write.csv(data.frame(index=seq_along(singular),singular_value=singular),file.path(out,'xdec/reference_singular_values.csv'),row.names=FALSE)
technical<-data.frame(check=c('source_parse','informative_set_matches_reference','reference_numerical_rank','reference_singular_value_ratio','max_absolute_profile_correlation','target_matrix_NA','raw_target_RPM_max','reference_value_max','XDec_main_executed','patient_proportions_estimated'),value=c(TRUE,TRUE,qr(as.matrix(ref),tol=1e-8)$rank,max(singular)/min(singular),max(abs(cr),na.rm=TRUE),sum(!is.finite(as.matrix(target))),max(target),max(ref),FALSE,FALSE))
# Deterministic synthetic fixture: only establishes transformation code behavior and background dependence.
toy<-matrix(c(0,1,8,2,20,4,1,0,2,7,3,30,0,0,1,0,9,3,6,3,0,7,0,4,2,9,1,3,6,0),nrow=5)
rownames(toy)<-paste0('toyRNA',1:5);colnames(toy)<-paste0('synthetic',1:6)
transformOnly<-function(m){q<-as.data.frame(t(env$quantile_normalisation(t(m))));z<-t(apply(t(q),2,env$logistic));z/max(z)}
a<-transformOnly(toy)
added<-cbind(toy,c(40,0,2,70,1),c(2,80,1,0,9))
b<-transformOnly(added)[,1:6]
technical<-rbind(technical,data.frame(check=c('synthetic_6sample_transform_finite','synthetic_fixture_background_change_max_abs','synthetic_test_is_patient_analysis'),value=c(all(is.finite(a)),max(abs(a-b)),FALSE)))
write.csv(technical,file.path(out,'xdec/technical_test_results.csv'),row.names=FALSE)
installed<-sapply(c('optparse','EDec','gplots','ggplot2','reshape2'),requireNamespace,quietly=TRUE)
write.csv(data.frame(package=names(installed),loadable=installed),file.path(out,'xdec/XDec_dependency_check.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'logs/XDec_technical_sessionInfo.txt'))
print(technical)
