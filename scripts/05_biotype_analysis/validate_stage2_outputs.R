options(stringsAsFactors=FALSE,warn=1)
out<-"revision_stage2"
checks<-data.frame(check=character(),status=character(),evidence=character())
add<-function(n,ok,e){checks<<-rbind(checks,data.frame(check=n,status=ifelse(ok,"PASS","FAIL"),evidence=as.character(e)))}

h<-read.csv(file.path(out,"STAGE2_HIGH_ABUNDANCE_RNA_RESULTS.csv"),check.names=FALSE)
b<-read.csv(file.path(out,"STAGE2_BIOTYPE_COMPOSITION.csv"),check.names=FALSE)
g<-read.csv(file.path(out,"STAGE2_GLOBAL_SHIFT_ANALYSIS.csv"),check.names=FALSE)
c<-read.csv(file.path(out,"STAGE2_CORE_FEATURE_ASSESSMENT.csv"),check.names=FALSE)
hm<-read.csv(file.path(out,"intermediate","STAGE2_HIGH_ABUNDANCE_MODELS.csv"),check.names=FALSE)
bm<-read.csv(file.path(out,"intermediate","STAGE2_BIOTYPE_TEMPORAL_MODELS.csv"),check.names=FALSE)

add("high_abundance_rows",nrow(h)==308,nrow(h))
add("high_abundance_fixed_features",length(unique(h$feature))==11,length(unique(h$feature)))
add("high_abundance_full_samples_by_label",all(h$n[h$dataset=="FULL"]>0),paste(range(h$n[h$dataset=="FULL"]),collapse=".."))
add("biotype_rows",nrow(b)==4116,nrow(b))
add("biotype_categories",setequal(unique(b$biotype),c("miRNA","piRNA","tRNA","Y_RNA","lincRNA","snoRNA","snRNA")),paste(sort(unique(b$biotype)),collapse=";"))
bf<-b[b$dataset=="FULL",]
ps<-aggregate(selected_7type_proportion~expression_id,bf,sum)
add("biotype_proportion_sum",max(abs(ps$selected_7type_proportion-1))<1e-12,sprintf("max_error=%.3g",max(abs(ps$selected_7type_proportion-1))))
add("biotype_primary_samples",length(unique(bf$expression_id))==306,length(unique(bf$expression_id)))
add("biotype_qc_samples",length(unique(b$expression_id[b$dataset=="QC_PASS_ONLY"]))==282,length(unique(b$expression_id[b$dataset=="QC_PASS_ONLY"])))
add("global_record_count",nrow(g)==764,nrow(g))
gp<-g[g$record_type=="PAIR",]
add("global_pair_records",nrow(gp)==588,nrow(gp))
add("global_CLR_identity",max(abs(gp$identity_error),na.rm=TRUE)<1e-10,sprintf("max_error=%.3g",max(abs(gp$identity_error),na.rm=TRUE)))
add("global_full_pairs_per_layer",all(table(gp$layer[gp$dataset=="FULL"])==153),paste(table(gp$layer[gp$dataset=="FULL"]),collapse=";"))
add("global_qc_pairs_per_layer",all(table(gp$layer[gp$dataset=="QC_PASS_ONLY"])==141),paste(table(gp$layer[gp$dataset=="QC_PASS_ONLY"]),collapse=";"))
add("core_rows",nrow(c)==49,nrow(c))
add("core_layer_counts",sum(c$layer=="miRNA")==16 && sum(c$layer=="GENCODE")==33,sprintf("miRNA=%d;GENCODE=%d",sum(c$layer=="miRNA"),sum(c$layer=="GENCODE")))
add("core_all_not_cargo_adjusted",all(c$cargo_adjusted_status=="NOT_EVALUATED_XDEC_NO_GO"),paste(unique(c$cargo_adjusted_status),collapse=";"))

for(ds in c("FULL","QC_PASS_ONLY")) for(fl in c("CSF","plasma")){
  q<-hm$dataset==ds & hm$biofluid==fl
  ex<-p.adjust(hm$P[q],"BH",n=sum(q))
  add(paste("high_BH",ds,fl,sep="_"),max(abs(ex-hm$FDR[q]),na.rm=TRUE)<1e-14,sprintf("n=%d",sum(q)))
  q<-bm$dataset==ds & bm$biofluid==fl
  ex<-p.adjust(bm$CLR_P[q],"BH",n=7)
  add(paste("biotype_CLR_BH",ds,fl,sep="_"),max(abs(ex-bm$CLR_FDR[q]),na.rm=TRUE)<1e-14,sprintf("n=%d",sum(q)))
}

required<-c("STAGE2_FROZEN_ANALYSIS_PLAN.md","XDEC_TECHNICAL_VALIDATION.md","XDEC_GO_NO_GO.md",
  "STAGE2_HIGH_ABUNDANCE_RNA_RESULTS.csv","STAGE2_BIOTYPE_COMPOSITION.csv","STAGE2_GLOBAL_SHIFT_ANALYSIS.csv",
  "STAGE2_QC_SENSITIVITY.csv","STAGE2_CORE_FEATURE_ASSESSMENT.csv","STAGE2_FINAL_REPORT.md","STAGE3_RECOMMENDATION.md")
for(f in required)add(paste0("required_",f),file.exists(file.path(out,f)),f)

xdec_forbidden<-c("STAGE2_XDEC_SAMPLE_PROPORTIONS.csv","STAGE2_XDEC_MODEL_DIAGNOSTICS.csv",
  "STAGE2_XDEC_TEMPORAL_RESULTS.csv","STAGE2_CARGO_ADJUSTED_RESULTS.csv")
for(f in xdec_forbidden)add(paste0("xdec_not_created_",f),!file.exists(file.path(out,f)),f)

plots<-list.files(file.path(out,"plots"),pattern="\\.png$",full.names=TRUE)
add("diagnostic_plots",length(plots)==5 && all(file.info(plots)$size>10000),paste(basename(plots),file.info(plots)$size,collapse=";"))

write.csv(checks,file.path(out,"STAGE2_OUTPUT_VALIDATION.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
if(any(checks$status!="PASS")){print(checks[checks$status!="PASS",]);stop("Stage 2 validation failed")}
cat(sprintf("STAGE2_VALIDATION_PASS %d checks\n",nrow(checks)))
