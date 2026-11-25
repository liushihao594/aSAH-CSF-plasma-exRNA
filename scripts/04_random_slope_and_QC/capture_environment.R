options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(limma)
  library(statmod)
  library(lme4)
  library(lmerTest)
  library(pbkrtest)
})
out <- "revision_stage1/logs/STAGE1_sessionInfo.txt"
capture.output({
  cat("Stage 1 environment captured: 2026-11-21 Asia/Shanghai\n")
  cat("Analysis package library search path:\n")
  print(.libPaths())
  cat("\nKey package versions:\n")
  for (p in c("limma", "statmod", "lme4", "lmerTest", "pbkrtest")) {
    cat(p, as.character(packageVersion(p)), "\n")
  }
  cat("\n")
  print(sessionInfo())
}, file = out)
