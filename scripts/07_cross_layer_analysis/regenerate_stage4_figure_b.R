options(stringsAsFactors = FALSE)
.libPaths(c("revision_stage1/R_library", .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
root <- "revision_stage4"
scores <- read.csv(file.path(root, "STAGE4_RAW_CLR_CROSSFITTED_SCORES.csv"), check.names = FALSE)
patient_slopes <- function(sc) do.call(rbind, lapply(unique(sc$subject), function(id) {
  z <- sc[sc$subject == id, , drop = FALSE]
  dx <- z$day - mean(z$day)
  data.frame(original_subject = id,
             miRNA_slope = sum(dx * z$miRNA_score) / sum(dx^2),
             GENCODE_slope = sum(dx * z$GENCODE_score) / sum(dx^2))
}))
slope_plot <- do.call(rbind, lapply(c("RAW", "CLR"), function(scale) {
  z <- patient_slopes(scores[scores$dataset == "FULL" & scores$scale == scale, ])
  z$scale <- scale
  z
}))
p <- ggplot(slope_plot, aes(miRNA_slope, GENCODE_slope)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  geom_point(color = "#0072B2", size = 2) +
  geom_text(data = slope_plot[slope_plot$original_subject == "EXR-KJENS1Sample12-DO", ],
            aes(label = "subject12"), hjust = 0, vjust = -1, nudge_x = .01, color = "#D55E00") +
  geom_smooth(method = "lm", se = FALSE, color = "grey30") + facet_wrap(~scale, scales = "free") +
  theme_bw() + labs(title = "Subject-specific cross-layer score slopes",
                    x = "miRNA slope", y = "GENCODE slope",
                    caption = "30 subjects; each slope uses the subject's 3-7 observed collection labels")
ggsave(file.path(root, "plots", "STAGE4_FIGURE_B_PATIENT_SLOPES.pdf"), p,
       width = 9, height = 5, device = cairo_pdf)
ggsave(file.path(root, "plots", "STAGE4_FIGURE_B_PATIENT_SLOPES.svg"), p,
       width = 9, height = 5)
