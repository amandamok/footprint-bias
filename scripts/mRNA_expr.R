rm(list=ls())

library(here)
library(ggplot2)

# load mRNA RPKM from Weinberg et al. (2016)
mrna_expr <- read.table(file.path(here(), "reference_data", "GSE75897_RiboZero_RPKMs.txt"),
                        stringsAsFactors=F, col.names=c("transcript", "RPKM"))

# load regression model
load(file.path(here(), "expts", "weinberg_20cds20_f5_2_f3_3", "weinberg_fit_150.Rda"))
model_coef <- data.frame(summary(weinberg_fit_150)$coefficients)
model_coef <- subset(model_coef, grepl("transcript", rownames(model_coef)))
rownames(model_coef) <- sub("transcript", "", rownames(model_coef))

# plot
common_transcripts <- rownames(model_coef)[rownames(model_coef) %in% mrna_expr$transcript]
model_coef <- model_coef[match(common_transcripts, rownames(model_coef)),]
plot_data <- data.frame(transcript = common_transcripts,
                        RPKM = mrna_expr$RPKM[match(common_transcripts, mrna_expr$transcript)],
                        FC = with(model_coef, exp(Estimate)),
                        FC_lower = with(model_coef, exp(Estimate - Std..Error)),
                        FC_upper = with(model_coef, exp(Estimate + Std..Error)))
mrna_text <- data.frame(cor = paste("rho =", signif(cor(plot_data$RPKM, plot_data$FC),2)),
                        pvalue = paste("p =", signif(cor.test(plot_data$RPKM, plot_data$FC)$p.value, 2)))
mrna_plot <- ggplot(plot_data, aes(x=RPKM, y=FC)) + geom_point() + theme_classic() +
  geom_smooth(method="lm", formula=y~x) + geom_errorbar(aes(ymin=FC_lower, ymax=FC_upper)) +
  ggtitle("Weinberg", subtitle=paste(nrow(plot_data), "transcripts")) +
  ylab(expr(paste("exp(", beta, ")"))) +
  geom_text(data=mrna_text, mapping=aes(x=500, y=max(plot_data$FC_upper), label=cor)) +
  geom_text(data=mrna_text, mapping=aes(x=500, y=0.95*max(plot_data$FC_upper), label=pvalue))

load(file.path(here(), "expts", "weinberg_20cds20_f5_2_f3_3", "weinberg_bam.Rda"))
transcript_cts <- aggregate(count~transcript, data=weinberg_bam, FUN=sum)
plot_data$RPF_ct <- transcript_cts$count[match(common_transcripts, transcript_cts$transcript)]
rpf_text <- data.frame(cor = paste("rho =", signif(cor(plot_data$RPF_ct, plot_data$FC),2)),
                        pvalue = paste("p =", signif(cor.test(plot_data$RPF_ct, plot_data$FC)$p.value, 2)))
rpf_plot <- ggplot(plot_data, aes(x=RPF_ct, y=FC)) + geom_point() + theme_classic() +
  geom_smooth(method="lm", formula=y~x) + geom_errorbar(aes(ymin=FC_lower, ymax=FC_upper)) +
  ylab(expr(paste("exp(", beta, ")"))) + xlab("footprint count") +
  geom_text(data=rpf_text, mapping=aes(x=75000, y=0.95, label=cor)) +
  geom_text(data=rpf_text, mapping=aes(x=75000, y=0.9, label=pvalue))

mrna_plot / rpf_plot