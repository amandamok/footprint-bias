rm(list=ls())

library(ggplot2)
library(here)

# load qPCR data ----------------------------------------------------------

se <- function(x) { sd(x)/sqrt(length(x)) }

oligos <- c("ATA", "TCC", "CCA", "CGT", "GAC", "GGG")
names(oligos) <- c(1:3, 7:9)

# read in qPCR data
qpcr_fname <- file.path(here(), "reference_data", "circligase_qpcr.csv")
raw_data <- read.csv(qpcr_fname, stringsAsFactors=F)
raw_data <- subset(raw_data, as.character(raw_data$Oligo) %in% names(oligos))
raw_data$end_seq <- oligos[as.character(raw_data$Oligo)]
raw_data$primer_type <- ifelse(raw_data$PRIMERS == "NM827_NM828", "circ", "cont")

# calculate efficiency per primer
eff <- aggregate(eff ~ primer_type, data=raw_data, FUN=mean)

# calcaulate amplification
raw_data$primer_eff <- eff$eff[match(raw_data$primer_type, eff$primer)]
raw_data$amp <- raw_data$primer_eff^raw_data$Cy0

# average across replicates
qpcr_data <- aggregate(amp ~ end_seq + primer_type + CL, data=raw_data, FUN=mean)
qpcr_data$amp_se <- aggregate(amp ~ end_seq + primer_type + CL, data=raw_data, FUN=se)$amp
qpcr_data <- reshape(qpcr_data, v.names=c("amp", "amp_se"), idvar=c("end_seq", "CL"),
                     timevar="primer_type", direction="wide")

# calculate ratios
qpcr_data$ratio <- with(qpcr_data, amp.cont/amp.circ)
qpcr_data$ratio <- qpcr_data$ratio / max(qpcr_data$ratio)

# error propagation
qpcr_data$se_percent <- with(qpcr_data,
                             sqrt( (amp_se.circ / amp.circ)^2 + (amp_se.cont / amp.cont)^2 ))
qpcr_data$se_abs <- with(qpcr_data, ratio * se_percent)

# add regression coefficients ---------------------------------------------

plot_data <- lapply(c("green", "lareau", "weinberg"),
                    function(expt) {
                      expt_dir <- paste0(expt, "_20cds20_f5_2_f3_3")
                      # load model fit
                      model_obj <- paste0(expt, "_fit_150")
                      load(file.path(here(), "expts", expt_dir, paste0(model_obj, ".Rda")))
                      # check if 28mer frame 0 is baseline subset
                      d5_ref <- get(model_obj)$xlevels$d5[1]
                      d3_ref <- get(model_obj)$xlevels$d3[1]
                      f5_header <- ifelse(d5_ref=="15" & sum(as.numeric(d5_ref), as.numeric(d3_ref), 3)==28,
                                          "genome_f5", "d515:genome_f5")
                      # pull coefficients from model
                      model_coef <- data.frame(summary(get(model_obj))$coefficients)
                      model_coef <- model_coef[match(paste0(f5_header, substr(oligos, 1, 2)),
                                                     rownames(model_coef)),]
                      model_coef$oligo <- oligos
                      # compute fold-change and upper/lower bounds
                      model_coef$FC <- with(model_coef, exp(Estimate))
                      model_coef$FC_lower <- with(model_coef, exp(Estimate - Std..Error))
                      model_coef$FC_upper <- with(model_coef, exp(Estimate + Std..Error))
                      # return qpcr_data with added coefficients
                      row_index <- match(qpcr_data$end_seq, model_coef$oligo)
                      return(cbind(qpcr_data,
                                   model_coef[row_index, grep("FC", colnames(model_coef))],
                                   expt=expt))
                    })
plot_data <- do.call(rbind, plot_data)
plot_data$CL <- paste("circligase", plot_data$CL)
levels(plot_data$expt) <- c("Green", "Lareau", "Weinberg")

corr_text <- unique(plot_data[, c("CL", "expt")])
corr_text$corr <- sapply(seq(nrow(corr_text)),
                         function(x) {
                           paste("cor =",
                                 round(with(subset(plot_data, CL==corr_text$CL[x] & expt==corr_text$expt[x]),
                                            cor(ratio, FC)),2))
                         })
corr_text$pvalue <- sapply(seq(nrow(corr_text)),
                           function(x) {
                             paste("p =",
                                   signif(with(subset(plot_data, CL==corr_text$CL[x] & expt==corr_text$expt[x]),
                                               cor.test(ratio, FC)$p.value),2))
                           })
corr_text$x <- min(plot_data$ratio) + 0.2*(max(plot_data$ratio)-min(plot_data$ratio))
corr_text$y_corr <- 1.3*max(plot_data$FC)
corr_text$y_pvalue <- 1.2*max(plot_data$FC)

ggplot(plot_data, aes(x=ratio, y=FC)) + geom_point() +
  facet_grid(CL ~ expt) +
  geom_errorbar(aes(ymin=FC_lower, ymax=FC_upper)) +
  geom_errorbarh(aes(xmin=ratio-se_abs, xmax=ratio+se_abs)) +
  geom_smooth(method="lm", formula="y ~ x") +
  geom_text(data=corr_text, mapping=aes(x=x, y=y_corr, label=corr)) +
  geom_text(data=corr_text, mapping=aes(x=x, y=y_pvalue, label=pvalue)) +
  theme_classic() + theme(panel.background = element_rect(fill = NA, color = "black")) +
  xlab("in vitro ligation efficiency") + ylab(expression("exp("*beta*")"))
