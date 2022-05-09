rm(list=ls())

library(here)
library(ggplot2)

expt_dir <- file.path(here(), "expts", "cheng_2018")
comparisons_dir <- file.path(expt_dir, "timecourse_comparison")

codon_properties <- read.table(file.path(here(), "reference_data", 
                                         "yeast_codon_properties.txt"),
                               header=T)

timepoints <- c(paste0(c(1.5, 3, 4.5, 6, 8, 10), "hr"), "spore")

# load all regression coefficients
all_coef <- lapply(timepoints,
                   function(x) {
                     tmp_comparison <- paste0("timepoint_0hr_timepoint_", x)
                     load(file.path(comparisons_dir, tmp_comparison,
                                    paste0(tmp_comparison, "_coef.Rda")))
                     within(interxn_coef, 
                            timepoint <- x)
                   })
all_coef <- do.call(rbind, all_coef)

# plot A site interaction term 
coef_A <- subset(all_coef, grepl("^timepoint.*:A", all_coef$term))
coef_A$codon <- sub("timepoint.*:A", "", coef_A$term)
coef_A$aa <- codon_properties$aa[match(coef_A$codon, codon_properties$codon)]
coef_A$timepoint[coef_A$timepoint=="spore"] <- "24hr"
coef_A$timepoint <- as.numeric(sub("hr", "", coef_A$timepoint))
coef_A$pvalue <- ifelse(p.adjust(coef_A$Pr...z.., method="fdr") < 0.05,
                        "p_adj < 0.05", "p_adj >= 0.05")
## timecourse
ggplot(coef_A, aes(x=timepoint, y=Estimate,
                   ymin=Estimate + qnorm(0.025)*Std..Error,
                   ymax=Estimate + qnorm(0.975)*Std..Error,
                   group=codon, col=pvalue)) + 
  geom_point() + geom_line(alpha=0.5) + geom_errorbar(alpha=0.5) + 
  theme_bw() + facet_wrap(~aa, scales="free_y") + labs(col="") + 
  scale_color_manual(values=c("p_adj < 0.05"="red", "p_adj >= 0.05"="black")) + 
  ggtitle("A site ribosome occupancy", subtitle="relative to 0hr timepoint") + 
  xlab("time (hr)") + ylab("logFC")
## volcano plots
ggplot(coef_A, aes(x=Estimate, y=-log10(Pr...z..), col=codon)) + 
  geom_point() + theme_bw() + guides(col="none") + 
  facet_wrap(~timepoint, scales="free_y") + 
  xlab(expr(beta)) + ylab("-log10(p)") + ggtitle("A site timecourse")

# plot P site interaction term
coef_P <- subset(all_coef, grepl("^timepoint.*:P", all_coef$term))
coef_P$codon <- sub("timepoint.*:P", "", coef_P$term)
coef_P$aa <- codon_properties$aa[match(coef_P$codon, codon_properties$codon)]
coef_P$timepoint[coef_P$timepoint=="spore"] <- "24hr"
coef_P$timepoint <- as.numeric(sub("hr", "", coef_P$timepoint))
coef_P$pvalue <- ifelse(p.adjust(coef_P$Pr...z.., method="fdr") < 0.05,
                        "p_adj < 0.05", "p_adj >= 0.05")
ggplot(coef_P, aes(x=timepoint, y=Estimate,
                   ymin=Estimate + qnorm(0.025)*Std..Error,
                   ymax=Estimate + qnorm(0.975)*Std..Error,
                   group=codon, col=pvalue)) + 
  geom_point() + geom_line(alpha=0.5) + geom_errorbar(alpha=0.5) + 
  theme_bw() + facet_wrap(~aa, scales="free_y") + labs(col="") + 
  scale_color_manual(values=c("p_adj < 0.05"="red", "p_adj >= 0.05"="black")) + 
  ggtitle("P site ribosome occupancy", subtitle="relative to 0hr timepoint") + 
  xlab("time (hr)") + ylab("logFC")
## volcano plots
ggplot(coef_P, aes(x=Estimate, y=-log10(Pr...z..), col=codon)) + 
  geom_point() + theme_bw() + guides(col="none") + 
  facet_wrap(~timepoint, scales="free_y") + 
  xlab(expr(beta)) + ylab("-log10(p)") + ggtitle("P site timecourse")

