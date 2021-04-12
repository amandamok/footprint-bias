rm(list=ls())


# load simulation values --------------------------------------------------

simriboseq_dir <- "~/simRiboSeq"
ref_dir <- file.path(simriboseq_dir, "refData")

codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")

## green expt: bias scores for ligate() & circularize()
# /mnt/lareaulab/rtunney/iXnos/results/green/s28_cod_n5p4_nt_n15p14/epoch30/codon_scores.tsv
bias_seq_2nt <- unique(substr(codons, start=1, stop=2))
green_biasFile <- "codon_scores.tsv"
green_biasScores <- read.table(file.path(ref_dir, green_biasFile))
colnames(green_biasScores) <- as.character(seq.int(from=-5, length.out=ncol(green_biasScores)))
# 5' bias probabilities
green_p5bias <- exp(green_biasScores[,"-5"])
green_p5bias <- green_p5bias/max(green_p5bias, na.rm=T)
green_p5bias[is.na(green_p5bias)] <- 0
names(green_p5bias) <- sort(codons)
green_p5bias_2nt <- sapply(bias_seq_2nt,
                           function(x) {
                             mean(green_p5bias[substr(names(green_p5bias), start=1, stop=2)==x])
                           })
green_p5bias_2nt[green_p5bias_2nt==0] <- 0.01
# 3' bias probabilities
green_n3bias <- exp(green_biasScores[,"3"])
green_n3bias <- green_n3bias/max(green_n3bias, na.rm=T)
green_n3bias[is.na(green_n3bias)] <- 0
names(green_n3bias) <- sort(codons)
green_n3bias_3nt <- green_n3bias
green_n3bias_3nt[green_n3bias_3nt==0] <- 0.01
# no bias (uniform bias)
no_p5bias <- rep(1, length(green_p5bias_2nt))
names(no_p5bias) <- names(green_p5bias_2nt)
no_n3bias <- rep(1, length(codons))
names(no_n3bias) <- codons

# load regression coefficients --------------------------------------------

expt_dir <- file.path(here(), "expts", "simulated_data")
expts <- c("noBias", "n3Bias", "p5Bias", "bothBias")
f5_coefs <- lapply(expts,
                   function(expt) {
                     load(file.path(expt_dir, expt, "model_fit_100genes.Rda"))
                     model_coef <- data.frame(summary(model_fit_100)$coefficients)
                     model_coef <- subset(model_coef,
                                          grepl("^genome_f5", rownames(model_coef)))
                     model_coef$seq <- sub("genome_f5", "", rownames(model_coef))
                     model_coef$expt <- expt
                     return(model_coef)
                   })
names(f5_coefs) <- expts
f3_coefs <- lapply(expts,
                   function(expt) {
                     load(file.path(expt_dir, expt, "model_fit_100genes.Rda"))
                     model_coef <- data.frame(summary(model_fit_100)$coefficients)
                     model_coef <- subset(model_coef,
                                          grepl("^genome_f3", rownames(model_coef)))
                     model_coef$seq <- sub("genome_f3", "", rownames(model_coef))
                     model_coef$expt <- expt
                     return(model_coef)
                   })
names(f3_coefs) <- expts

# add simulated values ----------------------------------------------------

no_p5bias <- no_p5bias / no_p5bias["AA"]
green_p5bias_2nt <- green_p5bias_2nt / green_p5bias_2nt["AA"]
f5_coefs$noBias$sim <- no_p5bias[f5_coefs$noBias$seq]
f5_coefs$n3Bias$sim <- no_p5bias[f5_coefs$n3Bias$seq]
f5_coefs$p5Bias$sim <- green_p5bias_2nt[f5_coefs$p5Bias$seq]
f5_coefs$bothBias$sim <- green_p5bias_2nt[f5_coefs$bothBias$seq]

no_n3bias <- no_n3bias / no_n3bias["AAA"]
green_n3bias_3nt <- green_n3bias_3nt / green_n3bias_3nt["AAA"]
f3_coefs$noBias$sim <- no_n3bias[f3_coefs$noBias$seq]
f3_coefs$n3Bias$sim <- green_n3bias_3nt[f3_coefs$n3Bias$seq]
f3_coefs$p5Bias$sim <- no_n3bias[f3_coefs$p5Bias$seq]
f3_coefs$bothBias$sim <- green_n3bias_3nt[f3_coefs$bothBias$seq]

# make plots --------------------------------------------------------------

noBias_f5 <- ggplot(f5_coefs$noBias, aes(x=sim, y=exp(Estimate))) +
  geom_point() + theme_bw() + geom_vline(xintercept=1) +
  xlab("simulation FC") + ylab(expr(paste("exp(", beta, ")"))) +
  xlim(0.40, 2.75) + ylim(0.40, 2.75)
noBias_f3 <- ggplot(f3_coefs$noBias, aes(x=sim, y=exp(Estimate))) +
  geom_point() + theme_bw() + geom_vline(xintercept=1) +
  xlab("simulation FC") + ylab(expr(paste("exp(", beta, ")"))) +
  xlim(0.03, 3.75) + ylim(0.03, 3.75)

n3Bias_f5 <- ggplot(f5_coefs$n3Bias, aes(x=sim, y=exp(Estimate))) +
  geom_point() + theme_bw() + geom_vline(xintercept=1) +
  xlab("simulation FC") + ylab(expr(paste("exp(", beta, ")"))) +
  xlim(0.40, 2.75) + ylim(0.40, 2.75)
n3Bias_f3 <- ggplot(f3_coefs$n3Bias, aes(x=sim, y=exp(Estimate))) +
  geom_point() + theme_bw() + geom_abline(intercept=0, slope=1) +
  xlab("simulation FC") + ylab(expr(paste("exp(", beta, ")"))) +
  xlim(0.03, 3.75) + ylim(0.03, 3.75)

p5Bias_f5 <- ggplot(f5_coefs$p5Bias, aes(x=sim, y=exp(Estimate))) +
  geom_point() + theme_bw() + geom_abline(intercept=0, slope=1) +
  xlab("simulation FC") + ylab(expr(paste("exp(", beta, ")"))) +
  xlim(0.40, 2.75) + ylim(0.40, 2.75)
p5Bias_f3 <- ggplot(f3_coefs$p5Bias, aes(x=sim, y=exp(Estimate))) +
  geom_point() + theme_bw() + geom_vline(xintercept=1) +
  xlab("simulation FC") + ylab(expr(paste("exp(", beta, ")"))) +
  xlim(0.03, 3.75) + ylim(0.03, 3.75)

bothBias_f5 <- ggplot(f5_coefs$bothBias, aes(x=sim, y=exp(Estimate))) +
  geom_point() + theme_bw() + geom_abline(intercept=0, slope=1) +
  xlab("simulation FC") + ylab(expr(paste("exp(", beta, ")"))) +
  xlim(0.40, 2.75) + ylim(0.40, 2.75)
bothBias_f3 <- ggplot(f3_coefs$bothBias, aes(x=sim, y=exp(Estimate))) +
  geom_point() + theme_bw() + geom_abline(intercept=0, slope=1) +
  xlab("simulation FC") + ylab(expr(paste("exp(", beta, ")"))) +
  xlim(0.03, 3.75) + ylim(0.03, 3.75)

noBias_f5 + n3Bias_f5 + p5Bias_f5 + bothBias_f5 +
  noBias_f3 + n3Bias_f3 + p5Bias_f3 + bothBias_f3 +
  plot_layout(nrow=2, byrow=T, heights=c(1,1))

ggsave(filename="simulation_f5_f3_coefs.pdf")

# aggregate plot ----------------------------------------------------------

f5_coefs <- do.call(rbind, f5_coefs)
f5_coefs$end_nt <- substr(f5_coefs$seq, 1, 1)
f5_coefs$region <- "5' bias"

f3_coefs <- do.call(rbind, f3_coefs)
f3_coefs$end_nt <- substr(f3_coefs$seq, 3, 3)
f3_coefs$region <- "3' bias"

all_coefs <- rbind(f5_coefs, f3_coefs)
all_coefs$region <- factor(all_coefs$region, levels=c("5' bias", "3' bias"))
all_coefs$expt <- factor(all_coefs$expt, levels=c("noBias", "n3Bias", "p5Bias", "bothBias"))
levels(all_coefs$expt) <- c("no bias", "only 3' bias", "only 5' bias", "5' and 3' bias")

ggplot(all_coefs, aes(x=sim, y=exp(Estimate), col=end_nt)) +
  geom_point(alpha=0.5, size=3) + theme_bw() + facet_grid(expt ~ region) +
  xlab("simulated parameter") + ylab("recovered parameter") +
  theme(legend.position="none")
