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
                     model_coef$f5 <- sub("genome_f5", "", rownames(model_coef))
                     return(model_coef)
                   })
names(f5_coefs) <- expts
f3_coefs <- lapply(expts,
                   function(expt) {
                     load(file.path(expt_dir, expt, "model_fit_100genes.Rda"))
                     model_coef <- data.frame(summary(model_fit_100)$coefficients)
                     model_coef <- subset(model_coef,
                                          grepl("^genome_f3", rownames(model_coef)))
                     model_coef$f3 <- sub("genome_f3", "", rownames(model_coef))
                     return(model_coef)
                   })
names(f3_coefs) <- expts

# add simulated values ----------------------------------------------------

no_p5bias <- no_p5bias / no_p5bias["AA"]
green_p5bias_2nt <- green_p5bias_2nt / green_p5bias_2nt["AA"]
f5_coefs$noBias$sim <- no_p5bias[f5_coefs$noBias$f5]
f5_coefs$n3Bias$sim <- no_p5bias[f5_coefs$n3Bias$f5]
f5_coefs$p5Bias$sim <- green_p5bias_2nt[f5_coefs$p5Bias$f5]
f5_coefs$bothBias$sim <- green_p5bias_2nt[f5_coefs$bothBias$f5]

no_n3bias <- no_n3bias / no_n3bias["AAA"]
green_n3bias_3nt <- green_n3bias_3nt / green_n3bias_3nt["AAA"]
f3_coefs$noBias$sim <- no_n3bias[f3_coefs$noBias$f3]
f3_coefs$n3Bias$sim <- green_n3bias_3nt[f3_coefs$n3Bias$f3]
f3_coefs$p5Bias$sim <- no_n3bias[f3_coefs$p5Bias$f3]
f3_coefs$bothBias$sim <- green_n3bias_3nt[f3_coefs$bothBias$f3]

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
