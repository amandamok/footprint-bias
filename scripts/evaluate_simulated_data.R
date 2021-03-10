rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

## green expt: bias scores for ligate() & circularize()
codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")
# /mnt/lareaulab/rtunney/iXnos/results/green/s28_cod_n5p4_nt_n15p14/epoch30/codon_scores.tsv
bias_seq_2nt <- unique(substr(codons, start=1, stop=2))
green_biasScores <- read.table(file.path("~/simRiboSeq", "refData", "codon_scores.tsv"))
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

expt_dir <- file.path(here(), "expts", "simulated_data")
expts <- c("noBias", "n3Bias", "p5Bias", "bothBias")
plot_data <- lapply(expts,
                    function(expt) {
                      # load model fit
                      model_obj <- paste0(expt, "_fit_200")
                      load(file.path(expts_dir, expt, paste0(model_obj, ".Rda")))
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
                      # return model parameters
                      model_coef$expt <- expt
                      return(model_coef)
                    })