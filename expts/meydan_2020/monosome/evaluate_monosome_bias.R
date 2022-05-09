rm(list=ls())

library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")

paper <- "meydan_2020"
expt <- "monosome"
results_dir <- file.path(here(), "expts", paper, expt)
num_genes <- 250

training_obj <- paste0(expt, "_training")
load(file.path(results_dir, paste0(training_obj, ".Rda")))

to_evaluate <- c("count",
                 paste0("correct_", num_genes, "_gc"))

codon_corr_obj <- paste0(expt, "_codon_corr")
assign(codon_corr_obj,
       lapply(to_evaluate,
              function(x) {
                evaluate_bias(get(training_obj), which_column=x,
                              transcript_fa_fname, transcript_length_fname,
                              type="codon")
              }))

nt_corr_obj <- paste0(expt, "_nt_corr")
assign(nt_corr_obj,
       lapply(to_evaluate,
              function(x) {
                evaluate_bias(get(training_obj), which_column=x,
                              transcript_fa_fname, transcript_length_fname,
                              type="nt")
              }))

save(list=c(codon_corr_obj, nt_corr_obj),
     file=file.path(results_dir, paste0(expt, "_corr.Rda")))

codon_corr_max <- max(unlist(get(codon_corr_obj)))
nt_corr_max <- max(unlist(get(nt_corr_obj)))

corr_plots <- (plot_bias(get(codon_corr_obj)[[1]]) +
                 coord_cartesian(ylim=c(0, codon_corr_max)) +
                 ggtitle(expt, subtitle="raw counts")) +
  (plot_bias(get(codon_corr_obj)[[2]]) +
     coord_cartesian(ylim=c(0, codon_corr_max)) +
     ggtitle("", subtitle="correct end nt + %gc")) +
  (plot_bias(get(nt_corr_obj)[[1]], type="nt") +
     coord_cartesian(ylim=c(0, nt_corr_max))) +
  (plot_bias(get(nt_corr_obj)[[2]], type="nt") +
     coord_cartesian(ylim=c(0, nt_corr_max))) +

  save(corr_plots, file=file.path(results_dir, "corr_plots.Rda"))

q(save="no")
