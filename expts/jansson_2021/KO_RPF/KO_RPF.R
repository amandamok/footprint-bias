rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "grch38.transcripts.fa")
transcript_length_fname <- file.path(ref_dir, "grch38.transcripts.lengths.tsv")
transcript_lengths <- load_lengths(transcript_length_fname)
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 250

regression_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

results_dir <- file.path(here(), "expts", "jansson_2021", "KO_RPF")

# read in footprint alignments --------------------------------------------

raw_bam_fname <- file.path(here(), "expts", "jansson_2021", "raw_data", "KO_RPF",
                           "KO_RPF_trim3_footprints.transcript.bam")

bam_fname <- file.path(results_dir, "KO_bam.Rda")
if(!file.exists(bam_fname)) {
  KO_bam <- load_bam(raw_bam_fname, transcript_fa_fname, transcript_length_fname,
                     offsets_fname, f5_length=f5_length, f3_length=f3_length)
  save(KO_bam, file=bam_fname)
} else {
  load(bam_fname)
}

# compute size/frame subsets ----------------------------------------------

subsets_fname <- file.path(results_dir, "KO_subsets.Rda")
if(!file.exists(subsets_fname)) {
  KO_d5_d3 <- count_d5_d3(KO_bam)
  KO_subsets <- KO_d5_d3$counts[1:(which(KO_d5_d3$counts$proportion>min_prop)[1]),
                                c("d5", "d3")]
  KO_subset_names <- sapply(seq(nrow(KO_subsets)),
                            function(x) {
                              paste("d5", KO_subsets$d5[x], "d3", KO_subsets$d3[x], sep="_")
                            })
  save(KO_d5_d3, KO_subsets, KO_subset_names, file=subsets_fname)
} else {
  load(subsets_fname)
}

# establish training data -------------------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=KO_bam, FUN=sum)
transcript_counts <- transcript_counts[order(transcript_counts$count, decreasing=T),]

training_set <- as.character(transcript_counts$transcript[1:num_genes])

# initialize training data for regression ---------------------------------

training_fname <- file.path(results_dir, "KO_training.Rda")
if(!file.exists(training_fname)) {
  KO_training <- init_data(transcript_fa_fname, transcript_length_fname,
                           d5_d3_subsets=KO_subsets, f5_length=f5_length, f3_length=f3_length,
                           which_transcripts=training_set)
  KO_training$transcript <- relevel(KO_training$transcript, ref=training_set[1])
  KO_training$count <- count_footprints(KO_bam, KO_training, "count")
  save(KO_training, file=training_fname)
} else {
  load(training_fname)
}

# compute regression ------------------------------------------------------

fit_fname <- file.path(results_dir, "KO_fit_250.Rda")
if(!file.exists(fit_fname)) {
  KO_fit_250 <- MASS::glm.nb(regression_model, data=KO_training, model=F)
  KO_coef_250 <- parse_coefs(KO_fit_250)
  save(KO_fit_250, file=fit_fname)
  save(KO_coef_250, file=file.path(results_dir, "KO_coef_250.Rda"))
} else {
  load(fit_fname)
}

# correct counts ----------------------------------------------------------

KO_bam$correct_250 <- correct_bias(KO_bam, KO_fit_250)
KO_training$correct_250 <- count_footprints(KO_bam, KO_training, "correct_250")

save(KO_bam, bam_fname)
save(KO_training, training_fname)

# evaluate bias -----------------------------------------------------------

to_evaluate <- c("count", "correct_250")
KO_codon_corr <- lapply(to_evaluate,
                        function(x) {
                          evaluate_bias(KO_training, which_column=x,
                                        transcript_fa_fname, transcript_length_fname,
                                        type="codon")
                        })
names(KO_codon_corr) <- to_evaluate
KO_nt_corr <- lapply(to_evaluate,
                     function(x) {
                       evaluate_bias(KO_training, which_column=x,
                                     transcript_fa_fname, transcript_length_fname,
                                     type="nt")
                     })
names(KO_nt_corr) <- to_evaluate

codon_max <- max(unlist(KO_codon_corr))
nt_max <- max(unlist(KO_nt_corr))

KO_plot <- (plot_bias(KO_codon_corr$count) +
              coord_cartesian(ylim=c(0, codon_max)) +
              ggtitle("KO", subtitle="raw")) +
  (plot_bias(KO_codon_corr$correct_250) +
     coord_cartesian(ylim=c(0, codon_max)) +
     ggtitle("", subtitle="corrected")) +
  (plot_bias(KO_nt_corr$count, type="nt") + coord_cartesian(ylim=c(0, nt_max))) +
  (plot_bias(KO_nt_corr$correct_250, type="nt") + coord_cartesian(ylim=c(0, nt_max)))

save(KO_codon_corr, KO_nt_corr, KO_plot, file="KO_corr.Rda")
