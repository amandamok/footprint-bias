rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "WBcel235.transcripts.fa")
transcript_length_fname <- file.path(ref_dir, "WBcel235.transcripts.lengths.txt")
transcript_lengths <- load_lengths(transcript_length_fname)
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 250

regression_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

results_dir <- file.path(here(), "expts", "aeschimann_2015", "GCop")

# read in footprint alignments --------------------------------------------

raw_bam_fname <- file.path(here(), "expts", "aeschimann_2015", "raw_data", "GCop",
                           "GCop_trim3_footprints.transcript.bam")

bam_fname <- file.path(results_dir, "GCop_bam.Rda")
if(!file.exists(bam_fname)) {
  GCop_bam <- load_bam(raw_bam_fname, transcript_fa_fname, transcript_length_fname,
                        offsets_fname, f5_length=f5_length, f3_length=f3_length)
  save(GCop_bam, file=bam_fname)
} else {
  load(bam_fname)
}

# compute size/frame subsets ----------------------------------------------

subsets_fname <- file.path(results_dir, "GCop_subsets.Rda")
if(!file.exists(subsets_fname)) {
  GCop_d5_d3 <- count_d5_d3(GCop_bam)
  GCop_subsets <- GCop_d5_d3$counts[1:(which(GCop_d5_d3$counts$proportion>min_prop)[1]),
                                      c("d5", "d3")]
  GCop_subset_names <- sapply(seq(nrow(GCop_subsets)),
                               function(x) {
                                 paste("d5", GCop_subsets$d5[x], "d3", GCop_subsets$d3[x], sep="_")
                               })
  save(GCop_d5_d3, GCop_subsets, GCop_subset_names, file=subsets_fname)
} else {
  load(subsets_fname)
}

# establish training data -------------------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=GCop_bam, FUN=sum)
transcript_counts <- transcript_counts[order(transcript_counts$count, decreasing=T),]

training_set <- as.character(transcript_counts$transcript[1:num_genes])

# initialize training data for regression ---------------------------------

training_fname <- file.path(results_dir, "GCop_training.Rda")
if(!file.exists(training_fname)) {
  GCop_training <- init_data(transcript_fa_fname, transcript_length_fname,
                              d5_d3_subsets=GCop_subsets, f5_length=f5_length, f3_length=f3_length,
                              which_transcripts=training_set)
  GCop_training$transcript <- relevel(GCop_training$transcript, ref=training_set[1])
  GCop_training$count <- count_footprints(GCop_bam, GCop_training, "count")
  save(GCop_training, file=training_fname)
} else {
  load(training_fname)
}

# compute regression ------------------------------------------------------

fit_fname <- file.path(results_dir, "GCop_fit_250.Rda")
if(!file.exists(fit_fname)) {
  GCop_fit_250 <- MASS::glm.nb(regression_model, data=GCop_training, model=F)
  GCop_coef_250 <- parse_coefs(GCop_fit_250)
  save(GCop_fit_250, file=fit_fname)
  save(GCop_coef_250, file=file.path(results_dir, "GCop_coef_250.Rda"))
} else {
  load(fit_fname)
}

# correct counts ----------------------------------------------------------

GCop_bam$correct_250 <- correct_bias(GCop_bam, GCop_fit_250)
GCop_training$correct_250 <- count_footprints(GCop_bam, GCop_training, "correct_250")

save(GCop_bam, file=bam_fname)
save(GCop_training, file=training_fname)

# evaluate bias -----------------------------------------------------------

to_evaluate <- c("count", "correct_250")
GCop_codon_corr <- lapply(to_evaluate,
                           function(x) {
                             evaluate_bias(GCop_training, which_column=x,
                                           transcript_fa_fname, transcript_length_fname,
                                           type="codon")
                           })
names(GCop_codon_corr) <- to_evaluate
GCop_nt_corr <- lapply(to_evaluate,
                        function(x) {
                          evaluate_bias(GCop_training, which_column=x,
                                        transcript_fa_fname, transcript_length_fname,
                                        type="nt")
                        })
names(GCop_nt_corr) <- to_evaluate

codon_max <- max(unlist(GCop_codon_corr))
nt_max <- max(unlist(GCop_nt_corr))

GCop_plot <- (plot_bias(GCop_codon_corr$count) +
                 coord_cartesian(ylim=c(0, codon_max)) +
                 ggtitle("WT", subtitle="raw")) +
  (plot_bias(GCop_codon_corr$correct_250) +
     coord_cartesian(ylim=c(0, codon_max)) +
     ggtitle("", subtitle="corrected")) +
  (plot_bias(GCop_nt_corr$count, type="nt") + coord_cartesian(ylim=c(0, nt_max))) +
  (plot_bias(GCop_nt_corr$correct_250, type="nt") + coord_cartesian(ylim=c(0, nt_max)))

save(GCop_codon_corr, GCop_nt_corr, GCop_plot, file="GCop_corr.Rda")

