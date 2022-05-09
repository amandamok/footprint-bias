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

results_dir <- file.path(here(), "expts", "aeschimann_2015", "GCkit")

# read in footprint alignments --------------------------------------------

raw_bam_fname <- file.path(here(), "expts", "aeschimann_2015", "raw_data", "GCkit",
                           "GCkit_trim3_footprints.transcript.bam")

bam_fname <- file.path(results_dir, "GCkit_bam.Rda")
if(!file.exists(bam_fname)) {
  GCkit_bam <- load_bam(raw_bam_fname, transcript_fa_fname, transcript_length_fname,
                        offsets_fname, f5_length=f5_length, f3_length=f3_length)
  save(GCkit_bam, file=bam_fname)
} else {
  load(bam_fname)
}

# compute size/frame subsets ----------------------------------------------

subsets_fname <- file.path(results_dir, "GCkit_subsets.Rda")
if(!file.exists(subsets_fname)) {
  GCkit_d5_d3 <- count_d5_d3(GCkit_bam)
  GCkit_subsets <- GCkit_d5_d3$counts[1:(which(GCkit_d5_d3$counts$proportion>min_prop)[1]),
                                      c("d5", "d3")]
  GCkit_subset_names <- sapply(seq(nrow(GCkit_subsets)),
                               function(x) {
                                 paste("d5", GCkit_subsets$d5[x], "d3", GCkit_subsets$d3[x], sep="_")
                               })
  save(GCkit_d5_d3, GCkit_subsets, GCkit_subset_names, file=subsets_fname)
} else {
  load(subsets_fname)
}

# establish training data -------------------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=GCkit_bam, FUN=sum)
transcript_counts <- transcript_counts[order(transcript_counts$count, decreasing=T),]

training_set <- as.character(transcript_counts$transcript[1:num_genes])

# initialize training data for regression ---------------------------------

training_fname <- file.path(results_dir, "GCkit_training.Rda")
if(!file.exists(training_fname)) {
  GCkit_training <- init_data(transcript_fa_fname, transcript_length_fname,
                              d5_d3_subsets=GCkit_subsets, f5_length=f5_length, f3_length=f3_length,
                              which_transcripts=training_set)
  GCkit_training$transcript <- relevel(GCkit_training$transcript, ref=training_set[1])
  GCkit_training$count <- count_footprints(GCkit_bam, GCkit_training, "count")
  save(GCkit_training, file=training_fname)
} else {
  load(training_fname)
}

# compute regression ------------------------------------------------------

fit_fname <- file.path(results_dir, "GCkit_fit_250.Rda")
if(!file.exists(fit_fname)) {
  GCkit_fit_250 <- MASS::glm.nb(regression_model, data=GCkit_training, model=F)
  GCkit_coef_250 <- parse_coefs(GCkit_fit_250)
  save(GCkit_fit_250, file=fit_fname)
  save(GCkit_coef_250, file=file.path(results_dir, "GCkit_coef_250.Rda"))
} else {
  load(fit_fname)
}

# correct counts ----------------------------------------------------------

GCkit_bam$correct_250 <- correct_bias(GCkit_bam, GCkit_fit_250)
GCkit_training$correct_250 <- count_footprints(GCkit_bam, GCkit_training, "correct_250")

save(GCkit_bam, file=bam_fname)
save(GCkit_training, file=training_fname)

# evaluate bias -----------------------------------------------------------

to_evaluate <- c("count", "correct_250")
GCkit_codon_corr <- lapply(to_evaluate,
                           function(x) {
                             evaluate_bias(GCkit_training, which_column=x,
                                           transcript_fa_fname, transcript_length_fname,
                                           type="codon")
                           })
names(GCkit_codon_corr) <- to_evaluate
GCkit_nt_corr <- lapply(to_evaluate,
                        function(x) {
                          evaluate_bias(GCkit_training, which_column=x,
                                        transcript_fa_fname, transcript_length_fname,
                                        type="nt")
                        })
names(GCkit_nt_corr) <- to_evaluate

codon_max <- max(unlist(GCkit_codon_corr))
nt_max <- max(unlist(GCkit_nt_corr))

GCkit_plot <- (plot_bias(GCkit_codon_corr$count) +
                 coord_cartesian(ylim=c(0, codon_max)) +
                 ggtitle("WT", subtitle="raw")) +
  (plot_bias(GCkit_codon_corr$correct_250) +
     coord_cartesian(ylim=c(0, codon_max)) +
     ggtitle("", subtitle="corrected")) +
  (plot_bias(GCkit_nt_corr$count, type="nt") + coord_cartesian(ylim=c(0, nt_max))) +
  (plot_bias(GCkit_nt_corr$correct_250, type="nt") + coord_cartesian(ylim=c(0, nt_max)))

save(GCkit_codon_corr, GCkit_nt_corr, GCkit_plot, file="GCkit_corr.Rda")

