rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

ref_dir <- file.path(here(), "expts", "jansson_2021", "reference_data")
transcript_fa_fname <- file.path(ref_dir, "grch38transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "grch38.transcripts.20cds20.lengths.tsv")
transcript_lengths <- load_lengths(transcript_length_fname)
offsets_fname <- file.path(here(), "reference_data", "Asite_rules.txt")

min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 250

regression_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

results_dir <- file.path(here(), "expts", "jansson_2021", "WT_RPF", "new_transcriptome")

# read in footprint alignments --------------------------------------------

raw_bam_fname <- "WT_RPF_trim3_footprints.transcript.bam"

bam_fname <- file.path(results_dir, "WT_bam.Rda")
if(!file.exists(bam_fname)) {
  WT_bam <- load_bam(raw_bam_fname, transcript_fa_fname, transcript_length_fname,
                     offsets_fname, f5_length=f5_length, f3_length=f3_length)
  save(WT_bam, file=bam_fname)
} else {
  load(bam_fname)
}

# compute size/frame subsets ----------------------------------------------

subsets_fname <- file.path(results_dir, "WT_subsets.Rda")
if(!file.exists(subsets_fname)) {
  WT_d5_d3 <- count_d5_d3(WT_bam)
  WT_subsets <- WT_d5_d3$counts[1:(which(WT_d5_d3$counts$proportion>min_prop)[1]),
                                c("d5", "d3")]
  WT_subset_names <- sapply(seq(nrow(WT_subsets)),
                            function(x) {
                              paste("d5", WT_subsets$d5[x], "d3", WT_subsets$d3[x], sep="_")
                            })
  save(WT_d5_d3, WT_subsets, WT_subset_names, file=subsets_fname)
} else {
  load(subsets_fname)
}

# establish training data -------------------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=WT_bam, FUN=sum)
transcript_counts <- transcript_counts[order(transcript_counts$count, decreasing=T),]

training_set <- as.character(transcript_counts$transcript[1:num_genes])

# initialize training data for regression ---------------------------------

training_fname <- file.path(results_dir, "WT_training.Rda")
if(!file.exists(training_fname)) {
  WT_training <- init_data(transcript_fa_fname, transcript_length_fname,
                           d5_d3_subsets=WT_subsets, f5_length=f5_length, f3_length=f3_length,
                           which_transcripts=training_set)
  WT_training$transcript <- relevel(WT_training$transcript, ref=training_set[1])
  WT_training$count <- count_footprints(WT_bam, WT_training, "count")
  save(WT_training, file=training_fname)
} else {
  load(training_fname)
}

# compute regression ------------------------------------------------------

fit_fname <- file.path(results_dir, "WT_fit_250.Rda")
if(!file.exists(fit_fname)) {
  WT_fit_250 <- MASS::glm.nb(regression_model, data=WT_training, model=F)
  WT_coef_250 <- parse_coefs(WT_fit_250)
  save(WT_fit_250, file=fit_fname)
  save(WT_coef_250, file=file.path(results_dir, "WT_coef_250.Rda"))
} else {
  load(fit_fname)
}

# correct counts ----------------------------------------------------------

WT_bam$correct_250 <- correct_bias(WT_bam, WT_fit_250)
WT_training$correct_250 <- count_footprints(WT_bam, WT_training, "correct_250")

save(WT_bam, file=bam_fname)
save(WT_training, file=training_fname)

# evaluate bias -----------------------------------------------------------

to_evaluate <- c("count", "correct_250")
WT_codon_corr <- lapply(to_evaluate,
                        function(x) {
                          evaluatrae_bias(WT_training, which_column=x,
                                        transcript_fa_fname, transcript_length_fname,
                                        type="codon")
                        })
names(WT_codon_corr) <- to_evaluate
WT_nt_corr <- lapply(to_evaluate,
                     function(x) {
                       evaluate_bias(WT_training, which_column=x,
                                     transcript_fa_fname, transcript_length_fname,
                                     type="nt")
                     })
names(WT_nt_corr) <- to_evaluate

codon_max <- max(unlist(WT_codon_corr))
nt_max <- max(unlist(WT_nt_corr))

WT_plot <- (plot_bias(WT_codon_corr$count) +
              coord_cartesian(ylim=c(0, codon_max)) +
              ggtitle("WT", subtitle="raw")) +
  (plot_bias(WT_codon_corr$correct_250) +
     coord_cartesian(ylim=c(0, codon_max)) +
     ggtitle("", subtitle="corrected")) +
  (plot_bias(WT_nt_corr$count, type="nt") + coord_cartesian(ylim=c(0, nt_max))) +
  (plot_bias(WT_nt_corr$correct_250, type="nt") + coord_cartesian(ylim=c(0, nt_max)))

save(WT_codon_corr, WT_nt_corr, WT_plot, file="WT_corr.Rda")
