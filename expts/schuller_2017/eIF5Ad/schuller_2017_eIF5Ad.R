#####
# Schuller 2017: eIF5A depletion (yCW33)

rm(list=ls())

library(ggplot2)
library(patchwork)
library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- 150

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)
utr5_length <- unique(transcript_lengths$utr5_length)
utr3_length <- unique(transcript_lengths$utr3_length)

expt <- "schuller_2017"

sra_ids <- c("SRR5008136", "SRR5008137")
raw_data_dir <- file.path(here(), "expts", expt, "raw_data")

results_dir <- file.path(here(), "expts", expt, "eIF5Ad")

# 1. read in footprint alignments -----------------------------------------

bam_fname <- file.path(results_dir, "eIF5Ad.bam")
if(!file.exists(bam_fname)) {
  system(paste("samtools merge", bam_fname,
               file.path(raw_data_dir, sra_ids, paste0(sra_ids, ".transcript.bam"))))
}

bam_dat_fname <- file.path(results_dir, "eIF5Ad_bam.Rda")
if(!file.exists(bam_dat_fname)) {
  eIF5Ad_bam <- load_bam(bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                         f5_length=f5_length, f3_length=f3_length)
  save(eIF5Ad_bam, file=bam_dat_fname)
} else {
  load(bam_dat_fname)
}

# 2. compute size/frame subsets -------------------------------------------

subsets_fname <- file.path(results_dir, "eIF5Ad_subsets.Rda")
if(!file.exists(subsets_fname)) {
  d5_d3 <- count_d5_d3(eIF5Ad_bam)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file=subsets_fname)
} else {
  load(subsets_fname)
}

# 3. establish training data ----------------------------------------------

# calulate average RPF density by codon, per transcript
transcript_density <- aggregate(count ~ transcript + cod_idx, data=eIF5Ad_bam, FUN=sum)
transcript_density <- aggregate(count ~ transcript, data=transcript_density, FUN=mean)
transcript_density <- transcript_density[order(transcript_density$count, decreasing=T),]

training_set <- as.character(transcript_density$transcript[1:num_genes])

# 4. initialize data.frame for regression ---------------------------------

training_fname <- file.path(results_dir, "eIF5Ad_training.Rda")
if(!file.exists(training_fname)) {
  eIF5Ad_training <- init_data(transcript_fa_fname, transcript_length_fname,
                               d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                               which_transcripts=training_set)
  eIF5Ad_training$transcript <- relevel(eIF5Ad_training$transcript, ref=training_set[1])
  eIF5Ad_training$count <- count_footprints(eIF5Ad_bam, eIF5Ad_training, "count")
  save(eIF5Ad_training, file=training_fname)
} else {
  load(training_fname)
}

# 5. compute regression ---------------------------------------------------

fit_fname <- file.path(results_dir, "eIF5Ad_fit_150.Rda")
if(!file.exists(fit_fname)) {
  eIF5Ad_fit_150 <- MASS::glm.nb(interxn_model, data=eIF5Ad_training, model=F)
  save(eIF5Ad_fit_150, file=fit_fname)
} else {
  load(fit_fname)
}

# 6. correct counts -------------------------------------------------------

if(!("correct_150" %in% colnames(eIF5Ad_bam))) {
  eIF5Ad_bam$correct_150 <- correct_bias(eIF5Ad_bam, eIF5Ad_fit_150)
  save(eIF5Ad_bam, file=bam_dat_fname)
}

if(!("correct_150" %in% colnames(eIF5Ad_training))) {
  eIF5Ad_training$correct_150 <- count_footprints(eIF5Ad_bam, eIF5Ad_training, "correct_150")
  save(eIF5Ad_training, file=training_fname)
}

# 7. evaluate and plot bias -----------------------------------------------

bias_eval_fname <- file.path(results_dir, "eIF5Ad_bias.Rda")
if(!file.exists(bias_eval_fname)) {
  eIF5Ad_codon_corr <- evaluate_bias(eIF5Ad_training, which_column="correct_150",
                                     transcript_fa_fname, transcript_length_fname,
                                     utr5=utr5_length, utr3=utr3_length,
                                     type="codon")
  eIF5Ad_bias_plot <- plot_bias(eIF5Ad_codon_corr, type="codon")
  save(eIF5Ad_codon_corr, eIF5Ad_bias_plot, file=bias_eval_fname)
}

q(save="no")

