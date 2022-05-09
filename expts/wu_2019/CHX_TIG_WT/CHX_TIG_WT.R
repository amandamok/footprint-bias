#####
# Wu 2019: WT CHX TIG

rm(list=ls())

library(ggplot2)
library(patchwork)
library(here)

scripts_dir <- file.path(here(), "scripts")
source(file.path(scripts_dir, "helper.R"))
source(file.path(scripts_dir, "prep_data.R"))
source(file.path(scripts_dir, "correct_bias.R"))
source(file.path(scripts_dir, "evaluate_bias.R"))

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules_rotated.txt")

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- 150

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)
utr5_length <- unique(transcript_lengths$utr5_length)
utr3_length <- unique(transcript_lengths$utr3_length)

sra_ids <- c("SRR7241919", "SRR7241920") # omit SRR8093858
raw_data_dir <- file.path(here(), "expts", "wu_2019", "raw_data")

results_dir <- file.path(here(), "expts", "wu_2019", "CHX_TIG_WT")

# 1. read in footprint alignments -----------------------------------------

bam_fname <- file.path(results_dir, "CHX_TIG_WT.bam")
if(!file.exists(bam_fname)) {
  system(paste("samtools merge", bam_fname,
               file.path(raw_data_dir, sra_ids, paste0(sra_ids, ".transcript.bam"))))
}

bam_dat_fname <- file.path(results_dir, "CHX_TIG_WT_bam.Rda")
if(!file.exists(bam_dat_fname)) {
  CHX_TIG_WT_bam <- load_bam(bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                             f5_length=f5_length, f3_length=f3_length)
  save(CHX_TIG_WT_bam, file=bam_dat_fname)
} else {
  load(bam_dat_fname)
}

# 2. compute size/frame subsets -------------------------------------------

subsets_fname <- file.path(results_dir, "CHX_TIG_WT_subsets.Rda")
if(!file.exists(subsets_fname)) {
  d5_d3 <- count_d5_d3(CHX_TIG_WT_bam)
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

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=CHX_TIG_WT_bam, FUN=sum)
transcript_counts <- transcript_counts[order(transcript_counts$count, decreasing=T),]

training_set <- as.character(transcript_counts$transcript[1:num_genes])

# 4. initialize data.frame for regression ---------------------------------

training_fname <- file.path(results_dir, "CHX_TIG_WT_training.Rda")
if(!file.exists(training_fname)) {
  CHX_TIG_WT_training <- init_data(transcript_fa_fname, transcript_length_fname,
                                   d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                                   which_transcripts=training_set)
  CHX_TIG_WT_training$transcript <- relevel(CHX_TIG_WT_training$transcript, ref=training_set[1])
  CHX_TIG_WT_training$count <- count_footprints(CHX_TIG_WT_bam, CHX_TIG_WT_training, "count")
  save(CHX_TIG_WT_training, file=training_fname)
} else {
  load(training_fname)
}

# 5. compute regression ---------------------------------------------------

fit_fname <- file.path(results_dir, "CHX_TIG_WT_fit_150.Rda")
if(!file.exists(fit_fname)) {
  CHX_TIG_WT_fit_150 <- MASS::glm.nb(interxn_model, data=CHX_TIG_WT_training, model=F)
  save(CHX_TIG_WT_fit_150, file=fit_fname)
} else {
  load(fit_fname)
}

