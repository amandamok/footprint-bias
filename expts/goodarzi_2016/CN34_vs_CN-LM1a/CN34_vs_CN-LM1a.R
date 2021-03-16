rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "grch38.transcripts.fa")
transcript_length_fname <- file.path(ref_dir, "grch38.transcripts.lengths.tsv")
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")
mapping_fname <- file.path(ref_dir, "")

results_dir <- file.path(here(), "expts", "goodarzi_2016", "CN34_vs_CN-LM1a")
raw_data_dir <- file.path(here(), "expts", "goodarzi_2016", "raw_data")

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- 200

sra_CN34 <- c("SRR3129148", "SRR3129149")
sra_LM1a <- c("SRR3129152", "SRR3129153")

# 1. load footprint alignments --------------------------------------------

# CN34
CN34_bam_fname <- file.path(results_dir, "CN34.bam")
if(!file.exists(CN34_bam_fname)) {
  system(paste("samtools merge", CN34_bam_fname,
               file.path(raw_data_dir, sra_CN34, paste0(sra_CN34, ".transcript.bam"))))
}
CN34_bam_dat_fname <- file.path(results_dir, "CN34_bam.Rda")
if(!file.exists(CN34_bam_dat_fname)) {
  CN34_bam <- load_bam(CN34_bam_fname, transcript_fa_fname,
                       transcript_length_fname, offsets_fname,
                       f5_length=f5_length, f3_length=f3_length)
  save(CN34_bam, file=CN34_bam_dat_fname)
} else {
  load(CN34_bam_dat_fname)
}

# LM1a
LM1a_bam_fname <- file.path(results_dir, "LM1a.bam")
if(!file.exists(LM1a_bam_fname)) {
  system(paste("samtools merge", LM1a_bam_fname,
               file.path(raw_data_dir, sra_LM1a, paste0(sra_LM1a, ".transcript.bam"))))
}
LM1a_bam_dat_fname <- file.path(results_dir, "LM1a_bam.Rda")
if(!file.exists(LM1a_bam_dat_fname)) {
  LM1a_bam <- load_bam(LM1a_bam_fname, transcript_fa_fname,
                       transcript_length_fname, offsets_fname,
                       f5_length=f5_length, f3_length=f3_length)
  save(LM1a_bam, file=LM1a_bam_dat_fname)
} else {
  load(LM1a_bam_dat_fname)
}

# 2. compute size/frame subsets -------------------------------------------

CN34_subsets_fname <- file.path(results_dir, "CN34_subsets.Rda")
if(!file.exists(CN34_subsets_fname)) {
  CN34_d5_d3 <- count_d5_d3(CN34_bam)
  CN34_d5_d3_subsets <- CN34_d5_d3$counts[1:(which(CN34_d5_d3$counts$proportion>min_prop)[1]),
                                          c("d5", "d3")]
  CN34_subset_names <- sapply(seq(nrow(CN34_d5_d3_subsets)),
                              function(x) {
                                paste("d5", CN34_d5_d3_subsets$d5[x],
                                      "d3", CN34_d5_d3_subsets$d3[x],
                                      sep="_")
                              })
  save(CN34_d5_d3, CN34_d5_d3_subsets, CN34_subset_names,
       file=CN34_subsets_fname)
} else {
  load(CN34_subsets_fname)
}

LM1a_subsets_fname <- file.path(results_dir, "LM1a_subsets.Rda")
if(!file.exists(LM1a_subsets_fname)) {
  LM1a_d5_d3 <- count_d5_d3(LM1a_bam)
  LM1a_d5_d3_subsets <- LM1a_d5_d3$counts[1:(which(LM1a_d5_d3$counts$proportion>min_prop)[1]),
                                          c("d5", "d3")]
  LM1a_subset_names <- sapply(seq(nrow(LM1a_d5_d3_subsets)),
                              function(x) {
                                paste("d5", LM1a_d5_d3_subsets$d5[x],
                                      "d3", LM1a_d5_d3_subsets$d3[x],
                                      sep="_")
                              })
  save(LM1a_d5_d3, LM1a_d5_d3_subsets, LM1a_subset_names,
       file=LM1a_subsets_fname)
} else {
  load(LM1a_subsets_fname)
}

# pick transcripts for training -------------------------------------------

