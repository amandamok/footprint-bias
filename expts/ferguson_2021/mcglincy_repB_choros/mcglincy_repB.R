rm(list=ls())

library(choros)
library(here)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

mcglincy_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "mcglincy_repB",
                                "mcglincy_trimmed_deduplicated_footprints.transcript.bam")

mcglincy_repB_v1 <- choros(mcglincy_bam_fname,
                           transcript_fa_fname,
                           transcript_length_fname,
                           offsets_fname)

save(mcglincy_repB_v1, file="mcglincy_repB_v1.Rda")

training_set <- readLines("training_transcripts.txt")

## v2: n = 150 transcripts for training
regression_fit_150 <- MASS::glm.nb(formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3),
                                   data=subset(training_data, transcript %in% training_set[1:150]),
                                   model=F)
regression_coefs_150 <- parse_coefs(regression_fit_150)
alignment_data$correct_150 <- correct_bias(alignment_data, regression_fit_150)
training_data$correct_150 <- count_footprints(alignment_data, training_data, "correct_150")
codon_corr <- lapply(c("count", "corrected_count", "correct_150"),
                     function(x) {
                       evaluate_bias(training_data, which_column=x,
                                     transcript_fa_fname, transcript_length_fname,
                                     type="codon")
                     })
codon_corr_plots <- lapply(codon_corr, plot_bias)
nt_corr <- lapply(c("count", "corrected_count", "correct_150"),
                  function(x) {
                    evaluate_bias(training_data, which_column=x,
                                  transcript_fa_fname, transcript_length_fname,
                                  type="nt")
                  })
nt_corr_plots <- lapply(nt_corr, plot_bias, type="nt")
mcglincy_repB_v2 <- list(bam_alignment = alignment_data,
                         regression_data = training_data,
                         model_fit = regression_fit_150,
                         regression_coefs = regression_coefs_150,
                         codon_bias = codon_corr,
                         nt_bias = nt_corr,
                         codon_plots = codon_corr_plots,
                         nt_plots = nt_corr_plots)
save(mcglincy_repB_v2, file="mcglincy_repB_v2.Rda")

## v3: n = 50 transcripts for training
regression_fit_50 <- MASS::glm.nb(formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3),
                                   data=subset(training_data, transcript %in% training_set[1:50]),
                                   model=F)
regression_coefs_50 <- parse_coefs(regression_fit_50)
alignment_data$correct_50 <- correct_bias(alignment_data, regression_fit_50)
training_data$correct_50 <- count_footprints(alignment_data, training_data, "correct_50")
codon_corr <- lapply(c("count", "corrected_count", "correct_50"),
                     function(x) {
                       evaluate_bias(training_data, which_column=x,
                                     transcript_fa_fname, transcript_length_fname,
                                     type="codon")
                     })
codon_corr_plots <- lapply(codon_corr, plot_bias)
nt_corr <- lapply(c("count", "corrected_count", "correct_50"),
                  function(x) {
                    evaluate_bias(training_data, which_column=x,
                                  transcript_fa_fname, transcript_length_fname,
                                  type="nt")
                  })
nt_corr_plots <- lapply(nt_corr, plot_bias, type="nt")
mcglincy_repB_v3 <- list(bam_alignment = alignment_data,
                         regression_data = training_data,
                         model_fit = regression_fit_50,
                         regression_coefs = regression_coefs_50,
                         codon_bias = codon_corr,
                         nt_bias = nt_corr,
                         codon_plots = codon_corr_plots,
                         nt_plots = nt_corr_plots)
save(mcglincy_repB_v3, file="mcglincy_repB_v3.Rda")


q(save="no")
