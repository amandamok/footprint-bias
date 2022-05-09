#####
# ottrUMI data

rm(list=ls())

library(patchwork)
library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules_ottr_v2.txt")

ottrUMI_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "WT_UMI_OTTR_repB", "no_dedup",
                               "ottrUMI_trimmed_deduplicated_footprints.transcript.bam")

min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 200

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)

# 1. read in footprint alignments -----------------------------------------

if(!file.exists("ottrUMI_f5_3_no_dedup_bam.Rda")) {
  ottrUMI_f5_3_no_dedup_bam <- load_bam(ottrUMI_bam_fname, transcript_fa_fname,
                               transcript_length_fname, offsets_fname,
                               f5_length=f5_length, f3_length=f3_length)
  save(ottrUMI_f5_3_no_dedup_bam, file="ottrUMI_f5_3_no_dedup_bam.Rda")
} else {
  load("ottrUMI_f5_3_no_dedup_bam.Rda")
}

# 2. compute size/frame subsets -------------------------------------------

if(!file.exists("ottrUMI_f5_3_no_dedup_subsets.Rda")) {
  d5_d3 <- count_d5_d3(ottrUMI_f5_3_no_dedup_bam)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file="ottrUMI_f5_3_no_dedup_subsets.Rda")
} else {
  load("ottrUMI_f5_3_no_dedup_subsets.Rda")
}

# 3. establish training and test sets -------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=ottrUMI_f5_3_no_dedup_bam, FUN=sum)
transcript_counts$length_aa <- transcript_lengths$cds_length[match(transcript_counts$transcript,
                                                                   transcript_lengths$transcript)] / 3
transcript_counts$TE <- with(transcript_counts, count / length_aa)
transcript_counts <- transcript_counts[order(transcript_counts$TE, decreasing=T),]

training_set <- as.character(transcript_counts$transcript)[1:num_genes]

# 4. initialize data frames for regression --------------------------------

if(!file.exists("ottrUMI_f5_3_no_dedup_training_data.Rda")) {
  ottrUMI_f5_3_no_dedup_training <- init_data(transcript_fa_fname, transcript_length_fname,
                                     d5_d3_subsets=d5_d3_subsets,
                                     f5_length=f5_length, f3_length=f3_length,
                                     which_transcripts=training_set)
  ottrUMI_f5_3_no_dedup_training$transcript <- relevel(ottrUMI_f5_3_no_dedup_training$transcript, ref=training_set[1])
  ottrUMI_f5_3_no_dedup_training$count <- count_footprints(ottrUMI_f5_3_no_dedup_bam, ottrUMI_f5_3_no_dedup_training, "count")
  save(ottrUMI_f5_3_no_dedup_training, file="ottrUMI_f5_3_no_dedup_training_data.Rda")
} else {
  load("ottrUMI_f5_3_no_dedup_training_data.Rda")
  training_set <- levels(ottrUMI_f5_3_no_dedup_training$transcript)
  training_set <- training_set[order(transcript_counts$TE[match(training_set,
                                                                transcript_counts$transcript)],
                                     decreasing=T)]
}

# 5. compute regressions --------------------------------------------------

if(!file.exists("ottrUMI_f5_3_no_dedup_fit_200.Rda")) {
  # check if reference level for d5 or d3 are multiples of 3
  if((as.numeric(levels(ottrUMI_f5_3_no_dedup_training$d5)[1]) %% 3) == 0) {
    ottrUMI_f5_3_no_dedup_training$d5 <- relevel(ottrUMI_f5_3_no_dedup_training$d5,
                                        ref=levels(ottrUMI_f5_3_no_dedup_training$d5)[2])
  }
  if((as.numeric(levels(ottrUMI_f5_3_no_dedup_training$d3)[1]) %% 3) == 0) {
    ottrUMI_f5_3_no_dedup_training$d3 <- relevel(ottrUMI_f5_3_no_dedup_training$d3,
                                        ref=levels(ottrUMI_f5_3_no_dedup_training$d3)[2])
  }
  ottrUMI_f5_3_no_dedup_fit_200 <- MASS::glm.nb(interxn_model,
                                       data=subset(ottrUMI_f5_3_no_dedup_training,
                                                   transcript %in% training_set[1:num_genes]))
  if(!("try-error" %in% class(ottrUMI_f5_3_no_dedup_fit_200))) {
    save(ottrUMI_f5_3_no_dedup_fit_200, file="ottrUMI_f5_3_no_dedup_fit_200.Rda")
    # save regression coefficients
    ottrUMI_f5_3_no_dedup_coef_200 <- coef(ottrUMI_f5_3_no_dedup_fit_200)
    save(ottrUMI_f5_3_no_dedup_coef_200, file="ottrUMI_f5_3_no_dedup_coef_200.Rda")
  }
} else {
  load("ottrUMI_f5_3_no_dedup_fit_200.Rda")
}

# 6. correct counts -------------------------------------------------------

ottrUMI_f5_3_no_dedup_bam$correct_200 <- correct_bias(ottrUMI_f5_3_no_dedup_bam, ottrUMI_f5_3_no_dedup_fit_200)
ottrUMI_f5_3_no_dedup_training$correct_200 <- count_footprints(ottrUMI_f5_3_no_dedup_bam,
                                                      ottrUMI_f5_3_no_dedup_training,
                                                      "correct_200")

save(ottrUMI_f5_3_no_dedup_bam, file="ottrUMI_f5_3_no_dedup_bam.Rda")
save(ottrUMI_f5_3_no_dedup_training, file="ottrUMI_f5_3_no_dedup_training_data.Rda")

# 7. evaluate bias --------------------------------------------------------

# training data: codons
ottrUMI_f5_3_no_dedup_training_codon_corr <- lapply(c("count", "correct_200"),
                                           function(x) {
                                             evaluate_bias(ottrUMI_f5_3_no_dedup_training,
                                                           which_column=x,
                                                           transcript_fa_fname,
                                                           transcript_length_fname,
                                                           type="codon")
                                           })

ottrUMI_f5_3_no_dedup_training_codon_corr <- do.call(cbind, ottrUMI_f5_3_no_dedup_training_codon_corr)
colnames(ottrUMI_f5_3_no_dedup_training_codon_corr) <- c("uncorrected", "corrected")
save(ottrUMI_f5_3_no_dedup_training_codon_corr, file="ottrUMI_f5_3_no_dedup_training_codon_corr.Rda")

# training data: nucleotides
ottrUMI_f5_3_no_dedup_training_nt_corr <- lapply(c("count", "correct_200"),
                                        function(x) {
                                          evaluate_bias(ottrUMI_f5_3_no_dedup_training,
                                                        which_column=x,
                                                        transcript_fa_fname,
                                                        transcript_length_fname,
                                                        type="nt")
                                        })
ottrUMI_f5_3_no_dedup_training_nt_corr <- do.call(cbind, ottrUMI_f5_3_no_dedup_training_nt_corr)
colnames(ottrUMI_f5_3_no_dedup_training_nt_corr) <- c("uncorrected", "corrected")
save(ottrUMI_f5_3_no_dedup_training_nt_corr, file="ottrUMI_f5_3_no_dedup_training_nt_corr.Rda")

# 8. make plots -----------------------------------------------------------

# training data
ottrUMI_f5_3_no_dedup_training_codon_plots <- lapply(c("uncorrected", "corrected"),
                                                      function(x) {
                                                        tmp <- ottrUMI_f5_3_no_dedup_training_codon_corr[, x]
                                                        names(tmp) <- rownames(ottrUMI_f5_3_no_dedup_training_codon_corr)
                                                        plot_bias(tmp, plot_subtitle=x, type="codon")
                                                      })
ottrUMI_f5_3_no_dedup_training_nt_plots <- lapply(c("uncorrected", "corrected"),
                                                   function(x) {
                                                     tmp <- ottrUMI_f5_3_no_dedup_training_nt_corr[, x]
                                                     names(tmp) <- rownames(ottrUMI_f5_3_no_dedup_training_nt_corr)
                                                     plot_bias(tmp, plot_subtitle=x, type="nt")
                                                   })
ottrUMI_f5_3_no_dedup_training_codon_plots[[1]] <- ottrUMI_f5_3_no_dedup_training_codon_plots[[1]] + ggtitle("ottrUMI (5' bias: 3nt): training data")
ottrUMI_f5_3_no_dedup_training_plots <- wrap_plots(ottrUMI_f5_3_no_dedup_training_codon_plots, nrow=1) / wrap_plots(ottrUMI_f5_3_no_dedup_training_nt_plots, nrow=1)
save(ottrUMI_f5_3_no_dedup_training_plots, file="ottrUMI_f5_3_no_dedup_training_plots.Rda")



(plot_bias(ottrUMI_f5_3_no_dedup_training_codon_corr[,1]) +
    ggtitle("ottrUMI (no deduplication)", subtitle="raw") + ylim(0, 0.07)) +
  (plot_bias(ottrUMI_f5_3_no_dedup_training_codon_corr[,2]) +
     ggtitle("", subtitle="5' 3nt / 3' 3nt") + ylim(0, 0.07)) +
  (plot_bias(ottrUMI_f5_3_no_dedup_training_nt_corr[,1], type="nt") + ylim(0, 0.04)) +
  (plot_bias(ottrUMI_f5_3_no_dedup_training_nt_corr[,2], type="nt") + ylim(0, 0.04))

q(save="no")
