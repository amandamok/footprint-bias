#####
# ottrUMI data

rm(list=ls())

library(patchwork)
library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

ottrUMI_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "WT_UMI_OTTR_repB",
                               "ottrUMI_trimmed_deduplicated_footprints.transcript.bam")

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- 200

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)

# 1. read in footprint alignments -----------------------------------------

if(!file.exists("ottrUMI_bam.Rda")) {
  ottrUMI_bam <- load_bam(ottrUMI_bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                          f5_length=f5_length, f3_length=f3_length)
  save(ottrUMI_bam, file="ottrUMI_bam.Rda")
} else {
  load("ottrUMI_bam.Rda")
}

# 2. compute size/frame subsets -------------------------------------------

if(!file.exists("ottrUMI_subsets.Rda")) {
  d5_d3 <- count_d5_d3(ottrUMI_bam)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file="ottrUMI_subsets.Rda")
} else {
  load("ottrUMI_subsets.Rda")
}

# 3. establish training and test sets -------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=ottrUMI_bam, FUN=sum)
transcript_counts$length_aa <- transcript_lengths$cds_length[match(transcript_counts$transcript,
                                                                   transcript_lengths$transcript)] / 3
transcript_counts$TE <- with(transcript_counts, count / length_aa)
transcript_counts <- transcript_counts[order(transcript_counts$TE, decreasing=T),]

# training set: random set of 100 genes of top 200 translated genes
# test set: remaining 100 genes in top 200 translated genes
seed <- as.integer(Sys.time()) %% 1e6
print(paste("Seed:", seed))
set.seed(seed)
training_set <- sort(sample.int(n=(2*num_genes), size=num_genes))
test_set <- c(1:(2*num_genes))[-training_set]
training_set <- as.character(transcript_counts$transcript)[training_set]
test_set <- as.character(transcript_counts$transcript)[test_set]

# 4. initialize data frames for regression --------------------------------

if(!file.exists("ottrUMI_training_data.Rda")) {
  ottrUMI_training <- init_data(transcript_fa_fname, transcript_length_fname,
                                d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                                which_transcripts=training_set)
  ottrUMI_training$transcript <- relevel(ottrUMI_training$transcript, ref=training_set[1])
  ottrUMI_training$count <- count_footprints(ottrUMI_bam, ottrUMI_training, "count")
  save(ottrUMI_training, file="ottrUMI_training_data.Rda")
} else {
  load("ottrUMI_training_data.Rda")
  training_set <- levels(ottrUMI_training$transcript)
  training_set <- training_set[order(transcript_counts$TE[match(training_set,
                                                                transcript_counts$transcript)],
                                     decreasing=T)]
}

if(!file.exists("ottrUMI_test_data.Rda")) {
  ottrUMI_test <- init_data(transcript_fa_fname, transcript_length_fname,
                            d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                            which_transcripts=test_set)
  ottrUMI_test$count <- count_footprints(ottrUMI_bam, ottrUMI_test, "count")
  save(ottrUMI_test, file="ottrUMI_test_data.Rda")
} else {
  load("ottrUMI_test_data.Rda")
  test_set <- levels(ottrUMI_test$transcript)
  test_set <- test_set[order(transcript_counts$TE[match(test_set,
                                                        transcript_counts$transcript)],
                             decreasing=T)]
}

# 5. compute regressions --------------------------------------------------

if(!file.exists("ottrUMI_fit_200.Rda")) {
  # check if reference level for d5 or d3 are multiples of 3
  if((as.numeric(levels(ottrUMI_training$d5)[1]) %% 3) == 0) {
    ottrUMI_training$d5 <- relevel(ottrUMI_training$d5, ref=levels(ottrUMI_training$d5)[2])
  }
  if((as.numeric(levels(ottrUMI_training$d3)[1]) %% 3) == 0) {
    ottrUMI_training$d3 <- relevel(ottrUMI_training$d3, ref=levels(ottrUMI_training$d3)[2])
  }
  ottrUMI_fit_200 <- MASS::glm.nb(interxn_model,
                                  data=subset(ottrUMI_training,
                                              transcript %in% training_set[1:num_genes]))
  if(!("try-error" %in% class(ottrUMI_fit_200))) {
    save(ottrUMI_fit_200, file="ottrUMI_fit_200.Rda")
    # save regression coefficients
    ottrUMI_coef_200 <- coef(ottrUMI_fit_200)
    save(ottrUMI_coef_200, file="ottrUMI_coef_200.Rda")
  }
} else {
  load("ottrUMI_fit_200.Rda")
}

# 6. correct counts -------------------------------------------------------

ottrUMI_bam$correct_200 <- correct_bias(ottrUMI_bam, ottrUMI_fit_200)
ottrUMI_training$correct_200 <- count_footprints(ottrUMI_bam,
                                                 ottrUMI_training,
                                                 "correct_200")
ottrUMI_test$correct_200 <- count_footprints(ottrUMI_bam,
                                             ottrUMI_test,
                                             "correct_200")

save(ottrUMI_bam, file="ottrUMI_bam.Rda")
save(ottrUMI_training, file="ottrUMI_training_data.Rda")
save(ottrUMI_test, file="ottrUMI_test_data.Rda")

# 7. evaluate bias --------------------------------------------------------

# training data: codons
ottrUMI_training_codon_corr <- lapply(c("count", "correct_200"),
                                      function(x) {
                                        evaluate_bias(ottrUMI_training, which_column=x,
                                                      transcript_fa_fname, transcript_length_fname,
                                                      type="codon")
                                      })

ottrUMI_training_codon_corr <- do.call(cbind, ottrUMI_training_codon_corr)
colnames(ottrUMI_training_codon_corr) <- c("uncorrected", "corrected")
save(ottrUMI_training_codon_corr, file="ottrUMI_training_codon_corr.Rda")

# training data: nucleotides
ottrUMI_training_nt_corr <- lapply(c("count", "correct_200"),
                                   function(x) {
                                     evaluate_bias(ottrUMI_training, which_column=x,
                                                   transcript_fa_fname, transcript_length_fname,
                                                   type="nt")
                                   })
ottrUMI_training_nt_corr <- do.call(cbind, ottrUMI_training_nt_corr)
colnames(ottrUMI_training_nt_corr) <- c("uncorrected", "corrected")
save(ottrUMI_training_nt_corr, file="ottrUMI_training_nt_corr.Rda")

# testing data: codons
ottrUMI_testing_codon_corr <- lapply(c("count", "correct_200"),
                                     function(x) {
                                       evaluate_bias(ottrUMI_test, which_column=x,
                                                     transcript_fa_fname, transcript_length_fname,
                                                     type="codon")
                                     })
ottrUMI_testing_codon_corr <- do.call(cbind, ottrUMI_testing_codon_corr)
colnames(ottrUMI_testing_codon_corr) <- c("uncorrected", "corrected")
save(ottrUMI_testing_codon_corr, file="ottrUMI_testing_codon_corr.Rda")

# testing data: nt
ottrUMI_testing_nt_corr <- lapply(c("count", "correct_200"),
                                  function(x) {
                                    evaluate_bias(ottrUMI_test, which_column=x,
                                                  transcript_fa_fname, transcript_length_fname,
                                                  type="nt")
                                  })
ottrUMI_testing_nt_corr <- do.call(cbind, ottrUMI_testing_nt_corr)
colnames(ottrUMI_testing_nt_corr) <- c("uncorrected", "corrected")
save(ottrUMI_testing_nt_corr, file="ottrUMI_testing_nt_corr.Rda")

# 8. make plots -----------------------------------------------------------

# training data
ottrUMI_training_codon_plots <- lapply(c("uncorrected", "corrected"),
                                       function(x) {
                                         tmp <- ottrUMI_training_codon_corr[, x]
                                         names(tmp) <- rownames(ottrUMI_training_codon_corr)
                                         plot_bias(tmp, plot_subtitle=x, type="codon")
                                       })
ottrUMI_training_nt_plots <- lapply(c("uncorrected", "corrected"),
                                    function(x) {
                                      tmp <- ottrUMI_training_nt_corr[, x]
                                      names(tmp) <- rownames(ottrUMI_training_nt_corr)
                                      plot_bias(tmp, plot_subtitle=x, type="nt")
                                    })

ottrUMI_training_codon_plots[[1]] <- ottrUMI_training_codon_plots[[1]] + ggtitle("ottrUMI: training data")
ottrUMI_training_plots <- wrap_plots(ottrUMI_training_codon_plots, nrow=1) / wrap_plots(ottrUMI_training_nt_plots, nrow=1)
save(ottrUMI_training_plots, file="ottrUMI_training_plots.Rda")

# test data
ottrUMI_test_codon_plots <- lapply(c("uncorrected", "corrected"),
                                   function(x) {
                                     tmp <- ottrUMI_testing_codon_corr[, x]
                                     names(tmp) <- rownames(ottrUMI_testing_codon_corr)
                                     plot_bias(tmp, plot_subtitle=x, type="codon")
                                   })
ottrUMI_test_nt_plots <- lapply(c("uncorrected", "corrected"),
                                function(x) {
                                  tmp <- ottrUMI_testing_nt_corr[, x]
                                  names(tmp) <- rownames(ottrUMI_testing_nt_corr)
                                  plot_bias(tmp, plot_subtitle=x, type="nt")
                                })
ottrUMI_test_codon_plots[[1]] <- ottrUMI_test_codon_plots[[1]] + ggtitle("ottrUMI: test data")
ottrUMI_test_plots <- wrap_plots(ottrUMI_test_codon_plots, nrow=1) / wrap_plots(ottrUMI_test_nt_plots, nrow=1)
save(ottrUMI_test_plots, file="ottrUMI_test_plots.Rda")






q(save="no")
