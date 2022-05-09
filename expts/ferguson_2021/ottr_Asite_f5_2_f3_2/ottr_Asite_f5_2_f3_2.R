#####
# ottr data

rm(list=ls())

library(patchwork)
library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules_ottr.txt")

ottr_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "WT_OTTR_repB",
                            "ottr_trimmed_footprints.transcript.bam")

min_prop <- 0.9
f5_length <- 2
f3_length <- 2
num_genes <- 200

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)

# 1. read in footprint alignments -----------------------------------------

if(!file.exists("ottr_bam.Rda")) {
  ottr_bam <- load_bam(ottr_bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                       f5_length=f5_length, f3_length=f3_length)
  save(ottr_bam, file="ottr_bam.Rda")
} else {
  load("ottr_bam.Rda")
}

# 2. compute size/frame subsets -------------------------------------------

if(!file.exists("ottr_subsets.Rda")) {
  d5_d3 <- count_d5_d3(ottr_bam)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file="ottr_subsets.Rda")
} else {
  load("ottr_subsets.Rda")
}

# 3. establish training and test sets -------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=ottr_bam, FUN=sum)
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

if(!file.exists("ottr_training_data.Rda")) {
  ottr_training <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                             which_transcripts=training_set)
  ottr_training$transcript <- relevel(ottr_training$transcript, ref=training_set[1])
  ottr_training$count <- count_footprints(ottr_bam, ottr_training, "count")
  save(ottr_training, file="ottr_training_data.Rda")
} else {
  load("ottr_training_data.Rda")
  training_set <- levels(ottr_training$transcript)
  training_set <- training_set[order(transcript_counts$TE[match(training_set,
                                                                transcript_counts$transcript)],
                                     decreasing=T)]
}

if(!file.exists("ottr_test_data.Rda")) {
  ottr_test <- init_data(transcript_fa_fname, transcript_length_fname,
                         d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                         which_transcripts=test_set)
  ottr_test$count <- count_footprints(ottr_bam, ottr_test, "count")
  save(ottr_test, file="ottr_test_data.Rda")
} else {
  load("ottr_test_data.Rda")
  test_set <- levels(ottr_test$transcript)
  test_set <- test_set[order(transcript_counts$TE[match(test_set,
                                                        transcript_counts$transcript)],
                             decreasing=T)]
}

# 5. compute regressions --------------------------------------------------

if(!file.exists("ottr_fit_200.Rda")) {
  # check if reference level for d5 or d3 are multiples of 3
  if((as.numeric(levels(ottr_training$d5)[1]) %% 3) == 0) {
    ottr_training$d5 <- relevel(ottr_training$d5, ref=levels(ottr_training$d5)[2])
  }
  if((as.numeric(levels(ottr_training$d3)[1]) %% 3) == 0) {
    ottr_training$d3 <- relevel(ottr_training$d3, ref=levels(ottr_training$d3)[2])
  }
  ottr_fit_200 <- MASS::glm.nb(interxn_model,
                               data=subset(ottr_training,
                                           transcript %in% training_set[1:num_genes]))
  if(!("try-error" %in% class(ottr_fit_200))) {
    save(ottr_fit_200, file="ottr_fit_200.Rda")
    # save regression coefficients
    ottr_coef_200 <- coef(ottr_fit_200)
    save(ottr_coef_200, file="ottr_coef_200.Rda")
  }
} else {
  load("ottr_fit_200.Rda")
}

# 6. correct counts -------------------------------------------------------

ottr_bam$correct_200 <- correct_bias(ottr_bam, ottr_fit_200)
ottr_training$correct_200 <- count_footprints(ottr_bam,
                                              ottr_training,
                                              "correct_200")
ottr_test$correct_200 <- count_footprints(ottr_bam,
                                          ottr_test,
                                          "correct_200")

save(ottr_bam, file="ottr_bam.Rda")
save(ottr_training, file="ottr_training_data.Rda")
save(ottr_test, file="ottr_test_data.Rda")

# 7. evaluate bias --------------------------------------------------------

# training data: codons
ottr_training_codon_corr <- lapply(c("count", "correct_200"),
                                   function(x) {
                                     evaluate_bias(ottr_training, which_column=x,
                                                   transcript_fa_fname, transcript_length_fname,
                                                   type="codon")
                                   })

ottr_training_codon_corr <- do.call(cbind, ottr_training_codon_corr)
colnames(ottr_training_codon_corr) <- c("uncorrected", "corrected")
save(ottr_training_codon_corr, file="ottr_training_codon_corr.Rda")

# training data: nucleotides
ottr_training_nt_corr <- lapply(c("count", "correct_200"),
                                function(x) {
                                  evaluate_bias(ottr_training, which_column=x,
                                                transcript_fa_fname, transcript_length_fname,
                                                type="nt")
                                })
ottr_training_nt_corr <- do.call(cbind, ottr_training_nt_corr)
colnames(ottr_training_nt_corr) <- c("uncorrected", "corrected")
save(ottr_training_nt_corr, file="ottr_training_nt_corr.Rda")

# testing data: codons
ottr_testing_codon_corr <- lapply(c("count", "correct_200"),
                                  function(x) {
                                    evaluate_bias(ottr_test, which_column=x,
                                                  transcript_fa_fname, transcript_length_fname,
                                                  type="codon")
                                  })
ottr_testing_codon_corr <- do.call(cbind, ottr_testing_codon_corr)
colnames(ottr_testing_codon_corr) <- c("uncorrected", "corrected")
save(ottr_testing_codon_corr, file="ottr_testing_codon_corr.Rda")

# testing data: nt
ottr_testing_nt_corr <- lapply(c("count", "correct_200"),
                               function(x) {
                                 evaluate_bias(ottr_test, which_column=x,
                                               transcript_fa_fname, transcript_length_fname,
                                               type="nt")
                               })
ottr_testing_nt_corr <- do.call(cbind, ottr_testing_nt_corr)
colnames(ottr_testing_nt_corr) <- c("uncorrected", "corrected")
save(ottr_testing_nt_corr, file="ottr_testing_nt_corr.Rda")

# 8. make plots -----------------------------------------------------------

# training data
ottr_training_codon_plots <- lapply(c("uncorrected", "corrected"),
                                    function(x) {
                                      tmp <- ottr_training_codon_corr[, x]
                                      names(tmp) <- rownames(ottr_training_codon_corr)
                                      plot_bias(tmp, plot_subtitle=x, type="codon")
                                    })
ottr_training_nt_plots <- lapply(c("uncorrected", "corrected"),
                                 function(x) {
                                   tmp <- ottr_training_nt_corr[, x]
                                   names(tmp) <- rownames(ottr_training_nt_corr)
                                   plot_bias(tmp, plot_subtitle=x, type="nt")
                                 })

ottr_training_codon_plots[[1]] <- ottr_training_codon_plots[[1]] + ggtitle("ottr: training data")
ottr_training_plots <- wrap_plots(ottr_training_codon_plots, nrow=1) / wrap_plots(ottr_training_nt_plots, nrow=1)
save(ottr_training_plots, file="ottr_training_plots.Rda")

# test data
ottr_test_codon_plots <- lapply(c("uncorrected", "corrected"),
                                function(x) {
                                  tmp <- ottr_testing_codon_corr[, x]
                                  names(tmp) <- rownames(ottr_testing_codon_corr)
                                  plot_bias(tmp, plot_subtitle=x, type="codon")
                                })
ottr_test_nt_plots <- lapply(c("uncorrected", "corrected"),
                             function(x) {
                               tmp <- ottr_testing_nt_corr[, x]
                               names(tmp) <- rownames(ottr_testing_nt_corr)
                               plot_bias(tmp, plot_subtitle=x, type="nt")
                             })
ottr_test_codon_plots[[1]] <- ottr_test_codon_plots[[1]] + ggtitle("ottr: test data")
ottr_test_plots <- wrap_plots(ottr_test_codon_plots, nrow=1) / wrap_plots(ottr_test_nt_plots, nrow=1)
save(ottr_test_plots, file="ottr_test_plots.Rda")






q(save="no")
