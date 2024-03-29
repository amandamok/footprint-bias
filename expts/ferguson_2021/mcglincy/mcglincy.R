#####
# mcglincy data

rm(list=ls())

library(patchwork)
library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

mcglincy_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "mcglincy",
                                "mcglincy_trimmed_deduplicated_footprints.transcript.bam")

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- 200

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)

# 1. read in footprint alignments -----------------------------------------

if(!file.exists("mcglincy_bam.Rda")) {
  mcglincy_bam <- load_bam(mcglincy_bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                           f5_length=f5_length, f3_length=f3_length)
  save(mcglincy_bam, file="mcglincy_bam.Rda")
} else {
  load("mcglincy_bam.Rda")
}

# 2. compute size/frame subsets -------------------------------------------

if(!file.exists("mcglincy_subsets.Rda")) {
  d5_d3 <- count_d5_d3(mcglincy_bam)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file="mcglincy_subsets.Rda")
} else {
  load("mcglincy_subsets.Rda")
}

# 3. establish training and test sets -------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=mcglincy_bam, FUN=sum)
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

if(!file.exists("mcglincy_training_data.Rda")) {
  mcglincy_training <- init_data(transcript_fa_fname, transcript_length_fname,
                                 d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                                 which_transcripts=training_set)
  mcglincy_training$transcript <- relevel(mcglincy_training$transcript, ref=training_set[1])
  mcglincy_training$count <- count_footprints(mcglincy_bam, mcglincy_training, "count")
  save(mcglincy_training, file="mcglincy_training_data.Rda")
} else {
  load("mcglincy_training_data.Rda")
  training_set <- levels(mcglincy_training$transcript)
  training_set <- training_set[order(transcript_counts$TE[match(training_set,
                                                                transcript_counts$transcript)],
                                     decreasing=T)]
}

if(!file.exists("mcglincy_test_data.Rda")) {
  mcglincy_test <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                             which_transcripts=test_set)
  mcglincy_test$count <- count_footprints(mcglincy_bam, mcglincy_test, "count")
  save(mcglincy_test, file="mcglincy_test_data.Rda")
} else {
  load("mcglincy_test_data.Rda")
  test_set <- levels(mcglincy_test$transcript)
  test_set <- test_set[order(transcript_counts$TE[match(test_set,
                                                        transcript_counts$transcript)],
                             decreasing=T)]
}

# 5. compute regressions --------------------------------------------------

if(!file.exists("mcglincy_fit_200.Rda")) {
  mcglincy_fit_200 <- MASS::glm.nb(interxn_model,
                                   data=subset(mcglincy_training,
                                               transcript %in% training_set[1:num_genes]))
  if(!("try-error" %in% class(mcglincy_fit_200))) {
    save(mcglincy_fit_200, file="mcglincy_fit_200.Rda")
    # save regression coefficients
    mcglincy_coef_200 <- coef(mcglincy_fit_200)
    save(mcglincy_coef_200, file="mcglincy_coef_200.Rda")
  }
} else {
  load("mcglincy_fit_200.Rda")
}

# 6. correct counts -------------------------------------------------------

mcglincy_bam$correct_200 <- correct_bias(mcglincy_bam, mcglincy_fit_200)
mcglincy_training$correct_200 <- count_footprints(mcglincy_bam,
                                                  mcglincy_training,
                                                  "correct_200")
mcglincy_test$correct_200 <- count_footprints(mcglincy_bam,
                                              mcglincy_test,
                                              "correct_200")

save(mcglincy_bam, file="mcglincy_bam.Rda")
save(mcglincy_training, file="mcglincy_training_data.Rda")
save(mcglincy_test, file="mcglincy_test_data.Rda")

# 7. evaluate bias --------------------------------------------------------

# training data: codons
mcglincy_training_codon_corr <- lapply(c("count", "correct_200"),
                                       function(x) {
                                         evaluate_bias(mcglincy_training, which_column=x,
                                                       transcript_fa_fname, transcript_length_fname,
                                                       type="codon")
                                       })

mcglincy_training_codon_corr <- do.call(cbind, mcglincy_training_codon_corr)
colnames(mcglincy_training_codon_corr) <- c("uncorrected", "corrected")
save(mcglincy_training_codon_corr, file="mcglincy_training_codon_corr.Rda")

# training data: nucleotides
mcglincy_training_nt_corr <- lapply(c("count", "correct_200"),
                                    function(x) {
                                      evaluate_bias(mcglincy_training, which_column=x,
                                                    transcript_fa_fname, transcript_length_fname,
                                                    type="nt")
                                    })
mcglincy_training_nt_corr <- do.call(cbind, mcglincy_training_nt_corr)
colnames(mcglincy_training_nt_corr) <- c("uncorrected", "corrected")
save(mcglincy_training_nt_corr, file="mcglincy_training_nt_corr.Rda")

# testing data: codons
mcglincy_testing_codon_corr <- lapply(c("count", "correct_200"),
                                      function(x) {
                                        evaluate_bias(mcglincy_test, which_column=x,
                                                      transcript_fa_fname, transcript_length_fname,
                                                      type="codon")
                                      })
mcglincy_testing_codon_corr <- do.call(cbind, mcglincy_testing_codon_corr)
colnames(mcglincy_testing_codon_corr) <- c("uncorrected", "corrected")
save(mcglincy_testing_codon_corr, file="mcglincy_testing_codon_corr.Rda")

# testing data: nt
mcglincy_testing_nt_corr <- lapply(c("count", "correct_200"),
                                   function(x) {
                                     evaluate_bias(mcglincy_test, which_column=x,
                                                   transcript_fa_fname, transcript_length_fname,
                                                   type="nt")
                                   })
mcglincy_testing_nt_corr <- do.call(cbind, mcglincy_testing_nt_corr)
colnames(mcglincy_testing_nt_corr) <- c("uncorrected", "corrected")
save(mcglincy_testing_nt_corr, file="mcglincy_testing_nt_corr.Rda")

# 8. make plots -----------------------------------------------------------

# training data
mcglincy_training_codon_plots <- lapply(c("uncorrected", "corrected"),
                                        function(x) {
                                          tmp <- mcglincy_training_codon_corr[, x]
                                          names(tmp) <- rownames(mcglincy_training_codon_corr)
                                          plot_bias(tmp, plot_subtitle=x, type="codon")
                                        })
mcglincy_training_nt_plots <- lapply(c("uncorrected", "corrected"),
                                     function(x) {
                                       tmp <- mcglincy_training_nt_corr[, x]
                                       names(tmp) <- rownames(mcglincy_training_nt_corr)
                                       plot_bias(tmp, plot_subtitle=x, type="nt")
                                     })

mcglincy_training_codon_plots[[1]] <- mcglincy_training_codon_plots[[1]] + ggtitle("mcglincy: training data")
mcglincy_training_plots <- wrap_plots(mcglincy_training_codon_plots, nrow=1) / wrap_plots(mcglincy_training_nt_plots, nrow=1)
save(mcglincy_training_plots, file="mcglincy_training_plots.Rda")

# test data
mcglincy_test_codon_plots <- lapply(c("uncorrected", "corrected"),
                                    function(x) {
                                      tmp <- mcglincy_testing_codon_corr[, x]
                                      names(tmp) <- rownames(mcglincy_testing_codon_corr)
                                      plot_bias(tmp, plot_subtitle=x, type="codon")
                                    })
mcglincy_test_nt_plots <- lapply(c("uncorrected", "corrected"),
                                 function(x) {
                                   tmp <- mcglincy_testing_nt_corr[, x]
                                   names(tmp) <- rownames(mcglincy_testing_nt_corr)
                                   plot_bias(tmp, plot_subtitle=x, type="nt")
                                 })
mcglincy_test_codon_plots[[1]] <- mcglincy_test_codon_plots[[1]] + ggtitle("mcglincy: test data")
mcglincy_test_plots <- wrap_plots(mcglincy_test_codon_plots, nrow=1) / wrap_plots(mcglincy_test_nt_plots, nrow=1)
save(mcglincy_test_plots, file="mcglincy_test_plots.Rda")






q(save="no")
