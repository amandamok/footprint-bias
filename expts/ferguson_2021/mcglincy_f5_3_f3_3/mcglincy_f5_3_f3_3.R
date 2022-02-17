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
f5_length <- 3
f3_length <- 3
num_genes <- 200

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)

# 1. read in footprint alignments -----------------------------------------

if(!file.exists("mcglincy_f5_3_bam.Rda")) {
  mcglincy_f5_3_bam <- load_bam(mcglincy_bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                                f5_length=f5_length, f3_length=f3_length)
  save(mcglincy_f5_3_bam, file="mcglincy_f5_3_bam.Rda")
} else {
  load("mcglincy_f5_3_bam.Rda")
}

# 2. compute size/frame subsets -------------------------------------------

if(!file.exists("mcglincy_f5_3_subsets.Rda")) {
  d5_d3 <- count_d5_d3(mcglincy_f5_3_bam)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file="mcglincy_f5_3_subsets.Rda")
} else {
  load("mcglincy_f5_3_subsets.Rda")
}

# 3. establish training and test sets -------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=mcglincy_f5_3_bam, FUN=sum)
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

if(!file.exists("mcglincy_f5_3_training_data.Rda")) {
  mcglincy_f5_3_training <- init_data(transcript_fa_fname, transcript_length_fname,
                                      d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                                      which_transcripts=training_set)
  mcglincy_f5_3_training$transcript <- relevel(mcglincy_f5_3_training$transcript, ref=training_set[1])
  mcglincy_f5_3_training$count <- count_footprints(mcglincy_f5_3_bam, mcglincy_f5_3_training, "count")
  save(mcglincy_f5_3_training, file="mcglincy_f5_3_training_data.Rda")
} else {
  load("mcglincy_f5_3_training_data.Rda")
  training_set <- levels(mcglincy_f5_3_training$transcript)
  training_set <- training_set[order(transcript_counts$TE[match(training_set,
                                                                transcript_counts$transcript)],
                                     decreasing=T)]
}

if(!file.exists("mcglincy_f5_3_test_data.Rda")) {
  mcglincy_f5_3_test <- init_data(transcript_fa_fname, transcript_length_fname,
                                  d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                                  which_transcripts=test_set)
  mcglincy_f5_3_test$count <- count_footprints(mcglincy_f5_3_bam, mcglincy_f5_3_test, "count")
  save(mcglincy_f5_3_test, file="mcglincy_f5_3_test_data.Rda")
} else {
  load("mcglincy_f5_3_test_data.Rda")
  test_set <- levels(mcglincy_f5_3_test$transcript)
  test_set <- test_set[order(transcript_counts$TE[match(test_set,
                                                        transcript_counts$transcript)],
                             decreasing=T)]
}

# 5. compute regressions --------------------------------------------------

if(!file.exists("mcglincy_f5_3_fit_200.Rda")) {
  mcglincy_f5_3_fit_200 <- MASS::glm.nb(interxn_model,
                                        data=subset(mcglincy_f5_3_training,
                                                    transcript %in% training_set[1:num_genes]))
  if(!("try-error" %in% class(mcglincy_f5_3_fit_200))) {
    save(mcglincy_f5_3_fit_200, file="mcglincy_f5_3_fit_200.Rda")
    # save regression coefficients
    mcglincy_f5_3_coef_200 <- coef(mcglincy_f5_3_fit_200)
    save(mcglincy_f5_3_coef_200, file="mcglincy_f5_3_coef_200.Rda")
  }
} else {
  load("mcglincy_f5_3_fit_200.Rda")
}

# 6. correct counts -------------------------------------------------------

mcglincy_f5_3_bam$correct_200 <- correct_bias(mcglincy_f5_3_bam, mcglincy_f5_3_fit_200)
mcglincy_f5_3_training$correct_200 <- count_footprints(mcglincy_f5_3_bam,
                                                       mcglincy_f5_3_training,
                                                       "correct_200")
mcglincy_f5_3_test$correct_200 <- count_footprints(mcglincy_f5_3_bam,
                                                   mcglincy_f5_3_test,
                                                   "correct_200")

save(mcglincy_f5_3_bam, file="mcglincy_f5_3_bam.Rda")
save(mcglincy_f5_3_training, file="mcglincy_f5_3_training_data.Rda")
save(mcglincy_f5_3_test, file="mcglincy_f5_3_test_data.Rda")

# 7. evaluate bias --------------------------------------------------------

# training data: codons
mcglincy_f5_3_training_codon_corr <- lapply(c("count", "correct_200"),
                                            function(x) {
                                              evaluate_bias(mcglincy_f5_3_training, which_column=x,
                                                            transcript_fa_fname, transcript_length_fname,
                                                            type="codon")
                                            })

mcglincy_f5_3_training_codon_corr <- do.call(cbind, mcglincy_f5_3_training_codon_corr)
colnames(mcglincy_f5_3_training_codon_corr) <- c("uncorrected", "corrected")
save(mcglincy_f5_3_training_codon_corr, file="mcglincy_f5_3_training_codon_corr.Rda")

# training data: nucleotides
mcglincy_f5_3_training_nt_corr <- lapply(c("count", "correct_200"),
                                         function(x) {
                                           evaluate_bias(mcglincy_f5_3_training, which_column=x,
                                                         transcript_fa_fname, transcript_length_fname,
                                                         type="nt")
                                         })
mcglincy_f5_3_training_nt_corr <- do.call(cbind, mcglincy_f5_3_training_nt_corr)
colnames(mcglincy_f5_3_training_nt_corr) <- c("uncorrected", "corrected")
save(mcglincy_f5_3_training_nt_corr, file="mcglincy_f5_3_training_nt_corr.Rda")

# testing data: codons
mcglincy_f5_3_testing_codon_corr <- lapply(c("count", "correct_200"),
                                           function(x) {
                                             evaluate_bias(mcglincy_f5_3_test, which_column=x,
                                                           transcript_fa_fname, transcript_length_fname,
                                                           type="codon")
                                           })
mcglincy_f5_3_testing_codon_corr <- do.call(cbind, mcglincy_f5_3_testing_codon_corr)
colnames(mcglincy_f5_3_testing_codon_corr) <- c("uncorrected", "corrected")
save(mcglincy_f5_3_testing_codon_corr, file="mcglincy_f5_3_testing_codon_corr.Rda")

# testing data: nt
mcglincy_f5_3_testing_nt_corr <- lapply(c("count", "correct_200"),
                                        function(x) {
                                          evaluate_bias(mcglincy_f5_3_test, which_column=x,
                                                        transcript_fa_fname, transcript_length_fname,
                                                        type="nt")
                                        })
mcglincy_f5_3_testing_nt_corr <- do.call(cbind, mcglincy_f5_3_testing_nt_corr)
colnames(mcglincy_f5_3_testing_nt_corr) <- c("uncorrected", "corrected")
save(mcglincy_f5_3_testing_nt_corr, file="mcglincy_f5_3_testing_nt_corr.Rda")

# 8. make plots -----------------------------------------------------------

# training data
mcglincy_f5_3_training_codon_plots <- lapply(c("uncorrected", "corrected"),
                                             function(x) {
                                               tmp <- mcglincy_f5_3_training_codon_corr[, x]
                                               names(tmp) <- rownames(mcglincy_f5_3_training_codon_corr)
                                               plot_bias(tmp, plot_subtitle=x, type="codon")
                                             })
mcglincy_f5_3_training_nt_plots <- lapply(c("uncorrected", "corrected"),
                                          function(x) {
                                            tmp <- mcglincy_f5_3_training_nt_corr[, x]
                                            names(tmp) <- rownames(mcglincy_f5_3_training_nt_corr)
                                            plot_bias(tmp, plot_subtitle=x, type="nt")
                                          })
mcglincy_f5_3_training_codon_plots[[1]] <- mcglincy_f5_3_training_codon_plots[[1]] + ggtitle("mcglincy (5' bias: 3nt): training data")
mcglincy_f5_3_training_plots <- wrap_plots(mcglincy_f5_3_training_codon_plots, nrow=1) / wrap_plots(mcglincy_f5_3_training_nt_plots, nrow=1)
save(mcglincy_f5_3_training_plots, file="mcglincy_f5_3_training_plots.Rda")

# test data
mcglincy_f5_3_test_codon_plots <- lapply(c("uncorrected", "corrected"),
                                         function(x) {
                                           tmp <- mcglincy_f5_3_testing_codon_corr[, x]
                                           names(tmp) <- rownames(mcglincy_f5_3_testing_codon_corr)
                                           plot_bias(tmp, plot_subtitle=x, type="codon")
                                         })
mcglincy_f5_3_test_nt_plots <- lapply(c("uncorrected", "corrected"),
                                      function(x) {
                                        tmp <- mcglincy_f5_3_testing_nt_corr[, x]
                                        names(tmp) <- rownames(mcglincy_f5_3_testing_nt_corr)
                                        plot_bias(tmp, plot_subtitle=x, type="nt")
                                      })
mcglincy_f5_3_test_codon_plots[[1]] <- mcglincy_f5_3_test_codon_plots[[1]] + ggtitle("mcglincy (5' bias: 3nt): test data")
mcglincy_f5_3_test_plots <- wrap_plots(mcglincy_f5_3_test_codon_plots, nrow=1) / wrap_plots(mcglincy_f5_3_test_nt_plots, nrow=1)
save(mcglincy_f5_3_test_plots, file="mcglincy_f5_3_test_plots.Rda")


# try new correction method -----------------------------------------------

# modify bam data
mcglincy_f5_3_bam$f5_1 <- substr(as.character(mcglincy_f5_3_bam$genome_f5), 1, 1)
mcglincy_f5_3_bam$f5_2 <- substr(as.character(mcglincy_f5_3_bam$genome_f5), 2, 2)
mcglincy_f5_3_bam$f5_3 <- substr(as.character(mcglincy_f5_3_bam$genome_f5), 3, 3)
mcglincy_f5_3_bam$f3_1 <- substr(as.character(mcglincy_f5_3_bam$genome_f3), 3, 3)
mcglincy_f5_3_bam$f3_2 <- substr(as.character(mcglincy_f5_3_bam$genome_f3), 2, 2)
mcglincy_f5_3_bam$f3_3 <- substr(as.character(mcglincy_f5_3_bam$genome_f3), 1, 1)

# modify training data
mcglincy_f5_3_training$f5_1 <- substr(as.character(mcglincy_f5_3_training$genome_f5), 1, 1)
mcglincy_f5_3_training$f5_2 <- substr(as.character(mcglincy_f5_3_training$genome_f5), 2, 2)
mcglincy_f5_3_training$f5_3 <- substr(as.character(mcglincy_f5_3_training$genome_f5), 3, 3)
mcglincy_f5_3_training$f3_1 <- substr(as.character(mcglincy_f5_3_training$genome_f3), 3, 3)
mcglincy_f5_3_training$f3_2 <- substr(as.character(mcglincy_f5_3_training$genome_f3), 2, 2)
mcglincy_f5_3_training$f3_3 <- substr(as.character(mcglincy_f5_3_training$genome_f3), 1, 1)

# 3nt 5' bias and 3nt 3' bias
fit_f5_3_f3_3 <- glm.nb(formula(count ~ transcript + A + P + E + d5 + d3 +
                                  f5_1 + f5_2 + f5_3 + f3_1 + f3_2 + f3_3),
                        mcglincy_f5_3_training, model=F)
mcglincy_f5_3_bam$correct_nt_f5_3_f3_3 <- correct_bias_nt(mcglincy_f5_3_bam,
                                                          fit_f5_3_f3_3)
mcglincy_f5_3_training$correct_nt_f5_3_f3_3 <- count_footprints(mcglincy_f5_3_bam,
                                                                mcglincy_f5_3_training,
                                                                which_column="correct_nt_f5_3_f3_3")
mcglincy_f5_3_nt_bias_codon_corr <- evaluate_bias(mcglincy_f5_3_training,
                                                  which_column="correct_nt_f5_3_f3_3",
                                                  transcript_fa_fname, transcript_length_fname,
                                                  type="codon")
mcglincy_f5_3_nt_bias_nt_corr <- evaluate_bias(mcglincy_f5_3_training,
                                               which_column="correct_nt_f5_3_f3_3",
                                               transcript_fa_fname, transcript_length_fname,
                                               type="nt")

fit_f5_2_f3_3 <- glm.nb(formula(count ~ transcript + A + P + E + d5 + d3 +
                                  f5_1 + f5_2 + f3_1 + f3_2 + f3_3),
                        mcglincy_f5_3_training, model=F)
mcglincy_f5_3_bam$correct_nt_f5_2_f3_3 <- correct_bias_nt(mcglincy_f5_3_bam,
                                                          fit_f5_2_f3_3,
                                                          which_bias=c("f5_1", "f5_2", "f3_1", "f3_2", "f3_3"))
mcglincy_f5_3_training$correct_nt_f5_2_f3_3 <- count_footprints(mcglincy_f5_3_bam,
                                                                mcglincy_f5_3_training,
                                                                which_column="correct_nt_f5_2_f3_3")
mcglincy_f5_3_nt_bias_codon_corr <- evaluate_bias(mcglincy_f5_3_training,
                                                  which_column="correct_nt_f5_3_f3_3",
                                                  transcript_fa_fname, transcript_length_fname,
                                                  type="codon")
mcglincy_f5_3_nt_bias_nt_corr <- evaluate_bias(mcglincy_f5_3_training,
                                               which_column="correct_nt_f5_3_f3_3",
                                               transcript_fa_fname, transcript_length_fname,
                                               type="nt")

# no interaction term -----------------------------------------------------

fit_no_interaction <- glm.nb(formula(count ~ transcript + A + P + E + d5 + d3 + genome_f5 + genome_f3),
                             mcglincy_f5_3_training, model=F)

# load correct_bias_noInteraction() from choros
mcglincy_f5_3_bam$correct_no_intrxn <- correct_bias_noInteraction(mcglincy_f5_3_bam,
                                                                  fit_no_interaction)

mcglincy_f5_3_training$correct_no_intrxn <- count_footprints(mcglincy_f5_3_bam,
                                                             mcglincy_f5_3_training,
                                                             which_column="correct_no_intrxn")

mcglincy_f5_3_no_intrxn_codon_corr <- evaluate_bias(mcglincy_f5_3_training,
                                                    which_column="correct_no_intrxn",
                                                    transcript_fa_fname, transcript_length_fname,
                                                    type="codon")
mcglincy_f5_3_no_intrxn_nt_corr <- evaluate_bias(mcglincy_f5_3_training,
                                                 which_column="correct_no_intrxn",
                                                 transcript_fa_fname, transcript_length_fname,
                                                 type="nt")

# 5' interaction ----------------------------------------------------------

fit_f5_interaction <- glm.nb(formula(count ~ transcript + A + P + E + d5*genome_f5 + d3 + genome_f3),
                             mcglincy_f5_3_training, model=F)

# load correct_bias_f5Interaction() from choros
mcglincy_f5_3_bam$correct_f5intrxn <- correct_bias_f5Interaction(mcglincy_f5_3_bam,
                                                                 fit_f5_interaction)
mcglincy_f5_3_training$correct_f5intrxn <- count_footprints(mcglincy_f5_3_bam,
                                                            mcglincy_f5_3_training,
                                                            which_column="correct_f5intrxn")
mcglincy_f5_3_f5intrxn_codon_corr <- evaluate_bias(mcglincy_f5_3_training,
                                                   which_column="correct_f5intrxn",
                                                   transcript_fa_fname, transcript_length_fname,
                                                   type="codon")
mcglincy_f5_3_f5intrxn_nt_corr <- evaluate_bias(mcglincy_f5_3_training,
                                                which_column="correct_f5intrxn",
                                                transcript_fa_fname, transcript_length_fname,
                                                type="nt")

# shrink extreme coefficients ---------------------------------------------

# evaluate corrected (interaction term; 3nt 5' bias / 3nt 3' bias; shrink extreme coefs)
# load version of correct_bias() that shrinks coefs
load("mcglincy_f5_3_fit_200.Rda")
ggplot(parse_coefs(mcglincy_f5_3_fit_200), aes(x=type, y=value, fill=type)) +
  geom_boxplot() + theme_classic() + theme(legend.position="none") +
  xlab("") + ylab("regression coefficient")
mcglincy_f5_3_bam$correct_200 <- correct_bias(mcglincy_f5_3_bam,
                                              mcglincy_f5_3_fit_200)
mcglincy_f5_3_training$correct_200 <- count_footprints(mcglincy_f5_3_bam,
                                                       mcglincy_f5_3_training,
                                                       which_column="correct_200")
mcglincy_f5_3_corrected_codon_corr <- evaluate_bias(mcglincy_f5_3_training,
                                                    which_column="correct_200",
                                                    transcript_fa_fname, transcript_length_fname,
                                                    type="codon")
mcglincy_f5_3_corrected_nt_corr <- evaluate_bias(mcglincy_f5_3_training,
                                                 which_column="correct_200",
                                                 transcript_fa_fname, transcript_length_fname,
                                                 type="nt")

# individual bias nt + interaction ----------------------------------------

mcglincy_f5_3_nt_interxn_fit <- MASS::glm.nb(count ~ transcript + A + P + E +
                                               f5_1*d5 + f5_2*d5 + f5_3*d5 +
                                               f3_1*d3 + f3_2*d3 + f3_3*d3,
                                             mcglincy_f5_3_training,
                                             model=F)
mcglincy_f5_3_bam$correct_1nt <- correct_bias(mcglincy_f5_3_bam,
                                              mcglincy_f5_3_nt_interxn_fit,
                                              which_column="count",
                                              which_f5="f5_1", which_f3="f3_1",
                                              fit_f5="f5_1", fit_f3="f3_1")
mcglincy_f5_3_bam$correct_2nt <- correct_bias(mcglincy_f5_3_bam,
                                              mcglincy_f5_3_nt_interxn_fit,
                                              which_column="correct_1nt",
                                              which_f5="f5_2", which_f3="f3_2",
                                              fit_f5="f5_2", fit_f3="f3_2")
mcglincy_f5_3_bam$correct_3nt <- correct_bias(mcglincy_f5_3_bam,
                                              mcglincy_f5_3_nt_interxn_fit,
                                              which_column="correct_2nt",
                                              which_f5="f5_3", which_f3="f3_3",
                                              fit_f5="f5_3", fit_f3="f3_3")
mcglincy_f5_3_bam$correct_nt2 <- correct_bias(mcglincy_f5_3_bam,
                                               mcglincy_f5_3_nt_interxn_fit,
                                               which_column="count",
                                               which_f5="f5_2", which_f3="f3_2",
                                               fit_f5="f5_2", fit_f3="f3_2")
mcglincy_f5_3_bam$correct_nt23 <- correct_bias(mcglincy_f5_3_bam,
                                               mcglincy_f5_3_nt_interxn_fit,
                                               which_column="correct_nt2",
                                               which_f5="f5_3", which_f3="f3_3",
                                               fit_f5="f5_3", fit_f3="f3_3")
mcglincy_f5_3_bam$rpf_f5_1 <- substr(as.character(mcglincy_f5_3_bam$rpf_f5), 1, 1)
mcglincy_f5_3_bam$correct_1nt_rpf_f5 <- correct_bias(mcglincy_f5_3_bam,
                                                       mcglincy_f5_3_nt_interxn_fit,
                                                       which_column="count",
                                                       which_f5="rpf_f5_1", which_f3="f3_1",
                                                       fit_f5="f5_1", fit_f3="f3_1")
mcglincy_f5_3_bam$correct_2nt_rpf_f5 <- correct_bias(mcglincy_f5_3_bam,
                                                     mcglincy_f5_3_nt_interxn_fit,
                                                     which_column="correct_1nt_rpf_f5",
                                                     which_f5="f5_2", which_f3="f3_2",
                                                     fit_f5="f5_2", fit_f3="f3_2")
mcglincy_f5_3_bam$correct_3nt_rpf_f5 <- correct_bias(mcglincy_f5_3_bam,
                                                     mcglincy_f5_3_nt_interxn_fit,
                                                     which_column="correct_2nt_rpf_f5",
                                                     which_f5="f5_3", which_f3="f3_3",
                                                     fit_f5="f5_3", fit_f3="f3_3")

mcglincy_f5_3_training$correct_1nt <- count_footprints(mcglincy_f5_3_bam,
                                                       mcglincy_f5_3_training,
                                                       which_column="correct_1nt")
mcglincy_f5_3_training$correct_2nt <- count_footprints(mcglincy_f5_3_bam,
                                                       mcglincy_f5_3_training,
                                                       which_column="correct_2nt")
mcglincy_f5_3_training$correct_3nt <- count_footprints(mcglincy_f5_3_bam,
                                                       mcglincy_f5_3_training,
                                                       which_column="correct_3nt")
mcglincy_f5_3_training$correct_nt23 <- count_footprints(mcglincy_f5_3_bam,
                                                        mcglincy_f5_3_training,
                                                        which_column="correct_nt23")
mcglincy_f5_3_training$correct_1nt_rpf_f5 <- count_footprints(mcglincy_f5_3_bam,
                                                              mcglincy_f5_3_training,
                                                              which_column="correct_1nt_rpf_f5")
mcglincy_f5_3_training$correct_2nt_rpf_f5 <- count_footprints(mcglincy_f5_3_bam,
                                                              mcglincy_f5_3_training,
                                                              which_column="correct_2nt_rpf_f5")
mcglincy_f5_3_training$correct_3nt_rpf_f5 <- count_footprints(mcglincy_f5_3_bam,
                                                              mcglincy_f5_3_training,
                                                              which_column="correct_3nt_rpf_f5")

nt_interxn <- c(paste0("correct_", 1:3, "nt"), "correct_nt23",
                paste0("correct_", 1:3, "nt_rpf_f5"))
mcglincy_f5_3_nt_interxn_codon_corr <- lapply(nt_interxn,
                                              function(x) {
                                                evaluate_bias(mcglincy_f5_3_training,
                                                              which_column=x,
                                                              transcript_fa_fname,
                                                              transcript_length_fname,
                                                              type="codon")
                                              })
names(mcglincy_f5_3_nt_interxn_codon_corr) <- nt_interxn
mcglincy_f5_3_nt_interxn_nt_corr <- lapply(nt_interxn,
                                              function(x) {
                                                evaluate_bias(mcglincy_f5_3_training,
                                                              which_column=x,
                                                              transcript_fa_fname,
                                                              transcript_length_fname,
                                                              type="nt")
                                              })
names(mcglincy_f5_3_nt_interxn_nt_corr) <- nt_interxn

save(mcglincy_f5_3_bam, mcglincy_f5_3_training,
     mcglincy_f5_3_nt_interxn_fit,
     mcglincy_f5_3_nt_interxn_codon_corr, mcglincy_f5_3_nt_interxn_nt_corr,
     file="mcglincy_f5_3_nt_interxn.Rda")

# plot --------------------------------------------------------------------

# evaluate raw codon & nt correlations
mcglincy_f5_3_raw_codon_corr <- evaluate_bias(mcglincy_f5_3_training,
                                              which_column="count",
                                              transcript_fa_fname, transcript_length_fname,
                                              type="codon")
mcglincy_f5_3_raw_nt_corr <- evaluate_bias(mcglincy_f5_3_training,
                                           which_column="count",
                                           transcript_fa_fname, transcript_length_fname,
                                           type="nt")

# evaluate bias correction: interaction model w/ 2nt 5' bias, 3nt 3' bias
load("../mcglincy/mcglincy_training_data.Rda")
mcglincy_corrected_codon_corr <- evaluate_bias(mcglincy_training,
                                               which_column="correct_200",
                                               transcript_fa_fname, transcript_length_fname,
                                               type="codon")
mcglincy_corrected_nt_corr <- evaluate_bias(mcglincy_training,
                                            which_column="correct_200",
                                            transcript_fa_fname, transcript_length_fname,
                                            type="nt")


(plot_bias(mcglincy_f5_3_raw_codon_corr) + ylim(0, 0.28) +
    ggtitle("McGlincy", subtitle="raw")) +
  (plot_bias(mcglincy_corrected_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="interaction: 2nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_no_intrxn_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="no interaction: 3nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_nt_bias_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="individual nt: 3nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_raw_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_corrected_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_no_intrxn_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_nt_bias_nt_corr, type="nt") + ylim(0, 0.14)) +
  plot_layout(nrow=2)

(plot_bias(mcglincy_f5_3_raw_codon_corr) + ylim(0, 0.28) +
    ggtitle("McGlincy", subtitle="raw")) +
  (plot_bias(mcglincy_corrected_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="interaction: 2nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_corrected_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="interaction ***: 3nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_f5intrxn_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="5' interaction ***: 3nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_raw_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_corrected_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_corrected_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_f5intrxn_nt_corr, type="nt") + ylim(0, 0.14)) +
  plot_layout(nrow=2)

(plot_bias(mcglincy_f5_3_raw_codon_corr) + ylim(0, 0.28) +
    ggtitle("McGlincy", subtitle="raw")) +
  (plot_bias(mcglincy_corrected_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="interaction: 2nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_corrected_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="interaction ***: 3nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_f5intrxn_codon_corr) + ylim(0, 0.28) +
     ggtitle("", subtitle="5' interaction ***: 3nt 5' bias / 3nt 3' bias")) +
  (plot_bias(mcglincy_f5_3_raw_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_corrected_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_corrected_nt_corr, type="nt") + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_f5intrxn_nt_corr, type="nt") + ylim(0, 0.14)) +
  plot_layout(nrow=2)

(plot_bias(mcglincy_f5_3_raw_codon_corr) + ylim(0, 0.28) +
    ggtitle("McGlincy", subtitle="raw")) +
  (plot_bias(mcglincy_f5_3_nt_interxn_codon_corr$correct_1nt) + ylim(0, 0.28) +
     ggtitle("", subtitle="nt interxn: 1nt")) +
  (plot_bias(mcglincy_f5_3_nt_interxn_codon_corr$correct_2nt) + ylim(0, 0.28) +
     ggtitle("", subtitle="nt interxn: 2nt")) +
  (plot_bias(mcglincy_f5_3_nt_interxn_codon_corr$correct_3nt) + ylim(0, 0.28) +
     ggtitle("", subtitle="nt interxn: 3nt")) +
  (plot_bias(mcglincy_f5_3_nt_interxn_codon_corr$correct_nt23) + ylim(0, 0.28) +
     ggtitle("", subtitle="nt interxn: nt23")) +
  (plot_bias(mcglincy_f5_3_raw_nt_corr) + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_nt_interxn_nt_corr$correct_1nt) + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_nt_interxn_nt_corr$correct_2nt) + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_nt_interxn_nt_corr$correct_3nt) + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_nt_interxn_nt_corr$correct_nt23) + ylim(0, 0.14)) +
  plot_layout(nrow=2)

(plot_bias(mcglincy_f5_3_raw_codon_corr) + ylim(0, 0.28) +
    ggtitle("McGlincy", subtitle="raw")) +
  (plot_bias(mcglincy_f5_3_nt_interxn_codon_corr$correct_1nt_rpf_f5) + ylim(0, 0.28) +
     ggtitle("", subtitle="nt interxn: 1nt (RPF 1st nt)")) +
  (plot_bias(mcglincy_f5_3_nt_interxn_codon_corr$correct_2nt_rpf_f5) + ylim(0, 0.28) +
     ggtitle("", subtitle="nt interxn: 2nt (RPF 1st nt)")) +
  (plot_bias(mcglincy_f5_3_nt_interxn_codon_corr$correct_3nt_rpf_f5) + ylim(0, 0.28) +
     ggtitle("", subtitle="nt interxn: 3nt (RPF 1st nt)")) +
  (plot_bias(mcglincy_f5_3_raw_nt_corr) + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_nt_interxn_nt_corr$correct_1nt_rpf_f5) + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_nt_interxn_nt_corr$correct_2nt_rpf_f5) + ylim(0, 0.14)) +
  (plot_bias(mcglincy_f5_3_nt_interxn_nt_corr$correct_3nt_rpf_f5) + ylim(0, 0.14)) +
  plot_layout(nrow=2)

(plot_bias(mcglincy_f5_3_raw_codon_corr) + ylim(0, 0.28) +
    ggtitle("McGlincy", subtitle="raw")) +
  (plot_bias(mcglincy_f5_3_nt_interxn_codon_corr$correct_1nt_rpf_f5))

save(list=c(grep("mcglincy", ls(), value=T),
            grep("fit", ls(), value=T)),
     file="all_fits.Rda")

q(save="no")
