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

ottrUMI_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "WT_UMI_OTTR_repB",
                               "ottrUMI_trimmed_deduplicated_footprints.transcript.bam")

min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 200

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)

# 1. read in footprint alignments -----------------------------------------

if(!file.exists("ottrUMI_f5_3_bam.Rda")) {
  ottrUMI_f5_3_bam <- load_bam(ottrUMI_bam_fname, transcript_fa_fname,
                               transcript_length_fname, offsets_fname,
                               f5_length=f5_length, f3_length=f3_length)
  save(ottrUMI_f5_3_bam, file="ottrUMI_f5_3_bam.Rda")
} else {
  load("ottrUMI_f5_3_bam.Rda")
}

# 2. compute size/frame subsets -------------------------------------------

if(!file.exists("ottrUMI_f5_3_subsets.Rda")) {
  d5_d3 <- count_d5_d3(ottrUMI_f5_3_bam)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file="ottrUMI_f5_3_subsets.Rda")
} else {
  load("ottrUMI_f5_3_subsets.Rda")
}

# 3. establish training and test sets -------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=ottrUMI_f5_3_bam, FUN=sum)
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

if(!file.exists("ottrUMI_f5_3_training_data.Rda")) {
  ottrUMI_f5_3_training <- init_data(transcript_fa_fname, transcript_length_fname,
                                     d5_d3_subsets=d5_d3_subsets,
                                     f5_length=f5_length, f3_length=f3_length,
                                     which_transcripts=training_set)
  ottrUMI_f5_3_training$transcript <- relevel(ottrUMI_f5_3_training$transcript, ref=training_set[1])
  ottrUMI_f5_3_training$count <- count_footprints(ottrUMI_f5_3_bam, ottrUMI_f5_3_training, "count")
  save(ottrUMI_f5_3_training, file="ottrUMI_f5_3_training_data.Rda")
} else {
  load("ottrUMI_f5_3_training_data.Rda")
  training_set <- levels(ottrUMI_f5_3_training$transcript)
  training_set <- training_set[order(transcript_counts$TE[match(training_set,
                                                                transcript_counts$transcript)],
                                     decreasing=T)]
}

if(!file.exists("ottrUMI_f5_3_test_data.Rda")) {
  ottrUMI_f5_3_test <- init_data(transcript_fa_fname, transcript_length_fname,
                                 d5_d3_subsets=d5_d3_subsets,
                                 f5_length=f5_length, f3_length=f3_length,
                                 which_transcripts=test_set)
  ottrUMI_f5_3_test$count <- count_footprints(ottrUMI_f5_3_bam, ottrUMI_f5_3_test, "count")
  save(ottrUMI_f5_3_test, file="ottrUMI_f5_3_test_data.Rda")
} else {
  load("ottrUMI_f5_3_test_data.Rda")
  test_set <- levels(ottrUMI_f5_3_test$transcript)
  test_set <- test_set[order(transcript_counts$TE[match(test_set,
                                                        transcript_counts$transcript)],
                             decreasing=T)]
}

# 5. compute regressions --------------------------------------------------

if(!file.exists("ottrUMI_f5_3_fit_200.Rda")) {
  # check if reference level for d5 or d3 are multiples of 3
  if((as.numeric(levels(ottrUMI_f5_3_training$d5)[1]) %% 3) == 0) {
    ottrUMI_f5_3_training$d5 <- relevel(ottrUMI_f5_3_training$d5,
                                        ref=levels(ottrUMI_f5_3_training$d5)[2])
  }
  if((as.numeric(levels(ottrUMI_f5_3_training$d3)[1]) %% 3) == 0) {
    ottrUMI_f5_3_training$d3 <- relevel(ottrUMI_f5_3_training$d3,
                                        ref=levels(ottrUMI_f5_3_training$d3)[2])
  }
  ottrUMI_f5_3_fit_200 <- MASS::glm.nb(interxn_model,
                                       data=subset(ottrUMI_f5_3_training,
                                                   transcript %in% training_set[1:num_genes]))
  if(!("try-error" %in% class(ottrUMI_f5_3_fit_200))) {
    save(ottrUMI_f5_3_fit_200, file="ottrUMI_f5_3_fit_200.Rda")
    # save regression coefficients
    ottrUMI_f5_3_coef_200 <- coef(ottrUMI_f5_3_fit_200)
    save(ottrUMI_f5_3_coef_200, file="ottrUMI_f5_3_coef_200.Rda")
  }
} else {
  load("ottrUMI_f5_3_fit_200.Rda")
}

# 6. correct counts -------------------------------------------------------

ottrUMI_f5_3_bam$correct_200 <- correct_bias(ottrUMI_f5_3_bam, ottrUMI_f5_3_fit_200)
ottrUMI_f5_3_training$correct_200 <- count_footprints(ottrUMI_f5_3_bam,
                                                      ottrUMI_f5_3_training,
                                                      "correct_200")
ottrUMI_f5_3_test$correct_200 <- count_footprints(ottrUMI_f5_3_bam,
                                                  ottrUMI_f5_3_test,
                                                  "correct_200")

save(ottrUMI_f5_3_bam, file="ottrUMI_f5_3_bam.Rda")
save(ottrUMI_f5_3_training, file="ottrUMI_f5_3_training_data.Rda")
save(ottrUMI_f5_3_test, file="ottrUMI_f5_3_test_data.Rda")

# 7. evaluate bias --------------------------------------------------------

# training data: codons
ottrUMI_f5_3_training_codon_corr <- lapply(c("count", "correct_200"),
                                           function(x) {
                                             evaluate_bias(ottrUMI_f5_3_training,
                                                           which_column=x,
                                                           transcript_fa_fname,
                                                           transcript_length_fname,
                                                           type="codon")
                                           })

ottrUMI_f5_3_training_codon_corr <- do.call(cbind, ottrUMI_f5_3_training_codon_corr)
colnames(ottrUMI_f5_3_training_codon_corr) <- c("uncorrected", "corrected")
save(ottrUMI_f5_3_training_codon_corr, file="ottrUMI_f5_3_training_codon_corr.Rda")

# training data: nucleotides
ottrUMI_f5_3_training_nt_corr <- lapply(c("count", "correct_200"),
                                        function(x) {
                                          evaluate_bias(ottrUMI_f5_3_training,
                                                        which_column=x,
                                                        transcript_fa_fname,
                                                        transcript_length_fname,
                                                        type="nt")
                                        })
ottrUMI_f5_3_training_nt_corr <- do.call(cbind, ottrUMI_f5_3_training_nt_corr)
colnames(ottrUMI_f5_3_training_nt_corr) <- c("uncorrected", "corrected")
save(ottrUMI_f5_3_training_nt_corr, file="ottrUMI_f5_3_training_nt_corr.Rda")

# 8. make plots -----------------------------------------------------------

load("../ottrUMI_Asite/ottrUMI_training_codon_corr.Rda")
ottrUMI_f5_2_f3_3_training_codon_corr <- ottrUMI_training_codon_corr
load("../ottrUMI_Asite/ottrUMI_training_nt_corr.Rda")
ottrUMI_f5_2_f3_3_training_nt_corr <- ottrUMI_training_nt_corr

load("../ottrUMI_Asite_f5_2_f3_2/ottrUMI_training_codon_corr.Rda")
ottrUMI_f5_2_f3_2_training_codon_corr <- ottrUMI_training_codon_corr
load("../ottrUMI_Asite_f5_2_f3_2/ottrUMI_training_nt_corr.Rda")
ottrUMI_f5_2_f3_2_training_nt_corr <- ottrUMI_training_nt_corr

(plot_bias(ottrUMI_f5_3_training_codon_corr[,1]) +
    ggtitle("ottrUMI", subtitle="raw") + ylim(0, 0.075)) +
  (plot_bias(ottrUMI_f5_2_f3_3_training_codon_corr[,2]) +
     ggtitle("", subtitle="5' 2nt / 3' 3nt") + ylim(0, 0.075)) +
  (plot_bias(ottrUMI_f5_3_training_codon_corr[,2]) +
     ggtitle("", subtitle="5' 3nt / 3' 3nt") + ylim(0, 0.075)) +
  (plot_bias(ottrUMI_f5_2_f3_2_training_codon_corr[,2]) +
     ggtitle("", subtitle="5' 2nt / 3' 2nt") + ylim(0, 0.075)) +
  (plot_bias(ottrUMI_f5_3_training_nt_corr[,1], type="nt") + ylim(0, 0.05)) +
  (plot_bias(ottrUMI_training_nt_corr[,2], type="nt") + ylim(0, 0.05)) +
  (plot_bias(ottrUMI_f5_3_training_nt_corr[,2], type="nt") + ylim(0, 0.05)) +
  (plot_bias(ottrUMI_f5_2_f3_2_training_nt_corr[,2], type="nt") + ylim(0, 0.05)) +
  plot_layout(nrow=2)


ottrUMI_f5_3_coefs <- parse_coefs(ottrUMI_f5_3_fit_200)
ggplot(ottrUMI_f5_3_coefs, aes(x=type, y=value, fill=type)) + geom_boxplot() + theme_classic()

q(save="no")
