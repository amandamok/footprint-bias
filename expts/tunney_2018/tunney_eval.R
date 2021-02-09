#####
# tunney data
# test consistency of regression coefficients
# test correction on genes not in training set

rm(list=ls())

if(!("choros" %in% rownames(installed.packages()))) {
  devtools::install_git("https://github.com/amandamok/choros")
}

library(patchwork)
library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

tunney_bam_fname <- "~/iXnos/expts/lareau_20cds20/process/lareau.transcript.bam"

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- c(50, 75, 100, 125, 150, 175, 200)
num_regression_genes <- max(num_genes)

interxn_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

transcript_lengths <- load_lengths(transcript_length_fname)
utr5_length <- unique(transcript_lengths$utr5_length)
utr3_length <- unique(transcript_lengths$utr3_length)

# 1. read in footprint alignments -----------------------------------------

if(!file.exists("tunney_bam.Rda")) {
  tunney_bam <- load_bam(tunney_bam_fname, transcript_fa_fname, transcript_length_fname, offsets_fname,
                         f5_length=f5_length, f3_length=f3_length)
  save(tunney_bam, file="tunney_bam.Rda")
} else {
  load("tunney_bam.Rda")
}

# 2. compute size/frame subsets -------------------------------------------

if(!file.exists("tunney_subsets.Rda")) {
  d5_d3 <- count_d5_d3(tunney_bam)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]), c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x], "d3", d5_d3_subsets$d3[x], sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names, file="tunney_subsets.Rda")
} else {
  load("tunney_subsets.Rda")
}

# 3. establish training and test sets -------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=tunney_bam, FUN=sum)
transcript_counts$length_aa <- transcript_lengths$cds_length[match(transcript_counts$transcript,
                                                                   transcript_lengths$transcript)] / 3
transcript_counts$TE <- with(transcript_counts, count / length_aa)
transcript_counts <- transcript_counts[order(transcript_counts$TE, decreasing=T),]

# training set: random set of 100 genes of top 200 translated genes
# test set: remaining 100 genes in top 200 translated genes
seed <- as.integer(Sys.time()) %% 1e6
print(paste("Seed:", seed))
set.seed(seed)
training_set <- sort(sample.int(n=(2*num_regression_genes), size=num_regression_genes))
test_set <- c(1:(2*num_regression_genes))[-training_set]
training_set <- as.character(transcript_counts$transcript)[training_set]
test_set <- as.character(transcript_counts$transcript)[test_set]

# 4. initialize data frames for regression --------------------------------

if(!file.exists("tunney_training_data.Rda")) {
  tunney_training <- init_data(transcript_fa_fname, transcript_length_fname,
                               d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                               which_transcripts=training_set)
  tunney_training$transcript <- relevel(tunney_training$transcript, ref=training_set[1])
  tunney_training$count <- count_footprints(tunney_bam, tunney_training, "count")
  save(tunney_training, file="tunney_training_data.Rda")
} else {
  load("tunney_training_data.Rda")
  training_set <- levels(tunney_training$transcript)
  training_set <- training_set[order(transcript_counts$TE[match(training_set,
                                                                transcript_counts$transcript)],
                                     decreasing=T)]
}

if(!file.exists("tunney_test_data.Rda")) {
  tunney_test <- init_data(transcript_fa_fname, transcript_length_fname,
                           d5_d3_subsets=d5_d3_subsets, f5_length=f5_length, f3_length=f3_length,
                           which_transcripts=test_set)
  tunney_test$count <- count_footprints(tunney_bam, tunney_test, "count")
  save(tunney_test, file="tunney_test_data.Rda")
} else {
  load("tunney_test_data.Rda")
  test_set <- levels(tunney_test$transcript)
  test_set <- test_set[order(transcript_counts$TE[match(test_set,
                                                           transcript_counts$transcript)],
                             decreasing=T)]
}

# 5. compute regressions --------------------------------------------------

for(x in num_genes) {
  print(paste(x, "genes"))
  fit_obj <- paste0("tunney_fit_", x)
  fit_fname <- paste0(fit_obj, ".Rda")
  if(!file.exists(fit_fname)) {
    assign(fit_obj,
           try(MASS::glm.nb(interxn_model,
                            data=subset(tunney_training, transcript %in% training_set[1:x]),
                            model=F)))
    if(!("try-error" %in% class(get(fit_obj)))) {
      save(list=fit_obj, file=fit_fname)
      # save regression coefficients
      coef_obj <- paste0("tunney_coef_", x)
      assign(coef_obj, coef(get(fit_obj)))
      save(list=coef_obj, file=paste0(coef_obj, ".Rda"))
    }
  } else {
    load(fit_fname)
  }
}

# figure out which regression computed
num_genes <- grep("tunney_fit_", list.files(), value=T)
num_genes <- sort(as.numeric(sub("tunney_fit_", "", sub("\\.Rda", "", num_genes))))
print(paste("regression fit for:", paste(num_genes, collapse=", "), "genes"))

# 6. correct counts -------------------------------------------------------

for(x in num_genes) {
  correction <- paste0("correct_", x)
  tunney_bam[, correction] <- correct_bias(tunney_bam, get(paste0("tunney_fit_", x)))
  tunney_training[, correction] <- count_footprints(tunney_bam, tunney_training, correction)
  tunney_test[, correction] <- count_footprints(tunney_bam, tunney_test, correction)
}
save(tunney_bam, file="tunney_bam.Rda")
save(tunney_training, file="tunney_training_data.Rda")
save(tunney_test, file="tunney_test_data.Rda")

# 7. evaluate bias --------------------------------------------------------

# training data: codons
tunney_training_codon_corr <- lapply(num_genes,
                                     function(x) {
                                       evaluate_bias(tunney_training, which_column=paste0("correct_", x),
                                                     transcript_fa_fname, transcript_length_fname,
                                                     utr5=utr5_length, utr3=utr3_length, type="codon")
                                     })
tunney_training_codon_corr <- data.frame(do.call(cbind, tunney_training_codon_corr))
colnames(tunney_training_codon_corr) <- paste0("correct_", num_genes)
save(tunney_training_codon_corr, file="tunney_training_codon_corr.Rda")

# training data: nucleotides
tunney_training_nt_corr <- lapply(num_genes,
                                  function(x) {
                                    evaluate_bias(tunney_training, which_column=paste0("correct_", x),
                                                  transcript_fa_fname, transcript_length_fname,
                                                  utr5=utr5_length, utr3=utr3_length, type="nt")
                                  })
tunney_training_nt_corr <- data.frame(do.call(cbind, tunney_training_nt_corr))
colnames(tunney_training_nt_corr) <- paste0("correct_", num_genes)
save(tunney_training_nt_corr, file="tunney_training_nt_corr.Rda")

# test data: codons
tunney_test_codon_corr <- lapply(num_genes,
                                 function(x) {
                                   evaluate_bias(tunney_test, which_column=paste0("correct_", x),
                                                 transcript_fa_fname, transcript_length_fname,
                                                 utr5=utr5_length, utr3=utr3_length, type="codon")
                                 })
tunney_test_codon_corr <- data.frame(do.call(cbind, tunney_test_codon_corr))
colnames(tunney_test_codon_corr) <- paste0("correct_", num_genes)
save(tunney_test_codon_corr, file="tunney_test_codon_corr.Rda")

# test data: nucleotides
tunney_test_nt_corr <- lapply(num_genes,
                              function(x) {
                                evaluate_bias(tunney_test, which_column=paste0("correct_", x),
                                              transcript_fa_fname, transcript_length_fname,
                                              utr5=utr5_length, utr3=utr3_length, type="nt")
                              })
tunney_test_nt_corr <- data.frame(do.call(cbind, tunney_test_nt_corr))
colnames(tunney_test_nt_corr) <- paste0("correct_", num_genes)
save(tunney_test_nt_corr, file="tunney_test_nt_corr.Rda")

# 8. make plots -----------------------------------------------------------

# training data
tunney_training_codon_plots <- lapply(num_genes,
                                      function(x) {
                                        tmp <- tunney_training_codon_corr[, paste0("correct_", x)]
                                        names(tmp) <- rownames(tunney_training_codon_corr)
                                        plot_bias(tmp, plot_subtitle=paste(x, "genes"), type="codon")
                                      })
tunney_training_nt_plots <- lapply(num_genes,
                                   function(x) {
                                     tmp <- tunney_training_nt_corr[, paste0("correct_", x)]
                                     names(tmp) <- rownames(tunney_training_nt_corr)
                                     plot_bias(tmp, plot_subtitle=paste(x, "genes"), type="nt")
                                   })
tunney_training_codon_plots[[1]] <- tunney_training_codon_plots[[1]] + ggtitle("tunney: training data")
tunney_training_plots <- wrap_plots(tunney_training_codon_plots, nrow=1) / wrap_plots(tunney_training_nt_plots, nrow=1)
save(tunney_training_plots, file="tunney_training_plots.Rda")

# test data
tunney_test_codon_plots <- lapply(num_genes,
                                  function(x) {
                                    tmp <- tunney_test_codon_corr[, paste0("correct_", x)]
                                    names(tmp) <- rownames(tunney_test_codon_corr)
                                    plot_bias(tmp, plot_subtitle=paste(x, "genes"), type="codon")
                                  })
tunney_test_nt_plots <- lapply(num_genes,
                               function(x) {
                                 tmp <- tunney_test_nt_corr[, paste0("correct_", x)]
                                 names(tmp) <- rownames(tunney_test_nt_corr)
                                 plot_bias(tmp, plot_subtitle=paste(x, "genes"), type="nt")
                               })
tunney_test_codon_plots[[1]] <- tunney_test_codon_plots[[1]] + ggtitle("tunney: test data")
tunney_test_plots <- wrap_plots(tunney_test_codon_plots, nrow=1) / wrap_plots(tunney_test_nt_plots, nrow=1)
save(tunney_test_plots, file="tunney_test_plots.Rda")

# knit markdown output ----------------------------------------------------

rmarkdown::render("tunney_eval.Rmd", output_format="html_document")

q(save="no")
