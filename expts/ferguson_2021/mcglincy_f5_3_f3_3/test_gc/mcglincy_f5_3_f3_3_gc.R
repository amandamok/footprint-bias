rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_length_fname)
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

transcript_seq <- load_fasta(transcript_fa_fname)
transcript_lengths <- load_lengths(transcript_length_fname)

regression_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3 + gc)

# functions ---------------------------------------------------------------

compute_rpf_gc <- function(dat, omit="A", transcript_seq, transcript_lengths) {
  # compute GC content in RPF, omitting A/P/E sites
  # dat: data.frame; contains columns c("transcript", "cod_idx", "d5", "d3")
  ## omit: character; which codon sites to omit, one of c("A", "AP", "APE")
  ## transcript_seq: named character vector; output from load_fasta()
  ## transcript_lengths: data.frame; output from load_lengths()
  utr5_lengths <- transcript_lengths$utr5_length
  names(utr5_lengths) <- transcript_lengths$transcript
  A_start <- utr5_lengths[as.character(dat$transcript)] + 1 + 3*(dat$cod_idx-1)
  num_omit_codons <- ifelse(grepl("E", omit), 2,
                            ifelse(grepl("P", omit), 1, 0))
  dat$d5 <- as.numeric(as.character(dat$d5))
  dat$d3 <- as.numeric(as.character(dat$d3))
  # 1. extract RPF 5' regions
  rpf_5_start <- A_start - dat$d5
  rpf_5_end <- A_start - 1 - 3*num_omit_codons
  rpf_5 <- mapply(substr, transcript_seq[as.character(dat$transcript)],
                  rpf_5_start, rpf_5_end)
  # 2. extract RPF 3' regions
  rpf_3_start <- A_start + 3
  rpf_3_end <- A_start + 2 + dat$d3
  rpf_3 <- mapply(substr, transcript_seq[as.character(dat$transcript)],
                  rpf_3_start, rpf_3_end)
  # 3. compute GC content
  rpf_regions <- paste0(rpf_5, rpf_3)
  gc_content <- strsplit(rpf_regions, split="")
  gc_content <- sapply(gc_content, function(x) sum(x %in% c("G", "C")))
  gc_content <- gc_content / with(dat, d5 + d3 + 3)
  return(gc_content)
}

# load data ---------------------------------------------------------------

load("../training_400/mcglincy_f5_3_bam.Rda")
load("../training_400/mcglincy_f5_3_training_data.Rda")

# annotate GC content -----------------------------------------------------

mcglincy_f5_3_bam$gc <- compute_rpf_gc(mcglincy_f5_3_bam, omit="APE",
                                       transcript_seq, transcript_lengths)
mcglincy_f5_3_training$gc <- compute_rpf_gc(mcglincy_f5_3_training, omit="APE",
                                            transcript_seq, transcript_lengths)

# compute regression ------------------------------------------------------

mcglincy_f5_3_gc_fit <- glm.nb(regression_model,
                               data=mcglincy_f5_3_training)
mcglincy_f5_3_gc_coef <- parse_coefs(mcglincy_f5_3_gc_fit)

save(mcglincy_f5_3_gc_fit, file="mcglincy_f5_3_gc_fit_400.Rda")
save(mcglincy_f5_3_gc_coef, file="mcglincy_f5_3_gc_coef_400.Rda")

# correct counts ----------------------------------------------------------

beta_gc <- coef(mcglincy_f5_3_gc_fit)
beta_gc <- beta_gc[names(beta_gc)=="gc"]
mean_gc <- with(mcglincy_f5_3_bam, sum(count*gc)/sum(count))

mcglincy_f5_3_bam$correct_400_gc <- correct_bias(mcglincy_f5_3_bam, mcglincy_f5_3_gc_fit)
mcglincy_f5_3_bam$correct_400_gc <- with(mcglincy_f5_3_bam,
                                         exp(log(correct_400_gc) + beta_gc*mean_gc -
                                               beta_gc*gc))
mcglincy_f5_3_bam$correct_400_gc <- with(mcglincy_f5_3_bam,
                                         correct_400_gc * sum(count) / sum(correct_400_gc, na.rm=T))

mcglincy_f5_3_training$correct_400_gc <- count_footprints(mcglincy_f5_3_bam,
                                                          mcglincy_f5_3_training,
                                                          "correct_400_gc")

save(mcglincy_f5_3_bam, file="mcglincy_f5_3_gc_bam.Rda")
save(mcglincy_f5_3_training, file="mcglincy_f5_3_gc_training_data.Rda")

# evaluate bias -----------------------------------------------------------

gc_codon_corr <- evaluate_bias(mcglincy_f5_3_training,
                               which_column="correct_400_gc",
                               transcript_fa_fname, transcript_length_fname,
                               type="codon")
gc_nt_corr <- evaluate_bias(mcglincy_f5_3_training,
                            which_column="correct_400_gc",
                            transcript_fa_fname, transcript_length_fname,
                            type="nt")

save(gc_codon_corr, gc_nt_corr, file="mcglincy_f5_3_gc_corr.Rda")
