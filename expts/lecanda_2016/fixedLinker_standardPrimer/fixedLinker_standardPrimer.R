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

min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 250

regression_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

paper <- "lecanda_2016"
sample <- "fixedLinker_standardPrimer"
results_dir <- file.path(here(), "expts", paper, sample)

# read in footprint alignments --------------------------------------------

raw_bam_fname <- file.path(here(), "expts", paper, "raw_data", sample,
                           paste0(sample, "_trim3_footprints.transcript.bam"))

bam_obj <- paste0(sample, "_bam")
bam_fname <- file.path(results_dir, paste0(bam_obj, ".Rda"))
if(!file.exists(bam_fname)) {
  assign(bam_obj,
         load_bam(raw_bam_fname, transcript_fa_fname, transcript_length_fname, 
                  offsets_fname, f5_length=f5_length, f3_length=f3_length))
  save(list=bam_obj, file=bam_fname)
} else {
  load(bam_fname)
}

# compute size/frame subsets ----------------------------------------------

d5_d3_obj <- paste0(sample, "_d5_d3")
subsets_obj <- paste0(sample, "_subsets")
subsets_fname <- file.path(results_dir, paste0(sample, "_subsets.Rda"))
if(!file.exists(subsets_fname)) {
  assign(d5_d3_obj, count_d5_d3(get(bam_obj)))
  assign(subsets_obj, 
         get(d5_d3_obj)$counts[1:(which(get(d5_d3_obj)$counts$proportion>min_prop)[1]),
                               c("d5", "d3")])
  save(list=c(d5_d3_obj, subsets_obj), file=subsets_fname)
} else {
  load(subsets_fname)
}

# establish training data -------------------------------------------------

# count footprints per transcript
transcript_counts <- aggregate(count ~ transcript, data=get(bam_obj), FUN=sum)
transcript_counts <- transcript_counts[order(transcript_counts$count, decreasing=T),]

training_set <- as.character(transcript_counts$transcript[1:num_genes])

# initialize training data for regression ---------------------------------

training_obj <- paste0(sample, "_training")
training_fname <- file.path(results_dir, paste0(training_obj, ".Rda"))
if(!file.exists(training_fname)) {
  assign(training_obj, 
         init_data(transcript_fa_fname, transcript_length_fname,
                   d5_d3_subsets=get(subsets_obj), f5_length=f5_length, f3_length=f3_length,
                   which_transcripts=training_set))
  assign(training_obj,
         within(get(training_obj), {
           transcript <- relevel(get(training_obj)$transcript, ref=training_set[1])
           count <- count_footprints(get(bam_obj), get(training_obj), "count")
         }))
  save(list=training_obj, file=training_fname)
} else {
  load(training_fname)
}

# compute regression ------------------------------------------------------

fit_obj <- paste0(sample, "_fit_", num_genes)
coef_obj <- paste0(sample, "_coef_", num_genes)
fit_fname <- file.path(results_dir, paste0(fit_obj, ".Rda"))
if(!file.exists(fit_fname)) {
  assign(fit_obj, MASS::glm.nb(regression_model, data=get(training_obj), model=F))
  assign(coef_obj, parse_coefs(get(fit_obj)))
  save(list=fit_obj, file=fit_fname)
  save(list=coef_obj, file=file.path(results_dir, paste0(coef_obj, ".Rda")))
} else {
  load(fit_fname)
}

# correct counts ----------------------------------------------------------

correction_name <- paste0("correct_", num_genes)
assign(bam_obj,
       within(get(bam_obj), {
         assign(correction_name,
                correct_bias(get(bam_obj), get(fit_obj)))
       }))
assign(training_obj,
       within(get(training_obj), {
         assign(correction_name,
                count_footprints(get(bam_obj), get(training_obj), correction_name))
       }))

save(list=bam_obj, file=bam_fname)
save(list=training_obj, file=training_fname)

# evaluate bias -----------------------------------------------------------

to_evaluate <- c("count", correction_name)

codon_corr_obj <- paste0(sample, "_codon_corr")
assign(codon_corr_obj,
       lapply(setNames(to_evaluate, to_evaluate),
              function(x) {
                evaluate_bias(get(training_obj), which_column=x,
                              transcript_fa_fname, transcript_length_fname,
                              type="codon")
              }))
nt_corr_obj <- paste0(sample, "_nt_corr")
assign(nt_corr_obj,
       lapply(setNames(to_evaluate, to_evaluate),
              function(x) {
                evaluate_bias(get(training_obj), which_column=x,
                              transcript_fa_fname, transcript_length_fname,
                              type="nt")
              }))

codon_max <- max(unlist(get(codon_corr_obj)))
nt_max <- max(unlist(get(nt_corr_obj)))

plot_obj <- paste0(sample, "_plot")
assign(plot_obj,
       (plot_bias(get(codon_corr_obj)$count) + 
          coord_cartesian(ylim=c(0, codon_max)) + 
          ggtitle(sample, subtitle="raw")) + 
         (plot_bias(get(codon_corr_obj)[[correction_name]]) + 
            coord_cartesian(ylim=c(0, codon_max)) + 
            ggtitle("", subtitle="corrected")) + 
         (plot_bias(get(nt_corr_obj)$count, type="nt") + 
            coord_cartesian(ylim=c(0, nt_max))) + 
         (plot_bias(get(nt_corr_obj)[[correction_name]], type="nt") + 
            coord_cartesian(ylim=c(0, nt_max))))

save(list=c(codon_corr_obj, nt_corr_obj, plot_obj),
     file=file.path(results_dir, paste0(sample, "_corr.Rda")))

q(save="no")
