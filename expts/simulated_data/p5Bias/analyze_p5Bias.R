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

expt <- "p5Bias"
results_dir <- file.path(here(), "expts", "simulated_data", expt)

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- c(100, 75, 50, 25)
exclude_codon5 <- 10
exclude_codon3 <- 10
minimum_cds_length <- exclude_codon5 + exclude_codon3 + 10

orig_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

# 1. load footprint alignments --------------------------------------------

bam_fname <- file.path(results_dir, paste0(expt, ".transcript.bam"))
bam_data_fname <- file.path(results_dir, paste0(expt, "_bam.Rda"))
if(!file.exists(bam_data_fname)) {
  bam_data <- load_bam(bam_fname, transcript_fa_fname, transcript_length_fname,
                       offsets_fname, f5_length=f5_length, f3_length=f3_length)
  save(bam_data, file=bam_data_fname)
} else {
  load(bam_data_fname)
}

# 2. compute size/frame subsets -------------------------------------------

subsets_fname <- file.path(results_dir, paste0(expt, "_subsets.Rda"))
if(!file.exists(subsets_fname)) {
  d5_d3 <- count_d5_d3(bam_data)
  d5_d3_subsets <- d5_d3$counts[1:(which(d5_d3$counts$proportion>min_prop)[1]),
                                c("d5", "d3")]
  subset_names <- sapply(seq(nrow(d5_d3_subsets)),
                         function(x) {
                           paste("d5", d5_d3_subsets$d5[x],
                                 "d3", d5_d3_subsets$d3[x],
                                 sep="_")
                         })
  save(d5_d3, d5_d3_subsets, subset_names,
       file=subsets_fname)
} else {
  load(subsets_fname)
}

# 3. pick transcripts for training ----------------------------------------

training_set_fname <- file.path(results_dir, "training_transcripts.txt")
if(!file.exists(training_set_fname)) {
  # calculate read count across transcript
  training_set <- aggregate(count ~ transcript, data=bam_data, FUN=sum)
  # calculate median codon density across transcript
  codon_density <- calculate_transcript_density(bam_data,
                                                transcript_length_fname,
                                                statistic="median")
  training_set$codon_density <- codon_density[match(training_set$transcript,
                                                    names(codon_density))]
  # filter out transcripts with 0 counts
  training_set <- subset(training_set, count > 0)
  # filter for transcripts with sufficient number of codons
  training_set$num_codons <- transcript_lengths$cds_length[match(training_set$transcript,
                                                                 transcript_lengths$transcript)] / 3
  training_set <- subset(training_set, num_codons > minimum_cds_length)
  # pick top transcripts by median codon density
  training_set <- training_set[order(training_set$codon_density, decreasing=T),]
  training_set <- as.character(training_set$transcript)[1:max(num_genes)]
  writeLines(training_set, training_set_fname)
} else {
  training_set <- readLines(training_set_fname)
}

# 4. initialize data for regression ---------------------------------------

training_data_fname <- file.path(results_dir, "training_data.Rda")
if(!file.exists(training_data_fname)) {
  # initialize training data
  training_data <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets=d5_d3_subsets,,
                             f5_length=f5_length, f3_length=f3_length,
                             which_transcripts=training_set)
  training_data$transcript <- relevel(training_data$transcript, ref=training_set[1])
  # count footprints
  training_data$count <- count_footprints(bam_data, training_data, "count")
  save(training_data, file=training_data_fname)
} else {
  load(training_data_fname)
}

# 5. compute model fit and correct bias -----------------------------------

for(x in num_genes) {
  print(paste(x, "genes"))
  # calculate model fit
  fit_name <- paste0("model_fit_", x)
  fit_fname <- file.path(results_dir, paste0("model_fit_", x, "genes.Rda"))
  if(!file.exists(fit_fname)) {
    print("... calculating regression fit")
    assign(fit_name,
           MASS::glm.nb(orig_model,
                        data=subset(training_data,
                                    transcript %in% training_set[1:x]),
                        model=F))
    save(list=fit_name, file=fit_fname)
  } else {
    load(fit_fname)
  }
  # correct read counts
  correction_name <- paste0("correct_", x)
  if(!(correction_name %in% colnames(bam_data))) {
    print("... correcting read counts")
    bam_data[correction_name] <- correct_bias(bam_data, get(fit_name))
    save(bam_data, file=bam_data_fname)
  }
  if(!(correction_name %in% colnames(training_data))) {
    print("... calculating corrected counts")
    training_data[correction_name] <- count_footprints(bam_data,
                                                       training_data,
                                                       correction_name)
    save(training_data, file=training_data_fname)
  }
}

# 6. evaluate bias --------------------------------------------------------

bias_codon_fname <- file.path(results_dir, "bias_codon.Rda")
if(!file.exists(bias_codon_fname)) {
  bias_codon <- lapply(c("count", paste0("correct_", num_genes)),
                       function(x) {
                         evaluate_bias(training_data, which_column=x,
                                       transcript_fa_fname, transcript_length_fname,
                                       type="codon")
                       })
  save(bias_codon, file=bias_codon_fname)
}

bias_nt_fname <- file.path(results_dir, "bias_nt.Rda")
if(!file.exists(bias_nt_fname)) {
  bias_nt <- lapply(c("count", paste0("correct_", num_genes)),
                    function(x) {
                      evaluate_bias(training_data, which_column=x,
                                    transcript_fa_fname, transcript_length_fname,
                                    type="nt")
                    })
  save(bias_nt, file=bias_nt_fname)
}

q(save="no")

