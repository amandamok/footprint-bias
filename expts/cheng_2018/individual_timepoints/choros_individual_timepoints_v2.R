rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)
library(choros)

# reference data
ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "sk1.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "sk1.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_length_fname)
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

# choros parameters
min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 250
regression_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)

# experiment variables
expt_dir <- file.path(here(), "expts", "cheng_2018")
# timepoints <- grep("txt", grep("timepoint", list.files(file.path(expt_dir, "raw_data")), value=T),
#                    value=T, invert=T)
timepoints <- paste0("timepoint_", c("8hr", "MATaa", "spore", "vegetative"))

for(timepoint in timepoints) {
  print(timepoint)
  
  results_dir <- file.path(expt_dir, "individual_timepoints", timepoint)
  if(!dir.exists(results_dir)) { dir.create(results_dir) }
  
  # 1. read in footprint alignments
  print("Loading bam file")
  bam_fname <- file.path(expt_dir, "raw_data", timepoint, 
                         paste0(timepoint, "_trim3_footprints.transcript.bam"))
  bam_obj <- paste0(timepoint, "_bam")
  assign(bam_obj,
         load_bam(bam_fname, transcript_fa_fname, transcript_length_fname,
                  offsets_fname, f5_length=f5_length, f3_length=f3_length))
  save(list=bam_obj, file=file.path(results_dir, paste0(bam_obj, ".Rda")))
  
  # 2. compute size/frame subsets
  print("Computing d5/d3 subsets")
  d5_d3_obj <- paste0(timepoint, "_d5_d3")
  assign(d5_d3_obj,
         count_d5_d3(get(bam_obj)))
  which_min_prop <- which(get(d5_d3_obj)$counts$proportion > min_prop)[1]
  subsets_obj <- paste0(timepoint, "_subsets")
  assign(subsets_obj,
         get(d5_d3_obj)$counts[1:which_min_prop, c("d5", "d3")])
  subset_names_obj <- paste0(timepoint, "_subset_names")
  assign(subset_names_obj,
         sapply(seq(nrow(get(subsets_obj))),
                function(x) {
                  paste("d5", get(subsets_obj)$d5[x],
                        "d3", get(subsets_obj)$d3[x], sep="_")
                }))
  save(list=c(d5_d3_obj, subsets_obj, subset_names_obj), 
       file=file.path(results_dir, paste0(timepoint, "_subsets.Rda")))
  
  # 3. establish training set
  ## count footprints per transcript
  print("Establishing training set")
  transcript_counts <- aggregate(count ~ transcript, data=get(bam_obj), FUN=sum)
  transcript_counts$length_aa <- transcript_lengths$cds_length[match(transcript_counts$transcript,
                                                                     transcript_lengths$transcript)]/ 3
  transcript_counts$TE <- with(transcript_counts, count / length_aa)
  transcript_counts <- transcript_counts[order(transcript_counts$TE, decreasing=T),]
  training_set <- as.character(transcript_counts$transcript)[1:num_genes]
  
  # 4. initialize data frames for regression
  print("Initializing training data")
  training_obj <- paste0(timepoint, "_training")
  assign(training_obj, 
         init_data(transcript_fa_fname, transcript_length_fname,
                   d5_d3_subsets=get(subsets_obj), f5_length=f5_length, f3_length=f3_length,
                   which_transcripts=training_set))
  assign(training_obj,
         within(get(training_obj), {
           transcript <- relevel(transcript, ref=training_set[1])
           count <- count_footprints(get(bam_obj), get(training_obj), "count")
         }))
  save(list=training_obj, file=file.path(results_dir, paste0(training_obj, ".Rda")))
  
  # 5. compute regression
  print("Computing regression")
  fit_obj <- paste0(timepoint, "_fit")
  coef_obj <- paste0(timepoint, "_coef")
  assign(fit_obj, MASS::glm.nb(regression_model, data=get(training_obj), model=F))
  save(list=fit_obj, file=file.path(results_dir, paste0(fit_obj, ".Rda")) )
  assign(coef_obj, parse_coefs(get(fit_obj)))
  save(list=coef_obj, file=file.path(results_dir, paste0(coef_obj, ".Rda")))
  
  # 6. correct counts
  print("Correcting counts")
  correct_name <- paste0("correct_", num_genes)
  assign(bam_obj,
         within(get(bam_obj), {
           assign(correct_name,
                  correct_bias(get(bam_obj), get(fit_obj)))
         }))
  assign(training_obj,
         within(get(training_obj), {
           assign(correct_name,
                  count_footprints(get(bam_obj), get(training_obj), "correct_250"))
         }))
  save(list=bam_obj, file=file.path(results_dir, paste0(bam_obj, ".Rda")))
  save(list=training_obj, file=file.path(results_dir, paste0(training_obj, ".Rda")))
  
  # 7. evaluate bias
  print("Evaluating bias")
  comparisons <- c("count", paste0("correct_", num_genes))
  ## codon correlations
  codon_corr_obj <- paste0(timepoint, "_codon_corr")
  assign(codon_corr_obj,
         within(list(), {
           for(x in comparisons) {
             assign(x, 
                    evaluate_bias(get(training_obj), which_column=x,
                                  transcript_fa_fname, transcript_length_fname,
                                  type="codon"))
             rm(x)
           }}))
  save(list=codon_corr_obj, file=file.path(results_dir, paste0(codon_corr_obj, ".Rda")))
  ## nt correlations
  nt_corr_obj <- paste0(timepoint, "_nt_corr")
  assign(nt_corr_obj,
         within(list(), {
           for(x in comparisons) {
             assign(x, 
                    evaluate_bias(get(training_obj), which_column=x,
                                  transcript_fa_fname, transcript_length_fname,
                                  type="nt"))
             rm(x)
           }}))
  save(list=nt_corr_obj, file=file.path(results_dir, paste0(nt_corr_obj, ".Rda")))
  
  # 8. make plots
  print("Generating plots")
  codon_plot_obj <- paste0(timepoint, "_codon_plots")
  assign(codon_plot_obj,
         within(list(), {
           for(x in comparisons) {
             assign(x, plot_bias(get(codon_corr_obj)[[x]], type="codon"))
             rm(x)
           }
         }))
  nt_plot_obj <- paste0(timepoint, "_nt_plots")
  assign(nt_plot_obj,
         within(list(), {
           for(x in comparisons) {
             assign(x, plot_bias(get(nt_corr_obj)[[x]], type="nt"))
             rm(x)
           }
         }))
  save(list=c(codon_plot_obj, nt_plot_obj),
       file=file.path(results_dir, paste0(timepoint, "_plots.Rda")))
  
  # 9. clean environment
  rm(list=grep(timepoint, ls(), value=T))
}

q(save="no")

