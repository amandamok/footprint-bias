rm(list=ls())

library(here)
library(choros)

# reference data
ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "sk1.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "sk1.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_length_fname)
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

expt_dir <- file.path(here(), "expts", "cheng_2018")

# choros parameters
f5_length <- 3
f3_length <- 3
num_genes <- 100
interxn_model <- formula(count ~ transcript*timepoint + A*timepoint + P*timepoint + E + d5*genome_f5 + d3*genome_f3)

timepoint_1 <- "timepoint_0hr"
timepoints <- paste0("timepoint_",
                     c(paste0(c(1.5, 3, 4.5, 6, 8, 10), "hr"),
                       "spore"))

for(timepoint_2 in timepoints) {
  print(timepoint_2)
  
  comparison <- paste(timepoint_1, timepoint_2, sep="_")
  results_dir <- file.path(expt_dir, "timecourse_comparison", comparison)
  if(!dir.exists(results_dir)) { dir.create(results_dir) }
  
  timepoint_1_dir <- file.path(expt_dir, "individual_timepoints", timepoint_1)
  timepoint_2_dir <- file.path(expt_dir, "individual_timepoints", timepoint_2)
  
  # 1. read in bam data
  print("... loading in bam data")
  load(file.path(timepoint_1_dir, paste0(timepoint_1, "_bam.Rda")))
  load(file.path(timepoint_2_dir, paste0(timepoint_2, "_bam.Rda")))
  
  # 2. choose d5/d3 subsets
  print("... choosing d5/d3 subsets")
  load(file.path(timepoint_1_dir, paste0(timepoint_1, "_subsets.Rda")))
  load(file.path(timepoint_2_dir, paste0(timepoint_2, "_subsets.Rda")))
  interxn_subsets <- unique(rbind(get(paste0(timepoint_1, "_subsets")),
                                  get(paste0(timepoint_2, "_subsets"))))
  save(interxn_subsets, file=file.path(results_dir, "interxn_subsets.Rda"))
  
  # 3. choose training set
  print("... choosing training set")
  assign(paste0(timepoint_1, "_density"),
         sort(calculate_transcript_density(get(paste0(timepoint_1, "_bam")),
                                           transcript_length_fname, statistic="median"),
              decreasing=T))
  assign(paste0(timepoint_2, "_density"),
         sort(calculate_transcript_density(get(paste0(timepoint_2, "_bam")),
                                           transcript_length_fname, statistic="median"),
              decreasing=T))
  ## choose genes in top (num_genes) for both data sets
  assign(paste0(timepoint_1, "_top_transcripts"),
         names(get(paste0(timepoint_1, "_density"))[1:num_genes]))
  assign(paste0(timepoint_2, "_top_transcripts"),
         names(get(paste0(timepoint_2, "_density"))[1:num_genes]))
  training_set <- get(paste0(timepoint_1, "_top_transcripts"))[get(paste0(timepoint_1, "_top_transcripts")) %in%
                                                                 get(paste0(timepoint_2, "_top_transcripts"))]
  ## sort remainder of genes by average density across data sets
  if(length(training_set) < num_genes) {
    top_transcripts <- unique(get(paste0(timepoint_1, "_top_transcripts")),
                              get(paste0(timepoint_2, "_top_transcripts")))
    remaining_transcripts <- top_transcripts[!(top_transcripts %in% training_set)]
    remaining_transcripts <- cbind(get(paste0(timepoint_1, "_density"))[match(remaining_transcripts, 
                                                                              names(get(paste0(timepoint_1, "_density"))))],
                                   get(paste0(timepoint_2, "_density"))[match(remaining_transcripts, 
                                                                              names(get(paste0(timepoint_2, "_density"))))])
    remaining_transcripts <- sort(rowMeans(remaining_transcripts), decreasing=T)
    training_set <- c(training_set, names(remaining_transcripts)[1:(num_genes-length(training_set))])
  }
  writeLines(training_set, 
             file.path(results_dir, paste0(comparison, "_training_set.txt")))
  
  # 4. initialize data frames for regression
  print("... generating training data")
  training_data <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets=interxn_subsets, f5_length=f5_length,
                             f3_length=f3_length, which_transcripts=training_set)
  assign(paste0(timepoint_1, "_training"),
         within(training_data, {
           count <- count_footprints(get(paste0(timepoint_1, "_bam")), 
                                     training_data, "count")
           timepoint <- timepoint_1
         }))
  assign(paste0(timepoint_2, "_training"),
         within(training_data, {
           count <- count_footprints(get(paste0(timepoint_2, "_bam")), 
                                     training_data, "count")
           timepoint <- timepoint_2
         }))
  training_data <- rbind(get(paste0(timepoint_1, "_training")),
                         get(paste0(timepoint_2, "_training")))
  save(training_data, file=file.path(results_dir, paste0(comparison, "_training_data.Rda")))
  
  # 5. compute regression fit
  print("... computing regression fit")
  interxn_fit <- MASS::glm.nb(interxn_model, training_data, model=F)
  save(interxn_fit, file=file.path(results_dir, paste0(comparison, "_fit.Rda")))
  
  rm(list=c(grep("_(bam|subsets|density|training)", ls(), value=T),
            "training_data", "interxn_fit"))
}
