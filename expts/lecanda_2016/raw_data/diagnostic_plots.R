rm(list=ls())

library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
scer_length_fname <- grep("scer.*lengths", list.files(ref_dir, full.names=T), value=T)

raw_data_dir <- file.path(here(), "expts", "lecanda_2016", "raw_data")
expts <- list.dirs(raw_data_dir, full.names=F, recursive=F)

for(expt in expts) {
  tmp_bam_fname <- file.path(raw_data_dir, expt,
                             paste0(expt, "_trim3_footprints.transcript.bam"))
  tmp_plot <- plot_diagnostic(tmp_bam_fname, scer_length_fname)
  save(tmp_plot, file=file.path(raw_data_dir, paste0(expt, "_diagnosticPlot.Rda")))
}

q(save="no")
