rm(list=ls())

library(choros)
library(here)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")

ottrUMI_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "WT_UMI_OTTR_repB",
                                "ottrUMI_trimmed_deduplicated_footprints.transcript.bam")

ottrUMI_repB <- choros(ottrUMI_bam_fname,
                           transcript_fa_fname,
                           transcript_length_fname,
                           offsets_fname)

save(ottrUMI_repB, file="ottrUMI_repB.Rda")


q(save="no")
