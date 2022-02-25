rm(list=ls())

library(here)
library(Biostrings)

ref_dir <- file.path(here(), "reference_data")

# sk1 transcripts
sk1_transcripts <- readDNAStringSet(file.path(ref_dir, "sk1.transcripts.20cds20.fa"))
sk1_cds <- lapply(sk1_transcripts, function(x) { subseq(x, 21, nchar(x)-20) })
sk1_aa <- lapply(sk1_cds, translate)
sk1_aa <- AAStringSet(sapply(sk1_aa, toString))
writeXStringSet(sk1_aa, file.path(ref_dir, "sk1.transcripts.aa.fa"))

# scer transcripts
scer_transcripts <- readDNAStringSet(file.path(ref_dir, "scer.transcripts.20cds20.fa"))
scer_cds <- lapply(scer_transcripts, function(x) { subseq(x, 21, nchar(x)-20 )})
scer_aa <- lapply(scer_cds, translate)
scer_aa <- AAStringSet(sapply(sk1_aa, toString))
writeXStringSet(sk1_aa, file.path(ref_dir, "scer.transcripts.aa.fa"))
