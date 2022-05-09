rm(list=ls())

library(here)
library(Biostrings)

# load suppl. table 6 from Jansson et al (2021)
## "A single canonical transcript representing each protein-coding gene was
## selected from the GRCh38, v.97 Ensembl annotation file)"
suppl_tables_fname <- "41594_2021_669_MOESM3_ESM.xlsx"
suppl_table_6 <- readxl::read_xlsx(suppl_tables_fname,
                                   sheet="ST6 Selected transcripts list")

# load transcript sequences from UCSC Genome Browser
## GENCODE V39 Genomic Sequence
## 20nt upstream + 5' UTR Exons + CDS Exons + 3' UTR Exons + 20nt downstream
## One FASTA record per gene
## CDS in upper case, UTR in lower case
all_transcripts <- readBStringSet("hgTables.txt")
all_transcripts <- unique(all_transcripts)
names(all_transcripts) <- sub("^hg38_knownGene_", "",
                        sub("\\.[[:digit:]] .*$", "",
                            names(all_transcripts)))

# subset all_transcripts to transcripts listed in canonical_transcripts
canonical_transcripts <- all_transcripts[which(names(all_transcripts) %in%
                                                 suppl_table_6$Transcript_stable_ID)]

# compute 5' UTR, CDS, 3' UTR lengths
transcript_seq <- lapply(as.character(canonical_transcripts),
                         function(x) {
                           nt_seq <- strsplit(x, split="")[[1]]
                           cds_positions <- which(nt_seq %in% LETTERS)
                           if(length(cds_positions) == 0) { return(NA) }
                           cds_start <- min(cds_positions)
                           cds_end <- max(cds_positions)
                           return(c(utr5=substr(x, cds_start-20, cds_start-1),
                                    cds=substr(x, cds_start, cds_end),
                                    utr3=substr(x, cds_end+1, cds_end+20)))
                         })
no_cds <- which(is.na(transcript_seq))
canonical_transcripts <- canonical_transcripts[-no_cds]
transcript_seq <- transcript_seq[-no_cds]
transcript_lengths <- nchar(do.call(rbind, transcript_seq))

# write transcript sequences to file
transcript_seq <- sapply(transcript_seq, paste, collapse="")
transcript_seq <- toupper(transcript_seq)
transcript_seq <- DNAStringSet(transcript_seq)
writeXStringSet(transcript_seq, "grch38.transcripts.20cds20.fa")

# write transcript lengths to file
write.table(transcript_lengths, "grch38.transcripts.20cds20.lengths.tsv",
            quote=F, row.names=T, col.names=F, sep="\t")
