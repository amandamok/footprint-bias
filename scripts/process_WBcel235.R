rm(list=ls())

library(here)
library(Biostrings)
library(foreach)

num_cores <- 20

ref_dir <- file.path(here(), "reference_data")

utr5_length <- 20
utr3_length <- 20

# load sequences + annotations --------------------------------------------

# genome sequence and GFF from NCBI
# https://www.ncbi.nlm.nih.gov/genome/41

# load genome sequences
genome_seq_fname <- file.path(ref_dir,
                              "GCF_000002985.6_WBcel235_genomic.fna")
genome_seq <- readDNAStringSet(genome_seq_fname)
genome_seq <- as.character(genome_seq)
names(genome_seq) <- sub(" .*", "", names(genome_seq))

# load genome gff
genome_gff_fname <- file.path(ref_dir,
                              "GCF_000002985.6_WBcel235_genomic.gff")
genome_gff <- read.table(genome_gff_fname, sep="\t", stringsAsFactors=F)
colnames(genome_gff) <- c("chr", "source", "feature", "start", "end",
                          "score", "strand", "frame", "attribute")

# function ----------------------------------------------------------------

reverse_transcribe <- function(x) {
  # x: character; sequence to be reverse transcribed
  x <- strsplit(x, split="")[[1]]
  rev_x <- rev(x)
  rt <- c("A", "T", "C", "G")
  names(rt) <- c("T", "A", "G", "C")
  rt_x <- rt[rev_x]
  return(paste(rt_x, collapse=""))
}

extract_sequences <- function(transcript_name, genome_gff,
                              utr5_length=20, utr3_length=20) {
  # transcript_name: character; transcript name
  # genome_gff: data.frame; corresponds to GFF annotation
  transcript_gff <- subset(genome_gff,
                           grepl(transcript_name, genome_gff$attribute,
                                 fixed=T))
  transcript_mRNA <- subset(transcript_gff, feature=="mRNA")
  transcript_strand <- transcript_mRNA$strand
  # 1. extract CDS sequences
  transcript_cds <- subset(transcript_gff, feature=="CDS")
  transcript_cds$sequence <- mapply(substr,
                                    genome_seq[transcript_cds$chr],
                                    transcript_cds$start,
                                    transcript_cds$end)
  if(transcript_strand == "+") {
    transcript_cds <- transcript_cds[order(transcript_cds$start, decreasing=F),]
  } else {
    transcript_cds <- transcript_cds[order(transcript_cds$start, decreasing=T),]
    transcript_cds$sequence <- sapply(transcript_cds$sequence, reverse_transcribe)
  }
  cds_seq <- paste(transcript_cds$sequence, collapse="")
  # 2. extract 5' UTR sequence
  if(transcript_strand == "+") {
    utr5_end <- transcript_cds$start[1] - 1
    utr5_start <- utr5_end - utr5_length + 1
    utr5_seq <- substr(genome_seq[transcript_mRNA$chr], utr5_start, utr5_end)
  } else { # (-) strand
    utr5_start <- transcript_cds$end[1] + 1
    utr5_end <- utr5_start + utr5_length - 1
    utr5_seq <- substr(genome_seq[transcript_mRNA$chr], utr5_start, utr5_end)
    utr5_seq <- reverse_transcribe(utr5_seq)
  }
  # 3. extract 3' UTR sequence
  if(transcript_strand == "+") {
    utr3_start <- transcript_cds$end[nrow(transcript_cds)] + 1
    utr3_end <- utr3_start + utr3_length - 1
    utr3_seq <- substr(genome_seq[transcript_mRNA$chr], utr3_start, utr3_end)
  } else { # (-) strand
    utr3_end <- transcript_cds$start[nrow(transcript_cds)] - 1
    utr3_start <- utr3_end - utr3_length + 1
    utr3_seq <- substr(genome_seq[transcript_mRNA$chr], utr3_start, utr3_end)
    utr3_seq <- reverse_transcribe(utr3_seq)
  }
  return(c("utr5"=utr5_seq, "cds"=cds_seq, "utr3"=utr3_seq))
}

# pull (spliced) transcripts ----------------------------------------------

gff_mRNA <- subset(genome_gff, feature=="mRNA")
gff_mRNA$transcript <- sapply(gff_mRNA$attribute,
                              function(x) {
                                sub("ID=rna-", "",
                                    strsplit(x, split=";")[[1]][1])
                              })

cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)
transcript_seq <- foreach(x=gff_mRNA$transcript, .combine='rbind') %dopar% {
  extract_sequences(x, genome_gff)
}
stopCluster(cl)
colnames(transcript_seq) <- c("utr5", "cds", "utr3")

# filter out non-unique CDS -----------------------------------------------

# for non-unique CDS: choose longest 5' UTR (else, 3' UTR)
which_duplicated <- which(duplicated(transcript_seq[, "cds"]))
remove_transcripts <- sapply(which_duplicated,
                             function(x) {
                               which_duplicated <- which(transcript_seq[, "cds"] ==
                                                           transcript_seq[x, "cds"])
                               tmp_gff <- gff_mRNA[which_duplicated,]
                               tmp_strand <- unique(tmp_gff$strand)
                               if(length(tmp_strand) > 1) { # duplicated CDS on different strands
                                 return(NA)
                               }
                               if(length(unique(tmp_gff$start)) > 1) { # pick longest 5' UTR
                                 if(tmp_strand == "+") {
                                   which_remove <- which_duplicated[which.min(tmp_gff$start)]
                                 } else {
                                   which_remove <- which_duplicated[which.max(tmp_gff$end)]
                                 }
                               } else { # pick longest 3' UTR
                                 if(tmp_strand == "+") {
                                   which_remove <- which_duplicated[which.max(tmp_gff$end)]
                                 } else {
                                   which_remove <- which_duplicated[which.min(tmp_gff$start)]
                                 }
                               }
                               return(which_remove)
                             })
remove_transcripts <- remove_transcripts[!is.na(remove_transcripts)]
transcript_seq <- transcript_seq[-remove_transcripts,]
gff_mRNA <- gff_mRNA[-remove_transcripts,]

# compute transcript lengths ----------------------------------------------

transcript_lengths <- nchar(transcript_seq)
transcript_lengths <- data.frame(gff_mRNA$transcript, transcript_lengths)
colnames(transcript_lengths) <- c("transcript", "utr5", "cds", "utr3")

# output transcriptome ----------------------------------------------------

transcript_fasta <- apply(transcript_seq, 1, paste, collapse="")
names(transcript_fasta) <- gff_mRNA$transcript
transcript_fasta <- DNAStringSet(transcript_fasta)
writeXStringSet(transcript_fasta, file.path(ref_dir, "WBcel235.transcripts.fa"))

write.table(transcript_lengths, row.names=F, col.names=F, quote=F, sep="\t",
            file=file.path(ref_dir, "WBcel235.transcripts.lengths.txt"))

# extract contaminant RNAs ------------------------------------------------

contaminant_types <- c("lnc_RNA", "ncRNA", "pseudogenic_tRNA",
                       "rRNA", "snoRNA", "snRNA", "tRNA")
contaminants <- subset(genome_gff, feature %in% contaminant_types)
contaminants$id <- sapply(contaminants$attribute,
                          function(x) {
                            sub("ID=rna-", "",
                                strsplit(x, split=";")[[1]][1])
                          })
contaminants$id <- sub(":.*$", "", contaminants$id)
contaminants_seq <- mapply(substr,
                           genome_seq[contaminants$chr],
                           contaminants$start,
                           contaminants$end)
contaminants_seq[contaminants$strand=="-"] <- reverse_transcribe(contaminants_seq[contaminants$strand=="-"])
names(contaminants_seq) <- contaminants$id
contaminants_seq <- DNAStringSet(contaminants_seq)
writeXStringSet(contaminants_seq, file.path(ref_dir, "WBcel235.contaminants.fa"))
