library(here)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(ggplot2)
library(foreach)

num_cores <- 20
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

min_utr5_length <- 18 + 2
min_utr3_length <- 12 + 2
min_cds_length <- 10

ref_dir <- file.path(here(), "reference_data")

combine_list <- function(list_obj) {
  # list_obj: list of Bioconductor List objects
  list_obj <- unlist(list_obj)
  names(list_obj) <- NULL
  list_obj <- do.call(c, list_obj)
  return(list_obj)
}

get_unique_cds <- function(gene_cds) {
  # gene_cds: DNAStringSet; CDS sequences of transcripts for same gene
  gene_cds <- split(gene_cds, as.character(gene_cds))
  gene_cds <- lapply(gene_cds,
                     function(x) {
                       if(length(x) == 1) {
                         return(x)
                       } else {
                         transcripts <- names(x)
                         utr5_lengths <- nchar(utr5_seq[match(transcripts, names(utr5_seq))])
                         utr3_lengths <- nchar(utr3_seq[match(transcripts, names(utr3_seq))])
                         longest_utr5 <- which.max(utr5_lengths)
                         longest_utr3 <- which.max(utr3_lengths)
                         if(!all(utr3_lengths >= min_utr3_length) &
                            all(utr5_lengths >= min_utr5_length)) {
                           return(x[longest_utr3])
                         } else {
                           return(x[longest_utr5])
                         }
                       }
                     })
  return(combine_list(gene_cds))
}

# load gencode gff annotations
gff_fname <- file.path(ref_dir, "gencode.v36.annotation.gff3")
gff <- rtracklayer::import.gff3(gff_fname)
gff_coding <- gff[gff$gene_type == "protein_coding" & gff$type=="transcript"]
txdb <- makeTxDbFromGFF(gff_fname)
tx <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
mapping <- mcols(tx)
mapping$gene_id <- as.character(mapping$gene_id)
rownames(mapping) <- mapping$tx_id

# protein-cds transcripts ----------------------------------------------

# 1. pull CDS sequences
cds <- cdsBy(txdb, by="tx", use.names=T)
cds_seq <- extractTranscriptSeqs(Hsapiens, cds)
# filter for CDS length is multiple of 3
cds_seq <- cds_seq[(width(cds_seq) %% 3) == 0]
# filter for protein-coding transcripts
cds_seq <- cds_seq[names(cds_seq) %in% unique(gff_coding$transcript_id)]
# remove CDS with length < 10 codons
cds_seq <- cds_seq[width(cds_seq) > 3*min_cds_length]

# 2. pull 5' UTR sequences
# pull 5' UTR entries
utr5 <- fiveUTRsByTranscript(txdb, use.names=T)
utr5 <- GenomicFeatures::extractTranscriptSeqs(Hsapiens, utr5)
utr5 <- utr5[names(utr5) %in% names(cds_seq)]
# add entries for transcripts without 5' UTRs
utr5_seq <- as.character(utr5)
utr5_seq <- c(utr5_seq, rep("", length(cds_seq) - length(utr5)))
names(utr5_seq)[names(utr5_seq)==""] <- names(cds_seq)[!(names(cds_seq) %in% names(utr5))]

# 3. pull 3' UTR sequences
# pull 3' UTR entries
utr3 <- threeUTRsByTranscript(txdb, use.names=T)
utr3 <- GenomicFeatures::extractTranscriptSeqs(Hsapiens, utr3)
utr3 <- utr3[names(utr3) %in% names(cds_seq)]
# add entries for transcripts without 3' UTRs
utr3_seq <- as.character(utr3)
utr3_seq <- c(utr3_seq, rep("", length(cds_seq) - length(utr3)))
names(utr3_seq)[names(utr3_seq)==""] <- names(cds_seq)[!(names(cds_seq) %in% names(utr3))]

# 4. subset to unique CDS sequences
cds_seq_by_gene <- split(cds_seq, mapping$gene_id[match(names(cds_seq), mapping$tx_name)])
unique_cds <- foreach(chunk = split(cds_seq_by_gene,
                                    cut(seq_along(cds_seq_by_gene),
                                        num_cores))) %dopar% {
                                          return(lapply(chunk, get_unique_cds))
                                        }
unique_cds <- combine_list(unique_cds)
unique_cds_seq <- as.character(unique_cds)

# 4. pad out 5' UTRs of insufficient lengths
utr5_seq <- utr5_seq[names(utr5_seq) %in% names(unique_cds_seq)]
good_utr5_seq <- utr5_seq[nchar(utr5_seq) >= min_utr5_length]
bad_utr5 <- names(utr5_seq)[nchar(utr5_seq) < min_utr5_length]
bad_utr5_seq <- tx[tx$tx_name %in% bad_utr5,]
bad_utr5_seq <- extractUpstreamSeqs(Hsapiens, bad_utr5_seq, width=min_utr5_length)
names(bad_utr5_seq) <- bad_utr5
utr5_seq <- c(good_utr5_seq, as.character(bad_utr5_seq))

# 5. pad out 3' UTRs of insufficient lengths
utr3_seq <- utr3_seq[names(utr3_seq) %in% names(unique_cds_seq)]
good_utr3_seq <- utr3_seq[nchar(utr3_seq) >= min_utr3_length]
bad_utr3 <- names(utr3_seq)[nchar(utr3_seq) < min_utr3_length]
bad_utr3_seq <- tx[tx$tx_name %in% bad_utr3,]
# bad_utr3_seq <- extractUpstreamSeqs(Hsapiens, bad_utr3_seq, width=min_utr3_length)
bad_utr3_seq <- trim(suppressWarnings(flank(bad_utr3_seq, width=min_utr3_length, start=F)))
bad_utr3_seq <- getSeq(Hsapiens, bad_utr3_seq)
names(bad_utr3_seq) <- bad_utr3
utr3_seq <- c(good_utr3_seq, as.character(bad_utr3_seq))

# 6. write transcript sequences to file
# concatenate 5' UTR, CDS, and 3' UTR
transcript_seqs <- data.frame(utr5 = utr5_seq[match(names(unique_cds_seq), names(utr5_seq))],
                              cds = unique_cds_seq,
                              utr3 = utr3_seq[match(names(unique_cds_seq), names(utr3_seq))],
                              stringsAsFactors=F)
transcripts_fasta <- apply(transcript_seqs, 1, paste, collapse="")
transcripts_fasta <- DNAStringSet(transcripts_fasta)
# write to file
writeXStringSet(transcripts_fasta,
                filepath=file.path(ref_dir, "grch38.transcripts.fa"))

# 7. write transcript lengths to file
# extract region lengths
transcript_lengths <- data.frame(transcript_id = rownames(transcript_seqs),
                                 utr5_length = nchar(transcript_seqs$utr5),
                                 cds_length = nchar(transcript_seqs$cds),
                                 utr3_length = nchar(transcript_seqs$utr3))
# write to file
write.table(transcript_lengths,
            file=file.path(ref_dir, "grch38.transcripts.lengths.tsv"),
            quote=F, sep="\t", row.names=F, col.names=F)

# 8. write transcript to gene conversion
mapping <- data.frame(subset(mapping, tx_name %in% transcript_lengths$transcript_id))
mapping <- mapping[, -1]
mapping$gene_name <- gff$gene_name[match(mapping$gene_id, gff$gene_id)]
transcripts_withUTRs <- intersect(names(good_utr5_seq), names(good_utr3_seq))
mapping$with_UTRs <- mapping$tx_name %in% transcripts_withUTRs
write.table(mapping,
            file=file.path(ref_dir, "grch38.transcripts.mappings.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T)

# 9. compute codon composition (counts) per transcript
codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")
transcript_codons <- sapply(as.character(transcript_seqs$cds),
                            function(x) {
                              # split CDS into codons
                              as_codon <- strsplit(x, split="")[[1]]
                              as_codon <- paste0(as_codon[c(T, F, F)],
                                                 as_codon[c(F, T, F)],
                                                 as_codon[c(F, F, T)])
                              # compute codon counts
                              sapply(codons, function(y) sum(as_codon==y))
                            })
transcript_codons <- data.frame(t(transcript_codons),
                                row.names=rownames(transcript_seq))
write.table(transcript_codons,
            file=file.path(ref_dir, "grch38.transcripts.codonCounts.tsv"),
            quote=F, sep="\t", row.names=T, col.names=T)

# # ncRNA -------------------------------------------------------------------
#
# # extract RNA entries
# ncRNA_types <- c("ribozyme", grep("RNA", unique(gff$gene_type), value=T))
# ncRNA_types <- ncRNA_types[ncRNA_types != "lncRNA"]
# ncRNA <- gff[gff$gene_type %in% ncRNA_types &
#                !grepl("PAR", gff$ID) &
#                gff$type == "transcript"]
#
# # pull sequences
# ncRNA_seq <- split(ncRNA, ncRNA$transcript_id)
# ncRNA_seq <- GenomicFeatures::extractTranscriptSeqs(Hsapiens, ncRNA_seq)
#
# # write to file
# writeXStringSet(ncRNA_seq,
#                 filepath=file.path(here(), "reference_data", "grch38.ncRNA.fa"))
#
# writeXStringSet(ncRNA_seq,
#                 filepath=file.path(here(), "reference_data", "grch38.tmp.fa"))

# close out ---------------------------------------------------------------

stopCluster(cl)
