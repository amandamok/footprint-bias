rm(list=ls())

library(here)
library(Biostrings)
library(IRanges)

ref_dir <- file.path(here(), "reference_data")
genome_fasta_fname <- file.path(ref_dir, "sk1_MvO_V1.fasta")
gff_fname <- file.path(ref_dir, "SK1_NCSL00000000_SGD.gene.MvO.curated.mRNA.20190128.gff")

utr5_length <- 20
utr3_length <- 20

# load genome fasta -------------------------------------------------------

genome_fasta <- Biostrings::readDNAStringSet(genome_fasta_fname)

# load SK1 gff file -------------------------------------------------------

gff <- read.table(gff_fname, stringsAsFactors=F)
colnames(gff) <- c("chr", "source", "type", "start", "end", "score", "strand", 
                   "phase", "attributes")
gff$ID <- sapply(gff$attributes,
                 function(x) {
                   sub("ID=", "", strsplit(x, split=";")[[1]][1])
                 })
gff$parent <- sapply(gff$attributes,
                     function(x) {
                       sub("Parent=", "", strsplit(x, split=";")[[1]][2])
                     })
gff$gene <- sapply(gff$attributes,
                   function(x) {
                     sub("gene=", "", strsplit(x, split=";")[[1]][3])
                   })
gff$width <- with(gff, end-start+1)

# extract mRNA entires from gff
mRNA <- subset(gff, type=="mRNA")

# remove entries too close to chromosome edges
too_close_to_chr_edge <- c(subset(mRNA, start-15 < 1)$parent,
                           subset(mRNA, end+15 > width(genome_fasta[mRNA$chr]))$parent)
## YML132W, YDL248W
mRNA <- subset(mRNA, !(parent %in% too_close_to_chr_edge))
gff <- subset(gff, !(parent %in% too_close_to_chr_edge))

# remove duplicate mRNA ranges
mRNA_ranges <- with(mRNA, paste(chr, start, end, sep="_"))
mRNA_ranges[duplicated(mRNA_ranges)]
## [1] "chr04_713726_714148" "chr10_548149_549333"
## subset(gff, parent %in% subset(mRNA, chr=="chr04" & start==713726 & end==714148)$parent)
## --> single exon genes: YDR133C and YDR134C --> keep YDR133C
## subset(gff, parent %in% subset(mRNA, chr=="chr10" & start==548149 & end==549333)$parent)
## --> single exon genes: YJR080C and YJR079W --> keep YJR080C
mRNA <- subset(mRNA, !(parent %in% c("YDR133C", "YJR080C")))
gff <- subset(gff, !(parent %in% c("YDR133C", "YJR080C")))

# split gff file by gene
gff_by_gene <- split(gff, gff$gene)
num_exons_by_gene <- sapply(gff_by_gene, function(x) sum(x$type=="CDS"))

# pull only CDS entries in gff_by_gene, label exons by order
cds_by_gene <- lapply(gff_by_gene, 
                      function(x) {
                        tmp_strand <- subset(x, type=="mRNA")$strand
                        tmp_cds <- subset(x, type=="CDS") 
                        if(tmp_strand == "+") {
                          tmp_cds$exon_order <- order(tmp_cds$start)
                        } else {
                          tmp_cds$exon_order <- rev(order(tmp_cds$start))
                        }
                        return(tmp_cds)
                      })

# sanity checks -----------------------------------------------------------

# check that mRNA start+end correspond with CDS/exon start+end
## mRNA start (all ok)
summary(sapply(gff_by_gene, 
               function(x) {
                 tmp_strand <- subset(x, type=="mRNA")$strand
                 mRNA_start <- ifelse(tmp_strand=="+",
                                      subset(x, type=="mRNA")$start,
                                      subset(x, type=="mRNA")$end)
                 CDS_start <- ifelse(tmp_strand=="+",
                                     min(subset(x, type=="CDS")$start),
                                     max(subset(x, type=="CDS")$end))
                 return(mRNA_start == CDS_start)
               }))
## mRNA end (all ok)
summary(sapply(gff_by_gene, 
               function(x) {
                 tmp_strand <- subset(x, type=="mRNA")$strand
                 mRNA_end <- ifelse(tmp_strand=="+",
                                    subset(x, type=="mRNA")$end,
                                    subset(x, type=="mRNA")$start)
                 CDS_end <- ifelse(tmp_strand=="+",
                                   max(subset(x, type=="CDS")$end),
                                   min(subset(x, type=="CDS")$start))
                 return(mRNA_end == CDS_end)
               }))

# check that mRNA entries start with ATG codon
summary(as.factor(sapply(seq(nrow(mRNA)),
                         function(x) {
                           tmp_strand <- mRNA$strand[x]
                           if(tmp_strand == "+") {
                             start_codon <- subseq(genome_fasta[mRNA$chr[x]],
                                                   mRNA$start[x],
                                                   mRNA$start[x]+2)
                           } else {
                             start_codon <- subseq(genome_fasta[mRNA$chr[x]],
                                                   mRNA$end[x]-2,
                                                   mRNA$end[x])
                             start_codon <- reverseComplement(start_codon)
                           }
                           return(as.character(start_codon))
                         })))
## AGC  AGT  ATA  ATG  ATT 
## 1    2    5 6206    1 

# check that CDS widths are multiple of 3
summary(as.factor(sapply(cds_by_gene, function(x) (sum(x$width)%%3))))
## n=73 with 1 extra nt
## n=102 with 1 missing nt

# extract 5' UTR sequence -------------------------------------------------

mRNA$utr5 <- sapply(seq(nrow(mRNA)),
                    function(x) {
                      tmp_strand <- mRNA$strand[x]
                      if(tmp_strand == "+") {
                        utr5 <- subseq(genome_fasta[mRNA$chr[x]],
                                       mRNA$start[x]-utr5_length,
                                       mRNA$start[x]-1)
                      } else {
                        utr5 <- subseq(genome_fasta[mRNA$chr[x]],
                                       mRNA$end[x]+1,
                                       mRNA$end[x]+utr5_length)
                        utr5 <- reverseComplement(utr5)
                      }
                      return(as.character(utr5))
                    })

# extract CDS sequence ----------------------------------------------------

mRNA$CDS <- sapply(mRNA$gene,
                   function(x) {
                     tmp_cds <- cds_by_gene[[x]]
                     tmp_strand <- tmp_cds$strand[1]
                     tmp_cds$exon_seq <- sapply(seq(nrow(tmp_cds)),
                                                function(y) {
                                                  tmp_seq <- subseq(genome_fasta[tmp_cds$chr[y]],
                                                                    tmp_cds$start[y],
                                                                    tmp_cds$end[y])
                                                  if(tmp_strand=="-") {
                                                    tmp_seq <- reverseComplement(tmp_seq)
                                                  }
                                                  return(as.character(tmp_seq))
                                                })
                     cds_seq <- tmp_cds$exon_seq[tmp_cds$exon_order]
                     return(paste(cds_seq, collapse=""))
                   })

# sanity check start codon
summary(as.factor(sapply(mRNA$CDS, function(x) substr(x, 1, 3))))
## ATA  ATG  ATT 
## 5 6208    2 

# sanity check CDS length
summary(as.factor(nchar(mRNA$CDS) %% 3))
## 0    1    2 
## 6040   73  102 

# extract 3' UTR sequence -------------------------------------------------

mRNA$utr3 <- sapply(seq(nrow(mRNA)),
                    function(x) {
                      tmp_strand <- mRNA$strand[x]
                      if(tmp_strand == "+") {
                        utr3 <- subseq(genome_fasta[mRNA$chr[x]],
                                       mRNA$end[x]+1,
                                       mRNA$end[x]+utr3_length)
                      } else {
                        utr3 <- subseq(genome_fasta[mRNA$chr[x]],
                                       mRNA$start[x]-utr3_length,
                                       mRNA$start[x]-1)
                        utr3 <- reverseComplement(utr3)
                      }
                      return(as.character(utr3))
                    })

# write transcriptome fasta and lengths to file ---------------------------

# transcriptome FASTA file
SK1_transcriptome <- with(mRNA, paste(utr5, CDS, utr3, sep=""))
names(SK1_transcriptome) <- mRNA$parent
SK1_transcriptome <- DNAStringSet(SK1_transcriptome)
writeXStringSet(SK1_transcriptome,
                file.path(ref_dir, "sk1.transcripts.20cds20.fa"))

# lengths file
SK1_lengths <- data.frame(mRNA$parent,
                          nchar(mRNA$utr5),
                          nchar(mRNA$CDS),
                          nchar(mRNA$utr3))
write.table(SK1_lengths, 
            file=file.path(ref_dir, "sk1.transcripts.20cds20.lengths.txt"),
            quote=F, row.names=F, col.names=F, sep="\t")
