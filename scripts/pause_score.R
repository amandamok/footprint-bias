rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

# define functions --------------------------------------------------------

get_A_codon <- function(transcript_name, cod_idx, transcript_seq, utr5_length=20) {
  # return A site codon corresponding to (transcript_name, cod_idx)
  ## transcript_name: character
  ## cod_idx: integer
  ## transcript_seq: list of character vectors; transcript (5'UTR, CDS, and 3'UTR) sequences
  ## utr5_length: integer; length of 5'UTR (assumed same for all transcripts)
  A_start <- utr5_length + 3*(cod_idx-1) + 1
  A_end <- utr5_length + 3*cod_idx
  A_codon <- substr(transcript_seq[transcript_name], A_start, A_end)
  return(A_codon)
}

load_fa <- function(transcript_fa_fname) {
  # load transcript sequences from genome .fa file
  ## transcripts_fa_fname: character; file path to transcriptome .fa file
  raw_text <- readLines(transcript_fa_fname)
  transcript_startLines <- grep(">", raw_text)
  num_transcripts <- length(transcript_startLines)
  transcript_names <- sapply(transcript_startLines,
                             function(x) {
                               gsub(">", "", strsplit(raw_text[x], split=" ")[[1]][1])
                             })
  transcript_startLines <- c(transcript_startLines, length(raw_text)+1) # add extra line for bookkeeping
  transcript_sequences <- sapply(1:num_transcripts,
                                 function(x) {
                                   startLine <- transcript_startLines[x]+1
                                   endLine <- transcript_startLines[x+1]-1
                                   transcriptSequence <- paste(raw_text[startLine:endLine], collapse="")
                                   return(transcriptSequence)
                                 })
  names(transcript_sequences) <- transcript_names
  return(transcript_sequences)
}

# load reference data -----------------------------------------------------

ref_dir <- file.path(here(), "reference_data")

transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_seq <- load_fa(transcript_fa_fname)

transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_length_fname)

codon_properties_fname <- file.path(ref_dir, "yeast_codon_properties.txt")
codon_properties <- read.delim(codon_properties_fname, stringsAsFactors=F)

# process data ------------------------------------------------------------

# 1. load bam files w/ raw & corrected counts
expts <- c("schuller_2017", "tunney_2018", "weinberg_2016")
for(expt in expts) {
  load(file.path(here(), "expts", expt,
                 paste0(strsplit(expt, split="_")[[1]][1], "_bam.Rda")))
}

groups <- expand.grid(expt=sapply(expts, function(x) strsplit(x, split="_")[[1]][1]),
                      type=c("raw", "corrected"))

cts_by_codon <- lapply(seq(nrow(groups)),
                       function(x) {
                         tmp_expt <- groups$expt[x]
                         tmp_type <- groups$type[x]
                         cat(paste(tmp_expt, tmp_type, "\n"))
                         # 2. aggregate counts by codon position
                         cat("... aggregating counts\n")
                         if(tmp_type=="raw") {
                           cts <- aggregate(count ~ transcript + cod_idx,
                                            data=get(paste0(tmp_expt, "_bam")),
                                            FUN=sum, na.rm=T)
                         } else {
                           cts <- aggregate(correct_200 ~ transcript + cod_idx,
                                            data=get(paste0(tmp_expt, "_bam")),
                                            FUN=sum, na.rm=T)
                           colnames(cts)[3] <- "count"
                         }
                         # 3. pull A site per codon position
                         cat("... pulling A site codon\n")
                         cts$A_site <- with(cts, get_A_codon(transcript, cod_idx,
                                                             transcript_seq=transcript_seq))
                         # 4. per transcript, compute mean density (exclude first and last 15nt)
                         cat("... calculating mean density per transcript\n")
                         buffer <- 15
                         transcript_mean <- subset(transcript_lengths, cds_length/3 > 2*buffer)
                         transcript_mean$mean_ct <- sapply(seq(nrow(transcript_mean)),
                                                           function(y) {
                                                             tmp_dat <- subset(cts,
                                                                               transcript==transcript_mean$transcript[y])
                                                             tmp_dat <- subset(tmp_dat,
                                                                               cod_idx > buffer &
                                                                                 cod_idx < (transcript_mean$cds_length[y]/3-buffer))
                                                             return(mean(tmp_dat$count))
                                                           })
                         # 5. per codon position, compute ratio of codon density to transcript mean
                         cat("... calculating normalized footprint density\n")
                         cts$norm <- with(cts, count / transcript_mean$mean_ct[match(transcript, transcript_mean$transcript)])
                         # 6. per codon: compute pause score
                         cat("... computing pause score\n")
                         pause_mean <- aggregate(norm ~ A_site, data=cts, FUN=mean)
                         pause_sd <- aggregate(norm ~ A_site, data=cts, FUN=sd)
                         # 7. compute correlation with 1/tAI
                         cat("... pulling tAI\n")
                         pause <- data.frame(codon=pause_mean$A_site,
                                             mean=pause_mean$norm,
                                             stddev=pause_sd$norm[match(pause_sd$A_site, pause_mean$A_site)],
                                             tAI=codon_properties$tAI[match(codon_properties$codon, pause_mean$A_site)],
                                             expt=tmp_expt,
                                             type=tmp_type)
                         return(pause)
                       })
cts_by_codon <- do.call(rbind, cts_by_codon)

# compute correlation -----------------------------------------------------

groups$spearman <- sapply(seq(nrow(groups)),
                          function(x) {
                            dat <- subset(cts_by_codon,
                                          expt==groups$expt[x] & type==groups$type[x])
                            dat <- subset(dat, (tAI != 0) & !is.na(tAI))
                            with(dat, cor(mean, 1/tAI, use="complete.obs", method="spearman"))
                          })
groups$pearson <- sapply(seq(nrow(groups)),
                         function(x) {
                           dat <- subset(cts_by_codon,
                                         expt==groups$expt[x] & type==groups$type[x])
                           dat <- subset(dat, (tAI != 0) & !is.na(tAI))
                           with(dat, cor(mean, 1/tAI, use="complete.obs", method="pearson"))
                         })

# plot --------------------------------------------------------------------

(ggplot(cts_by_codon,aes(x=1/tAI, y=mean, col=type)) +
  geom_point() + geom_smooth(method="lm", formula=y~x) +
  facet_grid(~expt) + theme_bw() +
  xlab("1/tAI") + ylab("mean normalized RPF density")) /
(ggplot(groups, aes(x=type, y=spearman, fill=type)) + geom_col() +
  facet_grid(~expt) + theme_bw() + theme(axis.text.x=element_text(angle=90)) +
  xlab("") + ylab("spearman correlation"))

save(cts_by_codon, groups, file="pause_score.Rda")
