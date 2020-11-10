rm(list=ls())

library(ggplot2)
library(patchwork)
library(reshape2)
library(prodlim)

source(file.path(here(), "scripts", "helper.R"))

# identify signal peptides in s. cerevisiae genome (SignalP) --------------

dna_fname <- file.path(here(), "reference_data", "scer.transcripts.20cds20.fa")
aa_fname <- file.path(here(), "reference_data", "scer.transcripts.aa.fa")
if(!file.exists(aa_fname)) { convert_to_aa(dna_fname, aa_fname,
                                           utr5_length=20, utr3_length=20) }

signalp_prefix <- file.path(here(), "reference_data", "scer.transcripts.aa")
if(!file.exists(paste0(signalp_prefix, "_summary.signalp5"))) {
  signalp_path <- "~/.local/bin/signalp"
  system(paste(signalp_path, "-fasta", aa_fname,
               "-org euk -gff3 -format short -prefix", signalp_prefix,
               "-tmp", file.path(here(), "reference_data")))
}

signal_peptides <- read.table(paste0(signalp_prefix, ".gff3"), stringsAsFactors=F,
                              col.names=c("seqid", "source", "type", "start",
                                          "end", "score", "strand",
                                          "phase", "attributes"))

# compare profiles before/after correction of signal peptides -------------

expts <- c("green", "lareau", "weinberg")

for(expt in expts) {
  load(file.path(here(), "expts", paste0(expt, "_20cds20_f5_2_f3_3"), paste0(expt, "_bam.Rda")))
}

transcript_lengths_fname <- file.path(here(), "reference_data", "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_lengths_fname)

sp_profiles <- lapply(seq(nrow(signal_peptides)),
                      function(x) {
                        tmp_transcript <- signal_peptides$seqid[x]
                        tmp_sp_length <- signal_peptides$end[x]
                        tmp_length <- transcript_lengths$cds_length[match(tmp_transcript,
                                                                          transcript_lengths$transcript)]/3
                        tmp_profile <- lapply(expts,
                                              function(expt) {
                                                first_codon <- 5
                                                # extract footprints assigned to transcript
                                                tmp_subset <- subset(get(paste0(expt, "_bam")),
                                                                     transcript==tmp_transcript)
                                                if(nrow(tmp_subset)== 0) {
                                                  output <- data.frame(transcript = tmp_transcript,
                                                                       cod_idx = rep(first_codon:(tmp_sp_length+5), times=4),
                                                                       expt = expt,
                                                                       count = 0,
                                                                       type = rep(c("raw", "raw (normalized)",
                                                                                    "corrected", "corrected (normalized)"),
                                                                                  times=length(first_codon:(tmp_sp_length+5))))
                                                } else {
                                                  # aggregate raw footprint counts
                                                  if(nrow(tmp_subset==1)) {
                                                    tmp_raw <- tmp_subset
                                                    tmp_correct <- tmp_subset
                                                  } else {
                                                    tmp_raw <- aggregate(count ~ transcript + cod_idx,
                                                                         data=tmp_subset, FUN=sum)
                                                    tmp_correct <- aggregate(correct_150 ~ transcript + cod_idx,
                                                                             data=tmp_subset, FUN=sum)
                                                  }
                                                  mean_raw <- sum(tmp_raw$count, na.rm=T) / tmp_length
                                                  mean_correct <- sum(tmp_correct$correct_150, na.rm=T) / tmp_length
                                                  indices_raw <- match(first_codon:(tmp_sp_length+5), tmp_raw$cod_idx)
                                                  indices_correct <- match(first_codon:(tmp_sp_length+5), tmp_correct$cod_idx)
                                                  # report raw, raw (normalized), corrected, corrected (normalized) counts
                                                  tmp <- data.frame(transcript = tmp_transcript,
                                                                    cod_idx = first_codon:(tmp_sp_length+5) - tmp_sp_length,
                                                                    expt = expt)
                                                  output <- rbind(within(tmp,
                                                                         {
                                                                           count <- tmp_raw$count[indices_raw]
                                                                           type <- "raw"
                                                                         }),
                                                                  within(tmp,
                                                                         {
                                                                           count <- tmp_raw$count[indices_raw] / mean_raw
                                                                           type <- "raw (normalized)"
                                                                         }),
                                                                  within(tmp,
                                                                         {
                                                                           count <- tmp_correct$correct_150[indices_correct]
                                                                           type <- "corrected"
                                                                         }),
                                                                  within(tmp,
                                                                         {
                                                                           count <- tmp_correct$correct_150[indices_correct] / mean_correct
                                                                           type <- "corrected (normalized)"
                                                                         }))
                                                }
                                                return(output)
                                              })
                        tmp_profile <- do.call(rbind, tmp_profile)
                        tmp_profile$count[is.na(tmp_profile$count)] <- 0
                        return(tmp_profile)
                      })
sp_profiles <- do.call(rbind, sp_profiles)
sp_profiles$type <- factor(sp_profiles$type,
                           levels=c("raw", "raw (normalized)", "corrected", "corrected (normalized)"))
sp_profiles$correction <- factor(ifelse(grepl("raw", as.character(sp_profiles$type)),
                                        "raw", "corrected"),
                                 levels=c("raw", "corrected"))
sp_profiles$normalized <- factor(ifelse(grepl("normalized", as.character(sp_profiles$type)),
                                        "+", "-"),
                                 levels=c("-", "+"))
levels(sp_profiles$expt) <- c("Green", "Lareau", "Weinberg")

sp_counts_raw <- aggregate(count ~ transcript + expt, data=subset(sp_profiles, type=="raw"), FUN=sum)
sp_counts_raw <- dcast(sp_counts_raw, transcript ~ expt)

transcripts_to_plot <- as.character(sp_counts_raw$transcript[sapply(seq(nrow(sp_counts_raw)),
                                                                function(x) { all(sp_counts_raw[x,][-1]>75) })])
sp_profiles_plots <- lapply(transcripts_to_plot,
                            function(x) {
                              tmp_subset <- subset(sp_profiles,
                                                   transcript == x & normalized == "-")
                              ggplot(tmp_subset, aes(x=cod_idx, y=count, col=correction)) +
                                geom_line() + facet_grid(expt ~ ., scales="free_y") + theme_classic() +
                                geom_hline(yintercept=0) + geom_vline(xintercept=0) +
                                ggtitle(x) +
                                xlab("distance to signal peptide cleavage site")
                            })
wrap_plots(sp_profiles_plots, ncol=5)

sp_counts_corrected <- aggregate(count ~ transcript + expt,
                                 data=subset(sp_profiles, type=="corrected"), FUN=sum)
comparison_data <- melt(sp_counts_raw, value.name="raw")
comparison_data$corrected <- sp_counts_corrected$count[row.match(sp_counts_corrected[, c("transcript", "expt")],
                                                                 comparison_data[, c("transcript", "variable")])]
ggplot(comparison_data, aes(x=raw, y=corrected)) + geom_point(alpha=0.5) +
  facet_wrap(~variable) + geom_abline(slope=1, intercept=0) + theme_classic() +
  xlab("raw counts") + ylab("corrected counts") +
  ggtitle("Signal peptide", subtitle=paste("n =", nrow(comparison_data)))
