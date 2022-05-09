rm(list=ls())

library(here)
library(choros)
library(Rsamtools)
library(prodlim)

bam_features <- c("rname", "pos", "seq", "qwidth")
bam_param <- Rsamtools::ScanBamParam(tag=c("ZW", "MD"), what=bam_features)

ref_dir <- file.path(here(), "reference_data")
transcript_lengths_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_lengths_fname)

# load mcglincy data
mcglincy_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "mcglincy",
                                "mcglincy_trimmed_deduplicated_footprints.transcript.bam")
mcglincy_bam_file <- data.frame(scanBam(BamFile(mcglincy_bam_fname), param=bam_param)[[1]])
load(file.path(here(), "expts", "ferguson_2021", "mcglincy", "mcglincy_bam.Rda"))

# load ottrUMI data
ottrUMI_bam_fname <- file.path(here(), "expts", "ferguson_2021", "raw_data", "WT_UMI_OTTR_repB",
                               "ottrUMI_trimmed_deduplicated_footprints.transcript.bam")
ottrUMI_bam_file <- data.frame(scanBam(BamFile(ottrUMI_bam_fname), param=bam_param)[[1]])
load(file.path(here(), "expts", "ferguson_2021", "ottrUMI_Asite", "ottrUMI_bam.Rda"))

# count up RPFs by codon --------------------------------------------------

mcglincy_codon_cts <- aggregate(cbind(count, correct_200) ~ transcript + cod_idx,
                                mcglincy_bam, FUN=sum, na.rm=T)
ottrUMI_codon_cts <- aggregate(cbind(count, correct_200) ~ transcript + cod_idx,
                               ottrUMI_bam, FUN=sum, na.rm=T)

mcglincy_transcript_cts <- aggregate(cbind(count, correct_200) ~ transcript,
                                     mcglincy_codon_cts, FUN=sum, na.rm=T)
mcglincy_transcript_cts[, c("count", "correct_200")] <- mcglincy_transcript_cts[, c("count", "correct_200")] / sum(mcglincy_bam$count)
ottrUMI_transcript_cts <- aggregate(cbind(count, correct_200) ~ transcript,
                                    ottrUMI_codon_cts, FUN=sum, na.rm=T)
ottrUMI_transcript_cts[, c("count", "correct_200")] <- ottrUMI_transcript_cts[, c("count", "correct_200")] / sum(ottrUMI_bam$count)
transcript_cts <- rowMeans(cbind(mcglincy_transcript_cts$count[match(transcript_lengths$transcript, mcglincy_transcript_cts$transcript)],
                                 mcglincy_transcript_cts$correct_200[match(transcript_lengths$transcript, mcglincy_transcript_cts$transcript)],
                                 ottrUMI_transcript_cts$count[match(transcript_lengths$transcript, ottrUMI_transcript_cts$transcript)],
                                 ottrUMI_transcript_cts$correct_200[match(transcript_lengths$transcript, ottrUMI_transcript_cts$transcript)]),
                           na.rm=T)
names(transcript_cts) <- transcript_lengths$transcript
transcript_cts <- sort(transcript_cts, decreasing=T)

codon_cts <- lapply(which(transcript_lengths$transcript %in% names(transcript_cts)[1:50]),
                    function(x) {
                      tmp_transcript <- transcript_lengths$transcript[x]
                      tmp_cds_length <- transcript_lengths$cds_length[x] / 3
                      tmp_codon_cts <- data.frame(transcript=tmp_transcript, cod_idx = seq(tmp_cds_length))
                      tmp_codon_cts$mcglincy_raw <- mcglincy_codon_cts$count[row.match(data.frame(tmp_transcript, tmp_codon_cts$cod_idx),
                                                                                       mcglincy_codon_cts[, c("transcript", "cod_idx")])]
                      tmp_codon_cts$mcglincy_corrected <- mcglincy_codon_cts$correct_200[row.match(data.frame(tmp_transcript, tmp_codon_cts$cod_idx),
                                                                                                   mcglincy_codon_cts[, c("transcript", "cod_idx")])]
                      tmp_codon_cts$ottrUMI_raw <- ottrUMI_codon_cts$count[row.match(data.frame(tmp_transcript, tmp_codon_cts$cod_idx),
                                                                                     ottrUMI_codon_cts[, c("transcript", "cod_idx")])]
                      tmp_codon_cts$ottrUMI_corrected <- ottrUMI_codon_cts$correct_200[row.match(data.frame(tmp_transcript, tmp_codon_cts$cod_idx),
                                                                                                 ottrUMI_codon_cts[, c("transcript", "cod_idx")])]
                      tmp_codon_cts[is.na(tmp_codon_cts)] <- 0
                      return(tmp_codon_cts)
                    })

transcript_profiles <- lapply(codon_cts,
                              function(x) {
                                plot_data <- melt(x, measure.vars=3:6)
                                plot_data$variable <- factor(plot_data$variable,
                                                             levels=c("mcglincy_raw", "mcglincy_corrected",
                                                                      "ottrUMI_raw", "ottrUMI_corrected"))
                                levels(plot_data$variable) <- c("McGlincy: raw", "McGlincy: corrected",
                                                                "OTTR: raw", "OTTR: corrected")
                                ggplot(plot_data, aes(x=cod_idx, y=value, fill=variable)) +
                                  geom_col() + facet_grid(variable ~.) + theme_bw() +
                                  theme(legend.position="none") + xlab("codon") + ylab("RPF count") +
                                  ggtitle(plot_data$transcript[1])
                              })
