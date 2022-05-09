rm(list=ls())

library(here)
library(choros)

expt_dir <- file.path(here(), "expts", "ferguson_2021")

transcript_lengths_fname <- file.path(here(), "reference_data",
                                      "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_lengths_fname)

# functions ---------------------------------------------------------------

log10_plus1 <- function(x) {
  # x: numeric
  log10(x + 1)
}

# mcglincy: aggregate counts ----------------------------------------------

mcglincy_bam_fname <- file.path(expt_dir, "mcglincy_f5_3_f3_3",
                                "test_gc", "mcglincy_f5_3_gc_bam.Rda")
load(mcglincy_bam_fname)

mcglincy_codon_cts <- aggregate(cbind(count, correct_400_gc) ~ transcript + cod_idx,
                                mcglincy_f5_3_bam, sum, na.action=NULL, na.rm=T)
colnames(mcglincy_codon_cts) <- c("transcript", "cod_idx", "count", "corrected")
mcglincy_transcript_cts <- aggregate(cbind(count, corrected) ~ transcript,
                                     mcglincy_codon_cts, sum, na.rm=T)

# ottrUMI: aggregate counts -----------------------------------------------

ottrUMI_bam_fname <- file.path(expt_dir, "ottrUMI_Asite_f5_3_f3_3",
                               "test_gc", "ottrUMI_f5_3_gc_bam.Rda")
load(ottrUMI_bam_fname)

ottrUMI_codon_cts <- aggregate(cbind(count, correct_200_gc) ~ transcript + cod_idx,
                               ottrUMI_f5_3_bam, sum, na.action=NULL, na.rm=T)
colnames(ottrUMI_codon_cts) <- c("transcript", "cod_idx", "count", "corrected")
ottrUMI_transcript_cts <- aggregate(cbind(count, corrected) ~ transcript,
                                    ottrUMI_codon_cts, sum, na.rm=T)

# compute correlations: by codon position ---------------------------------

all_codons <- data.frame(transcript=unlist(mapply(rep,
                                                  transcript_lengths$transcript,
                                                  transcript_lengths$cds_length/3)),
                         cod_idx=unlist(mapply(seq, transcript_lengths$cds_length/3)),
                         row.names=NULL)
all_codons$mcglincy_raw <- mcglincy_codon_cts$count[match_rows(all_codons, mcglincy_codon_cts, c("transcript", "cod_idx"))]
all_codons$mcglincy_corrected <- mcglincy_codon_cts$corrected[match_rows(all_codons, mcglincy_codon_cts, c("transcript", "cod_idx"))]
all_codons$ottrUMI_raw <- ottrUMI_codon_cts$count[match_rows(all_codons, ottrUMI_codon_cts, c("transcript", "cod_idx"))]
all_codons$ottrUMI_corrected <- ottrUMI_codon_cts$corrected[match_rows(all_codons, ottrUMI_codon_cts, c("transcript", "cod_idx"))]
all_codons[is.na(all_codons)] <- 0

all_codons_cor <- lapply(c("raw", "corrected"),
                         function(x) {
                           mcglincy <- all_codons[[paste0("mcglincy_", x)]]
                           ottrUMI <- all_codons[[paste0("ottrUMI_", x)]]
                           mcglincy[is.na(mcglincy)] <- 0
                           ottrUMI[is.na(ottrUMI)] <- 0
                           tmp_test <- cor.test(mcglincy, ottrUMI)
                           return(data.frame(estimate.cor=tmp_test$estimate,
                                             conf.int1=tmp_test$conf.int[1],
                                             conf.int2=tmp_test$conf.int[2],
                                             type=x))
                         })
all_codons_corr <- data.frame(do.call(rbind, all_codons_cor))

all_codons_plot <- (ggplot(all_codons,
                           aes(x=log10_plus1(mcglincy_raw), y=log10_plus1(ottrUMI_raw))) +
                      geom_hex() + geom_abline(slope=1, intercept=0) + theme_classic() +
                      scale_fill_gradient(low="lightgrey", high="blue") +
                      xlab("McGlincy") + ylab("OTTR w/ UMI") + ggtitle("", subtitle="raw counts")) +
  (ggplot(all_codons,
          aes(x=log10_plus1(mcglincy_corrected), y=log10_plus1(ottrUMI_corrected))) +
     geom_hex() + geom_abline(slope=1, intercept=0) + theme_classic() +
     scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("McGlincy") + ylab("OTTR w/ UMI") + ggtitle("", subtitle="corrected counts")) +
  (ggplot(all_codons_corr,
          aes(x=factor(type, levels=c("raw", "corrected")),
              y=estimate.cor, ymin=conf.int1, ymax=conf.int2, fill=type)) +
     geom_col() + geom_errorbar(width=0.5) + theme_classic() +
     xlab("") + ylab(expr(rho)) + guides(fill="none") +
     scale_fill_manual(values=c("raw"="red", "corrected"="blue")) +
     coord_cartesian(ylim=c(0.5, 0.7)))

# compute correlations: by transcript -------------------------------------

all_transcripts <- data.frame(transcript=transcript_lengths$transcript,
                              row.names=NULL)
all_transcripts$mcglincy_raw <- mcglincy_transcript_cts$count[match(all_transcripts$transcript, mcglincy_transcript_cts$transcript)]
all_transcripts$mcglincy_corrected <- mcglincy_transcript_cts$corrected[match(all_transcripts$transcript, mcglincy_transcript_cts$transcript)]
all_transcripts$ottrUMI_raw <- ottrUMI_transcript_cts$count[match(all_transcripts$transcript, ottrUMI_transcript_cts$transcript)]
all_transcripts$ottrUMI_corrected <- ottrUMI_transcript_cts$corrected[match(all_transcripts$transcript, ottrUMI_transcript_cts$transcript)]
all_transcripts[is.na(all_transcripts)] <- 0

all_transcripts_cor <- lapply(c("raw", "corrected"),
                              function(x) {
                                mcglincy <- all_transcripts[[paste0("mcglincy_", x)]]
                                ottrUMI <- all_transcripts[[paste0("ottrUMI_", x)]]
                                tmp_test <- cor.test(mcglincy, ottrUMI)
                                return(data.frame(estimate.cor=tmp_test$estimate,
                                                  conf.int1=tmp_test$conf.int[1],
                                                  conf.int2=tmp_test$conf.int[2],
                                                  type=x))
                              })
all_transcripts_corr <- data.frame(do.call(rbind, all_transcripts_cor))

all_transcripts_plot <- (ggplot(all_transcripts,
                                aes(x=log10_plus1(mcglincy_raw), y=log10_plus1(ottrUMI_raw))) +
                           geom_hex() + geom_abline(slope=1, intercept=0) + theme_classic() +
                           scale_fill_gradient(low="lightgrey", high="blue") +
                           xlab("McGlincy") + ylab("OTTR w/ UMI") + ggtitle("", subtitle="raw counts")) +
  (ggplot(all_transcripts,
          aes(x=log10_plus1(mcglincy_corrected), y=log10_plus1(ottrUMI_corrected))) +
     geom_hex() + geom_abline(slope=1, intercept=0) + theme_classic() +
     scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("McGlincy") + ylab("OTTR w/ UMI") + ggtitle("", subtitle="corrected counts")) +
  (ggplot(all_transcripts_corr,
          aes(x=factor(type, levels=c("raw", "corrected")),
              y=estimate.cor, ymin=conf.int1, ymax=conf.int2, fill=type)) +
     geom_col() + geom_errorbar(width=0.5) + theme_classic() +
     xlab("") + ylab(expr(rho)) + guides(fill="none") +
     scale_fill_manual(values=c("raw"="red", "corrected"="blue")) +
     coord_cartesian(ylim=c(0.8, 1)))

# compute correlations by transcript --------------------------------------

cor_by_transcript <- lapply(transcript_lengths$transcript,
                            function(x) {
                              tmp_codons <- subset(all_codons, transcript==x)
                              raw_cor <- cor.test(tmp_codons$mcglincy_raw,
                                                  tmp_codons$ottrUMI_raw)
                              corrected_cor <- cor.test(tmp_codons$mcglincy_corrected,
                                                        tmp_codons$ottrUMI_corrected)
                              return(data.frame(raw_estimate=raw_cor$estimate,
                                                raw_conf.int1=raw_cor$conf.int[1],
                                                raw_conf.int2=raw_cor$conf.int[2],
                                                corrected_estimate=corrected_cor$estimate,
                                                corrected_conf.int1=corrected_cor$conf.int[1],
                                                corrected_conf.int2=corrected_cor$conf.int[2]))
                            })
cor_by_transcript <- do.call(rbind, cor_by_transcript)
cor_by_transcript$avg_abundance <- rowMeans(within(all_transcripts,
                                                   {
                                                     mcglincy_raw <- mcglincy_raw/mean(mcglincy_raw)
                                                     mcglincy_corrected <- mcglincy_corrected/mean(mcglincy_corrected)
                                                     ottrUMI_raw <- ottrUMI_raw/mean(ottrUMI_raw)
                                                     ottrUMI_corrected <- ottrUMI_corrected/mean(ottrUMI_corrected)
                                                   })[, -1])[match(transcript_lengths$transcript,
                                                                   all_transcripts$transcript)]
cor_by_transcript$top_200 <- ifelse(cor_by_transcript$avg_abundance >=
                                      sort(cor_by_transcript$avg_abundance, decreasing=T)[200],
                                    "top 200 genes", "other")

ggplot(cor_by_transcript,
       aes(x=raw_estimate, y=corrected_estimate,
           size=avg_abundance, col=top_200)) +
  geom_point(alpha=0.25) + geom_abline(slope=1, intercept=0) + theme_classic() +
  geom_point(data=subset(cor_by_transcript, top_200=="top 200 genes"), alpha=0.5) +
  labs(size="average\nnormalized\nabundance", col="") +
  scale_color_manual(values=c("top 200 genes"="red", "other"="grey")) +
  xlab("cor(raw)") + ylab("cor(corrected)") + ggtitle("correlation by transcript")

