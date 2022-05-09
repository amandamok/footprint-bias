rm(list=ls())

library(here)
library(choros)
library(patchwork)

project_dir <- file.path(here(), "expts", "ferguson_2021")

transcript_lengths <- read.table("~/footprint-bias/reference_data/scer.transcripts.20cds20.lengths.txt",
                                 col.names=c("transcript", "utr5", "cds", "utr3"),
                                 stringsAsFactors=F)
transcript_lengths$num_codons <- transcript_lengths$cds / 3

# load data ---------------------------------------------------------------

# load mcglincy files
mcglincy_dir <- file.path(project_dir, "mcglincy_f5_3_f3_3", "no_dedup")
load(file.path(mcglincy_dir, "mcglincy_f5_3_no_dedup_bam.Rda"))
mcglincy_bam <- mcglincy_f5_3_no_dedup_bam
load(file.path(mcglincy_dir, "mcglincy_f5_3_no_dedup_fit_200.Rda"))
mcglincy_coef_ridge <- parse_coefs(mcglincy_f5_3_no_dedup_fit_200)
load(file.path(mcglincy_dir, "mcglincy_f5_3_no_dedup_training_codon_corr.Rda"))
mcglincy_testing_codon_corr <- mcglincy_f5_3_no_dedup_training_codon_corr
load(file.path(mcglincy_dir, "mcglincy_f5_3_no_dedup_training_nt_corr.Rda"))
mcglincy_testing_nt_corr <- mcglincy_f5_3_no_dedup_training_nt_corr
load(file.path(project_dir, "mcglincy_f5_3_f3_3", "mcglincy_f5_3_raw_corr.Rda"))
rm(list=grep("f5_3", ls(), value=T))

# load ottrUMI files
ottrUMI_dir <- file.path(project_dir, "ottrUMI_Asite_f5_3_f3_3", "no_dedup")
load(file.path(ottrUMI_dir, "ottrUMI_f5_3_no_dedup_bam.Rda"))
ottrUMI_bam <- ottrUMI_f5_3_no_dedup_bam
load(file.path(ottrUMI_dir, "ottrUMI_f5_3_no_dedup_fit_200.Rda"))
ottrUMI_coef_200 <- parse_coefs(ottrUMI_f5_3_no_dedup_fit_200)
load(file.path(ottrUMI_dir, "ottrUMI_f5_3_no_dedup_training_codon_corr.Rda"))
ottrUMI_testing_codon_corr <- ottrUMI_f5_3_no_dedup_training_codon_corr
load(file.path(ottrUMI_dir, "ottrUMI_f5_3_no_dedup_training_nt_corr.Rda"))
ottrUMI_testing_nt_corr <- ottrUMI_f5_3_no_dedup_training_nt_corr
rm(list=grep("f5_3", ls(), value=T))

# functions ---------------------------------------------------------------

corr_codons_per_transcript <- function(dat, group1, group2,
                                       cor_use = "na.or.complete",
                                       which_transcripts = NULL,
                                       zero_out = F,
                                       transcripts_lengths=transcript_lengths) {
  # dat: data.frame with columns "transcript", "cod_idx", and groups
  # group1: character; column in dat
  # group2: character; column in dat
  # cor_use: character; "use" argument for cor()
  # zero_out: logical; whether to pad out unobserved codon counts with 0
  # transcripts_lengths: data.frame;
  dat <- split(dat, dat$transcript)
  if(!is.null(which_transcripts)) {
    dat <- subset(dat, names(dat) %in% which_transcripts)
  }
  sapply(dat,
         function(x) {
           tmp_dat <- x[, c("transcript", "cod_idx", group1, group2)]
           if(zero_out) {
             tmp_transcript <- as.character(unique(tmp_dat$transcript))
             num_codons <- transcripts_lengths$num_codons[transcripts_lengths$transcript==tmp_transcript]
             tmp_dat <- tmp_dat[match(seq(num_codons), tmp_dat$cod_idx),]
             tmp_dat$transcript <- tmp_transcript
             tmp_dat$cod_idx <- seq(num_codons)
             tmp_dat[is.na(tmp_dat)] <- 0
           } else {
             tmp_dat <- subset(tmp_dat,
                               rowSums(is.na(tmp_dat)) == 0)
           }
           cor(tmp_dat[[group1]], tmp_dat[[group2]],
               use=cor_use)
         })
}

log10_plus1 <- function(x) {
  # x: numeric
  log10(x + 1)
}

plot_profiles <- function(which_transcript, transcripts_lengths=transcript_lengths,
                          dat=all_codons, which_samples=NULL,
                          expt_labels = NULL) {
  # which_transcript: character; transcript id
  # transcripts_lengths: data.frame; contains columns "transcript" and "num_codons"
  # dat: data.frame; columns "transcript", "codon", and sample names
  # which_samples: character vector; which samples to plot
  # expt_labels: named character vector; new labels for experiments
  # if NULL, plot all samples
  if(is.null(which_samples)) {
    which_samples <- colnames(dat)[!(colnames(dat) %in% c("transcript", "cod_idx"))]
  }
  num_codons <- transcripts_lengths$num_codons[transcripts_lengths$transcript==which_transcript]
  dat <- subset(dat,
                transcript == which_transcript,
                select = colnames(dat) %in% c("transcript", "cod_idx", which_samples))
  dat <- dat[match(seq(num_codons), dat$cod_idx),]
  dat$transcript <- which_transcript
  dat$cod_idx <- seq(num_codons)
  dat[is.na(dat)] <- 0
  plot_data <- lapply(which_samples,
                      function(x) {
                        data.frame(dat[, c("transcript", "cod_idx")],
                                   count = dat[[x]],
                                   expt = sub("_(..|)corrected", "", x),
                                   corrected = ifelse(grepl("uncorrected", x),
                                                      "uncorrected", "corrected"),
                                   row.names = NULL)
                      })
  plot_data <- do.call(rbind, plot_data)
  if(!is.null(expt_labels)) { levels(plot_data$expt) <- expt_labels[levels(plot_data$expt)] }
  ggplot(plot_data, aes(x=cod_idx, y=count)) +
    geom_col() + theme_bw() + facet_grid(expt ~ corrected)
}

# plot regression coefficients --------------------------------------------

all_coef <- rbind(data.frame(mcglincy_coef_ridge, expt="mcglincy"),
                  data.frame(ottrUMI_coef_200, expt="ottrUMI"))

codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")

# plot A site coefficients
coef_A <- subset(all_coef, type=="A")
coef_A <- reshape(coef_A, timevar="expt",
                  idvar=c("name", "type", "coef"), direction="wide")
A_mcglincy_umi <- ggplot(coef_A, aes(x=value.mcglincy, y=value.ottrUMI)) +
  geom_point(aes(col=coef)) + theme_classic() + labs(col="") +
  geom_abline(slope=1, intercept=0) + geom_smooth(method="lm") +
  xlab("McGlincy") + ylab("OTTR UMI") + theme(legend.position="none") +
  ggtitle("A", subtitle=paste0("cor = ", signif(cor(coef_A$value.mcglincy,
                                                   coef_A$value.ottrUMI,
                                                   use="na.or.complete",
                                                   method="spearman"), 3)))

# plot 5' bias coefficients
coef_f5 <- subset(all_coef, type=="genome_f5")
coef_f5 <- reshape(coef_f5, timevar="expt",
                   idvar=c("name", "type", "coef"), direction="wide")
f5_mcglincy_umi <- ggplot(coef_f5, aes(x=value.mcglincy, y=value.ottrUMI)) +
  geom_point(aes(col=coef)) + theme_classic() + labs(col="") +
  geom_abline(slope=1, intercept=0) + geom_smooth(method="lm") +
  xlab("McGlincy") + ylab("OTTR UMI") + theme(legend.position="none") +
  ggtitle("f5", subtitle=paste0("cor = ", signif(cor(coef_f5$value.mcglincy,
                                                     coef_f5$value.ottrUMI,
                                                     use="na.or.complete",
                                                     method="spearman"), 3)))

# plot 3' bias coefficients
coef_f3 <- subset(all_coef, type=="genome_f3")
coef_f3 <- reshape(coef_f3, timevar="expt",
                   idvar=c("name", "type", "coef"), direction="wide")
f3_mcglincy_umi <- ggplot(coef_f3, aes(x=value.mcglincy, y=value.ottrUMI)) +
  geom_point(aes(col=coef)) + theme_classic() + labs(col="") +
  geom_abline(slope=1, intercept=0) + geom_smooth(method="lm") +
  xlab("McGlincy") + ylab("OTTR UMI") + theme(legend.position="none") +
  ggtitle("f3", subtitle=paste0("cor = ", signif(cor(coef_f3$value.mcglincy,
                                                   coef_f3$value.ottrUMI,
                                                   use="na.or.complete",
                                                   method="spearman"), 3)))

# plot transcript and codon abundances ------------------------------------

# load mcglincy
mcglincy_codon_cts <- aggregate(cbind(count, correct_200) ~ transcript + cod_idx,
                                data=mcglincy_bam, FUN=sum, na.rm=T, na.action=NULL)
# mcglincy_codon_cts$correct_200_na <- aggregate(correct_200 ~ transcript + cod_idx,
#                                                data=mcglincy_bam, FUN=sum,
#                                                na.action=NULL, na.rm=F)$correct_200
mcglincy_transcript_cts <- aggregate(cbind(count, correct_200) ~ transcript,
                                     data=mcglincy_bam, FUN=sum, na.rm=T, na.action=NULL)
(plot_bias(mcglincy_testing_codon_corr[,1]) +
    ggtitle("McGlincy", subtitle="uncorrected") +
    coord_cartesian(ylim=c(0, 1.05*max(unlist(mcglincy_testing_codon_corr))))) +
  (plot_bias(mcglincy_testing_codon_corr[,2]) +
     ggtitle("", subtitle="corrected") +
     coord_cartesian(ylim=c(0, 1.05*max(unlist(mcglincy_testing_codon_corr))))) +
  (plot_bias(mcglincy_testing_nt_corr[,1]) +
     ggtitle("McGlincy", subtitle="uncorrected") +
     coord_cartesian(ylim=c(0, 1.05*max(unlist(mcglincy_testing_nt_corr))))) +
  (plot_bias(mcglincy_testing_nt_corr[,2]) +
     ggtitle("", subtitle="corrected") +
     coord_cartesian(ylim=c(0, 1.05*max(unlist(mcglincy_testing_nt_corr))))) +
  (ggplot(mcglincy_codon_cts, aes(x=log10_plus1(count), y=log10_plus1(correct_200))) +
     geom_hex() + geom_abline(slope=1, intercept=0) +
     theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("log10(uncorrected)") + ylab("log10(corrected)") +
     ggtitle("", subtitle="codon counts")) +
  (ggplot(mcglincy_transcript_cts, aes(x=log10_plus1(count), y=log10_plus1(correct_200))) +
     geom_hex() + geom_abline(slope=1, intercept=0) +
     theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("log10(uncorrected)") + ylab("log10(corrected)") +
     ggtitle("", subtitle="transcript counts")) +
  plot_layout(design="
              ABE
              CDF
              ")

# load ottrUMI
ottrUMI_codon_cts <- aggregate(cbind(count, correct_200) ~ transcript + cod_idx,
                               data=ottrUMI_bam, FUN=sum, na.rm=T, na.action=NULL)
# ottrUMI_codon_cts$correct_200_na <- aggregate(correct_200 ~ transcript + cod_idx,
#                                               data=ottrUMI_bam, FUN=sum,
#                                               na.action=NULL, na.rm=F)$correct_200
ottrUMI_transcript_cts <- aggregate(cbind(count, correct_200) ~ transcript,
                                    data=ottrUMI_bam, FUN=sum, na.rm=T, na.action=NULL)
(plot_bias(ottrUMI_testing_codon_corr[,1]) +
    ggtitle("OTTR-seq w/ UMI", subtitle="uncorrected") +
    coord_cartesian(ylim=c(0, 1.05*max(ottrUMI_testing_codon_corr)))) +
  (plot_bias(ottrUMI_testing_codon_corr[,2]) +
     ggtitle("", subtitle="corrected") +
     coord_cartesian(ylim=c(0, 1.05*max(ottrUMI_testing_codon_corr)))) +
  (plot_bias(ottrUMI_testing_nt_corr[,1], type="nt") +
     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
     coord_cartesian(ylim=c(0, 1.05*max(ottrUMI_testing_nt_corr)))) +
  (plot_bias(ottrUMI_testing_nt_corr[,2], type="nt") +
     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
     coord_cartesian(ylim=c(0, 1.05*max(ottrUMI_testing_nt_corr)))) +
  (ggplot(ottrUMI_codon_cts, aes(x=log10_plus1(count), y=log10_plus1(correct_200))) +
     geom_hex() + geom_abline(slope=1, intercept=0) +
     theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("log10(uncorrected)") + ylab("log10(corrected)") +
     ggtitle("", subtitle="codon counts")) +
  (ggplot(ottrUMI_transcript_cts, aes(x=log10_plus1(count), y=log10_plus1(correct_200))) +
     geom_hex() + geom_abline(slope=1, intercept=0) +
     theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("log10(uncorrected)") + ylab("log10(corrected)") +
     ggtitle("", subtitle="transcript counts")) +
  plot_layout(design="
              ABE
              CDF
              ")

# compare codon counts across protocols -----------------------------------

match_columns <- c("transcript", "cod_idx")

# aggregate data
colnames(mcglincy_codon_cts)[colnames(mcglincy_codon_cts)=="correct_200"] <- "corrected"
colnames(ottrUMI_codon_cts)[colnames(ottrUMI_codon_cts)=="correct_200"] <- "corrected"
observed_codons <- rbind(mcglincy_codon_cts, ottrUMI_codon_cts)
observed_codons <- unique(observed_codons[, match_columns])
observed_codons <- cbind(observed_codons,
                         mcglincy_codon_cts[row.match(observed_codons[, match_columns],
                                                      mcglincy_codon_cts[, match_columns]),
                                            c("count", "corrected")],
                         ottrUMI_codon_cts[row.match(observed_codons[, match_columns],
                                                     ottrUMI_codon_cts[, match_columns]),
                                           c("count", "corrected")],
                         stringsAsFactors=F)
colnames(observed_codons) <- c(match_columns,
                               "mcglincy_uncorrected", "mcglincy_corrected",
                               "ottrUMI_uncorrected", "ottrUMI_corrected")
# compute correlations
observed_codons_corr <- list(cor.test(observed_codons$mcglincy_uncorrected,
                                      observed_codons$ottrUMI_uncorrected),
                             cor.test(observed_codons$mcglincy_corrected,
                                      observed_codons$ottrUMI_corrected))
observed_codons_corr <- data.frame(comparison = "McGlincy\nOTTR w/ UMI",
                                   correction = c("uncorrected", "corrected"),
                                   correlation = sapply(observed_codons_corr, function(x) x$estimate),
                                   lower95 = sapply(observed_codons_corr, function(x) x$conf.int[1]),
                                   upper95 = sapply(observed_codons_corr, function(x) x$conf.int[2]))
observed_codons_corr$correction <- relevel(observed_codons_corr$correction,
                                           ref="uncorrected")

# generate plot
(ggplot(observed_codons, aes(x=log10_plus1(mcglincy_uncorrected), y=log10_plus1(ottrUMI_uncorrected))) +
    geom_hex() + geom_abline(slope=1, intercept=0) +
    theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
    xlab("McGlincy: log10(uncorrected)") + ylab("OTTR w/ UMI: log10(uncorrected)") +
    ggtitle("",
            paste0("cor = ",
                   signif(subset(observed_codons_corr,
                                 comparison=="McGlincy\nOTTR w/ UMI" &
                                   correction=="uncorrected")$correlation,
                          3)))) +
  (ggplot(observed_codons, aes(x=log10_plus1(mcglincy_corrected), y=log10_plus1(ottrUMI_corrected))) +
     geom_hex() + geom_abline(slope=1, intercept=0) +
     theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("McGlincy: log10(corrected)") + ylab("OTTR w/ UMI: log10(orrected)") +
     ggtitle("",
             paste0("cor = ",
                    signif(subset(observed_codons_corr,
                                  comparison=="McGlincy\nOTTR w/ UMI" &
                                    correction=="corrected")$correlation,
                           3)))) +
  (ggplot(observed_codons_corr, aes(x=comparison, y=correlation, fill=correction,
                                    ymin=lower95, ymax=upper95)) +
     geom_col(position=position_dodge()) +
     geom_errorbar(position=position_dodge(width=0.9), width=0.5) +
     theme_classic() + xlab("") + ylab("correlation") + labs(fill="") +
     scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Set1")))

# compare transcript counts across protocols ------------------------------

# aggregate data
all_transcripts <- unique(c(as.character(mcglincy_transcript_cts$transcript),
                            as.character(ottrUMI_transcript_cts$transcript)))
colnames(mcglincy_transcript_cts)[colnames(mcglincy_transcript_cts)=="correct_200"] <- "corrected"
colnames(ottrUMI_transcript_cts)[colnames(ottrUMI_transcript_cts)=="correct_200"] <- "corrected"
all_transcripts <- cbind(all_transcripts,
                         mcglincy_transcript_cts[match(all_transcripts,
                                                       mcglincy_transcript_cts$transcript),
                                                 c("count", "corrected")],
                         ottrUMI_transcript_cts[match(all_transcripts,
                                                      ottrUMI_transcript_cts$transcript),
                                                c("count", "corrected")],
                         stringsAsFactors=F)
colnames(all_transcripts) <- c("transcript",
                               "mcglincy_uncorrected", "mcglincy_corrected",
                               "ottrUMI_uncorrected", "ottrUMI_corrected")
all_transcripts <- subset(all_transcripts,
                          rowSums(is.na(all_transcripts[, -1])) != 6)

# compute correlations
all_transcripts_corr <- list(cor.test(all_transcripts$mcglincy_uncorrected,
                                      all_transcripts$ottrUMI_uncorrected),
                             cor.test(all_transcripts$mcglincy_corrected,
                                      all_transcripts$ottrUMI_corrected))
all_transcripts_corr <- data.frame(comparison = "McGlincy\nOTTR w/ UMI",
                                   correction = c("uncorrected", "corrected"),
                                   correlation = sapply(all_transcripts_corr, function(x) x$estimate),
                                   lower95 = sapply(all_transcripts_corr, function(x) x$conf.int[1]),
                                   upper95 = sapply(all_transcripts_corr, function(x) x$conf.int[2]))
all_transcripts_corr$correction <- relevel(all_transcripts_corr$correction,
                                           ref="uncorrected")

# generate plot
(ggplot(all_transcripts, aes(x=log10_plus1(mcglincy_uncorrected), y=log10_plus1(ottrUMI_uncorrected))) +
    geom_hex() + geom_abline(slope=1, intercept=0) +
    theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
    xlab("McGlincy: log10(uncorrected)") + ylab("OTTR w/ UMI: log10(uncorrected)") +
    ggtitle("",
            paste0("cor = ",
                   signif(subset(all_transcripts_corr,
                                 comparison=="McGlincy\nOTTR w/ UMI" &
                                   correction=="uncorrected")$correlation,
                          3)))) +
  (ggplot(all_transcripts, aes(x=log10_plus1(mcglincy_corrected), y=log10_plus1(ottrUMI_corrected))) +
     geom_hex() + geom_abline(slope=1, intercept=0) +
     theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("McGlincy: log10(corrected)") + ylab("OTTR w/ UMI: log10(corrected)") +
     ggtitle("",
             paste0("cor = ",
                    signif(subset(all_transcripts_corr,
                                  comparison=="McGlincy\nOTTR w/ UMI" &
                                    correction=="corrected")$correlation,
                           3)))) +
  (ggplot(all_transcripts_corr, aes(x=comparison, y=correlation, fill=correction,
                                    ymin=lower95, ymax=upper95)) +
     geom_col(position=position_dodge()) +
     geom_errorbar(position=position_dodge(width=0.9), width=0.5) +
     theme_classic() + xlab("") + ylab("correlation") + labs(fill="") +
     scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Set1")))

# codon correlation, per transcript ---------------------------------------

# per transcript, % codons with read counts
codons_with_cts <- sapply(seq(nrow(transcript_lengths)),
                          function(x) {
                            tmp_transcript <- transcript_lengths$transcript[x]
                            tmp_num_codons <- transcript_lengths$num_codons[x]
                            tmp_codons <- subset(observed_codons,
                                                 transcript == tmp_transcript,
                                                 select=grepl("corrected", colnames(observed_codons)))
                            colSums(!is.na(tmp_codons)) / tmp_num_codons
                          })
codons_with_cts <- data.frame(transcript = transcript_lengths$transcript,
                              t(codons_with_cts))

(ggplot(codons_with_cts, aes(x=mcglincy_uncorrected, y=ottrUMI_uncorrected)) +
    geom_hex() + theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
    geom_abline(slope=1, intercept=0) + xlab("McGlincy") + ylab("OTTR w/ UMI") +
    ggtitle("", subtitle="uncorrected")) +
  (ggplot(codons_with_cts, aes(x=mcglincy_corrected, y=ottrUMI_corrected)) +
     geom_hex() + theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
     geom_abline(slope=1, intercept=0) + xlab("McGlincy") + ylab("OTTR w/ UMI") +
     ggtitle("", subtitle="corrected"))

# per transcript, average read count per codon
codons_avg_cts <- sapply(seq(nrow(transcript_lengths)),
                         function(x) {
                           tmp_transcript <- transcript_lengths$transcript[x]
                           tmp_num_codons <- transcript_lengths$num_codons[x]
                           tmp_codons <- subset(observed_codons,
                                                transcript == tmp_transcript,
                                                select=grepl("corrected", colnames(observed_codons)))
                           colSums(tmp_codons, na.rm = T) / tmp_num_codons
                         })
codons_avg_cts <- data.frame(transcript = transcript_lengths$transcript,
                             t(codons_avg_cts))

(ggplot(codons_avg_cts, aes(x=log10(mcglincy_uncorrected+1), y=log10(ottrUMI_uncorrected+1))) +
    geom_hex() + theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
    geom_abline(slope=1, intercept=0) + xlab("log10(McGlincy)") + ylab("log10(OTTR w/ UMI)") +
    ggtitle("", subtitle="uncorrected")) +
  (ggplot(codons_avg_cts, aes(x=log10(mcglincy_corrected+1), y=log10(ottrUMI_corrected+1))) +
     geom_hex() + theme_classic() + scale_fill_gradient(low="lightgrey", high="blue") +
     geom_abline(slope=1, intercept=0) + xlab("log10(McGlincy)") + ylab("log10(OTTR w/ UMI)") +
     ggtitle("", subtitle="corrected"))

codon_corr_per_transcript <- data.frame(uncorrected = corr_codons_per_transcript(observed_codons,
                                                                                 "mcglincy_uncorrected",
                                                                                 "ottrUMI_uncorrected",
                                                                                 zero_out=T),
                                        corrected = corr_codons_per_transcript(observed_codons,
                                                                               "mcglincy_corrected",
                                                                               "ottrUMI_corrected",
                                                                               zero_out=T))
codon_corr_per_transcript$avg_codon_ct <- sapply(rownames(codon_corr_per_transcript),
                                                 function(x) {
                                                   tmp_cts <- subset(all_transcripts,
                                                                     transcript == x,
                                                                     select=grepl("_corrected",
                                                                                  colnames(all_transcripts)))
                                                   tmp_num_codons <- transcript_lengths$num_codons[transcript_lengths$transcript==x]
                                                   return(mean(unlist(tmp_cts), na.rm=T) / tmp_num_codons)
                                                 })
codon_corr_per_transcript <- subset(codon_corr_per_transcript, !is.na(avg_codon_ct)) # unobserved transcripts
codon_corr_per_transcript$transcript <- rownames(codon_corr_per_transcript)
codon_corr_per_transcript$num_codons <- transcript_lengths$num_codons[match(rownames(codon_corr_per_transcript),
                                                                            transcript_lengths$transcript)]
codon_corr_per_transcript$top100 <- (nrow(codon_corr_per_transcript) -
                                       rank(codon_corr_per_transcript$avg_codon_ct) + 1) <= 100

ggplot(codon_corr_per_transcript, aes(x=uncorrected,
                                      y=corrected,
                                      size=avg_codon_ct,
                                      col=top100)) +
  geom_point(shape=16, alpha=0.5) + theme_classic() +
  geom_point(data=subset(codon_corr_per_transcript, top100),
             shape=16, col="blue", alpha=0.25) +
  scale_color_manual(values=c("grey", "blue")) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  xlab("uncorrected") + ylab("corrected") +
  labs(col="top 100 translated\ntranscript", size="average read count\nper codon") +
  ggtitle("Per-transcript correlation of codon counts")

ggplot(data.frame(comparison=codon_corr_per_transcript$comparison,
                  corr=c(codon_corr_per_transcript$uncorrected,
                         codon_corr_per_transcript$corrected),
                  corrected=rep(c("uncorrected", "corrected"),
                                each=nrow(codon_corr_per_transcript)),
                  top100=ifelse(codon_corr_per_transcript$top100,
                                "top 100", "not top 100")),
       aes(x=comparison, y=corr, fill=corrected)) + theme_classic() +
  geom_violin(draw_quantiles=0.5) + theme_classic() +
  facet_grid(top100~., scales="free_y") + ylab("per-transcript correlation")

all_profiles <- lapply(levels(observed_codons$transcript),
                       function(x) {
                         plot_profiles(x)
                       })
names(all_profiles) <- levels(observed_codons$transcript)

# look at d5/d3 subsets ---------------------------------------------------

load(file.path(mcglincy_dir, "mcglincy_subsets.Rda"))
mcglincy_d5_d3 <- d5_d3
mcglincy_d5_d3_subsets <- d5_d3_subsets

load(file.path(ottr_dir, "OTTR_RepC_subsets.Rda"))
ottr_d5_d3 <- d5_d3
ottr_d5_d3_subsets <- d5_d3_subsets

load(file.path(ottrUMI_dir, "ottrUMI_RepC_subsets.Rda"))
ottrUMI_d5_d3 <- d5_d3
ottrUMI_d5_d3_subsets <- d5_d3_subsets

(mcglincy_d5_d3$plot + ggtitle("McGlincy")) +
  (ottr_d5_d3$plot + ggtitle("OTTR-seq")) +
  (ottrUMI_d5_d3$plot + ggtitle("OTTR-seq w/ UMI"))

(ggplot(mcglincy_bam, aes(x=count, y=is.na(correct_200), fill=is.na(correct_200))) +
    geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
    theme(legend.position="none", axis.text.y=element_blank()) +
    ggtitle("McGlincy") + xlab("read count") + ylab("")) +
  (ggplot(ottr_bam, aes(x=count, y=is.na(correct_200), fill=is.na(correct_200))) +
     geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
     theme(legend.position="none", axis.text.y=element_blank()) +
     ggtitle("OTTR-seq") + xlab("read count") + ylab("")) +
  (ggplot(ottrUMI_bam, aes(x=count, y=is.na(correct_200), fill=is.na(correct_200))) +
     geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
     theme(legend.position="none", axis.text.y=element_blank()) +
     ggtitle("OTTR-seq w/ UMI") + xlab("read count") + ylab("")) +
  (ggplot(mcglincy_codon_cts, aes(x=count, y=is.na(correct_200_na), fill=is.na(correct_200_na))) +
     geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
     theme(legend.position="none", axis.text.y=element_blank()) +
     xlab("codon count") + ylab("")) +
  (ggplot(ottr_codon_cts, aes(x=count, y=is.na(correct_200_na), fill=is.na(correct_200_na))) +
     geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
     theme(legend.position="none", axis.text.y=element_blank()) +
     xlab("codon count") + ylab("")) +
  (ggplot(ottrUMI_codon_cts, aes(x=count, y=is.na(correct_200_na), fill=is.na(correct_200_na))) +
     geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
     theme(legend.position="none", axis.text.y=element_blank()) +
     xlab("codon count") + ylab("")) +
  (ggplot(mcglincy_codon_cts, aes(x=correct_200, y=is.na(correct_200_na), fill=is.na(correct_200_na))) +
     geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
     theme(legend.position="none", axis.text.y=element_blank()) +
     xlab("corrected codon count") + ylab("")) +
  (ggplot(ottr_codon_cts, aes(x=correct_200, y=is.na(correct_200_na), fill=is.na(correct_200_na))) +
     geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
     theme(legend.position="bottom", axis.text.y=element_blank()) +
     xlab("corrected codon count") + ylab("") + labs(fill="corrected count\nis NA")) +
  (ggplot(ottrUMI_codon_cts, aes(x=correct_200, y=is.na(correct_200_na), fill=is.na(correct_200_na))) +
     geom_violin(alpha=0.5) + theme_classic() + coord_cartesian(xlim=c(0, 100)) +
     theme(legend.position="none", axis.text.y=element_blank()) +
     xlab("corrected codon count") + ylab(""))

