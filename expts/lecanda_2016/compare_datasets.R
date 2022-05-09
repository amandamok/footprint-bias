rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

paper <- "lecanda_2016"
paper_dir <- file.path(here(), "expts", paper)

expts <- list.dirs(paper_dir, full.names=F, recursive=F)
expts <- expts[expts != "raw_data"]

scer_lengths_fname <- file.path(here(), "reference_data", "scer.transcripts.20cds20.lengths.txt")
scer_lengths <- load_lengths(scer_lengths_fname)

# load data ---------------------------------------------------------------

for(expt in expts) {
  load(file.path(paper_dir, expt, "correct_gc", paste0(expt, "_gc_coef_250.Rda")))
  load(file.path(paper_dir, expt, "correct_gc", paste0(expt, "_bam.Rda")))
}

# compare regression coefficients -----------------------------------------

all_coefs <- lapply(expts,
                    function(expt) {
                      get(paste0(expt, "_gc_coef_250"))[, c("name", "group", "term")]
                    })
all_coefs <- unique(do.call(rbind, all_coefs))

all_coefs$fixed_standard_est <- fixedLinker_standardPrimer_gc_coef_250$estimate[match_rows(all_coefs,
                                                                                           fixedLinker_standardPrimer_gc_coef_250,
                                                                                           c("name", "group", "term"))]
all_coefs$fixed_standard_se <- fixedLinker_standardPrimer_gc_coef_250$std_error[match_rows(all_coefs,
                                                                                           fixedLinker_standardPrimer_gc_coef_250,
                                                                                           c("name", "group", "term"))]
all_coefs$random_standard_est <- randomLinker_standardPrimer_gc_coef_250$estimate[match_rows(all_coefs,
                                                                                             randomLinker_standardPrimer_gc_coef_250,
                                                                                             c("name", "group", "term"))]
all_coefs$random_standard_se <- randomLinker_standardPrimer_gc_coef_250$std_error[match_rows(all_coefs,
                                                                                             randomLinker_standardPrimer_gc_coef_250,
                                                                                             c("name", "group", "term"))]
all_coefs$random_random_est <- randomLinker_randomPrimer_gc_coef_250$estimate[match_rows(all_coefs,
                                                                                         randomLinker_randomPrimer_gc_coef_250,
                                                                                         c("name", "group", "term"))]
all_coefs$random_random_se <- randomLinker_randomPrimer_gc_coef_250$std_error[match_rows(all_coefs,
                                                                                         randomLinker_randomPrimer_gc_coef_250,
                                                                                         c("name", "group", "term"))]
all_coefs$group[is.na(all_coefs$group)] <- "intercept"
all_coefs$group <- factor(all_coefs$group,
                          levels=unique(all_coefs$group))

(ggplot(subset(all_coefs, !(group %in% c("transcript", "gc", "intercept"))),
        aes(x=fixed_standard_est, y=random_standard_est,
            xmin=fixed_standard_est + qnorm(0.025)*fixed_standard_se,
            xmax=fixed_standard_est + qnorm(0.975)*fixed_standard_se,
            ymin=random_standard_est + qnorm(0.025)*random_standard_se,
            ymax=random_standard_est + qnorm(0.975)*random_standard_se,
            col=group)) + theme_bw() +
    geom_abline(slope=1, intercept=0) + geom_point() +
    geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5) +
    xlab("fixed linker / standard primer") + ylab("random linker / standard primer") +
    facet_grid(~group) + guides(col="none")) +
  (ggplot(subset(all_coefs, !(group %in% c("transcript", "gc", "intercept"))),
          aes(x=fixed_standard_est, y=random_random_est,
              xmin=fixed_standard_est + qnorm(0.025)*fixed_standard_se,
              xmax=fixed_standard_est + qnorm(0.975)*fixed_standard_se,
              ymin=random_random_est + qnorm(0.025)*random_random_se,
              ymax=random_random_est + qnorm(0.975)*random_random_se,
              col=group)) + theme_bw() +
     geom_abline(slope=1, intercept=0) + geom_point() +
     geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5) +
     xlab("fixed linker / standard primer") + ylab("random linker / random primer") +
     facet_grid(~group) + guides(col="none")) +
  (ggplot(subset(all_coefs, !(group %in% c("transcript", "gc", "intercept"))),
          aes(x=random_standard_est, y=random_random_est,
              xmin=random_standard_est + qnorm(0.025)*random_standard_se,
              xmax=random_standard_est + qnorm(0.975)*random_standard_se,
              ymin=random_random_est + qnorm(0.025)*random_random_se,
              ymax=random_random_est + qnorm(0.975)*random_random_se,
              col=group)) + theme_bw() +
     geom_abline(slope=1, intercept=0) + geom_point() +
     geom_errorbar(alpha=0.5) + geom_errorbarh(alpha=0.5) +
     xlab("random linker / standard primer") + ylab("random linker / random primer") +
     facet_grid(~group) + guides(col="none")) +
  plot_layout(nrow=3)

# aggregate footprint counts ----------------------------------------------

scer_lengths$num_codons <- with(scer_lengths, cds_length/3)

cts_by_codon <- data.frame(transcript=unlist(mapply(rep, scer_lengths$transcript,
                                                    scer_lengths$num_codons,
                                                    simplify=F)),
                           cod_idx=unlist(lapply(scer_lengths$num_codons, seq.int)),
                           row.names=NULL)

expts_codon_cts <- lapply(expts,
                          function(expt) {
                            tmp <- aggregate(cbind(count, correct_250, correct_250_gc) ~ transcript + cod_idx,
                                             data=get(paste0(expt, "_bam")),
                                             FUN=sum, na.rm=T, na.action=na.pass)
                            colnames(tmp) <- c("transcript", "cod_idx",
                                               paste0(expt, "_raw"),
                                               paste0(expt, "_corrected"),
                                               paste0(expt, "_corrected_gc"))
                            return(tmp)
                          })
cts_by_codon <- cbind(cts_by_codon,
                      expts_codon_cts[[1]][match_rows(cts_by_codon,
                                                      expts_codon_cts[[1]],
                                                      c("transcript", "cod_idx")),
                                           c(3:5)],
                      expts_codon_cts[[2]][match_rows(cts_by_codon,
                                                      expts_codon_cts[[2]],
                                                      c("transcript", "cod_idx")),
                                           c(3:5)],
                      expts_codon_cts[[3]][match_rows(cts_by_codon,
                                                      expts_codon_cts[[3]],
                                                      c("transcript", "cod_idx")),
                                           c(3:5)])
cts_by_codon[is.na(cts_by_codon)] <- 0

rm(list=grep("_bam", ls(), value=T))

# correlate counts (across entire dataset) --------------------------------

plot_codon_cor <- function(dat, x, y, xlabel=x, ylabel=y) {
  # dat: data.frame; contains columns x and y
  # x: character; column in dat to plot as x axis
  # y: character; column in dat to plot as y axis
  # xlabel: character; label for x axis
  # ylabel: character; label for y axis
  ggplot(dat, aes_string(x=x, y=y)) +
    geom_hex() + geom_abline(slope=1, intercept=0) + theme_classic() +
    scale_fill_gradient(low="lightgrey", high="blue") +
    scale_x_log10() + scale_y_log10() + xlab(xlabel) + ylab(ylabel)
}

to_evaluate <- c("raw", "corrected", "corrected_gc")

corr_by_dataset <- sapply(to_evaluate,
                          function(x) {
                            cor_1 <- with(cts_by_codon,
                                          cor.test(cts_by_codon[, paste0(expts[1], "_", x)],
                                                   cts_by_codon[, paste0(expts[2], "_", x)]))
                            cor_2 <- with(cts_by_codon,
                                          cor.test(cts_by_codon[, paste0(expts[1], "_", x)],
                                                   cts_by_codon[, paste0(expts[3], "_", x)]))
                            cor_3 <- with(cts_by_codon,
                                          cor.test(cts_by_codon[, paste0(expts[2], "_", x)],
                                                   cts_by_codon[, paste0(expts[3], "_", x)]))
                            return(list(cor_1=cor_1, cor_2=cor_2, cor_3=cor_3))
                          }, simplify=F, USE.NAMES=T)

corr_by_dataset_plots <- sapply(to_evaluate,
                                function(x) {
                                  plot_1 <- plot_codon_cor(cts_by_codon[, -c(1,2)]+1,
                                                           paste(expts[1], x, sep="_"),
                                                           paste(expts[2], x, sep="_"),
                                                           expts[1], expts[2]) +
                                    ggtitle(x, subtitle=paste("cor =",
                                                              signif(corr_by_dataset[[x]]$cor_1$estimate,
                                                                     3)))
                                  plot_2 <- plot_codon_cor(cts_by_codon[, -c(1,2)]+1,
                                                           paste(expts[1], x, sep="_"),
                                                           paste(expts[3], x, sep="_"),
                                                           expts[1], expts[3]) +
                                    ggtitle("", subtitle=paste("cor =",
                                                               signif(corr_by_dataset[[x]]$cor_2$estimate,
                                                                      3)))
                                  plot_3 <- plot_codon_cor(cts_by_codon[, -c(1,2)]+1,
                                                           paste(expts[2], x, sep="_"),
                                                           paste(expts[3], x, sep="_"),
                                                           expts[2], expts[3]) +
                                    ggtitle("", subtitle=paste("cor =",
                                                               signif(corr_by_dataset[[x]]$cor_3$estimate,
                                                                      3)))
                                  return(plot_1 + plot_2 + plot_3)
                                }, simplify=F, USE.NAMES=T)


corr_by_dataset_barplot_dat <- unlist(corr_by_dataset, recursive=F)
corr_by_dataset_barplot_dat <- data.frame(estimate=sapply(corr_by_dataset_barplot_dat,
                                                          function(x) x$estimate),
                                          ci_lower=sapply(corr_by_dataset_barplot_dat,
                                                          function(x) x$conf.int[1]),
                                          ci_upper=sapply(corr_by_dataset_barplot_dat,
                                                          function(x) x$conf.int[2]),
                                          correction=factor(rep(to_evaluate, each=3),
                                                            levels=to_evaluate),
                                          comparison=rep(c(paste(expts[1], expts[2], sep="\n"),
                                                           paste(expts[1], expts[3], sep="\n"),
                                                           paste(expts[2], expts[3], sep="\n")),
                                                         times=length(to_evaluate)))

corr_by_dataset_plots$raw /
  corr_by_dataset_plots$corrected /
  corr_by_dataset_plots$corrected_gc/
  (ggplot(corr_by_dataset_barplot_dat,
          aes(x=correction, y=estimate, ymin=ci_lower, ymax=ci_upper, fill=correction)) +
     geom_col() + geom_errorbar(width=0.5) + theme_classic() + facet_grid(~comparison) +
     xlab("") + ylab("correlation") + guides(fill="none") +
     coord_cartesian(ylim=c(0.9, 1))) +
  plot_layout(heights=c(1, 1, 1, 1.5))

# correlate counts (by transcript) ----------------------------------------

cor_by_transcript <- lapply(split(cts_by_codon, cts_by_codon$transcript),
                            function(x) {
                              cor_by_eval <- sapply(to_evaluate,
                                                    function(y) {
                                                      cor_1 <- cor(x[, paste0(expts[1], "_", y)],
                                                                   x[, paste0(expts[2], "_", y)],
                                                                   method="spearman")
                                                      cor_2 <- cor(x[, paste0(expts[1], "_", y)],
                                                                   x[, paste0(expts[3], "_", y)],
                                                                   method="spearman")
                                                      cor_3 <- cor(x[, paste0(expts[2], "_", y)],
                                                                   x[, paste0(expts[3], "_", y)],
                                                                   method="spearman")
                                                      tmp <- list(cor_1, cor_2, cor_3)
                                                      names(tmp) <- paste(y, c("1v2", "1v3", "2v3"), sep="_")
                                                      return(tmp)
                                                    }, simplify=F, USE.NAMES=F)
                              return(unlist(cor_by_eval))
                            })
cor_by_transcript <- data.frame(do.call(rbind, cor_by_transcript))

cts_by_transcript <- aggregate(cbind(fixedLinker_standardPrimer_raw,
                                     randomLinker_standardPrimer_raw,
                                     randomLinker_randomPrimer_raw) ~ transcript,
                               data=cts_by_codon, FUN=sum)
cts_by_transcript <- cts_by_transcript[match(rownames(cor_by_transcript),
                                             cts_by_transcript$transcript),]
cor_by_transcript$mean_ct <- rowMeans(cts_by_transcript[, -1])
cor_by_transcript$top_250 <- cor_by_transcript$mean_ct >= sort(cor_by_transcript$mean_ct,
                                                               decreasing=T)[250]
cor_by_transcript$top_100 <- cor_by_transcript$mean_ct >= sort(cor_by_transcript$mean_ct,
                                                               decreasing=T)[100]

plot_cor_by_transcript <- function(dat, x, y, xlabel=x, ylabel=y, color=NULL) {
  # dat: data.frame; contains columns x, y, col
  # x: character; column in dat to plot as x axis
  # y: character; column in dat to plot as y axis
  # xlabel: character; label for x axis
  # ylabel: character; label for y axis
  # color: character; column in dat containing logicals
  if(is.null(color)) {
    ggplot(dat, aes_string(x=x, y=y)) + theme_classic() +
      geom_point(alpha=0.5, shape=16) + geom_abline(slope=1, intercept=0) +
      xlab(xlabel) + ylab(ylabel)
  } else {
    ggplot(dat, aes_string(x=x, y=y, col=color)) + theme_classic() +
      geom_point(alpha=0.25, shape=16) + geom_abline(slope=1, intercept=0) +
      geom_point(data=subset(dat, dat[[color]]), alpha=0.5, shape=16) +
      scale_color_manual(values=c("TRUE"="blue", "FALSE"="grey")) +
      xlab(xlabel) + ylab(ylabel)
  }
}

(plot_cor_by_transcript(cor_by_transcript, "raw_1v2", "corrected_1v2",
                       "raw", "corrected", "top_100") +
    ggtitle(paste0(expts[1], "\n", expts[2]))) +
  (plot_cor_by_transcript(cor_by_transcript, "raw_1v2", "corrected_gc_1v2",
                          "raw", "corrected_gc", "top_100") +
     ggtitle(paste0(expts[1], "\n", expts[2]))) +
  (plot_cor_by_transcript(cor_by_transcript, "raw_1v3", "corrected_1v3",
                          "raw", "corrected", "top_100") +
     ggtitle(paste0(expts[1], "\n", expts[3]))) +
  (plot_cor_by_transcript(cor_by_transcript, "raw_1v3", "corrected_gc_1v3",
                          "raw", "corrected_gc", "top_100")) +
  (plot_cor_by_transcript(cor_by_transcript, "raw_2v3", "corrected_2v3",
                          "raw", "corrected", "top_100") +
     ggtitle(paste0(expts[2], "\n", expts[3]))) +
  (plot_cor_by_transcript(cor_by_transcript, "raw_2v3", "corrected_gc_2v3",
                          "raw", "corrected_gc", "top_100")) +
  plot_layout(ncol=2)

(plot_cor_by_transcript(subset(cor_by_transcript, top_250),
                        "raw_1v2", "corrected_1v2", "raw", "corrected", "top_100") +
    ggtitle(paste0(expts[1], "\n", expts[2]))) +
  (plot_cor_by_transcript(subset(cor_by_transcript, top_250),
                          "raw_1v2", "corrected_gc_1v2", "raw", "corrected_gc", "top_100")) +
  (plot_cor_by_transcript(subset(cor_by_transcript, top_250),
                          "raw_1v3", "corrected_1v3", "raw", "corrected", "top_100") +
     ggtitle(paste0(expts[1], "\n", expts[3]))) +
  (plot_cor_by_transcript(subset(cor_by_transcript, top_250),
                          "raw_1v3", "corrected_gc_1v3", "raw", "corrected_gc", "top_100")) +
  (plot_cor_by_transcript(subset(cor_by_transcript, top_250),
                          "raw_2v3", "corrected_2v3", "raw", "corrected", "top_100") +
     ggtitle(paste0(expts[2], "\n", expts[3]))) +
  (plot_cor_by_transcript(subset(cor_by_transcript, top_250),
                          "raw_2v3", "corrected_gc_2v3", "raw", "corrected_gc", "top_100")) +
  plot_layout(ncol=2)
