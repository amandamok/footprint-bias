rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

sk1_scores <- read.table("~/iXnos/expts/sk1_vegetative/process/sk1_transcript_scores.txt")
sk1_scores$label <- " "
sk1_scores$label[rownames(sk1_scores)=="YPL119C"] <- "Dbp1"
sk1_scores$label[rownames(sk1_scores)=="YOR204W"] <- "Ded1"

sk1_lengths <- read.table(file.path(here(), "reference_data", 
                                    "sk1.transcripts.20cds20.lengths.txt"),
                          col.names=c("gene", "utr5", "cds", "utr3"))
sk1_scores$cds_length <- sk1_lengths$cds[match(rownames(sk1_scores), sk1_lengths$gene)]
sk1_scores$sk1_score_norm <- with(sk1_scores, sk1_score/cds_length)
sk1_scores$weinberg_score_norm <- with(sk1_scores, weinberg_score/cds_length)

# remove 2 outliers
sk1_scores <- sk1_scores[-which.max(sk1_scores$weinberg_score_norm),]
sk1_scores <- sk1_scores[-which.max(sk1_scores$weinberg_score_norm),]

scer_scores <- read.table("~/iXnos/expts/sk1_vegetative/process/scer_transcript_scores.txt")
scer_scores$label <- " "
scer_scores$label[rownames(scer_scores)=="YPL119C"] <- "Dbp1"
scer_scores$label[rownames(scer_scores)=="YOR204W"] <- "Ded1"

scer_lengths <- read.table(file.path(here(), "reference_data",
                                     "scer.transcripts.20cds20.lengths.txt"),
                           col.names=c("gene", "utr5", "cds", "utr3"))
scer_scores$cds_length <- scer_lengths$cds[match(rownames(scer_scores), scer_lengths$gene)]
scer_scores$sk1_score_norm <- with(scer_scores, sk1_score/cds_length)
scer_scores$weinberg_score_norm <- with(scer_scores, weinberg_score/cds_length)

(ggplot(sk1_scores, aes(x=weinberg_score_norm, y=sk1_score_norm, col=label)) + 
    geom_point(alpha=0.1) + geom_abline(slope=1, intercept=0, col="red") + 
    coord_fixed(ratio=1) + theme_bw() + labs(col="") + 
    geom_point(data=subset(sk1_scores, label != " ")) + 
    scale_color_manual(values=c(" "="black", "Dbp1"="blue", "Ded1"="yellow")) + 
    ggtitle("Vegetative", subtitle="SK1 transcripts") + 
    coord_cartesian(xlim=c(0.1, 0.6))) + 
  (ggplot(scer_scores, aes(x=weinberg_score_norm, y=sk1_score_norm, col=label)) + 
     geom_point(alpha=0.1) + geom_abline(slope=1, intercept=0, col="red") + 
     coord_fixed(ratio=1) + theme_bw() + labs(col="") + 
     geom_point(data=subset(sk1_scores, label != " ")) + 
     scale_color_manual(values=c(" "="black", "Dbp1"="blue", "Ded1"="yellow")) + 
     ggtitle("", subtitle="Scer transcripts") + 
     coord_cartesian(xlim=c(0.1, 0.6))) + 
  (ggplot(sk1_scores, aes(x=sk1_score_norm)) + 
     geom_histogram(binwidth=0.005, fill="lightgrey") + theme_bw() + 
     annotate(geom="point", y=0, col="blue", 
              x=sk1_scores$sk1_score_norm[sk1_scores$label=="Dbp1"]) + 
     annotate(geom="point", y=0, col="yellow", 
              x=sk1_scores$sk1_score_norm[sk1_scores$label=="Ded1"]) + 
     coord_cartesian(xlim=c(0.1, 0.6))) + 
  (ggplot(scer_scores, aes(x=sk1_score_norm)) + 
     geom_histogram(binwidth=0.005, fill="lightgrey") + theme_bw() + 
     annotate(geom="point", y=0, col="blue", 
              x=scer_scores$sk1_score_norm[scer_scores$label=="Dbp1"]) + 
     annotate(geom="point", y=0, col="yellow", 
              x=scer_scores$sk1_score_norm[scer_scores$label=="Ded1"]) + 
     coord_cartesian(xlim=c(0.1, 0.6))) +
  (ggplot(sk1_scores, aes(x=weinberg_score_norm)) + 
     geom_histogram(binwidth=0.005, fill="lightgrey") + theme_bw() + 
     annotate(geom="point", y=0, col="blue", 
              x=sk1_scores$weinberg_score_norm[sk1_scores$label=="Dbp1"]) + 
     annotate(geom="point", y=0, col="yellow", 
              x=sk1_scores$weinberg_score_norm[sk1_scores$label=="Ded1"]) + 
     coord_cartesian(xlim=c(0.1, 0.6))) + 
  (ggplot(scer_scores, aes(x=weinberg_score_norm)) + 
     geom_histogram(binwidth=0.005, fill="lightgrey") + theme_bw() + 
     annotate(geom="point", y=0, col="blue", 
              x=scer_scores$weinberg_score_norm[scer_scores$label=="Dbp1"]) + 
     annotate(geom="point", y=0, col="yellow", 
              x=scer_scores$weinberg_score_norm[scer_scores$label=="Ded1"]) + 
     coord_cartesian(xlim=c(0.1, 0.6))) + 
  plot_layout(ncol=2, byrow=T)

