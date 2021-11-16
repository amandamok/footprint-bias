library(here)
library(Rsamtools)
library(choros)
library(patchwork)

data_dir <- (file.path(here(), "expts", "ferguson_2021", "raw_data"))

ref_dir <- file.path(here(), "reference_data")
transcript_lengths_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_lengths_fname)

bam_features <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual", "qwidth")
bam_param <- ScanBamParam(tag=c("ZW", "MD"), what=bam_features)

mcglincy_bam_fname <- file.path(data_dir, "mcglincy",
                                "mcglincy_trimmed_deduplicated_footprints.transcript.bam")
ottr_bam_fname <- file.path(data_dir, "WT_OTTR_repB",
                            "ottr_trimmed_footprints.transcript.bam")
ottrUMI_bam_fname <- file.path(data_dir, "WT_UMI_OTTR_repB",
                               "ottrUMI_trimmed_deduplicated_footprints.transcript.bam")

# mcglincy ----------------------------------------------------------------

mcglincy_bam <- data.frame(scanBam(BamFile(mcglincy_bam_fname),
                                   param=bam_param))
mcglincy_bam <- subset(mcglincy_bam, flag == 0 & qwidth <= 40 & qwidth >= 20)
# 8 420 247 unique reads
# 65 168 856 alignments
# 7 127 550 RSEM counts
mcglincy_bam$cds_length <- transcript_lengths$cds_length[match(mcglincy_bam$rname,
                                                               transcript_lengths$transcript)]
mcglincy_bam$utr5_length <- transcript_lengths$utr5_length[match(mcglincy_bam$rname,
                                                                 transcript_lengths$transcript)]
mcglincy_bam$start_distance <- with(mcglincy_bam, pos-(utr5_length+1+3))
mcglincy_bam$stop_distance <- with(mcglincy_bam, (pos+qwidth-1)-(utr5_length+cds_length))
# start codon ribogrid
mcglincy_start_ribogrid <- subset(mcglincy_bam, start_distance < 10)
mcglincy_start_ribogrid <- aggregate(tag.ZW ~ start_distance + qwidth,
                                     mcglincy_start_ribogrid, sum)
mcglincy_start_ribogrid_plot <- ggplot(mcglincy_start_ribogrid,
                                       aes(x=start_distance, y=qwidth, fill=tag.ZW)) +
  geom_tile(color="black") + theme_classic() + labs(fill="weighted\ncount") +
  scale_fill_gradient(low="white", high="blue") + coord_fixed(ratio=1) +
  xlab("read 5' end to start codon (P site)") + ylab("RPF length") +
  scale_x_continuous(breaks=c(-4:2)*5) + scale_y_continuous(breaks=c(2:8)*5)
# stop codon ribogrid
mcglincy_stop_ribogrid <- subset(mcglincy_bam, stop_distance > -10)
mcglincy_stop_ribogrid <- aggregate(tag.ZW ~ stop_distance + qwidth,
                                    mcglincy_stop_ribogrid, sum)
mcglincy_stop_ribogrid_plot <- ggplot(mcglincy_stop_ribogrid,
                                      aes(x=stop_distance, y=qwidth, fill=tag.ZW)) +
  geom_tile(color="black") + theme_classic() + labs(fill="weighted\ncount") +
  scale_fill_gradient(low="white", high="blue") + coord_fixed(ratio=1) +
  xlab("read 3' end to stop codon") + ylab("RPF length") +
  scale_x_continuous(breaks=c(-2:4)*5) + scale_y_continuous(breaks=c(2:8)*5)
# length/frame plots
mcglincy_bam$frame <- with(mcglincy_bam, (pos - utr5_length - 1) %% 3)
mcglincy_length_frame_cts <- aggregate(tag.ZW ~ qwidth + frame, mcglincy_bam, sum)
mcglincy_length_plot <- ggplot(mcglincy_length_frame_cts,
                               aes(x=qwidth, y=tag.ZW/sum(tag.ZW), fill=as.factor(frame))) +
  geom_col() + theme_classic() + labs(fill="frame") + scale_x_continuous(breaks=20:40) +
  xlab("RPF length") + ylab("% of total RFP count (RSEM weighted)") +
  theme(legend.position="bottom", axis.text.x=element_text(angle=90, vjust=0.5))
mcglincy_length_frame_plot <- ggplot(mcglincy_length_frame_cts,
                                     aes(x=factor(frame,levels=c(2, 1, 0)), y=qwidth, fill=tag.ZW)) +
  geom_tile() + scale_fill_gradient(low="white", high="blue") + theme_classic() +
  coord_fixed(ratio=1) + xlab("frame") + ylab("RPF length") +
  scale_y_continuous(breaks=20:40) + labs(fill="weighted\ncount")
mcglincy_plot <- (mcglincy_length_plot + ggtitle("McGlincy rep B (deduplicated)")) +
  (mcglincy_start_ribogrid_plot / mcglincy_stop_ribogrid_plot) +
  mcglincy_length_frame_plot +
  plot_layout(widths=c(5, 5, 1))

# ottr --------------------------------------------------------------------

ottr_bam <- data.frame(scanBam(BamFile(ottr_bam_fname),
                               param=bam_param))
ottr_bam <- subset(ottr_bam, flag == 0 & qwidth <= 40)
# 10 450 547 unique reads
# 17 904 733 alignments
# 10 447 518 RSEM counts
ottr_bam$cds_length <- transcript_lengths$cds_length[match(ottr_bam$rname,
                                                           transcript_lengths$transcript)]
ottr_bam$utr5_length <- transcript_lengths$utr5_length[match(ottr_bam$rname,
                                                             transcript_lengths$transcript)]
ottr_bam$start_distance <- with(ottr_bam, pos-(utr5_length+1+3))
ottr_bam$stop_distance <- with(ottr_bam, (pos+qwidth-1)-(utr5_length+cds_length))
# start codon ribogrid
ottr_start_ribogrid <- subset(ottr_bam, start_distance < 10)
ottr_start_ribogrid <- aggregate(tag.ZW ~ start_distance + qwidth,
                                 ottr_start_ribogrid, sum)
ottr_start_ribogrid_plot <- ggplot(ottr_start_ribogrid,
                                   aes(x=start_distance, y=qwidth, fill=tag.ZW)) +
  geom_tile(color="black") + theme_classic() + labs(fill="weighted\ncount") +
  scale_fill_gradient(low="white", high="blue") + coord_fixed(ratio=1) +
  xlab("read 5' end to start codon (P site)") + ylab("RPF length") +
  scale_x_continuous(breaks=c(-4:2)*5) + scale_y_continuous(breaks=c(2:8)*5)
# stop codon ribogrid
ottr_stop_ribogrid <- subset(ottr_bam, stop_distance > -10)
ottr_stop_ribogrid <- aggregate(tag.ZW ~ stop_distance + qwidth,
                                ottr_stop_ribogrid, sum)
ottr_stop_ribogrid_plot <- ggplot(ottr_stop_ribogrid,
                                  aes(x=stop_distance, y=qwidth, fill=tag.ZW)) +
  geom_tile(color="black") + theme_classic() + labs(fill="weighted\ncount") +
  scale_fill_gradient(low="white", high="blue") + coord_fixed(ratio=1) +
  xlab("read 3' end to stop codon") + ylab("RPF length") +
  scale_x_continuous(breaks=c(-2:4)*5) + scale_y_continuous(breaks=c(2:8)*5)
# length/frame plots
ottr_bam$frame <- with(ottr_bam, (pos - utr5_length - 1) %% 3)
ottr_length_frame_cts <- aggregate(tag.ZW ~ qwidth + frame, ottr_bam, sum)
ottr_length_plot <- ggplot(ottr_length_frame_cts,
                           aes(x=qwidth, y=tag.ZW/sum(tag.ZW), fill=as.factor(frame))) +
  geom_col() + theme_classic() + labs(fill="frame") + scale_x_continuous(breaks=20:40) +
  xlab("RPF length") + ylab("% of total RFP count (RSEM weighted)") +
  theme(legend.position="bottom", axis.text.x=element_text(angle=90, vjust=0.5))
ottr_length_frame_plot <- ggplot(ottr_length_frame_cts,
                                 aes(x=factor(frame,levels=c(2, 1, 0)), y=qwidth, fill=tag.ZW)) +
  geom_tile() + scale_fill_gradient(low="white", high="blue") + theme_classic() +
  coord_fixed(ratio=1) + xlab("frame") + ylab("RPF length") +
  scale_y_continuous(breaks=20:40) + labs(fill="weighted\ncount")
ottr_plot <- (ottr_length_plot + ggtitle("OTTR rep B (deduplicated)")) +
  (ottr_start_ribogrid_plot / ottr_stop_ribogrid_plot) +
  ottr_length_frame_plot +
  plot_layout(widths=c(5, 5, 1))

# ottrUMI -----------------------------------------------------------------

ottrUMI_bam <- data.frame(scanBam(BamFile(ottrUMI_bam_fname),
                                  param=bam_param))
ottrUMI_bam <- subset(ottrUMI_bam, qwidth <= 40)
# 7 573 933 unique reads
# 12 641 132 alignments
# 7 573 748 RSEM counts
ottrUMI_bam$cds_length <- transcript_lengths$cds_length[match(ottrUMI_bam$rname,
                                                              transcript_lengths$transcript)]
ottrUMI_bam$utr5_length <- transcript_lengths$utr5_length[match(ottrUMI_bam$rname,
                                                                transcript_lengths$transcript)]
ottrUMI_bam$start_distance <- with(ottrUMI_bam, pos-(utr5_length+1+3))
ottrUMI_bam$stop_distance <- with(ottrUMI_bam, (pos+qwidth-1)-(utr5_length+cds_length))
# start codon ribogrid
ottrUMI_start_ribogrid <- subset(ottrUMI_bam, start_distance < 10)
ottrUMI_start_ribogrid <- aggregate(tag.ZW ~ start_distance + qwidth,
                                    ottrUMI_start_ribogrid, sum)
ottrUMI_start_ribogrid_plot <- ggplot(ottrUMI_start_ribogrid,
                                      aes(x=start_distance, y=qwidth, fill=tag.ZW)) +
  geom_tile(color="black") + theme_classic() + labs(fill="weighted\ncount") +
  scale_fill_gradient(low="white", high="blue") + coord_fixed(ratio=1) +
  xlab("read 5' end to start codon (P site)") + ylab("RPF length") +
  scale_x_continuous(breaks=c(-4:2)*5) + scale_y_continuous(breaks=c(2:8)*5)
# stop codon ribogrid
ottrUMI_stop_ribogrid <- subset(ottrUMI_bam, stop_distance > -10)
ottrUMI_stop_ribogrid <- aggregate(tag.ZW ~ stop_distance + qwidth,
                                   ottrUMI_stop_ribogrid, sum)
ottrUMI_stop_ribogrid_plot <- ggplot(ottrUMI_stop_ribogrid,
                                     aes(x=stop_distance, y=qwidth, fill=tag.ZW)) +
  geom_tile(color="black") + theme_classic() + labs(fill="weighted\ncount") +
  scale_fill_gradient(low="white", high="blue") + coord_fixed(ratio=1) +
  xlab("read 3' end to stop codon") + ylab("RPF length") +
  scale_x_continuous(breaks=c(-2:4)*5) + scale_y_continuous(breaks=c(2:8)*5)
# length/frame plots
ottrUMI_bam$frame <- with(ottrUMI_bam, (pos - utr5_length - 1) %% 3)
ottrUMI_length_frame_cts <- aggregate(tag.ZW ~ qwidth + frame, ottrUMI_bam, sum)
ottrUMI_length_plot <- ggplot(ottrUMI_length_frame_cts,
                              aes(x=qwidth, y=tag.ZW/sum(tag.ZW), fill=as.factor(frame))) +
  geom_col() + theme_classic() + labs(fill="frame") + scale_x_continuous(breaks=20:40) +
  xlab("RPF length") + ylab("% of total RFP count (RSEM weighted)") +
  theme(legend.position="bottom", axis.text.x=element_text(angle=90, vjust=0.5))
ottrUMI_length_frame_plot <- ggplot(ottrUMI_length_frame_cts,
                                    aes(x=factor(frame,levels=c(2, 1, 0)), y=qwidth, fill=tag.ZW)) +
  geom_tile() + scale_fill_gradient(low="white", high="blue") + theme_classic() +
  coord_fixed(ratio=1) + xlab("frame") + ylab("RPF length") +
  scale_y_continuous(breaks=20:40) + labs(fill="weighted\ncount")
ottrUMI_plot <- (ottrUMI_length_plot + ggtitle("OTTR w/ UMI rep B (deduplicated)")) +
  (ottrUMI_start_ribogrid_plot / ottrUMI_stop_ribogrid_plot) +
  ottrUMI_length_frame_plot +
  plot_layout(widths=c(5, 5, 1))

# save --------------------------------------------------------------------

save(mcglincy_bam, mcglincy_length_frame_cts,
     mcglincy_start_ribogrid, mcglincy_stop_ribogrid, mcglincy_plot,
     file=file.path(data_dir, "mcglincy_asite.Rda"))
save(ottr_bam, ottr_length_frame_cts,
     ottr_start_ribogrid, ottr_stop_ribogrid, ottr_plot,
     file=file.path(data_dir, "ottr_asite.Rda"))
save(ottrUMI_bam, ottrUMI_length_frame_cts,
     ottrUMI_start_ribogrid, ottrUMI_stop_ribogrid, ottrUMI_plot,
     file=file.path(data_dir, "ottrUMI_asite.Rda"))

q(save="no")
