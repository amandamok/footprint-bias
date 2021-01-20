rm(list=ls())

library(here)
library(Rsamtools)
library(ggplot2)
library(patchwork)

data_dir <- file.path(here(), "expts", "wu_2019", "raw_data")
bam_param <- Rsamtools::ScanBamParam(tag=c("ZW", "MD"), what=c("rname", "pos", "qwidth"))

# expts <- system(paste0("cat ", data_dir, "/*txt"), intern=T)

plot_dir <- file.path(here(), "expts", "wu_2019", "raw_data", "plots")

sra_table <- read.csv(file.path(data_dir, "GSE115162_SraRunTable.csv"), stringsAsFactors=F)
expts <- sra_table$Run[sra_table$Organism=="Saccharomyces cerevisiae" &
                         sra_table$molecule_subtype=="Ribosome Protected mRNA (15-34 nt)"]

lengths_fname <- file.path(here(), "reference_data", "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- read.table(lengths_fname, stringsAsFactors=F,
                                 col.names=c("transcript", "utr5", "cds", "utr3"))

for(x in expts) {
  print(x)

  pdf_fname <- file.path(plot_dir, paste0(x, ".pdf"))
  dat_fname <- file.path(plot_dir, paste0(x, ".Rda"))

  if(!file.exists(pdf_fname)) {
    if(!file.exists(dat_fname)) {
      # load bam alignment
      bam_fname <- file.path(data_dir, x, paste0(x, ".transcript.bam"))
      bam_file <- Rsamtools::BamFile(bam_fname)
      alignment <- data.frame(Rsamtools::scanBam(bam_file, param=bam_param)[[1]])

      # remove unaligned reads
      alignment <- subset(alignment, !is.na(rname))
      num_footprints <- nrow(alignment)

      # aggregate counts by read length, start position, and transcript
      size_by_pos <- rbind(data.frame(aggregate(tag.ZW ~ pos + qwidth + rname, data=alignment, FUN=length),
                                      type="read count"),
                           data.frame(aggregate(tag.ZW ~ pos + qwidth + rname, data=alignment, FUN=sum),
                                      type="RSEM weight"))
      colnames(size_by_pos) <- c("position", "length", "transcript", "count", "type")

      # restrict footprint sizes (15-34 nt)
      size_by_pos <- subset(size_by_pos, length >=15 & length<=34)

      # compute distance to start and stop codons
      size_by_pos$d5 <- size_by_pos$pos - (20 + 3 + 1)
      size_by_pos$d3 <- (size_by_pos$pos + size_by_pos$length - 1) -
        (transcript_lengths$cds[match(as.character(size_by_pos$transcript), transcript_lengths$transcript)] + 20 + 1)

      # subset to start codon
      start_codon <- subset(size_by_pos, d5 < 0)
      start_codon <- aggregate(count ~ d5 + length + type, data=start_codon, FUN=sum)

      # subset to stop codon
      stop_codon <- subset(size_by_pos, d3 > 0)
      stop_codon <- aggregate(count ~ d3 + length + type, data=stop_codon, FUN=sum)

      # calculate size distribution
      fp_sizes <- rbind(data.frame(aggregate(tag.ZW ~ qwidth, data=alignment, FUN=length),
                                   type="read count"),
                        data.frame(aggregate(tag.ZW ~ qwidth, data=alignment, FUN=sum),
                                   type="RSEM weight"))
      colnames(fp_sizes) <- c("length", "count", "type")
      fp_sizes$prop <- sapply(seq(nrow(fp_sizes)),
                              function(x) {
                                ifelse(fp_sizes$type[x]=="read count",
                                       fp_sizes$count[x] / sum(fp_sizes$count[fp_sizes$type=="read count"]),
                                       fp_sizes$count[x] / sum(fp_sizes$count[fp_sizes$type=="RSEM weight"]))
                              })

      save(size_by_pos, start_codon, stop_codon, fp_sizes,
           file=dat_fname)
    } else {
      load(dat_fname)
    }

    num_footprints <- sum(size_by_pos$count[size_by_pos$type=="read count"])

    # plot
    start_codon_plot <- ggplot(start_codon, aes(x=d5, y=length, fill=count)) + geom_tile(linetype=1, col=1) +
      theme_classic() + scale_fill_gradient(low="white", high="darkblue") + facet_grid(type~.) +
      scale_x_continuous(breaks=seq(min(start_codon$d5), max(start_codon$d5))) +
      scale_y_continuous(breaks=15:34) + xlab("5' digest length") +
      ggtitle(x, subtitle=paste(sra_table$X.GEO.label.[match(x, sra_table$Run)],
                                 "(", num_footprints, "alignments )")) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))
    stop_codon_plot <- ggplot(stop_codon, aes(x=d3, y=length, fill=count)) + geom_tile(linetype=1, col=1) +
      theme_classic() + scale_fill_gradient(low="white", high="darkblue") + facet_grid(type ~.) +
      scale_x_continuous(breaks=seq(min(stop_codon$d3), max(stop_codon$d3))) +
      scale_y_continuous(breaks=15:34) + xlab("3' digest length") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))
    sizes_plot <- ggplot(fp_sizes, aes(x=length, y=prop)) + geom_col() + facet_grid(type~.) + theme_classic() +
      scale_x_continuous(breaks=seq(min(fp_sizes$length), max(fp_sizes$length))) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5))

    ggsave(plot=wrap_plots(start_codon_plot, stop_codon_plot, sizes_plot, nrow=1),
           filename=pdf_fname, device="pdf", width=10, height=5, units="in")
  }
}

q(save="no")