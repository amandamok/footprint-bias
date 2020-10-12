## requires libraries: parallel, doParallel, foreach, reshape, prodlim, Rsamtools, ggplot2

## requires get_bias_seq() from helper.R

library(foreach)
library(ggplot2)

load_bam <- function(bam_fname, transcript_length_fname, offsets_fname,
                     f5_length=3, f3_length=3, full=F, nt_base=F) {
  # calculate proportion of footprints within each 5' and 3' digest length combination
  ## bam_fname: character; file.path to .bam alignment file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## f5_length: integer; length of 5' bias region
  ## f3_length: integer; length of 3' bias region
  ## full: logical; whether to import all fields in .bam alignment file
  ##
  # 1. read in footprints
  bam_file <- Rsamtools::BamFile(bam_fname)
  if(full) {
    features <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual", "qwidth")
  } else {
    features <- c("rname", "pos", "seq", "qwidth")
  }
  bam_param <- Rsamtools::ScanBamParam(tag=c("ZW", "MD"), what=features)
  alignment <- data.frame(Rsamtools::scanBam(bam_file, param=bam_param)[[1]])
  num_footprints <- nrow(alignment)
  print(paste("Read in", num_footprints, "total footprints"))
  print(paste("... Removing",
              sum(is.na(alignment$rname)),
              paste0("(", round(sum(is.na(alignment$rname)) / num_footprints * 100, 1), "%)"),
              "unaligned footprints"))
  alignment <- subset(alignment, !is.na(alignment$rname))
  # 2. assign 5' UTR lengths
  transcript_length <- load_lengths(transcript_length_fname)
  alignment$utr5_length <- transcript_length$utr5_length[match(alignment$rname,
                                                               transcript_length$transcript)]
  # 3. calculate frame
  alignment$frame <- (alignment$pos - alignment$utr5_length - 1) %% 3
  # 4. calculate 5' and 3' digest lengths
  offsets <- load_offsets(offsets_fname)
  alignment$d5 <- offsets$offset[prodlim::row.match(alignment[,c("frame", "qwidth")],
                                                    offsets[c("frame", "length")])]
  print(paste("... Removing",
              sum(is.na(alignment$d5)),
              paste0("(", round(sum(is.na(alignment$d5)) / num_footprints * 100, 1), "%)"),
              "footprints outside A site offset definitions"))
  alignment <- subset(alignment, !is.na(alignment$d5))
  # 5. calculate 3' digest lengths
  alignment$d3 <- with(alignment, qwidth - d5 - 3)
  # 6. calculate cod_idx, remove footprints mapped outside coding region
  alignment$cod_idx <- with(alignment, (pos + d5 - utr5_length + 2) / 3)
  alignment$cds_length <- transcript_length$cds_length[match(alignment$rname,
                                                             transcript_length$transcript)]/3
  outside_cds <- ((alignment$cod_idx <= 0) | (alignment$cod_idx > alignment$cds_length))
  print(paste("... Removing",
              sum(outside_cds),
              paste0("(", round(sum(outside_cds) / num_footprints * 100, 1), "%)"),
              "footprints outside CDS"))
  alignment <- subset(alignment, !outside_cds)
  # 7. pull bias sequences
  alignment$f5 <- substr(alignment$seq, 1, f5_length)
  alignment$f3 <- mapply(substr, alignment$seq, alignment$qwidth-f3_length+1, alignment$qwidth)
  invalid_bias_seq <- (grepl("N", alignment$f5) | grepl("N", alignment$f3))
  print(paste("... Removing",
              sum(invalid_bias_seq),
              paste0("(", round(sum(invalid_bias_seq) / num_footprints * 100, 1), "%)"),
              "footprints with N in bias region"))
  alignment <- subset(alignment, !invalid_bias_seq)
  alignment$f5 <- factor(alignment$f5)
  alignment$f3 <- factor(alignment$f3)
  # 8. return nt_base
  if(nt_base) {
    alignment$nt_base <- alignment$tag.MD
    levels(alignment$nt_base) <- sapply(levels(alignment$nt_base),
                                        function(x) {
                                          ifelse(grepl("^0(A|T|C|G)", x), substr(x, 2, 2), "-")
                                        })
    print(paste("...",
                sum(alignment$nt_base != "-"),
                paste0("(", round(sum(alignment$nt_base != "-") / nrow(alignment) * 100, 1), "%)"),
                "footprints with non-templated 5' base"))
    print(paste("... ... A:",
                sum(alignment$nt_base=="A"),
                paste0("(", round(sum(alignment$nt_base=="A") / sum(alignment$nt_base!="-") * 100, 1), "%)")))
    print(paste("... ... T:",
                sum(alignment$nt_base=="A"),
                paste0("(", round(sum(alignment$nt_base=="T") / sum(alignment$nt_base!="-") * 100, 1), "%)")))
    print(paste("... ... C:",
                sum(alignment$nt_base=="A"),
                paste0("(", round(sum(alignment$nt_base=="C") / sum(alignment$nt_base!="-") * 100, 1), "%)")))
    print(paste("... ... G:",
                sum(alignment$nt_base=="A"),
                paste0("(", round(sum(alignment$nt_base=="G") / sum(alignment$nt_base!="-") * 100, 1), "%)")))
  }
  # return data
  colnames(alignment)[colnames(alignment)=="rname"] <- "transcript"
  colnames(alignment)[colnames(alignment)=="tag.ZW"] <- "count"
  subset_features <- c("transcript", "cod_idx", "d5", "d3", "f5", "f3", "count")
  if(nt_base) { subset_features <- c(subset_features, "nt_base") }
  if(!full) { alignment <- alignment[, subset_features] }
  return(alignment)
}

init_data <- function(transcript_fa_fname, transcript_length_fname,
                      digest5_lengths=15:18, digest3_lengths=9:11, d5_d3_subsets=NULL,
                      f5_length=2, f3_length=2, num_cores=NULL, which_transcripts=NULL) {
  # initialize data.frame for downstream GLM
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## digest5_lengths: integer vector; legal 5' digest lengths
  ## digest3_lengths: integer vector; legal 3' digest lengths
  ## d5_d3_subsets: data.frame; columns "d5" and "d3" of d5/d3 subsets to initiate data over
  ## bias_length: integer; length of bias sequence
  ## num_cores: integer; number of cores to parallelize over
  ## which_transcripts: character vector; transcripts selected for regression
  transcript_seq <- load_fa(transcript_fa_fname)
  transcript_length <- load_lengths(transcript_length_fname)
  if(!is.null(which_transcripts)) {
    transcript_seq <- transcript_seq[which_transcripts]
    transcript_length <- subset(transcript_length, transcript %in% which_transcripts)
  }
  transcript <- unlist(mapply(rep, x=transcript_length$transcript, times=transcript_length$cds_length/3))
  cod_idx <- unlist(lapply(transcript_length$cds_length/3, seq))
  if(is.null(num_cores)) {
    num_cores <- parallel::detectCores()-8
  }
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  codons <- foreach(a=transcript, b=cod_idx,
                    c=transcript_length$utr5_len[match(transcript, transcript_length$transcript)],
                    .combine='rbind', .export=c("get_codons")) %dopar% {
                      get_codons(a, b, c, transcript_seq)
                    }
  if(!is.null(d5_d3_subsets)) {
    dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons),
                                   d5_d3_subsets)
  } else {
    dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons),
                                   expand.grid(d5=digest5_lengths, d3=digest3_lengths))
  }
  dat$f5 <- foreach(a=as.character(dat$transcript), b=dat$cod_idx, c=dat$d5,
                    d=transcript_length$utr5_length[match(dat$transcript, transcript_length$transcript)],
                    .combine='c', .export=c("get_bias_seq")) %dopar% {
                      get_bias_seq(a, b, c, d, transcript_seq, "f5", f5_length)
                    }
  dat$f3 <- foreach(a=as.character(dat$transcript), b=dat$cod_idx, c=dat$d3,
                    d=transcript_length$utr3_length[match(dat$transcript, transcript_length$transcript)],
                    .combine='c', .export=c("get_bias_seq")) %dopar% {
                      get_bias_seq(a, b, c, d, transcript_seq, "f3", f3_length)
                    }
  bases <- c("A", "C", "T", "G")
  parallel::stopCluster(cl)
  dat$d5 <- as.factor(dat$d5)
  dat$d3 <- as.factor(dat$d3)
  dat$f5 <- as.factor(dat$f5)
  dat$f3 <- as.factor(dat$f3)
  dat$count <- 0
  return(dat)
}

count_d5_d3 <- function(bam_dat, plot_title="") {
  # count number of footprints per d5/d3 combination
  ## bam_dat: data.frame; output from load_bam()
  ## plot_title: character; plot title for output heatmap
  subset_count <- aggregate(count ~ d5 + d3, data=bam_dat, FUN=sum)
  subset_count <- subset_count[order(subset_count$count, decreasing=T),]
  subset_count$proportion <- sapply(seq(nrow(subset_count)),
                                    function(x) {
                                      sum(subset_count$count[1:x])/sum(subset_count$count)
                                    })
  subset_count_plot <- ggplot(subset_count, aes(x=d5, y=d3, fill=count)) + geom_tile(col="black") +
    scale_fill_gradient(low="white", high="blue", name="Count") + theme_classic() +
    geom_text(aes(label=paste0(round(count/sum(count)*100, 1), "%"))) +
    ggtitle(plot_title) + xlab("5' digestion length") + ylab("3' digestion length")
  return(list(counts=subset_count, plot=subset_count_plot))
}

count_footprints <- function(bam_dat, regression_data, which_column="count") {
  # count up footprints by transcript, A site, and digest lengths
  ## bam_dat: data.frame; output from load_bam()
  ## regression_data: data.frame; output from init_data()
  ## which_column: character; name of column containing counts
  # count up footprints
  bam_dat <- subset(bam_dat, transcript %in% levels(regression_data$transcript))
  bam_dat <- aggregate(formula(paste(which_column, "~ transcript + cod_idx + d5 + d3")),
                          data=bam_dat, FUN=sum, na.rm=T)
  # add counts to regression data.frame
  features <- c("transcript", "cod_idx", "d5", "d3")
  match_rows <- prodlim::row.match(regression_data[, features], bam_dat[, features])
  counts <- bam_dat[match_rows, which_column]
  counts[is.na(counts)] <- 0
  # counts <- rep(0, nrow(regression_data))
  # counts[!is.na(match_rows)] <- bam_dat[match_rows[!is.na(match_rows)], which_column]
  return(counts)
}
