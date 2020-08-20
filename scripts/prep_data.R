## requires libraries: parallel, doParallel, foreach, reshape, prodlim, Rsamtools

library(foreach)

load_lengths <- function(lengths_fname) {
  # load transcript lengths table
  ## length_fname: character; file path to transcript lengths file
  transcript_lengths <- read.table(lengths_fname, stringsAsFactors=F,
                                   col.names=c("transcript", "utr5_length", "cds_length", "utr3_length"))
  return(transcript_lengths)
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

get_codons <- function(transcript_name, cod_idx, utr5_length, transcript_seq) {
  # return codons corresponding to A-, P-, and E-site 
  # for footprint originating from transcript_name and A site codon cod_idx
  ## transcript_name: character; correspond to names(transcript_seq)
  ## cod_idx: integer; codon index for A site codon
  ## utr5_legthn: integer; length of 5' utr region in transcript sequence
  ## transcript_seq: character vector; transcript (+ 5' and 3' UTR regions) sequences
  A_start <- utr5_length + 3*(cod_idx-1) + 1
  A_end <- utr5_length + 3*cod_idx
  A_codon <- substr(transcript_seq[transcript_name], A_start, A_end)
  P_start <- utr5_length + 3*(cod_idx-2) + 1
  P_end <- utr5_length + 3*(cod_idx-1)
  P_codon <- substr(transcript_seq[transcript_name], P_start, P_end)
  E_start <- utr5_length + 3*(cod_idx-3) + 1
  E_end <- utr5_length + 3*(cod_idx-2)
  E_codon <- substr(transcript_seq[transcript_name], E_start, E_end)
  codons <- c(A_codon, P_codon, E_codon)
  names(codons) <- c("A", "P", "E")
  return(codons)
}

get_bias_seq <- function(transcript_name, cod_idx, digest_length, utr5_length,
                         transcript_seq, bias_region, bias_length=2) {
  # get bias sequence at end of footprint
  ## transcript_name: character; transcript name, corresponds with item in names(transcript_seq)
  ## cod_idx: integer; index of A site codon
  ## digest_length: integer; d5 or d3 length between A site and footprint end
  ## utr5_length: integer; length of 5' UTR region (from lengths file)
  ## transcript_seq: character vector; transcript sequences (+ 5' and 3' UTR regions)
  ## bias_region: character; f5 or f3 (corresponding to 5' or 3' bias sequence)
  ## bias_length: integer; length of bias sequence
  if(bias_region=="f5") {
    seq_start <- utr5_length + 3*(cod_idx-1)+1 - digest_length
    seq_end <- seq_start + bias_length - 1
  } else {
    if(bias_region=="f3") {
      seq_end <- utr5_length + 3*cod_idx + digest_length
      seq_start <- seq_end - bias_length + 1
    }
  }
  bias_seq <- substr(transcript_seq[transcript_name], seq_start, seq_end)
  return(bias_seq)
}

init_data <- function(transcript_fa_fname, transcript_length_fname,
                      digest5_lengths=15:18, digest3_lengths=9:11, bias_length=2,
                      num_cores=NULL, which_transcripts=NULL, ntBase=F, jointBias=F) {
  # initialize data.frame for downstream GLM
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## digest5_lengths: integer vector; legal 5' digest lengths
  ## digest3_lengths: integer vector; legal 3' digest lengths
  ## bias_length: integer; length of bias sequence
  ## num_cores: integer; number of cores to parallelize over
  ## which_transcripts: character vector; transcripts selected for regression
  ## ntBase: logical ; whether to account for addition of non-templated base to 5' end
  ## jointBias: logical ; whether to model 5' and 3' biases jointly
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
  dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons),
                                 expand.grid(d5=digest5_lengths, d3=digest3_lengths))
  dat$f5 <- foreach(a=as.character(dat$transcript), b=dat$cod_idx, c=dat$d5, 
                    d=transcript_length$utr5_length[match(dat$transcript, transcript_length$transcript)],
                    .combine='c', .export=c("get_bias_seq")) %dopar% { 
                      get_bias_seq(a, b, c, d, transcript_seq, "f5", bias_length)
                    }
  dat$f3 <- foreach(a=as.character(dat$transcript), b=dat$cod_idx, c=dat$d3,
                    d=transcript_length$utr3_length[match(dat$transcript, transcript_length$transcript)],
                    .combine='c', .export=c("get_bias_seq")) %dopar% {
                      get_bias_seq(a, b, c, d, transcript_seq, "f3", bias_length)
                    }
  bases <- c("A", "C", "T", "G")
  if(ntBase) {
    dat <- foreach(x=seq.int(nrow(dat))) %dopar% {
      first_nt <- substr(dat$f5[x], start=1, stop=1)
      add_to <- substr(dat$f5[x], start=2, stop=bias_length)
      alt_f5s <- paste0(bases, add_to)
      nt_base <- bases
      nt_base[nt_base==first_nt] <- "-"
      return(data.frame(dat[x,], alt_f5=alt_f5s, nt_base=nt_base))
    }
    dat <- do.call(rbind, dat)
  }
  if(jointBias) {
    if(ntBase) {
      dat$joint_bias <- paste0(substr(dat$alt_f5, start=1, stop=2),
                               substr(dat$f3, start=bias_length, stop=bias_length))
    } else {
      dat$joint_bias <- paste0(substr(dat$f5, start=1, stop=2),
                               substr(dat$f3, start=bias_length, stop=bias_length))
    }
    dat$joint_bias <- as.factor(dat$joint_bias)
  }
  parallel::stopCluster(cl)
  dat$d5 <- as.factor(dat$d5)
  dat$d3 <- as.factor(dat$d3)
  dat$f5 <- as.factor(dat$f5)
  dat$f3 <- as.factor(dat$f3)
  dat$count <- 0
  return(dat)
}

load_offsets <- function(offsets_fname) {
  # load A site offset rules
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## rownames: footprint length
  ## colnames: frame (0, 1, 2)
  offsets <- read.table(offsets_fname, header=T)
  offsets <- data.frame(frame=as.vector(mapply(rep, 0:2, nrow(offsets))),
                        length=rep(as.numeric(rownames(offsets)), 3),
                        offset=c(offsets$frame_0, offsets$frame_1, offsets$frame_2))
  return(offsets)
}

load_bam <- function(bam_fname, transcript_length_fname, offsets_fname) {
  # calculate proportion of footprints within each 5' and 3' digest length combination
  ## bam_fname: character; file.path to .bam alignment file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  # read in footprints
  bam_file <- Rsamtools::BamFile(bam_fname)
  bam_param <- Rsamtools::ScanBamParam(tag=c("ZW", "MD"))
  alignment <- data.frame(Rsamtools::scanBam(bam_file, param=bam_param)[[1]])
  num_footprints <- nrow(alignment)
  print(paste("Read in", num_footprints, "total footprints"))
  print(paste("... Removing", 
              sum(is.na(alignment$rname)), 
              paste0("(", round(sum(is.na(alignment$rname)) / num_footprints * 100, 1), "%)"),
              "unaligned footprints"))
  alignment <- subset(alignment, !is.na(alignment$rname))
  # assign 5' UTR lengths
  transcript_length <- load_lengths(transcript_length_fname)
  alignment$utr5_length <- transcript_length$utr5_length[match(alignment$rname, 
                                                               transcript_length$transcript)]
  # calculate frame
  alignment$frame <- (alignment$pos - alignment$utr5_length - 1) %% 3
  # calculate 5' and 3' digest lengths
  offsets <- load_offsets(offsets_fname)
  alignment$d5 <- offsets$offset[prodlim::row.match(alignment[,c("frame", "qwidth")], 
                                                    offsets[c("frame", "length")])]
  print(paste("... Removing", 
              sum(is.na(alignment$d5)), 
              paste0("(", round(sum(is.na(alignment$d5)) / num_footprints * 100, 1), "%)"), 
              "footprints outside A site offset definitions"))
  alignment <- subset(alignment, !is.na(alignment$d5))
  # calculate 3' digest lengths
  alignment$d3 <- with(alignment, qwidth - d5 - 3)
  # calculate cod_idx, remove footprints mapped outside coding region
  alignment$cod_idx <- with(alignment, (pos + d5 - utr5_length + 2) / 3)
  alignment$cds_length <- transcript_length$cds_length[match(alignment$rname, 
                                                             transcript_length$transcript)]/3
  outside_cds <- ((alignment$cod_idx < 0) | (alignment$cod_idx > alignment$cds_length))
  print(paste("... Removing", 
              sum(outside_cds),
              paste0("(", round(sum(outside_cds) / num_footprints * 100, 1), "%)"), 
              "footprints outside CDS"))
  alignment <- subset(alignment, !outside_cds)
  # return data
  # alignment <- alignment[, c("rname", "cod_idx", "d5", "d3", "seq", "tag.ZW")]
  # colnames(alignment) <- c("transcript", "cod_idx", "d5", "d3", "seq", "count")
  return(alignment)
}

count_d5_d3 <- function(bam_dat, plot_title="") {
  # count number of footprints per d5/d3 combination
  ## bam_dat: data.frame; output from load_bam()
  ## plot_title: character; plot title for output heatmap
  subset_count <- aggregate(tag.ZW ~ d5 + d3, data=bam_dat, FUN=sum)
  subset_count <- subset_count[order(subset_count$tag.ZW, decreasing=T),]
  subset_count$proportion <- sapply(seq(nrow(subset_count)),
                                    function(x) {
                                      sum(subset_count$tag.ZW[1:x])/sum(subset_count$tag.ZW)
                                    })
  subset_count_plot <- ggplot(subset_count, aes(x=d5, y=d3, fill=tag.ZW)) + geom_tile(col="black") + 
    scale_fill_gradient(low="white", high="blue", name="Count") + theme_classic() + 
    geom_text(aes(label=paste0(round(tag.ZW/sum(tag.ZW)*100, 1), "%"))) + 
    ggtitle(plot_title) + xlab("5' digestion length") + ylab("3' digestion length") +
  return(list(counts=subset_count, plot=subset_count_plot))
}

count_footprints <- function(bam_dat, regression_data) {
  # count up footprints by transcript, A site, and digest lengths
  ## bam_dat: data.frame; output from load_bam()
  ## regression_data: data.frame; output from init_data()
  # count up footprints
  alignments <- aggregate(count ~ transcript + cod_idx + d5 + d3, data=bam_dat, FUN=sum)
  # add counts to regression data.frame
  match_rows <- prodlim::row.match(alignments[, c("rname", "cod_idx", "d5", "d3")], 
                                   regression_data[, c("transcript", "cod_idx", "d5", "d3")])
  alignments <- subset(alignments, !is.na(match_rows))
  match_rows <- match_rows[!is.na(match_rows)]
  regression_data$count[match_rows] <- alignments$count
  return(regression_data)
}


# archive -----------------------------------------------------------------

plot_profile <- function(regression_data, transcript_id, model_fit=NULL) {
  # plot ribosome profile (aggregated counts per codon position) for a transcript
  ## regression_data: data.frame; output by init_data() and count_footprints()
  ## transcript_id: character; name of transcript to be plotted
  ## model_fit: glm() object; output from performing regression
  # count up footprints per codon position
  data_subset <- subset(regression_data, transcript==transcript_id)
  data_subset_cts <- aggregate(count ~ transcript + cod_idx, data=data_subset, FUN=sum)
  data_subset_cts$type <- "data"
  # plot profile
  if(is.null(model_fit)) {
    profile_plot <- ggplot2::ggplot(data_subset_cts, aes(x=cod_idx, y=count)) + geom_line() +
      theme_bw() + xlab("codon position") + ylab("footprint count") + ggtitle(transcript_id)
  } else {
    data_pred <- predict(model_fit, newdata=data_subset, type="response")
    data_pred_cts <- cbind(data_subset, data_pred)
    data_pred_cts <- aggregate(data_pred ~ transcript + cod_idx, data=data_pred_cts, FUN=sum)
    data_pred_cts$type <-
      data_subset_cts$pred <- data_pred_cts$data_pred[prodlim::row.match(data_pred_cts[,c("transcript", "cod_idx")],
                                                                         data_subset_cts[,c("transcript", "cod_idx")])]
    data_subset_cts$type <- "model prediction"
    profile_plot <- ggplot2::ggplot(data_subset_cts) +
      geom_line(aes(x=cod_idx, y=count), col="black") +
      geom_line(aes(x=cod_idx, y=pred), col="red", alpha=0.5) +
      theme_bw() + xlab("codon position") + ylab("footprint count") + ggtitle(transcript_id) + labs(colour="")
  }
  return(profile_plot)
}

count_footprints_gene <- function(transcript_name, sam_fname, transcript_length_fname, 
                                  offsets_fname, regression_data, mapping_weights=F, ntBase=F) {
  # count up footprints by A site and digest lengths
  ## transcript_name: character; transcript to count footprints for
  ## sam_fname: character; file.path to .sam alignment file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset A site assignment rules .txt file
  ## regression_data: data.frame; output from init_data()
  ## mapping_weights: logical; whether .sam alignment file has mapping weights (i.e. output from RSEM)
  ## ntBase: logical ; whether to account for addition of non-templated base to 5' end
  rpf_sam <- system(paste("grep -w", transcript_name, sam_fname), intern=T)
  rpf_sam <- rpf_sam[!grepl("^@", rpf_sam)]
  rpf_sam <- strsplit(rpf_sam, split="\t")
  if(mapping_weights) {
    rpf_sam <- data.frame(t(sapply(rpf_sam, function(x) x[c(3, 4, 10, 13, 15)])), stringsAsFactors=F)
    colnames(rpf_sam) <- c("transcript", "pos", "rpf", "mismatch_seq", "wts")
    rpf_sam$wts <- as.numeric(sub("ZW:f:", "", rpf_sam$wts))
    rpf_sam$mismatch_1 <- grepl(":0[A|T|C|G]", rpf_sam$mismatch_seq) # mismatch in first nt of footprint
    rpf_sam$mismatch_2 <- grepl(":0[A|T|C|G]0", rpf_sam$mismatch_seq) # mismatch in first 2 nt of footprint
    rpf_sam <- subset(rpf_sam, !mismatch_2) # remove footprints with 2 non-templated bases
  } else {
    rpf_sam <- data.frame(t(sapply(rpf_sam, function(x) x[c(3, 4, 10)])), stringsAsFactors=F)
    colnames(rpf_sam) <- c("transcript", "pos", "rpf")
    rpf_sam$wts <- 1
  }
  rpf_sam <- rpf_sam[rpf_sam$transcript==transcript_name,]
  rpf_sam$pos <- as.numeric(as.character(rpf_sam$pos))
  # get 5' UTR lengths
  transcript_length <- load_lengths(transcript_length_fname)
  rpf_sam$utr5_length <- transcript_length$utr5_length[match(transcript_name, transcript_length$transcript)]
  # calculate 5' digest lengths, remove non-allowable footprints
  rpf_sam$frame <- (rpf_sam$pos - rpf_sam$utr5_length - 1) %% 3
  rpf_sam$length <- nchar(rpf_sam$rpf)
  offsets <- load_offsets(offsets_fname)
  rpf_sam$d5 <- offsets$offset[prodlim::row.match(rpf_sam[,c("frame", "length")], 
                                                  offsets[c("frame", "length")])]
  rpf_sam <- subset(rpf_sam, !is.na(rpf_sam$d5))
  # calculate 3' digest lengths
  rpf_sam$d3 <- with(rpf_sam, length-d5-3)
  # calculate cod_idx, remove footprints mapped outside coding region
  rpf_sam$cod_idx <- with(rpf_sam, (pos+d5-utr5_length+2)/3)
  rpf_sam$cds_length <- transcript_length$cds_length[match(transcript_name, transcript_length$transcript)]/3
  rpf_sam <- subset(rpf_sam, rpf_sam$cod_idx > 0)
  rpf_sam <- subset(rpf_sam, rpf_sam$cod_idx < rpf_sam$cds_length)
  # count up footprints
  if(ntBase) {
    rpf_sam$f5 <- substr(rpf_sam$rpf, start=1, stop=2)
    rpf_sam$d5[rpf_sam$mismatch_1] <- rpf_sam$d5[rpf_sam$mismatch_1] - 1
    rpf_sam <- aggregate(wts ~ cod_idx + d5 + d3 + f5, data=rpf_sam, FUN=sum)
  } else {
    rpf_sam <- aggregate(wts ~ cod_idx + d5 + d3, data=rpf_sam, FUN=sum)
  }
  # add counts to regression data.frame
  regression_data <- subset(regression_data, regression_data$transcript==transcript_name)
  variables <- c("cod_idx", "d5", "d3")
  if(ntBase) {
    rpf_rows <- prodlim::row.match(rpf_sam[, c(variables, "f5")], 
                                   regression_data[, c(variables, "alt_f5")])
  } else {
    rpf_rows <- prodlim::row.match(rpf_sam[, variables], regression_data[, variables])
  }
  # take care of mismatches
  rpf_sam <- subset(rpf_sam, !is.na(rpf_rows))
  rpf_rows <- rpf_rows[!is.na(rpf_rows)]
  # return data.frame with footprint counts
  regression_data$count[rpf_rows] <- rpf_sam$wts
  return(regression_data)
}
