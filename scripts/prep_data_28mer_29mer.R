## requires libraries: parallel, doParallel, foreach, reshape, prodlim

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

init_data <- function(transcript_fa_fname, transcript_length_fname, bias_length=2,
                      num_cores=NULL, which_transcripts=NULL) {
  # initialize data.frame for downstream GLM ; assume d5=15 and d3=10 for all footprints
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## bias_length: integer; length of bias sequence
  ## num_cores: integer; number of cores to parallelize over
  ## which_transcripts: character vector; transcripts selected for regression
  transcript_seq <- load_fa(transcript_fa_fname)
  transcript_length <- load_lengths(transcript_length_fname)
  if(!is.null(which_transcripts)) {
    transcript_seq <- transcript_seq[which_transcripts]
    transcript_length <- subset(transcript_length, transcript %in% which_transcripts)
  }
  # establish transcripts
  transcript <- unlist(mapply(rep, x=transcript_length$transcript, times=transcript_length$cds_length/3))
  # establish codon indices
  cod_idx <- unlist(lapply(transcript_length$cds_length/3, seq))
  # look up codons for A, P, and E sites
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
  # look up 5' and 3' bias sequences
  dat <- data.frame(transcript, cod_idx, codons, d5=15, d3=10)
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
  # add rows for non-templated base additions
  bases <- c("A", "C", "G", "T")
  dat <- do.call(rbind, 
                 foreach(x=seq.int(nrow(dat))) %dopar% {
                   suffix <- substr(dat$f5[x], start=1, stop=(bias_length-1))
                   tmp_dat <- dat[rep(x, each=5),]
                   tmp_dat$f5[2:5] <- paste0(bases, suffix)
                   tmp_dat$nt_base <- c("-", bases)
                   return(tmp_dat) })
  dat$f5 <- as.factor(dat$f5)
  dat$f3 <- as.factor(dat$f3)
  dat$nt_base <- as.factor(dat$nt_base)
  parallel::stopCluster(cl)
  # zero out counts
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

count_footprints_gene <- function(transcript_name, sam_fname, transcript_length_fname, 
                                  offsets_fname, regression_data, mapping_weights=F, prop_assign=F) {
  # count up footprints by A site and digest lengths
  ## transcript_name: character; transcript to count footprints for
  ## sam_fname: character; file.path to .sam alignment file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset A site assignment rules .txt file
  ## regression_data: data.frame; output from init_data()
  ## mapping_weights: logical; whether .sam alignment file has mapping weights (i.e. output from RSEM)
  ## prop_assign: either FALSE or data.frame [ output from compute_29mer_weight() ]
  rpf_sam <- system(paste("grep -w", transcript_name, sam_fname), intern=T)
  rpf_sam <- rpf_sam[!grepl("^@", rpf_sam)]
  rpf_sam <- strsplit(rpf_sam, split="\t")
  if(mapping_weights) {
    rpf_sam <- data.frame(t(sapply(rpf_sam, function(x) x[c(3, 4, 10, 13, 15)])), stringsAsFactors=F)
    colnames(rpf_sam) <- c("transcript", "pos", "rpf", "mismatch_seq", "wts")
    rpf_sam$wts <- as.numeric(sub("ZW:f:", "", rpf_sam$wts))
    rpf_sam$mismatch_1 <- grepl(":0[A|T|C|G]", rpf_sam$mismatch_seq) # mismatch in first nt of footprint
    # rpf_sam$mismatch_2 <- grepl(":0[A|T|C|G]0", rpf_sam$mismatch_seq) # mismatch in first 2 nt of footprint
    # rpf_sam <- subset(rpf_sam, !mismatch_2) # remove footprints with 2 non-templated bases
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
  rpf_sam$d5 <- offsets$offset[prodlim::row.match(rpf_sam[,c("frame", "length")], offsets[c("frame", "length")])]
  rpf_sam <- subset(rpf_sam, !is.na(rpf_sam$d5))
  # calculate 3' digest lengths
  rpf_sam$d3 <- with(rpf_sam, length-d5-3)
  # calculate cod_idx, remove footprints mapped outside coding region
  rpf_sam$cod_idx <- with(rpf_sam, (pos+d5-utr5_length+2)/3)
  rpf_sam$cds_length <- transcript_length$cds_length[match(transcript_name, transcript_length$transcript)]/3
  rpf_sam <- subset(rpf_sam, rpf_sam$cod_idx > 0)
  rpf_sam <- subset(rpf_sam, rpf_sam$cod_idx < rpf_sam$cds_length)
  # count up footprints
  if(class(prop_assign) == "data.frame") {
    rpf_sam_28mer <- within(subset(rpf_sam, length==28 & frame==0),
                            nt_base <- "-")
    rpf_sam_29mer_mismatch <- within(subset(rpf_sam, length==29 & frame==2), 
                                     { wts <- wts * prop_assign$weight[prop_assign$mismatch==T]
                                     nt_base <- substr(rpf, 1, 1)
                                     length <- length - 1
                                     frame <- 0
                                     })
    rpf_sam_29mer_noMismatch <- within(subset(rpf_sam, length==29 & frame==2), 
                                       { wts <- wts * prop_assign$weight[prop_assign$mismatch==F]
                                       nt_base <- "-"
                                       })
    rpf_sam <- rbind(rpf_sam_28mer, rpf_sam_29mer_mismatch, rpf_sam_29mer_noMismatch)
  } else {
    rpf_sam <- rbind(subset(rpf_sam, length==28 & frame==0), subset(rpf_sam, length==29 & frame==2))
    rpf_sam$nt_base <- "-"
    rpf_sam$nt_base[rpf_sam$length==29] <- substr(rpf_sam$rpf[rpf_sam$length==29], 1, 1)
  }
  rpf_sam <- subset(rpf_sam, nt_base != "N")
  rpf_sam <- aggregate(wts ~ cod_idx + nt_base, data=rpf_sam, FUN=sum)
  # add counts to regression data.frame
  regression_data <- subset(regression_data, regression_data$transcript==transcript_name)
  variables <- c("cod_idx", "nt_base")
  rpf_rows <- prodlim::row.match(rpf_sam[, variables], regression_data[, variables])
  # take care of mismatches
  rpf_sam <- subset(rpf_sam, !is.na(rpf_rows))
  rpf_rows <- rpf_rows[!is.na(rpf_rows)]
  # return data.frame with footprint counts
  regression_data$count[rpf_rows] <- rpf_sam$wts
  return(regression_data)
}

compute_29mer_weight <- function(sam_fname, transcript_length_fname, 
                                 which_transcripts, mapping_weights=F, start_codon=F) {
  # compute weighting of 29mer footprint to mismatch versus no mismatch
  ## sam_fname: character; file.path to .sam alignment file
  ## transcript_length_fname: character; file.path to transcriptome lengths file
  ## which_transcripts: character vector; transcript ids in 3rd column of .sam alignment file
  ## mapping_weights: logical; whether .sam alignment file has mapping weights (i.e. output from RSEM)
  utr5_length <- unique(load_lengths(transcript_length_fname)$utr5_length)
  num_cores <- min(parallel::detectCores()-8, length(which_transcripts))
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  dat <- foreach(a=which_transcripts, .export="count_29mer_mismatch_gene") %dopar% {
    count_29mer_mismatch_gene(a, sam_fname, utr5_length, mapping_weights, start_codon)
  }
  dat <- data.frame(mismatch = c(F, T), 
                    weight = rowSums(sapply(dat, function(x) x$weight)))
  dat$weight <- dat$weight / sum(dat$weight)
  parallel::stopCluster(cl)
  return(dat)
}

count_29mer_mismatch_gene <- function(transcript_name, sam_fname, utr5_length=20, mapping_weights=F, start_codon=F) {
  # per gene, count number of 29mers in frame 2 with/out mismatch
  ## transcript_name: character; transcript id in 3rd column of .sam alignment file
  ## sam_fname: character; file.path to .sam alignment file
  ## utr5_length: integer; number of nucleotides padding until start codon
  ## mapping_weights: logical; whether to .sam alignment file has mapping weights (i.e. output from RSEM)
  ## start_codon: logical; whether to only count footprints at start codon (i.e. pos=8)
  if(mapping_weights) {
    sam_dat <- system(paste("grep -w", transcript_name, sam_fname, "| grep -v ^@ | cut -f 4,10,13,15"), intern=T)
    sam_dat <- data.frame(matrix(unlist(strsplit(sam_dat, split="\t")), ncol=4, byrow=T), stringsAsFactors=F)
    colnames(sam_dat) <- c("pos", "rpf_seq", "mismatch", "weight")
    sam_dat$weight <- as.numeric(sub("ZW:f:", "", sam_dat$weight))
  } else {
    sam_dat <- system(paste("grep -w", transcript, sam_fname, "| grep -v ^@ | cut -f 4,10,13"), intern=T)
    sam_dat <- data.frame(matrix(unlist(strsplit(sam_dat, split="\t")), ncol=3, byrow=T), stringsAsFactors=F)
    colnames(sam_dat) <- c("pos", "rpf_seq", "mismatch")
    sam_dat$weight <- 1
  }
  sam_dat$pos <- as.numeric(sam_dat$pos)
  sam_dat$frame <- (sam_dat$pos - utr5_length - 1) %% 3
  sam_dat$rpf_length <- nchar(sam_dat$rpf_seq)
  sam_dat$mismatch <- grepl("MD:Z:0", sam_dat$mismatch)
  sam_dat <- subset(sam_dat, rpf_length==29 & frame==2)
  if(start_codon) { sam_dat <- subset(sam_dat, pos==8) }
  counts <- aggregate(weight ~ mismatch, data=sam_dat, FUN=sum)
  return(counts)
}

count_footprints <- function(sam_fname, transcript_length_fname, offsets_fname, 
                             regression_data, num_cores=NULL, mapping_weights=F, prop_assign=F,
                             start_codon=F) {
  # count up footprints by transcript, A site, and digest lengths
  ## sam_fname: character; file.path to .sam alignment file
  ## transcript_length_fname: character; file.path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## regression_data: data.frame; output from init_data()
  ## num_cores: integer; number of cores to parallelize over
  ## mapping_weights: logical; whether .sam alignment file has mapping weights (i.e. output from RSEM)
  ## prop_assign: logical; whether to use 29mer weights
  ## start_codon: logical; whether to limit 29mer weights to start codon footprints
  num_subsets <- length(levels(regression_data$transcript))
  # count footprints
  if(is.null(num_cores)) {
    num_cores <- parallel::detectCores()-8
  }
  num_cores <- ifelse(num_cores > num_subsets, num_subsets, num_cores)
  if(prop_assign) {
    prop_assign <- compute_29mer_weight(sam_fname, transcript_length_fname, 
                                        which_transcripts=levels(regression_data$transcript),
                                        mapping_weights, start_codon)
  }
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  regression_data <- foreach(a=levels(regression_data$transcript),
                             .export=c("count_footprints_gene", "load_offsets", "load_lengths"),
                             .combine=rbind) %dopar% {
                               tmp_data <- subset(regression_data, transcript==a)
                               count_footprints_gene(a, sam_fname, transcript_length_fname,
                                                     offsets_fname, tmp_data, 
                                                     mapping_weights, prop_assign)
                             }
  parallel::stopCluster(cl)
  return(regression_data)
}

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
