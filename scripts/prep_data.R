## requires libraries: parallel, doParallel, foreach, reshape, prodlim

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

init_data <- function(transcript_fa_fname, transcript_length_fname,
                      digest5_lengths=14:17, digest3_lengths=9:11, num_cores=NULL) {
  # initialize data.frame for downstream GLM
  ## transcript_fa_fname: character; file path to transcriptome .fa file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## digest5_lengths: integer vector; legal 5' digest lengths
  ## digest3_lengths: integer vector; legal 3' digest lengths
  transcript_seq <- load_fa(transcript_fa_fname)
  transcript_length <- load_lengths(transcript_length_fname)
  transcript <- unlist(mapply(rep, x=transcript_length$transcript, times=transcript_length$cds_length/3))
  cod_idx <- unlist(lapply(transcript_length$cds_length/3, seq))
  if(is.null(num_cores)) {
    num_cores <- parallel::detectCores()-8
  }
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  codons <- foreach(a=transcript, b=cod_idx, 
                    c=transcript_length$utr5_len[match(transcript, transcript_length$transcript)],
                    .combine='rbind', .export=c("get_codons")) %dopar% get_codons(a, b, c, transcript_seq)
  parallel::stopCluster(cl)
  dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons),
                                 expand.grid(d5=as.factor(digest5_lengths), d3=as.factor(digest3_lengths)))
  dat$count <- 0
  return(dat)
}

load_offsets <- function(offsets_fname) {
  # load A site offset rules
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## rownames: footprint length
  ## colnames: frame (0, 1, 2)
  offsets <- read.table(offsets_fname)
  offsets <- data.frame(frame=as.vector(mapply(rep, 0:2, nrow(offsets))),
                        length=rep(as.numeric(rownames(offsets)), 3),
                        offset=c(offsets$frame_0, offsets$frame_1, offsets$frame_2))
  return(offsets)
}

count_footprints_subset <- function(start_line, num_lines, 
                                    sam_fname, transcript_length_fname, offsets_fname, regression_data) {
  # count up footprints by transcript, A site, and digest lengths
  ## start_line: integer; first line of .sam file to load
  ## num_lines: integer; number of rows of .sam file to load
  ## sam_fname: character; file.path to .sam alignment file
  ## transcript_length_fname: character; file path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## regression_data: data.frame; output from init_data()
  colClasses <- c(rep("NULL", 2), "character", "integer", rep("NULL", 5), "character", rep("NULL", 4))
  rpf_sam <- read.table(sam_fname, skip=start_line-1, stringsAsFactors=F, 
                        colClasses=colClasses, nrows=num_lines)
  colnames(rpf_sam) <- c("transcript", "pos", "rpf")
  # calculate cod_idx and digest lengths
  rpf_sam$frame <- (rpf_sam$pos-1) %% 3
  rpf_sam$length <- nchar(rpf_sam$rpf)
  offsets <- load_offsets(offsets_fname)
  rpf_sam$d5 <- offsets$offset[prodlim::row.match(rpf_sam[,c("frame", "length")], 
                                                  offsets[c("frame", "length")])]
  rpf_sam$d3 <- with(rpf_sam, length-d5-3)
  transcript_length <- load_lengths(transcript_length_fname)
  rpf_sam$utr5_length <- transcript_length$utr5_length[match(rpf_sam$transcript, transcript_length$transcript)]
  rpf_sam$cod_idx <- with(rpf_sam, (pos+d5-utr5_length+2)/3)
  # count up footprints with same (transcript, cod_idx, d5, d3)
  variables <- c("transcript", "cod_idx", "d5", "d3")
  rpf_sam <- rpf_sam[, variables]
  rpf_sam <- aggregate(list(count=rep(1, nrow(rpf_sam))), rpf_sam, length)
  # add counts to regression data.frame
  rpf_rows <- prodlim::row.match(rpf_sam[, variables], regression_data[, variables])
  regression_data$count[rpf_rows] <- rpf_sam$count
  return(regression_data$count)
}

count_footprints <- function(sam_fname, transcript_length_fname, offsets_fname, 
                             regression_data, subset_size=1e5, num_cores=NULL) {
  # count up footprints by transcript, A site, and digest lengths
  ## sam_fname: character; file.path to .sam alignment file
  ## transcript_length_fname: character; file.path to transcriptome lengths file
  ## offsets_fname: character; file.path to offset / A site assignment rules .txt file
  ## regression_data: data.frame; output from init_data()
  # establish start_lines and num_lines for count_footprints_subset()
  n_sam_header_lines <- length(system(paste("grep @", sam_fname), intern=T))
  n_footprints <- as.numeric(strsplit(system(paste("wc -l", sam_fname), intern=T), split=" ")[[1]][1])
  n_footprints <- n_footprints - n_sam_header_lines
  num_subsets <- floor(n_footprints / subset_size)
  start_lines <- (1:num_subsets-1)*subset_size + 1 + n_sam_header_lines
  num_lines <- rep(subset_size, num_subsets)
  num_lines[num_subsets] <- (n_sam_header_lines + n_footprints) - start_lines[num_subsets] + 1
  # count footprints
  if(is.null(num_cores)) {
    num_cores <- parallel::detectCores()-8
  }
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  fp_counts <- foreach(a=start_lines, b=num_lines, 
                       .export=c("count_footprints_subset", "load_offsets", "load_lengths"),
                       .combine=cbind) %dopar% {
                         count_footprints_subset(a, b, sam_fname, transcript_length_fname, 
                                                 offsets_fname, regression_data)
                       }
  parallel::stopCluster(cl)
  fp_counts <- rowSums(fp_counts)
  regression_data$count <- fp_counts
  return(regression_data)
}

### testing functions
library(foreach)

transcript_fa_fname <- "../sample_data/riboViz_sampleData_transcripts_18cds18.fa"
transcript_length_fname <- "../sample_data/riboViz_sampleData_transcripts_18cds18_lengths.txt"
offsets_fname <- "../sample_data/Asite_rules.txt"
sam_fname <- "../sample_data/riboViz_sampleData_footprints.sam"

regression_data <- init_data(transcript_fa_fname, transcript_length_fname)
system.time({regression_data <- count_footprints(sam_fname, transcript_length_fname, offsets_fname, regression_data)})
save(regression_data, file="../sample_data/riboViz_sampleData_regression.Rda")
