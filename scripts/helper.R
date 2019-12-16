## requires libraries: parallel, doParallel, foreach, reshape

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

init_dat <- function(transcript_fa_fname, transcript_length_fname,
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
  stopCluster(cl)
  dat <- reshape::expand.grid.df(data.frame(transcript, cod_idx, codons),
                                 expand.grid(d5=as.factor(digest5_lengths), d3=as.factor(digest3_lengths)))
  dat$count <- 0
  return(dat)
}