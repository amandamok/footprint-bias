library(foreach)

correct_bias <- function(dat, fit) {
  # predict counts if bias sequences were set to reference level, allowing residuals
  ## dat: data.frame containing regression predictors
  ## fit: fitted model object
  # 1. establish f5 scaling factors
  f5_coefs <- coef(fit)[match(paste0("f5", levels(dat$f5)), names(coef(fit)))]
  names(f5_coefs) <- levels(dat$f5)
  f5_coefs[is.na(f5_coefs)] <- 0
  f5_coefs <- exp(f5_coefs)
  # 2. establish f3 scaling factors
  f3_coefs <- coef(fit)[match(paste0("f3", levels(dat$f3)), names(coef(fit)))]
  names(f3_coefs) <- levels(dat$f3)
  f3_coefs[is.na(f3_coefs)] <- 0
  f3_coefs <- exp(f3_coefs)
  # 3. calculate corrected counts
  dat$corrected_count <- dat$count / (f5_coefs[as.character(dat$f5)] * f3_coefs[as.character(dat$f3)])
  # 4. rescale predicted counts so they sum to original footprint count
  dat$corrected_count <- dat$corrected_count * sum(dat$count) / sum(dat$corrected_count)
  # 5. return dat with corrected counts
  return(dat)
}

correct_bias_predict <- function(dat, fit) {
  # predict counts if bias sequences were set to reference level, using model fit
  ## dat: data.frame containing regression predictors
  ## fit: fitted model object
  # 1. set all bias sequences to reference level
  dat_corrected <- dat
  dat_corrected$f5 <- factor(levels(dat_corrected$f5)[1], levels=levels(dat_corrected$f5))
  dat_corrected$f3 <- factor(levels(dat_corrected$f3)[1], levels=levels(dat_corrected$f3))
  # 2. predict counts for modified data, add to original data.frame
  dat$corrected_count <- predict(fit, newdata=dat_corrected, type="response")
  # 3. rescale predicted counts so they sum to original footprint count
  dat$corrected_count <- dat$corrected_count * sum(dat$count) / sum(dat$corrected_count)
  # 4. return dat with corrected counts
  return(dat)
}

correct_bias_sim <- function(dat, p5bias, n3bias) {
  # predict counts if bias sequences were set to reference level, using simulation parameters
  ## dat: data.frame containing regression predictors
  ## p5bias: named numeric vector; recovery probabilities from simulation
  ## n3bias: named numeric vector; recovery probabilities from simulation
  # 1. establish f5 scaling factors
  f5_ref_level <- levels(dat$f5)[1]
  f5_coefs <- p5bias / p5bias[f5_ref_level]
  f5_coefs[f5_coefs==0] <- 1
  # 2. establish f3 scaling factors
  f3_ref_level <- levels(dat$f3)[1]
  f3_coefs <- n3bias / n3bias[f3_ref_level]
  f3_coefs[f3_coefs==0] <- 1
  # 3. calculate corrected counts
  dat$corrected_count <- dat$count / (f5_coefs[as.character(dat$f5)] * f3_coefs[as.character(dat$f3)])
  # 4. rescale predicted counts so they sum to original footprint count
  dat$corrected_count <- dat$corrected_count * sum(dat$count) / sum(dat$corrected_count)
  # 5. return dat with corrected counts
  return(dat)
}

correct_bias_sam_chunk <- function(sam_fname, p5bias, n3bias, start_line, end_line, replace_weights=F, 
                                   tmp_dir, tmp_fname) {
  # add column to .sam alignment file (subset)
  ## sam_fname: character; file.path to .sam alignment file
  ## p5bias: named numeric vector; fold change values, relative to reference level
  ## n3bias: named numeric vector; fold change values, relative to reference level
  ## start_line: integer; first line in sam_fname to read
  ## end_line: integer; last line in sam_fname to read 
  ## replace_weights: logical; whether .sam alignment has existing weights (RSEM)
  ## tmp_dir: character; file.path to directory for temporary files
  ## tmp_fname: character; file name for output .sam file (subset)
  # extract footprint sequences and weights
  subset_command <- paste("head -n", as.integer(end_line), sam_fname, "| tail -n", as.integer(end_line-start_line+1))
  fp_seq <- system(paste(subset_command, "| cut -f 10"), intern=T)
  if(replace_weights==T) {
    fp_weights <- system(paste(subset_command, "| cut -f 15"), intern=T)
    fp_weights <- as.numeric(sub("ZW:f:", "", fp_weights))
  } else {
    fp_weights <- rep(1, length(fp_seq))
  }
  # compute correction weights
  f5_length <- nchar(names(p5bias)[1])
  f5_seq <- substr(fp_seq, start=1, stop=f5_length)
  f5_weights <- p5bias[match(f5_seq, names(p5bias))]
  f5_weights[is.na(f5_weights)] <- 1 # reads with f5 sequence not in names(p5bias)
  f3_length <- nchar(names(n3bias)[1])
  f3_seq <- substr(fp_seq, start=nchar(fp_seq)-f3_length+1, stop=nchar(fp_seq))
  f3_weights <- n3bias[match(f3_seq, names(n3bias))]
  f3_weights[is.na(f3_weights)] <- 1 # reads with f3 sequence not in names(n3bias)
  # correct counts
  fp_weights <- fp_weights / (f5_weights * f3_weights)
  # append column to .sam file
  system(paste(subset_command, ">", file.path(tmp_dir, paste0(tmp_fname, ".tmp"))))
  write.table(fp_weights, file=file.path(tmp_dir, paste0(tmp_fname, ".wts")), 
              quote=F, row.names=F, col.names=F)
  system(paste("paste", file.path(tmp_dir, paste0(tmp_fname, ".tmp")), 
               file.path(tmp_dir, paste0(tmp_fname, ".wts")), ">", file.path(tmp_dir, tmp_fname)))
  system(paste("rm", file.path(tmp_dir, paste0(tmp_fname, ".tmp")), file.path(tmp_dir, paste0(tmp_fname, ".wts"))))
}

correct_bias_sam <- function(sam_fname, p5bias, n3bias, replace_weights=F, tmp_dir, new_fname,
                             chunk_size=1e6, num_cores=NULL) {
  # add column to .sam alignment file
  ## sam_fname: character; file.path to .sam alignment file
  ## p5bias: named numeric vector; fold change values, relative to reference level
  ## n3bias: named numeric vector; fold change values, relative to reference level
  ## replace_weights: logical; whether .sam alignment has existing weights (RSEM)
  ## tmp_dir: character; file.path to directory for temporary files
  ## new_fname: character; file.path for corrected .sam file output
  ## chunk_size: integer; number of lines of .sam file to process at a time
  ## num_cores: integer; number of cores to parallelize over
  # create temporary directory
  if(!dir.exists(tmp_dir)) {dir.create(tmp_dir)}
  # define chunks to be processed by correct_bias_sam_chunk()
  num_lines <- as.numeric(strsplit(system(paste("wc -l", sam_fname), intern=T), split=" ")[[1]][1])
  num_header_lines <- as.numeric(strsplit(system(paste("grep ^@", sam_fname, "| wc -l"), intern=T), split=" ")[[1]][1])
  chunk_starts <- seq(from=num_header_lines+1, to=num_lines, by=chunk_size)
  chunk_ends <- chunk_starts + chunk_size - 1
  chunk_ends[length(chunk_ends)] <- num_lines
  chunk_names <- paste0("chunk_", seq.int(length(chunk_starts)), ".sam")
  # parallelize processing
  if(is.null(num_cores)) {
    num_cores <- parallel::detectCores()-8
  }
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  foreach(tmp_start=chunk_starts, tmp_end=chunk_ends, tmp_fname=chunk_names,
          .export=c("correct_bias_sam_chunk"), .inorder=F, .errorhandling='pass') %dopar% {
            correct_bias_sam_chunk(sam_fname, p5bias, n3bias, tmp_start, tmp_end, replace_weights,
                                   tmp_dir, tmp_fname)
            return(tmp_start)
          }
  parallel::stopCluster(cl)
  # merge .sam header and chunks
  system(paste("head -n", num_header_lines, sam_fname, ">", file.path(tmp_dir, "header.sam")))
  system(paste("cat", file.path(tmp_dir, "header.sam"), 
               paste0(sapply(chunk_names, function(x) file.path(tmp_dir, x)), collapse=" "), ">", new_fname))
  corrected_num_lines <- as.numeric(strsplit(system(paste("wc -l", new_fname), intern=T), split=" ")[[1]][1])
  if(corrected_num_lines != num_lines) {print("LINES MISSING")}
  system(paste("rm -r", tmp_dir))
}

count_by_codon <- function(dat) {
  # collapse data to single counts per codon
  ## dat: data.frame with colnames: transcript, cod_idx, count, corrected_count
  count <- aggregate(count ~ transcript + cod_idx, data=dat, FUN=sum)
  corrected <- aggregate(corrected_count ~ transcript + cod_idx, data=dat, FUN=sum)
  cts_by_codon <- unique(dat[,c("transcript", "cod_idx", "A", "P", "E")])
  cts_by_codon$count <- count$count[prodlim::row.match(cts_by_codon[,c("transcript", "cod_idx")],
                                                       count[,c("transcript", "cod_idx")])]
  cts_by_codon$corrected_count <- corrected$corrected_count[prodlim::row.match(cts_by_codon[,c("transcript", "cod_idx")],
                                                                               corrected[,c("transcript", "cod_idx")])]
  return(cts_by_codon)
}

readFAfile <- function(faFile, pad5, pad3) {
  ## read in .fa file where sequence broken up over multiple lines, convert to list of vectors of codons per transcript
  # faFile: character; path to .fa file of transcript sequences
  # pad5: numeric; number of nt padding 5' end of cds
  # pad3: numeric; number of nt padding 3' end of cds
  rawFile <- readLines(faFile)
  transcriptStartLines <- grep(">", rawFile)
  nTranscripts <- length(transcriptStartLines)
  transcriptNames <- sapply(transcriptStartLines,
                            function(x) {
                              gsub(">", "", strsplit(rawFile[x], split=" ")[[1]][1])
                            })
  transcriptStartLines <- c(transcriptStartLines, length(rawFile)+1, length(rawFile)+1) # add extra line for bookkeeping
  faList <- sapply(1:nTranscripts,
                   function(x) {
                     startLine <- transcriptStartLines[x]+1
                     endLine <- transcriptStartLines[x+1]-1
                     transcriptSequence <- paste(rawFile[startLine:endLine], collapse="")
                     nCodons <- floor((nchar(transcriptSequence))/3)
                     sequenceOffset <- pad5 %% 3
                     codonSequence <- substring(transcriptSequence, 
                                                first=(3*(1:nCodons-1)+1)+sequenceOffset, 
                                                last=(3*(1:nCodons))+sequenceOffset)
                     names(codonSequence) <- as.character(seq.int(from=-floor(pad5/3), length.out=nCodons)) # start codon is 0
                     return(codonSequence)
                   })
  names(faList) <- transcriptNames
  return(faList)
}

evaluate_bias <- function(dat, transcripts_fa_fname, utr5=18, utr3=15, trunc5=20, trunc3=20) {
  # perform iXnos regression and generate leave-one-out correlation plots
  ## dat: data.frame containing regression predictors and column of corrected counts
  ## transcripts_fa_fname: character; filepath to transcriptome .fa file
  ## utr5: integer; 5' UTR padding in transcriptome .fa file
  ## utr3: integer; 3' UTR padding in transcriptome .fa file
  ## trunc5: integer; number of 5' codons to leave out in codon correlation regression
  ## trunc3: integer; number of 3' codons to leave out in codon correlation regression
  # 1. aggregate counts by codon
  cts_by_codon <- count_by_codon(dat)
  # 2. remove truncated codons, scale footprint counts by transcript mean
  transcripts <- levels(dat$transcript)
  cts_by_codon_trunc <- lapply(transcripts,
                               function(x) {
                                 tmp_subset <- subset(cts_by_codon, transcript==x)
                                 # remove trunc5 and trunc3 codons
                                 num_codons <- max(tmp_subset$cod_idx)
                                 tmp_subset <- subset(tmp_subset,
                                                      cod_idx > trunc5 & cod_idx <= (num_codons-trunc3))
                                 # scale per-codon counts by transcript mean
                                 tmp_subset$count <- tmp_subset$count / mean(tmp_subset$count)
                                 tmp_subset$corrected_count <- tmp_subset$corrected_count / mean(tmp_subset$corrected_count)
                                 # return truncated & scaled counts
                                 return(tmp_subset)
                               })
  cts_by_codon <- do.call(rbind, cts_by_codon_trunc)
  # 3. make data.frame for iXnos regression
  transcripts_seq <- readFAfile(transcripts_fa_fname, utr5, utr3)
  codons <- data.frame(t(sapply(1:nrow(cts_by_codon),
                                function(x) {
                                  transcript <- as.character(cts_by_codon$transcript[x])
                                  cod_idx <- cts_by_codon$cod_idx[x] - 1 # transcripts_seq uses 0-based indexing
                                  Asite_index <- which(names(transcripts_seq[[transcript]])==as.character(cod_idx)) 
                                  codons <- transcripts_seq[[transcript]][(Asite_index-(floor(utr5/3))):(Asite_index+(floor(utr3/3)))]
                                  return(codons)
                                })))
  codon_positions <- c(paste0("p", floor(utr5/3):(3)), "E", "P", "A", paste0("n", 1:floor(utr3/3)))
  colnames(codons) <- codon_positions
  count_dat <- data.frame(count = cts_by_codon$count, codons)
  corrected_dat <- data.frame(count=cts_by_codon$corrected_count, codons)
  # 4. full model
  count_full_cor <- cor(count_dat$count,
                        predict(lm(count ~ ., data=count_dat)))
  corrected_full_cor <- cor(corrected_dat$count,
                            predict(lm(count ~ ., data=corrected_dat)))
  # 5. leave-one-out models
  count_loo_cor <- sapply(codon_positions,
                          function(x) {
                            cor(count_dat$count,
                                predict(lm(count ~ ., data=count_dat[,-which(colnames(count_dat)==x)])))
                          })
  corrected_loo_cor <- sapply(codon_positions,
                              function(x) {
                                cor(corrected_dat$count,
                                    predict(lm(count ~ ., data=corrected_dat[,-which(colnames(corrected_dat)==x)])))
                              })
  # 6. return model correlations
  model_cor <- data.frame(uncorrected = c(count_full_cor, count_loo_cor),
                          corrected = c(corrected_full_cor, corrected_loo_cor))
  rownames(model_cor) <- c("full", codon_positions)
  return(model_cor)
}

plot_bias <- function(model_cor, plotTitle="") {
  # make iXnos codon correlation plots
  ## model_cor: data.frame; output from evaluate_bias()
  ## plotTitle: character; 
  cor_diff <- data.frame(diff = c(model_cor$uncorrected[1] - model_cor$uncorrected[2:nrow(model_cor)],
                                  model_cor$corrected[1] - model_cor$corrected[2:nrow(model_cor)]),
                         site = rownames(model_cor)[-1],
                         type = rep(colnames(model_cor), each=nrow(model_cor)-1))
  cor_diff$site <- factor(cor_diff$site, levels=rownames(model_cor)[-1])
  cor_diff$type <- factor(cor_diff$type, levels=colnames(model_cor))
  ggplot(cor_diff, aes(x=site, y=diff)) + geom_bar(stat="identity") + facet_wrap(~type) +
    theme_bw() + xlab("position") + ylab(expression(Delta*" correlation")) + ggtitle(plotTitle)
}

codon_corr <- function(regression_data, model_fit, correction, transcripts_fa_fname,
                       utr5=20, utr3=20, trunc5=20, trunc3=20, p5bias=NULL, n3bias=NULL, plot_title="") {
  # wrapper function to correct data from regression fit, compute and plot codon correlations
  ## regression_data: data.frame containing regression predictors; output from init_data() and count_footprints()
  ## model_fit: fitted model object
  ## correction: character; which correction function to use
  ### one of "correct_bias", "correct_bias_predict", or "correct_bias_sim"
  ### correct_bias: takes regression coefficient from model_fit, applies correction coefficient as exp(beta)
  ### correct_bias_predict: takes prediction from model_fit (i.e. ignore residuals)
  ### correct_bias_sim: takes simulation probabilities as correction coefficients
  ## transcripts_fa_fname: character; file.path to transcriptome .fa file
  ## utr5: integer; 5' UTR padding in transcriptome .fa file
  ## utr3: integer; 3' UTR padding in transcriptome .fa file
  ## trunc5: integer; number of 5' codons to leave out in codon correlation regression
  ## trunc3: integer; number of 3' codons to leave out in codon correlation regression
  ## p5bias: named numeric vector; recovery probabilities from simulation for correct_bias_sim()
  ## n3bias: named numeric vector; recovery probabilities from simulation for correct_bias_sim()
  ## plot_title: character; title for codon correlation plot
  # 1. apply correction to regression_data
  if(correction %in% c("correct_bias", "correct_bias_predict")) {
    dat <- do.call(correction, args=list(dat=regression_data, fit=model_fit))
  } else {
    dat <- correct_bias_sim(regression_data, p5bias, n3bias)
  }
  # 2. evaluate bias
  loo_cors <- evaluate_bias(dat, transcripts_fa_fname, utr5, utr3, trunc5, trunc3)
  # 3. plot biases
  loo_plot <- plot_bias(loo_cors, plot_title)
  # 4. return model correlations and plot
  return(list(codon_corrs=loo_cors, codon_plot=loo_plot))
}
