correct_bias <- function(dat, fit) {
  # predict counts if bias sequences were set to reference level
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
  # predict counts if bias sequences were set to reference level
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

count_by_codon <- function(dat) {
  # collapse data to single counts per codon
  ## dat: data.frame with colnames: transcript, cod_idx, count, corrected_count
  require(prodlim)
  count <- aggregate(count ~ transcript + cod_idx, data=dat, FUN=sum)
  corrected <- aggregate(corrected_count ~ transcript + cod_idx, data=dat, FUN=sum)
  cts_by_codon <- unique(dat[,c("transcript", "cod_idx", "A", "P", "E")])
  cts_by_codon$count <- count$count[row.match(cts_by_codon[,c("transcript", "cod_idx")],
                                              count[,c("transcript", "cod_idx")])]
  cts_by_codon$corrected_count <- corrected$corrected_count[row.match(cts_by_codon[,c("transcript", "cod_idx")],
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
  transcriptStartLines <- c(transcriptStartLines, length(rawFile)+1)
  transcriptStartLines <- c(transcriptStartLines, length(rawFile)+1) # add extra line for bookkeeping
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
                     names(codonSequence) <- as.character(seq.int(from=-floor(pad5/3), length.out=nCodons))
                     return(codonSequence)
                   })
  names(faList) <- transcriptNames
  return(faList)
}

evaluate_bias <- function(dat, transcripts_fa_fname, utr5=18, utr3=15) {
  # perform iXnos regression and generate leave-one-out correlation plots
  ## dat: data.frame containing regression predictors and column of corrected counts
  ## transcripts_fa_fname: character; filepath to transcriptome .fa file
  # 1. aggregate counts by codon
  cts_by_codon <- count_by_codon(dat)
  # 2. scale footprint counts by transcript
  transcripts <- levels(dat$transcript)
  transcripts_count_mean <- sapply(transcripts, 
                                   function(x) {
                                     mean(subset(cts_by_codon, transcript==x)$count)
                                   })
  cts_by_codon$count <- cts_by_codon$count / transcripts_count_mean[as.character(cts_by_codon$transcript)]
  transcripts_corrected_mean <- sapply(transcripts, 
                                       function(x) {
                                         mean(subset(cts_by_codon, transcript==x)$corrected_count)
                                       })
  cts_by_codon$corrected_count <- cts_by_codon$corrected_count / transcripts_corrected_mean[as.character(cts_by_codon$transcript)]
  # 3. make data.frame for iXnos regression
  transcripts_seq <- readFAfile(transcripts_fa_fname, utr5, utr3)
  codons <- data.frame(t(sapply(1:nrow(cts_by_codon),
                                function(x) {
                                  transcript <- as.character(cts_by_codon$transcript[x])
                                  cod_idx <- cts_by_codon$cod_idx[x]
                                  Asite_index <- which(names(transcripts_seq[[transcript]])==as.character(cod_idx-1))
                                  codons <- transcripts_seq[[transcript]][(Asite_index-(floor(utr5/3))):(Asite_index+(floor(utr3/3)))]
                                  return(codons)
                                })))
  count_dat <- data.frame(count = cts_by_codon$count, codons)
  corrected_dat <- data.frame(count=cts_by_codon$corrected_count, codons)
  # 4. full model
  count_full_cor <- cor(count_dat$count,
                        predict(lm(count ~ ., data=count_dat)))
  corrected_full_cor <- cor(corrected_dat$count,
                            predict(lm(count ~ ., data=corrected_dat)))
  # 5. leave-one-out models
  count_loo_cor <- sapply(2:ncol(count_dat),
                          function(x) {
                            cor(count_dat$count,
                                predict(lm(count ~ ., data=count_dat[,-x])))
                          })
  corrected_loo_cor <- sapply(2:ncol(corrected_dat),
                              function(x) {
                                cor(corrected_dat$count,
                                    predict(lm(count ~ ., data=corrected_dat[,-x])))
                              })
  # 6. return model correlations
  model_cor <- data.frame(uncorrected = c(count_full_cor, count_loo_cor),
                          corrected = c(corrected_full_cor, corrected_loo_cor))
  rownames(model_cor) <- c("full", as.character((-floor(utr5/3)):(-3)), "E", "P", "A", as.character(1:(floor(utr3/3))))
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
