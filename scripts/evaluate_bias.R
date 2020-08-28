### requires libraries: prodlim, ggplot2

library(ggplot2)

evaluate_bias <- function(dat, which_column="count", transcripts_fa_fname, 
                          utr5=18, utr3=15, trunc5=20, trunc3=20, 
                          num_f5_codons=6, num_f3_codons=6, type="codon") {
  # perform iXnos regression and generate leave-one-out correlation plots
  ## dat: data.frame containing regression predictors and column of corrected counts
  ## which_column: character; name of column containing counts
  ## transcripts_fa_fname: character; filepath to transcriptome .fa file
  ## utr5: integer; 5' UTR padding in transcriptome .fa file
  ## utr3: integer; 3' UTR padding in transcriptome .fa file
  ## trunc5: integer; number of 5' codons to leave out in codon correlation regression
  ## trunc3: integer; number of 3' codons to leave out in codon correlation regression
  ## num_f5_codons: integer; number of codons to include 5' of A site in regression
  ## num_f3_codons: integer; number of codons to include 3' of A site in regression
  ## type: character; one of "codon" or "nt"
  # 1. aggregate counts by codon
  cts_by_codon <- aggregate(formula(paste(which_column, "~ transcript + cod_idx")), 
                            data=dat, FUN=sum, na.rm=T)
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
                                  codons <- transcripts_seq[[transcript]][(Asite_index-num_f5_codons):(Asite_index+num_f3_codons)]
                                  return(codons)
                                })))
  colnames(codons) <- c(paste0("n", num_f5_codons:3), 
                        "E", "P", "A", 
                        paste0("p", 1:num_f3_codons))
  if(type=="nt") {
    codons <- data.frame(lapply(codons, as.character), stringsAsFactors=F)
    codons <- sapply(seq(nrow(codons)), function(x) {paste0(codons[x,], collapse="")})
    codons <- data.frame(matrix(unlist(strsplit(codons, split="")), 
                                ncol=3*(num_f5_codons+num_f3_codons+1), byrow=T))
    colnames(codons) <- c(paste0("n", (3*num_f3_codons):7),
                          paste0(rep(c("E", "P", "A"), each=3), 0:2),
                          paste0("p", 1:(3*num_f5_codons)))
  }
  count_dat <- data.frame(count = cts_by_codon$count, codons)
  # 4. full model
  count_full_cor <- cor(count_dat$count, predict(lm(count ~ ., data=count_dat)))
  # 5. leave-one-out models
  count_loo_cor <- sapply(colnames(codons),
                          function(x) {
                            cor(count_dat$count,
                                predict(lm(count ~ ., data=count_dat[,-which(colnames(count_dat)==x)])))
                          })
  # 6. return model correlations
  model_cor <- c(count_full_cor, count_loo_cor)
  names(model_cor) <- c("full", colnames(codons))
  return(model_cor)
}

plot_bias <- function(model_cor, plot_title="", plot_subtitle="", type="codon") {
  # make iXnos codon correlation plots
  ## model_cor: numeric vector; output from evaluate_bias()
  ## plot_title: character; title for output plot
  ## plot_subtitle: character; subtitle for output plot
  ## type: character; one of "codon" or "nt"
  names(model_cor) <- sub("n", "-", rownames(model_cor))
  names(model_cor) <- sub("p", "", rownames(model_cor))
  cor_diff <- data.frame(diff=(model_cor - model_cor[1])[-1],
                         site=names(model_cor)[-1])
  cor_diff$site <- factor(cor_diff$site, levels=names(model_cor)[-1])
  bias_plot <- ggplot(cor_diff, aes(x=site, y=diff)) + geom_bar(stat="identity") + 
    theme_bw() + xlab("position") + ylab(expression(Delta*" correlation")) + ggtitle(plot_title)
  cor_diff <- data.frame(diff = c(model_cor$uncorrected[1] - model_cor$uncorrected[2:nrow(model_cor)],
                                  model_cor$corrected[1] - model_cor$corrected[2:nrow(model_cor)]),
                         site = rownames(model_cor)[-1],
                         type = rep(colnames(model_cor), each=nrow(model_cor)-1))
  cor_diff$site <- factor(cor_diff$site, levels=rownames(model_cor)[-1])
  cor_diff$type <- factor(cor_diff$type, levels=colnames(model_cor))
  bias_plot <- ggplot(cor_diff, aes(x=site, y=diff)) + geom_bar(stat="identity") + facet_wrap(~type) +
    theme_bw() + xlab("position") + ylab(expression(Delta*" correlation")) + 
    ggtitle(plot_title, subtitle=plot_subtitle)
  if(type=="nt") {
    bias_plot <- bias_plot + theme(axis.text.x=element_text(angle=90, vjust=0.5))
  }
  ### TODO: add color to bars according to position
  return(bias_plot)
}