rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "grch38.transcripts.fa")
transcript_length_fname <- file.path(ref_dir, "grch38.transcripts.lengths.tsv")
transcript_lengths <- load_lengths(transcript_length_fname)
offsets_fname <- file.path(ref_dir, "Asite_rules.txt")
mapping_fname <- file.path(ref_dir, "grch38.transcripts.mappings.tsv")
mapping <- read.table(mapping_fname, header=T, stringsAsFactors=F)
transcript_codons_fname <- file.path(ref_dir, "grch38.transcripts.codonCounts.tsv")
transcript_codons <- read.table(transcript_codons_fname)

results_dir <- file.path(here(), "expts", "goodarzi_2016", "CN34_vs_CN-LM1a")
raw_data_dir <- file.path(here(), "expts", "goodarzi_2016", "raw_data")

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- 100
exclude_codon5 <- 10
exclude_codon3 <- 10
minimum_cds_length <- exclude_codon5 + exclude_codon3 + 10

sra_CN34 <- c("SRR3129936", "SRR3129937")
sra_LM1a <- c("SRR3129940", "SRR3129941")

orig_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)
interxn_model <- formula(count ~ transcript*expt + A*expt + P + E + d5*genome_f5 + d3*genome_f3)

# 1. load footprint alignments --------------------------------------------

# CN34
CN34_bam_fname <- file.path(results_dir, "CN34.bam")
if(!file.exists(CN34_bam_fname)) {
  system(paste("samtools merge", CN34_bam_fname,
               file.path(raw_data_dir, sra_CN34, paste0(sra_CN34, ".transcript.bam"))))
}
CN34_bam_dat_fname <- file.path(results_dir, "CN34_bam.Rda")
if(!file.exists(CN34_bam_dat_fname)) {
  CN34_bam <- load_bam(CN34_bam_fname, transcript_fa_fname,
                      transcript_length_fname, offsets_fname,
                      f5_length=f5_length, f3_length=f3_length)
  save(CN34_bam, file=CN34_bam_dat_fname)
} else {
  load(CN34_bam_dat_fname)
}

# LM1a
LM1a_bam_fname <- file.path(results_dir, "LM1a.bam")
if(!file.exists(LM1a_bam_fname)) {
  system(paste("samtools merge", LM1a_bam_fname,
               file.path(raw_data_dir, sra_LM1a, paste0(sra_LM1a, ".transcript.bam"))))
}
LM1a_bam_dat_fname <- file.path(results_dir, "LM1a_bam.Rda")
if(!file.exists(LM1a_bam_dat_fname)) {
  LM1a_bam <- load_bam(LM1a_bam_fname, transcript_fa_fname,
                      transcript_length_fname, offsets_fname,
                      f5_length=f5_length, f3_length=f3_length)
  save(LM1a_bam, file=LM1a_bam_dat_fname)
} else {
  load(LM1a_bam_dat_fname)
}

# 2. compute size/frame subsets -------------------------------------------

CN34_subsets_fname <- file.path(results_dir, "CN34_subsets.Rda")
if(!file.exists(CN34_subsets_fname)) {
  CN34_d5_d3 <- count_d5_d3(CN34_bam)
  CN34_d5_d3_subsets <- CN34_d5_d3$counts[1:(which(CN34_d5_d3$counts$proportion>min_prop)[1]),
                                        c("d5", "d3")]
  CN34_subset_names <- sapply(seq(nrow(CN34_d5_d3_subsets)),
                             function(x) {
                               paste("d5", CN34_d5_d3_subsets$d5[x],
                                     "d3", CN34_d5_d3_subsets$d3[x],
                                     sep="_")
                             })
  save(CN34_d5_d3, CN34_d5_d3_subsets, CN34_subset_names,
       file=CN34_subsets_fname)
} else {
  load(CN34_subsets_fname)
}

LM1a_subsets_fname <- file.path(results_dir, "LM1a_subsets.Rda")
if(!file.exists(LM1a_subsets_fname)) {
  LM1a_d5_d3 <- count_d5_d3(LM1a_bam)
  LM1a_d5_d3_subsets <- LM1a_d5_d3$counts[1:(which(LM1a_d5_d3$counts$proportion>min_prop)[1]),
                                        c("d5", "d3")]
  LM1a_subset_names <- sapply(seq(nrow(LM1a_d5_d3_subsets)),
                             function(x) {
                               paste("d5", LM1a_d5_d3_subsets$d5[x],
                                     "d3", LM1a_d5_d3_subsets$d3[x],
                                     sep="_")
                             })
  save(LM1a_d5_d3, LM1a_d5_d3_subsets, LM1a_subset_names,
       file=LM1a_subsets_fname)
} else {
  load(LM1a_subsets_fname)
}

interxn_subsets <- unique(rbind(CN34_d5_d3_subsets, LM1a_d5_d3_subsets))
save(interxn_subsets,
     file=file.path(results_dir, "interxn_subsets.Rda"))

# 3. pick transcripts for training ----------------------------------------

training_set_fname <- file.path(results_dir, "training_transcripts.txt")
if(!file.exists(training_set_fname)) {
  # calculate
  ## transcript count
  CN34_ct <- aggregate(count ~ transcript, data=CN34_bam, FUN=sum)
  mapping$CN34_ct <- CN34_ct$count[match(mapping$tx_name, CN34_ct$transcript)]
  LM1a_ct <- aggregate(count ~ transcript, data=LM1a_bam, FUN=sum)
  mapping$LM1a_ct <- LM1a_ct$count[match(mapping$tx_name, LM1a_ct$transcript)]
  ## median codon density
  CN34_density <- calculate_transcript_density(CN34_bam, transcript_length_fname,
                                              statistic="median")
  mapping$CN34_median <- CN34_density[match(mapping$tx_name, names(CN34_density))]
  LM1a_density <- calculate_transcript_density(LM1a_bam, transcript_length_fname,
                                              statistic="median")
  mapping$LM1a_median <- LM1a_density[match(mapping$tx_name, names(LM1a_density))]
  ## average of median codon density
  mapping[is.na(mapping)] <- 0
  mapping$median <- with(mapping, (CN34_median + LM1a_median)/2)
  # filter transcripts
  ## transcript counts > 0 for both conditions
  training_set <- subset(mapping, (CN34_ct > 0) & (LM1a_ct > 0)) # n=23,190
  ## only 1 stop codon
  num_stop_codons <- transcript_codons[match(training_set$tx_name,
                                             rownames(transcript_codons)),]
  num_stop_codons <- rowSums(num_stop_codons[, c("TAG", "TAA", "TGA")])
  training_set <- subset(training_set, num_stop_codons == 1) # n=21,359
  ## sufficient number of codons
  num_codons <- transcript_lengths$cds_length[match(training_set$tx_name,
                                                    transcript_lengths$transcript)] / 3
  training_set <- subset(training_set, num_codons > minimum_cds_length) # n=21,324
  # pick top representative transcript per gene: maximum average of median codon density
  training_set <- split(training_set, training_set$gene_id) # n=11,227
  training_set <- lapply(training_set,
                         function(x) {
                           x[which.max(x$median),]
                         })
  training_set <- do.call(rbind, training_set)
  training_set <- training_set[order(training_set$median, decreasing=T),]
  # diagnostic plot
  training_set_plot <- (ggplot(training_set[1:num_genes,],
                               aes(x=CN34_ct, y=LM1a_ct, col=seq(num_genes))) +
                          geom_point(alpha=0.5) + geom_abline(slope=1, intercept=0) +
                          theme_bw() + ggtitle("reads per transcript") +
                          xlab("CN34") + ylab("LM1a") + labs(col="order") +
                          scale_colour_continuous(low="blue", high="lightgrey")) /
    (ggplot(training_set[1:num_genes,], aes(x=CN34_median, y=LM1a_median, col=seq(num_genes))) +
       geom_point(alpha=0.5) + geom_abline(slope=1, intercept=0) + theme_bw() +
       ggtitle("median reads per codon") + xlab("CN34") + ylab("LM1a") + labs(col="order") +
       scale_colour_continuous(low="blue", high="lightgrey"))
  ggsave(training_set_plot,
         filename=file.path(results_dir, "CN34_LM1a_training_set.pdf"))
  # write results
  training_set <- as.character(training_set$tx_name)[1:num_genes]
  writeLines(training_set, training_set_fname)
  save(mapping, file=file.path(results_dir, "transcript_mapping.Rda"))
} else {
  training_set <- readLines(training_set_fname)
}

# 4. initialize data for regression ---------------------------------------

training_data_fname <- file.path(results_dir, "training_data.Rda")
if(!file.exists(training_data_fname)) {
  # initialize training data
  training_data <- init_data(transcript_fa_fname, transcript_length_fname,
                             d5_d3_subsets=interxn_subsets,
                             f5_length=f5_length, f3_length=f3_length,
                             which_transcripts=training_set)
  training_data$transcript <- relevel(training_data$transcript, ref=training_set[1])
  # count CN34 footprints
  CN34_training <- training_data
  CN34_training$count <- count_footprints(CN34_bam, training_data, "count")
  CN34_training$expt <- "CN34"
  # count LM1a footprints
  LM1a_training <- training_data
  LM1a_training$count <- count_footprints(LM1a_bam, training_data, "count")
  LM1a_training$expt <- "LM1a"
  # combine data for interaction regression
  interxn_training <- rbind(CN34_training, LM1a_training)
  interxn_training$expt <- factor(interxn_training$expt,
                                  levels=c("CN34", "LM1a"))
  save(interxn_training, file=training_data_fname)
} else {
  load(training_data_fname)
}

# 5. compute interaction regression ---------------------------------------

interxn_model_fname <- file.path(results_dir,
                                 paste0("interxn_model_fit_", num_genes, "genes.Rda"))
if(!file.exists(interxn_model_fname)) {
  interxn_model_fit <- MASS::glm.nb(interxn_model,
                                    data=interxn_training, model=F)
  save(interxn_model_fit, file=interxn_model_fname)
} else {
  load(interxn_model_fname)
}

# volcano plot for interaction terms
codon_properties <- read.table(file.path(ref_dir, "yeast_codon_properties.txt"),
                               header=T)
model_coef <- data.frame(coefficients(summary(interxn_model_fit)))
interxn_coef <- model_coef[grep("exptLM1a:", rownames(model_coef)),]
interxn_coef$codon <- substr(rownames(interxn_coef),
                             nchar(rownames(interxn_coef))-2,
                             nchar(rownames(interxn_coef)))
interxn_coef$aa <- codon_properties$aa[match(interxn_coef$codon,
                                             codon_properties$codon)]
interxn_coef$p_log10 <- -log10(interxn_coef$Pr...z..)
volcano_plot <- ggplot(interxn_coef, aes(x=Estimate, y=p_log10)) +
  geom_point(aes(col=aa)) + geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  xlab(expr(beta)) + ylab("-log10(p value)") + labs(color="amino acid") +
  ggtitle("CN34 vs. CN34-LM1a") + theme_classic() +
  geom_text(data=subset(interxn_coef, p.adjust(Pr...z.., "bonferroni") < 0.05),
            aes(x=Estimate, y=p_log10, label=codon),
            check_overlap=T) +
  geom_text(data=subset(interxn_coef, codon %in% c("CGG", "GAA", "GAG")),
            aes(x=Estimate, y=p_log10, label=codon),
            check_overlap=T, col="red") +
  geom_hline(yintercept=-log10(0.05/nrow(interxn_coef)),
             col="red", linetype="dashed")
ggsave(volcano_plot,
       filename=file.path(results_dir, "CN34_LM1a_volcano.pdf"))

rm(interxn_model_fit)

# 6. evaluate bias correction ---------------------------------------------

# CN34
## fit regression
CN34_training <- subset(interxn_training, expt=="CN34")
CN34_training <- subset(CN34_training,
                        with(CN34_training, paste0("d5_", d5, "_d3_", d3)) %in% CN34_subset_names)
CN34_model_fname <- file.path(results_dir, paste0("CN34_model_fit_", num_genes, "genes.Rda"))
if(!file.exists(CN34_model_fname)) {
  CN34_model_fit <- MASS::glm.nb(orig_model, data=CN34_training, model=F)
  save(CN34_model_fit, file=CN34_model_fname)
} else {
  load(CN34_model_fname)
}
## correct counts
CN34_bam$corrected_count <- correct_bias(CN34_bam, CN34_model_fit)
save(CN34_bam, file=file.path(results_dir, "CN34_bam.Rda"))
CN34_training$corrected_count <- count_footprints(CN34_bam, CN34_training, "corrected_count")
save(CN34_training, file=file.path(results_dir, "CN34_training.Rda"))
## evaluate bias
CN34_bias_fname <- file.path(results_dir, "CN34_bias.Rda")
if(!file.exists(CN34_bias_fname)) {
  CN34_bias_raw_codon <- evaluate_bias(CN34_training, which_column="count",
                                      transcript_fa_fname, transcript_length_fname,
                                      type="codon")
  CN34_bias_raw_nt <- evaluate_bias(CN34_training, which_column="count",
                                   transcript_fa_fname, transcript_length_fname,
                                   type="nt")
  CN34_bias_corrected_codon <- evaluate_bias(CN34_training, which_column="corrected_count",
                                            transcript_fa_fname, transcript_length_fname,
                                            type="codon")
  CN34_bias_corrected_nt <- evaluate_bias(CN34_training, which_column="corrected_count",
                                         transcript_fa_fname, transcript_length_fname,
                                         type="nt")
  save(CN34_bias_raw_codon, CN34_bias_raw_nt,
       CN34_bias_corrected_codon, CN34_bias_corrected_nt,
       file=file.path(results_dir, "CN34_bias.Rda"))
  CN34_bias_plot <- (plot_bias(CN34_bias_raw_codon)+ggtitle("CN34: uncorrected") + ylim(0, 0.055)) +
    (plot_bias(CN34_bias_raw_nt) + ylim(0, 0.02)) +
    (plot_bias(CN34_bias_corrected_codon)+ggtitle("CN34: corrected") + ylim(0, 0.055)) +
    (plot_bias(CN34_bias_corrected_nt) + ylim(0, 0.02)) +
    plot_layout(ncol=2, byrow=F)
  ggsave(CN34_bias_plot,
         filename=file.path(results_dir, "CN34_bias_plot.pdf"))
}

# LM1a
## fit regression
LM1a_training <- subset(interxn_training, expt=="LM1a")
LM1a_training <- subset(LM1a_training,
                        with(LM1a_training, paste0("d5_", d5, "_d3_", d3)) %in% LM1a_subset_names)
LM1a_model_fname <- file.path(results_dir, paste0("LM1a_model_fit_", num_genes, "genes.Rda"))
if(!file.exists(LM1a_model_fname)) {
  LM1a_model_fit <- MASS::glm.nb(orig_model, data=LM1a_training, model=F)
  save(LM1a_model_fit, file=LM1a_model_fname)
} else {
  load(LM1a_model_fname)
}
## correct counts
LM1a_bam$corrected_count <- correct_bias(LM1a_bam, LM1a_model_fit)
save(LM1a_bam, file=file.path(results_dir, "LM1a_bam.Rda"))
LM1a_training$corrected_count <- count_footprints(LM1a_bam, LM1a_training, "corrected_count")
save(LM1a_training, file=file.path(results_dir, "LM1a_training.Rda"))
## evaluate bias
LM1a_bias_fname <- file.path(results_dir, "LM1a_bias.Rda")
if(!file.exists(LM1a_bias_fname)) {
  LM1a_bias_raw_codon <- evaluate_bias(LM1a_training, which_column="count",
                                      transcript_fa_fname, transcript_length_fname,
                                      type="codon")
  LM1a_bias_raw_nt <- evaluate_bias(LM1a_training, which_column="count",
                                   transcript_fa_fname, transcript_length_fname,
                                   type="nt")
  LM1a_bias_corrected_codon <- evaluate_bias(LM1a_training, which_column="corrected_count",
                                            transcript_fa_fname, transcript_length_fname,
                                            type="codon")
  LM1a_bias_corrected_nt <- evaluate_bias(LM1a_training, which_column="corrected_count",
                                         transcript_fa_fname, transcript_length_fname,
                                         type="nt")
  save(LM1a_bias_raw_codon, LM1a_bias_raw_nt,
       LM1a_bias_corrected_codon, LM1a_bias_corrected_nt,
       file=file.path(results_dir, "LM1a_bias.Rda"))
  LM1a_bias_plot <- (plot_bias(LM1a_bias_raw_codon)+ggtitle("LM1a: uncorrected") + ylim(0, 0.05)) +
    (plot_bias(LM1a_bias_raw_nt) + ylim(0, 0.03)) +
    (plot_bias(LM1a_bias_corrected_codon)+ggtitle("LM1a: corrected") + ylim(0, 0.05)) +
    (plot_bias(LM1a_bias_corrected_nt) + ylim(0, 0.03)) +
    plot_layout(ncol=2, byrow=F)
  ggsave(LM1a_bias_plot,
         filename=file.path(results_dir, "LM1a_bias_plot.pdf"))
}

q(save="no")