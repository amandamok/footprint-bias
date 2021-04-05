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

results_dir <- file.path(here(), "expts", "goodarzi_2016", "MDA_vs_MDA-LM2")
raw_data_dir <- file.path(here(), "expts", "goodarzi_2016", "raw_data")

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- 100
exclude_codon5 <- 10
exclude_codon3 <- 10
minimum_cds_length <- exclude_codon5 + exclude_codon3 + 10

sra_MDA <- c("SRR3129936", "SRR3129937")
sra_LM2 <- c("SRR3129940", "SRR3129941")

orig_model <- formula(count ~ transcript + A + P + E + d5*genome_f5 + d3*genome_f3)
interxn_model <- formula(count ~ transcript*expt + A*expt + P + E + d5*genome_f5 + d3*genome_f3)

# 1. load footprint alignments --------------------------------------------

# MDA
MDA_bam_fname <- file.path(results_dir, "MDA.bam")
if(!file.exists(MDA_bam_fname)) {
  system(paste("samtools merge", MDA_bam_fname,
               file.path(raw_data_dir, sra_MDA, paste0(sra_MDA, ".transcript.bam"))))
}
MDA_bam_dat_fname <- file.path(results_dir, "MDA_bam.Rda")
if(!file.exists(MDA_bam_dat_fname)) {
  MDA_bam <- load_bam(MDA_bam_fname, transcript_fa_fname,
                      transcript_length_fname, offsets_fname,
                      f5_length=f5_length, f3_length=f3_length)
  save(MDA_bam, file=MDA_bam_dat_fname)
} else {
  load(MDA_bam_dat_fname)
}

# LM2
LM2_bam_fname <- file.path(results_dir, "LM2.bam")
if(!file.exists(LM2_bam_fname)) {
  system(paste("samtools merge", LM2_bam_fname,
               file.path(raw_data_dir, sra_LM2, paste0(sra_LM2, ".transcript.bam"))))
}
LM2_bam_dat_fname <- file.path(results_dir, "LM2_bam.Rda")
if(!file.exists(LM2_bam_dat_fname)) {
  LM2_bam <- load_bam(LM2_bam_fname, transcript_fa_fname,
                      transcript_length_fname, offsets_fname,
                      f5_length=f5_length, f3_length=f3_length)
  save(LM2_bam, file=LM2_bam_dat_fname)
} else {
  load(LM2_bam_dat_fname)
}

# 2. compute size/frame subsets -------------------------------------------

MDA_subsets_fname <- file.path(results_dir, "MDA_subsets.Rda")
if(!file.exists(MDA_subsets_fname)) {
  MDA_d5_d3 <- count_d5_d3(MDA_bam)
  MDA_d5_d3_subsets <- MDA_d5_d3$counts[1:(which(MDA_d5_d3$counts$proportion>min_prop)[1]),
                                        c("d5", "d3")]
  MDA_subset_names <- sapply(seq(nrow(MDA_d5_d3_subsets)),
                             function(x) {
                               paste("d5", MDA_d5_d3_subsets$d5[x],
                                     "d3", MDA_d5_d3_subsets$d3[x],
                                     sep="_")
                             })
  save(MDA_d5_d3, MDA_d5_d3_subsets, MDA_subset_names,
       file=MDA_subsets_fname)
} else {
  load(MDA_subsets_fname)
}

LM2_subsets_fname <- file.path(results_dir, "LM2_subsets.Rda")
if(!file.exists(LM2_subsets_fname)) {
  LM2_d5_d3 <- count_d5_d3(LM2_bam)
  LM2_d5_d3_subsets <- LM2_d5_d3$counts[1:(which(LM2_d5_d3$counts$proportion>min_prop)[1]),
                                        c("d5", "d3")]
  LM2_subset_names <- sapply(seq(nrow(LM2_d5_d3_subsets)),
                             function(x) {
                               paste("d5", LM2_d5_d3_subsets$d5[x],
                                     "d3", LM2_d5_d3_subsets$d3[x],
                                     sep="_")
                             })
  save(LM2_d5_d3, LM2_d5_d3_subsets, LM2_subset_names,
       file=LM2_subsets_fname)
} else {
  load(LM2_subsets_fname)
}

interxn_subsets <- unique(rbind(MDA_d5_d3_subsets, LM2_d5_d3_subsets))
save(interxn_subsets,
     file=file.path(results_dir, "interxn_subsets.Rda"))

# 3. pick transcripts for training ----------------------------------------

training_set_fname <- file.path(results_dir, "training_transcripts.txt")
if(!file.exists(training_set_fname)) {
  # calculate
  ## transcript count
  MDA_ct <- aggregate(count ~ transcript, data=MDA_bam, FUN=sum)
  mapping$MDA_ct <- MDA_ct$count[match(mapping$tx_name, MDA_ct$transcript)]
  LM2_ct <- aggregate(count ~ transcript, data=LM2_bam, FUN=sum)
  mapping$LM2_ct <- LM2_ct$count[match(mapping$tx_name, LM2_ct$transcript)]
  ## median codon density
  MDA_density <- calculate_transcript_density(MDA_bam, transcript_length_fname,
                                              statistic="median")
  mapping$MDA_median <- MDA_density[match(mapping$tx_name, names(MDA_density))]
  LM2_density <- calculate_transcript_density(LM2_bam, transcript_length_fname,
                                              statistic="median")
  mapping$LM2_median <- LM2_density[match(mapping$tx_name, names(LM2_density))]
  ## average of median codon density
  mapping[is.na(mapping)] <- 0
  mapping$median <- with(mapping, (MDA_median + LM2_median)/2)
  # filter transcripts
  ## transcript counts > 0 for both conditions
  training_set <- subset(mapping, (MDA_ct > 0) & (LM2_ct > 0)) # n=22,405
  ## only 1 stop codon
  num_stop_codons <- transcript_codons[match(training_set$tx_name,
                                             rownames(transcript_codons)),]
  num_stop_codons <- rowSums(num_stop_codons[, c("TAG", "TAA", "TGA")])
  training_set <- subset(training_set, num_stop_codons == 1) # n=20,714
  ## sufficient number of codons
  num_codons <- transcript_lengths$cds_length[match(training_set$tx_name,
                                                    transcript_lengths$transcript)] / 3
  training_set <- subset(training_set, num_codons > minimum_cds_length) # n=20,685
  # pick top representative transcript per gene: maximum average of median codon density
  training_set <- split(training_set, training_set$gene_id) # n=11,112
  training_set <- lapply(training_set,
                         function(x) {
                           x[which.max(x$median),]
                         })
  training_set <- do.call(rbind, training_set)
  training_set <- training_set[order(training_set$median, decreasing=T),]
  # diagnostic plot
  training_set_plot <- (ggplot(training_set[1:num_genes,],
                               aes(x=MDA_ct, y=LM2_ct, col=seq(num_genes))) +
                          geom_point(alpha=0.5) + geom_abline(slope=1, intercept=0) +
                          theme_bw() + ggtitle("reads per transcript") +
                          xlab("MDA") + ylab("LM2") + labs(col="order") +
                          scale_colour_continuous(low="blue", high="lightgrey")) /
    (ggplot(training_set[1:num_genes,], aes(x=MDA_median, y=LM2_median, col=seq(num_genes))) +
       geom_point(alpha=0.5) + geom_abline(slope=1, intercept=0) + theme_bw() +
       ggtitle("median reads per codon") + xlab("MDA") + ylab("LM2") + labs(col="order") +
       scale_colour_continuous(low="blue", high="lightgrey"))
  ggsave(training_set_plot,
         filename=file.path(results_dir, "MDA_LM2_training_set.pdf"))
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
  # count MDA footprints
  MDA_training <- training_data
  MDA_training$count <- count_footprints(MDA_bam, training_data, "count")
  MDA_training$expt <- "MDA"
  # count LM2 footprints
  LM2_training <- training_data
  LM2_training$count <- count_footprints(LM2_bam, training_data, "count")
  LM2_training$expt <- "LM2"
  # combine data for interaction regression
  interxn_training <- rbind(MDA_training, LM2_training)
  interxn_training$expt <- factor(interxn_training$expt,
                                  levels=c("MDA", "LM2"))
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
interxn_coef <- model_coef[grep("exptLM2:", rownames(model_coef)),]
interxn_coef$codon <- substr(rownames(interxn_coef),
                             nchar(rownames(interxn_coef))-2,
                             nchar(rownames(interxn_coef)))
interxn_coef$aa <- codon_properties$aa[match(interxn_coef$codon,
                                             codon_properties$codon)]
interxn_coef$p_log10 <- -log10(interxn_coef$Pr...z..)
volcano_plot <- ggplot(interxn_coef, aes(x=Estimate, y=p_log10)) +
  geom_point(aes(col=aa)) + geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  xlab(expr(beta)) + ylab("-log10(p value)") + labs(color="amino acid") +
  ggtitle("MDA vs. MDA-LM2") + theme_classic() +
  geom_text(data=subset(interxn_coef, p.adjust(Pr...z.., "bonferroni") < 0.05),
            aes(x=Estimate, y=p_log10, label=codon),
            check_overlap=T) +
  geom_text(data=subset(interxn_coef, codon %in% c("CGG", "GAA", "GAG")),
            aes(x=Estimate, y=p_log10, label=codon),
            check_overlap=T, col="red") +
  geom_hline(yintercept=-log10(0.05/nrow(interxn_coef)),
             col="red", linetype="dashed")
ggsave(volcano_plot,
       filename=file.path(results_dir, "MDA_LM2_volcano.pdf"))

# 6. evaluate bias correction ---------------------------------------------

# MDA
## fit regression
MDA_training <- subset(interxn_training, expt=="MDA")
MDA_training <- subset(MDA_training,
                       with(MDA_training, paste0("d5_", d5, "_d3_", d3)) %in% MDA_subset_names)
MDA_model_fname <- file.path(results_dir, paste0("MDA_model_fit_", num_genes, "genes.Rda"))
if(!file.exists(MDA_model_fname)) {
  MDA_model_fit <- MASS::glm.nb(orig_model, data=MDA_training, model=F)
  save(MDA_model_fit, file=MDA_model_fname)
} else {
  load(MDA_model_fname)
}
## correct counts
MDA_bam$corrected_count <- correct_bias(MDA_bam, MDA_model_fit)
save(MDA_bam, file=file.path(results_dir, "MDA_bam.Rda"))
MDA_training$corrected_count <- count_footprints(MDA_bam, MDA_training, "corrected_count")
save(MDA_training, file=file.path(results_dir, "MDA_training.Rda"))
## evaluate bias
MDA_bias_fname <- file.path(results_dir, "MDA_bias.Rda")
if(!file.exists(MDA_bias_fname)) {
  MDA_bias_raw_codon <- evaluate_bias(MDA_training, which_column="count",
                                      transcript_fa_fname, transcript_length_fname,
                                      type="codon")
  MDA_bias_raw_nt <- evaluate_bias(MDA_training, which_column="count",
                                   transcript_fa_fname, transcript_length_fname,
                                   type="nt")
  MDA_bias_corrected_codon <- evaluate_bias(MDA_training, which_column="corrected_count",
                                            transcript_fa_fname, transcript_length_fname,
                                            type="codon")
  MDA_bias_corrected_nt <- evaluate_bias(MDA_training, which_column="corrected_count",
                                         transcript_fa_fname, transcript_length_fname,
                                         type="nt")
  save(MDA_bias_raw_codon, MDA_bias_raw_nt,
       MDA_bias_corrected_codon, MDA_bias_corrected_nt,
       file=file.path(results_dir, "MDA_bias.Rda"))
  MDA_bias_plot <- (plot_bias(MDA_bias_raw_codon)+ggtitle("MDA: uncorrected") + ylim(0, 0.05)) +
    (plot_bias(MDA_bias_raw_nt) + ylim(0, 0.02)) +
    (plot_bias(MDA_bias_corrected_codon)+ggtitle("MDA: corrected") + ylim(0, 0.05)) +
    (plot_bias(MDA_bias_corrected_nt) + ylim(0, 0.02)) +
    plot_layout(ncol=2, byrow=F)
  ggsave(MDA_bias_plot,
         filename=file.path(results_dir, "MDA_bias_plot.pdf"))
}

# LM2
## fit regression
LM2_training <- subset(interxn_training, expt=="LM2")
LM2_training <- subset(LM2_training,
                       with(LM2_training, paste0("d5_", d5, "_d3_", d3)) %in% LM2_subset_names)
LM2_model_fname <- file.path(results_dir, paste0("LM2_model_fit_", num_genes, "genes.Rda"))
if(!file.exists(LM2_model_fname)) {
  LM2_model_fit <- MASS::glm.nb(orig_model, data=LM2_training, model=F)
  save(LM2_model_fit, file=LM2_model_fname)
} else {
  load(LM2_model_fname)
}
## correct counts
LM2_bam$corrected_count <- correct_bias(LM2_bam, LM2_model_fit)
save(LM2_bam, file=file.path(results_dir, "LM2_bam.Rda"))
LM2_training$corrected_count <- count_footprints(LM2_bam, LM2_training, "corrected_count")
save(LM2_training, file=file.path(results_dir, "LM2_training.Rda"))
## evaluate bias
LM2_bias_fname <- file.path(results_dir, "LM2_bias.Rda")
if(!file.exists(LM2_bias_fname)) {
  LM2_bias_raw_codon <- evaluate_bias(LM2_training, which_column="count",
                                      transcript_fa_fname, transcript_length_fname,
                                      type="codon")
  LM2_bias_raw_nt <- evaluate_bias(LM2_training, which_column="count",
                                   transcript_fa_fname, transcript_length_fname,
                                   type="nt")
  LM2_bias_corrected_codon <- evaluate_bias(LM2_training, which_column="corrected_count",
                                            transcript_fa_fname, transcript_length_fname,
                                            type="codon")
  LM2_bias_corrected_nt <- evaluate_bias(LM2_training, which_column="corrected_count",
                                         transcript_fa_fname, transcript_length_fname,
                                         type="nt")
  save(LM2_bias_raw_codon, LM2_bias_raw_nt,
       LM2_bias_corrected_codon, LM2_bias_corrected_nt,
       file=file.path(results_dir, "LM2_bias.Rda"))
  LM2_bias_plot <- (plot_bias(LM2_bias_raw_codon)+ggtitle("LM2: uncorrected") + ylim(0, 0.045)) +
    (plot_bias(LM2_bias_raw_nt) + ylim(0, 0.03)) +
    (plot_bias(LM2_bias_corrected_codon)+ggtitle("LM2: corrected") + ylim(0, 0.045)) +
    (plot_bias(LM2_bias_corrected_nt) + ylim(0, 0.03)) +
    plot_layout(ncol=2, byrow=F)
  ggsave(LM2_bias_plot,
         filename=file.path(results_dir, "LM2_bias_plot.pdf"))
}

q(save="no")