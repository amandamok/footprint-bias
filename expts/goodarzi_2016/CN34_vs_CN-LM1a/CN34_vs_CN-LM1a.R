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

results_dir <- file.path(here(), "expts", "goodarzi_2016", "CN34_vs_CN34-LM1a")
raw_data_dir <- file.path(here(), "expts", "goodarzi_2016", "raw_data")

min_prop <- 0.9
f5_length <- 2
f3_length <- 3
num_genes <- 100
exclude_codon5 <- 10
exclude_codon3 <- 10

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
  training_set <- subset(mapping, (CN34_ct > 0) & (LM1a_ct > 0)) # n=22,405
  ## only 1 stop codon
  num_stop_codons <- transcript_codons[match(training_set$tx_name,
                                             rownames(transcript_codons)),]
  num_stop_codons <- rowSums(num_stop_codons[, c("TAG", "TAA", "TGA")])
  training_set <- subset(training_set, num_stop_codons == 1) # n=20,714
  # pick top representative transcript per gene: maximum average of median codon density
  training_set <- split(training_set, training_set$gene_id) # n=11,113
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
  writeLines(as.character(training_set$tx_name), training_set_fname)
  save(mapping, file=file.path(results_dir, "transcript_mapping.tsv"))
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

interxn_model_fname <- file.path(results_dir, "interxn_model_fit_200genes.Rda")
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
  geom_text(data=subset(interxn_coef, p_log10>5),
            aes(x=Estimate, y=p_log10, label=codon),
            nudge_y=0.25, check_overlap=T) +
  geom_text(data=subset(interxn_coef, codon %in% c("CGG", "GAA", "GAG")),
            aes(x=Estimate, y=p_log10, label=codon),
            nudge_y=0.25, check_overlap=T, col="red")
ggsave(volcano_plot,
       filename=file.path(results_dir, "CN34_LM1a_volcano.pdf"))

# 6. evaluate bias correction ---------------------------------------------

# CN34
## fit regression
CN34_model_fname <- file.path(results_dir, "CN34_model_fit_200genes.Rda")
CN34_training <- subset(interxn_training, expt=="CN34")
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
       # CN34_bias_corrected_codon, CN34_bias_corrected_nt,
       file=file.path(results_dir, "CN34_bias.Rda"))
}

# # debugging
# CN34_transcript_count <- aggregate(count ~ transcript,
#                                   data=CN34_bam,
#                                   FUN=sum)
# CN34_transcript_count$corrected_count <- aggregate(corrected_count ~ transcript,
#                                                   data=CN34_bam,
#                                                   FUN=sum)$corrected_count
# CN34_transcript_count$gene_name <- mapping$gene_name[match(CN34_transcript_count$transcript,
#                                                           mapping$tx_name)]
# CN34_training_transcripts <- levels(CN34_training$transcript)
# CN34_training_genes <- mapping$gene_name[match(CN34_training_transcripts, mapping$tx_name)]
# CN34_transcript_count_subset <- subset(CN34_transcript_count, gene_name %in% CN34_training_genes)
# (ggplot(CN34_transcript_count_subset, aes(x=count, y=corrected_count, col=gene_name)) +
#     geom_point(alpha=0.5, size=2) + geom_abline(slope=1, intercept=0) + theme_bw() +
#     theme(legend.position="none") + ggtitle("CN34")) +
#   ((ggplot(CN34_transcript_count_subset, aes(x=count)) +
#       geom_density() + theme_bw() + ylab("")) /
#      (ggplot(CN34_transcript_count_subset, aes(x=corrected_count)) +
#         geom_density() + theme_bw() + ylab("")))
# CN34_coef <- data.frame(summary(CN34_model_fit)$coefficients)
# CN34_coef$coef <- NA
# CN34_coef$coef[grepl("^transcript", rownames(CN34_coef))] <- "transcript"
# CN34_coef$coef[grepl("^A", rownames(CN34_coef))] <- "A"
# CN34_coef$coef[grepl("^P", rownames(CN34_coef))] <- "P"
# CN34_coef$coef[grepl("^E", rownames(CN34_coef))] <- "E"
# CN34_coef$coef[grepl("^genome_f5", rownames(CN34_coef))] <- "f5"
# CN34_coef$coef[grepl("^genome_f3", rownames(CN34_coef))] <- "f3"
# CN34_coef$coef[grepl(":genome_f5", rownames(CN34_coef))] <- "f5_interxn"
# CN34_coef$coef[grepl(":genome_f3", rownames(CN34_coef))] <- "f3_interxn"
# CN34_coef$coef[grepl("^d5", rownames(CN34_coef)) & !grepl(":", rownames(CN34_coef))] <- "d5"
# CN34_coef$coef[grepl("^d3", rownames(CN34_coef)) & !grepl(":", rownames(CN34_coef))] <- "d3"
# CN34_coef$p_log10 <- -log10(CN34_coef$Pr...z..)
# ggplot(subset(CN34_coef, coef %in% c("f5", "f3")),
#        aes(x=Estimate, y=p_log10, col=coef)) +
#   geom_point() + theme_bw()

# LM1a
## fit regression
LM1a_model_fname <- file.path(results_dir, "LM1a_model_fit_200genes.Rda")
LM1a_training <- subset(interxn_training, expt=="LM1a")
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
       # LM1a_bias_corrected_codon, LM1a_bias_corrected_nt,
       file=file.path(results_dir, "LM1a_bias.Rda"))
}

q(save="no")