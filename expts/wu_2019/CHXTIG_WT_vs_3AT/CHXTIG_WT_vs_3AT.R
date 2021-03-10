rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

results_dir <- file.path(here(), "expts", "wu_2019", "CHXTIG_WT_vs_3AT")

# individual regressions --------------------------------------------------

load(file.path(here(), "expts", "wu_2019", "CHX_TIG_WT", "CHX_TIG_WT_fit_150.Rda"))
load(file.path(here(), "expts", "wu_2019", "CHX_TIG_3AT", "CHX_TIG_3AT_fit_150.Rda"))

CHX_TIG_WT <- data.frame(summary(CHX_TIG_WT_fit_150)$coefficients)
CHX_TIG_3AT <- data.frame(summary(CHX_TIG_3AT_fit_150)$coefficients)

A_coefs <- CHX_TIG_WT_fit_150$xlevels$A[-1]
A_coefs <- data.frame(A_site = A_coefs,
                      CHX_TIG_WT = CHX_TIG_WT$Estimate[match(paste0("A", A_coefs), rownames(CHX_TIG_WT))],
                      CHX_TIG_WT_se = CHX_TIG_WT$Std..Error[match(paste0("A", A_coefs), rownames(CHX_TIG_WT))],
                      CHX_TIG_3AT = CHX_TIG_3AT$Estimate[match(paste0("A", A_coefs), rownames(CHX_TIG_3AT))],
                      CHX_TIG_3AT_se = CHX_TIG_3AT$Std..Error[match(paste0("A", A_coefs), rownames(CHX_TIG_3AT))])
A_coefs$CHX_TIG_WT_est <- exp(A_coefs$CHX_TIG_WT)
A_coefs$CHX_TIG_WT_upper <- exp(A_coefs$CHX_TIG_WT + A_coefs$CHX_TIG_WT_se)
A_coefs$CHX_TIG_WT_lower <- exp(A_coefs$CHX_TIG_WT - A_coefs$CHX_TIG_WT_se)
A_coefs$CHX_TIG_3AT_est <- exp(A_coefs$CHX_TIG_3AT)
A_coefs$CHX_TIG_3AT_upper <- exp(A_coefs$CHX_TIG_3AT + A_coefs$CHX_TIG_3AT_se)
A_coefs$CHX_TIG_3AT_lower <- exp(A_coefs$CHX_TIG_3AT - A_coefs$CHX_TIG_3AT_se)

codon_table <- read.table(file.path(here(), "reference_data", "yeast_codon_properties.txt"),
                          header=T, stringsAsFactors=F)
A_coefs$aa <- codon_table$aa[match(A_coefs$A_site, codon_table$codon)]

# plot
marginal_plot <- ggplot(A_coefs, aes(x=CHX_TIG_WT_est, y=CHX_TIG_3AT_est)) + geom_point(aes(col=aa)) +
  geom_errorbar(aes(xmin=CHX_TIG_WT_lower, xmax=CHX_TIG_WT_upper,
                    ymin=CHX_TIG_3AT_lower, ymax=CHX_TIG_3AT_upper)) +
  geom_abline(slope=1, intercept=0) + theme_classic() +
  ggtitle("Wu 2019: CHX TIG", subtitle="trained on 19-23mers in top translated 150 genes (separate per dataset)") +
  xlab(expr(paste("WT: exp(", beta, ")"))) + ylab(expr(paste("3AT: exp(", beta, ")"))) +
  geom_text(data=subset(A_coefs, aa=="H"), aes(label=A_site, x=CHX_TIG_WT_est+0.2, y=CHX_TIG_3AT_est+0.2))
ggsave(marginal_plot, filename=file.path(results_dir, "marginal_coefficients.pdf"),
       device="pdf", width=10, height=5, units="in")

# regression with interaction term ----------------------------------------

scripts_dir <- file.path(here(), "scripts")
source(file.path(scripts_dir, "helper.R"))
source(file.path(scripts_dir, "prep_data.R"))
source(file.path(scripts_dir, "correct_bias.R"))
source(file.path(scripts_dir, "evaluate_bias.R"))

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
f5_length <- 2
f3_length <- 3

# load bam alignments
load(file.path(here(), "expts", "wu_2019", "CHX_TIG_WT", "CHX_TIG_WT_bam.Rda"))
load(file.path(here(), "expts", "wu_2019", "CHX_TIG_3AT", "CHX_TIG_3AT_bam.Rda"))

# compute top translated transcripts
transcripts_WT <- aggregate(count ~ transcript, data=CHX_TIG_WT_bam, FUN=sum)
transcripts_WT <- transcripts_WT[order(transcripts_WT$count, decreasing=T),]
transcripts_3AT <- aggregate(count ~ transcript, data=CHX_TIG_3AT_bam, FUN=sum)
transcripts_3AT <- transcripts_3AT[order(transcripts_3AT$count, decreasing=T),]
which_transcripts <- as.character(transcripts_WT$transcript[1:250])
which_transcripts <- which_transcripts[which_transcripts %in% transcripts_3AT$transcript[1:250]]

# decide d5/d3 subsets
load(file.path(here(), "expts", "wu_2019", "CHX_TIG_WT", "CHX_TIG_WT_subsets.Rda"))
subsets_WT <- d5_d3_subsets
load(file.path(here(), "expts", "wu_2019", "CHX_TIG_3AT", "CHX_TIG_3AT_subsets.Rda"))
subsets_3AT <- d5_d3_subsets
which_subsets <- unique(rbind(subsets_WT, subsets_3AT))

# initiate data.frame for regression
training_data <- init_data(transcript_fa_fname, transcript_length_fname,
                           d5_d3_subsets=which_subsets, f5_length=f5_length, f3_length=f3_length,
                           which_transcripts=which_transcripts)
training_WT <- within(training_data, {
  count <- count_footprints(CHX_TIG_WT_bam, training_data, "count")
  expt <- "WT"
})
training_3AT <- within(training_data, {
  count <- count_footprints(CHX_TIG_3AT_bam, training_data, "count")
  expt <- "3AT"
})
training_data <- rbind(training_WT, training_3AT)

# compute regression model
model <- formula(count ~ transcript + A*expt + P + E + d5*genome_f5 + d3*genome_f3)
WT_3AT_fit <- MASS::glm.nb(model, data=training_data, model=F)
save(WT_3AT_fit, file=file.path(results_dir, "WT_3AT_fit.Rda"))

