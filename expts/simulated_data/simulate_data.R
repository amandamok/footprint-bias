rm(list=ls())

library(parallel)
library(here)

# load scripts ------------------------------------------------------------

simriboseq_dir <- "~/simRiboSeq"

script_dir <- file.path(simriboseq_dir, "scripts")
output_dir <- file.path(here(), "expts", "scer_simASite_biasDelta_f5_2_f3_3")
ref_dir <- file.path(simriboseq_dir, "refData")

scripts <- c("helper.R", "simTranscriptome.R", "simRibosomeDist.R", "simFootprints.R")
for(x in scripts) {source(file.path(script_dir, x))}

# reference data ----------------------------------------------------------

## general params:
codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")
delta5_bias <- c(0.4, 0.3, 0.2)
names(delta5_bias) <- as.character(15:17)
delta3_bias <- c(0.35, 0.4, 0.25)
names(delta3_bias) <- as.character(9:11)
minSize <- 27
maxSize <- 31
nReads <- 6e7
partSize <- 1e6
nParts <- nReads/partSize

## weinberg expt: gene lengths and abundances for simRho()
weinberg_file <- "cts_by_codon.size.27.31.txt"
weinberg_data <- readRawProfiles(file.path(ref_dir, weinberg_file))
weinberg_lengths <- lengths(weinberg_data)
weinberg_abundances <- sapply(weinberg_data, sum)
weinberg_nRibosomes <- sum(weinberg_abundances)
rm(weinberg_data)

## weinberg expt: codon TE scores for simPi()
weinberg_codonTEfile <- "tunney_supp_table_2_codon_scores.csv"
weinberg_codonScores <- read.csv(file.path(ref_dir, weinberg_codonTEfile), header=T, stringsAsFactors=F)
weinberg_codonTE <- weinberg_codonScores$X0
weinberg_codonTE <- weinberg_codonTE - min(weinberg_codonTE) # scale up so min = 0.1
weinberg_stopCodons <- codons[which(!(codons %in% weinberg_codonScores$codon))]
weinberg_codonTE <- c(weinberg_codonTE, rep(0, length(weinberg_stopCodons))) # add TE for stop codons
names(weinberg_codonTE) <- c(weinberg_codonScores$codon, weinberg_stopCodons)

## green expt: bias scores for ligate() & circularize()
# /mnt/lareaulab/rtunney/iXnos/results/green/s28_cod_n5p4_nt_n15p14/epoch30/codon_scores.tsv
bias_seq_2nt <- unique(substr(codons, start=1, stop=2))
green_biasFile <- "codon_scores.tsv"
green_biasScores <- read.table(file.path(ref_dir, green_biasFile))
colnames(green_biasScores) <- as.character(seq.int(from=-5, length.out=ncol(green_biasScores)))
# 5' bias probabilities
green_p5bias <- exp(green_biasScores[,"-5"])
green_p5bias <- green_p5bias/max(green_p5bias, na.rm=T)
green_p5bias[is.na(green_p5bias)] <- 0
names(green_p5bias) <- sort(codons)
green_p5bias_2nt <- sapply(bias_seq_2nt,
                           function(x) {
                             mean(green_p5bias[substr(names(green_p5bias), start=1, stop=2)==x])
                           })
green_p5bias_2nt[green_p5bias_2nt==0] <- 0.01
# 3' bias probabilities
green_n3bias <- exp(green_biasScores[,"3"])
green_n3bias <- green_n3bias/max(green_n3bias, na.rm=T)
green_n3bias[is.na(green_n3bias)] <- 0
names(green_n3bias) <- sort(codons)
green_n3bias_3nt <- green_n3bias
green_n3bias_3nt[green_n3bias_3nt==0] <- 0.01
# no bias (uniform bias)
no_p5bias <- rep(1, length(green_p5bias_2nt))
names(no_p5bias) <- names(green_p5bias_2nt)
no_n3bias <- rep(1, length(codons))
names(no_n3bias) <- codons

## non-templated base
noRTbias <- c(1, 0, 0, 0, 0)
names(noRTbias) <- c("", "A", "T", "C", "G")


# generate ribosome profiles ----------------------------------------------

# generate transcriptome: s. cerevisiae
scerFile <- "scer.transcripts.20cds20.fa"
scer <- readFAfile(file.path(ref_dir, scerFile), pad5=20, pad3=20)

# simulate ribosome distributions
## weinberg data for expt lengths, abundances, codon TE
scer_rho <- simRho(scer, exptLengths=weinberg_lengths, exptAbundances=weinberg_abundances)
scer_pi <- simPi(scer, codonTE=weinberg_codonTE)
save(scer_rho, scer_pi, file=file.path(output_dir, "rho_pi.Rda"))

# simulate data: no bias --------------------------------------------------

# biased delta5, delta3
# minSize=27, maxSize=31
# 5' bias: no_p5bias
# 3' bias: no_n3bias
# no extra base: noRTbias

simulation_name <- "noBias"
if(!dir.exists(file.path(output_dir, simulation_name))) {
  dir.create(file.path(output_dir, simulation_name))
}

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0(simulation_name, "_part", i)
  part_filename <- paste0(simulation_name, "_60Mreads_part", i)
  assign(partName,
         value=simFootprints(scer, nRibosomes=partSize,
                             rhos=scer_rho, pis=scer_pi,
                             delta5=delta5_bias, delta3=delta3_bias,
                             ligBias=no_n3bias, RTBias=noRTbias, circBias=no_p5bias,
                             digest_transcript=digest_transcript))
  writeFootprintsFQ(get(partName),
                    file.path(output_dir, simulation_name, paste0(part_filename, ".fq")))
  save(list=partName,
       file=file.path(output_dir, simulation_name, paste0(part_filename, ".Rda")))
  rm(list=partName)
}


# simulate data: 3' bias --------------------------------------------------

# biased delta5, delta3
# minSize=27, maxSize=31
# 5' bias: no_p5bias
# 3' bias: green_n3bias_3nt
# no extra base: noRTbias

simulation_name <- "n3Bias"
if(!dir.exists(file.path(output_dir, simulation_name))) {
  dir.create(file.path(output_dir, simulation_name))
}

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0(simulation_name, "_part", i)
  part_filename <- paste0(simulation_name, "_60Mreads_part", i)
  assign(partName,
         value=simFootprints(scer, nRibosomes=partSize,
                             rhos=scer_rho, pis=scer_pi,
                             delta5=delta5_bias, delta3=delta3_bias,
                             ligBias=green_n3bias_3nt, RTBias=noRTbias, circBias=no_p5bias,
                             digest_transcript=digest_transcript))
  writeFootprintsFQ(get(partName),
                    file.path(output_dir, simulation_name, paste0(part_filename, ".fq")))
  save(list=partName,
       file=file.path(output_dir, simulation_name, paste0(part_filename, ".Rda")))
  rm(list=partName)
}

# simulate data: 5' bias --------------------------------------------------

# biased delta5, delta3
# minSize=27, maxSize=31
# 5' bias: green_p5bias_2nt
# 3' bias: no_n3bias
# no extra base: noRTbias

simulation_name <- "p5Bias"
if(!dir.exists(file.path(output_dir, simulation_name))) {
  dir.create(file.path(output_dir, simulation_name))
}

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0(simulation_name, "_part", i)
  part_filename <- paste0(simulation_name, "_60Mreads_part", i)
  assign(partName,
         value=simFootprints(scer, nRibosomes=partSize,
                             rhos=scer_rho, pis=scer_pi,
                             delta5=delta5_bias, delta3=delta3_bias,
                             ligBias=no_n3bias, RTBias=noRTbias, circBias=green_p5bias_2nt,
                             digest_transcript=digest_transcript))
  writeFootprintsFQ(get(partName),
                    file.path(output_dir, simulation_name, paste0(part_filename, ".fq")))
  save(list=partName,
       file=file.path(output_dir, simulation_name, paste0(part_filename, ".Rda")))
  rm(list=partName)
}

# simulate data: 5' and 3' bias -------------------------------------------

# biased delta5, delta3
# minSize=27, maxSize=31
# 5' bias: green_p5bias_2nt
# 3' bias: green_n3bias_3nt
# no extra base: noRTbias

simulation_name <- "bothBias"
if(!dir.exists(file.path(output_dir, simulation_name))) {
  dir.create(file.path(output_dir, simulation_name))
}

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0(simulation_name, "_part", i)
  part_filename <- paste0(simulation_name, "_60Mreads_part", i)
  assign(partName,
         value=simFootprints(scer, nRibosomes=partSize,
                             rhos=scer_rho, pis=scer_pi,
                             delta5=delta5_bias, delta3=delta3_bias,
                             ligBias=green_n3bias_3nt, RTBias=noRTbias, circBias=green_p5bias_2nt,
                             digest_transcript=digest_transcript))
  writeFootprintsFQ(get(partName),
                    file.path(output_dir, simulation_name, paste0(part_filename, ".fq")))
  save(list=partName,
       file=file.path(output_dir, simulation_name, paste0(part_filename, ".Rda")))
  rm(list=partName)
}

# quit --------------------------------------------------------------------

quit(save="no")