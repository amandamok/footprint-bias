rm(list=ls())

# set up ------------------------------------------------------------------

outputDir <- "~/footprint-bias/scalability"
refDir <- "~/simRiboSeq/refData"

source("~/simRiboSeq/scripts/helper.R")
source("~/footprint-bias/scripts/prep_data.R")

if(!("scalability" %in% list.files("~/footprint-bias"))) {
  system("mkdir ~/footprint-bias/scalability")
}

weinberg_sam_fname <- "~/iXnos/expts/weinberg_20cds20/process/weinberg.footprints.sam"
offsets_fname <- "~/footprint-bias/sample_data/Asite_rules.txt"

regression <- formula(count ~ transcript + A + P + E + d5 + d3 + f5 + f3)

# load data ---------------------------------------------------------------

# transcriptome
scer_transcriptome <- readFAfile(file.path(refDir, "scer.transcripts.20cds20.fa"),
                                 pad5=20, pad3=20)
scer_lengths <- read.table(file.path(refDir, "scer.transcripts.20cds20.lengths.txt"),
                           col.names=c("transcript", "utr5_length", "cds_length", "utr3_length"))

# weinberg cts_by_codon data
weinberg_cts_by_codon <- readRawProfiles(file.path(refDir, "cts_by_codon.size.27.31.txt"))
weinberg_cts_by_transcript <- sapply(weinberg_cts_by_codon, sum)
transcripts <- sort(weinberg_cts_by_transcript, decreasing=T)

# 3 genes -----------------------------------------------------------------

n <- 3
expt <- paste0("genes_", n)
if(!(expt %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, expt)))
}

# select top 3 translated transcripts
transcripts_3 <- transcripts[1:n]
transcripts_seq_3 <- scer_transcriptome[match(names(transcripts_3), names(scer_transcriptome))]

# write transcriptome .fa file
transcripts_nCodons_3 <- scer_lengths$cds_length[match(names(transcripts_seq_3), scer_lengths$transcript)]/3
transcripts_seq_3 <- lapply(1:n,
                            function(x) {
                              # 6x 5' UTR codons + 5x 3' UTR codons
                              transcripts_seq_3[[x]][1:(transcripts_nCodons_3[x]+6+5)] 
                            })
names(transcripts_seq_3) <- names(transcripts_3)
transcript_fa_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15.fa"))
writeTranscriptomeFA(transcripts_seq_3, transcript_fa_fname)

# write transcriptome lengths file
transcripts_lengths_3 <- scer_lengths[match(names(transcripts_3), scer_lengths$transcript),]
transcripts_lengths_3$utr5_length <- 18
transcripts_lengths_3$utr3_length <- 15
transcript_length_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15_lengths.txt"))
write.table(transcripts_lengths_3, transcript_length_fname, quote=F)

# create footprints FASTA file
transcripts_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_names.txt"))
sam_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_footprints.sam"))
writeLines(names(transcripts_3), con=transcripts_fname)
system(paste("grep -wFf", transcripts_fname, weinberg_sam_fname, ">", sam_fname))

# create regression_data
print("3 genes: init data")
time_3_init <- bench::mark(regression_data_3 <- init_data(transcript_fa_fname, 
                                                          transcript_length_fname, bias_length=2),
                           iterations=1)
regression_data_3$transcript <- relevel(regression_data_3$transcript, ref=names(transcripts)[1])
print(paste("Time: ", time_3_init$min, " ; Memory: ", time_3_init$mem_alloc))
print(paste("Number of rows:", nrow(regression_data_3)))

# count footprints
print("3 genes: count footprints")
time_3_count <- bench::mark(regression_data_3 <- count_footprints(sam_fname, transcript_length_fname, 
                                                                  offsets_fname, regression_data_3),
                            iterations=1)
print(paste("Time: ", time_3_count$min, " ; Memory: ", time_3_init$mem_alloc))

# perform regression
print("3 genes: regression")
time_3_regression <- bench::mark(nb_fit_3 <- MASS::glm.nb(regression, data=regression_data_3),
                                 iterations=1)
print(paste("Time: ", time_3_regression$min, " ; Memory: ", time_3_regression$mem_alloc))

# save data
save(regression_data_3, file=file.path(outputDir, expt, "regression_data_3.Rda"))
save(nb_fit_3, file=file.path(outputDir, expt, "negbin_fit_3.Rda"))
save(regression_data_3, nb_fit_3, 
     time_3_init, time_3_count, time_3_regression,
     file=file.path(outputDir, expt, "data_3genes.Rda"))
rm(regression_data_3, nb_fit_3)

# 10 genes ----------------------------------------------------------------

n <- 10
expt <- paste0("genes_", n)
if(!(expt %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, expt)))
}

# select top 10 translated transcripts
transcripts_10 <- transcripts[1:n]
transcripts_seq_10 <- scer_transcriptome[match(names(transcripts_10), names(scer_transcriptome))]

# write transcriptome .fa file
transcripts_nCodons_10 <- scer_lengths$cds_length[match(names(transcripts_seq_10), scer_lengths$transcript)]/3
transcripts_seq_10 <- lapply(1:n,
                             function(x) {
                               # 6x 5' UTR codons + 5x 10' UTR codons
                               transcripts_seq_10[[x]][1:(transcripts_nCodons_10[x]+6+5)] 
                             })
names(transcripts_seq_10) <- names(transcripts_10)
transcript_fa_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15.fa"))
writeTranscriptomeFA(transcripts_seq_10, transcript_fa_fname)

# write transcriptome lengths file
transcripts_lengths_10 <- scer_lengths[match(names(transcripts_10), scer_lengths$transcript),]
transcripts_lengths_10$utr5_length <- 18
transcripts_lengths_10$utr3_length <- 15
transcript_length_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15_lengths.txt"))
write.table(transcripts_lengths_10, transcript_length_fname, quote=F)

# create footprints FASTA file
transcripts_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_names.txt"))
sam_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_footprints.sam"))
writeLines(names(transcripts_10), con=transcripts_fname)
system(paste("grep -wFf", transcripts_fname, weinberg_sam_fname, ">", sam_fname))

# create regression_data
print("10 genes: init data")
time_10_init <- bench::mark(regression_data_10 <- init_data(transcript_fa_fname, 
                                                            transcript_length_fname, bias_length=2),
                            iterations=1)
print(paste("Time: ", time_10_init$min, " ; Memory: ", time_10_init$mem_alloc))
print(paste("Number of rows:", nrow(regression_data_10)))

# count footprints
print("10 genes: count footprints")
time_10_count <- bench::mark(regression_data_10 <- count_footprints(sam_fname, transcript_length_fname, 
                                                                    offsets_fname, regression_data_10),
                             iterations=1)
print(paste("Time: ", time_10_count$min, " ; Memory: ", time_10_count$mem_alloc))

# perform regression
print("10 genes: regression")
time_10_regression <- bench::mark(nb_fit_10 <- MASS::glm.nb(regression, data=regression_data_10),
                                  iterations=1)
print(paste("Time: ", time_10_regression$min, " ; Memory: ", time_10_regression$mem_alloc))

# save data
save(regression_data_10, file=file.path(outputDir, expt, "regression_data_10.Rda"))
save(nb_fit_10, file=file.path(outputDir, expt, "negbin_fit_10.Rda"))
save(regression_data_10, nb_fit_10, 
     time_10_init, time_10_count, time_10_regression,
     file=file.path(outputDir, expt, "data_10genes.Rda"))
rm(regression_data_10, nb_fit_10)

# 20 genes ----------------------------------------------------------------

n <- 20
expt <- paste0("genes_", n)
if(!(expt %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, expt)))
}

# select top 20 translated transcripts
transcripts_20 <- transcripts[1:n]
transcripts_seq_20 <- scer_transcriptome[match(names(transcripts_20), names(scer_transcriptome))]

# write transcriptome .fa file
transcripts_nCodons_20 <- scer_lengths$cds_length[match(names(transcripts_seq_20), scer_lengths$transcript)]/3
transcripts_seq_20 <- lapply(1:n,
                             function(x) {
                               # 6x 5' UTR codons + 5x 20' UTR codons
                               transcripts_seq_20[[x]][1:(transcripts_nCodons_20[x]+6+5)] 
                             })
names(transcripts_seq_20) <- names(transcripts_20)
transcript_fa_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15.fa"))
writeTranscriptomeFA(transcripts_seq_20, transcript_fa_fname)

# write transcriptome lengths file
transcripts_lengths_20 <- scer_lengths[match(names(transcripts_20), scer_lengths$transcript),]
transcripts_lengths_20$utr5_length <- 18
transcripts_lengths_20$utr3_length <- 15
transcript_length_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15_lengths.txt"))
write.table(transcripts_lengths_20, transcript_length_fname, quote=F)

# create footprints FASTA file
transcripts_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_names.txt"))
sam_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_footprints.sam"))
writeLines(names(transcripts_20), con=transcripts_fname)
system(paste("grep -wFf", transcripts_fname, weinberg_sam_fname, ">", sam_fname))

# create regression_data
print("20 genes: init data")
time_20_init <- bench::mark(regression_data_20 <- init_data(transcript_fa_fname, 
                                                            transcript_length_fname, bias_length=2),
                            iterations=1)
print(paste("Time: ", time_20_init$min, " ; Memory: ", time_20_init$mem_alloc))

# count footprints
print("20 genes: count footprints")
time_20_count <- bench::mark(regression_data_20 <- count_footprints(sam_fname, transcript_length_fname, 
                                                                    offsets_fname, regression_data_20),
                             iterations=1)
print(paste("Time: ", time_20_count$min, " ; Memory: ", time_20_count$mem_alloc))

# perform regression
print("20 genes: regression")
time_20_regression <- bench::mark(nb_fit_20 <- MASS::glm.nb(regression, data=regression_data_20),
                                  iterations=1)
print(paste("Time: ", time_20_regression$min, " ; Memory: ", time_20_regression$mem_alloc))

# save data
save(regression_data_20, file=file.path(outputDir, expt, "regression_data_20.Rda"))
save(nb_fit_20, file=file.path(outputDir, expt, "negbin_fit_20.Rda"))
save(regression_data_20, nb_fit_20, 
     time_20_init, time_20_count, time_20_regression,
     file=file.path(outputDir, expt, "data_20genes.Rda"))
rm(regression_data_20, nb_fit_20)

# 50 genes ----------------------------------------------------------------

n <- 50
expt <- paste0("genes_", n)
if(!(expt %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, expt)))
}

# select top 50 translated transcripts
transcripts_50 <- transcripts[1:n]
transcripts_seq_50 <- scer_transcriptome[match(names(transcripts_50), names(scer_transcriptome))]

# write transcriptome .fa file
transcripts_nCodons_50 <- scer_lengths$cds_length[match(names(transcripts_seq_50), scer_lengths$transcript)]/3
transcripts_seq_50 <- lapply(1:n,
                             function(x) {
                               # 6x 5' UTR codons + 5x 50' UTR codons
                               transcripts_seq_50[[x]][1:(transcripts_nCodons_50[x]+6+5)]
                             })
names(transcripts_seq_50) <- names(transcripts_50)
transcript_fa_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15.fa"))
writeTranscriptomeFA(transcripts_seq_50, transcript_fa_fname)

# write transcriptome lengths file
transcripts_lengths_50 <- scer_lengths[match(names(transcripts_50), scer_lengths$transcript),]
transcripts_lengths_50$utr5_length <- 18
transcripts_lengths_50$utr3_length <- 15
transcript_length_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15_lengths.txt"))
write.table(transcripts_lengths_50, transcript_length_fname, quote=F)

# create footprints FASTA file
transcripts_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_names.txt"))
sam_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_footprints.sam"))
writeLines(names(transcripts_50), con=transcripts_fname)
system(paste("grep -wFf", transcripts_fname, weinberg_sam_fname, ">", sam_fname))

# create regression_data
print("50 genes: init data")
time_50_init <- bench::mark(regression_data_50 <- init_data(transcript_fa_fname, 
                                                            transcript_length_fname, bias_length=2),
                            iterations=1)
print(paste("Time: ", time_50_init$min, " ; Memory: ", time_50_init$mem_alloc))
print(paste("Number of rows:", nrow(regression_data_50)))

# count footprints
print("50 genes: count footprints")
time_50_count <- bench::mark(regression_data_50 <- count_footprints(sam_fname, transcript_length_fname, 
                                                                    offsets_fname, regression_data_50),
                             iterations=1)
print(paste("Time: ", time_50_count$min, " ; Memory: ", time_50_count$mem_alloc))

# perform regression
print("50 genes: regression")
time_50_regression <- bench::mark(nb_fit_50 <- MASS::glm.nb(regression, data=regression_data_50),
                                  iterations=1)
print(paste("Time: ", time_50_regression$min, " ; Memory: ", time_50_regression$mem_alloc))

# save data
save(regression_data_50, file=file.path(outputDir, expt, "regression_data_50.Rda"))
save(nb_fit_50, file=file.path(outputDir, expt, "negbin_fit_50.Rda"))
save(regression_data_50, nb_fit_50, 
     time_50_init, time_50_count, time_50_regression,
     file=file.path(outputDir, expt, "data_50genes.Rda"))
rm(regression_data_50, nb_fit_50)

# 100 genes ---------------------------------------------------------------

n <- 100
expt <- paste0("genes_", n)
if(!(expt %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, expt)))
}

# select top 100 translated transcripts
transcripts_100 <- transcripts[1:n]
transcripts_seq_100 <- scer_transcriptome[match(names(transcripts_100), names(scer_transcriptome))]

# write transcriptome .fa file
transcripts_nCodons_100 <- scer_lengths$cds_length[match(names(transcripts_seq_100), scer_lengths$transcript)]/3
transcripts_seq_100 <- lapply(1:n,
                              function(x) {
                                # 6x 5' UTR codons + 5x 100' UTR codons
                                transcripts_seq_100[[x]][1:(transcripts_nCodons_100[x]+6+5)]
                              })
names(transcripts_seq_100) <- names(transcripts_100)
transcript_fa_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15.fa"))
writeTranscriptomeFA(transcripts_seq_100, transcript_fa_fname)

# write transcriptome lengths file
transcripts_lengths_100 <- scer_lengths[match(names(transcripts_100), scer_lengths$transcript),]
transcripts_lengths_100$utr5_length <- 18
transcripts_lengths_100$utr3_length <- 15
transcript_length_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15_lengths.txt"))
write.table(transcripts_lengths_100, transcript_length_fname, quote=F)

# create footprints FASTA file
transcripts_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_names.txt"))
sam_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_footprints.sam"))
writeLines(names(transcripts_100), con=transcripts_fname)
system(paste("grep -wFf", transcripts_fname, weinberg_sam_fname, ">", sam_fname))

# create regression_data
print("100 genes: init data")
time_100_init <- bench::mark(regression_data_100 <- init_data(transcript_fa_fname, 
                                                              transcript_length_fname, bias_length=2),
                             iterations=1)
print(paste("Time: ", time_100_init$min, " ; Memory: ", time_100_init$mem_alloc))

# count footprints
print("100 genes: count footprints")
time_100_count <- bench::mark(regression_data_100 <- count_footprints(sam_fname, transcript_length_fname, 
                                                                      offsets_fname, regression_data_100),
                              iterations=1)
print(paste("Time: ", time_100_count$min, " ; Memory: ", time_100_count$mem_alloc))
print(paste("Number of rows:", nrow(regression_data_100)))

# perform regression
print("100 genes: regression")
time_100_regression <- bench::mark(nb_fit_100 <- MASS::glm.nb(regression, data=regression_data_100),
                                   iterations=1)
print(paste("Time: ", time_100_regression$min, " ; Memory: ", time_100_regression$mem_alloc))

# save data
save(regression_data_100, file=file.path(outputDir, expt, "regression_data_100.Rda"))
save(nb_fit_100, file=file.path(outputDir, expt, "negbin_fit_100.Rda"))
save(regression_data_100, nb_fit_100, 
     time_100_init, time_100_count, time_100_regression,
     file=file.path(outputDir, expt, "data_100genes.Rda"))
rm(regression_data_100, nb_fit_100)

# 250 genes ---------------------------------------------------------------

n <- 250
expt <- paste0("genes_", n)
if(!(expt %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, expt)))
}

# select top 250 translated transcripts
transcripts_250 <- transcripts[1:n]
transcripts_seq_250 <- scer_transcriptome[match(names(transcripts_250), names(scer_transcriptome))]

# write transcriptome .fa file
transcripts_nCodons_250 <- scer_lengths$cds_length[match(names(transcripts_seq_250), scer_lengths$transcript)]/3
transcripts_seq_250 <- lapply(1:n,
                              function(x) {
                                # 6x 5' UTR codons + 5x 250' UTR codons
                                transcripts_seq_250[[x]][1:(transcripts_nCodons_250[x]+6+5)]
                              })
names(transcripts_seq_250) <- names(transcripts_250)
transcript_fa_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15.fa"))
writeTranscriptomeFA(transcripts_seq_250, transcript_fa_fname)

# write transcriptome lengths file
transcripts_lengths_250 <- scer_lengths[match(names(transcripts_250), scer_lengths$transcript),]
transcripts_lengths_250$utr5_length <- 18
transcripts_lengths_250$utr3_length <- 15
transcript_length_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_18cds15_lengths.txt"))
write.table(transcripts_lengths_250, transcript_length_fname, quote=F)

# create footprints FASTA file
transcripts_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_names.txt"))
sam_fname <- file.path(outputDir, expt, paste0("scer_", n, "_transcripts_footprints.sam"))
writeLines(names(transcripts_250), con=transcripts_fname)
system(paste("grep -wFf", transcripts_fname, weinberg_sam_fname, ">", sam_fname))

# create regression_data
print("250 genes: init data")
time_250_init <- bench::mark(regression_data_250 <- init_data(transcript_fa_fname, 
                                                              transcript_length_fname, bias_length=2),
                             iterations=1)
print(paste("Time: ", time_250_init$min, " ; Memory: ", time_250_init$mem_alloc))

# count footprints
print("250 genes: count footprints")
time_250_count <- bench::mark(regression_data_250 <- count_footprints(sam_fname, transcript_length_fname, 
                                                                      offsets_fname, regression_data_250),
                              iterations=1)
print(paste("Time: ", time_250_count$min, " ; Memory: ", time_250_count$mem_alloc))
print(paste("Number of rows:", nrow(regression_data_250)))

# perform regression
print("250 genes: regression")
time_250_regression <- bench::mark(nb_fit_250 <- MASS::glm.nb(regression, data=regression_data_250),
                                   iterations=1)
print(paste("Time: ", time_250_regression$min, " ; Memory: ", time_250_regression$mem_alloc))

# save data
save(regression_data_250, file=file.path(outputDir, expt, "regression_data_250.Rda"))
save(nb_fit_250, file=file.path(outputDir, expt, "negbin_fit_250.Rda"))
save(regression_data_250, nb_fit_250, 
     time_250_init, time_250_count, time_250_regression,
     file=file.path(outputDir, expt, "data_250genes.Rda"))
rm(regression_data_250, nb_fit_250)


# evaluation --------------------------------------------------------------

# knit markdown file

rmarkdown::render("scalability_eval.Rmd", output_format = "html_document")
