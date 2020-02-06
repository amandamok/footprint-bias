rm(list=ls())

source("~/footprint-bias/scripts/prep_data.R")

data_dir <- "~/footprint-bias/sample_data"

rpf_fa_fname <- file.path(data_dir, "riboViz_simASite.fa")
rpf_sam_fname <- file.path(data_dir, "riboViz_simASite_footprints.sam")
transcripts_fa_fname <- file.path(data_dir, "riboViz_transcripts_18cds15.fa")
transcripts_length_fname <- file.path(data_dir, "riboViz_transcripts_18cds15_lengths.txt")
Asite_rules_fname <- file.path(data_dir, "Asite_rules.txt")

# initialize data frame ---------------------------------------------------

Asite_rules <- read.table(Asite_rules_fname)
d5_range <- range(unlist(Asite_rules), na.rm=T)
rpf_lengths <- range(as.numeric(rownames(Asite_rules)))
d3_range <- range(t(as.numeric(rownames(Asite_rules))) - Asite_rules - 3, na.rm=T)
rpf_counts <- init_data(transcripts_fa_fname, transcripts_length_fname,
                        digest5_lengths=d5_range[1]:d5_range[2],
                        digest3_lengths=d3_range[1]:d3_range[2],
                        bias_length=2)

# count footprints (by annotation) ----------------------------------------

# load read ids (contain annotation information)
rpf_ids <- readLines(rpf_fa_fname)
rpf_ids <- rpf_ids[(1:(length(rpf_ids)/2))*2-1]
rpf_ids <- sub(">", "", rpf_ids)

# split read ids
rpf_anno <- data.frame(do.call(rbind, strsplit(rpf_ids, split="_"))[,1:4])
colnames(rpf_anno) <- c("transcript", "cod_idx", "d5", "d3")

# update cod_idx (rpf_counts$cod_idx are 1-based; simulations are 0-based)
rpf_anno$cod_idx <- as.numeric(as.character(rpf_anno$cod_idx)) + 1

# count up footprints with same (transcript, cod_idx, d5, d3)
rpf_anno <- aggregate(list(count=rep(1, nrow(rpf_anno))), rpf_anno, length)

# add counts to initialized data frame
rpf_rows <- prodlim::row.match(rpf_anno[, c("transcript", "cod_idx", "d5", "d3")],
                               rpf_counts[, c("transcript", "cod_idx", "d5", "d3")])
rpf_counts_true <- rpf_counts
rpf_counts_true$count[rpf_rows] <- rpf_anno$count

# count of footprints (mapped) --------------------------------------------

rpf_counts_mapped <- count_footprints(rpf_sam_fname, transcripts_length_fname,
                                      Asite_rules_fname, rpf_counts)


# check if annotations and mapped counts match ----------------------------

true_mapped_index <- prodlim::row.match(rpf_counts_true[, c("transcript", "cod_idx", "d5", "d3")],
                                        rpf_counts_mapped[, c("transcript", "cod_idx", "d5", "d3")])
summary(rpf_counts_true$count == rpf_counts_mapped$count[true_mapped_index])

# write counts data frame -------------------------------------------------

rpf_counts_codon <- aggregate(count ~ transcript + cod_idx, rpf_counts_mapped, sum)

write.table(rpf_counts_mapped, 
            file=file.path(data_dir, "riboViz_simASite_countsTable.txt"),
            quote=F, row.names=F, col.names=T)
write.table(rpf_counts_codon,
            file=file.path(data_dir, "riboViz_simASite_countsPerPosition.txt"),
            quote=F, row.names=F, col.names=T)
