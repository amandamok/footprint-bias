rm(list=ls())

fasta_file <- "sk1_revised_100331_chr_to_use.fa"

fasta_num_lines <- as.numeric(strsplit(system(paste("wc -l", fasta_file), intern=T), 
                                       split=" ")[[1]][1])

fasta_lines <- system(paste('grep -n ^">"', fasta_file), intern=T)
fasta_starts <- as.numeric(sapply(fasta_lines, function(x) { strsplit(x, split=":")[[1]][1] }))
names(fasta_starts) <- sapply(fasta_lines, function(x) strsplit(x, split=">")[[1]][2])
num_chr <- length(fasta_starts)
fasta_starts <- c(fasta_starts, fasta_num_lines+1)

for(x in seq(num_chr)) {
  tmp_start <- fasta_starts[x]
  tmp_stop <- fasta_starts[x+1] - 1
  tmp_name <- paste0("sk1.", names(fasta_starts)[x], ".fa.gz")
  system(paste0("sed -n '", tmp_start, ",", tmp_stop, "'p ", fasta_file,
                " | gzip > ", tmp_name))
}

q(save="no")
