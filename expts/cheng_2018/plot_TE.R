rm(list=ls())

# process RPF file
rpf <- readLines("GSE108778_timecourse_replicate_2_fps.txt")
rpf <- rpf[-length(rpf)] # last line is empty
rpf <- strsplit(rpf, split="\t")
rpf <- data.frame(do.call(rbind, rpf))
colnames(rpf) <- rpf[1, ]
colnames(rpf) <- sapply(colnames(rpf), function(x) strsplit(x, split=" ")[[1]][1])
rownames(rpf) <- rpf[, 1]
rpf <- rpf[-1, ]
rpf <- rpf[, -1]
for(x in colnames(rpf)) { rpf[[x]] <- as.numeric(rpf[[x]]) }

# process mRNA file
mrna <- readLines("GSE108778_timecourse_replicate_2_mRNA.txt")
mrna <- mrna[-length(mrna)] # last line is empty
mrna <- strsplit(mrna, split="\t")
mrna <- data.frame(do.call(rbind, mrna))
colnames(mrna) <- mrna[1, ]
colnames(mrna) <- sapply(colnames(mrna), function(x) strsplit(x, split="(-| )")[[1]][1])
colnames(mrna) <- sub("repl1", "", colnames(mrna))
colnames(mrna)[colnames(mrna)=="1.5"] <- "1.5hr"
rownames(mrna) <- mrna[, 1]
mrna <- mrna[-1, ]
mrna <- mrna[, -1]
for(x in colnames(mrna)) { mrna[[x]] <- as.numeric(mrna[[x]]) }
mrna <- mrna[match(rownames(rpf), rownames(mrna)), ]
mrna <- mrna[, match(colnames(rpf), colnames(mrna))]

# calculate TE
TE <- lapply(seq(ncol(mrna)), function(x) rpf[, x] / mrna[, x])
TE <- data.frame(do.call(cbind, TE))
rownames(TE) <- rownames(mrna)
colnames(TE) <- colnames(mrna)

genes_TE <- data.frame(timepoint=rep(colnames(TE), times=2),
                       TE=c(unlist(TE[rownames(TE)=="YPL119C",]),
                            unlist(TE[rownames(TE)=="YOR204W",])),
                       gene=rep(c("Dbp1", "Ded1"), each=ncol(TE)))
genes_TE$timepoint <- factor(genes_TE$timepoint, levels=colnames(TE))
(ggplot(genes_TE, aes(x=timepoint, y=TE, fill=gene)) + 
  geom_col() + theme_bw() + facet_grid(~gene) + xlab("") + guides(fill="none")) / 
  (ggplot(subset(genes_TE, timepoint %in% c("4.5hr", "6hr")), aes(x=gene, y=TE, fill=gene)) + 
     geom_col() + theme_bw() + facet_grid(~timepoint) + xlab("") + guides(fill="none"))

