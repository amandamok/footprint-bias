dssp <- readLines("~/footprint-bias/expts/ferguson_2021/YKL060C.dssp")
dssp <- subset(dssp, seq(length(dssp)) > grep("#", dssp))
dssp <- sub("^( )+", "", dssp)
dssp <- lapply(dssp, function(x) { strsplit(x, split="( )+")[[1]][1:5] })
dssp <- data.frame(do.call(rbind, dssp), stringsAsFactors=F)
colnames(dssp) <- c("row_number", "index", "chain", "aa", "structure")
structure_alphabet <- c(H="alpha helix",
                        B="beta bridge",
                        E="strand",
                        G="helix 3",
                        I="helix 5",
                        P="helix PPII",
                        "T"="turn",
                        S="bend")
dssp$structure[!(dssp$structure %in% names(structure_alphabet))] <- "L"
structure_alphabet <- c(structure_alphabet, L="loop")
dssp$structure_label <- structure_alphabet[dssp$structure]
dssp$helix_number <- NA
helix_number <- 1
for(x in seq(nrow(dssp))) {
  if(dssp$structure_label[x] == "alpha helix") {
    if(dssp$structure_label[x-1] != "alpha helix") {
      helix_number <- helix_number + 1
    }
    dssp$helix_number[x] <- paste0("helix_", helix_number)
  } else {
    dssp$helix_number[x] <- "other"
  }
}

dss <- read.table("~/footprint-bias/expts/ferguson_2021/YKL060C.dss",
                  sep=":", col.names=c("index", "structure"))

YKL060C$structure <- dssp$structure[match(YKL060C$cod_idx, dssp$index)]
YKL060C$structure_label <- dssp$structure_label[match(YKL060C$cod_idx, dssp$index)]
YKL060C$structure_label[is.na(YKL060C$structure_label)] <- "stop codon"
YKL060C$helix_number <- dssp$helix_number[match(YKL060C$cod_idx, dssp$index)]
YKL060C$helix_number[is.na(YKL060C$helix_number)] <- "other"
levels(YKL060C$variable) <- c("raw", "corrected")
YKL060C$dss_structure <- dss$structure[match(YKL060C$cod_idx, dss$index)]
YKL060C$dss_structure_label <- c(H="helix", L="loop", S="sheet")[YKL060C$dss_structure]

ggplot(YKL060C, aes(x=cod_idx, y=value, fill=dss_structure_label)) +
  geom_col() + theme_bw() + facet_grid(variable ~.) + labs(fill="dss structure") +
  xlab("codon") + ylab("RPF count") + ggtitle("YKL060C", subtitle="McGlincy") +
  # scale_fill_discrete(type=c(RColorBrewer::brewer.pal(12, "Paired"), "grey")) +
  geom_hline(data=aggregate(value ~ variable, YKL060C, median),
             aes(yintercept=value), linetype="dashed")


xx <- data.frame(common = unlist(mapply(rep, names(xx), lengths(xx))),
                 ORF = unlist(xx))
