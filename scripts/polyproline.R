rm(list=ls())

library(here)
library(choros)

source(file.path(here(), "scripts", "helper.R"))

transcripts_aa_fname <- file.path(here(), "reference_data", "scer.transcripts.aa.fa")
transcripts_aa <- load_fa(transcripts_aa_fname)

window_width <- 10


# identify polyproline (>3) regions ---------------------------------------

# subset to transcripts with ≥3 consecutive prolines
polyproline <- transcripts_aa[sapply(transcripts_aa, function(x) grepl("PPP", x))]

# identify codon positions of poly-proline tracts
polyproline_regions <- lapply(seq_along(polyproline),
                              function(gene) {
                                aa <- strsplit(polyproline[gene], split="")[[1]]
                                in_PPP <- which(sapply(seq.int(length(aa)-2),
                                                       function(pos) {
                                                         all(aa[pos:(pos+2)]=="P")
                                                       }))
                                # report polyproline regions
                                PPP <- data.frame(transcript=names(polyproline)[gene],
                                                  start=in_PPP[1],
                                                  end=in_PPP[1]+2)
                                if(length(in_PPP)==1) {
                                  return(PPP)
                                } else {
                                  for(i in 2:length(in_PPP)) {
                                    if(in_PPP[i] <= PPP$end[nrow(PPP)]) {
                                      PPP$end[nrow(PPP)] <- in_PPP[i]+2
                                    } else {
                                      PPP <- rbind(PPP,
                                                   data.frame(transcript=names(polyproline)[gene],
                                                              start=in_PPP[i],
                                                              end=in_PPP[i]+2))
                                    }
                                  }
                                  return(PPP)
                                }
                              })
polyproline_regions <- do.call(rbind, polyproline_regions)
polyproline_regions$length <- with(polyproline_regions, end-start+1)
polyproline_regions$transcript <- as.character(polyproline_regions$transcript)


# load experimental data --------------------------------------------------

# green
load(file.path(here(), "expts", "allModels", "green_20cds20_f5_2_f3_3", "green_bam.Rda"))
green_bam <- subset(green_bam, transcript %in% names(polyproline))
green_dat_raw <- data.frame(aggregate(count ~ transcript + cod_idx,
                                      data=green_bam, FUN=sum),
                            type="raw")
green_dat_corrected <- data.frame(aggregate(correct_150 ~ transcript + cod_idx,
                                            data=green_bam, FUN=sum),
                                  type="corrected")
colnames(green_dat_corrected)[colnames(green_dat_corrected)=="correct_150"] <- "count"
green_dat <- data.frame(rbind(green_dat_raw, green_dat_corrected),
                        expt="green")
rm(green_bam, green_dat_raw, green_dat_corrected)

# lareau
load(file.path(here(), "expts", "allModels", "lareau_20cds20_f5_2_f3_3", "lareau_bam.Rda"))
lareau_bam <- subset(lareau_bam, transcript %in% names(polyproline))
lareau_dat_raw <- data.frame(aggregate(count ~ transcript + cod_idx,
                                       data=lareau_bam, FUN=sum),
                             type="raw")
lareau_dat_corrected <- data.frame(aggregate(correct_150 ~ transcript + cod_idx,
                                             data=lareau_bam, FUN=sum),
                                   type="corrected")
colnames(lareau_dat_corrected)[colnames(lareau_dat_corrected)=="correct_150"] <- "count"
lareau_dat <- data.frame(rbind(lareau_dat_raw, lareau_dat_corrected),
                         expt="lareau")
rm(lareau_bam, lareau_dat_raw, lareau_dat_corrected)

# weinberg
load(file.path(here(), "expts", "allModels", "weinberg_20cds20_f5_2_f3_3", "weinberg_bam.Rda"))
weinberg_bam <- subset(weinberg_bam, transcript %in% names(polyproline))
weinberg_dat_raw <- data.frame(aggregate(count ~ transcript + cod_idx,
                                         data=weinberg_bam, FUN=sum),
                               type="raw")
weinberg_dat_corrected <- data.frame(aggregate(correct_150 ~ transcript + cod_idx,
                                               data=weinberg_bam, FUN=sum),
                                     type="corrected")
colnames(weinberg_dat_corrected)[colnames(weinberg_dat_corrected)=="correct_150"] <- "count"
weinberg_dat <- data.frame(rbind(weinberg_dat_raw, weinberg_dat_corrected),
                           expt="weinberg")
rm(weinberg_bam, weinberg_dat_raw, weinberg_dat_corrected)

# aggregate
all_dat <- rbind(green_dat, lareau_dat, weinberg_dat)
all_dat$transcript <- as.character(all_dat$transcript)
rm(green_dat, lareau_dat, weinberg_dat)


# count reads per region --------------------------------------------------

region_dat <- polyproline_regions
for(x in expts) {
  for(y in c("raw", "corrected")) {
    print(paste(x, y))
    region_dat[paste0(x, "_", y)] <- sapply(seq(nrow(polyproline_regions)),
                                            function(i) {
                                              tmp <- subset(all_dat, expt==x &
                                                              type==y &
                                                              transcript == polyproline_regions$transcript[i] &
                                                              cod_idx >= polyproline_regions$start[i] &
                                                              cod_idx <= polyproline_regions$end[i])
                                              return(sum(tmp$count))
                                            })
  }
}
window_dat <- polyproline_regions
for(x in expts) {
  for(y in c("raw", "corrected")) {
    print(paste(x, y))
    window_dat[paste0(x, "_", y)] <- sapply(seq(nrow(polyproline_regions)),
                                            function(i) {
                                              tmp <- subset(all_dat, expt==x &
                                                              type==y &
                                                              transcript == polyproline_regions$transcript[i] &
                                                              cod_idx >= polyproline_regions$start[i]  - window_width&
                                                              cod_idx <= polyproline_regions$end[i] + window_width)
                                              return(sum(tmp$count))
                                            })
  }
}
polyproline_regions$region_avg <- rowMeans(region_dat[,5:ncol(region_dat)])
polyproline_regions$region_diff <- rowMeans(data.frame(with(region_dat, abs(green_corrected-green_raw)),
                                                       with(region_dat, abs(lareau_corrected-lareau_raw)),
                                                       with(region_dat, abs(weinberg_corrected-weinberg_raw))))
polyproline_regions$window_avg <- rowMeans(window_dat[,5:ncol(window_dat)])
polyproline_regions$window_diff <- rowMeans(data.frame(with(window_dat, abs(green_corrected-green_raw)),
                                                       with(window_dat, abs(lareau_corrected-lareau_raw)),
                                                       with(window_dat, abs(weinberg_corrected-weinberg_raw))))
polyproline_regions <- polyproline_regions[with(polyproline_regions, order(region_diff/region_avg, decreasing=T)),]

# plot individual regions -------------------------------------------------

best_PPP <- subset(polyproline_regions, region_avg>150)
best_PPP_plots <- lapply(31:45,
                         function(which_PPP) {
                           ggplot(subset(all_dat, transcript==best_PPP$transcript[which_PPP] &
                                           cod_idx >= best_PPP$start[which_PPP]-window_width &
                                           cod_idx <= best_PPP$end[which_PPP]+window_width),
                                  aes(x=cod_idx, y=count, col=type)) +
                             geom_line() + facet_grid(expt ~ .) +
                             geom_vline(xintercept=best_PPP$start[which_PPP], linetype=2) +
                             geom_vline(xintercept=best_PPP$end[which_PPP], linetype=2) +
                             theme_classic() + scale_color_manual(values=1:2) +
                             ggtitle(paste0(best_PPP$transcript[which_PPP], ": ",
                                            best_PPP$start[which_PPP], "-", best_PPP$end[which_PPP]))
                         })
wrap_plots(best_PPP_plots, nrow=3)

# metagene plot -----------------------------------------------------------

metagene <- rbind(data.frame(pos=seq(-10, 10),
                             count=colSums(t(sapply(seq(nrow(polyproline_regions)),
                                                    function(x) {
                                                      tmp_transcript <- polyproline_regions$transcript[x]
                                                      tmp_end <- polyproline_regions$end[x]
                                                      tmp_window <- seq(tmp_end-window_width,
                                                                        tmp_end+window_width)
                                                      tmp_dat <- subset(all_dat,
                                                                        transcript == tmp_transcript &
                                                                          cod_idx %in% tmp_window &
                                                                          type == "raw" & expt == "green")
                                                      tmp_cts <- tmp_dat$count[match(tmp_window, tmp_dat$cod_idx)]
                                                      return(tmp_cts)
                                                    })), na.rm=T),
                             expt="green", type="raw"),
                  data.frame(pos=seq(-10, 10),
                             count=colSums(t(sapply(seq(nrow(polyproline_regions)),
                                                    function(x) {
                                                      tmp_transcript <- polyproline_regions$transcript[x]
                                                      tmp_end <- polyproline_regions$end[x]
                                                      tmp_window <- seq(tmp_end-window_width,
                                                                        tmp_end+window_width)
                                                      tmp_dat <- subset(all_dat,
                                                                        transcript == tmp_transcript &
                                                                          cod_idx %in% tmp_window &
                                                                          type == "corrected" & expt == "green")
                                                      tmp_cts <- tmp_dat$count[match(tmp_window, tmp_dat$cod_idx)]
                                                      return(tmp_cts)
                                                    })), na.rm=T),
                             expt="green", type="corrected"),
                  data.frame(pos=seq(-10, 10),
                             count=colSums(t(sapply(seq(nrow(polyproline_regions)),
                                                    function(x) {
                                                      tmp_transcript <- polyproline_regions$transcript[x]
                                                      tmp_end <- polyproline_regions$end[x]
                                                      tmp_window <- seq(tmp_end-window_width,
                                                                        tmp_end+window_width)
                                                      tmp_dat <- subset(all_dat,
                                                                        transcript == tmp_transcript &
                                                                          cod_idx %in% tmp_window &
                                                                          type == "raw" & expt == "lareau")
                                                      tmp_cts <- tmp_dat$count[match(tmp_window, tmp_dat$cod_idx)]
                                                      return(tmp_cts)
                                                    })), na.rm=T),
                             expt="lareau", type="raw"),
                  data.frame(pos=seq(-10, 10),
                             count=colSums(t(sapply(seq(nrow(polyproline_regions)),
                                                    function(x) {
                                                      tmp_transcript <- polyproline_regions$transcript[x]
                                                      tmp_end <- polyproline_regions$end[x]
                                                      tmp_window <- seq(tmp_end-window_width,
                                                                        tmp_end+window_width)
                                                      tmp_dat <- subset(all_dat,
                                                                        transcript == tmp_transcript &
                                                                          cod_idx %in% tmp_window &
                                                                          type == "corrected" & expt == "lareau")
                                                      tmp_cts <- tmp_dat$count[match(tmp_window, tmp_dat$cod_idx)]
                                                      return(tmp_cts)
                                                    })), na.rm=T),
                             expt="lareau", type="corrected"),
                  data.frame(pos=seq(-10, 10),
                             count=colSums(t(sapply(seq(nrow(polyproline_regions)),
                                                    function(x) {
                                                      tmp_transcript <- polyproline_regions$transcript[x]
                                                      tmp_end <- polyproline_regions$end[x]
                                                      tmp_window <- seq(tmp_end-window_width,
                                                                        tmp_end+window_width)
                                                      tmp_dat <- subset(all_dat,
                                                                        transcript == tmp_transcript &
                                                                          cod_idx %in% tmp_window &
                                                                          type == "raw" & expt == "weinberg")
                                                      tmp_cts <- tmp_dat$count[match(tmp_window, tmp_dat$cod_idx)]
                                                      return(tmp_cts)
                                                    })), na.rm=T),
                             expt="weinberg", type="raw"),
                  data.frame(pos=seq(-10, 10),
                             count=colSums(t(sapply(seq(nrow(polyproline_regions)),
                                                    function(x) {
                                                      tmp_transcript <- polyproline_regions$transcript[x]
                                                      tmp_end <- polyproline_regions$end[x]
                                                      tmp_window <- seq(tmp_end-window_width,
                                                                        tmp_end+window_width)
                                                      tmp_dat <- subset(all_dat,
                                                                        transcript == tmp_transcript &
                                                                          cod_idx %in% tmp_window &
                                                                          type == "corrected" & expt == "weinberg")
                                                      tmp_cts <- tmp_dat$count[match(tmp_window, tmp_dat$cod_idx)]
                                                      return(tmp_cts)
                                                    })), na.rm=T),
                             expt="weinberg", type="corrected"))
ggplot(metagene, aes(x=pos, y=count, col=type)) +
  geom_line() + facet_grid(expt ~ .) + geom_vline(xintercept=0, linetype=2) +
  theme_classic() + scale_color_manual(values=1:2) +
  ggtitle("metagene plot of polyproline (≥3) regions") +
  xlab("codon position, relative to 3'-most proline residue")

# save data ---------------------------------------------------------------

save(all_dat, file=file.path(here(), "scripts", "polyproline_all_dat.Rda"))
save(metagene, file=file.path(here(), "scripts", "polyproline_metagene.Rda"))
save(polyproline_regions, file=file.path(here(), "scripts", "polyproline_regions.Rda"))
