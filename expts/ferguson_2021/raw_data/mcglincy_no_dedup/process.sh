
# input: trimmed reads, filtered out rRNA and ncRNA

# align to transcriptome
bowtie -v 2 -p 10 -S --norc -a \
~/footprint-bias/reference_data/scer.transcripts.20cds20 \
../mcglincy/mcglincy_trimmed_not_rrna_ncrna.fq > \
mcglincy_trimmed_footprints.sam 2> \
mcglincy_trimmed_footprints.bowtiestats
## 17 437 864 input reads
## 8 434 340 reads w/ ≥1 alignment
## 65 260 589 output alignments

# compute RSEM weights
rsem-calculate-expression --seed-length 15 --sam \
mcglincy_trimmed_footprints.sam \
~/footprint-bias/reference_data/scer.transcripts.20cds20 \
mcglincy_trimmed_footprints > \
mcglincy_trimmed_footprints.rsem.stdout 2> \
mcglincy_trimmed_footprints.rsem.stderr
