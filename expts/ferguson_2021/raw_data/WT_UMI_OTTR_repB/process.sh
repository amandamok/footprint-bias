#!/bin/bash

# cutadapt version 3.4
# umi_tools version 1.0.0
# bowtie version 1.0.0
# samtools version 1.7
## note: unset PYTHONPATH before running cutadapt

ref_dir="$HOME/footprint-bias/reference_data"
unset PYTHONPATH

# 1. trim constant adaptor
cat R1CY5Nb_round1.fastq.gz R1CY5Nb_round2.fastq.gz | zcat | \
cutadapt -j 10 -m 20 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--untrimmed-output ottrUMI_untrimmed.fq -o ottrUMI_trimmed_adaptor.fq - > \
ottrUMI_trimmed_adaptor.cutadapt

# 2. deduplicate: extract 6nt from 5' end as UMI and trim
umi_tools extract --extract-method=regex --bc-pattern="^(?P<umi_1>.{6}).+" \
-I ottrUMI_trimmed_adaptor.fq \
-S ottrUMI_trimmed_adaptor_umi.fq \
-L ottrUMI_trimmed_adaptor_umi.umitools

# 3. trim 1nt from 5' end --> NTA
cutadapt -j 10 -m 20 -u 1 --rename='{id} NTA={cut_prefix}' -o \
ottrUMI_trimmed_adaptor_umi_nta.fq ottrUMI_trimmed_adaptor_umi.fq > \
ottrUMI_trimmed_adaptor_umi_nta.cutadapt
## A: 3 330 267 (6.1%)
## C: 41 027 846 (75.8%)
## G: 4 361 232 (8.0%)
## T: 5 441 549 (10.0%)

# 4. trim 1nt from 3' end --> TPRT
cutadapt -j 10 -m 20 -u -1 --rename='{id} {comment} TPRT={cut_suffix}' -o \
ottrUMI_trimmed_adaptor_umi_nta_tprt.fq ottrUMI_trimmed_adaptor_umi_nta.fq > \
ottrUMI_trimmed_adaptor_umi_nta_tprt.cutadapt
## A: 39 535 664 (76.8%)
## C: 8 358 (0.0%)
## G: 11 917 248 (23.2%)
## T: 13 169 (0.0%)
## N: 192 (0.0%)

# 3. filter reads by length
cutadapt -m 27 -M 31 -j 10 -o ottrUMI_trimmed_adaptor_umi_nta_tprt_25to29.fq \
ottrUMI_trimmed_adaptor_umi_nta_tprt.fq > \
ottrUMI_trimmed_adaptor_umi_nta_tprt_25to29.cutadapt
## 51 474 631 input reads
## 27 979 235 (54.4%) reads too short
## 9 616 407 (18.7%) reads too long
## 13 878 989 (27.0%) reads written

# 5. remove contaminant reads: rRNA
bowtie -v 2 -p 10 -S --un ottrUMI_trimmed_not_rrna.fq \
$ref_dir/ScerRRNA ottrUMI_trimmed_adaptor_umi_nta_tprt_25to29.fq > \
ottrUMI_trimmed_rrna.sam 2> ottrUMI_trimmed_rrna.bowtiestats
## 13 878 989 input reads
## 9 692 641 (69.84%) reads with ≥1 alignment
## 4 186 348 (30.16%) unaligned reads

# 6. remove contaminant reads: ncRNA
bowtie -v 2 -p 10 -S --un ottrUMI_trimmed_not_rrna_ncrna.fq \
$ref_dir/rna_coding ottrUMI_trimmed_not_rrna.fq > \
ottrUMI_trimmed_ncrna.sam 2> ottrUMI_trimmed_ncrna.bowtiestats
## 4 186 348 input reads
## 36 350 (0.87%) reads with ≥1 alignment
## 4 149 998 (99.13%) unaligned reads

# 7. deduplicate: report all best alignments (-a --best --strata)
bowtie -v 2 -p 10 -S --norc -a --best --strata \
--un ottrUMI_trimmed_unaligned.fq \
$ref_dir/scer.transcripts.20cds20 \
ottrUMI_trimmed_not_rrna_ncrna.fq > \
ottrUMI_trimmed_not_rrna_ncrna_best.sam 2> \
ottrUMI_trimmed_not_rrna_ncrna_best.bowtiestats
## 3 933 044 (94.77%) reads with ≥1 alignment
## 216 954 (5.23%) unaligned reads
## 5 679 106 reported alignments

# 8. deduplicate: sort by read name
samtools view -h -F 0x04 ottrUMI_trimmed_not_rrna_ncrna_best.sam | \
samtools sort -n -O SAM -o ottrUMI_trimmed_not_rrna_ncrna_best_sorted.sam

# 9. deduplicate: choose positionally-first alignment per read
Rscript $HOME/footprint-bias/scripts/filter_sam.R \
-i ottrUMI_trimmed_not_rrna_ncrna_best_sorted.sam \
-o ottrUMI_trimmed_not_rrna_ncrna_best_sorted_filtered.sam

# 10. deduplicate: sort and index
samtools sort ottrUMI_trimmed_not_rrna_ncrna_best_sorted_filtered.sam \
-o ottrUMI_trimmed_not_rrna_ncrna_unique.bam -O BAM
samtools index ottrUMI_trimmed_not_rrna_ncrna_unique.bam

# 11. deduplicate: deduplicate UMIs
umi_tools dedup --output-stats=deduplicated --out-sam --read-length \
-I ottrUMI_trimmed_not_rrna_ncrna_unique.bam \
-S ottrUMI_trimmed_deduplicated.sam \
-L ottrUMI_trimmed_deduplicated.umitools \
-E ottrUMI_trimmed_deduplicated.error
## 3 933 044 reads in
## 808 567 positions deduplicated
## 2 283 876 reads out

# 12. deduplicate: convert into fastq
samtools fastq ottrUMI_trimmed_deduplicated.sam > ottrUMI_trimmed_deduplicated.fq

# 13. align to transcriptome
bowtie -v 2 -p 10 -S --norc -a $ref_dir/scer.transcripts.20cds20 \
ottrUMI_trimmed_deduplicated.fq > \
ottrUMI_trimmed_deduplicated_footprints.sam 2> \
ottrUMI_trimmed_deduplicated_footprints.bowtiestats
## 2 283 876 unique reads
## 3 620 934 alignments

# 14. compute multimapping weights
rsem-calculate-expression --sam \
ottrUMI_trimmed_deduplicated_footprints.sam \
$ref_dir/scer.transcripts.20cds20 \
ottrUMI_trimmed_deduplicated_footprints > \
ottrUMI_trimmed_deduplicated_footprints.rsem.stdout 2> \
ottrUMI_trimmed_deduplicated_footprints.rsem.stderr
