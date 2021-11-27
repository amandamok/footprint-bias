#!/bin/bash

# cutadapt version 3.4
# umi_tools version 1.0.0
# bowtie version 1.0.0
# samtools version 1.7
## note: unset PYTHONPATH before running cutadapt 

ref_dir="$HOME/footprint-bias/reference_data"
unset PYTHONPATH

# 1. preprocess: trim constant 3' adaptor
cat McGlincy_round1.fastq.gz McGlincy_round2.fastq.gz | zcat | \
cutadapt -a CGTAAAGATCGGAAGAGCAC -j 10 -m 20 \
--untrimmed-output mcglincy_untrimmed.fq \
-o mcglincy_trimmed_adaptor.fq - > mcglincy_trimmed_adaptor.cutadapt

# 2. preprocess: extract 2nt from 5' end and 5nt from 3' end as UMI; trim
umi_tools extract --extract-method=regex \
--bc-pattern="(?P<umi_1>.{2})(.+)(?P<umi_2>.{5})$" \
-I mcglincy_trimmed_adaptor.fq \
-S mcglincy_trimmed_adaptor_umi.fq \
-L mcglincy_trimmed_adaptor_umi.umitools

# 3. remove contaminant reads: rRNA
bowtie -v 2 -p 10 -S --un mcglincy_trimmed_not_rrna.fq \
$ref_dir/ScerRRNA mcglincy_trimmed_adaptor_umi.fq > \
mcglincy_trimmed_rrna.sam 2> mcglincy_trimmed_rrna.bowtiestats
## 78 342 763 input reads
## 59 351 550 (75.76%) reads with ≥1 alignment
## 18 991 213 (24.24%) unaligned reads

# 4. remove contaminant reads: ncRNA
bowtie -v 2 -p 10 -S --un mcglincy_trimmed_not_rrna_ncrna.fq \
$ref_dir/rna_coding mcglincy_trimmed_not_rrna.fq > \
mcglincy_trimmed_ncrna.sam 2> mcglincy_trimmed_ncrna.bowtiestats
## 18 991 213 input reads
## 1 579 426 (8.32%) reads with ≥1 alignment
## 17 411 787 (91.68%) unaligned reads

# 5. deduplicate: report all best alignments (-a --best --strata)
bowtie -v 2 -p 10 -S --norc -a --best --strata \
--un mcglincy_trimmed_unaligned.fq \
$ref_dir/scer.transcripts.20cds20 \
mcglincy_trimmed_not_rrna_ncrna.fq > \
mcglincy_trimmed_not_rrna_ncrna_best.sam 2> \
mcglincy_trimmed_not_rrna_ncrna_best.bowtiestats
## 17 411 787 input reads
## 8 408 263 (48.29%) reads with ≥1 alignment
## 9 003 524 (51.71%) unaligned reads
## 16 642 462 reported alignments

# 6. deduplicate: sort by read name
samtools view -h -F 0x04 mcglincy_trimmed_not_rrna_ncrna_best.sam | \
samtools sort -n -O SAM -o mcglincy_trimmed_not_rrna_ncrna_best_sorted.sam

# 7. deduplicate: choose positionally-first alignment per read
Rscript $HOME/footprint-bias/scripts/filter_sam.R \
-i mcglincy_trimmed_not_rrna_ncrna_best_sorted.sam \
-o mcglincy_trimmed_not_rrna_ncrna_best_sorted_filtered.sam

# 8. deduplicate: sort and index
samtools sort mcglincy_trimmed_not_rrna_ncrna_best_sorted_filtered.sam \
-o mcglincy_trimmed_not_rrna_ncrna_unique.bam -O BAM
samtools index mcglincy_trimmed_not_rrna_ncrna_unique.bam

# 9. deduplicate: deduplicate UMIs
umi_tools dedup --output-stats=deduplicated --in-sam --out-sam --read-length \
-I mcglincy_trimmed_not_rrna_ncrna_unique.bam \
-S mcglincy_trimmed_deduplicated.sam \
-L mcglincy_trimmed_deduplicated.umitools \
-E mcglincy_trimmed_deduplicated.error
## 8 408 263 reads in
## 1 940 063 positions deduplicated
## 3 906 538 reads out

# 10. deduplicate: convert into fastq
samtools fastq mcglincy_trimmed_deduplicated.sam > mcglincy_trimmed_deduplicated.fq

# 11. align to transcriptome
bowtie -v 2 -p 10 -S --norc -a $ref_dir/scer.transcripts.20cds20 \
mcglincy_trimmed_deduplicated.fq > \
mcglincy_trimmed_deduplicated_footprints.sam 2> \
mcglincy_trimmed_deduplicated_footprints.bowtiestats
## 3 906 538 unique reads
## 27 819 829 alignments

# 12. compute multimapping weights
rsem-calculate-expression --seed-length 15 --sam \
mcglincy_trimmed_deduplicated_footprints.sam \
$ref_dir/scer.transcripts.20cds20 \
mcglincy_trimmed_deduplicated_footprints > \
mcglincy_trimmed_deduplicated_footprints.rsem.stdout 2> \
mcglincy_trimmed_deduplicated_footprints.rsem.stderr

