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
cutadapt -a CGTAAAGATCGGAAGAGCAC -j 10 -m 15 \
--untrimmed-output mcglincy_untrimmed.fq \
-o mcglincy_trimmed_adaptor.fq - > mcglincy_trimmed_adaptor.cutadapt

# 2. preprocess: trim 2nt from 5' end (from cDNA ligation)
cutadapt -j 10 -u 2 -o mcglincy_trimmed_adaptor_cDNA.fq \
mcglincy_trimmed_adaptor.fq > mcglincy_trimmed_adaptor_cDNA.cutadapt

# 3. deduplicate: extract 5nt from 3' end as UMI and trim
umi_tools extract --extract-method=regex --3prime --bc-pattern=".+(?P<umi_1>.{5})$" \
-I mcglincy_trimmed_adaptor_cDNA.fq \
-S mcglincy_trimmed_adaptor_cDNA_umi.fq \
-L mcglincy_trimmed_adaptor_cDNA_umi.umitools

# 4. remove contaminant reads: rRNA
bowtie -v 2 -p 10 -S --un mcglincy_trimmed_not_rrna.fq \
$ref_dir/ScerRRNA mcglincy_trimmed_adaptor_cDNA_umi.fq > \
mcglincy_trimmed_rrna.sam 2> mcglincy_trimmed_rrna.bowtiestats
## 98 433 619 input reads
## 79 245 242 (80.51%) reads with ≥1 alignment
## 19 188 377 (19.49%) unaligned reads

# 5. remove contaminant reads: ncRNA
bowtie -v 2 -p 10 -S --un mcglincy_trimmed_not_rrna_ncrna.fq \
$ref_dir/rna_coding mcglincy_trimmed_not_rrna.fq > \
mcglincy_trimmed_ncrna.sam 2> mcglincy_trimmed_ncrna.bowtiestats
## 19 188 377 input reads
## 1 750 513 (9.12%) reads with ≥1 alignment
## 17 437 864 (90.88%) unaligned reads

# 6. deduplicate: report all best alignments (-a --best --strata)
bowtie -v 2 -p 10 -S --norc -a --best --strata \
--un mcglincy_trimmed_unaligned.fq \
$ref_dir/scer.transcripts.20cds20 \
mcglincy_trimmed_not_rrna_ncrna.fq > \
mcglincy_trimmed_not_rrna_ncrna_best.sam 2> \
mcglincy_trimmed_not_rrna_ncrna_best.bowtiestats
## 8 434 340 (48.37%) reads with ≥1 alignment
## 9 003 524 (51.63%) unaligned reads
## 16 832 967 reported alignments

# 7. deduplicate: sort by read name
samtools view -h -F 0x04 mcglincy_trimmed_not_rrna_ncrna_best.sam | \
samtools sort -n -O SAM -o mcglincy_trimmed_not_rrna_ncrna_best_sorted.sam

# 8. deduplicate: choose positionally-first alignment per read
Rscript $HOME/footprint-bias/scripts/filter_sam.R \
-i mcglincy_trimmed_not_rrna_ncrna_best_sorted.sam \
-o mcglincy_trimmed_not_rrna_ncrna_best_sorted_filtered.sam

# 9. deduplicate: sort and index
samtools sort mcglincy_trimmed_not_rrna_ncrna_best_sorted_filtered.sam \
-o mcglincy_trimmed_not_rrna_ncrna_unique.bam -O BAM
samtools index mcglincy_trimmed_not_rrna_ncrna_unique.bam

# 10. deduplicate: deduplicate UMIs
umi_tools dedup --output-stats=deduplicated --in-sam --out-sam --read-length \
-I mcglincy_trimmed_not_rrna_ncrna_unique.bam \
-S mcglincy_trimmed_deduplicated.sam \
-L mcglincy_trimmed_deduplicated.umitools \
-E mcglincy_trimmed_deduplicated.error
## 8 434 340 reads in
## 1 947 395 positions deduplicated
## 3 479 630 reads out

# 11. deduplicate: convert into fastq
samtools fastq mcglincy_trimmed_deduplicated.sam > mcglincy_trimmed_deduplicated.fq

# 12. align to transcriptome
bowtie -v 2 -p 10 -S --norc -a $ref_dir/scer.transcripts.20cds20 \
mcglincy_trimmed_deduplicated.fq > \
mcglincy_trimmed_deduplicated_footprints.sam 2> \
mcglincy_trimmed_deduplicated_footprints.bowtiestats
## 3 479 630 unique reads
## 28 757 047 alignments

# 13. compute multimapping weights
rsem-calculate-expression --seed-length 15 --sam \
mcglincy_trimmed_deduplicated_footprints.sam \
$ref_dir/scer.transcripts.20cds20 \
mcglincy_trimmed_deduplicated_footprints > \
mcglincy_trimmed_deduplicated_footprints.rsem.stdout 2> \
mcglincy_trimmed_deduplicated_footprints.rsem.stderr
