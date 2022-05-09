#!/bin/bash

unset PYTHONPATH

ref_dir="$HOME/footprint-bias/reference_data"

# 1. merge replicates
cat SRR*fastq > disome_raw.fq

# 2. trim UMIs
umi_tools extract --extract-method=regex \
--bc-pattern="(?P<umi_1>.{2})(.+)(?P<umi_2>.{5})$" \
-I disome_raw.fq -S disome_trimUMI.fq -L disome_trimUMI.umitools

# 4. remove contaminant reads: rRNA
bowtie -v 2 -p 10 -S --un disome_trimUMI_not_rrna.fq ${ref_dir}/ScerRRNA \
disome_trimUMI.fq > disome_trimUMI_rrna.sam 2> \
disome_trimUMI_rrna.bowtiestats

# 5. remove contaminant reads: ncRNA
bowtie -v 2 -p 10 --un disome_trimUMI_not_rrna_ncrna.fq ${ref_dir}/rna_coding \
disome_trimUMI_not_rrna.fq > disome_trimUMI_ncrna.sam 2> \
disome_trimUMI_ncrna.bowtiestats

# 6. deduplicate: report all best alignments (-a --best --strata)
bowtie -v 2 -p 10 -S --norc -a --best --strata --un \
disome_trimUMI_unaligned.fq ${ref_dir}/scer.transcripts.20cds20 \
disome_trimUMI_not_rrna_ncrna.fq > disome_trimUMI_not_rrna_ncrna_best.sam 2> \
disome_trimUMI_not_rrna_ncrna_best.bowtiestats

# 7. deduplicate: sort by read name
samtools view -h -F 0x04 disome_trimUMI_not_rrna_ncrna_best.sam | \
samtools sort -n -O SAM -o disome_trimUMI_not_rrna_ncrna_best_sorted.sam

# 8. deduplicate: choose positionally-first alignment per reads
Rscript ${HOME}/footprint-bias/scripts/filter_sam.R \
-i disome_trimUMI_not_rrna_ncrna_best_sorted.sam \
-o disome_trimUMI_not_rrna_ncrna_best_sorted_filtered.sam

# 9. deduplicate: sort and index
samtools sort disome_trimUMI_not_rrna_ncrna_best_sorted_filtered.sam \
-o disome_trimUMI_not_rrna_ncrna_unique.bam -O BAM

samtools index disome_trimUMI_not_rrna_ncrna_unique.bam

# 10. deduplicate: deduplicate UMIs
umi_tools dedup --output-stats=deduplicated --in-sam --out-sam --read-length \
-I disome_trimUMI_not_rrna_ncrna_unique.bam \
-S disome_trimUMI_deduplicated.sam \
-L disome_trimUMI_deduplicated.umitools \
-E disome_trimUMI_deduplicated.error

# 11. deduplicate: convert into fastq
samtools fastq disome_trimUMI_deduplicated.sam > disome_trimUMI_deduplicated.fq

# 6. align to transcriptome
bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/scer.transcripts.20cds20 \
disome_trimUMI_deduplicated.fq > disome_trimUMI_footprints.sam 2> \
disome_trimUMI_footprints.bowtiestats

# 7. compute multimapping weights
rsem-calculate-expression --sam disome_trimUMI_footprints.sam \
${ref_dir}/scer.transcripts.20cds20 disome_trimUMI_footprints > \
disome_trimUMI_footprints.rsem.stdout 2> disome_trimUMI_footprints.rsem.stderr
