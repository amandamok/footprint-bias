#!/bin/bash

ref_index="../../reference_data/grch38.transcripts.20cds20"

# re-align to transcriptome
echo "WT: aligning to transcriptome"
bowtie -v 2 -p 10 -S --norc -a ${ref_index} \
../../raw_data/WT_RPF/WT_RPF_trim3_not_rrna_trna.fq > \
WT_RPF_trim3_footprints.sam 2> WT_RPF_trim3_footprints.bowtiestats

# re-compute multimapping weights
echo "WT: computing multimapping weights"
rsem-calculate-expression --seed-length 22 --sam WT_RPF_trim3_footprints.sam ${ref_index} WT_RPF_trim3_footprints > WT_RPF_trim3_footprints.rsem.stdout 2> WT_RPF_trim3_footprints.rsem.stderr

# convert .bam to .sam
echo "WT: converting from bam to sam"
samtools view -h WT_RPF_trim3_footprints.transcript.bam > \
WT_RPF_trim3_footprints.transcripts.sam
