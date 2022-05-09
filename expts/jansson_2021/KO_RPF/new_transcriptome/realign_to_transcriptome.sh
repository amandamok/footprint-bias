#!/bin/bash

ref_index="../../reference_data/grch38.transcripts.20cds20"

# re-align to transcriptome
echo "KO: aligning to transcriptome"
bowtie -v 2 -p 10 -S --norc -a ${ref_index} \
../../raw_data/KO_RPF/KO_RPF_trim3_not_rrna_trna.fq > \
KO_RPF_trim3_footprints.sam 2> KO_RPF_trim3_footprints.bowtiestats

# re-compute multimapping weights
echo "KO: computing multimapping weights"
rsem-calculate-expression --seed-length 22 --sam KO_RPF_trim3_footprints.sam ${ref_index} KO_RPF_trim3_footprints > KO_RPF_trim3_footprints.rsem.stdout 2> KO_RPF_trim3_footprints.rsem.stderr

# convert .bam to .sam
echo "KO: converting from bam to sam"
samtools view -h KO_RPF_trim3_footprints.transcript.bam > \
KO_RPF_trim3_footprints.transcripts.sam
