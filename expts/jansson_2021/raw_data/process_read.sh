#!/bin/bash

unset PYTHONPATH

ref_dir="$HOME/footprint-bias/reference_data"
min_length=22

for expt in $(ls -d *RPF)
do
  echo ${expt}
  cd ${expt}
  
  # 1. pool replicates
  echo "... pooling replicates"
  cat SRR*fastq > ${expt}_raw.fq

  # 2. trim 3' adapter
  echo "... trimming 3' adapter"
  cutadapt -a AGATCGGAAGAGCACACGTCT -m ${min_length} --trimmed-only -o \
  ${expt}_trim3.fq ${expt}_raw.fq > ${expt}_trim3.cutadapt

  # 3. align against rRNA, tRNA; capture unaligned reads
  echo "... removing rRNA, tRNA"
  bowtie -v 2 -p 10 -S --un ${expt}_trim3_not_rrna_trna.fq \
  ${ref_dir}/human_rrna_trna ${expt}_trim3.fq > ${expt}_trim3_rrna_trna.sam 2> \
  ${expt}_trim3_rrna_trna.bowtiestats

  # 4. align against transcriptome
  echo "... aligning to transcriptome"
  bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/grch38.transcripts \
  ${expt}_trim3_not_rrna_trna.fq > ${expt}_trim3_footprints.sam 2> \
  ${expt}_trim3_footprints.bowtiestats

  # 5. computing multimapping weights
  echo "... calculating multimapping weights"
  rsem-calculate-expression --seed-length ${min_length} --sam ${expt}_trim3_footprints.sam ${ref_dir}/grch38.transcripts ${expt}_trim3_footprints > ${expt}_trim3_footprints.rsem.stdout 2> ${expt}_trim3_footprints.rsem.stderr

  # 6. convert .bam to .sam
  echo "... converting from bam to sam"
  samtools view -h ${expt}_trim3_footprints.transcript.bam > \
  ${expt}_trim3_footprints.transcripts.sam

  cd ..
  
done
