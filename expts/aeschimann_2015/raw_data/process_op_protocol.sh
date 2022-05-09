#!/bin/bash

ref_dir="$HOME/footprint-bias/reference_data"
min_length=22

unset PYTHONPATH

# adapter sequences from UC Davis
# https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina-adapter-sequences_1000000002694-00.pdf

for expt in $(ls -d *op)
do
  echo ${expt}
  cd ${expt}

  # 1. pool replicates
  echo "... pooling replicates"
  cat SRR*fastq > ${expt}_raw.fq

  # 2. trim 3' adapter
  echo "... trimming adapters"
  cutadapt -a TGGAATTCTCGGGTGCCAAGG -m ${min_length} \
  --trimmed-only -o ${expt}_trim3.fq ${expt}_raw.fq > ${expt}_trim3.cutadapt

  # 3. align against RNA contaminants; capture unaligned reads
  echo "... removing RNA contaminants"
  bowtie -v 2 -p 10 -S --un ${expt}_trim3_not_contaminants.fq \
  ${ref_dir}/WBcel235.contaminants ${expt}_trim3.fq > \
  ${expt}_trim3_contaminants.sam 2> ${expt}_trim3_contaminants.bowtiestats

  # 4. align against transcriptome
  echo "... aligning to transcriptome"
  bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/WBcel235.transcripts \
  ${expt}_trim3_not_contaminants.fq > ${expt}_trim3_footprints.sam 2> \
  ${expt}_trim3_footprints.bowtiestats

  # 5. computing multimapping weights
  echo "... calculating multimapping weights"
  rsem-calculate-expression --seed-length ${min_length} --sam \
  ${expt}_trim3_footprints.sam ${ref_dir}/WBcel235.transcripts \
  ${expt}_trim3_footprints > ${expt}_trim3_footprints.rsem.stdout 2> \
  ${expt}_trim3_footprints.rsem.stderr

  # 6. convert .bam to .sam
  echo "... converting from bam to sam"
  samtools view -h ${expt}_trim3_footprints.transcript.bam > \
  ${expt}_trim3_footprints.transcripts.sam

  cd ..
done

