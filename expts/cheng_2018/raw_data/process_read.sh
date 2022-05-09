#!/bin/bash

unset PYTHONPATH

ref_dir="$HOME/footprint-bias/reference_data"
min_length=22

for timepoint_file in $(ls timepoint*txt)
do
  timepoint=$(basename ${timepoint_file} .txt)
  echo ${timepoint}
  cd ${timepoint}

  # 1. pool replicates
  echo "... pooling replicates"
  cat *fastq > ${timepoint}_raw.fq

  # 2. trim polyA from 3' end
  echo "... trimming 3' polyA"
  cutadapt -a AAAAAAAAAAAAAAAAAAAACT -m ${min_length} --trimmed-only -o \
    ${timepoint}_trim3.fq ${timepoint}_raw.fq > ${timepoint}_trim3.cutadapt

  # 3. align against rRNA, capture unaligned reads
  echo "... removing rRNA"
  bowtie -v 2 -p 10 -S --un ${timepoint}_trim3_not_rrna.fq ${ref_dir}/sk1.rrna \
    ${timepoint}_trim3.fq > ${timepoint}_trim3_rrna.sam 2> \
    ${timepoint}_trim3_rrna.bowtiestats

  # 4. align against ncRNA, capture unaligned reads
  echo "... removing ncRNA"
  bowtie -v 2 -p 10 -S --un ${timepoint}_trim3_not_rrna_ncrna.fq \
    ${ref_dir}/rna_coding ${timepoint}_trim3_not_rrna.fq > \
    ${timepoint}_trim3_ncrna.sam 2> ${timepoint}_trim3_ncrna.bowtiestats

  # 5. align against transcriptome
  echo "... aligning to transcriptome"
  bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/sk1.transcripts.20cds20 \
    ${timepoint}_trim3_not_rrna_ncrna.fq > \
    ${timepoint}_trim3_footprints.sam 2> \
    ${timepoint}_trim3_footprints.bowtiestats

  # 6. calculate multi-mapping weights
  echo "... calculating multimapping weights"
  rsem-calculate-expression --seed-length ${min_length} --sam \
    ${timepoint}_trim3_footprints.sam ${ref_dir}/sk1.transcripts.20cds20 \
    ${timepoint}_trim3_footprints > \
    ${timepoint}_trim3_footprints.rsem.stdout 2> \
    ${timepoint}_trim3_footprints.rsem.stderr

  # 7. convert .bam to .sam
  samtools view -h ${timepoint}_trim3_footprints.transcript.bam > \
	  ${timepoint}_trim3_footprints.transcripts.sam

  cd ..
done
