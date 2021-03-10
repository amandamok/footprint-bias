#!/bin/bash

R CMD BATCH simulate_data.R

for expt in noBias n3Bias p5Bias bothBias; do
    echo ${expt}
    cd ${expt}

    # 1. combine all fastq files
    echo "... concatenating fastq files"
    cat parts/*.fq > ${expt}.fq

    # 2. trim adaptor
    echo "... trimming adaptor"
    cutadapt -a CTGTAGGCACCATCAAT -m 18 -o ${expt}_trimmed.fq ${expt}.fq > \
        ${expt}_trimmed.cutadapt

    # 3. align to transcriptome
    echo "... aligning to transcriptome"
    bowtie -a --norc -v 2 -p 20 -S ~/iXnos/genome_data/scer.transcripts.20cds20 \
        ${expt}_trimmed.fq > ${expt}_footprints.sam \
        2> ${expt}_footprints.bowtiestats

    # 4. compute multi-mapping weights
    echo "... computing multi-mapping weights"
    rsem-calculate-expression --seed-length 18 --sam ${expt}_footprints.sam \
        ~/iXnos/genome_data/scer.transcripts.20cds20 ${expt} > \
        ${expt}.rsem.stdout 2> ${expt}.rsem.stderr

    # 5. convert from .bam to .sam
    echo "... converting to sam format"
    samtools view -h ${expt}.transcript.bam > ${expt}.transcript.sam

    cd ..
done
