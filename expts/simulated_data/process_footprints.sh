#!/bin/bash

for expt in noBias n3Bias p5Bias bothBias; do
    echo ${expt}
    cd ${expt}
    
    # 1. change file extension from fasta to fastq
    for file in *.fa; do
        mv -- "$file" "${file%.fa}.fq"
    done
    
    # 2. combine all fastq files
    cat *.fq > ${expt}.fq
    
    # 3. trim adaptor
    cutadapt -a CTGTAGGCACCATCAAT -m 18 -o ${expt}_trimmed.fq ${expt}.fq > \ 
        ${expt}_trimmed.cutadapt
        
    # 4. align to transcriptome
    bowtie -a --norc -v 2 -p 20 -S ~/iXnos/genome_data/scer.transcripts.20cds20 \
        ${expt}_trimmed.fq > ${expt}_footprints.sam \ 
        2> ${expt}_footprints.bowtiestats
        
    # 5. compute multi-mapping weights
    rsem-calculate-expression --seed-length 18 --sam ${expt}_footprints.sam \
        ~ /iXnos/genome_data/scer.transcripts.20cds20 ${expt} > \ 
        ${expt}.rsem.stdout 2> ${expt}.rsem.stderr
        
    # 6. convert from .bam to .sam
    samtools view -h ${expt}.transcript.bam > ${expt}.transcript.sam

    cd ..
done
