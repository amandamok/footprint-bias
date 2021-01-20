#!/bin/bash

for expt in $(cat adaptor_1.txt)
do
    echo ${expt}

    # 1. create directory for dataset
    mkdir ${expt}
    cd ${expt}

    # 2. download .fastq file
    echo "... downloading fastq file"
    fastq-dump ${expt}

    # 3. trim 3' adaptor
    echo "... trimming 3' adaptor"
    cutadapt -a CTGTAGGCACCATCAAT -m 15 -o ${expt}_trim3.fq ${expt}.fastq > \
        ${expt}_trim3.cutadapt

    # 4. trim 5' adaptor
    echo "... triming 5' adaptor"
    cutadapt -g ^NNNY --match-read-wildcards -m 15 -o ${expt}_trim3_trim5.fq \
        ${expt}_trim3.fq > ${expt}_trim3_trim5.cutadapt

    # 5. remove rRNA
    echo "... removing rRNA"
    bowtie -v 2 -p 20 -S --un ${expt}_not_rRNA.fq ~/iXnos/genome_data/ScerRRNA \
        ${expt}_trim3_trim5.fq > ${expt}_rRNA.sam 2> ${expt}_rRNA.bowtiestats

    # 6. remove ncRNA
    echo "... removing ncRNA"
    bowtie -v 2 -p 20 -S --un ${expt}_not_ncRNA.fq ~/iXnos/genome_data/rna_coding \
        ${expt}_not_rRNA.fq > ${expt}_ncRNA.sam 2> ${expt}_ncRNA.bowtiestats

    # 7. align to transcriptome
    echo "... aligning to transcriptome"
    bowtie -a --norc -v 2 -p 20 -S ~/iXnos/genome_data/scer.transcripts.20cds20 \
        ${expt}_not_ncRNA.fq > ${expt}_footprints.sam 2> ${expt}_footprints.bowtiestats

    # 8. compute multi-mapping weights
    echo "... computing multi-mapping weights"
    rsem-calculate-expression --seed-length 18 --sam ${expt}_footprints.sam \
        ~/iXnos/genome_data/scer.transcripts.20cds20 ${expt} > ${expt}.rsem.stdout \
        2> ${expt}.rsem.stderr

    # 9. convert from .bam to .sam
    echo "... converting to .sam"
    samtools view -h ${expt}.transcript.bam > ${expt}.transcript.sam

    cd ..
done

for expt in $(cat adaptor_2.txt)
do
    echo ${expt}

    # 1. create directory for dataset
    mkdir ${expt}
    cd ${expt}

    # 2. download .fastq file
    echo "... downloading fastq file"
    fastq-dump ${expt}

    # 3. trim 3' adaptor
    echo "... trimming 3' adaptor"
    cutadapt -a NNNNNNCACTCGGGCACCAAGGA --match-read-wildcards -m 15 -o \
        ${expt}_trim3.fq ${expt}.fastq > ${expt}_trim3.cutadapt

    # 4. trim 5' adaptor
    echo "... triming 5' adaptor"
    cutadapt -g ^NNNY --match-read-wildcards -m 15 -o ${expt}_trim3_trim5.fq \
        ${expt}_trim3.fq > ${expt}_trim3_trim5.cutadapt

    # 5. remove rRNA
    echo "... removing rRNA"
    bowtie -v 2 -p 20 -S --un ${expt}_not_rRNA.fq ~/iXnos/genome_data/ScerRRNA \
        ${expt}_trim3_trim5.fq > ${expt}_rRNA.sam 2> ${expt}_rRNA.bowtiestats

    # 6. remove ncRNA
    echo "... removing ncRNA"
    bowtie -v 2 -p 20 -S --un ${expt}_not_ncRNA.fq ~/iXnos/genome_data/rna_coding \
        ${expt}_not_rRNA.fq > ${expt}_ncRNA.sam 2> ${expt}_ncRNA.bowtiestats

    # 7. align to transcriptome
    echo "... aligning to transcriptome"
    bowtie -a --norc -v 2 -p 20 -S ~/iXnos/genome_data/scer.transcripts.20cds20 \
        ${expt}_not_ncRNA.fq > ${expt}_footprints.sam 2> ${expt}_footprints.bowtiestats

    # 8. compute multi-mapping weights
    echo "... computing multi-mapping weights"
    rsem-calculate-expression --seed-length 18 --sam ${expt}_footprints.sam \
        ~/iXnos/genome_data/scer.transcripts.20cds20 ${expt} > ${expt}.rsem.stdout \
        2> ${expt}.rsem.stderr

    # 9. convert from .bam to .sam
    echo "... converting to .sam"
    samtools view -h ${expt}.transcript.bam > ${expt}.transcript.sam

    cd ..
done
