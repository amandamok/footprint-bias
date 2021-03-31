#!/bin/bash
ref_dir="../../../../reference_data"
min_length="15"

# download reads and trim adaptor
for expt in $(cut -f1 -d"," sra_metadata.csv | grep SRR)
do
    echo ${expt}

    # 1 create directory
    mkdir ${expt}
    cd ${expt}

    # 2. download .fastq file
    echo "... downloading fastq file"
    fasterq-dump ${expt}

    # 3. trim 3' adaptor
    echo "... trimming 3' adaptor"
    cutadapt -a AGATCGGAAGAGCACACGTCT -m ${min_length} \
        -o ${expt}_trim3.fq ${expt}.fastq > ${expt}_trim3.cutadapt

     # 4. remove rRNA
    echo "... removing rRNA"
    bowtie -v 2 -p 20 -S --un ${expt}_not_ncRNA.fq ${ref_dir}/grch38.ncRNA \
        ${expt}_trim3.fq > ${expt}_ncRNA.sam 2> ${expt}_ncRNA.bowtiestats

    # 5. align to transcriptome
    echo "... aligning to transcriptome"
    bowtie -a --norc -v 2 -p 20 -S ${ref_dir}/grch38.transcripts \
        ${expt}_not_ncRNA.fq > ${expt}_footprints.sam 2> ${expt}_footprints.bowtiestats

    # 6. compute multi-mapping weights
    echo "... computing multi-mapping weights"
    rsem-calculate-expression --seed-length ${min_length} -p 20 --alignments \
        ${expt}_footprints.sam ${ref_dir}/grch38.transcripts ${expt} > \
        ${expt}.rsem.stdout 2> ${expt}.rsem.stderr

    cd ..
done
