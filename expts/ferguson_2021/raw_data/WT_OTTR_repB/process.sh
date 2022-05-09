#!/bin/bash

# cutadapt version 3.4
# _tools version 1.0.0
## note: unset PYTHONPATH before running cutadapt 

ref_dir="$HOME/footprint-bias/reference_data"
unset PYTHONPATH

# 1. trim constant adaptor
cat R1b_round1.fastq.gz R1b_round2.fastq.gz | zcat | \
    cutadapt -j 10 -m 20 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --untrimmed-output ottr_untrimmed.fq -o ottr_trimmed_adaptor.fq - > \ 
    ottr_trimmed_adaptor.cutadapt

# 2. trim 1nt from 5' end --> NTA
cutadapt -j 10 -m 20 -u 1 --rename='{id} NTA={cut_prefix}' -o \
    ottr_trimmed_adaptor_nta.fq ottr_trimmed_adaptor.fq > 
    ottr_trimmed_adaptor_nta.cutadapt
## A: 5 398 707 (7.1%)
## C: 60 347 729 (79.3%)
## G: 4 116 198 (5.4%)
## T: 6 196 084 (8.1%)
## N: 24 555 (0.0%)

# 3. trim 1nt from 3' end --> TPRT
cutadapt -j 10 -m 20 -u -1 --rename='{id} {comment} TPRT={cut_suffix}' -o \
    ottr_trimmed_adaptor_nta_tprt.fq ottr_trimmed_adaptor_nta.fq > \
    ottr_trimmed_adaptor_nta_tprt.cutadapt
## A: 54 929 083 (75.2%)
## C: 10 043 (0.0%)
## G: 18 122 951 (24.8%)
## T: 17 473 (0.0%)
## N: 233 (0.0%)

# 4. remove contaminant reads: rRNA
bowtie -v 2 -p 10 -S --un ottr_trimmed_not_rrna.fq \
    $ref_dir/ScerRRNA ottr_trimmed_adaptor_nta_tprt.fq > \
    ottr_trimmed_rrna.sam 2> ottr_trimmed_rrna.bowtiestats

# 5. remove contaminant reads: ncRNA
bowtie -v 2 -p 10 -S --un ottr_trimmed_not_rrna_ncrna.fq \
    $ref_dir/rna_coding ottr_trimmed_not_rrna.fq > \
    ottr_trimmed_ncrna.sam 2> ottr_trimmed_ncrna.bowtiestats

# 6. align to transcriptome
bowtie -v 2 -p 10 -S --norc -a $ref_dir/scer.transcripts.20cds20 \ 
    ottr_trimmed_not_rrna_ncrna.fq > ottr_trimmed_footprints.sam 2> \ 
    ottr_trimmed_footprints.bowtiestats
## 10 450 600 reads
## 17 904 807 alignments

# 7. compute multimapping weights
rsem-calculate-expression --seed-length 20 --sam \ 
    ottr_trimmed_footprints.sam \
    $ref_dir/scer.transcripts.20cds20 \ 
    ottr_trimmed_footprints > \
    ottr_trimmed_footprints.rsem.stdout 2> \
    ottr_trimmed_footprints.rsem.stderr
