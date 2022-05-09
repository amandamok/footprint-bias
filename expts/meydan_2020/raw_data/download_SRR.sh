#!/bin/bash

# sample name in Table S3
# SRR accession from SraRunTable

mkdir monosome
cd monosome
fastq-dump SRR10302098
fastq-dump SRR10302100
cd ..

mkdir disome
cd disome
fastq-dump SRR10302102
fastq-dump SRR10302104
fastq-dump SRR10302108
