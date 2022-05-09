#!/bin/bash

grep ribosome SraRunTable.txt | grep 'wild type' | cut -f1 -d',' > wt_srr
grep ribosome SraRunTable.txt | grep 'HelaS3-SNORD45C-KO' | cut -f1 -d',' > ko_srr

echo "downloading WT files"
mkdir WT_RPF
cd WT_RPF
for SRR in $(cat ../wt_srr)
do
  echo "... downloading fastq file" ${SRR}
  fastq-dump ${SRR}
done
cd ..

echo "downloading KO files"
mkdir KO_RPF
cd KO_RPF
for SRR in $(cat ../ko_srr)
do
  echo "... downloading fastq file" ${SRR}
  fastq-dump ${SRR}
done
