#!/bin/bash

for x in {1..4}
do
  expt=GSM161159$x
  echo "downloading " ${expt} " files"
  mkdir ${expt}
  cd ${expt}
  for SRR in $(grep ${expt} ../SraRunTable.txt | cut -f1 -d',')
  do
    echo "... downloading fastq file" ${SRR}
    fastq-dump ${SRR}
  done
  cd ..
done

mv GSM1611591 SGop
mv GSM1611592 GCop
mv GSM1611593 SGkit
mv GSM1611594 GCkit
