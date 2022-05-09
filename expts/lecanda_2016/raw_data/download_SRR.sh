#!/bin/bash

for expt in $(cut -f2 GEO_table.txt | grep -v mRNA | grep -v DualLigation | cut -f1 -d'_' | sort | uniq) 
do
  echo "downloading" ${expt}
  mkdir ${expt}
  cd ${expt}
  for GSM in $(grep ${expt} ../GEO_table.txt | grep -v mRNA | cut -f1)
  do
    SRR=$(grep ${GSM} ../SraRunTable.txt | cut -f1 -d',')
    echo "... downloading fastq file" ${SRR}
    fastq-dump ${SRR}
  done
  cd ..
done

mv 4+3N randomLinker_randomPrimer
mv 4N randomLinker_standardPrimer
mv Nonrandom fixedLinker_standardPrimer
