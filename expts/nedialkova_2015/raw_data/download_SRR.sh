#!/bin/bash

grep WT_ribo_YPD SraRunTable.txt | grep -v noCHX | cut -f1 -d',' > wt.txt
grep ncs2?_ribo_YPD SraRunTable.txt | grep -v noCHX | cut -f1 -d',' > ncs2del.txt
grep WT_ribo_YPD_noCHX SraRunTable.txt | cut -f1 -d',' > wt_noCHX.txt
grep ncs2?_ribo_YPD_noCHX SraRunTable.txt | cut -f1 -d',' > ncs2del_noCHX.txt
grep ncs2?elp6?_ribo_YPD SraRunTable.txt | grep -v hbs1 | grep -v dom34 | cut -f1 -d',' > ncs2del_elp6del.txt

for expt_file in $(ls *txt | grep -v SraRunTable)
do
  expt_name=$(echo ${expt_file} | cut -f1 -d'.')
  echo "downloading" ${expt_name}
  mkdir ${expt_name}
  cd ${expt_name}
  for SRR in $(cat ../${expt_file})
  do
    echo "... downloading fastq file" ${SRR}
    fastq-dump ${SRR}
  done
  cd ..
done

