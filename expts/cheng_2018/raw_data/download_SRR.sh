#!/bin/bash

for timepoint_file in $(ls timepoint*txt)
do 
  timepoint=$(basename ${timepoint_file} .txt)
  echo ${timepoint}

  # create directory for datasets
  mkdir ${timepoint} 

  # download .fastq files
  for SRR in $(cat ${timepoint_file})
  do
    echo "... downloading fastq file" ${SRR}
    cd ${timepoint}
    fastq-dump ${SRR}
    cd ..
  done
done

