#!/bin/bash

ref_dir="/mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised"

# generate genome sequence from bowtie2 index
bowtie2-inspect ${ref_dir}/sk1_revised_100331_chr_to_use > sk1_revised_100331_chr_to_use.fa

# split genome fasta into individual chromosomes 
Rscript split_genome_by_chr.R

# generate BSgenome package for SK1
if [ ! -f BSgenome.ScerSK1.brarLab.ScerSK1_0.0.1.tar.gz ]
then
  Rscript -e 'BSgenome::forgeBSgenomeDataPkg("BSgenome.ScerSK1.brarLab.ScerSK1-seed")'
  R CMD build BSgenome.ScerSK1.brarLab.ScerSK1 --no-manual --no-build-vignettes
  R CMD check BSgenome.ScerSK1.brarLab.ScerSK1_0.0.1.tar.gz --no-vignettes --no-build-vignettes --no-manual
  R CMD INSTALL BSgenome.ScerSK1.brarLab.ScerSK1_0.0.1.tar.gz --no-docs
else 
  echo "Package already installed"
fi
