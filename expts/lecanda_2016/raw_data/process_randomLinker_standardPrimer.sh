#!/bin/bash

unset PYTHONPATH

ref_dir="$HOME/footprint-bias/reference_data"
min_length=22
expt="randomLinker_standardPrimer"

cd ${expt}

# 1. pool replicates
echo "... pooling replicates"
cat SRR*fastq > ${expt}_raw.fq

# 2. trim 3' adapter
echo "... trimming 3' adapter"
cutadapt -a CTGTAGGCACCATCAAT -m ${min_length} --trimmed-only -o \
${expt}_trim3.fq ${expt}_raw.fq > ${expt}_trim3.cutadapt

# 3. trim random linker
echo "... trimming random linker"
umi_tools extract --extract-method=regex --bc-pattern="(.+)(?P<umi_1>.{4})$" \
-I ${expt}_trim3.fq -S ${expt}_trim3_trimUMI.fq -L ${expt}_trim3_trimUMI.umitools

# 3. remove contaminant reads: rRNA
echo "... removing rRNA"
bowtie -v 2 -p 10 -S --un ${expt}_trim3_not_rrna.fq ${ref_dir}/ScerRRNA \
${expt}_trim3_trimUMI.fq > ${expt}_trim3_rrna.sam 2> \
${expt}_trim3_rrna.bowtiestats

# 4. remove contaminant reads: ncRNA
echo "... removing ncRNA"
bowtie -v 2 -p 10 -S --un ${expt}_trim3_not_rrna_ncrna.fq \
${ref_dir}/rna_coding ${expt}_trim3_not_rrna.fq > \
${expt}_trim3_ncrna.sam 2> ${expt}_trim3_ncrna.bowtiestats

# 5. align to transcriptome 
echo "... aligning to transcriptome"
bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/scer.transcripts.20cds20 \
${expt}_trim3_not_rrna_ncrna.fq > ${expt}_trim3_footprints.sam 2> \
${expt}_trim3_footprints.bowtiestats

# 6. compute multimapping weights
echo "... compute multimapping weights"
rsem-calculate-expression --sam ${expt}_trim3_footprints.sam \
${ref_dir}/scer.transcripts.20cds20 ${expt}_trim3_footprints > \
${expt}_trim3_footprints.rsem.stdout 2> ${expt}_trim3_footprints.rsem.stderr
