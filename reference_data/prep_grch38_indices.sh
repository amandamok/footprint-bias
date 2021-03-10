### download human references (GRCh38.p13)
# gff3 of gene annotations from ENSEMBL
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gff3.gz

gunzip *gz

R CMD BATCH ../scripts/process_grch38.R

# bowtie-build grch38.ncRNA.fa grch38.ncRNA
bowtie-build grch38.transcripts.fa grch38.transcripts
rsem-prepare-reference grch38.transcripts.fa grch38.transcripts