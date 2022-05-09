# download and process reads

`GSE108778_series_matrix.txt`
downloaded from GEO
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108778

`SraRunTable.txt`
downloaded from SRA Run Selector (Total Metadata)
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA428526&o=acc_s%3Aa

1. identify SRR accession numbers per timepoint: `parse_SRR.R`

2. download FASTQ files: `download_SRR.sh`

3. process reads: `process_read.sh`
- trim polyA tail (`cutadapt`)
- filter out rRNA reads (`bowtie`)
- filter ncRNA reads (`bowtie`)
- align to transcriptome (`bowtie`)
- calcualte multi-mapping weights (`RSEM`)
