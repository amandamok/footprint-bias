#!/bin/bash


# gb-june2016-fp_tc-tubeA_idx15_combined

# decompress file, clip the polyA

zcat /mnt/brarlab/zecheng/Data_160718/gb-june2016-fp_tc-tubeA_idx15_combined.fastq.gz | fastx_clipper -Q33 -a AAAAAAAAAAAAAAA -l22 -c -n -v 2>Statistics/gb-june2016-fp_tc-tubeA_idx15_combined_clip.txt > gb-june2016-fp_tc-tubeA_idx15_combined_clipped.fq

# align against rRNA reference file, capture the unaligned reads

bowtie2 -p36 --solexa-quals --un gb-june2016-fp_tc-tubeA_idx15_combined_norrna.fq -x /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_rrna -U gb-june2016-fp_tc-tubeA_idx15_combined_clipped.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx15_combined_rrna.txt > /dev/null

# align against yeast genomic reference file

tophat -p36 --solexa-quals --GTF /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.gtf --no-novel-juncs --output-dir gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_chr_to_use gb-june2016-fp_tc-tubeA_idx15_combined_norrna.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx15_combined_genome.txt

# exclude multiple alignment reads - optional

samtools view -b -q50 gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome/accepted_hits.bam > gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome_unique.bam

# preform quality control, find A site offsets

samtools index gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome_unique.bam

fp-framing -o Statistics/gb-june2016-fp_tc-tubeA_idx15_combined_unique -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome_unique.bam


# count reads for each gene, generate qexpr.txt

fp-count -o gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome_unique_qexpr.txt -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome_unique.bam


# make wiggle-track files

wiggle-track -o gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome_unique.wig -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx15_combined_vs_genome_unique.bam



# gb-june2016-fp_tc-tubeA_idx14_combined

# decompress file, clip the polyA

zcat /mnt/brarlab/zecheng/Data_160718/gb-june2016-fp_tc-tubeA_idx14_combined.fastq.gz | fastx_clipper -Q33 -a AAAAAAAAAAAAAAA -l22 -c -n -v 2>Statistics/gb-june2016-fp_tc-tubeA_idx14_combined_clip.txt > gb-june2016-fp_tc-tubeA_idx14_combined_clipped.fq

# align against rRNA reference file, capture the unaligned reads

bowtie2 -p36 --solexa-quals --un gb-june2016-fp_tc-tubeA_idx14_combined_norrna.fq -x /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_rrna -U gb-june2016-fp_tc-tubeA_idx14_combined_clipped.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx14_combined_rrna.txt > /dev/null

# align against yeast genomic reference file

tophat -p36 --solexa-quals --GTF /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.gtf --no-novel-juncs --output-dir gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_chr_to_use gb-june2016-fp_tc-tubeA_idx14_combined_norrna.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx14_combined_genome.txt

# exclude multiple alignment reads - optional

samtools view -b -q50 gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome/accepted_hits.bam > gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome_unique.bam

# preform quality control, find A site offsets

samtools index gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome_unique.bam

fp-framing -o Statistics/gb-june2016-fp_tc-tubeA_idx14_combined_unique -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome_unique.bam


# count reads for each gene, generate qexpr.txt

fp-count -o gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome_unique_qexpr.txt -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome_unique.bam


# make wiggle-track files

wiggle-track -o gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome_unique.wig -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx14_combined_vs_genome_unique.bam






# gb-june2016-fp_tc-tubeA_idx12_combined

# decompress file, clip the polyA

zcat /mnt/brarlab/zecheng/Data_160718/gb-june2016-fp_tc-tubeA_idx12_combined.fastq.gz | fastx_clipper -Q33 -a AAAAAAAAAAAAAAA -l22 -c -n -v 2>Statistics/gb-june2016-fp_tc-tubeA_idx12_combined_clip.txt > gb-june2016-fp_tc-tubeA_idx12_combined_clipped.fq

# align against rRNA reference file, capture the unaligned reads

bowtie2 -p36 --solexa-quals --un gb-june2016-fp_tc-tubeA_idx12_combined_norrna.fq -x /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_rrna -U gb-june2016-fp_tc-tubeA_idx12_combined_clipped.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx12_combined_rrna.txt > /dev/null

# align against yeast genomic reference file

tophat -p36 --solexa-quals --GTF /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.gtf --no-novel-juncs --output-dir gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_chr_to_use gb-june2016-fp_tc-tubeA_idx12_combined_norrna.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx12_combined_genome.txt

# exclude multiple alignment reads - optional

samtools view -b -q50 gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome/accepted_hits.bam > gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome_unique.bam

# preform quality control, find A site offsets

samtools index gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome_unique.bam

fp-framing -o Statistics/gb-june2016-fp_tc-tubeA_idx12_combined_unique -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome_unique.bam


# count reads for each gene, generate qexpr.txt

fp-count -o gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome_unique_qexpr.txt -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome_unique.bam


# make wiggle-track files

wiggle-track -o gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome_unique.wig -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx12_combined_vs_genome_unique.bam







# gb-june2016-fp_tc-tubeA_idx11_combined

# decompress file, clip the polyA

zcat /mnt/brarlab/zecheng/Data_160718/gb-june2016-fp_tc-tubeA_idx11_combined.fastq.gz | fastx_clipper -Q33 -a AAAAAAAAAAAAAAA -l22 -c -n -v 2>Statistics/gb-june2016-fp_tc-tubeA_idx11_combined_clip.txt > gb-june2016-fp_tc-tubeA_idx11_combined_clipped.fq

# align against rRNA reference file, capture the unaligned reads

bowtie2 -p36 --solexa-quals --un gb-june2016-fp_tc-tubeA_idx11_combined_norrna.fq -x /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_rrna -U gb-june2016-fp_tc-tubeA_idx11_combined_clipped.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx11_combined_rrna.txt > /dev/null

# align against yeast genomic reference file

tophat -p36 --solexa-quals --GTF /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.gtf --no-novel-juncs --output-dir gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_chr_to_use gb-june2016-fp_tc-tubeA_idx11_combined_norrna.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx11_combined_genome.txt

# exclude multiple alignment reads - optional

samtools view -b -q50 gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome/accepted_hits.bam > gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome_unique.bam

# preform quality control, find A site offsets

samtools index gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome_unique.bam

fp-framing -o Statistics/gb-june2016-fp_tc-tubeA_idx11_combined_unique -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome_unique.bam


# count reads for each gene, generate qexpr.txt

fp-count -o gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome_unique_qexpr.txt -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome_unique.bam


# make wiggle-track files

wiggle-track -o gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome_unique.wig -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx11_combined_vs_genome_unique.bam






# gb-june2016-fp_tc-tubeA_idx10_combined

# decompress file, clip the polyA

zcat /mnt/brarlab/zecheng/Data_160718/gb-june2016-fp_tc-tubeA_idx10_combined.fastq.gz | fastx_clipper -Q33 -a AAAAAAAAAAAAAAA -l22 -c -n -v 2>Statistics/gb-june2016-fp_tc-tubeA_idx10_combined_clip.txt > gb-june2016-fp_tc-tubeA_idx10_combined_clipped.fq

# align against rRNA reference file, capture the unaligned reads

bowtie2 -p36 --solexa-quals --un gb-june2016-fp_tc-tubeA_idx10_combined_norrna.fq -x /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_rrna -U gb-june2016-fp_tc-tubeA_idx10_combined_clipped.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx10_combined_rrna.txt > /dev/null

# align against yeast genomic reference file

tophat -p36 --solexa-quals --GTF /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.gtf --no-novel-juncs --output-dir gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_chr_to_use gb-june2016-fp_tc-tubeA_idx10_combined_norrna.fq 2>Statistics/gb-june2016-fp_tc-tubeA_idx10_combined_genome.txt

# exclude multiple alignment reads - optional

samtools view -b -q50 gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome/accepted_hits.bam > gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome_unique.bam

# preform quality control, find A site offsets

samtools index gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome_unique.bam

fp-framing -o Statistics/gb-june2016-fp_tc-tubeA_idx10_combined_unique -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome_unique.bam


# count reads for each gene, generate qexpr.txt

fp-count -o gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome_unique_qexpr.txt -b /mnt/brarlab/gloriabrar/Libraries/SangerSGRP/SK1Revised/bowtie2_sk1_revised/sk1_revised_100331_cleanest.bed -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome_unique.bam


# make wiggle-track files

wiggle-track -o gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome_unique.wig -a Statistics/asites_fp.txt gb-june2016-fp_tc-tubeA_idx10_combined_vs_genome_unique.bam


