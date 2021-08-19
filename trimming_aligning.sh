#!/bin/bash

#Set up conda environment
conda create -y -n chikv_amplicon
#Activate 
conda activate chikv_amplicon

#Install necessary packages
conda install cutadapt
conda install ivar

#Unzip files
gunzip *.fastq.gz

#Trim forward and reverse primers from paired reads 
#Test on one file
#cutadapt -g ^GCCATCATTAAATATGCAGCCAGC -G ^GGTGTGTCTCTTAGGGGACACATATACC -o trimmed_R1.fastq -p trimmed_R2.fastq DFMO1_S3_L001_R1_001.fastq DFMO1_S3_L001_R2_001.fastq

#Loop to trim primers from all paired read files 
for file in *R1_001.fastq;
do
    prefix=${file%_R1_001.fastq} #Get file prefix
    cutadapt -g ^GCCATCATTAAATATGCAGCCAGC -G ^GGTGTGTCTCTTAGGGGACACATATACC -o ${prefix}_R1.fastq -p ${prefix}_R2.fastq ${prefix}_R1_001.fastq ${prefix}_R2_001.fastq > ${prefix}_cutadapt_report.txt 
done

## bowtie2 was used to index the reference sequence 
bowtie2-build chikv_gene_reference.fasta chikv_gene_ref

## reads were mapped to the reference sequence with bowtie 2

for file in *_R1.fastq;
do
    prefix=${file%_R1.fastq} #Get file prefix
    bowtie2 -x chikv_gene_ref -1 ${prefix}_R1.fastq -2 ${prefix}_R2.fastq -S ${prefix}.mapped.sam > ${prefix}_bowtie_report.txt
done

#Convert *sam to *bam
for i in *.sam; 
do bn=$(basename $i .sam);
samtools view -b ${bn}.sam -o ${bn}.bam
done

#sort and index bam
for i in *.bam;  
do bn=$(basename $i .bam); 
samtools sort -o ${bn}.sorted.bam ${bn}.bam && samtools index ${bn}.sorted.bam; 
done

#Create consensus using mpileup and ivar
for i in *.sorted.bam;  
do bn=$(basename $i .sorted.bam); 
samtools mpileup -aa -A -d 0 -Q 20 ${bn}.sorted.bam | ivar consensus -t 0 -m 20 -q 20 -p "${bn}.consensus"; 
done

#Create a TSV with the SNPs
for i in *.sorted.bam;  
do bn=$(basename $i .sorted.bam); 
samtools mpileup -aa -A -d 0 -B -Q 20 --reference chikv_gene_reference.fasta ${bn}.sorted.bam | ivar variants -t 0 -m 20 -p "${bn}.variants" -r chikv_gene_reference.fasta; 
done