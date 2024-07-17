# Hisat2_featureCount_RNAseq_analysis

## Data

The bulk RNA sequence data of sorted monocytes, derived from six patiens with COVID-19 and three healty donars from this study.
https://www.nature.com/articles/s41590-020-00832-x 

This will download the SRA files using parallel-fastq-dump and convert them in fastq.gz 
```
parallel-fastq-dump --sra-id SRR12926701 SRR12926702 SRR12926703 SRR12926704 SRR12926705 SRR12926706 SRR12926707 SRR12926708 SRR12926709 --threads 4 --outdir output/ --split-files --gzip
```

## Alignment

Alignemnt was conducted using Hisat2.

Reference geneome download and indexing 

```
#This will download the human reference genome and annotation file from UCSC website
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz
gunzip *.gz

#Index the ref genome
hisat2-build hg38.fa hg38_index
```

This will align the samples with indexed ref genome using hisat2 and convert the sam file to bam file using samtools.
Later sam files and bam files were removed and sorted bam files were kept for later analysis.
```
#!/bin/bash

# List of sample names
samples=("SRR12926701" "SRR12926702" "SRR12926703" "SRR12926704" "SRR12926705" "SRR12926706" "SRR12926707" "SRR12926708" "SRR12926709")


# Loop through each sample ids

for sample in "${samples[@]}"; do
  # Input file paths for forward and reverse reads
  forward_read="Rawdata/${sample}_1.fastq.gz"
  reverse_read="Rawdata/${sample}_2.fastq.gz"

  # Unzip input FASTQ files as hisat2 run on FASTQ files
  gunzip -c "$forward_read" > Rawdata/"${sample}_1.fastq"
  gunzip -c "$reverse_read" > Rawdata/"${sample}_2.fastq"

  # Output file paths for SAM, BAM, and sorted BAM files
  sam_file="${sample}.sam"
  bam_file="${sample}.bam"
  sorted_bam_file="${sample}.sorted.bam"

  # Align reads to the reference genome using HISAT2
  hisat2 -x  reference/hg38_index -1 Rawdata/"${sample}_1.fastq" -2 Rawdata/"${sample}_2.fastq" -S "$sam_file" || { echo "HISAT2 alignment failed for $sample"; exit 1; }

  # Load samtools module
  module load samtools/1.15
  # Convert SAM file to BAM file using samtools view
  samtools view -bS "$sam_file" > "$bam_file" || { echo "SAM to BAM conversion failed for $sample"; exit 1; }

  # Sort the BAM file using samtools sort
  samtools sort "$bam_file" -o "$sorted_bam_file" || { echo "BAM sorting failed for $sample"; exit 1; }

  # Create an index for the sorted BAM file using samtools index
  samtools index "$sorted_bam_file" || { echo "BAM indexing failed for $sample"; exit 1; }

  # Clean up intermediate files
  
  rm "$sam_file" "$bam_file" Rawdata/"${sample}_1.fastq" Rawdata/"${sample}_2.fastq"

done
```

## Counting reads

The featureCounts is used for counting the number of reads that align to genomic features such as genes, exons, or other genomic regions. It is a part of the Subread package, developed by the WEHI Bioinformatics team.
It takes bam files as input and Uses a gene annotation file (GTF/GFF format) that defines the genomic features to be counted (e.g., gene exons).

```
featureCounts -p -T 8 -a reference/gh38.ensGene.gtf -t exon -g gene_id -o counts_matrix.txt bam_files/*sorted.bam 
```