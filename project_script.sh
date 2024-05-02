#!/bin/bash

export LANG=en_US.UTF-8

# check for arguments
project_number=$1

if [ -d Projekt ]; then
	rm -r Projekt
fi

mkdir Projekt
mkdir -p Projekt/reads/pair_end
mkdir -p Projekt/reads/single_end
mkdir -p Projekt/ref
mkdir -p Projekt/fastqc/raw
mkdir -p Projekt/fastqc/trimmed
mkdir -p Projekt/trimmed_reads/pair_end
mkdir -p Projekt/trimmed_reads/single_end
mkdir -p Projekt/bams/pair_end
mkdir -p Projekt/bams/single_end
mkdir -p Projekt/DE_output

# create hisat2 index
echo "Creating hisat2 index..."

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz -O Projekt/ref/hg19.fa.gz
gunzip Projekt/ref/hg19.fa.gz
hisat2-build -p 2 Projekt/ref/hg19.fa hg19_index
rm Projekt/ref/hg19.fa

echo "Created hisat2 index"

# download reads and split them into SE and PE
echo "Looking for project..."

esearch -db sra -query "$project_number" \
	| efetch -format runinfo -mode xml \
	| xtract -pattern SraRunInfo -element Run \
	| awk '{for (i=1; i<=NF; i++) print $i}' > SRR_list.txt

echo "Downloading reads..."


while IFS= read -r run_number; do
	fastq-dump -X 2000 "$run_number" --split-files -O Projekt/reads
done < SRR_list.txt

echo "Splitting reads into PE and SE..."




for file in ./Projekt/reads/*_1.fastq; do
    if [ -e "$file" ]; then
        read_name=$(basename "$file" "_1.fastq")
        if [ -e "Projekt/reads/${read_name}_2.fastq" ]; then
            mv "$file" "Projekt/reads/pair_end/"
            mv "Projekt/reads/${read_name}_2.fastq" "Projekt/reads/pair_end/"
        else
            mv "$file" "Projekt/reads/single_end/"
        fi
    fi
done

echo "Reads are downloaded and split"



# fastqc before trimming
echo "Commencing quality control before trimming..."

for raw_read in Projekt/reads/pair_end/*.fastq; do
    fastqc -o Projekt/fastqc/raw "$raw_read"
done


for raw_read in Projekt/reads/single_end/*.fastq; do
    fastqc -o Projekt/fastqc/raw "$raw_read"
done

 
# trimmomatic
echo "Trimming pair-end reads"

for raw_read in Projekt/reads/pair_end/*.fastq; do
	read_name=$(basename "$raw_read" _1.fastq)
	java -jar /bioapp/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
		"Projekt/reads/pair_end/${read_name}_1.fastq" "Projekt/reads/pair_end/${read_name}_2.fastq" \
		"Projekt/trimmed_reads/pair_end/${read_name}_trimmed_1.fastq" "Projekt/trimmed_reads/pair_end/${read_name}_1_unpaired_1.fa" \
		"Projekt/trimmed_reads/pair_end/${read_name}_trimmed_2.fastq" "Projekt/trimmed_reads/pair_end/${read_name}_2_unpaired.fa" \
		ILLUMINACLIP:pe_adapter.fa:5:20:5 SLIDINGWINDOW:4:20 
done


echo "Trimming single-end reads"

for raw_read in Projekt/reads/single_end/*.fastq; do
	read_name=$(basename "$raw_read" _1.fastq)
	java -jar /bioapp/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
		"Projekt/reads/single_end/${read_name}_1.fastq" \
		"Projekt/trimmed_reads/single_end/${read_name}_trimmed_1.fastq" \
		ILLUMINACLIP:se_adapter.fa:2:20:5 SLIDINGWINDOW:4:20 
done 


# fastqc after trimming

for trimmed_read in Projekt/trimmed_reads/pair_end/*.fastq; do
    fastqc -o Projekt/fastqc/trimmed "$trimmed_read"
done

for trimmed_read in Projekt/trimmed_reads/single_end/*.fastq; do
    fastqc -o Projekt/fastqc/trimmed "$trimmed_read"
done


# hisat2 alignment and bam creation
echo "Commencing alignment for pair-end reads..."

for read in Projekt/trimmed_reads/pair_end/*_1.fastq; do
    read_name=$(basename "$read" _trimmed_1.fastq)
    echo "hisat2 alignment for $read_name started"
    hisat2 -q -x hg19_index \
        -1 "Projekt/trimmed_reads/pair_end/${read_name}_trimmed_1.fastq" \
        -2 "Projekt/trimmed_reads/pair_end/${read_name}_trimmed_2.fastq" \
        -S "Projekt/bams/pair_end/${read_name}.sam"

    samtools view -Sb -@ 2 "Projekt/bams/pair_end/${read_name}.sam" > "Projekt/bams/pair_end/${read_name}.bam"
    rm "Projekt/bams/pair_end/${read_name}.sam"
done


for read in Projekt/trimmed_reads/pair_end/*_1.fastq; do
    read_name=$(basename "$read" _trimmed_1.fastq)
    echo "hisat2 alignment for $read_name started"
    hisat2 -q -x hg19_index \
        -U "Projekt/trimmed_reads/single_end/${read_name}_trimmed_1.fastq" \
        -S "Projekt/bams/single_end/${read_name}.sam"

    samtools view -Sb -@ 2 "Projekt/bams/single_end/${read_name}.sam" > "Projekt/bams/single_end/${read_name}.bam"
    rm "Projekt/bams/single_end/${read_name}.sam"
done


# featureCounts
echo "Commencing feature count for pair-end reads..."

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz -O Projekt/ref/hg19.gtf.gz
gunzip Projekt/ref/hg19.gtf.gz

featureCounts -a Projekt/ref/hg19.gtf \
 	-g "gene_name" \
	-o "Projekt/MiSeq_counts.txt" \
    Projekt/bams/pair_end/*.bam

echo "Commencing feature count for single-end reads..."
featureCounts -a Projekt/ref/hg19.gtf \
 	-g "gene_name" \
	-o "Projekt/NextSeq_counts.txt" \
    Projekt/bams/single_end/*.bam

echo "features counted"


# DE in R
echo "Commencing DE analysis in R..."

Rscript de_analysis.R "Projekt/MiSeq_counts.txt"

Rscript de_analysis.R "Projekt/NextSeq_counts.txt"