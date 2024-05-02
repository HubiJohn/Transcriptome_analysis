This project has been done as a part of the _Transcriptome analysis_ (pl. Analiza transkryptomu) course during my studies in Wroclaw University of Environmental and Life Sciences. 
The aim was to construct a semi-automated pipeline that would requeire minimal command line input from the user.

The RNA-Seq analysis pipeline consists of the following steps:
    - Quality control (FASTQC)
    - Sequence data cleaning if required (Trimmomatic)
    - Alignment; mapping to a reference genome (hisat2, STAR) usually preceded by ref genome indexing
    - Counting reads for each individual gene (featureCounts)
    - Differential expression analysis, including normalization and p-value adjustment (DESeq2 in R)

This project.
The executable `project_script.sh` takes one argument, the PRJNA bioproject number. Then, the following is done:
    - download all associated fastq files using `fastq-dump`.
    - split them into pair-end and single-end files in order to apply the appriopriate parameters in further steps
    - conduct RNA-Seq analysis steps
    - run an R script to create graphical interpretations (heatmaps, PCA) 

The PL_project_summary.pdf file contains highlights of the code used in the project 
and the results obtained using the bioproject PRJNA313294.
This bioproject consists of 4 single-end SRR files and 4 pair-end SRR files.
Each set has Zika-infected and control cells transcriptomes.
