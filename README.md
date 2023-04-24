#  Identification of targets for the microRNA miR-4728-3p in ribosome and polysome sequencing data

## Data

The data for this project consists of ribosome and polysome profiling data for SKBR3 (HER2-positive breast cancer cell line) treated with antisense oligonucleotide (ASO) to block the miRNA or a control oligonucleotide, with matched total RNA-seq. There are three replicates per condition.

## Analysis environment

A conda environment was created for this analysis and the environment file is attached.

- samtools (v)
- gfftools (v)
- htseq (v)
- fastqc (v)

The following tools were not available through Conda and were installed manually following the installation steps provided in their respective documentations:

- HISAT2 (v2.2.1)
- Novoalign (v)

The following R packages are required:
- anota2seq
- clusterProfiler
- annotationdbi


## 0. Quality control

Both datasets were tested using FastQC to ensure that all sample files were suitable for use in downstream analyses.

### 0.1. Polysome profiling data

```shell
cat files_list.txt | while read file; do fastqc -o 0_fastqc --noextract $file
```

### 0.2. Ribosome profiling data

```shell
cat files_list.txt | while read file; do fastqc -o 0_fastqc --noextract $file
```

## 1. Alignment of reads

### 1.1. Polysome profiling data

#### 1.1.1. Building index

The GRCh38.p13 assembly of the human genome from the Genome Reference Consortium and the Gencode 43 release of the human transcriptome annotation (Ensembl release 109) was used to build a reference index for the alignment, following the [protocol](http://daehwankimlab.github.io/hisat2/howto/) provided in the HISAT2 documentation.

Unfortunately, the indexing of a human genome with transcript annotations required over 160GB of memory and could not be completed. Therefore, a prebuilt index used previously for a ScanB project was acquired and used. The index was created using hg38 assembly, gencode v41 annotation for the transcriptome annotation, and the dbSNP build 155 for the SNP data.

#### 1.1.2. Alignment

The reads from the polysome profiling were aligned to the reference index, using [HISAT2](http://daehwankimlab.github.io/hisat2/manual/). The alignment output was directly piped into samtools to sort and convert them from SAM to BAM formats.

```shell
cat files_list.txt | while read line;do file_dir=$(echo ../data/$line); file_name=$(echo $line | cut -d "/" -f 2); file_folder=$(echo $line | cut -d "/" -f 1); mate_1=$(echo ${file_dir}_R1_001.fastq.gz); mate_2=$(echo ${file_dir}_R2_001.fastq.gz); echo $file_name; hisat2 -p 15 -q --fr --new-summary --summary-file 1_alignments/summary/$file_name.sam.summary --dta --rna-strandness RF --non-deterministic --max-intronlen 2000000 -x ../reference/hisat2/genome_snp_tran -1 $mate_1 -2 $mate2 | samtools sort -o 1_alignments/$file_folder/$file_name.bam; done;
```

#### 1.1.3. Summary statistics from alignment

A list of summary files and their directories was generated using bash.

```shell
ls summary_files | while read file; do echo summary_files/$file; done > sample_list.txt
```

A custom Python script was used to parse the summary files with the list and visualized using R. The two sets of alignments, with the prebuilt and scanb index, were compared. No significant differences were identified, therefore the alignment with the most up-to-date reference index was used for the subsequent analyses.

### 1.2. Ribosome profiling data

#### 1.2.1. Building the index

While Novoalign includes the protocol for generating an index with known transcripts, the protocol was extremely long and poorly documented. It was decided that the benefits of attempting to create an index with transcriptome annotation was not great enough to outweigh the time and resources it would require. Therefore, an index was created with only the genome assembly. The same GRCh38 genome assembly used for the HISAT2 index was used.

```shell
# Index is created
novoindex ~/reference/novoalign/GRCh38_no_alt_maskedGRC.nix ~/reference/raw/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta

# A soft link to the index was created in the directory where novoalign will be run
ln -s ~/reference/novoalign/GRCh38_no_alt_maskedGRC.nix ~/ribosome/index_link.nix
```

In addition to the reference genome, another index was created for the complete human ribosomal DNA repeating units (U13369.1).

```shell
# Index is created
novoindex ~/reference/novoalign/rRNA.nix ~/reference/raw/human_complete_rRNA.fasta

# A soft link to the index is created
ln -s ~/reference/novoalign/rRNA.nix ~/ribosome/rRNA_link.nix
```

#### 1.2.2. Alignment with Novoalign

Novoalign was run using the index file created above. The alignment parameters were set to be extra sensitive compared to regular RNA seq alignments, due to the short read lengths of the RPFs. Novoalign also allows the clipping of adapter sequences that may be present in the 3' ends of the reads. Therefore, the adapter sequence used for the Ribo-seq procedure was provided to the aligner.

```shell
ls data | while read file; do sample_name=$(echo $file | sed "s/.fastq.gz//"); echo $file; novoalign -c 15 -d index_link.nix -f data/$file -F STDFQ -a AGATCGGAAGAGCACACGTCT -l 17 -h -1 -1 -t 90 -g 50 -x 15 -o SAM -o FullNW -r All 51 -e 51 2> 1_alignments/novoalign/genome/$file.summary | samtools sort -o 1_alignments/novoalign/genome/$sample_name.bam; done;
```

The ribosome dataset was aligned again to the rDNA repeats to determine the approximate levels of rRNA still remaining in the data. It was determined that while it is probably not necessary to remove such reads mapping as rRNA, it would still be beneficial to see how successful the rRNA depletion step was. '

```shell
ls data | while read file; do sample_name=$(echo $file | sed "s/.fastq.gz//"); echo $file; novoalign -c 15 -d rRNA_link.nix -f data/$file -F STDFQ -a AGATCGGAAGAGCACACGTCT -l 17 -h -1 -1 -t 90 -g 50 -x 15 -o SAM -o FullNW -r All 51 -e 51 2> 1_alignments/novoalign/rRNA/$file.summary | samtools sort -o 1_alignments/novoalign/rRNA/$sample_name.bam; done;
```

#### 1.2.3. Alignment with HISAT2

Another set of alignments were produced using HISAT2 in order to compare the two mapping softwares. It was believed that Novoalign would be better for the alignment of the ribo-seq dataset due to the short read lengths.

#### Trimming

Before the alignment could be performed, the adapter sequences were trimmed as HISAT2 is not capable of adapter trimming on its own, unlike Novoalign.

```shell
ls data | while read file; do sample_name=$(echo $file|sed "s/.1.fastq.gz//"); echo $file; trimmomatic SE -threads 15 -phred33 -summary 0_trimming/$sample_name.trim.summary data/$file 0_trimming/$sample_name.trimmed.fastq ILLUMINACLIP:$HOME/bin/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:25 MINLEN 17; done;
```

The trimmed sequences were passed to FastQC to evaluate the trimming results.

```shell
ls 0_trimming/ | grep -v ".summary" | while read file; do echo $file; fastqc -o 0_fastqc/trimmed/ --noextract 0_trimming/$file; done;
```

#### Alignment

The trimmed sequences were then aligned to the reference genome.

```shell
ls data | while read file; do sample_name=$(echo $file | sed "s/.1.fastq.gz//"); echo $file; hisat2 -p 15 -q --phred33 --new-summary --summary-file 1_alignments/hisat2/summary/$sample_name.sam.summary --dta --rna-strandness R --non-deterministic --max-intronlen 2000000 -x ../reference/hisat2/genome_snp_tran -U data/$file | samtools sort -o 1_alignments/hisat2/$file_folder/$sample_name.bam; done;
```

***NOTE***: The alignments with HISAT2 were extremely poor (<1% alignment).

## 2. Gene counts

### 2.1. Polysome data

#### 2.1.1. Read counts of polysome data with Stringtie

While it is possible to assemble the transcripts from the alignments to produce a new annotation GTF file and detecting novel transcripts, this is not necessary for this project as we are only interested in the differential expression and translation of known genes between the control and case groups. Therefore, Stringtie was used directly with the original Gencode annotation file as the reference to produce read counts.

```shell
cat sample_list.txt | while read line; do echo $line; stringtie -G /raidset/reference/scanb/hg38.giab_gencode41_snp155/processed/gencode.v41.primary_assmebly.annotation.ucsc.filtered.gtf -e -B -A stringtie/counts/${line}_gene_abund.tab -o stringtie/counts/$line.gtf HISAT2/genome_scanb/$line.bam; done;
```

Since Stringtie can only produce FPKM and TPM values, a Python script (prepDE.py3) provided in the Stringtie package was used to calculate hypothetical read counts using the following formula:

Hypothetical read count = Coverage * Transcript length / Read length

The script can be run in the directory with all of the subdirectories where the GTF files resulting from the last step are stored. In this case "stringtie/counts" directory.

```shell
prepDE.py3
```

#### 2.1.1. Read counts of polysome data with HTSeq

Because stringtie can only calculate hypothetical read counts instead of producing raw read counts, HTSeq, which is capable of producing raw counts, was used to quantify the trascripts once again, in order to compare the count values produced.

```shell
files=$(cat alignments_list.txt | awk '{print}' ORS=" "); htseq-count -f bam -r pos -s reverse -t exon -m intersection-strict --nonunique=all $files ../reference/raw/gencode.v41.primary_assembly.annotation.ucsc.filtered.gtf
```

### 2.2. Ribosome data

#### 2.2.1. Read counts of ribosome data with HTSeq

```shell
files=$(cat alignments_list.txt | awk '{print}' ORS=" "); htseq-count -f bam -r pos -s reverse -t exon -m intersection-strict --nonunique=all $files ../reference/raw/gencode.v41.primary_assembly.annotation.ucsc.filtered.gtf
```

## 3. Differential expression analysis

### 3.1. Polysome data

The R package Anota2seq was used for the differential expression analysis of the polysome data, generating log2FC values for the total mRNA and translated (polysome associated) mRNA.

### 3.2. Ribosome data

Anota2seq could not be used for the analysis of the Ribo-seq data, as it requires at least 3 replicates per condition and one of the samples was removed due to a significantly low total number of reads. Therefore, similar analysis packages: Xtail and DESeq2 were used.

## 4. Gene set enrichment analysis (GSEA)

The R package clusterProfiler was used.

