#  Identification of targets for the microRNA miR-4728-3p in ribosome and polysome sequencing data

## Data

The data for this project consists of ribosome and polysome profiling data for SKBR3 (HER2-positive breast cancer cell line) treated with antisense oligonucleotide (ASO) to block the miRNA or a control oligonucleotide, with matched total RNA-seq. There are three replicates per condition.

## Analysis environment

A conda environment was created for this analysis and the environment file is attached.

- samtools (v)
- gfftools (v)
- htseq (v)

The following tools were not available and were installed manually following the installation steps provided in their respective documentations:

- HISAT2 (v2.2.1)
- stringtie (v2.2.1)

## 1. Alignment of reads

### 1.1. Polysome profiling data

#### 1.1.1. Building index

The GRCh38.p13 assembly of the human genome from the Genome Reference Consortium and the Gencode 43 release of the human transcriptome annotation (Ensembl release 109) was used to build a reference index for the alignment, following the [protocol](http://daehwankimlab.github.io/hisat2/howto/) provided in the HISAT2 documentation.

Unfortunately, the indexing of a human genome with transcript annotations required over 160GB of memory and could not be completed. Therefore, a prebuilt index provided by [HISAT2](http://daehwankimlab.github.io/hisat2/download/) was used for the alignment.

Due to the lack of documentation on the version of the genome annotation used for the prebuilt index files, a comparison between the exons and splice sites included in the prebuilt index to the ones in the Gencode 43 release of the human transcriptome annotation was made.

#### Extracting the exons and splice sites from prebuilt index

The HISAT2 team provides a script to inspect and extract information from a set of index files generated using their suite.

```shell
hisat2-inspect --ss-all grch38_tran/genome_tran > prebuilt.ss
hisat2-inspect --exon grch38_tran/genome_tran > prebuilt.exons
```

#### Extracting the exons and splice sites from Gencode annotation

Similarly, the tools designed for building index files from genome assemblies and annotation files can be used to generate the same files for the Gencode annotation.

```shell
hisat2_extract_splice_sites.py gencode.v43.annotation.gtf > gencode.ss
hisat2_extract_splice_sites.py gencode.v43.annotation.gtf > gencode.exons
```

Using the generated ***.exons*** and ***.ss*** files, a simple comparison between the number of lines can be made to identify the difference in the number of splice sites and exons.

```shell
wc -l file_name.exons
wc -l file_name.ss
```

| Type | Exons | Splice sites |
| ---- | ---- | ---- |
| Prebuilt index | 308795 | 347295 |
| Gencode v43 Primary assembly | 328613 | 402486 |

Unique scaffolds/chromosomes with annotated exons from the two ***.exon*** files were identified and compared.

```shell
cat prebuilt.exons| cut -f 1 | cut -d " " -f 1 | sort -V | uniq > prebuilt_exons_scaffolds.txt

cat gencode.exons | cut -f 1 | sed "s/chrM/MT/" | sed "s/chr//" |sort -V | uniq > gencode_exons_scaffolds.txt

diff gencode_exons_scaffolds.txt prebuilt_exons_scaffolds.txt
```

There were 12 scaffolds present in the prebuilt index that were not found in the Gencode annotation.

#### 1.1.2. Alignment

The reads from the polysome profiling were aligned to the GRCh38 Human Genome as the reference, using [HISAT2](http://daehwankimlab.github.io/hisat2/manual/). The reference [index](http://daehwankimlab.github.io/hisat2/download/#h-sapiens) provided by HISAT2. The max intron length was set at 2Mb (REFERENCE). The alignment output was directly piped into samtools to sort and convert them from SAM to BAM formats.

```shell
cat files_list.txt | while read line;do file_dir=$(echo ../data/$line); file_name=$(echo $line | cut -d "/" -f 2); file_folder=$(echo $line | cut -d "/" -f 1); mate_1=$(echo ${file_dir}_R1_001.fastq.gz); mate_2=$(echo ${file_dir}_R2_001.fastq.gz); echo $file_name; hisat2 -p 15 -q --fr --new-summary --summary-file genome_tran/summary/$file_name.sam.summary --dta --rna-strandness RF --non-deterministic --max-intronlen 2000000 -x ../../reference/grch38_tran/genome_tran -1 $mate_1 -2 $mate2 | samtools sort -o genome_tran/$file_folder/$file_name.bam; done;
```

A different index file, with a more up-to-date transcript annotation and SNP data was obtained from ScanB and was used to generate a new set of alignments.

```shell
cat files_list.txt | while read line;do file_dir=$(echo ../data/$line); file_name=$(echo $line | cut -d "/" -f 2); file_folder=$(echo $line | cut -d "/" -f 1); mate_1=$(echo ${file_dir}_R1_001.fastq.gz); mate_2=$(echo ${file_dir}_R2_001.fastq.gz); echo $file_name; hisat2 -p 15 -q --fr --new-summary --summary-file genome_tran/summary/$file_name.sam.summary --dta --rna-strandness RF --non-deterministic --max-intronlen 2000000 -x ../../reference/grch38_tran/genome_tran -1 $mate_1 -2 $mate2 | samtools sort -o genome_tran/$file_folder/$file_name.bam; done;
```

#### 1.1.3. Summary statistics from alignment

A list of summary files and their directories was generated using bash.

```shell
ls summary_files | while read file; do echo summary_files/$file; done > sample_list.txt
```

A custom Python script was used to parse the summary files with the list and visualized using R. The two sets of alignments, with the prebuilt and scanb index, were compared. No significant differences were identified, therefore the alignment with the most up-to-date reference index was used for the subsequent analyses.

#### 1.1.4. Transcriptome assembly and gene feature quantification with Stringtie

While it is possible to skip the transcript assembly steps, since the goal of the project is to identify targets of the micro RNA in question, it would be beneficial to identify any novel transcripts even if the human genome has already been extensively annotated. Therefore, the assembly step was included in the analysis. The aligned reads were used to assemble the transcripts for each sample.

```shell
cat sample_list.txt | while read line; do echo $line; stringtie -G /raidset/reference/scanb/hg38.giab_gencode41_snp155/processed/gencode.v41.primary_assembly.annotation.ucsc.filtered.gtf -o stringtie/assembly/$line.gtf HISAT2/genome_scanb/$line.bam; done;
```

Then, the assembled transcripts were merged to create a comprehensive and consistent annotation of all of the gene structures found in the all of the samples, so that transcripts can be compared across the samples in subsequent analyses.

```shell
# First a list of GTF files to be merged has to be made
cat sample_list | while read line; do echo stringtie/assembly/$line.gtf; done > merge_list.txt

# The list is passed to stringtie to merge
stringtie --merge -G /raidset/reference/scanb/hg38.giab_gencode41_snp155/gencode.v41.primary_assembly.annotation.ucsc.filtered.gtf --rf -o stringtie/stringtie_merged.gtf merge_list.txt
```

The merged file was compared against the original reference annotation to identify novel transcripts.

```shell
gffcompare -r /raidset/reference/scanb/hg38.giab_gencode41_snp155/gencode.v41.primary_assembly.annotation.ucsc.filtered.gtf -o merged stringtie/stringtie_merged.gtf
```

Finally, the merged GTF file was used as a new reference to quantify the transcripts from the aligned reads.

```shell
cat sample_list.txt | while read line; do echo $line; stringtie -G stringtie/stringtie_merged.gtf -e -B -A stringtie/counts/$line_gene_abund.tab -o stringtie/counts/$line.gtf HISAT2/genome_scanb/$line.bam; done;
```

As stringtie is only capable of producing FPKM and TPM values for gene counts instead of the raw counts necessary for downstream analysis with DESeq2, a Python script (prepDE.py3) provided in the stringtie package was used to calculate hypothetical transcript read counts based on the following formula:

    Hypothetical read count = Coverage * Transcript length / Read length



#### 1.1.5. Gene feature quantification with HTSeq

Because stringtie can only calculate hypothetical read counts instead of producing raw read counts, HTSeq, which is capable of producing raw counts, was used to quantify the trascripts once again, in order to compare the count values produced.

