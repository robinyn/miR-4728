#  Identification of targets for the microRNA miR-4728-3p in ribosome and polysome sequencing data

## Data

The data for this project consists of ribosome and polysome profiling data for SKBR3 (HER2-positive breast cancer cell line) treated with antisense oligonucleotide (ASO) to block the miRNA or a control oligonucleotide, with matched total RNA-seq. There are three replicates per condition.

## Analysis environment

A conda environment was created for this analysis and the environment file is attached.

- samtools (v)
- gfftools (v)

Due to dependency issues, a separate conda environment was created for HTseq. The environment file is attached.

- htseq (v)

The following tools were not available through Conda and were installed manually following the installation steps provided in their respective documentations:

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

#### 1.1.4. Read counts with Stringtie

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

#### 1.1.5. Gene feature quantification with HTSeq

Because stringtie can only calculate hypothetical read counts instead of producing raw read counts, HTSeq, which is capable of producing raw counts, was used to quantify the trascripts once again, in order to compare the count values produced.

```shell
cat sample_list.txt | while read line; do sample_name=$(echo $line | cut -d "/" -f 1); echo $line; htseq-count -f bam -r pos -s reverse -t exon -m intersection-nonempty  --nonunique=none --addtional-attr=gene_name HISAT2/genome_scanb/$line.bam /raidset/reference/scanb/hg38.giab_gencode41_snp155/processed/gencode.v41.primary_assembly.annotation.ucsc.filtered.gtf > htseq/$sample_name.count; done;
```

### 1.2. Ribosome profiling data

#### 1.2.1. Building the index

In order to build an index file for Novoalign with known transcripts and splice sites, the protocol outlined in the [Novoalign documentation](https://www.novocraft.com/documentation/novoalign-2/novoalign-user-guide/rnaseq-analysis-mrna-and-the-spliceosome/) was used.

#### Preparing genome and transcript annotation files

The included Novoindex software, for generating the index file for Novoalign to use as the reference while mapping, requires three input FASTA files:

    1. Splice sites file
    2. Transcripts file
    3. Masked genome file (where genes are masked as Ns)

To generate the first two, the **MakeTranscriptome** program from the **Useq** package is required. The following [documentation](https://useq.sourceforge.net/usageRNASeq.html) was used to generate the FASTA files.

