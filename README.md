# cage-dev
CRISPR KO Analysis based on Genomic Editing data (development)

## Introduction
A CRISPR-cas9 based Genome Editing data analysis pipeline, for the
analysis of indels and microhomology patterns from CRISPR-Cas9
Knock-Out NGS data.

## Implementation
* Python >= 2.7
* Numpy >= 1.9.2
* Scipy >= 0.15.1
* Pandas >= 0.16.0
* scikit-learn >= 0.16.1
* lxml >= 3.4.4
* pyfasta >= 0.5.2
* bwa >= 0.7.12
* samtools >= 0.1.19
* bedtools >= 2.23.0

## Presetting
Make sure to perform this presetting carefully. Because reference setting
is **very important**.

For the sake of simplicity, we use hg19 as the example.

1. Download the hg19 genome(fasta file) from UCSC, put it in certain
   directory, name it `hg19.fa` and set the directory path as
   `$FASTADB`.

2. Generate `bwa` index files from `hg19.fa`, put them in certain directory
   and set the directory path as `$BWADB`.
   
3. (Optional) Download the hg19 gene annotation files from UCSC, convert it to
   `bed-6` format with the 4th column being the gene name, put them in
   certain directory and set the directory path as `$BEDDB`. Here are the
   renamed file:
>>
File | Standard | Requirement
------|-----|-----
hg19ref.bed | Refseq |required
hg19ucsc.bed | UCSC Gene | optional
hg19gencode.bed | GENCODE | optional



## Installation
```
git clone https://github.com/bm2-lab/cage-dev.git
```

## Usage
```bash
python cage.py <command> [option] ...
```

## Command
1. `sg` for sgRNA processing
2. `prep` for preprocessing
3. `mh` for microhomology detection
4. `indel` for indel frameshifting paradigm using **LASSO**
5. `las` for general data source using **LASSO**
6. `vis` for result visualization

## sgRNA Processing
```bash
python cage.py sg -s <sgRNA.fq>
	              -d <target directory>
                  -g <reference genome>
				  -t <bwa threads>
```
For more detail on the options, see `python cage.py sg -h`.

## Preprocessing
* Single-end
```bash
python cage.py prep -s <sg file>
	                -f <reads.fq>
	                -d <target directory>
                    -g <reference genome>
					-t <bwa threads>
```

* Paired-end
```bash
python cage.py prep -s <sg file>
                    -f <reads_1.fq>
					-r <reads_2.fq>
					-d <target directory>
					-g <reference genome>
					-t <bwa threads>
```
For more detail on the options, see `python cage.py prep -h`.

## Microhomology Detection
```bash
python cage.py mh -i <samind file>
                  -d <target directory>
	              -g <reference genome>
```
For more detail on the options, see `python cage.py mh -h`.

## Indel Analysis
```bash
python cage.py indel -i <samind file>
                     -s <sg file>
                     -d <target directory>
	                 -g <reference genome>
```
For more detail on the options, see `python cage.py indel -h`.

## General Lasso Feature Selection
```bash
python cage.py las -i <label file>
                   -s <sg file>
                   -d <target directory>
	               -g <reference genome>
```
For more detail on the options, see `python cage.py las -h`.

## Visualization
```bash
python cage.py vis -f <feature report file>
                   -d <target directory>
```
For more detail on the options, see `python cage.py vis -h`.				
