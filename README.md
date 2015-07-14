# cage-dev
**C**RISPR KO **A**nalysis based on **G**enomic **E**diting data (development)

## Introduction
A CRISPR-cas9 based Genome Editing data analysis pipeline, for the analysis of indels and microhomology patterns, the identification of personalized features correlated to sgRNA KO efficiency on heterogeneous experimental conditions, and the evaluation of the sgRNA KO efficiency based on the CRISPR-Cas9
Knock-Out NGS data or the sgRNA KO assay data.

The ultimate goals of CAGE are (1) CAGE provides a standard CROWDSOURCING platform for the users to share the CRISPR-Cas9 based gene KO data, (2) CAGE provides an efficient interface to analysis and visualize the CRISPR-based KO NGS data, (3) CAGE provides a robust learning pipeline to derive the sequence determinants from heterogeneous genome editing data for different cell types and organisms, and (4) CAGE provides an personalized scoring framework for on-target sgRNA design based on the derived sequence determinants for specific cell types or organisms.

Currently CAGE records the optimal sgRNA KO efficiency prediction models and the personalized score functions in sgRNA design for the following XXX cell types. The optimal results for new cell types as well as the the current ones will be updated timely.  


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
* LaTeX (for visualization)

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
1. `sg`    Process sgRNA sequences into sgRNA information table
2. `prep`  Process NGS data into sgRNA-Indel Table
3. `mh`    Microhomology Detection
4. `indel` Feature selection and model prediction on sgRNA OTF ratio based on NGS data
5. `fs`    Feature selection and model prediction on clearly defined sgRNA KO efficiency
6. `eval`  sgRNA KO efficiency evaluation
7. `vis`   Visualization of feature selection result

## sgRNA processing
```bash
python cage.py sg -s <sgRNA.fq>
	              -o <output directory>
                  -g <reference genome>
				  -t <bwa threads>
```
For more detail on the options, see `python cage.py sg -h`.

## NGS data preprocessing
* Single-end
```bash
python cage.py prep -s <sg file>
	                -f <reads.fq>
	                -o <output directory>
                    -g <reference genome>
					-t <bwa threads>
```

* Paired-end
```bash
python cage.py prep -s <sg file>
                    -f <reads_1.fq>
					-r <reads_2.fq>
					-o <output directory>
					-g <reference genome>
					-t <bwa threads>
```
For more detail on the options, see `python cage.py prep -h`.

## Microhomology detection
```bash
python cage.py mh -i <samind file>
                  -o <output directory>
	              -g <reference genome>
```
For more detail on the options, see `python cage.py mh -h`.

## Feature selection and model prediction on sgRNA OTF Ratio based on NGS data
```bash
python cage.py indel -i <samind file>
                     -s <sg file>
                     -o <output directory>
	                 -g <reference genome>
```
For more detail on the options, see `python cage.py indel -h`.

## Feature selection and model prediction on clearly defined sgRNA KO efficiency 
```bash
python cage.py fs -i <label file>
                  -s <sg file>
                  -o <output directory>
	              -g <reference genome>
				  -m <lasso|logit>
```
For more detail on the options, see `python cage.py fs -h`.

## sgRNA KO efficiency evaluation
```bash
python cage.py eval -s <sg file>
                    -f <score function file>
                    -o <output directory>
					-g <reference genome>
```
For more detail on the options, see `python cage.py eval -h`.

## Visualization
```bash
python cage.py vis -f <feature report file>
                   -o <output directory>
```
For more detail on the options, see `python cage.py vis -h`.				

## Test
For commands testing, `cd test` first, then execute the following
commands.

* Testing `sg`: `sh test.sh sg`
* Testing single-end `prep`: `sh test.sh prep_se`
* Testing pair-end `prep`: `sh test.sh prep_pe`
* Testing `mh`: `sh test.sh mh`
* Testing `indel` without auto detection: `sh test.sh indel`
* Testing `indel` with auto detection: `sh test.sh indel_a`
* Testing `fs` using *LASSO* without auto detection: `sh test.sh
fs_las`
* Testing `fs` using *LASSO* with auto detection: `sh test.sh fs_las_a`
* Testing `fs` using *Logistic Regression* without auto detection: `sh test.sh
fs_log`
* Testing `fs` using *Logistic Regression* with auto detection: `sh test.sh
fs_log_a`
* Testing `eval`: `sh test.sh eval`
* Testing `vis`: `sh test.sh vis`
