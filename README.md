# CAGE
**C**RISPR KO **A**nalysis based on **G**enomic **E**diting data —— Toward personalized sgRNA design in heterogeneous experimental conditions

## Introduction
A CRISPR-cas9 based genome-editing data analysis resource and platform, for the analysis of indels and microhomology patterns, the identification of personalized features correlated to sgRNA KO efficiency on heterogeneous experimental conditions, and the evaluation of the sgRNA KO efficiency based on the CRISPR-Cas9
KO NGS data or the sgRNA KO assay data.

The ultimate goals of CAGE:

1. CAGE provides a standard CROWDSOURCING platform for users to
share the CRISPR-Cas9 based gene KO data.

2. CAGE provides an efficient interface to analyze and visualize the
CRISPR-based KO NGS data.

3. CAGE provides a robust learning pipeline to derive the sequence
determinants from heterogeneous genome editing data.

4. CAGE provides an personalized evaluation framework for on-target sgRNA design based on the derived sequence determinants torward specific cell types or organisms.

Currently CAGE records the optimal sgRNA KO efficiency prediction
models and the personalized evaluation models in sgRNA design for the
following different cell types. The optimal results for new cell types as well
as the the current ones will be updated regularly. Users can select the existing evaluation model for a specific cell type for sgRNA design, or they can use their own sgRNA KO data to generate a new personalized evaluation model for a new cell type for furthur sgRNA KO efficiency evaluation.

Evaluation Model| Species | Cell Type | KO Efficiency Measurement | Data Type | Learning Model | Performance | Actual sgRNA Library Size | Accession | Time Stamp
---------|-----|-----|-------|---------|-------|------|-------|----|------
a375_1 | Homo sapiens | A375 | See [Doench et al.][1] | numerical | LASSO | r2=0.64 | 1248 | - | 2015-7-16
el4_1 | Mus musculus | EL4 | See [Doench et al.][1] | numerical | LASSO | r2=0.65 | 858 | - | 2015-7-16
mesc_1 | Mus musculus | mESC | OTF Ratio | numerical | LASSO | r2=0.72 | 99 | [ERP003292][2] | 2015-7-16
rn2c_1 | Mus musculus | RN2c | OTF Ratio | numerical | LASSO | r2=0.89| 26 | [SRP057117][3] | 2015-7-16
hela_1 | Homo sapiens | Hela | OTF Ratio | numerical | LASSO | r2=0.87| 68 | [SRP042061][4] | 2015-7-16
dr_1 | Danio rerio | *AB/Tu | OTF Ratio | numerical | LASSO | r2=0.91| 47 | [SRP052749][5] | 2015-7-16
hl60_nonribo | Homo sapiens | HL60 | See [Xu et al.][6] | categorical|Logistic | AUC=0.76 | 908 | - | 2015-7-16
hl60_ribo | Homo sapiens | HL60 | See [Xu et al.][6] | categorical| Logistic | AUC=0.79 | 373 | - | 2015-7-16
mesc_2 | Mus musculus | mESC | See [Xu et al.][6] | categorical | Logistic | AUC=0.81 | 1028 | - | 2015-7-16
hek293t_1 | Homo sapiens | HEK293T | See [Chari et al.][7] | categorical| Logistic | AUC=0.77 | 279 | - | 2015-7-23

[1]: http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html

[2]: http://www.ebi.ac.uk/ena/data/view/ERP003292

[3]: http://www.ncbi.nlm.nih.gov/sra/?term=SRP057117

[4]: http://www.ncbi.nlm.nih.gov/sra/?term=SRP042061

[5]: http://www.ncbi.nlm.nih.gov/sra/?term=SRP052749

[6]: http://genome.cshlp.org/content/early/2015/07/01/gr.191452.115

[7]: http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3473.html


## Implementation
* Python 2.7
* Numpy >= 1.9.2
* Scipy >= 0.15.1
* Pandas >= 0.16.0
* scikit-learn >= 0.16.1
* lxml >= 3.4.4
* pyfasta >= 0.5.2
* bwa >= 0.7.12
* samtools >= 0.1.19
* bedtools >= 2.23.0
* pyslep (for multi-task group lasso)
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

##### Install pyslep
For performing multi-task group lasso, pyslep is necessary.

1. `cd pyslep`
2. `sh setup_pyslep.sh`

## Usage
```
python cage.py <command> [option] ...
```

## Command
1. `sg`    Process sgRNA sequences into sgRNA information table
2. `prep`  Process NGS data
3. `mh`    Microhomology Detection
4. `indel` Indel frameshifting paradigm analysis
5. `fs`    Feature selection and model prediction on clearly defined sgRNA KO efficiency
6. `mt`    Feature selection with multi-task group LASSO on clearly defined sgRNA KO efficiency for cross-platform data
7. `eval`  sgRNA KO efficiency evaluation and the scanning of a given genome region for sgRNA design
8. `vis`   Visualization of feature selection result

## Data File Format
File Type | Suffix | Usage
----|----|----
sg file | .sg | sgRNA information table
samind file | .samind | reads mapping result
mnst file | .mnst | microhomology information table
iost file | .iost | sgRNA-indel table
seq file | .seq | original sequence feature table
fesrep file | _fesrep.xml | feature selection and model prediction report
pkl file | .pkl | evaluation model file
st file | .st | evaluation result file
label file | (arbitrary) | user-customized evaluation file for feature selection and model prediction

#####Note
For label file, the first and the last column will be regarded as
sgID and score respectively. File header should exist and the header
of the first column must be sgID. See the following example.

sgID | ... | score
----|----|----
sg1 | ... | 0.1
sg2 | ... | 0.2

## sgRNA processing
Generate sgRNA Information Table (sg file)
```
python cage.py sg -s <sgRNA.fq>
	              -o <output directory>
                  -g <reference genome> (e.g. hg19)
				  -t <bwa threads> (default 1)
				  -a (annotate, optional)
```
For more detail on the options, see `python cage.py sg -h`.

## NGS data preprocessing
* Single-end
```
python cage.py prep -s <sg file>
	                -f <reads.fq>
	                -o <output directory>
                    -g <reference genome>
					-t <bwa threads>
```

* Paired-end
```
python cage.py prep -s <sg file>
                    -f <reads_1.fq>
					-r <reads_2.fq>
					-o <output directory>
					-g <reference genome>
					-t <bwa threads>
```
For more detail on the options, see `python cage.py prep -h`.

## Microhomology detection
```
python cage.py mh -i <samind file>
                  -o <output directory>
	              -g <reference genome>
```
For more detail on the options, see `python cage.py mh -h`.

## Indel frameshifting paradigm analysis
Generate sgRNA-indel table (iost file)
```
python cage.py indel -i <samind file>
                     -s <sg file>
                     -o <output directory>
	                 -g <reference genome>
					 -t <read-count cutoff> (default 0)
```
For more detail on the options, see `python cage.py indel -h`.

## Feature selection and model prediction on clearly defined sgRNA KO efficiency
* Manual
```
python cage.py fs -i <label file>
                  -s <sg file>
                  -o <output directory>
	              -g <reference genome>
				  -t <reads cutoff> (default 0)
				  -u <upstream region length> (default 30)
				  -w <downstream region length> (without PAM, default 27)
				  -c <cross-validation folds> (default 5)
				  -j <number of CPU cores used> (default 1)
				  -m <lasso|logit>
```

* Auto
```
python cage.py fs -i <label file>
                  -s <sg file>
                  -o <output directory>
	              -g <reference genome>
				  -t <reads cutoff> (default 0)
				  -a (auto detection for sequence region)
				  --init-radius <init radius> (default 0)
				  -r <radius> (default 200)
				  --step <detection step> (default 5)
				  -c <cross-validation folds> (default 5)
				  -j <number of CPU cores used> (default 1)
				  -m <lasso|logit> (method)
```
For more detail on the options, see `python cage.py fs -h`.

## Feature selection with multi-task group LASSO on clearly defined sgRNA KO efficiency for cross-platform data
```
python cage.py mt -i [<label file> [<label file> ...]]
                  -s [<sg file> [<sg file> ...]]
                  -o <output directory>
	              -g [<ref genome> [<ref genome> ...]]
				  -u <upstream region length> (default 30)
				  -w <downstream region length> (without PAM, default 27)
				  -d <selection strictness>
```

## sgRNA KO efficiency evaluation and the scanning of a given genome region for sgRNA design
* Evaluation with Genome Scanner
```
python cage.py eval -c <target chromosome>
                    -b <start coordinate>
                    -e <end coordinate>
                    -m <pkl file>
                    -o <output directory>
					-g <reference genome>
					-d <two-sided|pos|neg> (scan direction)
					-t <bwa threads>
```

* Evaluation with sgRNA Information Table
```
python cage.py eval -s <sg file>
                    -m <pkl file>
                    -o <output directory>
					-g <reference genome>
```
For more detail on the options, see `python cage.py eval -h`.

## Visualization
```
python cage.py vis -f <feature report file>
                   -o <output directory>
```
For more detail on the options, see `python cage.py vis -h`.				

## Example
To run examples, `cd example` first, then execute the following commands.

* `sg`: `sh exam.sh sg`
* `prep` for single-end: `sh exam.sh prep_se`
* `prep` for pair-end: `sh exam.sh prep_pe`
* `mh`: `sh exam.sh mh`
* `indel`: `sh exam.sh indel`
* `fs` using *LASSO* without auto detection: `sh exam.sh fs_las`
* `fs` using *LASSO* with auto detection: `sh exam.sh fs_las_a`
* `fs` using *Logistic Regression* without auto detection: `sh exam.sh
  fs_log`
* `fs` using *Logistic Regression* with auto detection: `sh exam.sh
fs_log_a`
* `mt`: `sh exam.sh mt`
* `eval` using genome scanner: `sh exam.sh eval`
* `eval` using sgRNA information table: `sh exam.sh eval_sg`
* `vis`: `sh exam.sh vis`
