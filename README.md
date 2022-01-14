# DeBreak

A SV caller for long-read single-molecular sequencing data.

Author: Maggi Chen

Email: maggic@uab.edu

Draft date: Jan. 12, 2022

## Quick Start
```sh
git clone https://github.com/Maggi-Chen/DeBreak.git
cd DeBreak/
./debreak -h

# quick SV calling with sorted bam
debreak --bam merged.sort.bam --outpath debreak_out/

# Accurate SV calling with sorted bam (reference genome needed)
debreak --bam merged.sort.bam --outpath debreak_out/ --rescue_large_ins --rescue_dup --poa --ref hg38.fa 

# SV discovery in cancer/complex genome
debreak --bam merged.sort.bam --outpath debreak_out/ --rescue_large_ins --poa --ref hg38.fa --tumor

```



## Description

DeBreak is a tool for SV discovery with long read data. The input should be a sorted BAM (PacBio CLR, PacBio HiFi, Oxford Nanopore, or mixed platform). The reference genome is also needed when enabling full functions. The output is a standard VCF file containing all SVs that have passed filters. This program was tested on a x86_64 Linux system with a 128GB physical memory.


## Depencency

Dependencies for DeBreak:

* python 2.7  
* pysam  (tested with version 0.17.0)
* minimap2  (tested with version 2.10 and 2.15)
* wtdbg2  (tested with version 2.5)



## Installation

```
git clone https://github.com/Maggi-Chen/DeBreak.git
```
Then, please also add this directory to your PATH:
```
export PATH=$PWD/DeBreak/:$PATH
```


To simplify the environment setup process, Anaconda2 (https://www.anaconda.com/) is recommended.
To create an environment with conda:
```
conda create --name deb python=2.7
conda activate deb
conda install -c bioconda minimap2=2.15
conda install -c bioconda samtools=1.9
conda install -c bioconda pysam=0.16.0.1
conda install -c bioconda wtdbg=2.5

```

A test dataset is available to verify successful installation:
```
cd DeBreak/
debreak --bam testdata/test_read.bam -o test_out/ --poa --rescue_large_ins --rescue_dup --ref testdata/test_ref.fa
```
(The DeBreak SV discovery on test dataset should finish within several minutes with 4 CPUs and 400MB memory.)


## General usage


```
debreak [-h] --bam <sort.bam>

SV caller for long-read sequencing data

optional arguments:
  -h, --help                       Show this help message and exit
  -v, --version                    Show program's version number and exit
  --bam BAM                        Input sorted bam. index required
  --samlist SAMLIST                A list of SAM files of same sample
  -o, --outpath OUTPATH            Output directory
  -p, --prefix PREFIX              Prefix of output
  --min_size MIN_SIZE              Minimal size of detected SV
  --max_size MAX_SIZE              Maxminal size of detected SV
  -d, --depth DEPTH                Sequencing depth of this dataset
  -m, --min_support MIN_SUPPORT    Minimal number of supporting reads for one event
  --min_quality MIN_QUALITY        Minimal mapping quality of reads
  --aligner ALIGNER                Aligner used to generate BAM/SAM
  -t, --thread THREAD       Number of threads
  --rescue_dup                     Rescue DUP from INS calls. minimap2,ref required
  --rescue_large_ins               Rescue large INS. wtdbg2,minimap2,ref required
  --poa                            POA for accurate breakpoint. wtdbg2,minimap2,ref required
  --no_genotype                    Disable genotyping
  -r, --ref REF                Reference genome. Should be same with SAM/BAM
  --maxcov MAXCOV                  Maximal coverage for a SV. Suggested maxcov as 2 times mean depth
  --skip_detect                    Skip SV raw signal detection
  --tumor                          Allow clustered SV breakpoints during raw SV signal detection

```

## Use cases
DeBreak requires a input of read alignment results in BAM format. If you start with sequencing reads (Fasta or Fastq format), you may use minimap2 and samtools to map them to a reference genome before you can apply DeBreak:
```
minimap2 -a reference.fa  movie1.fastq | samtools sort -o movie1.bam
minimap2 -a reference.fa  movie2.fastq | samtools sort -o movie2.bam
...
samtools merge merged.sort.bam  movie1.bam movie2.bam
samtools index merged.sort.bam
```

DeBreak can be applied with full function (with accurate SV breakpoints):
```
debreak --bam merged.sort.bam -o debreak_out/ --rescue_large_ins --rescue_dup --poa --ref hg38.fa
```
Or with only basic functions to quickly call SVs (without accurate SV breakpoints, may miss some large insertions):
```
debreak --bam merged.sort.bam -o debreak_out/
```

### Options of Debreak
#### 1. --min_support, minimal number of supporting reads
min_support is the most important argument for SV filtering of DeBreak. It should be adjusted according to sequencing depth of input BAM.
Without given min_support/depth information, DeBreak estimates the depth of input dataset to assign a resonable min_support.
```
debreak --bam merged.sort.bam -o debreak_out/ --min_support 5 
```
Suggested min_support at each depth:

| Depth    | Suggested min_supp   |
| -------- |:--------------------:|
| 10       | 3                    |
| 20       | 4                    |
| 30       | 5                    |
| 40       | 6                    |
| 50       | 7                    |
| 60       | 8                    |
| 70       | 9                    |
| 80       | 10                   |

If you specify coverage of dataset with --depth, DeBreak will calculate min_support according to the table.
```
debreak --bam merged.sort.bam -o debreak_out/ --depth 70
```

#### 2. --poa, partial order alignment for SV breakpoint refinement
This option calls POA module of DeBreak. minimap2, wtdbg2, and reference genome are required.
After clustering and filtering SV calls, DeBreak collects all SV-containing reads for each SV candidate and performs POA with wtdbg2 to generate a consensus sequence with higher base accuracy. It then maps consensus sequences to reference genome and infers precise SV breakpoint positions.
```
debreak --bam merged.sort.bam -o debreak_out/ --poa --ref reference.fa
```

#### 3. --rescue_large_ins, detect large insertions using local de novo assembly
This option calls rescue_large_insertion module of DeBreak to identify insertions that are longer than the sequencing reads. minimap2, wtdbg2, and reference genome are required. 
DeBreak scans whole genome for read alignments with clipped end. It identifies candidate insertion breakpoints with enriched clipped alignments, and performs local de novo assembly at each candidate insertion site.
minimap2 is used to map all assembled contigs to the reference. DeBreak then detects ultra-large insertions from the contig alignment results.
```
debreak --bam merged.sort.bam -o debreak_out/ --rescue_large_ins --ref reference.fa
```

#### 4. --rescue_dup, rescue duplications that are falsely considered as insertions
This option calls rescue_duplication module of DeBreak. minimap2 and reference genome are required.
This module maps inserted sequence of each insertion call back to local region near insertion breakpoint. If the inserted sequence can be properly mapped, it suggests that this is a duplication instead of novel insertion.
```
debreak --bam merged.sort.bam -o debreak_out/ --rescue_dup --ref reference.fa
```


#### 5. --tumor, tumor mode for cancer genomes
This option sets looser criteria during SV raw signal detection, allowing identification of potentially closer SV breakpoints in complex rearrangements. 
```
debreak --bam tumor.sort.bam -o debreak_out/ --tumor
```


## Output files
The output directory includes:
```
debreak.vcf                   Standard VCF file of SV calls. The chromosome, coordinates, size, type, number of supporting reads, mapping quality, genotype, and multi-allelic information are recorded for each SV call.
debreak-allsv-merged-final    Tab-delimited SV list, containing the name of reads supporting each SV call.
sv_raw_calls/                 Includes all SV raw signals on each chromosome.
debreak_poa_workspace/        Temporary files during SV breakpoint refinement with POA. For debug purpose.
debreak_ins_workspace/        Temporary files during ultra-large insertion detection. For debug purpose.
map_depth/                    Temporary files during sequencing depth estimation. For debug purpose.
```



