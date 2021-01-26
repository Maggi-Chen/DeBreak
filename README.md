# DeBreak

A SV caller for long-read sigle-molecular sequencing reads.

Author: Maggi Chen

Email: maggic@uab.edu

Draft date: Oct. 25, 2018

## Quick Start
```sh
git clone https://github.com/Maggi-Chen/DeBreak
cd DeBreak/
./debreak -h
# Detect SVs from sorted bam
./debreak --bam merged.sort.bam --outpath debreak_out/
# Detect SVs form a list of sam files
ls movie1.sam movie2.sam movie3.sam > all_sam_list
./debreak --samlist all_sam_list --outpath debreak_out_sam/
# Optional modules
./debreak --bam merged.sort.bam --outpath debreak_out/ --no_genotype --rescue_dup --poa 

```



## Description

DeBreak is a tool for Structural Variant calling from long-read sequencing data. The input should be a sorted bam file, or a list of SAM files. The out put will be a vcf file, or a list of SVs in bed format with option --output_bed.

The program was tested on a x86_64 Linux system with a 8GB physical memory.

## Depencency:

python 2.7  
python package:  pysam  

For full function of DeBreak, some tools are required for optional module of DeBreak:
* minimap2  (tested with version 2.10 and 2.15)
* wtdbg2    (tested with version 2.1)
* samtools  (tested with version 1.9)


## Installation

```
git clone https://github.com/Maggi-Chen/DeBreak
```
Then, please also add this directory to your PATH:
```
export PATH=$PWD/DeBreak/:$PATH
```

## General usage


```
debreak --bam <sort.bam> --samlist <list_of_sam> [options]
  required arguments:
  --bam BAM             input sorted BAM. index required
  --samlist SAMLIST     a list of SAM files of same sample

  optional arguments:
  -h, --help            show this help message and exit
  --aligner             aligner used to generate BAM/SAM
  --coverage            coverage of this dataset
  --min_size            minimal size of detected SV
  --max_size            maxminal size of detected SV
  --min_support         minimal number of supporting reads for one event
  --min_quality         minimal mapping quality of reads
  --no_genotype         disable genotyping
  --outpath             output directory
  --prefix              prefix of output
  --ref                 reference genome. Should be same ref with SAM/BAM
  --rescue_dup          rescue duplication from insertion calls. minimap2, required
  --rescue_long_ins     rescue large insertions. wtdbg2, minimap2, ref required
  --poa                 POA for accurate breakpoint. wtdbg2, minimap2, ref required.
  --thread              number of threads

```

## Use cases
DeBreak requires a input of alignment results in SAM or BAM format. If you start with sequencing reads (fasta or fastq format), you may use minimap2 to map them to a reference genome before you can apply DeBreak.
```
minimap2 reference.fa  movie1.fa -o movie1.sam
minimap2 reference.fa  movie2.fa -o movie2.sam
...
```
You can directly use a list of SAM as input:
```
ls movie2.sam movie2.sam > sam_list
debreak --samlist sam_list 
```
Or you can merge and sort all SAM files to generate a sorted BAM file using samtools:
```
samtools merge - movie1.fa movie2.fa | samtools sort -o merged.sort.bam
samtools index merged.sort.bam
debreak --bam merged.sort.bam
```


### Options of Debreak
#### --min_support, minimal number of supporting reads
min_support is the most important argument in filter of DeBreak. It should be adjusted according to depth of input BAM/SAM.
With no input of min_support, DeBreak estimates the depth of input BAM/SAM to give a resonable min_support.

Or you can input one value to adjust the output result.
```
debreak --bam merged.sort.bam --min_support 5
debreak --samlist sam_list --min_support 5
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

If you give coverage of dataset with --depth, DeBreak will calculate min_support according to it.
```
debreak --bam merged.sort.bam --depth 30
debreak --samlist sam_list --depth 30
```

#### --rescue_dup, rescue duplication from insertion calls
This option calls rescue_duplication module of DeBreak. minimap2, reference genome are required.
This module maps inserted sequence of each insertion call back to local region near insertion breakpoint. If the inserted sequence can be mapped, it suggests that this is a duplication instead of novel-sequence insertion.
```
debreak --bam merged.sort.bam --rescue_dup --ref reference.fa
debreak --samlist sam_list --rescue_dup --ref reference.fa
```

#### --rescue_long_ins, rescue large insertions using local de novo assembly
This option calls rescue_large_insertion module of DeBreak. minimap2, wtdbg2, and reference genome are required. Input must be sorted BAM.
DeBreak scans whole genome for insertion breakpoint candidates. It collects all reads mapped near to the candidate, and collects unmapped reads that overlap with these mapped reads.
wtdbg2 is used to generate contigs from both mapped and unmapped reads near breakpoint candidate, and minimap2 is used to map contig to reference. DeBreak then detects large insertions from the contig alignment results.
```
debreak --bam merged.sort.bam --rescue_long_ins --ref reference.fa
```


#### --poa, partial order alignment for accurate breakpoint
This option calls POA module of DeBreak. minimap2, wtdbg2, and reference genome are required.
After calling and filtering SV calls, DeBreak collects all supporting reads for each SV and calls wtdbg2 to generate a consensus sequence with less errors. It then maps consensus sequence to reference genome and detects SV breakpoint at a base resolution.
```
debreak --bam merged.sort.bam --poa --ref reference.fa
debreak --samlist --poa --ref reference.fa
```

