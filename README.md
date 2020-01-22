![mitgard_logo](/mitgard_logo.png)

MITGARD
=======
<!---[![Latest GitHub release](https://img.shields.io/github/release/pedronachtigall/MITGARD.svg)](https://github.com/pedronachtigall/MITGARD/releases/latest) -->
<!---[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3403273.svg)](https://doi.org/10.5281/zenodo.3403273) -->
<!---[![Published in Genome Biology](https://img.shields.io/badge/published%20in-Genome%20Biology-blue.svg)](https://doi.org/10.1101/gr.214270.116) -->

MITGARD (**Mit**ochondrial **G**enome **A**ssembly from **R**NA-seq **D**ata) is a computational tool designed to recover the mitochondrial genome from RNA-seq data of any Eukaryote species.

Getting Started
=================

# Installation

Download the master folder and follow the steps below:
```
unzip MITGARD-master.zip
export PATH=$PATH:path/to/MITGARD-master/bin/
```


# Requirements

- [Python](https://www.python.org/)
- [Samtools](http://quinlanlab.org/tutorials/samtools/samtools.html) (v1.9)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (v2.3.5.1)
- [Minimap2](https://github.com/lh3/minimap2) (v2.17)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) (v2.8.5)
- [SPAdes](http://cab.spbu.ru/software/spades/) (v3.13.1)
- [MitoZ](https://github.com/linzhi2013/MitoZ) (v2.4)

Ensure that all requirements are working properly.
If you need help on installing all requirements to run MITGARD, check the ["Installing_dependencies"](https://github.com/pedronachtigall/MITGARD/blob/master/installing_dependencies.md) file.

# Usage

```
Usage: MITGARD.py [options]

Options:
  -h, --help            show this help message and exit
  -s string, --sample=string
                        Mandatory - sample ID to be used in the output files
                        and final mitogenome assembly
  -S fq, --single_end=fq
                        Mandatory - input single-end fq file (FASTQ format),
                        /path/to/single_end.fq ; the fq file can be in .gz the
                        compressed format
  -1 fq, --paired_read1=fq
                        Mandatory - input paired-end read 1 fq file (FASTQ
                        format), /path/to/paired_read1.fq ; the fq file can be
                        in .gz the compressed format
  -2 fq, --paired_read2=fq
                        Mandatory - input paired-end read 2 fq file (FASTQ
                        format), /path/to/paired_read2.fq ; the fq file can be
                        in .gz the compressed format
  -R fasta, --reference=fasta
                        Mandatory - input mitogenome in FASTA format to be
                        used as reference, /path/to/reference.fa
  -L string, --low_coverage=string
                        Optional - this parameter decide what to do with low
                        coverage regions. Use "N" to assign N's in the low
                        coverage regions. Use "R" to assign the Reference
                        nucleotides in the low coverage regions. [default=N]
  -c int, --cpu=int     Optional - number of threads to be used in each step
                        [default=1]
  -M string, --memory=string
                        Optional - Max memory usage to be passed to Trinity
                        assembler [default=4G], use the same format as stated
                        by Trinity assembler

```

Usage in PAIRED-END mode:
```
MITGARD.py -s sample_id -1 paired_R1.fq -2 paired_R2.fq -R reference.fa
```
Usage in SINGLED-END mode:
```
MITGARD.py -s sample_id -S single_end.fq -R reference.fa
```

Check our [tutorial](https://github.com/pedronachtigall/MITGARD/blob/master/TUTORIAL.md) to learn how to use MITGARD.

**Warning:** Avoid to have space and/or special characters at the reference header used to bait reads for the mitochondrial genome assembly. For instance, if the header presents this format ```>NC_010972.2 Anolis carolinensis mitochondrion, complete genome```, replace the "spaces" for "underscore" (i.e., ```>NC_010972.2_Anolis_carolinensis_mitochondrion,_complete_genome```) or remove the description after the first "space" (i.e., ```>NC_010972.2```). Be careful with sequences exported from bioinformatics software (e.g., Geneious) and/or opening the reference file in Windows or macOS, cause it can add special characters (such as ```\r``` at the end of each line instead of the common ```\n```, which can lead to errors during MITGARD process).


Pipeline workflow
=================
```
|====================================================|
|>Map reads to a reference mtDNA                     |
|    >De novo Assembly - Trinity, SPAdes and MitoZ   |
|    >Genome-guided Assembly - Trinity               |
|>Mix all contigs                                    |
|>Map contigs to reference mtDNA                     |
|>Convert contigs to mitogenome assembly             |
|====================================================|
```

Reference
=========

If you use or discuss MITGARD, please cite:

Nachtigall et al., under review

License
=======

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)

Contact
=======
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com


Frequently Asked Questions (FAQ)
================================
<!---
Why the name of the software is MITGARD?
- The tool is named MITGARD, due to the similarity with [MIDGARD](https://en.wikipedia.org/wiki/Midgard). MIDGARD is the name of Earth in the Norse mythology and protected by a serpent. MITGARD was firstly designed to reconstruct the mitochondrial genome os brazilian snakes. 
-->

How can I trimm and perform quality filter on my fasq files?
- We use [Trim_Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to perform quality and adapter trimming by using the command below:
```
trim_galore --paired --phred33 --length 75 -q 5 --stringency 1 -e 0.1 -o trim_galore_out Sample_R1.fastq.gz Sample_R2.fastq.gz
```

Can I use merged reads instead of my paired-end reads?
- Yes, we recommend to use [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) to merge paire-end files by using the command below:
```
pear -k -j 32 -f R1.fastq.gz -r R2.fastq.gz -o output
```
