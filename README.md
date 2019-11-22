<!--- ![mitgard_logo](/mitgard_logo.png) -->

MITGARD
=======
<!---[![Latest GitHub release](https://img.shields.io/github/release/pedronachtigall/MITGARD.svg)](https://github.com/pedronachtigall/MITGARD/releases/latest) -->
<!---[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3403273.svg)](https://doi.org/10.5281/zenodo.3403273) -->
<!---[![Published in Genome Biology](https://img.shields.io/badge/published%20in-Genome%20Biology-blue.svg)](https://doi.org/10.1101/gr.214270.116) -->

MITGARD () ...

Getting Started
=================

# Installation

Decompress the tar.gz file:
```

```


# Requirements
<!---
- [Python3](https://www.python.org/)
- [Samtools](http://quinlanlab.org/tutorials/samtools/samtools.html) (v1.9)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (v2.3.5.1)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) (v2.8.5)
- [SPAdes](http://cab.spbu.ru/software/spades/) (v3.13.1)
- [MitoZ](https://github.com/linzhi2013/MitoZ) (v2.4)
-->

Check the ["Installing the dependencies"](https://github.com/pedronachtigall/MITGARD/installing_dependencies.md) file to get help on installing all requirements to run MITGARD.

# Usage

```

```

Avoid to have space and/or special characters at the reference header used to bait reads for the mitochondrial genome assembly.


Pipeline workflow
=================

<!--- add a figure with the pipeline -->

Reference
=========

If you use or discuss MITGARD, please cite:

Nachtigall et al., under review

License
=======

<!---[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) -->

Contact
=======

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

<!---
Frequently Asked Questions (FAQ)
================================

Why the name of the software is MITGARD?
- The tool is named MITGARD, due to the similarity with [MIDGARD](https://en.wikipedia.org/wiki/Midgard). MIDGARD is the name of Earth in the Norse mythology and protected by a serpent. MITGARD was firstly designed to reconstruct the mitochondrial genome os brazilian snakes. 

How can I trimm and perform quality filter on my fasq files?
- We normally use [Trim_Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to perform quality and adapter trimming by using the following command:
```
trim_galore --paired --phred33 --length 75 -q 5 --stringency 1 -e 0.1 -o trim_galore_out Sample_R1.fastq.gz Sample_R2.fastq.gz
```

Can I use merged reads instead of my paire-end reads?
- Yes, we recommend to use [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) to merge paire-end files by using the following command:
```
pear -k -j 32 -f R1.fastq.gz -r R2.fastq.gz -o output
```


-->
