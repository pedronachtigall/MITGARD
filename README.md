![mitgard_logo](/mitgard_logo.png)

MITGARD
=======
<!---[![Latest GitHub release](https://img.shields.io/github/release/pedronachtigall/MITGARD.svg)](https://github.com/pedronachtigall/MITGARD/releases/latest) -->
<!---[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3403273.svg)](https://doi.org/10.5281/zenodo.3403273) -->
[![Published in Briefings in Bioinformatics](https://img.shields.io/badge/published%20in-Briefings%20in%20Bioinformatics-blue)](https://doi.org/10.1093/bib/bbaa429)

MITGARD (**Mit**ochondrial **G**enome **A**ssembly from **R**NA-seq **D**ata) is a computational tool designed to recover the mitochondrial genome from RNA-seq data of any Eukaryote species.

Getting Started
=================

# Installation

Download the master folder and follow the steps below:
```
unzip MITGARD-master.zip
export PATH=$PATH:path/to/MITGARD-master/bin/
```
Or git clone the MITGARD respository and add the bin folder into your PATH:
```
git clone https://github.com/pedronachtigall/MITGARD.git
export PATH=$PATH:path/to/MITGARD/bin/
```


# Requirements

- [Python](https://www.python.org/) and [Biopython](https://biopython.org/)
- [Samtools](http://quinlanlab.org/tutorials/samtools/samtools.html) (v1.9)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (v2.3.5.1)
- [Minimap2](https://github.com/lh3/minimap2) (v2.17)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) (v2.8.5)
- [SPAdes](http://cab.spbu.ru/software/spades/) (v3.13.1)
- [MitoZ](https://github.com/linzhi2013/MitoZ) (v2.4)

Ensure that all requirements are working properly.

:warning:**Installing MITGARD and dependencies - alternative 1**

If you need help on installing all requirements manually to run MITGARD, check the ["Installing_dependencies"](https://github.com/pedronachtigall/MITGARD/blob/master/installing_dependencies.md) file.

:warning:**Installing MITGARD and dependencies - alternative 2**

If you want to install MITGARD and all dependencies using [Conda environment](https://docs.conda.io/projects/conda/en/latest/index.html), follow the steps below:

1 - Ensure you have both conda-forge and bioconda channels added to the conda (it is important to add them in this order so that the priority is set correctly):
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

2 - Set the environment for MITGARD with all dependencies by following the commands below:
```
conda create -n mitgard_env samtools=1.9 bowtie2=2.3.5 minimap2=2.17 trinity=2.8.5 spades=3.13.1 libgd=2.2.4 python=3.6 biopython=1.69 ete3 perl-list-moreutils perl-params-validate perl-clone circos=0.69 perl-bioperl blast=2.2.31 hmmer=3.1b2 bwa=0.7.12 infernal=1.1.1 tbl2asn openjdk

git clone https://github.com/pedronachtigall/MITGARD.git
echo "export PATH=$PATH:$(pwd)/MITGARD/bin/" >> ~/.bash_profile

cd MITGARD/bin/
git clone https://github.com/linzhi2013/MitoZ.git
tar -jxvf MitoZ/version_2.4-alpha/release_MitoZ_v2.4-alpha.tar.bz2
echo "export PATH=$PATH:$(pwd)/release_MitoZ_v2.4-alpha/" >> ~/.bash_profile

source ~/.bash_profile

conda activate mitgard_env

python3
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
quit()
```
 - If necessary, change the permissions of all files in MITGARD/bin/: ```chmod 777 MITGARD/bin/*```
 - Then, run MITGARD with paired-end or single-end mode.
 - After, running MITGARD deactivate the conda environment ```conda deactivate```.

:warning: **Installing MITGARD and dependencies - alternative 3**

There is two scripts available to perform a semi-automated installation of MITGARD and its dependencies using Conda, just follow the commands below:
```
git clone https://github.com/pedronachtigall/MITGARD.git
cd MITGARD/
bash install_MITGARD.sh
conda activate mitgard_env
python install_NCBITaxa.py
```
 - Then, run MITGARD with paired-end or single-end mode.
 - If necessary, change the permissions of all files in MITGARD/bin/: ```chmod 777 MITGARD/bin/*```

:warning: **Installing MITGARD and dependencies - alternative 4**

Install all dependencies through conda with an ".yml" file with few manual steps. Ensure that the conda channels are added correctly (see "Installing MITGARD and dependencies - alternative 2") and follow the commands below:
```
git clone https://github.com/pedronachtigall/MITGARD.git
cd MITGARD/
conda env create -f mitgard_env.yml
conda activate mitgard_env
python install_NCBITaxa.py
cd bin/
git clone https://github.com/linzhi2013/MitoZ.git
tar -jxvf MitoZ/version_2.4-alpha/release_MitoZ_v2.4-alpha.tar.bz2
cd ..
echo "export PATH=$PATH:$(pwd)/bin/" >> ~/.bash_profile
echo "export PATH=$PATH:$(pwd)/bin/release_MitoZ_v2.4-alpha/" >> ~/.bash_profile
source ~/.bash_profile
```
 - Then, run MITGARD with paired-end or single-end mode.
 - If necessary, change the permissions of all files in MITGARD/bin/: ```chmod 777 MITGARD/bin/*```

Pipeline workflow
=================
```
|=======================================================|
|>Map reads to a reference mtDNA                        |
|    >De novo Assembly - Trinity, rnaSPAdes and MitoZ   |
|    >Genome-guided Assembly - Trinity                  |
|>Mix all contigs                                       |
|>Map contigs to reference mtDNA                        |
|>Convert contigs to mitogenome assembly                |
|=======================================================|
```

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
  -r string, --rearrangement=string
                        Optional - this parameter turn on/off an additional
                        step to check for rearrangements using the
                        mitochondrial contigs generated. To turn on set
                        "True". [Default=False]
  -g int, --genetic_code=int
                        Optional - set the genetic code to be used in the "--
                        genetic_code" parameter of MitoZ assembler [Default=2]
  -C string, --clade=string
                        Optional - set the taxonomic clade to be used in the "
                        --clade" parameter of MitoZ assembler
                        [Default=Chordata]
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

:warning:**Warning:**

- Avoid to have space and/or special characters at the reference header used to bait reads for the mitochondrial genome assembly. For instance, if the header presents this format ```>NC_010972.2 Anolis carolinensis mitochondrion, complete genome```, replace the "spaces" for "underscore" (i.e., ```>NC_010972.2_Anolis_carolinensis_mitochondrion,_complete_genome```) or remove the description after the first "space" (i.e., ```>NC_010972.2```).
- Be careful with sequences exported from bioinformatics software (e.g., Geneious) and/or opening the reference file in Windows or macOS, cause it can add special characters (such as ```\r``` at the end of each line instead of the common ```\n```, which can lead to errors during MITGARD process).

:warning:**Tips:**

**[T1]** MITGARD has two optional parameters that will be passed directly to MitoZ, which is ```--clade``` and ```--genetic_code```, which by default is set to ```--clade Chordata``` and ```--genetic_code 2```. If you are using data from any species representative of other clade, please specify it by using these parameters. Please refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more details. Specifying the best "clade" and "genetic_code" may improve the contigs generated by MitoZ, which may also improve the final mitogenome assembly generated by MITGARD.

**[T2]** MITGARD can perform an additional step to check for rearrangements. Use the options ```-r/--rearrangement True``` to perform this check step. If the user wants to run MITGARD with its default parameter and then performs the rearrangement check step later, just use the following command ```RearrangementCheck.py sample_mitogenome.fa sample_contigs.fa output_folder```.
   - After running the "RearrangementCheck" step, annotate the mitogenomes (e.g., using MITOS2 web server or the "annotate" module from Mitoz; see **Q4** in the FAQ section) and perform a manual check and correction of the gene order if needed.
   - However, this implementation still has some limitations regarding the size of contigs generated to assembly the mitogenome. It is only able to perform the analysis when at least one of the contigs presents a similar size to the reference used (i.e., >= 90%). In our tests, most of the datasets were able to reach this range of size, which allowed to check for rearrangements. If the dataset generates more fragmented contigs, the user should perform a manual check on the contigs generated (i.e., a file with the contigs is kept during MITGARD processing, ```sample_contigs.fasta```), which will help to identify rearrangements and correct the mitogenome assembly.

Output
======
By default, MITGARD outputs the mitogenome assembled and also keeps the files used to generate the final assembly. If the user performs the rearrangement check step (i.e., option ```--rearrangement True```), MITGARD will also create a directory named "RearrangementCheck".
```
MITGARD_output/
├──  sample_mitogenome.fa
├──  sample_contigs.fasta
├──  align.sam
├──  consensus.mfa.fasta
├──  mapped/
│   ├──  sample_mapped.fasta
│   ├──  sample_mapped.fq
│   └──  sample_mapped_sorted.bam
└──  RearrangementCheck/
    ├──  alignN.sam
    ├──  consensusN.mfa.fasta
    ├──  contigsN
    ├──  newrefN
    └──  newrefN_mitogenome.fa
```

Reference
=========

If you use or discuss MITGARD, please cite:

Nachtigall et al. (2021) MITGARD: an automated pipeline for mitochondrial genome assembly in eukaryotic species using RNA-seq data. Briefings in Bioinformatics. DOI:[https://doi.org/10.1093/bib/bbaa429](https://doi.org/10.1093/bib/bbaa429)


License
=======

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)

Contact
=======
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com


Frequently Asked Questions (FAQ)
================================

**[Q1]** How can I trimm and perform quality filter on my fasq files?
- We use [Trim_Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to perform quality and adapter trimming by using the command below:
```
trim_galore --paired --phred33 --length 75 -q 5 --stringency 1 -e 0.1 -o trim_galore_out Sample_R1.fastq.gz Sample_R2.fastq.gz
```

**[Q2]** Can I use merged reads instead of my paired-end reads?
- Yes, we recommend to use [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) to merge paire-end files by using the command below:
```
pear -k -j 32 -f R1.fastq.gz -r R2.fastq.gz -o output
```

**[Q3]** What OS do I need to use MITGARD?
- We tested CodAn in Ubuntu 16 and 18. But, we believe that MITGARD should work on any UNIX OS able to have all dependencies necessary to run MITGARD.

**[Q4]** How can I annotate the mitogenome assembled by MITGARD?
- There is several tools and pipelines available to annotate mitochondrial genomes. However, we recommend using the [MITOS2 webserver](http://mitos2.bioinf.uni-leipzig.de/index.py) OR using the "annotate" module from [MitoZ](https://github.com/linzhi2013/MitoZ) with the command below:
```
MitoZ.py annotate --genetic_code auto --clade Chordata --outprefix annotation_output --thread_number N --fastafile mitogenome.fa
```
**[Q5]** The script ```install_NCBITaxa.py``` is returning errors or not properly installing and/or updating the NCBITaxa properly. What to do?
- Try to update ete3 module through ```pip install --upgrade ete3```.
- If updating "ete3" module not works, it may be due to an unstable connection to NCBI. Then, download the ```taxdump.tar.gz``` file by yourself by following the instructions below.
```
wget -c http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
python
from ete3 import NCBITaxa
NCBITaxa(taxdump_file='/path/to/taxdump.tar.gz')
quit()
```
