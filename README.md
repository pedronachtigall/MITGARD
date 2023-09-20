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
- [Samtools](http://quinlanlab.org/tutorials/samtools/samtools.html) (>=v1.9)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (>=v2.3.5.1)
- [Minimap2](https://github.com/lh3/minimap2) (>=v2.17)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) (v2.8.5)
- [SPAdes](http://cab.spbu.ru/software/spades/) (v3.13.1)
- [MitoZ](https://github.com/linzhi2013/MitoZ) (v2.4) - Optional
- [Hifiasm](https://github.com/chhylp123/hifiasm) - Optional (MITGARD-LR)
- [Canu](https://github.com/marbl/canu) - Optional (MITGARD-LR)
- [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) - Optional (MITGARD-LR)

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

2 - Set the environment for MITGARD with all dependencies by following the commands below (please, change "/.bash_profile" to "/.bashrc" if it is your bash source):
```
conda create -n mitgard_env samtools=1.9 bowtie2=2.3.5 minimap2=2.17 trinity=2.8.5 spades=3.13.1 libgd=2.2.4 python=3.6 biopython ete3 perl-list-moreutils perl-params-validate perl-clone circos=0.69 perl-bioperl blast=2.2.31 hmmer=3.1b2 bwa=0.7.12 infernal=1.1.1 tbl2asn openjdk

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

Install all dependencies through conda with an ".yml" file with few manual steps. Ensure that the conda channels are added correctly (see "Installing MITGARD and dependencies - alternative 2") and follow the commands below (please, change "/.bash_profile" to "/.bashrc" if it is your bash source):
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

:warning: **Conda installation**

[![Install with conda](https://img.shields.io/badge/Install%20with-conda-success)](https://anaconda.org/bioconda/mitgard)

MITGARD can be installed with Conda by using the command: `conda install -c bioconda mitgard`

The user can also create an environment with the command: `conda create -n mitgard_env -c bioconda mitgard`. Then, activate the environment `conda activate mitgard_env` and run MITGARD with paired-end or single-end mode.

- Please, notice that the Conda installation of MITGARD does not install MitoZ. If you want to also use MitoZ in the assembling step, follow one of MITGARD's installation alternatives mentioned above or adjust the "Conda create" command line based on your expertise.

- If you are a MacOS user, notice that the Trinity assembler may present errors when installed through Conda in MacOS. In this case, check how to properly install Trinity in MacOS following Trinity's [documentation](https://github.com/trinityrnaseq/trinityrnaseq/wiki) and follow one of MITGARD's installation alternatives mentioned above by adjusting it.

- We noticed that few systems may break the Bowtie2 tool (and, consequently, MITGARD) due to a lack of the [```tbb```](https://en.wikipedia.org/wiki/Threading_Building_Blocks) library (error message: ```error while loading shared libraries: libtbb.so.2```). This issue may be simply solved by installing this library using ```sudo apt-get install libtbb-dev``` to install the library in the system or ```conda install -c conda-forge tbb``` to install the library in the activated Conda environment.

:warning: **Conda & Mamba installation**

If Conda is taking too long to solve the environment, you can take advantgae of [Mamba](https://github.com/mamba-org/mamba), which is a reimplementation of the Conda package manager in C++. It usually solves the environment quickly and helps to debug it if needed.

```
conda create -n mitgard_env mamba
conda activate mitgard_env
mamba install "python >=3.6,<3.9" biopython bowtie2 samtools trinity=2.8.5 spades=3.13.1 hifiasm minimap2 canu
git clone https://github.com/pedronachtigall/MITGARD.git
echo "export PATH=$PATH:$(pwd)/MITGARD/bin/" >> ~/.bash_profile
```
 - change ```~/.bash_profile``` to ```~/.bashrc``` if needed.

:warning: **Docker installation**

[![Docker build](https://img.shields.io/badge/Docker-build-blue)](https://hub.docker.com/repository/docker/pedronachtigall/mitgard)

If the user takes advantage of [Docker](https://docs.docker.com/) in its system, we have a pre-built Dockerfile that allows an easy build and containerization of MITGARD. Just follow the steps below:
- Git clone MITGARD repository (`git clone https://github.com/pedronachtigall/MITGARD.git`) and change to MITGARD directory (`cd MITGARD`)
- Build the container: `docker build -t mitgard:v1.0 .` (It may take a few minutes)
- In your working directory (the reads and reference files should be in there), enter in the container shell: `docker run -v $PWD:/project --rm -it mitgard:v1.0`
- Run MITGARD with paired-end or single-end mode

The user may also pull MITGARD container direct from the Docker repository following the steps below:
- Pull MITGARD container: `docker pull pedronachtigall/mitgard:latest`
- Run MITGARD container: `docker run -v $PWD:/project --rm -it pedronachtigall/mitgard:latest`
    - Please, notice that you should be in the folder containing your reads and reference files
- Run MITGARD with paired-end or single-end mode

Pipeline workflow
=================
```
|==========================================================|
|>Map reads to a reference mtDNA                           |
|    >De novo Assembly - Trinity, rnaSPAdes and/or MitoZ   |
|    >Genome-guided Assembly - Trinity                     |
|>Mix all contigs                                          |
|>Map contigs to reference mtDNA                           |
|>Convert contigs to mitogenome assembly                   |
|==========================================================|
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
  -a string, --assembler=string
                        Optional - this parameter sets the assembler tools to
                        be used to generate the mitochondrial contigs. Choice
                        one from trinity, spades and mitoz or make combination
                        between them (e.g., "trinity,spades,mitoz", or
                        "trinity,mitoz", etc). By default, MITGARD uses
                        trinity and spades, which were shown to recover the
                        entire mitochondrial genome using RNA-seq data.
                        [default=trinity,spades]
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

Check our [quick tutorial](https://github.com/pedronachtigall/MITGARD/tree/master/Tutorial) designed with simulated data to learn how to use MITGARD. Check our [SRA tutorial](https://github.com/pedronachtigall/MITGARD/blob/master/TUTORIAL_SRA.md) to use MITGARD with real dataset available at SRA.

:warning:**Warning:**

- Avoid to have space and/or special characters at the reference header used to bait reads for the mitochondrial genome assembly. For instance, if the header presents this format ```>NC_010972.2 Anolis carolinensis mitochondrion, complete genome```, replace the "spaces" for "underscore" (i.e., ```>NC_010972.2_Anolis_carolinensis_mitochondrion,_complete_genome```) or remove the description after the first "space" (i.e., ```>NC_010972.2```).
- Be careful with sequences exported from bioinformatics software (e.g., Geneious) and/or opening the reference file in Windows or macOS, cause it can add special characters (such as ```\r``` at the end of each line instead of the common ```\n```, which can lead to errors during MITGARD process).

:warning:**Tips:**

**[T1]** MITGARD has a parameter to set the assemblers used in the pipeline (i.e., `-a/--assembler` option). By default MITGARD uses Trinity and rnaSPAdes (i.e. `-a trinity,spades`), but the user can set to use any of the three assemblers that MITGARD was designed to work with (e.g., Trinity, SPAdes and/or MitoZ). For instance, if the user wants to use all three assemblers just set `-a trinity,spades,mitoz`, to use only Trinity set `-a trinity`, to use Trinity and MitoZ set `-a trinity,mitoz` and so on. Our tests revealed that Trinity and rnaSPAdes works well to recover the entire mitochondrial genome and we strongly recommend to always use both assemblers, which is the default (i.e., `-a trinity,spades`).

**[T2]** If setting MitoZ as assembler, MITGARD has two optional parameters that will be passed directly to MitoZ, which is ```--clade``` and ```--genetic_code```, which by default is set to ```--clade Chordata``` and ```--genetic_code 2```. If you are using data from any species representative of other clade, please specify it by using these parameters. Please refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for more details. Specifying the best "clade" and "genetic_code" may improve the contigs generated by MitoZ, which may also improve the final mitogenome assembly generated by MITGARD.

**[T3]** MITGARD can perform an additional step to check for rearrangements. Use the options ```-r/--rearrangement True``` to perform this check step. If the user wants to run MITGARD with its default parameter and then performs the rearrangement check step later, just use the following command ```RearrangementCheck.py sample_mitogenome.fa sample_contigs.fa output_folder```.
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

MITGARD-LR (Long-Read mode)
===========================
We designed a mode to perform mitochondrial genome assembly using **L**ong-**R**eads data (**MITGARD-LR**).

It is an experimental pipeline, which is returning satisfactory results. However, it was not thoroughly tested yet.

**Pipeline Workflow**
```
|==========================================================|
|>Map reads to a reference mtDNA                           |
|>Retrieve mitochondrial reads                             |
|    >De novo assembly - hifiasm or canu                   |
|>Map reference mtDNA to contigs                           |
|>Convert best match to mitogenome assembly                |
|==========================================================|
```

**Usage**
```
Usage: MITGARD-LR.py [options]

Options:
  -h, --help            show this help message and exit
  -s string, --sample=string
                        Mandatory - sample ID to be used in the output files
                        and final mitogenome assembly
  -r fq, --reads=fq     Mandatory - input long-reads fq file (FASTQ format),
                        /path/to/long_read.fq ; the fq file can be in .gz the
                        compressed format (e.g. long_read.fq.gz).
  -R fasta, --reference=fasta
                        Mandatory - input mitogenome in FASTA format to be
                        used as reference, /path/to/reference.fa
  -m string, --method=string
                        Optional - this parameter indicates the type of long-
                        read data being used (e.g., "pacbio_hifi",
                        "pacbio_clr", or "nanopore"). If not set, the
                        "pacbio_hifi" will be considered.
                        [default=pacbio_hifi]
  -l int, --length=int  Optional - this parameter indicates the estimated size
                        of the final mitochondiral genome (in bp; e.g., use
                        17000 for 17Kb). If not set, the estimated size will
                        be considered similar to the reference being used.
                        [default=length(reference)]
  -c int, --cpu=int     Optional - number of threads to be used in each step
                        [default=1]
```

Usage with PACBIO-HIFI data:
```
MITGARD-LR.py -s sample_id -m pacbio_hifi -r pacbio_hifi.fq(.gz) -R reference.fa
```

Usage with PACBIO-CRL data:
```
MITGARD-LR.py -s sample_id -m pacbio_clr -r pacbio_clr.fq(.gz) -R reference.fa
```

Usage with NANOPORE data:
```
MITGARD-LR.py -s sample_id -m nanopore -r nanopore.fq(.gz) -R reference.fa
```

Check our [quick tutorial](https://github.com/pedronachtigall/MITGARD/tree/master/Tutorial) to learn how to use MITGARD-LR and to test it.

**Output**
```
MITGARD_output/
├──  sample_mitogenome.fasta
├──  mitogenome_contigs.fasta ("mito_genome.p_ctg.fa" OR "mito_canu/mito_genome.contigs.fasta")
├──  mito_mapped.sam
├──  match_report.txt
└──  filtered_reads.fastq
```

 - If using `pacbio_clr` or `nanopore` data, you should consider to polish the final mitogenome (`sample_mitogenome.fasta`) to correct errors resulted from their known high error rates. You can use the mitochondrial reads (`filtered_reads.fastq`) to polish the mitogenome using tools like [Inspector](https://github.com/Maggi-Chen/Inspector), [Racon](https://github.com/isovic/racon), [BlockPolish](https://github.com/huangnengCSU/BlockPolish), or other.

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
