A quick tutorial to install all requirements to run MITGARD.
<!---
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (v2.3.5.1)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) (v2.8.5)
- [SPAdes](http://cab.spbu.ru/software/spades/) (v3.13.1)
- [MitoZ](https://github.com/linzhi2013/MitoZ) (v2.4)
-->

Samtools
========
Download the lastest release of [Samtools](http://quinlanlab.org/tutorials/samtools/samtools.html)
```
git clone https://github.com/samtools/htslib
git clone https://github.com/samtools/samtools
cd samtools
make
cp samtools ~/bin
```

Bowtie2
=======
Download Bowtie2 v2.3.5.1 [here](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/) and follow the steps below:
```
unzip bowtie2-2.3.5.1-linux-x86_64.zip
export PATH=$PATH:path/to/bowtie2-2.3.5.1-linux-x86_64/
```

Trinity
=======

SPAdes
======

MitoZ
=====
