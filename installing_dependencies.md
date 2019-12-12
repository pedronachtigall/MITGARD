A quick tutorial to install all requirements to run MITGARD.
<!---
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
Download the lastest release of Trinity [here](https://github.com/trinityrnaseq/trinityrnaseq/releases)
```
tar -xf trinityrnaseq.tar.gz
cd trinityrnaseq
make
export PATH=$PATH:path/to/trinityrnaseq/
```

SPAdes
======
Download SPAdes [here](http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz).
```
tar -xf SPAdes-3.13.0-Linux.tar.gz
export PATH=$PATH:path/to/SPAdes-3.13.0-Linux/bin/
```

MitoZ
=====
