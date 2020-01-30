A quick tutorial to install all requirements to run MITGARD.

Python and Biopython
====================
```
sudo apt install python3 pip
pip install biopython ete3
```

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

Minimap2
========
Download the master folder of minimap2 at [github](https://github.com/lh3/minimap2).
```
unzip minimap2-master.zip
cd minimap2-master
make
export PATH=$PATH:path/to/minimap2-master/
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
Download the MitoZ release 2.4 [here](https://github.com/linzhi2013/MitoZ/blob/master/version_2.4-alpha/release_MitoZ_v2.4-alpha.tar.bz2).
```
tar -jxvf release_MitoZ_v2.4-alpha.tar.bz2
export PATH=$PATH:path/to/release_MitoZ_v2.4-alpha/
```
The, enter in the python shell to download a dependecy from NCBI to run MitoZ:
```
#enter in the python shell:
python3
#then follow the step below:
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
quit()
```
