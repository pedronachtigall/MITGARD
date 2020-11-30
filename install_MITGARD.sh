#!/bin/bash

echo "Starting MITGARD dependencies installation, it may take several minutes..."
echo "Setting conda environment..."
##add conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
##install env
#$CONDAPATH = $(conda info | grep -i 'base environment' | sed 's/       base environment : //g' | sed 's/  (writable)//g')
conda create -n mitgard_env samtools=1.9 bowtie2=2.3.5 minimap2=2.17 trinity=2.8.5 spades=3.13.1 libgd=2.2.4 python=3.6 biopython=1.69 ete3=3.0.0b35 perl-list-moreutils perl-params-validate perl-clone circos=0.69 perl-bioperl blast=2.2.31 hmmer=3.1b2 bwa=0.7.12 infernal=1.1.1 tbl2asn openjdk
conda deactivate
#source $CONDAPATH"/etc/profile.d/conda.sh"
#conda activate mitgard_env
##install ete3
#echo "Installing NCBITaxa, it may take several minutes..."
#python install_NCBITaxa.py
##install MitoZ
cd bin/
git clone https://github.com/linzhi2013/MitoZ.git
tar -jxvf MitoZ/version_2.4-alpha/release_MitoZ_v2.4-alpha.tar.bz2
cd ..
echo "Adding MITGARD and MitoZ to your PATH."
#cwd = $(pwd)
#echo $(pwd)"/bin/"
#echo $(pwd)"/bin/release_MitoZ_v2.4-alpha/"
#export PATH=$PATH:$(pwd)"/bin/"
echo "export PATH=$PATH:$(pwd)/bin/" >> ~/.bash_profile
#export PATH=$PATH:$(pwd)"/bin/release_MitoZ_v2.4-alpha/"
echo "export PATH=$PATH:$(pwd)/bin/release_MitoZ_v2.4-alpha/" >> ~/.bash_profile

source ~/.bash_profile

#echo "*************************"
#echo "Important!"
#echo "Add the path/to/MITGARD/bin/ and path/to/MITGARD/bin/release_MitoZ_v2.4-alpha/ in your ~/.bash_profile or ~/.bashrc files to have it permanently in your PATH."
#echo "If necessary, change the permissions of all files in MITGARD/bin: chmod 777 MITGARD/bin/*.py"
#echo "*************************"

echo "MITGARD isntallation is almost done. Proceed to the next steps:"
echo " - conda activate mitgard_env"
echo " - python install_NCBITaxa.py"

#END
