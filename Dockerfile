FROM ubuntu:16.04

MAINTAINER Pedro G Nachtigall <pedronachtigall[AT]gmail[DOT]com>

ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update
RUN apt-get install -y wget bzip2 unzip git openjdk-8-jdk python3 python3-biopython gcc make cmake build-essential libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev

RUN apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# create a dir to save all programs
RUN mkdir /app

#install bowtie2
RUN wget https://ufpr.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip && unzip bowtie2-2.4.4-linux-x86_64.zip && mv bowtie2-2.4.4-linux-x86_64 /app

# install trinity
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/refs/tags/Trinity-v2.8.5.tar.gz && tar -xvzf Trinity-v2.8.5.tar.gz && cd /trinityrnaseq-Trinity-v2.8.5 && make && cd .. && mv trinityrnaseq-Trinity-v2.8.5 /app

# install spades
RUN wget http://cab.spbu.ru/files/release3.15.2/SPAdes-3.15.2-Linux.tar.gz && tar -xzf SPAdes-3.15.2-Linux.tar.gz && mv SPAdes-3.15.2-Linux /app

# install samtools
# samtools requirements: gcc make libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev
# htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && tar -vxjf htslib-1.9.tar.bz2 && mv htslib-1.9 /app && cd /app/htslib-1.9 && make && cd ../..
# BCFTools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && tar -vxjf bcftools-1.9.tar.bz2 && mv bcftools-1.9 /app && cd /app/bcftools-1.9 && make && cd ../..
# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar -vxjf samtools-1.9.tar.bz2 && mv samtools-1.9 /app && cd /app/samtools-1.9 && make

# install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.20/minimap2-2.20_x64-linux.tar.bz2 && tar -vxjf minimap2-2.20_x64-linux.tar.bz2 && mv minimap2-2.20_x64-linux /app

# install MITGARD
RUN git clone https://github.com/pedronachtigall/MITGARD.git && mv MITGARD /app && chmod +x /app/MITGARD/bin/*.py

# set PATH
ENV LC_ALL=C
ENV PATH=/app/bowtie2-2.4.4-linux-x86_64:$PATH
ENV PATH=/app/trinityrnaseq-Trinity-v2.8.5:$PATH
ENV PATH=/app/bcftools-1.9:$PATH
ENV PATH=/app/samtools-1.9:$PATH
ENV PATH=/app/htslib-1.9:$PATH
ENV PATH=/app/minimap2-2.20_x64-linux:$PATH
ENV PATH=/app/SPAdes-3.15.2-Linux/bin:$PATH
ENV PATH=/app/MITGARD/bin:$PATH

VOLUME /project
WORKDIR /project

CMD ["/bin/bash"]
