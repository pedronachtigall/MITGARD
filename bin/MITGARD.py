#!/usr/bin/env python3

#MITGARD - Mitochondrial Genome Assembly from RNA-seq Data
#Author: Pedro Gabriel Nachtigall - pedronachtigall@gmail.com

##modules
import os
import sys
import datetime as dt
from optparse import OptionParser
try:
    from Bio import SeqIO
except:
    print("Biopython module is not installed.\nInstall it using \"pip install biopython\" or check https://biopython.org for more details.")

##functions
#parse fasta files to mix paired-end into a single file
def _ParseFasta_(out, fasta1, fasta2):
    OUT = open(out, "w")
    count = 0
    for record in SeqIO.parse(fasta1, "fasta"):
        SEQ = str(record.seq)
        count += 1
        OUT.write(">seq"+str(count)+"\n"+SEQ+"\n")
    for record in SeqIO.parse(fasta2, "fasta"):
        SEQ = str(record.seq)
        count += 1
        OUT.write(">seq"+str(count)+"\n"+SEQ+"\n")
    OUT.close()

#parse fastq files to mix paired-end into a single file
def _ParseFq_(out, fq1, fq2):
    OUT = open(out, "w")
    count = 0
    for record in SeqIO.parse(fq1, "fastq"):
        SEQ = str(record.seq)
        PHRED = "G"*len(SEQ)
        count += 1
        OUT.write("@seq"+str(count)+"\n"+SEQ+"\n+\n"+PHRED+"\n")
    for record in SeqIO.parse(fq2, "fastq"):
        SEQ = str(record.seq)
        PHRED = "G"*len(SEQ)
        count += 1
        OUT.write("@seq"+str(count)+"\n"+SEQ+"\n+\n"+PHRED+"\n")
    OUT.close()

#run minimap2 on contigs and reference
def _minimap2_(sample, contigs, reference, lowcoverage):
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> generating final assembly...\n")
    os.system("minimap2 -ax splice "+reference+" "+contigs+" > align.sam")

    #convert sam to msa
    os.system("sam2msa.py "+reference+" align.sam consensus.mfa.fasta")
    #statinfo = os.stat("consensus.mfa.fasta")
    #if statinfo.st_size == 0:
    #    print("Something is wrong with the reference header, check for special characters and eliminate them from the reference file!")
    #convert msa to consensus
    os.system("msa2consensus.py "+sample+"_mitogenome consensus.mfa.fasta "+sample+"_mitogenome.fa "+lowcoverage)
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Final Assembly Ready!")

#run rearrangement check
def _RearrangementCheck_(mitogenome, contigs, outF):
    os.system("RearrangementCheck.py "+mitogenome+" "+contigs+" "+outF)

#run pipeline on single-end mode
def _SingleEnd_(sample,single,reference,N,M,lowcoverage,gc,clade,assembler):
    #index reference
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> indexing reference...\n")
    os.system("mkdir bowtie_index")
    os.system("bowtie2-build "+reference+" bowtie_index/REF")

    #mapping reads to reference
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> mapping reads to reference with LOCAL and ENDTOEND mode...\n")
    #if not os.path.exists("mapped/"):
    #    os.system("makedir mapped")
    os.system("mkdir mapped")
    bowtieLOCAL = "bowtie2 --local --threads "+N+" -x bowtie_index/REF -U "+single+" -S mapped/"+sample+"_LOCAL.sam"
    os.system(bowtieLOCAL)
    bowtieENDTOEND="bowtie2 --end-to-end --threads "+N+" -x bowtie_index/REF -U "+single+" -S mapped/"+sample+"_ENDTOEND.sam"
    os.system(bowtieENDTOEND)

    #convert SAM to BAM using samtools and delete .sam files
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> converting files to perform assemblies...\n")
    os.system("samtools view -@ "+N+" -b -S -o mapped/"+sample+"_LOCAL.bam mapped/"+sample+"_LOCAL.sam")
    os.system("rm mapped/"+sample+"_LOCAL.sam")
    os.system("samtools view -@ "+N+" -b -S -o mapped/"+sample+"_ENDTOEND.bam mapped/"+sample+"_ENDTOEND.sam")
    os.system("rm mapped/"+sample+"_ENDTOEND.sam")
    #extract only mapped reads
    os.system("samtools view -@ "+N+" -b -F 4 mapped/"+sample+"_LOCAL.bam > mapped/"+sample+"_LOCAL_mapped.bam")
    os.system("rm mapped/"+sample+"_LOCAL.bam")
    os.system("samtools view -@ "+N+" -b -F 4 mapped/"+sample+"_ENDTOEND.bam > mapped/"+sample+"_ENDTOEND_mapped.bam")
    os.system("rm mapped/"+sample+"_ENDTOEND.bam")
    #sort bam file
    os.system("samtools sort -@ "+N+" mapped/"+sample+"_LOCAL_mapped.bam -o mapped/"+sample+"_LOCAL_mapped_sorted.bam")
    os.system("rm mapped/"+sample+"_LOCAL_mapped.bam")
    os.system("samtools sort -@ "+N+" mapped/"+sample+"_ENDTOEND_mapped.bam -o mapped/"+sample+"_ENDTOEND_mapped_sorted.bam")
    os.system("rm mapped/"+sample+"_ENDTOEND_mapped.bam")
    #convert bam to fasta/q
    os.system("samtools fasta -@ "+N+" -0 mapped/"+sample+"_LOCAL_mapped.fasta mapped/"+sample+"_LOCAL_mapped_sorted.bam")
    os.system("samtools fasta -@ "+N+" -0 mapped/"+sample+"_ENDTOEND_mapped.fasta mapped/"+sample+"_ENDTOEND_mapped_sorted.bam")
    os.system("samtools fastq -@ "+N+" -0 mapped/"+sample+"_LOCAL_mapped.fq mapped/"+sample+"_LOCAL_mapped_sorted.bam")
    os.system("samtools fastq -@ "+N+" -0 mapped/"+sample+"_ENDTOEND_mapped.fq mapped/"+sample+"_ENDTOEND_mapped_sorted.bam")

    ##assemblers
    if "trinity" in assembler or "Trinity" in assembler:
        #Genome Guided Trinity assembly
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running Trinity Genome-Guided assemblies...\n")
        os.system("Trinity --genome_guided_bam mapped/"+sample+"_LOCAL_mapped_sorted.bam --max_memory "+M+" --genome_guided_max_intron 1000 --CPU "+N+" --output "+sample+"_trinityGG_LOCAL/ --full_cleanup")
        os.system("Trinity --genome_guided_bam mapped/"+sample+"_ENDTOEND_mapped_sorted.bam --max_memory "+M+" --genome_guided_max_intron 1000 --CPU "+N+" --output "+sample+"_trinityGG_ENDTOEND/ --full_cleanup")

        #De Novo Trinity assembly
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running Trinity De Novo assemblies...\n")
        os.system("Trinity --seqType fa --single mapped/"+sample+"_LOCAL_mapped.fasta --max_memory "+M+" --CPU "+N+" --output "+sample+"_trinitydenovo_LOCAL/ --full_cleanup")
        os.system("Trinity --seqType fa --single mapped/"+sample+"_ENDTOEND_mapped.fasta --max_memory "+M+" --CPU "+N+" --output "+sample+"_trinitydenovo_ENDTOEND/ --full_cleanup")

    if "spades" in assembler or "SPAdes" in assembler:
        #SPAdes assembly
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running rnaSPAdes assemblies...\n")
        os.system("rnaspades.py --threads "+N+" --phred-offset 33 -s mapped/"+sample+"_LOCAL_mapped.fq -o "+sample+"_spades_LOCAL/")
        os.system("rnaspades.py --threads "+N+" --phred-offset 33 -s mapped/"+sample+"_ENDTOEND_mapped.fq -o "+sample+"_spades_ENTOEND/")

    if "mitoz" in assembler or "MitoZ" in assembler:
        #MitoZ assembly
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running MitoZ assemblies...\n")
        os.system("MitoZ.py assemble --genetic_code "+gc+" --clade "+clade+" --thread_number "+N+" --outprefix "+sample+"_mitoz_LOCAL --fastq1 mapped/"+sample+"_LOCAL_mapped.fq")
        os.system("MitoZ.py assemble --genetic_code "+gc+" --clade "+clade+" --thread_number "+N+" --outprefix "+sample+"_mitoz_ENDTOEND --fastq1 mapped/"+sample+"_ENDTOEND_mapped.fq")

    #mix all contigs generated by each assembler
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> mixing all assemblies...\n")
    #os.system("cat "+sample+"_trinityGG_ENDTOEND/Trinity-GG.fasta "+sample+"_trinityGG_LOCAL/Trinity-GG.fasta "+sample+"_trinitydenovo_LOCAL.Trinity.fasta "+sample+"_trinitydenovo_ENDTOEND.Trinity.fasta "+sample+"_spades_LOCAL/transcripts.fasta "+sample+"_spades_ENTOEND/transcripts.fasta "+sample+"_mitoz_LOCAL.result/work71.mitogenome.fa "+sample+"_mitoz_ENDTOEND.result/work71.mitogenome.fa > "+sample+"_contigs.fasta")
    if "trinity" in assembler or "Trinity" in assembler:
        os.system("cat "+sample+"_trinityGG_ENDTOEND/Trinity-GG.fasta "+sample+"_trinityGG_LOCAL/Trinity-GG.fasta "+sample+"_trinitydenovo_LOCAL.Trinity.fasta "+sample+"_trinitydenovo_ENDTOEND.Trinity.fasta >> "+sample+"_contigs.fasta")
    if "spades" in assembler or "SPAdes" in assembler:
        os.system("cat "+sample+"_spades_LOCAL/transcripts.fasta "+sample+"_spades_ENTOEND/transcripts.fasta >> "+sample+"_contigs.fasta")
    if "mitoz" in assembler or "MitoZ" in assembler:
        os.system("cat "+sample+"_mitoz_LOCAL.result/work71.mitogenome.fa "+sample+"_mitoz_ENDTOEND.result/work71.mitogenome.fa >> "+sample+"_contigs.fasta")
    contigs = sample+"_contigs.fasta"
    contigs = sample+"_contigs.fasta"
    _minimap2_(sample, contigs, reference,lowcoverage)

#run pipeline on paired-end mode
def _PairedEnd_(sample,pairONE,pairTWO,reference,N,M,lowcoverage,gc,clade,assembler):
    #index reference
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> indexing reference...\n")
    os.system("mkdir bowtie_index")
    os.system("bowtie2-build "+reference+" bowtie_index/REF")

    #mapping reads to reference
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> mapping reads to reference with LOCAL and ENDTOEND mode...\n")
    #if not os.path.exists("mapped/"):
    #    os.system("makedir mapped")
    os.system("mkdir mapped")
    bowtieLOCAL = "bowtie2 --local --threads "+N+" -x bowtie_index/REF -1 "+pairONE+" -2 "+pairTWO+" -S mapped/"+sample+"_LOCAL.sam"
    os.system(bowtieLOCAL)
    bowtieENDTOEND="bowtie2 --end-to-end --threads "+N+" -x bowtie_index/REF --1 "+pairONE+" -2 "+pairTWO+" -S mapped/"+sample+"_ENDTOEND.sam"
    os.system(bowtieENDTOEND)

    #convert SAM to BAM using samtools and delete .sam files
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> converting files to perform assemblies...\n")
    os.system("samtools view -@ "+N+" -b -S -o mapped/"+sample+"_LOCAL.bam mapped/"+sample+"_LOCAL.sam")
    os.system("rm mapped/"+sample+"_LOCAL.sam")
    os.system("samtools view -@ "+N+" -b -S -o mapped/"+sample+"_ENDTOEND.bam mapped/"+sample+"_ENDTOEND.sam")
    os.system("rm mapped/"+sample+"_ENDTOEND.sam")
    #extract only mapped reads
    os.system("samtools view -@ "+N+" -b -F 4 mapped/"+sample+"_LOCAL.bam > mapped/"+sample+"_LOCAL_mapped.bam")
    os.system("rm mapped/"+sample+"_LOCAL.bam")
    os.system("samtools view -@ "+N+" -b -F 4 mapped/"+sample+"_ENDTOEND.bam > mapped/"+sample+"_ENDTOEND_mapped.bam")
    os.system("rm mapped/"+sample+"_ENDTOEND.bam")
    #sort bam file
    os.system("samtools sort -@ "+N+" mapped/"+sample+"_LOCAL_mapped.bam -o mapped/"+sample+"_LOCAL_mapped_sorted.bam")
    os.system("rm mapped/"+sample+"_LOCAL_mapped.bam")
    os.system("samtools sort -@ "+N+" mapped/"+sample+"_ENDTOEND_mapped.bam -o mapped/"+sample+"_ENDTOEND_mapped_sorted.bam")
    os.system("rm mapped/"+sample+"_ENDTOEND_mapped.bam")
    #convert bam to fasta/q
    os.system("samtools fasta -@ "+N+" -1 mapped/"+sample+"_LOCAL_mapped_1.fasta -2 mapped/"+sample+"_LOCAL_mapped_2.fasta mapped/"+sample+"_LOCAL_mapped_sorted.bam")
    os.system("samtools fasta -@ "+N+" -1 mapped/"+sample+"_ENDTOEND_mapped_1.fasta -2 mapped/"+sample+"_ENDTOEND_mapped_2.fasta mapped/"+sample+"_ENDTOEND_mapped_sorted.bam")
    os.system("samtools fastq -@ "+N+" -1 mapped/"+sample+"_LOCAL_mapped_1.fq -2 mapped/"+sample+"_LOCAL_mapped_2.fq mapped/"+sample+"_LOCAL_mapped_sorted.bam")
    os.system("samtools fastq -@ "+N+" -1 mapped/"+sample+"_ENDTOEND_mapped_1.fq -2 mapped/"+sample+"_ENDTOEND_mapped_2.fq mapped/"+sample+"_ENDTOEND_mapped_sorted.bam")
    #convert fasta/q from mapped reads into single-end
    _ParseFasta_("mapped/"+sample+"_LOCAL_mapped.fasta", "mapped/"+sample+"_LOCAL_mapped_1.fasta", "mapped/"+sample+"_LOCAL_mapped_2.fasta")
    _ParseFasta_("mapped/"+sample+"_ENDTOEND_mapped.fasta", "mapped/"+sample+"_ENDTOEND_mapped_1.fasta", "mapped/"+sample+"_ENDTOEND_mapped_1.fasta")
    _ParseFq_("mapped/"+sample+"_LOCAL_mapped.fq", "mapped/"+sample+"_LOCAL_mapped_1.fq", "mapped/"+sample+"_LOCAL_mapped_2.fq")
    _ParseFq_("mapped/"+sample+"_ENDTOEND_mapped.fq", "mapped/"+sample+"_ENDTOEND_mapped_1.fq", "mapped/"+sample+"_ENDTOEND_mapped_2.fq")

    ##assemblers
    if "trinity" in assembler or "Trinity" in assembler:
        #Genome Guided Trinity assembly
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running Trinity Genome-Guided assemblies...\n")
        os.system("Trinity --genome_guided_bam mapped/"+sample+"_LOCAL_mapped_sorted.bam --max_memory "+M+" --genome_guided_max_intron 1000 --CPU "+N+" --output "+sample+"_trinityGG_LOCAL/ --full_cleanup")
        os.system("Trinity --genome_guided_bam mapped/"+sample+"_ENDTOEND_mapped_sorted.bam --max_memory "+M+" --genome_guided_max_intron 1000 --CPU "+N+" --output "+sample+"_trinityGG_ENDTOEND/ --full_cleanup")

        #De Novo Trinity assembly
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running Trinity De Novo assemblies...\n")
        os.system("Trinity --seqType fa --single mapped/"+sample+"_LOCAL_mapped.fasta --max_memory "+M+" --CPU "+N+" --output "+sample+"_trinitydenovo_LOCAL/ --full_cleanup")
        os.system("Trinity --seqType fa --single mapped/"+sample+"_ENDTOEND_mapped.fasta --max_memory "+M+" --CPU "+N+" --output "+sample+"_trinitydenovo_ENDTOEND/ --full_cleanup")

    if "spades" in assembler or "SPAdes" in assembler:
        #SPAdes assembly
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running rnaSPAdes assemblies...\n")
        os.system("rnaspades.py --threads "+N+" --phred-offset 33 -s mapped/"+sample+"_LOCAL_mapped.fq -o "+sample+"_spades_LOCAL/")
        os.system("rnaspades.py --threads "+N+" --phred-offset 33 -s mapped/"+sample+"_ENDTOEND_mapped.fq -o "+sample+"_spades_ENTOEND/")

    if "mitoz" in assembler or "MitoZ" in assembler:
        #MitoZ assembly
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running MitoZ assemblies...\n")
        os.system("MitoZ.py assemble --genetic_code "+gc+" --clade "+clade+" --thread_number "+N+" --outprefix "+sample+"_mitoz_LOCAL --fastq1 mapped/"+sample+"_LOCAL_mapped.fq")
        os.system("MitoZ.py assemble --genetic_code "+gc+" --clade "+clade+" --thread_number "+N+" --outprefix "+sample+"_mitoz_ENDTOEND --fastq1 mapped/"+sample+"_ENDTOEND_mapped.fq")

    #mix all contigs generated by each assembler
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> mixing all assemblies...\n")
    #os.system("cat "+sample+"_trinityGG_ENDTOEND/Trinity-GG.fasta "+sample+"_trinityGG_LOCAL/Trinity-GG.fasta "+sample+"_trinitydenovo_LOCAL.Trinity.fasta "+sample+"_trinitydenovo_ENDTOEND.Trinity.fasta "+sample+"_spades_LOCAL/transcripts.fasta "+sample+"_spades_ENTOEND/transcripts.fasta "+sample+"_mitoz_LOCAL.result/work71.mitogenome.fa "+sample+"_mitoz_ENDTOEND.result/work71.mitogenome.fa > "+sample+"_contigs.fasta")
    if "trinity" in assembler or "Trinity" in assembler:
        os.system("cat "+sample+"_trinityGG_ENDTOEND/Trinity-GG.fasta "+sample+"_trinityGG_LOCAL/Trinity-GG.fasta "+sample+"_trinitydenovo_LOCAL.Trinity.fasta "+sample+"_trinitydenovo_ENDTOEND.Trinity.fasta >> "+sample+"_contigs.fasta")
    if "spades" in assembler or "SPAdes" in assembler:
        os.system("cat "+sample+"_spades_LOCAL/transcripts.fasta "+sample+"_spades_ENTOEND/transcripts.fasta >> "+sample+"_contigs.fasta")
    if "mitoz" in assembler or "MitoZ" in assembler:
        os.system("cat "+sample+"_mitoz_LOCAL.result/work71.mitogenome.fa "+sample+"_mitoz_ENDTOEND.result/work71.mitogenome.fa >> "+sample+"_contigs.fasta")
    contigs = sample+"_contigs.fasta"
    _minimap2_(sample, contigs, reference,lowcoverage)

##>>>>Options
def __main__():
    parser = OptionParser()
    parser.add_option("-s", "--sample", dest="sample", help="Mandatory - sample ID to be used in the output files and final mitogenome assembly", metavar="string", default=None)
    parser.add_option("-S", "--single_end", dest="single", help="Mandatory - input single-end fq file (FASTQ format), /path/to/single_end.fq ; the fq file can be in .gz the compressed format", metavar="fq", default=None)
    parser.add_option("-1", "--paired_read1", dest="pairONE", help="Mandatory - input paired-end read 1 fq file (FASTQ format), /path/to/paired_read1.fq ; the fq file can be in .gz the compressed format", metavar="fq", default=None)
    parser.add_option("-2", "--paired_read2", dest="pairTWO", help="Mandatory - input paired-end read 2 fq file (FASTQ format), /path/to/paired_read2.fq ; the fq file can be in .gz the compressed format", metavar="fq", default=None)
    parser.add_option("-R", "--reference", dest="reference", help="Mandatory - input mitogenome in FASTA format to be used as reference, /path/to/reference.fa", metavar="fasta", default=None)
    parser.add_option("-a", "--assembler", dest="assembler", help="Optional - this parameter sets the assembler tools to be used to generate the mitochondrial contigs. Choice one from trinity, spades and mitoz or make combination between them (e.g., \"trinity,spades,mitoz\", or \"trinity,mitoz\", etc). By default, MITGARD uses trinity and spades, which were shown to recover the entire mitochondrial genome using RNA-seq data. [default=trinity,spades]", metavar="string", default="trinity,spades")
    parser.add_option("-L", "--low_coverage", dest="lowcoverage", help="Optional - this parameter decide what to do with low coverage regions. Use \"N\" to assign N\'s in the low coverage regions. Use \"R\" to assign the Reference nucleotides in the low coverage regions. [default=N]", metavar="string", default="N")
    parser.add_option("-r", "--rearrangement", dest="rearrangement", help="Optional - this parameter turn on/off an additional step to check for rearrangements using the mitochondrial contigs generated. To turn on set \"True\". [Default=False]", metavar="string", default="False")
    parser.add_option("-g", "--genetic_code", dest="gc", help="Optional - set the genetic code to be used in the \"--genetic_code\" parameter of MitoZ assembler [Default=2]", metavar="int", default="2")
    parser.add_option("-C", "--clade", dest="clade", help="Optional - set the taxonomic clade to be used in the \"--clade\" parameter of MitoZ assembler [Default=Chordata]", metavar="string", default="Chordata")
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used in each step [default=1]", metavar="int", default="1")
    parser.add_option("-M", "--memory", dest="M", help="Optional - Max memory usage to be passed to Trinity assembler [default=4G], use the same format as stated by Trinity assembler", metavar="string", default="4G")

    (options, args) = parser.parse_args()

    if options.sample == None and (options.single == None or options.pairONE == None or options.pairTWO == None) and options.reference == None:
        print(
        """

####################################################
#  ___  ________ _____ _____   ___  ____________   #
#  |  \/  |_   _|_   _|  __ \ / _ \ | ___ \  _  \  #
#  | .  . | | |   | | | |  \// /_\ \| |_/ / | | |  #
#  | |\/| | | |   | | | | __ |  _  ||    /| | | |  #
#  | |  | |_| |_  | | | |_\ \| | | || |\ \| |/ /   #
#  \_|  |_/\___/  \_/  \____/\_| |_/\_| \_|___/    #
#                                                  #
####################################################

>>>> MITGARD v1.0 November 2019 <<<<
      ****Use -h for help!****

USAGE PAIRED-END MODE:
MITGARD.py -s sample_id -1 paired_R1.fq -2 paired_R2.fq -R reference.fa

USAGE SINGLED-END MODE:
MITGARD.py -s sample_id -S single_end.fq -R reference.fa

***Pay attention: the output files will be generated at the current directory!!!
***Be carefull on the number of threads and memory being used!!!
    > this will dramatically change the processing time!!!
        """)
        quit()

    if options.single != None and options.reference != None and options.sample != None:

        print("""

	####################################################
	#  ___  ________ _____ _____   ___  ____________   #
	#  |  \/  |_   _|_   _|  __ \ / _ \ | ___ \  _  \  #
	#  | .  . | | |   | | | |  \// /_\ \| |_/ / | | |  #
	#  | |\/| | | |   | | | | __ |  _  ||    /| | | |  #
	#  | |  | |_| |_  | | | |_\ \| | | || |\ \| |/ /   #
	#  \_|  |_/\___/  \_/  \____/\_| |_/\_| \_|___/    #
	#                                                  #
	####################################################

        """)
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting MITGARD (v1.0 November 2019) on Single-End Mode...")
        CWD = os.getcwd()
        print("\tSingle-end file -> "+options.single)
        print("\tReference -> "+options.reference)
        print("\tAssembler -> "+options.assembler)
        print("\tLow Coverage regions -> "+options.lowcoverage)
        print("\tNumber of threads -> "+options.cpu)
        print("\tNumber of memory -> "+options.M)
        print("\tThe output files and folders will be generated at the current directory -> "+CWD+"\n")

        #check the reference file
        countR = 0
        refOUT = open(options.reference+"_REFERENCE","w")
        for record in SeqIO.parse(options.reference, "fasta"):
            countR += 1
            seq = str(record.seq)
            SEQ = "\n".join([seq[n:n+100] for n in range(0, len(seq), 100)])
            refOUT.write(">"+str(record.id)+"\n"+SEQ+"\n")
        refOUT.close()
        if countR > 1:
            print("Error:\nThe reference used is a Multi Fasta file (i.e., more than one sequence is present in the reference file), please provide only one sequence in the fasta file to be used as reference (-option \"-R\").")
            quit()
        if countR == 1:
            options.reference += "_REFERENCE"

        #os.mkdir("MITGARD_output_"+options.sample)
        #os.chdir("MITGARD_output_"+options.sample)

        _SingleEnd_(options.sample,
                    options.single,
                    options.reference,
                    options.cpu,
                    options.M,
                    options.lowcoverage,
                    str(options.gc),
                    options.clade,
                    options.assembler)

        print("\tFinal assembly -> "+options.sample+"_mitogenome.fa")

    if options.pairONE != None and options.pairTWO != None and options.reference != None and options.sample != None:

        print("""

	####################################################
	#  ___  ________ _____ _____   ___  ____________   #
	#  |  \/  |_   _|_   _|  __ \ / _ \ | ___ \  _  \  #
	#  | .  . | | |   | | | |  \// /_\ \| |_/ / | | |  #
	#  | |\/| | | |   | | | | __ |  _  ||    /| | | |  #
	#  | |  | |_| |_  | | | |_\ \| | | || |\ \| |/ /   #
	#  \_|  |_/\___/  \_/  \____/\_| |_/\_| \_|___/    #
	#                                                  #
	####################################################

        """)
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting MITGARD (v1.0 November 2019) on Paired-End Mode...")
        CWD = os.getcwd()
        print("\tPaired-end files -> "+options.pairONE+" & "+options.pairTWO)
        print("\tReference -> "+options.reference)
        print("\tAssembler -> "+options.assembler)
        print("\tLow Coverage regions -> "+options.lowcoverage)
        print("\tNumber of threads -> "+options.cpu)
        print("\tNumber of memory -> "+options.M)
        print("\tThe output files and folders will be generated at the current directory -> "+CWD+"\n")

        #check the reference file
        countR = 0
        refOUT = open(options.reference+"_REFERENCE","w")
        for record in SeqIO.parse(options.reference, "fasta"):
            countR += 1
            seq = str(record.seq)
            SEQ = "\n".join([seq[n:n+100] for n in range(0, len(seq), 100)])
            refOUT.write(">"+str(record.id)+"\n"+SEQ+"\n")
        refOUT.close()
        if countR > 1:
            print("Error:\nThe reference used is a Multi Fasta file (i.e., more than one sequence is present in the reference file), please provide only one sequence in the fasta file to be used as reference (-option \"-R\").")
            quit()
        if countR == 1:
            options.reference += "_REFERENCE"

        #os.mkdir("MITGARD_output_"+options.sample)
        #os.chdir("MITGARD_output_"+options.sample)

        _PairedEnd_(options.sample,
                    options.pairONE,options.pairTWO,
                    options.reference,
                    options.cpu,
                    options.M,
                    options.lowcoverage,
                    str(options.gc),
                    options.clade,
                    options.assembler)

        print("\tFinal assembly -> "+options.sample+"_mitogenome.fa")

    if options.rearrangement != "False":
        print("\n\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting rearrangement check (additional step)...")
        CWD = os.getcwd()+"/"
        outF = CWD+"RearrangementCheck/"
        os.mkdir(outF)
        contigs = options.sample+"_contigs.fasta"
        mitogenome = options.sample+"_mitogenome.fa"
        _RearrangementCheck_(mitogenome, contigs, outF)
        print("\n\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Rearrangement check finished!!!\n\tCheck the results in "+CWD+"RearrangementCheck/")

__main__()

#END
