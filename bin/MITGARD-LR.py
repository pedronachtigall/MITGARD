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
#get the size of reference to use in further analysis
def _GetRefSize_(ref):
    R = SeqIO.to_dict(SeqIO.parse(ref, 'fasta'))
    for k in R.keys():
        RefSize = len(R[k].seq)
    return RefSize

#get only mapped reads
def _ParseSAM_(sam):
    final = {}
    a = open(sam,"r")
    for line in a:
        if not line.startswith("@"):
            line1 = line.strip().split()
            readID = line1[0]
            readSEQ = line1[9]
            readQUAL = line1[10]
            final[readID] = [readSEQ, readQUAL]
    a.close()
    return final

#parse blast to generate the final mitogenome
#ref contig perc_id alength mismatch gap stref endref stcontig endcontig evalue score
def _ParseBLAST_(blast, estsize):
    mitos = {}
    bestmatch = {"mitogenome":["none", 0, 0, 0]}
    a = open(blast,"r")
    for line in a:
        line1 = line.strip().split("\t")
        contig = line1[1]
        alength = line1[3]
        st = line1[8]
        end = line1[9]
        count = 1
        if int(alength) <= int(estsize)+2000 and int(alength) >= int(estsize)-2000:
            mitos["matchregion_"+str(count)] = [contig, alength, st, end]
            count += 1
    for k in mitos.keys():
        bestsize = int(bestmatch["mitogenome"][1])
        matchsize = int(mitos[k][1])
        #the match with lower difference from the ref size is kept as the best match
        if abs(matchsize-int(estsize)) < abs(bestsize-int(estsize)):
            bestmatch["mitogenome"] = mitos[k]
            bestmatch["mitogenome"].append(k)
    a.close()
    return bestmatch, mitos

def _MITGARD_(sample, lr, reference, method, estsize, threads):
    #if hifi
    if method == "pacbio_hifi":
        #detect mito reads
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> mapping reads to reference...\n")
        os.system("minimap2 -t "+threads+" --secondary=no -ax map-hifi "+reference+" "+lr+" > mito.sam")
        os.system("samtools view -@ "+threads+" -F4 -F 0x800 -o mito_mapped.sam mito.sam")
        os.system("rm mito.sam")
        #filter reads by size (to avoid NUMTS)
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> retrieving mitochondrial reads and filtering them by size...\n")
        S = _ParseSAM_("mito_mapped.sam")
        out = open("filtered_reads.fastq", "w")
        for k in S.keys():
            if len(S[k][0]) >= int((int(estsize)+2000)/5) and len(S[k][0]) <= int(estsize)+2000:
                out.write("@"+k+"\n"+S[k][0]+"\n+\n"+S[k][1]+"\n")
        out.close()
        #assembly mito contigs
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running hifiasm assembly...\n")
        os.system("hifiasm --primary -t "+threads+" -o mito_genome filtered_reads.fastq")
        os.system("awk \'/^S/{print \">\"$2;print $3}\' mito_genome.p_ctg.gfa > mito_genome.p_ctg.fa")
        #BLAST search against reference to retrieve the final mitogenome
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> detecting the final mitogenome sequence...\n")
        os.system("makeblastdb -in mito_genome.p_ctg.fa -out blastDB/mitogenome -dbtype nucl")
        os.system("blastn -query "+reference+" -db blastDB/mitogenome -out blast.out -evalue 1E-6 -qcov_hsp_perc 40 -num_threads "+threads+" -outfmt 6")
        best, match = _ParseBLAST_("blast.out", estsize)
        #analyze blast output to check if we have a confident final mitogenome assembly
        #otherwise, print a message to state that the user must check outputs to confirm the final mitogenome
        if best["mitogenome"][0] == "none" and len(match) == 0:
            print("\n\tThe mitochondrial contigs may be too fragmented, please check the files \"mito_genome.p_ctg.fa\" and \"blast.out\". It may reveal some partial mitochondrial sequences that you can use in further analysis.")
            print("\n\tConsider to run MITGARD-LR setting one of the other methods, such as \"-m pacbio_clr\" or \"-m nanopore\". It may help to retrieve more mitochondrial reads to the assembly step due to a less stringency mapping.")
            quit()
        if best["mitogenome"][0] == "none" and len(match) > 0:
            print("\n\tMITGARD-LR was not able to retrieve the best mitogenome sequence among all possibilities. Please check the files \"mito_genome.p_ctg.fa\", \"blast.out\", and \"match_report.txt\" to retrieve the final mitogenome. It may reveal a complete or partial mitochondrial sequences that you can use in further analysis.")
            print("\n\tConsider to run MITGARD-LR setting one of the other methods, such as \"-m pacbio_clr\" or \"-m nanopore\". It may help to retrieve more mitochondrial reads to the assembly step due to a less stringency mapping.")
            MREPORT = open("match_report.txt","w")
            for k in match.keys():
                MREPORT.write(k+"\t"+"\t".join(match[k])+"\n")
            MREPORT.close()
            quit()
        if best["mitogenome"][0] != "none" and len(match) > 0:
            MREPORT = open("match_report.txt","w")
            for k in match.keys():
                MREPORT.write(k+"\t"+"\t".join(match[k])+"\n")
            MREPORT.close()
            MITOGENOME = open(sample+"_mitogenome.fasta","w")
            contigs = SeqIO.to_dict(SeqIO.parse("mito_genome.p_ctg.fa", 'fasta'))
            id = best["mitogenome"][0]
            start = best["mitogenome"][2]
            end = best["mitogenome"][3]
            matchID = best["mitogenome"][4]
            region = str(contigs[id].seq[int(start):int(end)])
            SEQ = "\n".join([region[n:n+100] for n in range(0,len(region),100)])
            ID = sample+"_mitogenome"
            MITOGENOME.write(">"+ID+"\n"+SEQ+"\n")
            MITOGENOME.close()
            print("\n\tMITGARD-LR considered the "+matchID+" in the \"match_report.txt\" as the best representative for the final mitogenome.")
            print("\n\tPlease check the final mitogenome assembly (\"SAMPLE_mitogenome.fasta\") and the \"match_report.txt\" file to ensure MITGARD-LR retrieved a confident mitogenome for your sample!!!")

    #if pb-clr or ont
    if method == "pacbio_clr" or method == "nanopore":
        #detect mito reads
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> mapping reads to reference...\n")
        os.system("minimap2 -t "+threads+" --secondary=no -ax map-pb "+reference+" "+lr+" > mito.sam")
        os.system("samtools view -@ "+threads+" -F4 -F 0x800 -o mito_mapped.sam mito.sam")
        os.system("rm mito.sam")
        #filter reads by size (to avoid NUMTS)
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> retrieving mitochondrial reads and filtering them by size...\n")
        S = _ParseSAM_("mito_mapped.sam")
        out = open("filtered_reads.fastq", "w")
        for k in S.keys():
            if len(S[k][0]) >= int((int(estsize)+2000)/5) and len(S[k][0]) <= int(estsize)+2000:
                out.write("@"+k+"\n"+S[k][0]+"\n+\n"+S[k][1]+"\n")
        out.close()
        #assembly mito contigs
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> running canu assembly...\n")
        if method == "pacbio_clr":
            os.system("canu -p mito_genome -d mito_canu genomeSize="+str(estsize)+" minInputCoverage=0.01 stopOnLowCoverage=0.01 -pacbio filtered_reads.fastq")
            #output with contigs: mito_canu/mito_genome.contigs.fasta
        if method == "nanopore":
            os.system("canu -p mito_genome -d mito_canu genomeSize="+str(estsize)+" minInputCoverage=0.01 stopOnLowCoverage=0.01 -nanopore filtered_reads.fastq")
            #output with contigs: mito_canu/mito_genome.contigs.fasta
        #BLAST search against reference to retrieve the final mitogenome
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> detecting the final mitogenome sequence...\n")
        os.system("makeblastdb -in mito_canu/mito_genome.contigs.fasta -out blastDB/mitogenome -dbtype nucl")
        os.system("blastn -query "+reference+" -db blastDB/mitogenome -out blast.out -evalue 1E-6 -qcov_hsp_perc 40 -num_threads "+threads+" -outfmt 6")
        best, match = _ParseBLAST_("blast.out", estsize)
        #analyze blast output to check if we have a confident final mitogenome assembly
        #otherwise, print a message to state that the user must check outputs to confirm the final mitogenome
        if best["mitogenome"][0] == "none" and len(match) == 0:
            print("\n\tThe mitochondrial contigs may be too fragmented, please check the files \"mito_canu/mito_genome.contigs.fasta\" and \"blast.out\". It may reveal some partial mitochondrial sequences that you can use in further analysis.")
            quit()
        if best["mitogenome"][0] == "none" and len(match) > 0:
            print("\n\tMITGARD-LR was not able to retrieve the best mitogenome sequence among all possibilities. Please check the files \"mito_canu/mito_genome.contigs.fasta\", \"blast.out\", and \"match_report.txt\" to retrieve the final mitogenome. It may reveal a complete or partial mitochondrial sequences that you can use in further analysis.")
            MREPORT = open("match_report.txt","w")
            for k in match.keys():
                MREPORT.write(k+"\t"+"\t".join(match[k])+"\n")
            MREPORT.close()
            quit()
        if best["mitogenome"][0] != "none" and len(match) > 0:
            MREPORT = open("match_report.txt","w")
            for k in match.keys():
                MREPORT.write(k+"\t"+"\t".join(match[k])+"\n")
            MREPORT.close()
            MITOGENOME = open(sample+"_mitogenome.fasta","w")
            contigs = SeqIO.to_dict(SeqIO.parse("mito_canu/mito_genome.contigs.fasta", 'fasta'))
            id = best["mitogenome"][0]
            start = best["mitogenome"][2]
            end = best["mitogenome"][3]
            matchID = best["mitogenome"][4]
            region = str(contigs[id].seq[int(start):int(end)])
            SEQ = "\n".join([region[n:n+100] for n in range(0,len(region),100)])
            ID = sample+"_mitogenome"
            MITOGENOME.write(">"+ID+"\n"+SEQ+"\n")
            MITOGENOME.close()
            print("\n\tMITGARD-LR considered the "+matchID+" in the \"match_report.txt\" as the best representative for the final mitogenome.")
            print("\n\tPlease check the final mitogenome assembly (\"SAMPLE_mitogenome.fasta\") and the \"match_report.txt\" file to ensure MITGARD-LR retrieved a confident mitogenome for your sample!!!")

##>>>>Options
def __main__():
    parser = OptionParser()
    parser.add_option("-s", "--sample", dest="sample", help="Mandatory - sample ID to be used in the output files and final mitogenome assembly", metavar="string", default=None)
    parser.add_option("-r", "--reads", dest="reads", help="Mandatory - input long-reads fq file (FASTQ format), /path/to/long_read.fq ; the fq file can be in .gz the compressed format (e.g. long_read.fq.gz).", metavar="fq", default=None)
    parser.add_option("-R", "--reference", dest="reference", help="Mandatory - input mitogenome in FASTA format to be used as reference, /path/to/reference.fa", metavar="fasta", default=None)
    parser.add_option("-m", "--method", dest="method", help="Optional - this parameter indicates the type of long-read data being used (e.g., \"pacbio_hifi\", \"pacbio_clr\", or \"nanopore\"). If not set, the \"pacbio_hifi\" will be considered. [default=pacbio_hifi]", metavar="string", default="pacbio_hifi")
    parser.add_option("-l", "--length", dest="length", help="Optional - this parameter indicates the estimated size of the final mitochondiral genome (in bp; e.g., use 17000 for 17Kb). If not set, the estimated size will be considered similar to the reference being used. [default=length(reference)]", metavar="int", default=None)
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used in each step [default=1]", metavar="int", default="1")

    (options, args) = parser.parse_args()

    if options.sample == None and options.reads == None and options.reference == None:
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

>>>> MITGARD v1.1 July 2023 <<<<
      ****Use -h for help!****

>>>> MITGARD LONG-READ MODE <<<<

USAGE - PACBIO-HIFI DATA:
MITGARD-LR.py -s sample_id -m pacbio_hifi -r pacbio_hifi.fq(.gz) -R reference.fa

USAGE - PACBIO-CRL DATA:
MITGARD-LR.py -s sample_id -m pacbio_clr -r pacbio_clr.fq(.gz) -R reference.fa

USAGE - NANOPORE DATA:
MITGARD-LR.py -s sample_id -m nanopore -r nanopore.fq(.gz) -R reference.fa

        """)
        quit()

    if options.reads != None and options.reference != None and options.sample != None:

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
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting MITGARD (v1.1 July 2023) on Long-Read Mode...")
        CWD = os.getcwd()
        print("\tLong-read file -> "+options.reads)
        print("\tReference -> "+options.reference)
        print("\tMethod -> "+options.method)
        if options.length != None:
            print("\tEstimated size -> "+options.length)
        print("\tNumber of threads -> "+options.cpu)
        print("\tThe output files will be generated at the current directory -> "+CWD+"\n")

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
            print("Error:\nThe reference used is a Multi Fasta file (i.e., more than one sequence is present in the reference file). Please, provide only one sequence in the fasta file to be used as reference (option \"-R\").")
            quit()
        if countR == 1:
            options.reference += "_REFERENCE"

        #estimating mitogenome size
        if options.length == None:
            options.length = _GetRefSize_(options.reference)
            print("\tThe estimated size was calculated using the reference ->", options.length)

        #os.mkdir("MITGARD_output_"+options.sample)
        #os.chdir("MITGARD_output_"+options.sample)

        _MITGARD_(options.sample,
                    options.reads,
                    options.reference,
                    options.method,
                    options.length,
                    str(options.cpu))

        print("\n\tFinal assembly -> "+options.sample+"_mitogenome.fa\n")

__main__()

#END