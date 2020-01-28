#!/usr/bin/python

#script to convert msa files to consensus in fasta format

import sys

def _MostFrequent_(List):
    counter = 0
    num = List[0]
    for i in List:
        curr_frequency = List.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i
    return num

def _ParseMSA_(fasta):
    final = []
    a = open(fasta).read()
    a1 = a.split(">")
    #to get the REF
    a2 = a1[1].split("\n")
    IDref = a2[0]
    SEQref = a2[1]
    final.append(SEQref)
    #to get the other seqs
    for n in range(2,len(a1)):
        a2 = a1[n].split("\n")
        ID = a2[0]
        SEQ = a2[1]
        final.append(SEQ)
    return final

def _GetConsensusR_(fasta):
    SEQ = ""
    a = _ParseMSA_(fasta)
    for n in range(0, len(a[0])):
        NUC = ""
        temp = []
        for m in range(1, len(a)):
            if a[m][n] == ".":
                NUC = a[0][n]
            if a[m][n] in ["A","T","C","G","-"]:
                temp.append(a[m][n])
        if len(temp) > 0:
            NUC = _MostFrequent_(temp)
            if NUC == "-":
                NUC = ""
                SEQ += NUC
            else:
                SEQ += NUC
        else:
            SEQ += NUC
    return SEQ

def _GetConsensusN_(fasta):
    SEQ = ""
    a = _ParseMSA_(fasta)
    for n in range(0, len(a[0])):
        NUC = ""
        temp = []
        for m in range(1, len(a)):
            if a[m][n] == ".":
                NUC = "N"
            if a[m][n] in ["A","T","C","G","-"]:
                temp.append(a[m][n])
        if len(temp) > 0:
            NUC = _MostFrequent_(temp)
            if NUC == "-":
                NUC = ""
                SEQ += NUC
            else:
                SEQ += NUC
        else:
            SEQ += NUC
    return SEQ

def _GenOutput_(header, fasta, output, lowcoverage):
    if lowcoverage == "N":
        SEQ = _GetConsensusN_(fasta)
    if lowcoverage == "R":
        SEQ = _GetConsensusR_(fasta)
    OUT = open(output,"w")
    OUT.write(">"+header+"\n"+"\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])+"\n")
    OUT.close()

def _main_():

    if len (sys.argv) != 5:
        print("Basic usage: msa2consensus.py header msa.fa out.fa lowcoverage")
        print("\t> header: header for the consensus sequence")
        print("\t> msa.fa: input alignment file in fasta format")
        print("\t> out.fa: output consensus sequence")
        print("\t> lowcoverage: use N, to complete with \"N\", or \"R\", to complete with the REF")
        quit()

    header = sys.argv[1]
    fasta = sys.argv[2]
    output = sys.argv[3]
    lowcoverage = sys.argv[4]
    _GenOutput_(header, fasta, output, lowcoverage)

_main_()

#END
