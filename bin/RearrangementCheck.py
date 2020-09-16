#!/usr/bin/env python

#script designed to check rearrangments in the mitogenome assembled by MITGARD
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import	os
import	sys
from Bio import SeqIO

def _ParseFasta_(fasta):
	final = {}
	for record in SeqIO.parse(fasta, "fasta"):
		ID = str(record.id)
		SEQ = str(record.seq)
		final[ID] = SEQ
	return final

def _GetRefLen_(ref):
	sequence = ""
	a = open(ref,"r")
	for line in a:
		if not line.startswith(">"):
			sequence += line.strip()
	return len(sequence), sequence

#run rearrangement check
def _RearrangementCheck_(mitogenome, contigs , outF):
	rl, rs = _GetRefLen_(mitogenome)
	newref = {}
	C = _ParseFasta_(contigs)
	count = 0
	for k in C.keys():
		if len(C[k]) >= rl-(rl*0.90) and len(C[k]) <= rl+(rl*0.90):
			newref["contig_"+str(count)] = C[k]
			count += 1
	if len(newref.keys()) > 0:
		os.system("mkdir "+outF)
		con = list(newref.keys())
		#print(con)
		for n in range(0, len(con)):
			NR = open(outF+"newref"+str(n), "w")
			NR.write(">"+con[n]+"\n"+newref[con[n]])
			NR.close()
			CONTIGS = open(outF+"contigs"+str(n), "w")
			for m in range(0, len(con)):
				if m != n:
					CONTIGS.write(">"+con[m]+"\n"+newref[con[m]]+"\n")
			CONTIGS.close()
			reference = outF+"newref"+str(n)
			CON = outF+"contigs"+str(n)
			os.system("minimap2 -ax splice "+reference+" "+CON+" > "+outF+"align"+str(n)+".sam")
			os.system("sam2msa.py "+reference+" "+outF+"align"+str(n)+".sam "+outF+"consensus"+str(n)+".mfa.fasta")
			os.system("msa2consensus.py "+"newref"+str(n)+"_mitogenome "+outF+"consensus"+str(n)+".mfa.fasta "+outF+"newref"+str(n)+"_mitogenome.fa N")

		print("\n\n>>>> You should annotate the mitogenome assembled and all mitogenomes generated without the reference presented at the folder \""+outF+"\". You should annotate the final mitogenome (from MITGARD) and all mitogenomes from the RearrangementCheck step (i.e., mitogenomes named as \"newrefX_mitogenome.fa\"). Then compare the gene order, if you notice rearrangements, then consider the the order presented by the \"newrefX_mitogenome.fa\". Perform a manual curation in the mitogenome assembled! If they present similar gene order, then consider the final mitogenome!")
	else:
		print("\n\n>>>> The contigs were shorter than the mitogenome used as reference/assembled. The RearrangementCheck step is not able to perform any other mitogenome assembly. You should annotate the mitogenome and all contigs to manually check if the genes are presenting a similar arrangement. If they present a similar order, then consider the final mitogenome. If they are rearranged in contigs, then you should consider the order presented in the contigs and perform a manual curation in the mitogenome assembled!")

def	_main_():

	if len(sys.argv) != 4:
		print("Basic usage: RearrangementCheck.py mitogenome.fa contigs.fa output_folder")
		print("where:")
		print("\t mitogenome.fa: mitogenome generated by MITGARD")
		print("\t contigs.fa: contigs generated by MITGARD")
		print("\t output_folder: folder to output results")
		quit()

	mitogenome = sys.argv[1]
	contigs = sys.argv[2]
	outF = sys.argv[3]

	if not outF.endswith("/"):
		outF += "/"

	_RearrangementCheck_(mitogenome, contigs , outF)

_main_()

#END