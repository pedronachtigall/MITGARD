#!/usr/bin/env python

#script designed to convert sam files to msa fasta file for downstream analysis

import	os
import	sys

g_ref_name = ""
g_ref_seq = ""
g_ref_len = 0

g_padding_seq = ""
g_count_error = 0

g_file_out = None

def	append_sequence_to_output(name, seq):

	global	g_file_out

	g_file_out.write(">" + name + "\n")
	g_file_out.write(seq + "\n")
	return

def	on_ref(name, seq):

	global	g_ref_name
	global	g_ref_seq
	global	g_ref_len
	global	g_padding_seq

	if (g_ref_len != 0):
		g_count_error += 1
		sys.stderr.write ("Error: reference has more than one sequence. Only one reference sequence is supported\n")
		return

	g_ref_name = name
	g_ref_seq = seq
	g_ref_len = len (g_ref_seq)
	append_sequence_to_output (name, seq)

	g_padding_seq = ""
	for pos in range (len(g_ref_seq)):
		g_padding_seq += "."

	return

def	sam_cigar_split(sam_cigar):

	list_digit = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
	list_alphabet = ["M", "I", "D", "N", "S", "H", "P", "=", "X"]

	list_split= []

	number_current = 0

	if (len(sam_cigar) == 0):
		return list_split

	for character in sam_cigar:
		if (character in list_digit):
			digit = int (character)
			number_current = (10 * number_current) + digit
		elif (character in list_alphabet):
			if (number_current > 0):
				list_split.append((number_current, character))
			number_current = 0

	return list_split

def	on_sam(list_sam_field):

	global	g_ref_name
	global	g_ref_len

	if len(list_sam_field) < 11:
		return

	sam_qname = list_sam_field[0]
	sam_flag = int(list_sam_field[1])
	sam_rname = list_sam_field[2]
	sam_pos = int(list_sam_field[3])
	sam_mapq = int(list_sam_field[4])
	sam_cigar = list_sam_field[5]
	sam_rnext = list_sam_field[6]
	sam_pnext = int(list_sam_field[7])
	sam_tlen = int(list_sam_field[8])
	sam_seq = list_sam_field[9]
	sam_qual = list_sam_field[10]

	if sam_pos < 1:
		# Unmapped read
		return

	if sam_rname != g_ref_name:
		print("Reference sequence name does not match")
		print("Sequence name used in sam file: " + sam_rname)
		print("Sequence name given in reference file: "  + g_ref_name)
		return

	index_seq_current = 0
	seq_current = ""

	is_first_cigar = True
	list_cigar = sam_cigar_split (sam_cigar)
	for tuple_cigar in list_cigar:
		cigar_number = tuple_cigar[0]
		cigar_character = tuple_cigar[1]

		if (cigar_character == "M"):
			# Matched
			seq_matched = sam_seq[index_seq_current:index_seq_current + cigar_number]
			seq_current += seq_matched
			index_seq_current = index_seq_current + cigar_number
		elif (cigar_character == "I"):
			# Does not support insertion
			index_seq_current = index_seq_current + cigar_number
		elif (cigar_character == "D"):
			# Deletion
			for deletion in range (cigar_number):
				seq_current += "-"
		elif (cigar_character == "N"):
			# Skipped region
			for skipped in range (cigar_number):
				seq_current += "."
		elif (cigar_character == "S"):
			# Soft clipping, clipped sequence present in seq
			if (is_first_cigar == True):
				sam_pos -= cigar_number
			seq_clip = sam_seq[index_seq_current:index_seq_current + cigar_number]
			seq_current += seq_clip
			index_seq_current = index_seq_current + cigar_number
		elif (cigar_character == "H"):
			# Hard clipping, clipped sequence not present in seq
			for clip in range (cigar_number):
				seq_current += "#"
		elif (cigar_character == "P"):
			# Padding
			print ("Error, padding", list_fields)
			g_count_error += 1
		else:
			print ("Error, unknown", list_fields)
			g_count_error += 1

		is_first_cigar = False

	# Left padding or clipping
	if (sam_pos < 1):
		# clipping
		seq_current = seq_current[1 - sam_pos:]
		sam_pos = 1

	padding_left = g_padding_seq[0:sam_pos-1]


	# Right padding
	padding_right = g_padding_seq[sam_pos - 1 + len (seq_current):g_ref_len]
	seq_padded = padding_left + seq_current + padding_right
	if (len(seq_padded) > g_ref_len):
		seq_padded = seq_padded [0:g_ref_len]
	append_sequence_to_output (sam_qname, seq_padded)
	return

def	_ParseSamFile_(file_in, callback):

	MANDATORY_NUMBER_OF_SAM_FIELDS = 11

	f = open (file_in, "r")
	for line in f:
		line = line.rstrip ("\n\r")
		if (len (line) == 0):
			continue
		if (line.startswith ("@")):
			continue
		list_sam_field = line.split ()
		if (len (list_sam_field) < MANDATORY_NUMBER_OF_SAM_FIELDS):
			continue

		callback (list_sam_field)

	f.close ()
	return

def	_ParseFastaFile_(file_in, callback):

	name = ""
	seq = ""

	f = open(file_in, "r")
	for line in f:
		line = line.rstrip("\n\r")
		if len(line) == 0:
			continue

		if line.startswith(">"):
			if seq != "":
				callback(name, seq)
			name = line[1:]
			seq = ""
		else:
			line = line.replace(" ", "")
			line = line.upper()
			seq += line
	if seq != "":
		callback(name, seq)

	f.close()
	return

def	_main_():

	global	g_file_out
	global	g_count_error

	if len(sys.argv) != 4:
		print("Basic usage: sam2msa.py REF.fa SORTED.sam MSA_OUT.fa")
		print("where:")
		print("\t ref.fasta: reference fasta file to read")
		print("\t in.sam: input sam file to read")
		print("\t out.fasta: output sequence file in aligned fasta format")
		quit()

	file_in_ref = sys.argv[1]
	file_in_sam = sys.argv[2]
	file_out_fasta = sys.argv[3]

	g_count_error = 0

	g_file_out = open(file_out_fasta, 'w')

	_ParseFastaFile_(file_in_ref, on_ref)

	if g_count_error != 0:
		sys.stderr.write("Error, in input files\n")
		quit()

	_ParseSamFile_(file_in_sam, on_sam)

	g_file_out.close()

_main_()

#END
