#!/usr/bin/python
"Converts a fasta file to a list of sequences, with one sequence per line"
from Bio import SeqIO
import argparse
import os, sys

def grabFasta(fastaFileName):
    if not isinstance(fastaFileName,str):
        return (j for i in fastaFileName for j in grabFasta(i))
    if not os.path.isdir(fastaFileName):
        return (i for i in SeqIO.parse(fastaFileName,"fasta"))
    else:
        fastas = os.listdir(fastaFileName)
        return (j for i in fastas for j in grabFasta(os.path.join(fastaFileName,i)))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("infile", nargs="+")
	args = parser.parse_args()
	fastas = grabFasta(args.infile)
	for i in fastas:
		print i.seq.tostring()

