#!/usr/bin/python
import argparse
import Data
import sys, os

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("infile")
        parser.add_argument("-o", "--output")
        parser.add_argument("-p","--percentile", nargs=2, type=float)
        parser.add_argument("-r","--rng", nargs=2, type=float)
        args = parser.parse_args()
	x = Data.ChipData.fromFileName(args.infile)
	try:
		os.mkdir(args.output)
	except OSError:
		pass
	os.chdir(args.output)
	ranges = open("dyn_range.txt","w")
	for ix,i in enumerate(x):
		sample = i.samples
		vals = zip(i.values.transpose()[0],i.peptides)
		vals = sorted(vals)# sorted smallest to largest
		vals.reverse()# now largest to smallest
		if args.percentile:
			npeps = len(vals)
			lower = args.percentile[0]*npeps
			upper = args.percentile[1]*npeps
		elif args.rng:
			lower,upper = args.rng
		vals,peps = zip(*vals[int(lower):int(upper)])
		pepfile = open("{0}_{1}.txt".format(sample,ix),"w")
		pepfile.write("\n".join([str(i) for i in peps]))
		ranges.write("{0}\t{1}\n".format(min(vals),max(vals)))
		pepfile.close()
	open("peptides.txt","w").write("\n".join(x.peptides))
