#!/usr/bin/python
import Data
from Searcher import SubSearcher,pepMap
import pickle,os
import sys
import numpy as np
from multiprocessing import Pool, Value
"""
Performs indexing and searching for standard sequence n-mer analysis.
Very rough.
"""
"standardSeqs.py inFile fuzzy"

infile = sys.argv[1]
fuzzy = int(sys.argv[2])

if not os.path.exists("barf"):
	x = Data.ChipData.fromFileName(infile)

	sys.stderr.write("preparing peptide map\n")
	y = pepMap(x.peptides, (3,7))
	os.mkdir("barf")
	os.chdir("barf")
	
	sys.stderr.write("preparing barf directory\n")
	pickle.dump(y,open("pmap.barf","w"))
	for i,v in enumerate(x):
		v.samples = [v.samples]
		pickle.dump(v,open(str(i),"w"))
	samps = [(i,v.samples) for i,v in enumerate(x)]
	pickle.dump(samps,open("sample_list.barf","w"))
else:
	os.chdir("barf")

#parallel part
sys.stderr.write("running sequence calcs in parllel\n")
def parfunc(args):
	colNum = args[0]
	fuzzy = args[1]
	import numpy as np
	pepmap = pickle.load(open("pmap.barf"))
	data = pickle.load(open(str(colNum)))
	ss = SubSearcher(pepmap,data,rm=True)
	res = ss.searchTop(1000, np.median, fuz=fuzzy, minProb=0.001)
	pickle.dump(res,"{0}.out".format(colNum))
pool = Pool(processes=10)

samps = pickle.load(open("sample_list.barf"))

pool.map(parfunc, [(i,fuzzy) for i,v in enumerate(samps)])

sys.stderr.write("finished running sequence calcs in parllel\n")

sys.stderr.write("collecting results\n")

res = [(v[1],i,pickle.load("{0}.out".format(i))) for i,v in enumerate(samps)]

pickle.dump(open("final.out","w"),res)


