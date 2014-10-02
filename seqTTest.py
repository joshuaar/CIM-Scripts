#!/usr/bin/python
import Data,argparse,sys
import CalcEnrichments as c
import numpy as np

def FDR(x):
    """
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R  
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
    q = sum([1.0/i for i in xrange(1,len(x)+1)])
    l = [q*len(x)/i*x[j] for i,j in zip(reversed(xrange(1,len(x)+1)),o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return l

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
        parser.add_argument("infile")
        parser.add_argument("-c","--classes")
        parser.add_argument("-p","--pvalue", type=float)
	parser.add_argument("-l","--lower_pval", type=float)
	parser.add_argument("-g","--group",required=True)
	parser.add_argument("-a","--adjust",action="store_true")
	parser.add_argument("-n","--normalize",action="store_true")
	parser.add_argument("-d","--downs",action="store_true")
	parser.add_argument("-u","--ups",action="store_true")
	parser.add_argument("-e","--enrich",type=str,nargs=2)
        args = parser.parse_args()
	
	x = Data.ChipData.fromFileName(args.infile)
	if args.normalize:
		x = x.medianNormalize()
	if args.classes:
		classes = open(args.classes).read().split()
		x.samples = classes
	if args.ups or args.downs:
		x1=x.selectCols(args.group)
		x2=x.selectCols(args.group,True)
		vals1=[np.mean(i) for i in x1.values]
		vals2=[np.mean(i) for i in x2.values]
		ups = [j>i for i,j in zip(vals1,vals2)]
		if args.downs:
			ups = [j<i for i,j in zip(vals1,vals2)]
		x.values = x.values[np.array(ups),:]
		x.peptides = [i for i,j in zip(x.peptides,ups) if j]
	selected = x.select(args.group)
	pvals = selected.significance
	if args.adjust:
		pvals = FDR(pvals)
	out = zip(x.peptides,pvals)
	if args.pvalue:
		out = [i for i in out if i[1] < args.pvalue]
	if args.lower_pval:
		out = [i for i in out if i[1] > args.lower_pval]
	if args.enrich:
		frm = int(args.enrich[0])
		to = int(args.enrich[1])
		libseqs = x.peptides 
		winfnc = lambda peps: [k for i in peps for j in range(frm,to) for k in c.window(i,j)]
		library = c.PeptideLibrary(libseqs,winfnc)
		chosen = zip(*out)[0]
		c.calcEnrichments(library,libseqs,chosen,frm,to)
	else:
		for i,j in out:
			print "{0}\t{1}".format(i,j)
	
