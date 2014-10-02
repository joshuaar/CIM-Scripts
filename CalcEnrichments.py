#!/usr/bin/python
import sys

"""
To run from command line:
    Windows:
        CalcEnrichments.py librarySeqs.txt selectedSeqs.txt 3 7
"""

from collections import Counter
from scipy.stats import hypergeom,binom

def window(string,length):
    if isinstance(string,list): # If its a list of things instead of one thing
        out = map(lambda x: window(x,length), string) # window all of them and put them together
        return denest(out)
    out = []
    for i in range(len(string) - length + 1):
        out = out + [string[i:i+length]]
    return out

defaultwin = lambda peps: [k for i in peps for j in range(3,7) for k in window(i,j) ]

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

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

def calcEnrichments(library,libseqs,chosen,frm,to):
	for i in chosen:
		if not i in libseqs:
			sys.stderr.write("WARN: Sequence {0} not in library\n".format(i))
	subseqs = library.windowfunc(chosen)
	n_chosen = len(subseqs)
	subseqs = Counter(subseqs)
	res = [(i, subseqs[i], library.pepcounts[i], 1-binom.cdf(subseqs[i],n_chosen,library.pepcounts[i]*1.0/library.N_pop), subseqs[i]*1.0/library.pepcounts[i]) for i in subseqs]
	
	res_t = zip(*res)
	pvals = res_t[3]
	pvals_adj = FDR(pvals)
	res_t += [pvals_adj]
	res = zip(*res_t)
	print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format("seq","n_enriched","n_library","p_binom","frac","p_adjusted")
	for seq,n_enriched,n_lib,p_binom,frac,p_adj in sorted(res,key=lambda x:x[3]):
		if n_lib > 1 and n_enriched > 1 and not "GSG" in seq:
			print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(seq,n_enriched,n_lib,p_binom,frac,p_adj)

class PeptideLibrary(object):

	def __init__(self,peps,windowfunc=defaultwin):
		self.peps = peps
		self.windowfunc = windowfunc
		subseqs = windowfunc(peps)
		self.N_pop = len(subseqs)
		self.pepcounts = Counter(subseqs)

	def p_hypergeom(self,N_pop,n_chosen,K_pop,k_success):
		return hypergeom.cdf(k_success,N_pop,K_pop,n_chosen)

	def test_seq(self,seq,k_success,n_chosen):
		K_pop = self.pepcounts.get(seq,0)
		return self.p_hypergeom(self.N_pop,n_chosen,K_pop,k_success)

	def calc_enrichments(self,subseqs,n_chosen):
		return [(i,self.test_seq(i,subseqs[i],n_chosen)) for i in subseqs]
		
if __name__ == "__main__":
	import sys,csv
	def read(arg):
		with open(arg) as f:
			return [i[0] for i in csv.reader(f) if len(i) > 0]
	frm = int(sys.argv[3]) # first length range
	to = int(sys.argv[4]) # second length range
	winfnc = lambda peps: [k for i in peps for j in range(frm,to) for k in window(i,j) ]
	libseqs = read(sys.argv[1])
	library = PeptideLibrary(libseqs,winfnc)
	chosen = read(sys.argv[2])
	calcEnrichments(library,libseqs,chosen,frm,to)

