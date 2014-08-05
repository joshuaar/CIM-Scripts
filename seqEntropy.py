#!/usr/bin/python
"""
a program for computing properties of sequence spaces.
Takes a list of sequences and produces statistics as 
directed by a command line interface
"""
from collections import Counter
import sys
saved = sys.stdout
sys.stdout= open("/dev/null","w")
from Indexer import window
import itertools as it
import math, random
import argparse
sys.stdout = saved


class SeqSpace(object):
	"Abstract class for sequence space"
	def __init__(self, seqs):
		self.dist = self.counts(self.window(seqs))
		self.total = sum(zip(*self.dist.items())[1])
	def window(self, seqs):
		raise NotImplementedError( "Abstract class" )
	def counts(self, windows):
		return Counter(windows)
	def entropy(self):
		try:
			total = sum(zip(*self.dist.items())[1])
		except IndexError:
			return 0
		probs = ( (value*1.0/total)  for key,value in self.dist.items() )
		logprobs = (math.log(i,2) for i in probs)
		entropyTerms = (p * logp for p,logp in zip(probs, logprobs) )
		entropy = -sum(entropyTerms)
		return entropy
	def entropyNorm(self):
		ent=self.entropy()
		return ent / math.log(len(self.elements()),2)
	def elements(self):
		return self.dist.keys()
	def incommon(self,otherSeqSpace):
		return set(self.elements()).intersection(otherSeqSpace.elements())
	def p(self,element):
		"Gets the probability of an element under this sequence space and distribution"
		return self.dist.get(element, 1) * 1.0 / self.total
	def n(self,element):
		return self.dist[element]
	def kldiv(self, otherSeqSpace,sym=True):
		"""
		computes the KL divergence between two distributions defined under the same sequence space
		sym argument speficies if this should be forced to be symmetric (this metric normally is not)
		"""
	
		elements = set(self.elements()) # + otherSeqSpace.elements())
		switch = False
		if sym and len(set(otherSeqSpace.elements())) < len(elements):
			elements = set(otherSeqSpace.elements())
			switch = True
		if not switch:
			thiscounts = [self.n(i) for i in elements]
			othercounts = [ otherSeqSpace.n(i)+1 for i in elements ] # the counts for each element in p from q
		else:
			thiscounts = [otherSeqSpace.n(i) for i in elements]
			othercounts = [ self.n(i)+1 for i in elements ] # the counts for each element in p from q
			
		thistotal = sum(thiscounts)
		othertotal = sum(othercounts)
		thisp = ( cnt*1.0 / thistotal for cnt in thiscounts)
		otherp = ( cnt*1.0 / othertotal for cnt in othercounts)
		klterms = ( (p/q) * p for p,q in zip(thisp,otherp) )
		return sum(klterms) 
		
		

		
		
class pentamerSpace(SeqSpace):
	""" Linear space of all 5 mers """
	def window(self, seqs):
		return window(seqs, 5)	
class linearSpace(SeqSpace):
	""" Linear space of all 3 to 7 mers """
	def window(self, seqs):
		return [j for i in range(3,7) for j in window(seqs, i)]	

class residueSpace(SeqSpace):
	def window(self, seqs):
		return window(seqs, 1)

class gappedSpace(SeqSpace):
	"Gapped 4 to 7 mers"
	def window(self, seqs):
		return [j for i in seqs for j in self.getcombs(i,8,4)]
	def getcombs(self,pep,n1,n2):
		"Gets combinations where n1 is the window length and n2 is number of amino acids"
		getsmallcombs = lambda pep,n:[[(j,pep[j]) for j in i] for i in it.combinations(range(len(pep)),n)]
		if len(pep) <= n1:
			pepwin = [pep]
		else:
			pepwin = window(pep,n1)
		rawcombs = [j for i in pepwin for j in getsmallcombs(i,n2)]
		normalizecomb = lambda comb: [(i-comb[0][0],j) for i,j in comb]
		normalizedcombs = [normalizecomb(comb) for comb in rawcombs]
		joincomb = lambda comb: ".".join([str(i)+j for i,j in comb])
		return [joincomb(i) for i in normalizedcombs]

class cysteinGaps(SeqSpace):
	"Looks only at pairs of cysteins and their spacings"
	def window(self,seqs):
		return [j for i in seqs for j in self.getaagaps(i,"C")]
	
	def getaagaps(self, seq, aa):
		ixs = [i for i,v in enumerate(seq) if v == aa]
		return [ixs[i+1]-ixs[i] for i in range(len(ixs)-1)]

def checkSeqSpaceArg(arg):
	g = globals().copy()
	for name, obj in g.iteritems():
		if name == arg:
			return True
	return False

def construct(seqSpaceName):
	g = globals()
	return g[seqSpaceName]

def shuffle(seqs):
	return ["".join([random.choice(i) for j in i]) for i in seqs]
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("infile", nargs="?" ,type=argparse.FileType("r"), default=sys.stdin)
	parser.add_argument("--space", required=True)
	parser.add_argument("--dist", action="store_true")
	parser.add_argument("-s", "--shuffle", action="store_true")
	parser.add_argument("-kl", "--kldiv", type=argparse.FileType("r"))
	args = parser.parse_args()
	if not checkSeqSpaceArg(args.space):
		parser.error("Not a valid sequence space")
	peps = args.infile.read().split()
	if args.shuffle:
		sspace = construct(args.space)(shuffle(peps))
	else:
		sspace = construct(args.space)(peps)
	
	#run program
	if args.dist:
		for i,j in sspace.dist.items():
			print i,j
	elif args.kldiv:
		otherPeps = args.kldiv.read().split()
		sspace2 = construct(args.space)(otherPeps)
		common = len(sspace.incommon(sspace2))
		total = len(sspace.elements())
		print "Elements in common:\t", common
		print "Total elements:\t", total
		print "Ratio:\t", common*1.0/total
		print "KL-Divergence:\t", sspace.kldiv(sspace2)
	else:
		print "Entropy\t",sspace.entropy()
		print "Normalized Entropy\t",sspace.entropyNorm()
		

