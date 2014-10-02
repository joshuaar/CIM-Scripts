"""
Computes properties of peptide collections, estimates errors. Relies on BioPython
Author: Josh Richer, 2014
"""

from Bio.SeqUtils import ProtParam as p
import numpy as np
import scipy as sp
import scipy.stats

def calcErr(data, confidence=0.95):
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
	return m, m-h, m+h
	
def calcProps(seqs):
	objs = [p.ProteinAnalysis(i) for i in seqs]
	pIs = calcErr(map(lambda x: x.isoelectric_point(),objs))#
	gravys = calcErr(map(lambda x: x.gravy(),objs))#
	aromat = calcErr(map(lambda x: x.aromaticity(),objs))#
	instabl = calcErr(map(lambda x: x.instability_index(),objs))#
	#flexi = calcErr(map(lambda x: x.flexibility(),objs))#x
	#sec = calcErr(map(lambda x: x.secondary_structure_fraction(),objs))#x
	return {"pI":pIs, "gravy":gravys, "aromaticity":aromat, "instability":instabl}




