#!/usr/bin/python
"""
Impliments linear regression based on AA composition.
From http://www.biomedcentral.com/1471-2164/13/79
"""
import sys
saved = sys.stdout
sys.stdout = open("/dev/null","w")
import Data
from collections import Counter
import argparse, os, csv, random
import statsmodels.api as sm
import numpy as np
sys.stdout = saved

def calcAAMatrix(peptides):
	"Calculates amino acid matrix from peptide list"	
	alphabet = sorted(list(set([j for i in peptides for j in i])))
	out = []
	for peptide in peptides:
		counts = Counter(peptide)
		row = [counts.get(i,0) for i in alphabet]
		out.append(row)
	return out,alphabet

def calcLinearRegression(aaMatrix, yvals):
	"Performs OLS Linear Regression to predict intensity values from aa composition"
	ols = sm.OLS(yvals, aaMatrix)
	res = ols.fit()
	return res

def preprocess(data):
	data.peptides = data.rmEnd()
	return data

def fitModel(data, coeff=None, r2=None, residuals=None, limit=None, intercept=True):
	stdErrs = open("coeff_err.txt","w")
	allresid = []
	for cnt,i in enumerate(data):
		y = i.values[:,0]
		if limit:
			lower, upper = limit
			y,X=zip(*[k for k in zip(y,aa[0]) if k[0] < upper and k[0] > lower]) # filter on intensity limit
		else:
			X=aa[0]
		if intercept:
			X = sm.add_constant(X) # fit a y intercept
		linreg = calcLinearRegression(X,y)
		summ = open("summary_"+i.samples.replace("|","")+"_"+str(cnt)+".txt","w")
		summ.write(linreg.summary().as_text())
		summ.close() # finished writing summary
		if coeff:
			coeff.write("\t".join([str(j) for j in linreg.params])+"\n")
			stdErrs.write("\t".join([str(j) for j in linreg.HC0_se])+"\n")
		if r2:
			r2.write("{0}\t{1}".format(i.samples,linreg.rsquared)+"\n")
		if residuals:
			resid = linreg.resid
			allresid.append(resid)
	if residuals:
		allresid = np.array(allresid).transpose()
		data.values = allresid
		data.toFile(residuals)

def getInflationFactors(aa,alphabet,intercept):
	"gets variance inflation factors of X matrix"
	aa_ = np.array(aa)
	vifs = []
	for i in range(aa_.shape[1]):
		X = np.delete(aa_,i,1)
		if intercept:
			X = sm.add_constant(X)
		print X
		y = aa_[:,i].transpose()
		linreg = calcLinearRegression(X,y)
		r2 = linreg.rsquared
		vif = 1.0/(1.0 - r2)
		vifs.append(vif)
	return zip(vifs,alphabet)
			
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
        parser.add_argument("infile")
        parser.add_argument("-c", "--coeff", action="store_true")
        parser.add_argument("-a", "--alphabet", action="store_true")
        parser.add_argument("-i", "--intercept", action="store_true", default=False)
        parser.add_argument("-m", "--remove_end", action="store_true")
        parser.add_argument("-s", "--shuffle", action="store_true")
        parser.add_argument("-v", "--variance_inflation", action="store_true")
        parser.add_argument("-o", "--output")
        parser.add_argument("-t", "--log_transform", action="store_true")
        parser.add_argument("-l","--limit", nargs=2, type=float)
	parser.add_argument("--col", type=int)
        args = parser.parse_args()
	data = Data.ChipData.fromFileName(args.infile)
	if args.remove_end:
		data = preprocess(data)
	if args.log_transform:
		#data.values = data.values+1
		data.values = np.log(data.values+1)
	if args.shuffle:
		random.shuffle(data.peptides)
	if args.col:
		data = data[args.col]
	aa = calcAAMatrix(data.peptides)
	if args.variance_inflation:
		print "\n".join(["{0}\t{1}".format(*i) for i in getInflationFactors(*aa,intercept=args.intercept)])
		sys.exit(0)
	if args.alphabet:
		print "\t".join(aa[1])
		sys.exit(0)
	elif args.output:
		try:
			os.mkdir(args.output)	
		except OSError:
			pass
		os.chdir(args.output)
		open("alphabet.txt","w").write("\t".join(aa[1]))
		coeff = open("coeff.txt","w")
		r2 = open("r2.txt","w")
		residuals = "residuals.csv"
		if args.limit:
			l = args.limit
		else:
			l = None
		fitModel(data, coeff, r2, residuals,limit=l,intercept=args.intercept)
	else:
		fitModel(data, r2=sys.stdout,intercept=args.intercept)
