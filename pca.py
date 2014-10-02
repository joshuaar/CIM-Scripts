#!/usr/bin/python
"""
Creates a PCA plot from a matrix and user input
Requires matplotlib, numpy, scikit learn libraries
Usage: python pca.py -c CLASSLABELS INPUTDATA
"""
import csv, argparse
from matplotlib.mlab import PCA
import numpy as np
from matplotlib import pyplot as plt
from sklearn import decomposition
import numpy as np
from matplotlib.font_manager import FontProperties


def getPoints(dataTable):
	pca = decomposition.PCA(2)
	pca.fit(dataTable)
	return pca.components_.transpose()

markers = "ov^<>8sp*hH+"
colors = "bgrcmyk"

nmarkers = len(markers)
ncolors = len(colors)

def clsIX(classVector):
	"Finds indexes for each class"
	out = {}
	for i,v in enumerate(classVector):
		lst = out.get(v,[])
		lst.append(i)
		out[v]=lst
	return out

def plotPoints(points,classVector,title,markerSize=8):
	points = np.array(points)
	classmarkers = {v: markers[i % nmarkers] for i,v in enumerate(set(classVector))}
	classcolors = {v: colors[i % ncolors] for i,v in enumerate(set(classVector))}
	clsixs = clsIX(classVector)# a dictionary providing which indexes in points to find data for each class
	for cls in sorted(clsixs.keys()):
		ixs = clsixs[cls]# indexes for the current class
		plt.plot(points[ixs,0],points[ixs,1],classmarkers[cls],markersize=markerSize,color=classcolors[cls],alpha=0.5,label=cls)
	plt.xlabel("PC1")
	plt.ylabel("PC2")
	xs = points.transpose()[0]
	ys = points.transpose()[1]
	
	plt.xlim = [min(xs),max(xs)]
	plt.ylim = [min(ys),max(ys)]
	plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
	plt.title(title)
	plt.show()
	

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
	parser.add_argument("infile", type=argparse.FileType("r"))
	parser.add_argument("-c","--cls", required=True, type=argparse.FileType("r"))
	parser.add_argument("-d","--delimiter",default="\t")
	parser.add_argument("-p","--points",default=False,action="store_true")
	parser.add_argument("-t","--title",default="PCA Plot")
	parser.add_argument("-s","--markersize",default=7,type=int)
        args = parser.parse_args()
	dataTable = [[float(j) for j in i] for i in csv.reader(args.infile,delimiter=args.delimiter)]
	dataTable = np.array(dataTable).transpose()
	pcapoints = getPoints(dataTable)
	if args.points:
		print "\n".join(["\t".join([str(j) for j in i]) for i in pcapoints])
	else:
		plotPoints(pcapoints, args.cls.read().split(), args.title, args.markersize)


