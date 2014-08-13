#!/usr/bin/python
'''Second attempt at disease classification'''
# set global variables in windows to be
# classif2.py -i 0.25 -p 0.05 -n 3 -s -o test "E:\csv"
import numpy as np
import Data
from Query import denest
from sklearn import svm
import argparse,random,csv,os
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
def sampleReshuffle(dat):
	newsamples = random.sample(dat.samples,len(dat.samples))
	return newsamples
def crossValIteration(dat,classes,cutoff,prop=0.9,reshuffle=False):
	if reshuffle:
		dat.samples = sampleReshuffle(dat)
	saved_samples = [i for i in dat.samples]
	dat.samples = ["{0}_$$_{1}".format(i,v) for i,v in enumerate(dat.samples)]
	train,test=dat.splitTraining(prop, classes)
	selectedSampleIndicies = [int(i.split("_$$_")[0]) for i in test.samples]
	dat.samples = saved_samples
	test.samples = [i.split("_$$_")[1] for i in test.samples]
	train.samples = [i.split("_$$_")[1] for i in train.samples]
	print "Training set has {0} samples from classes: {1}".format(len(train.samples),",".join(set(train.samples)))
	print "Test set has {0} samples from classes: {1}".format(len(test.samples),",".join(set(test.samples)))
	print "Selecting data..."
	# select features for each disease
	downSelected = [train.select(i) for i in classes if not i=="Normal"]#runs T test for each class except normals
	downSelectedPVals = [i.significance for i in downSelected]
	downSelectedIndex = [i.getSigIndex(cutoff) for i in downSelected] #
	print "Number of selections made for each class:"
	dsRet = downSelectedIndex
	downSelectedIndex=denest(downSelectedIndex)
	trainSel = train.filterIndex(downSelectedIndex)
	testSel = test.filterIndex(downSelectedIndex) # filter to select only desired peptides

	print "Setting up SVM..."
	Xtrain = trainSel.values.transpose()
	Ytrain = trainSel.samples
	Ytrain = label_binarize(Ytrain,classes=classes)
	n_classes = len(classes)

	clf=OneVsRestClassifier(svm.SVC(kernel='linear',probability=True))
	Xtest = testSel.values.transpose()

	Ytest = label_binarize(testSel.samples,classes=classes)

	yscore = clf.fit(Xtrain,Ytrain).decision_function(Xtest)

	fpr = dict()
	tpr = dict()
	roc_auc = dict()
	if len(classes) == 2:
		k = 1
		fpr[k],tpr[k], _ = roc_curve(Ytest,yscore)
		roc_auc[k] = auc(fpr[k],tpr[k])
	
	else:
		for i in range(len(classes)):
			fpr[i], tpr[i], _ = roc_curve(Ytest[:, i], yscore[:, i])
			roc_auc[i] = auc(fpr[i], tpr[i])
	fpr["micro"], tpr["micro"], _ = roc_curve(Ytest.ravel(), yscore.ravel())
	roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])


	plt.figure()
	plt.plot(fpr["micro"], tpr["micro"], label='micro-average (area = {0:0.2f})'.format(roc_auc["micro"]))
	print fpr.keys()
	for i in [k for k in fpr if isinstance(k,int)]:
		plt.plot(fpr[i], tpr[i], label='class {0} (area = {1:0.2f})'.format(classes[i], roc_auc[i]))
	plt.plot([0, 1], [0, 1], 'k--')
	plt.xlim([-0.01, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('')
	plt.legend(loc="lower right")
	plt.show()

def dictize(result): # result is a bunch of acc's for iterations of the script above
    out ={}
    for i in result:
        try:
            out[i[0]]
        except KeyError:
            out[i[0]]=[]
        out[i[0]].append(i[1])
    return out



def sensitivity(resDict):
    if not isinstance(resDict,dict):
        resDict = dictize(resDict)
    out = {}
    for i in resDict:
        sens = len([j for j in resDict[i] if j==i])*1.0/len(resDict[i])
        out[i] = sens
    return out
        

def specificity(resDict):
    if not isinstance(resDict,dict):
        resDict = dictize(resDict)
    out = {}
    outDenom = {}
    for i in resDict:
        out[i]=0
        outDenom[i] = 0
        for j in resDict:
            if not i == j:
                out[i]+=len([k for k in resDict[j] if not k == i])
                outDenom[i]+=len(resDict[j])
    for i in out:
        out[i]=out[i]*1./outDenom[i]
    return out

def compress(resDicts):
    '''For parsing sensitivity and specifictiy results'''
    mkLists = lambda x: {i:[x[i]]for i in x}
    resDicts = map(mkLists,resDicts)
    comp = lambda x,y: {i:x[i]+y[i] for i in x}
    return reduce(comp,resDicts)

def makeTable(x):
    "Turn a list of dicts with same keys into a matrix"
    columns = list(set([j for i in x for j in i.keys()]))
    out = []
    for row in x:
        row_build = []
        for col in columns:
            row_build.append(row[col])
        out.append(row_build)
    out = [columns] + out
    return out

def doStats(res,format=True):
    specif = map(specificity,res)
    sensid = map(sensitivity,res)
    spec = compress(specif)
    sens = compress(sensid)
    print sensid,specif
    meansd = lambda x: (np.mean(x),np.std(x))
    specsd = [meansd(spec[i]) for i in spec]
    senssd = [meansd(sens[i]) for i in sens]
    out = zip(spec.keys(),specsd,senssd)
    if not format:
        return out
    else:
        txtOut = "Disease\tSpecificity(Mean)\tSpecificity(SD)\tSensitivity(Mean)\tSensitivity(SD)\n"
        txtOut+="\n".join(["{0}\t{1:.3}\t{2:.3}\t{3:.3}\t{4:.3}".format(i[0],i[1][0],i[1][1],i[2][0],i[2][1]) for i in out])
        return txtOut,makeTable(specif),makeTable(sensid)

#expects input of peptides for each class, calculates overlap between them
def pepsInCommon(selPeps):
	f = lambda x,y: len([i for i in x if i in y])
	return [[f(j) for j in selPeps] for i in selPeps]

def pvalTable(pvals,classes):
	l = len(pvals[0])
	pv2 = [[j[i] for j in pvals] for i in range(l)] # by class (transpose)t
	print pv2
	out = {cls:dat for cls,dat in zip(classes,pv2)}
	return out
		    
if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("infile")
        parser.add_argument("-i", "--proportion",type=float,required=True)
        parser.add_argument("-p", "--p_value",type=float,required=True)
        parser.add_argument("-s", "--shuffle",action="store_true")
        args = parser.parse_args()
	print "Importing Data..."
	dat=Data.ChipData.fromFileName(args.infile) # read and normalize raw data?
	print "Setting classes..."
	classes = set(dat.samples)
	classes = [i for i in classes]
	print "Splitting data..."
	crossValIteration(dat,classes,args.p_value,args.proportion,reshuffle=args.shuffle)
