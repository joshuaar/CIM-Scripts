"""
Old code for manipulating 2-d labeled arrays.
It is useful and performant, but clunky and unmaintainable. 
Josh Richer, 2010
"""

import numpy as np
import re,sys
from collections import Counter
try:
    import Selection as sel
    defaultTTest = sel.tTest
except ImportError:
    sys.stderr.write("WARN: Feature selection could not be imported, do not use statistical feature selection.")
    sel = None
    defaultTTest=None
import random
import Result

def uniquify(chipData):
    peps = chipData.peptides
    vals = chipData.values
    pepdict = {}
    for i,v in enumerate(peps):
        pepdict[v] = pepdict.get(v,[]) + [i]
    out = []
    for i in pepdict:
        rows = vals[pepdict[i],]
        average = np.apply_along_axis(np.mean,0,rows)
        out.append((i,average))
    out = sorted(out)
    outpeps,outvals = zip(*out)
    outvals = np.vstack(outvals)
    chipData.values = outvals
    chipData.peptides = outpeps
       
def correlate(chipData1, chipData2):
    "Correlate all the columns in cd1 with those in cd2" 
    cd1 = chipData1.values
    cd2 = chipData2.values
    return [[(np.corrcoef(i,j)[0][1],chipData1.samples[i1],chipData2.samples[i2]) for i2,j in enumerate(cd2.transpose())] for i1,i in enumerate(cd1.transpose())]


class ChipData(object):
    def __init__(self,samples,peptide,values,rmEnd=False):
        self.samples = samples
        self.peptides = self.pepProcess(peptide)
        self.values = values
        try: # reshape 1d arrays
            self.values.shape[1]
        except IndexError:
            d1 = values.shape[0]
            d2 = 1
            self.values=self.values.reshape((d1,d2))
            
    @classmethod
    def fromFileName(cls, fileName,col=None):
        f=open(fileName)
        samples=f.readline()
        samples = samples.replace("\t",",")
        samples = samples.strip().split(",")[1:]
        if col:#Single col read
            samples = [samples[col]]
        peptides = []
        values = []
        for i in f:
            i = i.replace("\t",",")
            row = i.strip().split(",")
            peptides.append(row[0])
            row = [float(i) for i in row[1:]]
            if col:#Single col read
                row = [row[col]]
            values.append(row)
        values=np.array(values)
        f.close()
        return cls(samples,peptides,values)
    def renameCols(self,pattern,newName=None):
        if newName == None:
            newName = pattern
        cols = self.selectCols(pattern)
        notcols = self.selectCols(pattern,opposite=True)
        cols.samples = [newName for i in range(len(cols.samples))] # rename all samples to the new name
        return cols.join(notcols)
    
    def toFile(self,fileName):
        header = self.samples
        rows = [["Peptides"]+header]
        for index,row in enumerate(self.peptides):
            peptide = [row]
            self.values[index,:]
            numbers = [j for j in self.values[index,:]]
            rows.append(peptide+numbers)
        if isinstance(fileName,str):
            with open(fileName,mode="w") as f:
                csvWrite(f,rows)
        else:
            csvWrite(fileName,rows)
            
    def rmEnd(self,n=3):
        '''Find the linker of the peptide (CSG,GSG,GSC, CRH...) and remove it'''
        dctmax = lambda dct: dct[max(dct,key=lambda x:dct[x])]
        front=Counter([i[0:n] for i in self.peptides])
        back=Counter([i[-n:] for i in self.peptides])
        if dctmax(front) > dctmax(back):
            peptides = [i[n:] for i in self.peptides]
        else:
            peptides = [i[0:-n] for i in self.peptides]
        return peptides
            
    def filter(self,pepList):
        '''Downselect so only the features in peplist remain as data rows'''
        indicies,peps = zip(*[(i,v) for i,v in enumerate(self.peptides) if v in pepList])
        vals = self.values[indicies,:]
        return ChipData(self.samples,peps,vals)

    def sortOnCol(self,colNum):
        vals = self.values[:,colNum]
        peps = np.array(self.peptides)
        both = np.vstack((vals,peps,range(len(self.peptides)))).transpose()
        return both[both[:,0].argsort(),:]
    
    def getSeqsAbove(self,colNum,n):
        res = self.sortOnCol(colNum)
        vals = [i for i,v in enumerate(res[:,0]) if float(v) > n]
        return res[vals,1]

    def getTopSeqs(self,colNum,n,indicies=False):
        vals = self.values[:,colNum]
        peps = np.array(self.peptides)
        both = np.vstack((vals,peps)).transpose()
        return both[both[:,0].argsort(),1][-n:]
    
    def filterIndex(self,indicies):
        peps = [self.peptides[i] for i in indicies]
        vals = self.values[indicies,:]
        return ChipData(self.samples,peps,vals)
    
    #Just removes single quotes right now
    def pepProcess(self,peptide,rmLinker=True):
        out = ["".join(i.split("'")) for i in peptide]
        return out
    
    def select(self,group1,group2=None,f=defaultTTest,chipData=True):
        '''
        Selects 2 groups based on regular expression patterns, and applies a significance test given by function f
        Returns a chipdata object, or just the significance vector, given by the chipData flag
        If group2 is not given, then it is assumed group2 is everything not in group1
        '''
        g1 = self.selectCols(group1)
        if group2 == None:
            g2 = self.selectCols(group1,opposite=True)
        else:
            g2 = self.selectCols(group2)
        significance = [f(g1.values[i,:],g2.values[i,:]) for i in range(len(self.peptides))]
        if chipData:
            return ChipDataSelected(g1.join(g2),significance)
        else:
            return significance
    
    def selectInd(self,caseGrp,controlGrp,f=defaultTTest, chipData=True):
        '''
        Does feature selection using 1 sample tests individual (cases) v group (controls)
        '''
        cases = self.selectCols(caseGrp)
        controls = self.selectCols(controlGrp)
        out = []
        for i in cases:
            significance = [f(i.values[j],controls.values[j,:]) for j in range(len(self.peptides))]
            if chipData:
                out.append(ChipDataSelected(i,significance))
            else:
                out.append(significance)
        return ChipDataSelectedInd(out)
        
    
    def medianNormalize(self): # returns a median normalized ChipData
        newVals = self.values.copy()
        ncol = len(self.values[0,])
        for i in range(ncol):
            col = newVals[:,i]
            mCol = np.median(col)
            newVals[:,i] = col/mCol
        return ChipData(self.samples,self.peptides,newVals)

    def split(self,percent):
        '''Returns two chip data randomly selected: Split %, not Split %'''
        n = len(self.samples)
        nSel = int(percent * n)
        indicies = random.sample(range(n),nSel)
        oIndicies = [i for i in range(n) if not i in indicies]
	if len(indicies) == 0:
		indicies.append(oindicies.pop())
	elif len(oIndicies) == 0:
		oIndicies.append(indicies.pop())
	#above ensures no split has 0 elements
        
        peptides = self.peptides
        values = self.values[:,indicies]
        samples = [v for i,v in enumerate(self.samples) if i in indicies] # select the proper samples
        
        #the opposite values (!indicies)
        oPeptides = self.peptides
        oValues = self.values[:,oIndicies]
        oSamples = [v for i,v in enumerate(self.samples) if not i in indicies]
        
        return ChipData(samples,peptides,values),ChipData(oSamples,oPeptides, oValues)
    
    def splitTraining(self,percent,classes):
        '''Splits the data into training and test sets, taking a percentage of each class'''
        tr=[]
        te=[]
        for i in classes:
            tri,tei=self.selectCols(i).split(percent)
            print i,tei.samples
            tr.append(tri)
            te.append(tei)
        trOut = tr[0]
        teOut = te[0]
        for i in tr[1:]:
            trOut=trOut.join(i)
        for i in te[1:]:
            teOut=teOut.join(i)
        return trOut,teOut
	#generates a training/test split for each iteration of a leave one out scheme
	def splitLOOCV(self):
		for i in range(len(self.samples)):
			yield self[0:i].join(self[i+1:]),self[i]
			
    def selectCols(self,pattern,opposite=False): # returns columns matching a certain regex pattern
        selection = match(pattern,self.samples,opposite)
        samples = [self.samples[i] for i in selection]
        values = self.values[:,selection]
        return ChipData(samples,self.peptides,values)
    
    def mergeCols(self,pattern="",name=None,opposite=False,f=np.mean):
        selection = self.selectCols(pattern,opposite)
        if len(selection.samples)==0:
            raise Exception("Cannot merge, no samples selected")
        notSelection = self.selectCols(pattern,not opposite)
        selection.values = f(selection.values,axis=1)
        if opposite:
            pattern = "NOT_"+pattern
        if name == None:
            name = pattern
        selection = ChipData([name],selection.peptides,selection.values)
        return selection.join(notSelection)
    
    def join(self,otherData,ignore=False):
        if not self.peptides == otherData.peptides and not ignore:
            raise Exception("Peptide indicies not the same.")
        samples = self.samples + otherData.samples
        values = np.hstack((self.values, otherData.values))
        peptides = self.peptides
        return ChipData(samples,peptides,values)
    
    def innerJoin(self,otherData):
        samples = self.samples + otherData.samples
        
        otherPepDict = {v:i for i,v in enumerate(otherData.peptides)}
        selfInd = []
        otherInd = []
        for i,v in enumerate(self.peptides):
            try:
                otherInd.append(otherPepDict[v])
                selfInd.append(i)
            except KeyError:
                pass
            
        selfVals = self.values[selfInd,:]
        otherVals = otherData.values[otherInd,:]
        
        selfPeps = [self.peptides[i] for i in selfInd]
        filterPeps = [otherData.peptides[i] for i in otherInd]
        
        if not selfPeps == filterPeps:
            raise Exception("Peptide indicies not the same, this should never happen")
        values = np.hstack((selfVals,otherVals))
        peptides = selfPeps
        return ChipData(samples,peptides,values)    

    # g is for groups, i is for individuals
    # This gets the fold change for groups~not groups or individuals~individuals not in group
    def foldChange(self,groups,mode="g",f=None):
        if isinstance(groups,str):
            groups = [groups]
        values = []
        samples = []
        peptides = self.peptides
        if mode == "g":
            for i in groups:
                sel = self.selectCols(i)
                nsel = self.selectCols(i,opposite=True)
                if f==None:
                    sel = sel.mergeCols()
                    nsel = sel.mergeCols()
                    values.append(sel.values/nsel.values)
                else:
                    values.append(f(sel.values,nsel.values))
                samples.append(i)
        if mode == "i":
            for i in groups:
                sel = self.selectCols(i)
                nsel = self.selectCols(i,opposite=True)
                for j,v in enumerate(sel.samples):
                    sel_i = sel[j] # should be 1 column
                    if f==None:
                        nsel = nsel.mergeCols()
                        values.append(sel_i.values/nsel.values) # the values for sample i
                    else:
                        values.append(f(sel_i.values,nsel.values))
                    samples.append(v) # the sample name for sample i
        return ChipData(samples,peptides,np.hstack(values))

    def getMedians(self,indicies,ax=0):
        return np.median(self.values[indicies,:],ax)
    #Get # of peptides each column above a certain value
    def n_above(self,n):
        out = []
        for i in self:
            out.append(len([j for j in i.values if j > n]))
        return out
    
    def find(self,pattern):
        valid=[i for i,v in enumerate(self.peptides) if not re.match(pattern,v) == None]
        peps = [self.peptides[i] for i in valid]
        vals = self.values[valid,:]
        samp = self.samples
        return ChipData(samp,peps,vals)
    
    def tuple(self):
        '''Returns a tuple with peptide names to the right, values to the left'''
        return [(v,self.values[i]) for i,v in enumerate(self.peptides)]
    
    #returns the inverse quantile of an input q
    #Used if you want to find the intensity value corresponding to the quantile q
    def quant_inverse(self,q):
        out = []
        for i in self:
            cur=np.sort(i.values,0)
            ind = int(q*(len(cur)-1))
            out.append(cur[ind][0])
        return out
    
    def getFastas(self):
        peps = self.peptides
        out = ''
        for i in peps:
            out+=">{0}\n{1}\n".format(i,i)
        return out[0:-2]
    
    def __getitem__(self,index):
        peptides = self.peptides
        values = self.values[:,index]
        samples = self.samples[index]
        return ChipData(samples,peptides,values)

# Creates a mask of true false bits for selecting columns in data or elements in a list
def match(pattern,array,opposite=False):
    singleMatch = lambda s:not re.search(pattern,s) == None
    selection = map(singleMatch,array)
    if opposite:
        return [i for i,v in enumerate(selection) if not v]
    return [i for i,v in enumerate(selection) if v]

#Takes the results of Data.searchmotifs and finds the # of motifs in common between results
def resultSimilarity(results,percent=False):
    seqs = [[j[0] for j in i[3]] for i in results]
    out = []
    for i in seqs:
        row = [len([k for k in j if k in i]) for j in seqs]
        if percent:
            l=len(i)*1.0
            row =[round(j/l,2) for j in row]
        out.append(row)
    return out

def csvWrite(file,data):
    for i in data:
        row=[str(j) for j in i]
        row=",".join(row)+"\n"
        file.write(row)

class ChipDataSelected(ChipData):
    def __init__(self,chipData, significance=None):
        samples = chipData.samples
        peptide = chipData.peptides
        values = chipData.values
        super(ChipDataSelected,self).__init__(samples,peptide,values)
        self.significance = np.array(significance)
    
    @classmethod
    def fromRaw(cls,samples,peptide,values, significance=None):
        return cls(ChipData(samples,peptide,values),significance)
        
    def downSel(self,cutoff,lower=True):
        if cutoff == "Phil":
            cutoff = 1./len(self.peptides) # Phil significant
        if cutoff == "Bon":
            cutoff = 0.05/len(self.peptides) # Bonferroni significant
        samp = self.samples
        pep = self.peptides
        sig = self.significance
        if lower:
            cutOnes = [i for i in sig if i < cutoff]
        else:
            cutOnes = [i for i in sig if i >= cutoff]
        if len(cutOnes) == 0:
            return ChipDataSelected.fromRaw(samp,[],np.array(()),np.array(()))
        if lower:
            pep,val,sig=zip(*[(pep[i],i,sig[i]) for i in range(len(sig)) if sig[i] < cutoff])
        else:
            pep,val,sig=zip(*[(pep[i],i,sig[i]) for i in range(len(sig)) if sig[i] >= cutoff])
        print "{0} peptides selected".format(len(pep))
        val = self.values[val,:]
        return ChipDataSelected.fromRaw(samp,pep,val,sig)
    def getSigIndex(self,cutoff="Phil",lower=True):
        if cutoff == "Phil":
            cutoff = 1.0/len(self.peptides)
        if cutoff == "Bon":
            cutoff = 0.05/len(self.peptides)
        if lower:
            return [i for i,v in enumerate(self.significance) if v < cutoff]
        else:
            return [i for i,v in enumerate(self.significance) if v >= cutoff]

class ChipDataSelectedInd(object):
    '''For personalized feature selection'''
    def __init__(self,chipDataInds):
        '''Defines a collection of independant chipData objects'''
        self.chipData = chipDataInds
    def downSel(self,cutoff,lower=True):
        '''If lower is true, it keeps everything below a certain value, else it keeps everything above it'''
        return ChipDataSelectedInd([i.downSel(cutoff,lower) for i in self.chipData])
