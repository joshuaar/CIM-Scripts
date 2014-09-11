'''Implements feature selection methods for tabular data (array data)'''
import scipy.stats as s
import numpy as np
def tTest(list1,list2):
    f1 = len(list1) <= 1
    f2 = len(list2) <= 1
    if not f1 and not f2:# if they both have multiple values, use the 2 sample t test
        return s.ttest_ind(list1,list2)[1]
    if f1 and not f2:# if list 2 has multiple values and list 1 has single value, use 1 sample ttest
        return s.ttest_1samp(list2,list1[0])[1]
    if not f1 and f2:# list 2 has single values, list 1 has multiple values
        return s.ttest_1samp(list1,list2[0])[1]


def foldChange(list1,list2):
    a=np.mean(list1)
    b=np.mean(list2)
    return a/b
