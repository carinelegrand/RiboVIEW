#!/usr/local/bin/python
# -*- coding: utf-8 -*-
#
# Libraries:
import numpy as np
import sys

# Test for significant eigenvalues
def TracyWidom(variances, nv, no, outfile, tol=0.1):
    # Tracy-Widom test for eigenvalues as described in: 
    # - "The distributions of random matrix theory and their applications"
    #   Tracy CA and Widom H 2008
    # - "Population Structure and Eigenanalysis", 
    #     Nick Patterson, Alkes L Price, David Reich, PLOS Genetics 2006.
    # 
    # Inputs :
    # - Variances summing to 1, corresponding to eigenvalues from a PCA or SVD
    # - number of variables (equal to number of variances)
    # - number of observations in matrix from which pca and variances were computed.
    # - name of output file
    # - optionally, tol = tolerance for sum of variances

    # -- Init, checks
    if nv != len(variances) :
        print "nv should be equal to length of variances. Please check the inputs."
        #sys.exit()
    
    if abs(1 - sum([x for x in variances])) > tol :
        print sum([x for x in variances])
        print "Eigenvalues (variances) should add up to 1. Please check the input or,"
        print " if due to limited accuracy, you can try modifying the 'tol' parameter (default is 0.1)."
        #sys.exit()
    
    TWpercentiles=np.array([[0.05,0.01,0.001],[0.9794,2.0236,3.2730]])
    TW=np.array(10*[float('Nan')])
    cList=np.array(10*[float('Nan')])
    
    # -- Calculation TW (asymptotically) distributed x for the 10 largest eigenvalues
    for i in range(nv):
        va=variances[i:]
        nVar=np.shape(va)[0]
        n,m=no,nv
        mu=((np.sqrt(n-1)+np.sqrt(m))**2) /n
        sigma=(np.sqrt(n-1)+np.sqrt(m))/np.sqrt(n) * (1/np.sqrt(n-1)+1/np.sqrt(m))**(1/3)
        L1=m*va[0]/np.sum(va)
        TW[i]=(L1-mu)/sigma
        cList[i] = 1.000
        if TW[i]>TWpercentiles[1,0] : 
            cList[i] = 0.05
        if TW[i]>TWpercentiles[1,1] : 
            cList[i] = 0.01
        if TW[i]>TWpercentiles[1,2] : 
            cList[i] = 0.001
    
    # -- Write to file
    with open(outfile, 'w') as fout :
        fout.write("\tTW\tsignif-at-level\n")
        #for row in table.T: 
        for i in range(nv) :
            fout.write( str(TW[i]) + '\t' + str(cList[i]) + '\n')
          
    """ End """

    
    """
    # -- Print TW
    print "Significance levels:    0.05     0.01     0.001 "
    print "Significance threshold: 0.9794   2.0236   3.2730"
    print "PC1 to PC10:"
    print TW
    print "Signif at 5% level ?:"
    print TW>TWpercentiles[1,0]
    print "Signif at 1% level ?:"
    print TW>TWpercentiles[1,1]
    print "Signif at 0.1% level ?:"
    print TW>TWpercentiles[1,2]
    """
