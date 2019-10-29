#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Author : Carine Legrand
# Date : 10.04.2018
#
# Purpose :
#   Read n files containing a table mRNA/position/coverage,
#   merge to one table with each mRNA/position and corresponding
#   coverage from each file, or 0 if the mRNA/position pair is
#   absent from one or more input files
#

def unionMat(listInputFiles,outFile):
    
    verbose=False
    
    listFichiers = listInputFiles.split(';')
    
    n = len(listFichiers)
    
    if verbose : print "n = ",n
    
    # get list of idenfifier mRNA;pos
    
    if verbose : print "get list of idenfifier mRNA;pos"
    
    listIDpre=[]    
    for i in range(n) :    
        with open(listFichiers[i],'rb') as fichier:  
            for row in fichier:
                elements = row.split('\t')
                mRNA = elements[0]
                pos  = elements[1]
                ID=mRNA+";"+pos
                if mRNA != "mRNA" :
                    listIDpre.append(ID)

    listID = sorted(set(listIDpre))
    
    #with open(outFile, 'w') as fout :
    #    for ID in listID :
    #        fout.write( ID + '\n')
    
    
    # Initialize coverages at 0
    
    if verbose : print "Initialize coverages at 0"

    mRNA_pos_coverages = {}
    for ID in listID : 
        mRNA_pos_coverages[ID] = {}
        for i in range(n) :
            mRNA_pos_coverages[ID][i]=0
     
    # Update coverage wherever there is a value
    
    if verbose : print "Update coverage wherever there is a value"

    for i in range(n) :    
        with open(listFichiers[i],'rb') as fichier:  
            for row in fichier:
                elements = row.split('\t')
                mRNA = elements[0]
                pos  = elements[1]
                ID=mRNA+";"+pos
                if mRNA != "mRNA" :
                    cov = elements[2].split('\n')[0]
                    mRNA_pos_coverages[ID][i] = cov
                    
    # Write to file
    
    if verbose : print "Write to file"

    with open(outFile, 'w') as fout :
    
        title_pre = '\t'.join([ 'file'+str(i+1) for i in range(n) ])
        
        fout.write('mRNA;pos\t'+title_pre+'\n')
        
        for ID in listID :
            row_pre = '\t'.join([ str(mRNA_pos_coverages[ID][i]) for i in range(n) ])
            fout.write(ID + '\t' + row_pre + '\n')
            
    # End
    
    if verbose : print "End"
      
   












