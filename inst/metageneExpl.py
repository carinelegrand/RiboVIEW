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

import numpy as np


def metageneExpl(inFile, res1, res2):
    
    # Parameters
    verbose = False
    t_AUG = range(-20,51)                # = [[-20 ; 50]]
    t_res1 = np.arange(-1,2,res1)        # = {-1,   -1+res1,   ..., 2-res1  }
    t_res2 = np.arange(-0.1, 0.3, res2)  # = {-0.1, -0.1+res2, ..., 0.3-res2}

    # Initialisations
    metagene          = {}
    metagene_nonSTART = 0
    metagene_START    = {}
    leakage_START     = {}
    leakage_STOP      = {}
    cov_5UTR          = {}
    cov_CDS           = {}
    cov_3UTR          = {}
    cov_5UTR_all      = 0
    cov_CDS_all       = 0
    cov_3UTR_all      = 0
    
    for elt in t_AUG :
        metagene_START[elt] = 0
    
    for t in t_res2 :
        leakage_START[t] = 0
        leakage_STOP[t] = 0

    ### Nb footprints in UTRS and CDS, as well as in small intervals along a standard transcripts
    # Reading metagene input file
    with open(inFile,'rb') as fichier:
        # Read mRNA, positions and coverage row by row
        for row in fichier:
            elements = row.split('\t')
            mRNA     = elements[0]
            
            if mRNA != "mRNA" :
                
                posNtInAsite  = int(elements[1])
                posMetagene   = float(elements[2])
                cov           = int(elements[3].split('\n')[0])
                
                
                # for inflation around AUG
                if posNtInAsite >= 0 and posNtInAsite <= 15  :
                    metagene_START[posNtInAsite] += cov
                else :
                    metagene_nonSTART += cov
                
                # for leakage at START and STOP codons
                for t in t_res2 :
                    if posMetagene > t and posMetagene <= (t+res2) :
                        leakage_START[t] += cov
                    if posMetagene > (1+t) and posMetagene <= (1+t+res2) :
                        leakage_STOP[t] += cov
                
                # Init metagene[mRNA] if not yet done
                if verbose : print "init metagene"
                if mRNA not in metagene :
                    metagene[mRNA] = {}
                    cov_5UTR[mRNA] = 0
                    cov_CDS[mRNA]  = 0
                    cov_3UTR[mRNA] = 0
                    for t in t_res1 :
                        metagene[mRNA][t] = 0
                
                # Update metagene[mRNA]
                if verbose : print "update metagene"
                for t in t_res1 :
                    if posMetagene >= t and posMetagene < (t+res1) :
                        metagene[mRNA][t] += cov
                
                # coverage in 5'UTR, CDS, 3'UTR
                if verbose : print "coverage in UTRs and CDS"
                if posMetagene < 0 :
                    cov_5UTR[mRNA] += cov
                    cov_5UTR_all   += cov
                elif posMetagene >=0 and posMetagene < 1 :
                    cov_CDS[mRNA]  += cov
                    cov_CDS_all    += cov
                elif posMetagene >=1 and posMetagene < 2 :
                    cov_3UTR[mRNA] += cov
                    cov_3UTR_all   += cov
    
    """ Exploit results by mRNA and prepare output """
    if verbose : print "End of loop, continue calculations"
    
    # initialisations
    list_nUTRdivnCDS                = []
    nb_mRNA_with_nonvoidCDS         = 0
    nb_mRNA_with_nonvoidCDS_voidUTR = 0
    metaG                           = {}
    
    # aggregate by or over mRNAs
    for mRNA in cov_CDS :
        #if cov_CDS[mRNA] > 0 :
        if cov_CDS[mRNA] > 0 :
            if verbose : print "for mRNA in cov_CDS and nonvoid CDS, mRNA : ", mRNA
            sum_metagene        = sum([metagene[mRNA][t] for t in t_res1])
            if verbose : print "metaG[mRNA], float(sum_metagene) = ", float(sum_metagene)
            metaG[mRNA]         = [ (float(metagene[mRNA][t]) / float(sum_metagene)) for t in t_res1]
            cov_CDSnonvoid_5UTR = cov_5UTR[mRNA] #for list.5.nonvoidCDS.3
            cov_CDSnonvoid_CDS  = cov_CDS[mRNA] 
            cov_CDSnonvoid_3UTR = cov_3UTR[mRNA]
            #
            if verbose : print "list_nUTRdivnCDS.append, float(cov_CDS[mRNA]) = ", float(cov_CDS[mRNA])
            list_nUTRdivnCDS.append( float(cov_5UTR[mRNA] + cov_3UTR[mRNA]) / float(cov_CDS[mRNA]) )
            #
            nb_mRNA_with_nonvoidCDS += 1
            if (cov_5UTR[mRNA] + cov_3UTR[mRNA]) == 0 :
                nb_mRNA_with_nonvoidCDS_voidUTR += 1
    
    # fraction coverage in CDS / UTRs
    if verbose : print "fraction coverage in CDS / UTRs"
    fx_0covUTR   = float(nb_mRNA_with_nonvoidCDS_voidUTR) / float(nb_mRNA_with_nonvoidCDS)
    fxCovUTR_med = np.median(list_nUTRdivnCDS)
    fq75, fq25 = np.percentile(list_nUTRdivnCDS, [75 ,25])
    fxCovUTR_iqr = fq75 - fq25
    
    
    # Metagene median and inter-quartile range
    if verbose : print "Metagene median and inter-quartile range"
    metaG_m_pre = []   #metaG.m  <- apply(metaG, 2, function(x) {median(x, na.rm = TRUE)})
    metaG_iqr_pre = []
    for i in range(len(t_res1)) :
        ll = [metaG[mRNA][i] for mRNA in metaG]
        #
        metaG_m_pre.append(np.mean(ll))    #metaG_m.append(np.median(ll)) 
        #
        q75, q25 = np.percentile(ll, [75 ,25])
        metaG_iqr_pre.append( (q75 - q25) )
    
    if verbose : print "metaG_m_pre :"
    if verbose : print metaG_m_pre
    
    metaG_sum = sum(metaG_m_pre)
    if verbose : print "sum :"
    if verbose : print metaG_sum

    metaG_m   = [float(x)/float(metaG_sum) for x in metaG_m_pre]
    metaG_iqr = [float(x)/float(metaG_sum) for x in metaG_iqr_pre]
    if verbose : print "metaG_iqr :"
    if verbose : print metaG_iqr
    if verbose : print "metaG_m :"
    if verbose : print metaG_m
    
    # Inflation around AUG 
    if verbose : print "Inflation around AUG"
    #infl <- c(infl, (metagene.start.all[[i]][21] / (sum(metagene.start.all[[i]])+ metagene.nonstart.all[[i]])))
    infl = float(metagene_START[0]) / float( sum([metagene_START[elt] for elt in metagene_START]) + metagene_nonSTART )
    
    # % Reads in UTR
    if verbose : print "% Reads in UTR"
    covUTRpc = 100. * float(cov_5UTR_all + cov_3UTR_all) / float(cov_5UTR_all + cov_CDS_all + cov_3UTR_all)
    
    # Turn dictionaries into lists
    metaG_start = [metagene_START[t] for t in t_AUG]
    leakStart   = [leakage_START[t] for t in t_res2]  #leakage.start.all[[i]]    <- leakage.START
    leakStop    = [leakage_STOP[t] for t in t_res2]   #leakage.stop.all[[i]]     <- leakage.STOP
    covUTR      = [fx_0covUTR, fxCovUTR_med, fxCovUTR_iqr, covUTRpc]
    
    ## Optional outputs :
    """
    outFile=inFile+'.cov5UTRCDS3UTR'
    with open(outFile, 'w') as ficresu :
        ficresu.write('mRNA\tcov_5UTR\tcov_CDS\tcov_3UTR\n')
        for mRNA in cov_CDS :
            ficresu.write(mRNA + '\t' + str(cov_5UTR[mRNA]) + '\t' + str(cov_CDS[mRNA]) + '\t' + str(cov_3UTR[mRNA]) + '\n')   
    """
    """
    outFile=inFile+'.metaG'
    with open(outFile, 'w') as ficresu :
        titre = '\t'.join(map(str,t_res1))
        ficresu.write('mRNA\t'+titre+'\n')
        for mRNA in metaG :
            metaGstr = '\t'.join(map(str,metaG[mRNA]))
            ficresu.write(mRNA + '\t' + metaGstr + '\n')   
    """
    """
    # for inclusion in liste :
    metaG_table = np.array([metaG[mRNA] for mRNA in metaG])
    """
    
    """ Return results """
    if verbose : print "return results"
        
    noms  = ["covUTR" "infl", "metaG_m", "metaG_iqr", "metaG_start", "leakStart", "leakStop"]
    
    liste = [covUTR, infl, metaG_m, metaG_iqr, metaG_start, leakStart, leakStop, noms]
    
    return liste

