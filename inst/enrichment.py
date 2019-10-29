#!/usr/bin/python

### Content : 
# Functions - write_dico_double_to_file 
#           - enrichment

import pysam 
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys
import numpy as np

def write_dico_double_to_file(dico_double, keys1, keys2, outFile, suffix, opposite=False) :
  """
  
  Write a double dictionary as a labelled table in an output file
  
  Inputs :
  
    dico_double     type=dict   aligned reads (BAM file) 
    keys1           type=*      CDS annotation (output of bed2table)
    keys2           type=*      reference transcript sequences (FASTA)
    outFile         type=str    minimum read length considered
    suffix          type=str    maximum read length considered
    opposite        type=bool   if True, the opposite of keys2 (columns) is taken
  
  Outputs :
  
    direct output   : 
        none
  
    indirect output : 
        file outFile+suffix written
  
  """ 
  
  with open(outFile+suffix, 'w') as fout :
    
    # Write header containing column names (keys2)
    fout.write("\t"+"\t".join(map(str, keys2))+"\n")
    
    # Write row name (keys2) in first column, then corresponding elements
    for key1 in keys1 :
      # Define row
      if not opposite :
        row       = [dico_double[key1][key2] for key2 in keys2]
      else :
        # Note : [::-1] reverses indices
        row       = [dico_double[key1][key2] for key2 in keys2[::-1]]
      
      # Write row name, followed by values, separator = tabulation
      fout.write( str(key1) + '\t')
      for elt in row :
        fout.write( str(elt) + '\t' )
      
      fout.write( '\n' )
      
      # END of write_dico_double


def enrichment(bam,refCDS,refSeq,mini,maxi,outFile, versionStrip=False, messages=True):
  """
  Enrichment and bulk codon occupancy, rpkm, single gene's codon occupancy
  
  in input :
  
  bam     type=str  aligned reads (BAM file) 
  refCDS  type=str  CDS annotation (output of bed2table)
  refSeq  type=str  reference transcript sequences (FASTA)
  mini    type=int  minimum read length considered
  maxi    type=int  maximum read length considered
  versionStrip type=bool indicate if .1, .2, etc should be removed from transcript names                         defaults to False
  
  in output :
  
  outFile type=str  bulk codon occupancy
  
  
  Some local variables :
  
    iCodon          type=int   iCodon is in [-ex3 ; +ex5], by default [-90 ; +90]
                                 A codon at iCodon offset of current position is
                                 at (pos - iCodon) since "pos" is in 5'-3' direction
                                 while "iCodon" runs from 3' to 5' (by default -90
                                 at 3' furthest end and +90 at 5' furthest end).
  
  """ 
  
  #Parameters
  offset  = 1 #allows 1 clutching base at start of read
  verbose = False
  ex3 = 90
  ex5 = 90
  
  if verbose : print "Read CDS annotation file (tab-separated table)"
  
  #################################
  # Database of transcript's structure : local start, local stop, strand
  # as dictionary containing gene ID as key with local start and end sites as tuple.
  
  localDic = {}
  
  with open(refCDS, mode='r') as f :
    header=f.readline()
    
    for line in f:
      elements=line.strip('\n').split('\t')
      if verbose : print elements 
      #localDic[ID] =  (localStart, localEnd, strandDir, exonLength)
      if not elements[0] in localDic :
        localDic[elements[0]] = (int(elements[1]), int(elements[2]), elements[3], int(elements[4]))
      elif messages :
        print "Annotation for mRNA transcript ",elements[0],"was found several times in ",refCDS,"."
        print " -> Only the first occurence is used here. "
        print " -> Please edit ",refCDS," if you wish to select another occurence."
  
  
  if verbose : print "Read alignment file"
  
  #samfile = pysam.Samfile("/home/nader/Documents/thesis/analysis/mapFiles/merged/mergedBM_ST_DNMT2.bam", "rb")
  samfile = pysam.Samfile(bam, "rb")
  mappedReads = 0. #number of mapped reads
  
  
  if verbose : print "Read reference transcript file: Transcript ID should be as in refCDS and as in bam"
  
  transcript = SeqIO.to_dict(SeqIO.parse(open(refSeq), "fasta")) #function seqIO.to_dict() turns a SeqRecord iterator (list) into a dictionary (in memory). 
  
  
  if verbose : print "Main calculations block"
  
  """ Read codon lists in transcript file """
  codonLists = {}
  
  for mRNA in transcript:
    sequence = str(transcript[mRNA].seq) #sequence of the transcript for each mRNA (ucsc id) in a form of a 'seq' object --> change to a string
    sequence = sequence.lower().replace('u','t') #all in lowercase, and U's replaced by T's.
    """
    Getting a list of codons for each CDS: {1st codon= codonList[0], 2nd codon= codonList[1], codonList[2], ... , last codon= codonList[-1]}
    First codon should be the START codon
    Last codon should be the STOP codon
    """  
    codonList = []
    if mRNA in localDic:
      for i in range(localDic[mRNA][0], localDic[mRNA][1], 3): #here, I need to make sure that scope of the variables are correctly taken into account each time. From start to the end of the mRNA (using local coordinates) in increments of 3.
        #  {localStart of mRNA, localStart of mRNA +3, localStart of mRNA +6, ..., localEnd of mRNA -3}
        if (i + 3) < localDic[mRNA][1]:
          codon = sequence[i:i+3]
        else:
          codon = sequence[i:localDic[mRNA][1]]
        
        codon=codon.replace('t', 'u')
        codonList.append(codon)
      
      codonLists[mRNA] = codonList
  
  
  """ Initialisations for codon coverage """
  codonCover_A_site = defaultdict(lambda:defaultdict(int)) #mRNA ID: codon position calculated by considering the local position of the assigned nucleotide and dividing by 3.
  codonCover_P_site = defaultdict(lambda:defaultdict(int))
  codonCover_E_site = defaultdict(lambda:defaultdict(int))
  codonCover_minus1 = defaultdict(lambda:defaultdict(int))
  codonCover_minus2 = defaultdict(lambda:defaultdict(int))
  codonCover_minus3 = defaultdict(lambda:defaultdict(int))
  codonCover_plus1 = defaultdict(lambda:defaultdict(int))
  codonCover_plus2 = defaultdict(lambda:defaultdict(int))
  codonCover_plus3 = defaultdict(lambda:defaultdict(int))
  codonCover_1first = defaultdict(lambda:defaultdict(int))
  codonCover_2first = defaultdict(lambda:defaultdict(int))
  codonCover_3first = defaultdict(lambda:defaultdict(int))
  codonCover_1last = defaultdict(lambda:defaultdict(int))
  codonCover_2last = defaultdict(lambda:defaultdict(int))
  codonCover_3last = defaultdict(lambda:defaultdict(int))
  
  """ Dictionary containing read cound per mRNA """
  dico_mRNA_nbReads = defaultdict(int) #dico_mRNA_nbReads[mRNAID] = nbReads
  dico_mRNA_POS_nbReads = defaultdict(lambda:defaultdict(int)) #dico_mRNA_POS_nbReads[refid][alignread.pos] += 1
  
  """ Initialisations for nucleotide frequencies at boundaries (first 9, last 9) """
  ntCount_1 = defaultdict(int)
  ntCount_2 = defaultdict(int)
  ntCount_3 = defaultdict(int)
  ntCount_4 = defaultdict(int)
  ntCount_5 = defaultdict(int)
  ntCount_6 = defaultdict(int)
  ntCount_7 = defaultdict(int)
  ntCount_8 = defaultdict(int)
  ntCount_9 = defaultdict(int)
  #
  ntCount_minus1 = defaultdict(int)
  ntCount_minus2 = defaultdict(int)
  ntCount_minus3 = defaultdict(int)
  ntCount_minus4 = defaultdict(int)
  ntCount_minus5 = defaultdict(int)
  ntCount_minus6 = defaultdict(int)
  ntCount_minus7 = defaultdict(int)
  ntCount_minus8 = defaultdict(int)
  ntCount_minus9 = defaultdict(int)
  
  
  """
  Assigning a specific A-site nucleotide based on the length of the fragments and 
  converting the nucleotide position to a codon position that would correspond 
  to a specific codon in the codon list. 
  """
  #Nota bene: - alignread.mapping_quality always 0 (unmapped) or 255 (mapped)
  #           - alignread.query_qualities always 14 or more <=> pvalue is lower or equal to 0.05
  for alignread in samfile:
    mappedReads+=1
    if alignread.tid >= 0 :
      if verbose : print "alignread.tid>0"
      # accession name of reference sequence
      refid = samfile.getrname(alignread.tid)
      
      # keep accession only, not the further information
      refid = refid.split('|')[0]
      
      if verbose : print "refid : "+refid
      
      # keep or remove version number
      if versionStrip:
        if verbose : print "versionStrip"
        #transcript version nb to rm in case transcript names in fasta reference contain it but not the annotation file.
        refid = refid.split('.')[0]
      
      if refid in localDic:
        if verbose : print "refid in localDic"
        # start and end of transcript, as read in reference file
        localStart = localDic[refid][0]
        localEnd = localDic[refid][1]
        orientation = localDic[refid][2]
        if (alignread.pos != None) and (25 <= alignread.alen <= 32): 
            
            # increment number of reads aligning to current transcript
            dico_mRNA_nbReads[refid] += 1
            
            # increment number of reads in the A-site conditioned on the existence of a sound main CDS
            if verbose: 
              print "refid             = ",refid
              print "localDic[refid]   = ",localDic[refid]
              print "codonLists[refid] = ",codonLists[refid]
            
            if refid in codonLists :
              # Ignore transcripts which do not have any coding sequence : 
              if codonLists[refid] != []:
                if ( ((localDic[refid][1] - localDic[refid][0]) % 3 == 0) 
                     and (codonLists[refid][0] == "aug") 
                     and ((codonLists[refid][-1] == 'uaa') or (codonLists[refid][-1] == 'uag') or (codonLists[refid][-1] == 'uga')) ):
                  lastPossibleReadStart = localDic[refid][3] - 25 #localDic[refid][3] is the transcript length (or else said, cumulative exons length) ; 25 is the minimum read length
                  if (alignread.pos>=0 and alignread.pos < lastPossibleReadStart) : #make double sure that codon is still inside transcript (5' cap to 3' polyA tail)
                    # A-site position is at +15/read start ; offset=1 because sometimes it is 16 ; integer division leads to codon position = 5 in both cases
                    codonPosition_A_site_metagene = (alignread.pos + 15 - localStart + offset) / 3
                    # increment number of reads
                    dico_mRNA_POS_nbReads[refid][codonPosition_A_site_metagene] += 1
                  
                  
                  if verbose: print "alignread.pos: ",alignread.pos
                  #for forward and reverse strands 
                  #alignread.alen gives length of the reads (between 26-32 nt)    
                  NucleotidePosition_A_site = 0 #position of the A-site on each alignment. 
                  NucleotidePosition_P_site = 0
                  NucleotidePosition_E_site = 0
                  NucleotidePosition_minus1 = 0
                  NucleotidePosition_minus2 = 0
                  NucleotidePosition_minus3 = 0
                  NucleotidePosition_plus1 = 0
                  NucleotidePosition_plus2 = 0
                  NucleotidePosition_plus3 = 0
                  
                  codonPosition_A_site = 0 #position of the A-site in respect to the codon list. 
                  codonPosition_P_site = 0
                  codonPosition_E_site = 0
                  codonPosition_minus1 = 0
                  codonPosition_minus2 = 0
                  codonPosition_minus3 = 0
                  codonPosition_plus1 = 0
                  codonPosition_plus2 = 0
                  codonPosition_plus3 = 0
                  
                  #modif CL/Nader : "if not alignread.is_reverse" is UNNECESSARY - REMOVED => but no 2ce as much coverage everywhere because strand bias in Bowtie.
                  #if 25 <= alignread.alen <= 32:
                  if mini <= alignread.alen <= maxi:   #28 <= alignread.alen <= 31:
                       # Positions in the ribosome are defined wrt the aligned read
                       # Positions of nucleotides
                       # From 5' strand to 3', downstream codons, E site comes first, then P, and A last
                       # Indeed, A accepts the next unread codon, from side 3'.
                       
                       """
                       
                       CODON POSITION RATIONALE
                       ------------------------
                       
                       Ingolia :  - A-site position is at read start + 15 bases if read length is 29, 30,
                                  - +14 possibly if read length is 28,
                                  - +16 if read length is 31, 32 or 33.
                       
                       In practice : 
                                  - recurrence plot should show synchro start at -12/P-site <=> -15/A-site
                                                       or sometimes:         at -13/P-site <=> -16/A-site (corrected if offset is set to 1)  
                                  - offset length sometimes exhibits a large peak at 0[3], but also a 2nd-large peak at 2[3] (Mice and Human)
                                    => peak at 0 means most reads start where a codon starts
                                    => peak at 2 means that many reads start one base before a codon start, because an isolate base clutches at 5'end (side where the coding sequence starts, E-site side of the ribosome)
                                    => A-site position should be +16 instead of +15 in this case
                                    => offset=1 given in input of this script
                                    => integer division (NucleotidePosition_xxxx - localStart + offset) /3
                                       So that A-site of reads beginning with codon start are correctly assigned (pos 0 + 15, int division : 5 ; with offset +1: 0 + 16, int division is still 5), 
                                       But also A-site of reads beginning with 1 base before codon start are correctly assigned (pos 0 + 15, int division : 4 which is false ; with offset +1: 0 + 16, int division is 5 which is correct.).
                                   
                                    => S.Pombe : peak at 1[3] often occurs (sometimes more often than 0[3])
                                    =>           means that one base (-12) was chomped away during the experiment
                                    =>           A-site position should be +14 instead of +15 in the case of S.Pombe
                                    =>           already taken care of by integer division, a +1 offset does not influence that.
                                   
                                   - other sites (P, E, etc.) should benefit also from this rationale. 
                                   
                       Case 1st 3, last 3 codons : no offset, raw position   
                       
                       """
                       NucleotidePosition_1last   = alignread.pos + alignread.alen
                       NucleotidePosition_2last = alignread.pos + alignread.alen - 3
                       NucleotidePosition_3last = alignread.pos + alignread.alen - 6
                       #
                       NucleotidePosition_minus3 = alignread.pos + 24
                       NucleotidePosition_minus2 = alignread.pos + 21
                       NucleotidePosition_minus1 = alignread.pos + 18
                       NucleotidePosition_A_site = alignread.pos + 15
                       NucleotidePosition_P_site = alignread.pos + 12
                       NucleotidePosition_E_site = alignread.pos + 9
                       NucleotidePosition_plus1  = alignread.pos + 6
                       NucleotidePosition_plus2  = alignread.pos + 3
                       NucleotidePosition_plus3  = alignread.pos
                       #
                       # Positions of codons, conditioned to their lying in coding sequence
                       # from which the first 15 aminoacids are excluded: [start+45 ; end]
                       #   And: update coverage . for corresponding reference sequence (refid)
                       #                        . at codon position
                       #
                       #if (NucleotidePosition_minus3 < localEnd-15) and (NucleotidePosition_minus3 > (localStart+45)):                  # change v10/v9: exclude 5 last codons
                       # boundaries codons relate to ligation bias, which isn't affected by incomplete-digestion-realted offset
                       codonPosition_1last  = (NucleotidePosition_1last - localStart) / 3
                       codonPosition_2last  = (NucleotidePosition_2last - localStart) / 3
                       codonPosition_3last  = (NucleotidePosition_3last - localStart) / 3
                       #
                       codonPosition_minus3 = (NucleotidePosition_minus3 - localStart + offset) / 3 
                       codonPosition_minus2 = (NucleotidePosition_minus2 - localStart + offset) / 3
                       codonPosition_minus1 = (NucleotidePosition_minus1 - localStart + offset) / 3
                       codonPosition_A_site = (NucleotidePosition_A_site - localStart + offset) / 3
                       codonPosition_P_site = (NucleotidePosition_P_site - localStart + offset) / 3
                       codonPosition_E_site = (NucleotidePosition_E_site - localStart + offset) / 3
                       codonPosition_plus1  = (NucleotidePosition_plus1  - localStart + offset) / 3
                       codonPosition_plus2  = (NucleotidePosition_plus2  - localStart + offset) / 3
                       codonPosition_plus3  = (NucleotidePosition_plus3  - localStart + offset) / 3
                       #
                       codonPosition_3first  = (NucleotidePosition_plus1 - localStart) / 3
                       codonPosition_2first  = (NucleotidePosition_plus2 - localStart) / 3
                       codonPosition_1first  = (NucleotidePosition_plus3 - localStart) / 3
                       
                       
                       # boundaries codons relate to ligation bias, which is potentially affected (different codon usage) in CDS start,
                       # thus the same rationale is used than for A-site, etc codons 
                       
                       if ( (codonPosition_1last >= 15) and (codonPosition_1last < (len(codonLists[refid])-5)) ) :
                         codonCover_1last[refid][codonPosition_1last] += 1
                       
                       if ( (codonPosition_2last >= 15) and (codonPosition_2last < (len(codonLists[refid])-5)) ) :
                         codonCover_2last[refid][codonPosition_2last] += 1
                       
                       if ( (codonPosition_3last >= 15) and (codonPosition_3last < (len(codonLists[refid])-5)) ) :
                         codonCover_3last[refid][codonPosition_3last] += 1
                       
                       if ( (codonPosition_minus3 >= 15) and (codonPosition_minus3 < (len(codonLists[refid])-5)) ) :
                         codonCover_minus3[refid][codonPosition_minus3] += 1
                       
                       if ( (codonPosition_minus2 >= 15) and (codonPosition_minus2 < (len(codonLists[refid])-5)) ) :
                         codonCover_minus2[refid][codonPosition_minus2] += 1
                       
                       if ( (codonPosition_minus1 >= 15) and (codonPosition_minus1 < (len(codonLists[refid])-5)) ) :
                         codonCover_minus1[refid][codonPosition_minus1] += 1
                       
                       if ( (codonPosition_A_site >= 15) and (codonPosition_A_site < (len(codonLists[refid])-5)) ) :
                         codonCover_A_site[refid][codonPosition_A_site] += 1
                       
                       if ( (codonPosition_P_site >= 15) and (codonPosition_P_site < (len(codonLists[refid])-5)) ) :
                         codonCover_P_site[refid][codonPosition_P_site] += 1
                       
                       if ( (codonPosition_E_site >= 15) and (codonPosition_E_site < (len(codonLists[refid])-5)) ) :
                         codonCover_E_site[refid][codonPosition_E_site] += 1
                       
                       if ( (codonPosition_plus1 >= 15) and (codonPosition_plus1 < (len(codonLists[refid])-5)) ) :
                         codonCover_plus1[refid][codonPosition_plus1] += 1
                       
                       if ( (codonPosition_plus2 >= 15) and (codonPosition_plus2 < (len(codonLists[refid])-5)) ) :
                         codonCover_plus2[refid][codonPosition_plus2] += 1
                       
                       if ( (codonPosition_plus3 >= 15) and (codonPosition_plus3 < (len(codonLists[refid])-5)) ) :
                         codonCover_plus3[refid][codonPosition_plus3] += 1
                       
                       if ( (codonPosition_3first >= 15) and (codonPosition_3first < (len(codonLists[refid])-5)) ) :
                         codonCover_3first[refid][codonPosition_3first] += 1
                       
                       if ( (codonPosition_2first >= 15) and (codonPosition_2first < (len(codonLists[refid])-5)) ) :
                         codonCover_2first[refid][codonPosition_2first] += 1
                       
                       if ( (codonPosition_1first >= 15) and (codonPosition_1first < (len(codonLists[refid])-5)) ) :
                         codonCover_1first[refid][codonPosition_1first] += 1
                       
                       
                       """
                       For nt boundaries, consider a smaller subset, simplified procedure, where the 
                         full footprint should be far from start and stop.
                       """
                       if ( (codonPosition_1first >= 15) and (codonPosition_1last < (len(codonLists[refid])-5)) ) :
                          #
                          ntCount_1[alignread.seq[0]] += 1
                          ntCount_2[alignread.seq[1]] += 1
                          ntCount_3[alignread.seq[2]] += 1
                          ntCount_4[alignread.seq[3]] += 1
                          ntCount_5[alignread.seq[4]] += 1
                          ntCount_6[alignread.seq[5]] += 1
                          ntCount_7[alignread.seq[6]] += 1
                          ntCount_8[alignread.seq[7]] += 1
                          ntCount_9[alignread.seq[8]] += 1
                          #
                          ntCount_minus1[alignread.seq[-1]] += 1
                          ntCount_minus2[alignread.seq[-2]] += 1
                          ntCount_minus3[alignread.seq[-3]] += 1
                          ntCount_minus4[alignread.seq[-4]] += 1
                          ntCount_minus5[alignread.seq[-5]] += 1
                          ntCount_minus6[alignread.seq[-6]] += 1
                          ntCount_minus7[alignread.seq[-7]] += 1
                          ntCount_minus8[alignread.seq[-8]] += 1
                          ntCount_minus9[alignread.seq[-9]] += 1
  
  
  """ million mapped reads """
  millionMapped = float(mappedReads) / float(1e6)
  
  codonUsage_A_site = defaultdict(lambda:defaultdict(int)) #codon and its count associated with it for each gene.
  codonUsage_P_site = defaultdict(lambda:defaultdict(int)) 
  codonUsage_E_site = defaultdict(lambda:defaultdict(int)) 
  codonUsage_minus1 = defaultdict(lambda:defaultdict(int)) 
  codonUsage_minus2 = defaultdict(lambda:defaultdict(int)) 
  codonUsage_minus3 = defaultdict(lambda:defaultdict(int)) 
  codonUsage_plus1 = defaultdict(lambda:defaultdict(int)) 
  codonUsage_plus2 = defaultdict(lambda:defaultdict(int)) 
  codonUsage_plus3 = defaultdict(lambda:defaultdict(int)) 
  codonUsage_1last = defaultdict(lambda:defaultdict(int)) 
  codonUsage_2last = defaultdict(lambda:defaultdict(int)) 
  codonUsage_3last = defaultdict(lambda:defaultdict(int)) 
  codonUsage_1first = defaultdict(lambda:defaultdict(int)) 
  codonUsage_2first = defaultdict(lambda:defaultdict(int)) 
  codonUsage_3first = defaultdict(lambda:defaultdict(int)) 
  
  # A dictionary for total counts of each codon
  codonCount_A_site = defaultdict(int)
  codonCount_P_site = defaultdict(int)
  codonCount_E_site = defaultdict(int)
  codonCount_minus1 = defaultdict(int)
  codonCount_minus2 = defaultdict(int)
  codonCount_minus3 = defaultdict(int)
  codonCount_plus1 = defaultdict(int)
  codonCount_plus2 = defaultdict(int)
  codonCount_plus3 = defaultdict(int)
  codonCount_3last = defaultdict(int)
  codonCount_2last = defaultdict(int)
  codonCount_1last = defaultdict(int)
  codonCount_3first = defaultdict(int)
  codonCount_2first = defaultdict(int)
  codonCount_1first = defaultdict(int)
  
  # Initialize for each codon
  fullcodonlist=['aaa','aac','aag','aau','aca','acc','acg','acu','aga','agc','agg','agu','aua',
                 'auc','aug','auu','caa','cac','cag','cau','cca','ccc','ccg','ccu','cga','cgc',
                 'cgg','cgu','cua','cuc','cug','cuu','gaa','gac','gag','gau','gca','gcc','gcg',
                 'gcu','gga','ggc','ggg','ggu','gua','guc','gug','guu','uaa','uac','uag','uau',
                 'uca','ucc','ucg','ucu','uga','ugc','ugg','ugu','uua','uuc','uug','uuu']
  for item in fullcodonlist :
    codonCount_A_site[item]=0
    codonCount_P_site[item]=0
    codonCount_E_site[item]=0
    codonCount_minus1[item]=0
    codonCount_minus2[item]=0
    codonCount_minus3[item]=0
    codonCount_plus1[item]=0
    codonCount_plus2[item]=0
    codonCount_plus3[item]=0
    codonCount_3last[item]=0
    codonCount_2last[item]=0
    codonCount_1last[item]=0
    codonCount_3first[item]=0
    codonCount_2first[item]=0
    codonCount_1first[item]=0
  
  # A dictionary for mRNA-level counts of each codon
  dico_mRNA_codon_count = {}
  
  # A dictionary for occupancy at codon resolution
  SingleCodonOccupancy = defaultdict(lambda:defaultdict(int))
  
  # A dictionary for strong paused sites
  lists_mRNA_paused_sites = []
  lists_mRNA_paused_sites_names=['mRNA','codon-pos-rel-CDSstart','codon','nb-footprint-per-gene-median','sequenceContext']
  lists_mRNA_paused_sites.append(lists_mRNA_paused_sites_names)
  
  # Intialize intermediate variables for Hussmann calculation of mean and sd of codon enrichment
  #   over transcripts
  # For Hussmann enrichment, mean and variance
  
  #codon_i_seq_count   = defaultdict(lambda:defaultdict(int))
  s_IF                = defaultdict(lambda:defaultdict(int))
  s_IF_H              = defaultdict(lambda:defaultdict(int))
  sum_e_gi            = defaultdict(lambda:defaultdict(float))
  sum2_e_gi           = defaultdict(lambda:defaultdict(float))
  sum_codonEnrichmnt  = defaultdict(lambda:defaultdict(float))
  sum2_codonEnrichmnt = defaultdict(lambda:defaultdict(float))
  T_ci                = defaultdict(lambda:defaultdict(float))
  U_ci                = defaultdict(lambda:defaultdict(float))
  W_i                 = defaultdict(int)
  codon_usage_glob    = defaultdict(float)
  
  # Initialize to 0 or 0.
  for codonSeq in fullcodonlist :
    for iCodon in range(-ex3,ex5+1) :
      #codon_i_seq_count[codonSeq][iCodon]   = 0
      s_IF[codonSeq][iCodon]                = 0
      s_IF_H[codonSeq][iCodon]              = 0
      sum_e_gi[codonSeq][iCodon]            = 0.
      sum2_e_gi[codonSeq][iCodon]           = 0.
      sum_codonEnrichmnt[codonSeq][iCodon]  = 0.
      sum2_codonEnrichmnt[codonSeq][iCodon] = 0.
      T_ci[codonSeq][iCodon]                = 0.
      U_ci[codonSeq][iCodon]                = 0.
      W_i[iCodon]                           = 0
      codon_usage_glob[codonSeq]            = 0.
  
  #####
  """Codon usage for each mRNA, summary over all mRNAs"""
  compteur=0
  # Nb of genes where actual calculations take place
  n_g     =0
  Total=len(transcript)
  for mRNA in transcript:
    if verbose :
      compteur += 1
      if compteur % 1000 == 0 or compteur == 1 or compteur == Total :
        print "mRNA "+str(compteur)+" out of "+str(Total)
    
    """
    if sensible coding sequences (3*n long, beginning with AUG, ending with TAA, TAG or TGA):
       codon Usage at each site for this CDS/target/mRNA is recorded
    """
    if mRNA in codonLists :
      if codonLists[mRNA] != []:
        if ((localDic[mRNA][1] - localDic[mRNA][0]) % 3 == 0) and (codonLists[mRNA][0] == "aug"): #check if the transcript is divisible by 3 and starts with an AUG. 
          if (codonLists[mRNA][-1] == 'uaa') or (codonLists[mRNA][-1] == 'uag') or (codonLists[mRNA][-1] == 'uga'): #only those CDS that end with a stop codon
            """
            codon usage at transcript level 
            = count(each codon, in current mRNA)        
            #codons in CDS=  dico_mRNA_codon_count[mRNA][codon]                          
            """
            #set(codonLists[mRNA]) #unique elements
            # dictionary giving count of each codon, in current mRNA.
            dico_mRNA_codon_count[mRNA] = dict((x,codonLists[mRNA].count(x)) for x in set(codonLists[mRNA]))
            # Codon usage in CDS restricted to AUG + 15 codons to STOP - 5 codons
            codon_usage_list_g = codonLists[mRNA][15:-5]
            #
            if len(codon_usage_list_g) > 0 :
              #codon_usage_g = dict((x,1.) for x in fullcodonlist)
              codon_usage_g = dict((x, float(codon_usage_list_g.count(x))/float(len(codon_usage_list_g)) ) for x in set(codon_usage_list_g))
            else :
              if verbose : 
                print "Warning : len(codon_usage_list_g) is non positive, calculation of codon_usage_g is impacted."
                print "          len(codon_usage_list_g) = "+str(len(codon_usage_list_g))
            
            """
            codon usage
            and, for A-site: single codon occupancy and paused sites  
            """
            # Initialize sum_r_gi and codonEnrichmnt
            sum_r_gi       = defaultdict(lambda:defaultdict(int))
            sum2_r_gi      = defaultdict(lambda:defaultdict(int))
            n_cig          = defaultdict(lambda:defaultdict(int))
            t_cig          = defaultdict(lambda:defaultdict(float))
            u_cig          = defaultdict(lambda:defaultdict(float))
            w_ig           = defaultdict(int)
            #codonEnrichmnt = defaultdict(lambda:defaultdict(float))
            seq_iCodon_TF  = defaultdict(lambda:defaultdict(bool))
            s_IF_thisgene  = defaultdict(lambda:defaultdict(int))
            s_IF_thisgene_H  = defaultdict(lambda:defaultdict(int))
            
            # Local codon usage of this gene, for readable formula of s_IF
            a = codonLists[mRNA]
            
            # Set initial value to 0 or false
            for codonSeq in fullcodonlist: 
              for iCodon in range(-ex3,ex5+1) :
                sum_r_gi[codonSeq][iCodon]       = 0
                sum2_r_gi[codonSeq][iCodon]      = 0
                n_cig[codonSeq][iCodon]          = 0
                t_cig[codonSeq][iCodon]          = 0.
                u_cig[codonSeq][iCodon]          = 0.
                w_ig[iCodon]                     = 0
                #codonEnrichmnt[codonSeq][iCodon] = 0.
                seq_iCodon_TF[codonSeq][iCodon]  = False
                s_IF_thisgene[codonSeq][iCodon]  = 0
                s_IF_thisgene_H[codonSeq][iCodon]  = 0
            
            # CodonList and codon usage of this mRNA for Hussmann calculations
            codonList                   = codonLists[mRNA]
            listdico                    = [codonCover_A_site[mRNA][position] for position in codonCover_A_site[mRNA] if position >= ex5 and position < (len(codonList)-ex3)]
            #listdico                    = [codonCover_A_site[mRNA][position] for position in codonCover_A_site[mRNA]]
            #
            # Calculate s_IF_thisgene
            CDSlen=len(codonList)
            #print "mRNA="+mRNA
            for pos in range(15, (CDSlen-5)) :
              for iCodon in range(-ex3,ex5+1) :
                # Eligible position ?
                if pos >= max(15, 15+iCodon) and pos < min((CDSlen-5), (CDSlen-5+iCodon)) :
                  codonSeq = codonList[pos-iCodon]
                  s_IF_thisgene[codonSeq][iCodon] += 1
                  ## Cardinal of set s_IF = {(g,i) : d<i<l_g-d, c_g,i-F = I}
                  ##     with I = codon identity / sequence
                  ##          F = offset / iCodon
                  ## = Codon usage stratified by codonSeq and iCodon
                
                if pos >= ex5 and pos < (CDSlen-ex3) :
                  codonSeq = codonList[pos-iCodon]
                  s_IF_thisgene_H[codonSeq][iCodon] += 1
            
            
            # strictly Hussmann
            read_count_wo_excl_dom      = 0
            if len(listdico) > 0 :
              meanCoverage_Asite_thismRNA      = float(sum(listdico)) / float(len(listdico)) 
            
            # Codon coverage and expected nb reads
            for position in codonCover_A_site[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_A_site[mRNA][codonLists[mRNA][position]] += codonCover_A_site[mRNA][position] 
                #if codonCover_A_site[mRNA][position] > 3:
                #"""
                #Single Codon Occupancy at A site
                #"""
                # nb of codons from CDS start to CDS end.
                nbCodons          =  ( localDic[mRNA][1] - localDic[mRNA][0] ) / 3  #- 20 
                # nb of codons from this specific sort in this specific mRNA
                codonCountInCDS   =  dico_mRNA_codon_count[mRNA][codonLists[mRNA][position]]
                ## if collapsed over identical codons from different positions along gene
                #ExpectedmRNACodon =  dico_mRNA_nbReads[mRNA] / float(nbCodons * codonCountInCDS)
                # as in NFC or SCO, at codon resolution
                if nbCodons >0 :
                  ExpectedmRNACodon =  dico_mRNA_nbReads[mRNA] / float(nbCodons)
                elif messages :
                  print "Warning : nbCodons is non positive, calculation of ExpectedmRNACodon is impacted."
                  print "          nbCodons = "+str(nbCodons)
                
                if ExpectedmRNACodon > 0 :
                  SingleCodonOccupancy[mRNA][position] = codonCover_A_site[mRNA][position] / ExpectedmRNACodon
                elif messages :
                  print "Warning : ExpectedmRNACodon is non positive, calculation of SingleCodonOccupancy["+mRNA+"]["+str(position)+"] is impacted."
                  print "          ExpectedmRNACodon = "+str(ExpectedmRNACodon)
              
              
              ### paused sites in 5-codons window
              #   positions 0-14 excluded because 15 first codons are biased by residual initiation and by AUG-related structure, or else,
              #   positions 15-24 also excluded because they do not have 10 valid codons on the 5' side,
              #   positions CDSend - 5 codons excluded because 5 codons could be biased by leakage near stop codon,
              #             CDSend - 15-10 codons excluded because they would not have 10 valid codons on the 3' side.
              if position > 24 and position < (nbCodons - 15):
                # median coverage at A site +/- 2 codons
                if (position-2) in codonCover_A_site[mRNA]: covmoins2 = codonCover_A_site[mRNA][position-2]
                else:                                       covmoins2 = 0
                if (position-1) in codonCover_A_site[mRNA]: covmoins1 = codonCover_A_site[mRNA][position-1]
                else:                                       covmoins1 = 0
                if (position+1) in codonCover_A_site[mRNA]: covplus1 = codonCover_A_site[mRNA][position+1]
                else:                                       covplus1 = 0
                if (position+2) in codonCover_A_site[mRNA]: covplus2 = codonCover_A_site[mRNA][position+2]
                else:                                       covplus2 = 0
                covA = codonCover_A_site[mRNA][position]
                mediancov_5codons = np.median( [ covmoins2, 
                                                 covmoins1,
                                                 covA,
                                                 covplus1,
                                                 covplus2 ] )
                # consider only where median coverage is 2 or more
                if mediancov_5codons >= 2:
                  # identify paused sites
                  if ExpectedmRNACodon > 0 :
                    is_peak = codonCover_A_site[mRNA][position] / ExpectedmRNACodon
                  elif messages :
                    print "Warning : ExpectedmRNACodon is non positive, calculation of is_peak is impacted."
                    print "          ExpectedmRNACodon = "+str(ExpectedmRNACodon)
                  
                  if is_peak > 25:
                    # store ribosome pause site information ; store -10 -> +10 codons relative to the paused site
                    lists_mRNA_paused_sites.append( [ mRNA, 
                                                    position, 
                                                    codonLists[mRNA][position], 
                                                    is_peak,
                                                    ''.join([codonLists[mRNA][x] for x in range(position-10,position+11)]) ] )
              
              
              #### Hussmann and own - update enrichment at each position, each codon
              #   
              # iCodon from -ex5 to +ex3 [-90 ; 90] unless in excluded domain
              #if position >= ex5 and position < (len(codonList) - ex3) :
              if position >= 15 and position < (len(codonList) - 5) :
                ## exclude domain [0;ex5]
                #iCodon_min = -ex5
                #if position < (2*ex5) :
                #  iCodon_min = - (position - ex5)
                # 
                ## exclude domain [len-90;len]
                #iCodon_max = ex3
                #if position > (len(codonList) - (2*ex3)) :
                #  iCodon_max = len(codonList) - ex3 - position
                #
                # Update Enrichment
                #for iCodon in range(iCodon_min, (iCodon_max+1)) :
                for iCodon in range(-ex3,ex5+1) :
                  #if position >= max(15, 15-iCodon) and position < min((CDSlen-5), (CDSlen-5-iCodon)) :
                  # (pos >=15 and pos-iCodon >= 15) and (pos < CDSlen-5 and pos-iCodon < CDSlen-5)
                  if position >= max(15, 15+iCodon) and position < min((CDSlen-5), (CDSlen-5+iCodon)) :
                    # local utilities
                    codonSeq                         = codonList[(position-iCodon)]
                    covHere                          = codonCover_A_site[mRNA][position]
                    # update
                    # Strictly Hussmann - sum(r_gi)_i=d,i=l_g-d, stratified by codon sequence (identity) and iCodon,
                    #   in this gene
                    if position >= ex5 and position < (len(codonList) - ex3) :
                      seq_iCodon_TF[codonSeq][iCodon]  = True
                      sum_r_gi[codonSeq][iCodon]  += covHere
                      sum2_r_gi[codonSeq][iCodon] += covHere**2
                    
                    # Codon enrichment own
                    n_cig[codonSeq][iCodon] += float(covHere)
                    # codon usage global et local : if codon_usage_g[codonSeq] > 0. :
                    # codon_usage_s_IF : if s_IF_thisgene[codonSeq][iCodon] > 0 :
                    if codon_usage_g[codonSeq] > 0. :
                      # codon usage local and iCodon-dependent (codon_usage_s_IF is in [0 ; 1])
                      # codon_usage_s_IF : codon_usage_s_IF = float(s_IF_thisgene[codonSeq][iCodon]) / float(sum([s_IF_thisgene[I][iCodon] for I in s_IF_thisgene]))
                      # codon_usage_s_IF : t_cig[codonSeq][iCodon] += float(covHere) / float(codon_usage_s_IF)
                      # codon usage global : 
                      t_cig[codonSeq][iCodon] += float(covHere) #/ float(codon_usage_g[codonSeq])
                      # codonUsage local   : t_cig[codonSeq][iCodon] += float(covHere) / float(codon_usage_g[codonSeq])
                    elif verbose :
                      # codon_usage_s_IF : print "Warning : s_IF_thisgene[codonSeq][iCodon] is non positive, calculation of t_cig[codonSeq][iCodon] is impacted."
                      # print "          s_IF_thisgene[codonSeq][iCodon] = "+str(s_IF_thisgene[codonSeq][iCodon])
                      # codon usage global et local :
                      print "Warning : codon_usage_g[codonSeq] is non positive, calculation of t_cig[codonSeq][iCodon] is impacted."
                      print "          codon_usage_g[codonSeq] = "+str(codon_usage_g[codonSeq])
                      print "          g = mRNA = "+mRNA
                      print "          codonSeq = "+codonSeq
                      print "          iCodon   = "+str(iCodon)
            
            
            # Hussmann - this is M_g - stratified by iCodon since sum_r_gi is itself stratified by iCodon.
            M_g         = {}
            #len_support = len(codonList) - ex3 - ex5 + 1
            len_support = max(1,len(codonList) - ex3 - ex5 + 1)
            for iCodon in range(-ex3, ex5+1) :
               sum_r_gi_by_iCodon = sum([sum_r_gi[codonSeq][iCodon] for codonSeq in sum_r_gi]) 
               # M_g is not calculated in case support is null.
               if len_support > 0 :
                 #M_g[iCodon]        = float(sum_r_gi_by_iCodon) / float(len_support)
                 M_g[iCodon]        = max(0.1, float(sum_r_gi_by_iCodon) / float(len_support))
                 
                 # Own, but keep filter on M_g
                 if M_g[iCodon] >= 0.1 :
                   w_ig[iCodon] = sum([n_cig[codonSeq][iCodon] for codonSeq in n_cig])
                   W_i[iCodon]  += w_ig[iCodon]
                   # count of mRNAs, for simplified standard error calculation
                   n_g          += 1
                   
                   for codonSeq in t_cig :
                     # Test on t_cig's positiveness is necessary because this cumul 
                     #   is always initialized at 0, even if non defined for local 
                     #   codonSeq and iCodon.
                     if t_cig[codonSeq][iCodon] > 0. :
                       T_ci[codonSeq][iCodon] += t_cig[codonSeq][iCodon]
                       if w_ig[iCodon] > 0 :
                         u_cig[codonSeq][iCodon] = (t_cig[codonSeq][iCodon])**2 / w_ig[iCodon]
                         U_ci[codonSeq][iCodon] += u_cig[codonSeq][iCodon]
                       elif verbose :
                         print "Warning : w_ig[iCodon] is non positive, calculation of u_cig[codonSeq][iCodon] is impacted."
                         print "          w_ig[iCodon] = "+str(w_ig[iCodon])                   
                         print "          g = mRNA     = "+mRNA
                         print "          codonSeq     = "+codonSeq
                         print "          iCodon       = "+str(iCodon)
                     
                     # codon_usage_list_glob
                     if iCodon == 0 and codonSeq in codon_usage_g :
                       codon_usage_glob[codonSeq] += w_ig[0] * codon_usage_g[codonSeq]
            
            
            #### Hussmann - Summarize over transcripts
            # (this cannot be done before, since a same codon, same iCodon, can be found in multiple [-ex5 ; +ex3] windows around AUG
            
            for codonSeq in sum_r_gi :
              
              for iCodon in sum_r_gi[codonSeq] :
                if iCodon in M_g :
                  if M_g[iCodon] >= 0.1 and seq_iCodon_TF[codonSeq][iCodon] :
                    # Count of pairs (g=mRNA,i=position) corresponding to (codonSeq, iCodon)
                    #codon_i_seq_count[codonSeq][iCodon]   +=1
                    s_IF[codonSeq][iCodon]                += s_IF_thisgene[codonSeq][iCodon]
                    s_IF_H[codonSeq][iCodon]              += s_IF_thisgene_H[codonSeq][iCodon]
                    #
                    # Relative enrichment e_gi by codon identity I and offset F (resp. codonSeq and iCodon)
                    # 
                    # Strictly Hussmann : 
                    #    (sum_r_gi[codonSeq][iCodon] / M_g[iCodon]) which is 
                    #    - the sum of e_gi over all (g=this gene, i) in s_IF_thisgene,
                    #    - for this codon identity I (sequence) and
                    #    - for this F=iCodon
                    #    then, variable "sum_e_gi" collects e_gi by gene over all genes, but is still stratified
                    #    - by I and F (identity/sequence, offset/iCodon)
                    #
                    # for sum and ultimately mean
                    sum_e_gi[codonSeq][iCodon]            += sum_r_gi[codonSeq][iCodon] / M_g[iCodon]
                    #sum_codonEnrichmnt[codonSeq][iCodon]  += codonEnrichmnt[codonSeq][iCodon]
                    # for variance
                    sum2_e_gi[codonSeq][iCodon]           += sum2_r_gi[codonSeq][iCodon] / (M_g[iCodon])**2
                    #sum2_codonEnrichmnt[codonSeq][iCodon] += (codonEnrichmnt[codonSeq][iCodon])**2
            
            
            # rest of positions
            
            for position in codonCover_P_site[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_P_site[mRNA][codonLists[mRNA][position]] += codonCover_P_site[mRNA][position] 
            
            for position in codonCover_E_site[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_E_site[mRNA][codonLists[mRNA][position]] += codonCover_E_site[mRNA][position] 
            
            for position in codonCover_minus1[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_minus1[mRNA][codonLists[mRNA][position]] += codonCover_minus1[mRNA][position] 
            
            for position in codonCover_minus2[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_minus2[mRNA][codonLists[mRNA][position]] += codonCover_minus2[mRNA][position] 
            
            for position in codonCover_minus3[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_minus3[mRNA][codonLists[mRNA][position]] += codonCover_minus3[mRNA][position] 
            
            for position in codonCover_plus1[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_plus1[mRNA][codonLists[mRNA][position]] += codonCover_plus1[mRNA][position] 
            
            for position in codonCover_plus2[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_plus2[mRNA][codonLists[mRNA][position]] += codonCover_plus2[mRNA][position] 
            
            for position in codonCover_plus3[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_plus3[mRNA][codonLists[mRNA][position]] += codonCover_plus3[mRNA][position]
            
            for position in codonCover_1first[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_1first[mRNA][codonLists[mRNA][position]] += codonCover_1first[mRNA][position]
            
            for position in codonCover_2first[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_2first[mRNA][codonLists[mRNA][position]] += codonCover_2first[mRNA][position]
            
            for position in codonCover_3first[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_3first[mRNA][codonLists[mRNA][position]] += codonCover_3first[mRNA][position]
            
            for position in codonCover_1last[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_1last[mRNA][codonLists[mRNA][position]] += codonCover_1last[mRNA][position]
            
            for position in codonCover_2last[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_2last[mRNA][codonLists[mRNA][position]] += codonCover_2last[mRNA][position]
            
            for position in codonCover_3last[mRNA]: 
              if position != len(codonLists[mRNA]): 
                codonUsage_3last[mRNA][codonLists[mRNA][position]] += codonCover_3last[mRNA][position]
  
  """
  Reproduce Hussmann et al (2015) PLoS Genetics 11(12):
       Extension to ex3(90) codons upstream (3'  "minus") / A P E / ex5(90) codons downstream (5' "plus") 
  - get mean and sd of codon coverage and codon enrichment     
  """
  
  if verbose : print "Calculations for sum_r_gi and codonEnrichmnt"
  
  meancodonEnrichmnt_H  = defaultdict(lambda:defaultdict(float))
  meancodonEnrichmnt    = defaultdict(lambda:defaultdict(float))
  sdevcodonEnrichmnt_H  = defaultdict(lambda:defaultdict(float))
  sdevcodonEnrichmnt    = defaultdict(lambda:defaultdict(float))
  serrcodonEnrichmnt    = defaultdict(lambda:defaultdict(float))
  
  # Normalize codon_usage_glob[codonSeq]
  for codonSeq in fullcodonlist :
    codon_usage_glob[codonSeq] = codon_usage_glob[codonSeq] / W_i[0]
  
  # Mean and dispersion (sd) of codon enrichment
  for codonSeq in fullcodonlist :
    for iCodon in range(-ex3,ex5+1) :
      # Own
      if W_i[iCodon] > 0 and codon_usage_glob[codonSeq] > 0 :
        # codon usage global :
        meancodonEnrichmnt[codonSeq][iCodon] = T_ci[codonSeq][iCodon] / float(W_i[iCodon]) / codon_usage_glob[codonSeq]
        sdevcodonEnrichmnt[codonSeq][iCodon] = (U_ci[codonSeq][iCodon] / float(W_i[iCodon]) / (codon_usage_glob[codonSeq])**2 - (meancodonEnrichmnt[codonSeq][iCodon])**2)**(0.5)
        serrcodonEnrichmnt[codonSeq][iCodon] = sdevcodonEnrichmnt[codonSeq][iCodon] / float(n_g)**0.5
        # codon usage local (including local-stratified using s_IF) : meancodonEnrichmnt[codonSeq][iCodon] = T_ci[codonSeq][iCodon] / float(W_i[iCodon])
        #sdevcodonEnrichmnt[codonSeq][iCodon] = (U_ci[codonSeq][iCodon] / float(W_i[iCodon]) - (meancodonEnrichmnt[codonSeq][iCodon])**2)**(0.5)
      else :
        if W_i[iCodon] <= 0 and verbose :
            print "Warning : W_i[iCodon] is non positive, meancodonEnrichmnt[codonSeq][iCodon] is set to NA."
            print "          W_i[iCodon]  = "+str(W_i[iCodon])                   
            print "          iCodon       = "+str(iCodon)
        
        meancodonEnrichmnt[codonSeq][iCodon] = "NA"
        sdevcodonEnrichmnt[codonSeq][iCodon] = "NA"
      
      # N
      #nf   = float(codon_i_seq_count[codonSeq][iCodon])
      nf_H = float(s_IF_H[codonSeq][iCodon])
      #if nf > 1. and nf_H > 1. :
      if nf_H > 1. :
        # Relative enrichment, mean
        meancodonEnrichmnt_H[codonSeq][iCodon] = float(sum_e_gi[codonSeq][iCodon]) / nf_H
        #meancodonEnrichmnt[codonSeq][iCodon]   = float(sum_codonEnrichmnt[codonSeq][iCodon]) / nf
        # Esperance X squared, for sd calculations below
        Expectation_Xsquare_codonEnrichmnt_H   = float(sum2_e_gi[codonSeq][iCodon]) / nf_H
        #Expectation_Xsquare_codonEnrichmnt     = float(sum2_codonEnrichmnt[codonSeq][iCodon]) / nf
        # Standard deviation, unbiased
        sqrt_n_nm1_H                           = ( nf_H / (nf_H - 1.) )**(0.5)
        # numerically, E[X]**2 and E[X**2] are sometimes equal with numerical precision leading to 
        #      the latter being smaller than the former -> 0. in this case
        EX_2_m_E_X2_H = Expectation_Xsquare_codonEnrichmnt_H - (meancodonEnrichmnt_H[codonSeq][iCodon])**2
        if EX_2_m_E_X2_H >= 0 :
          sdevcodonEnrichmnt_H[codonSeq][iCodon] = sqrt_n_nm1_H * (EX_2_m_E_X2_H)**(0.5)
        else :
          meancodonEnrichmnt_H[codonSeq][iCodon] = 1. 
          sdevcodonEnrichmnt_H[codonSeq][iCodon] = 0.
        
        ## Idem for alternative calculation method
        #sqrt_n_nm1                             = ( nf / (nf - 1.) )**(0.5)
        #EX_2_m_E_X2 = Expectation_Xsquare_codonEnrichmnt - (meancodonEnrichmnt[codonSeq][iCodon])**2
        #if EX_2_m_E_X2 >= 0 :
        #  sdevcodonEnrichmnt[codonSeq][iCodon]   = sqrt_n_nm1   * (EX_2_m_E_X2)**(0.5)
        #else :
        #  meancodonEnrichmnt[codonSeq][iCodon] = 1.
        #  sdevcodonEnrichmnt[codonSeq][iCodon] = 0.
  
  
  """
  
  END Calculations
  
  -
  
  START Output to files
  
  """
  
  if verbose : print "End of main calculation block - start output to files"
  
  """
  Codon Counts across all mRNAs
  """
  
  for mRNA in codonUsage_A_site:
    for codon in codonUsage_A_site[mRNA]:
      codonCount_A_site[codon] += codonUsage_A_site[mRNA][codon]
  
  for mRNA in codonUsage_P_site:
    for codon in codonUsage_P_site[mRNA]:
      codonCount_P_site[codon] += codonUsage_P_site[mRNA][codon]
  
  for mRNA in codonUsage_E_site:
    for codon in codonUsage_E_site[mRNA]:
      codonCount_E_site[codon] += codonUsage_E_site[mRNA][codon]
  
  for mRNA in codonUsage_minus1:
    for codon in codonUsage_minus1[mRNA]:
      codonCount_minus1[codon] += codonUsage_minus1[mRNA][codon]
  
  for mRNA in codonUsage_minus2:
    for codon in codonUsage_minus2[mRNA]:
      codonCount_minus2[codon] += codonUsage_minus2[mRNA][codon]
  
  for mRNA in codonUsage_minus3:
    for codon in codonUsage_minus3[mRNA]:
      codonCount_minus3[codon] += codonUsage_minus3[mRNA][codon]
  
  for mRNA in codonUsage_plus1:
    for codon in codonUsage_plus1[mRNA]:
      codonCount_plus1[codon] += codonUsage_plus1[mRNA][codon]
  
  for mRNA in codonUsage_plus2:
    for codon in codonUsage_plus2[mRNA]:
      codonCount_plus2[codon] += codonUsage_plus2[mRNA][codon]
  
  for mRNA in codonUsage_plus3:
    for codon in codonUsage_plus3[mRNA]:
      codonCount_plus3[codon] += codonUsage_plus3[mRNA][codon]
  
  for mRNA in codonUsage_1first:
    for codon in codonUsage_1first[mRNA]:
      codonCount_1first[codon] += codonUsage_1first[mRNA][codon]
  
  for mRNA in codonUsage_2first:
    for codon in codonUsage_2first[mRNA]:
      codonCount_2first[codon] += codonUsage_2first[mRNA][codon]
  
  for mRNA in codonUsage_3first:
    for codon in codonUsage_3first[mRNA]:
      codonCount_3first[codon] += codonUsage_3first[mRNA][codon]
  
  for mRNA in codonUsage_1last:
    for codon in codonUsage_1last[mRNA]:
      codonCount_1last[codon] += codonUsage_1last[mRNA][codon]
  
  for mRNA in codonUsage_2last:
    for codon in codonUsage_2last[mRNA]:
      codonCount_2last[codon] += codonUsage_2last[mRNA][codon]
  
  for mRNA in codonUsage_3last:
    for codon in codonUsage_3last[mRNA]:
      codonCount_3last[codon] += codonUsage_3last[mRNA][codon]
  
  
  
  
  """
  Table output of codon occupancy
  """
  table=np.array([
                 [codonCount_minus3[x] for x in fullcodonlist],
                 [codonCount_minus2[x] for x in fullcodonlist],
                 [codonCount_minus1[x] for x in fullcodonlist],
                 [codonCount_A_site[x] for x in fullcodonlist],
                 [codonCount_P_site[x] for x in fullcodonlist],
                 [codonCount_E_site[x] for x in fullcodonlist],
                 [codonCount_plus1[x] for x in fullcodonlist],
                 [codonCount_plus2[x] for x in fullcodonlist],
                 [codonCount_plus3[x] for x in fullcodonlist]
                 ])
  tableT=table.T
  
  with open(outFile+"_BCO", 'w') as fout :
    fout.write("\tup3\tup2\tup1\tA\tP\tE\tdo1\tdo2\tdo3\n")
    #for row in table.T: 
    for i in range(np.shape(table.T)[0]) :
      row=tableT[i,:]
      currentcodon=fullcodonlist[i]
      fout.write( currentcodon + '\t')
      for elt in row :
        fout.write( str(elt) + '\t' )
      
      fout.write( '\n' )
  
  """
  Table output of codon count at footprint start, footprint end
  """
  
  tabel=np.array([
                 [codonCount_1first[x] for x in fullcodonlist],
                 [codonCount_2first[x] for x in fullcodonlist],
                 [codonCount_3first[x] for x in fullcodonlist],
                 [codonCount_1last[x] for x in fullcodonlist],
                 [codonCount_2last[x] for x in fullcodonlist],
                 [codonCount_3last[x] for x in fullcodonlist]
                 ])
  tabelT=tabel.T
  
  with open(outFile+'.limCod', 'w') as fout :
    fout.write("\tcodon1\tcodon2\tcodon3\tcodon3tolast\tcodon2tolast\tlastcodon\n")
    #for row in tabel.T: 
    for i in range(np.shape(tabel.T)[0]) :
      row=tabelT[i,:]
      currentcodon=fullcodonlist[i]
      fout.write( currentcodon + '\t')
      for elt in row :
        fout.write( str(elt) + '\t' )
      
      fout.write( '\n' )
  
  
  """
  Table output of nt count at footprint start, footprint end
  """
  ntlist  = ['A','C','G','T']
  ntlistU = ['A','C','G','U']
  bleta=np.array([
                 [ntCount_1[x] for x in ntlist],
                 [ntCount_2[x] for x in ntlist],
                 [ntCount_3[x] for x in ntlist],
                 [ntCount_4[x] for x in ntlist],
                 [ntCount_5[x] for x in ntlist],
                 [ntCount_6[x] for x in ntlist],
                 [ntCount_7[x] for x in ntlist],
                 [ntCount_8[x] for x in ntlist],
                 [ntCount_9[x] for x in ntlist],
                 [ntCount_minus9[x] for x in ntlist],
                 [ntCount_minus8[x] for x in ntlist],
                 [ntCount_minus7[x] for x in ntlist],
                 [ntCount_minus6[x] for x in ntlist],
                 [ntCount_minus5[x] for x in ntlist],
                 [ntCount_minus4[x] for x in ntlist],
                 [ntCount_minus3[x] for x in ntlist],
                 [ntCount_minus2[x] for x in ntlist],
                 [ntCount_minus1[x] for x in ntlist]
                 ])
  bletaT=bleta.T
  
  with open(outFile+'.limNt', 'w') as fout :
    fout.write("\tnt1\tnt2\tnt3\tnt4\tnt5\tnt6\tnt7\tnt8\tnt9\tnt9thfromend\tnt8thfromend\tnt7thfromend\tnt6thfromend\tnt5thfromend\tnt4thfromend\tnt3rdfromend\tnt2ndfromend\tnt1stfromend\n")
    #for row in bleta.T: 
    for i in range(np.shape(bleta.T)[0]) :
      row=bletaT[i,:]
      currentNt=ntlistU[i]
      fout.write( currentNt + '\t')
      for elt in row :
        fout.write( str(elt) + '\t' )
      
      fout.write( '\n' )
  
  
  
  """ SCO and read counts files """                                                                 
  fileSCO = open(outFile+'.SCO', 'w')
  fileSCO.write("# (mRNA,codon) occupancy by mRNA ID, by increasing codonPosition\n")      
  fileSCO.write("mRNA\tcodonPosition\tcodon\tcoverage\tsingleCodonOccupancy\t%codingSequenceCovered\n")
  #
  fileNbreads = open(outFile+'.nbreads', 'w')
  fileNbreads.write("mRNA\tread.count\n")
  #
  for mRNA in SingleCodonOccupancy:
    #print mRNA
    positions = [position for position in SingleCodonOccupancy[mRNA]]
    positions.sort()
    # marker for filter : at least 1 one in 2 codons should be covered
    continuity = 100.*len(positions) / float(max(1,len(codonLists[mRNA])-15))
    # print infor for ordered positions
    for position in positions:
      fileSCO.write(mRNA+'\t'+
                    str(position)+'\t'+
                    str(codonLists[mRNA][position])+'\t'+
                    str(codonCover_A_site[mRNA][position])+'\t'+
                    str(round(SingleCodonOccupancy[mRNA][position],4))+'\t'+
                    str(round(continuity,2))+'\n')
    
    #
    nbreadsRNA = sum([codonCover_A_site[mRNA][position] for position in positions])
    fileNbreads.write(mRNA+"\t"+str(nbreadsRNA)+"\n")
  
  fileSCO.close()
  fileNbreads.close()
  
  
  """ RPKM and Metagene files """
  filerpkm = open(outFile+'.RPKM', 'w')
  filerpkm.write("mRNA\tRPKM\n")
  #
  filemeta = open(outFile+'.metagene', 'w')
  filemeta.write("mRNA\tpos-nt-in-A\tpos-metagene\tcoverage\n")
  #nt pos in A site -oo->0 5'UTR, 0->oo CDS and 3'UTR
  #standardized positions [-1-0[: 5'UTR, [0-1[: CDS, [1-2]: 3'UTR 
  
  for mRNA in dico_mRNA_POS_nbReads:
    positions = [position for position in dico_mRNA_POS_nbReads[mRNA]]
    positions.sort()
    # RPKM (Nature Methods 5:621-628, 2008, Motarzavi et al)
    coveragemRNA = sum([dico_mRNA_POS_nbReads[mRNA][position] for position in positions])
    localStart = localDic[mRNA][0]
    localEnd = localDic[mRNA][1]
    transcriptlength = localDic[mRNA][3]
    kilobaseExon = float(localEnd-localStart) / float(1e3)
    if (kilobaseExon * millionMapped) > 0. :
      rpkm = coveragemRNA / kilobaseExon / millionMapped
    elif verbose :
      print "Warning : kilobaseExon or millionMapped is non positive, calculation of rpkm is impacted."
      print "          kilobaseExon  = "+str(kilobaseExon)
      print "          millionMapped =  "+str(millionMapped)
    
    # print RPKM
    filerpkm.write(mRNA+'\t'+str(rpkm)+'\n')
    # print metagene information
    for i in range(len(positions)) :
      #
      pos = positions[i]
      #
      if pos<0:
        if localStart > 0 :
          pos_std=round(3*pos/float(localStart),4)
        elif verbose :
          print "Warning : localStart is non positive, calculation of pos_std (metagene) is impacted."
          print "          localStart  = "+str(localStart)
          print "          mRNA        = "+str(mRNA)
      elif pos<(localEnd-localStart):
        if (localEnd-localStart) > 0 :
          pos_std=round(3*pos/float(localEnd-localStart),4)
        elif verbose :
          print "Warning : (localEnd-localStart) is non positive, calculation of pos_std (metagene) is impacted."
          print "          (localEnd-localStart)  = "+str((localEnd-localStart))
          print "          mRNA                   = "+str(mRNA)
      else:
        if (transcriptlength-localEnd) > 0 :
          pos_std=round((3*pos-(localEnd-localStart))/float(transcriptlength-localEnd),4)
        elif verbose :
          print "Warning : (transcriptlength-localEnd) is non positive, calculation of pos_std (metagene) is impacted."
          print "          (transcriptlength-localEnd)  = "+str((transcriptlength-localEnd))
          print "          mRNA                         = "+str(mRNA)
      #
      pos_cov = dico_mRNA_POS_nbReads[mRNA][pos]
      #
      filemeta.write(str(mRNA) + '\t' + str(pos) + '\t' + str(pos_std) + '\t' + str(pos_cov) + '\n')
  
  
  
  filerpkm.close()
  #
  filemeta.close()
  
  
  """ Paused Sites """                                                     #change v10/v9: suppl calc and output
  np.savetxt(outFile+'.paused', np.array(lists_mRNA_paused_sites), fmt='%s', delimiter='\t')
  
  
  
  # print "Codon usage calculation and output (could last a little while)"
  # 
  # """
  # Codon Usage in translated CDS
  # - each mRNA counted once
  # - x nb_reads
  # # Counter of codons actually translated, quantitative
  # # codonUsageQuantitativeCDS[codon] = Sum_aligned.read-1-..-N Nb.codons.in.CDS[15:-5]
  # """
  # 
  # codonUsageQuantitativeCDS = Counter([])
  # 
  # for mRNA in codonUsage_A_site:
  #    # sum over positions of SingleCodonOccupancy (i.e. between AUG+15codons and Stop-5codons),
  #    # but coverage is contained in codonCover, not in SingleCodonOcc.
  #    nbCDS = sum([codonCover_A_site[mRNA][position] for position in SingleCodonOccupancy[mRNA]])
  #    codonUsageQuantitativeCDS = codonUsageQuantitativeCDS + Counter(nbCDS*codonLists[mRNA][15:-5])
  #    
  # 
  # """ Write Codon Usage, each mRNA counted only once to File """
  # array_codonUsageQuantitativeCDS = np.array([[key,val] for (key,val) in sorted(codonUsageQuantitativeCDS.iteritems())])  
  # #array_codonUsageQuantitativeCDS = np.array([[key,val] for (key,val) in codonUsageQuantitativeCDS.iteritems()])  
  # # save to file
  # np.savetxt(outFile+'.codonusage', array_codonUsageQuantitativeCDS, fmt='%s', delimiter='\t')
  
  
  
  """ Codon coverage and enrichment in -ex3 (3') to +ex5 (5') to repr. Hussman et al. 2015 and further. """
  """     OBSOLETE note : iCodon Above is -ex5 at 5' and +ex3 at 3'                                              """
  """     OBSOLETE       conversion to the opposite is done here via 'colKeys' and flag 'opposite=True'         """
  # arguments for function 'write_dico_double_to_file'
  rowkeys = fullcodonlist
  colkeys = range(-ex3,ex5+1)
  
  # Write codon coverage to file
  #   using function call : write_dico_double_to_file(dico_double, keys1, keys2, outFile, suffix) 
  write_dico_double_to_file(meancodonEnrichmnt_H,     rowkeys, colkeys, outFile, '_Hussmann_mean'    , opposite=False) #True) #prev. codonCover
  write_dico_double_to_file(sdevcodonEnrichmnt_H,     rowkeys, colkeys, outFile, '_Hussmann_sd'      , opposite=False) #True) #prev. codonCover
  
  # Write codon enrichment to file
  write_dico_double_to_file(meancodonEnrichmnt,     rowkeys, colkeys, outFile, '_Codon-enrichment-unbiased_mean', opposite=False) #True) #prev CodonEnrichmnt
  write_dico_double_to_file(sdevcodonEnrichmnt,     rowkeys, colkeys, outFile, '_Codon-enrichment-unbiased_sd'  , opposite=False) #True) #prev CodonEnrichmnt
  
  
  """
  The End
  """ 
  if verbose : print "End"


