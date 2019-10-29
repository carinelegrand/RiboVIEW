#!/usr/bin/python
import pysam 
from collections import defaultdict
from Bio import SeqIO
import numpy as np

def periodicity(bam,refCDS,refSeq,outFile, versionStrip=False, messages=True):
  """
  Metagene analysis at the start site to see if coverage has 3-nt periodicity.
  
  in input :
                          
  bam          type=str  aligned reads (BAM file) 
  refCDS       type=str  CDS annotation (output of bed2table)
  refSeq       type=str  reference transcript sequences (FASTA)
  versionStrip type=bool indicate if .1, .2, etc should be removed from transcript names                         defaults to False
  
  in output :
  
  outFile type=str  count of read starts per position around AUG, for each length
  
  """
  
  verbose=False
  
  if messages : print "Read alignment file (sam/bam)" 
  
  #aligned file ("rb" for .bam file and "r" for .sam file):
  samfile = pysam.Samfile(bam, "rb")
  
  if messages : print "Read CDS annotation file (tab-separated table)"
  
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
      localDic[elements[0]] = (int(elements[1]), int(elements[2]), elements[3], int(elements[4]))
    
  if messages : print "Read reference sequences (fasta)"     
 
  #################################
  # Database of reference sequences, from FASTA file, filtered for sound sequences
  #    - first codon of CDS must be ATG
  #    - last codon of CDS must be TAA, TAG or TGA

  #reference transcript file:
  transcript = SeqIO.to_dict(SeqIO.parse(open(refSeq), "fasta")) #function seqIO.to_dict() turns a SeqRecord iterator (list) into a dictionary (in memory). 

  codonLists = {}
  for mRNA in localDic:
      #print mRNA
      # sequence
      if mRNA in transcript:
        sequence = str(transcript[mRNA].seq)
        sequence = sequence.lower().replace('u','t')
        #print sequence
        codonList = []
        #if mRNA in localDic:
        for i in range(localDic[mRNA][0], localDic[mRNA][1], 3): #here, I need to make sure that scope of the variables are correctly taken into account each time. From start to the end of the mRNA (using local coordinates) in increments of 3.
          #  {localStart of mRNA, localStart of mRNA +3, localStart of mRNA +6, ..., localEnd of mRNA -3}
          if i + 3 < localDic[mRNA][1]:
            codon = sequence[i:i+3]
          else:
            codon = sequence[i:localDic[mRNA][1]]
          codonList.append(codon)
        #print codonList
        """
        if sensible coding sequences (3*n long, beginning with AUG, ending with TAA, TAG or TGA):
           codon Usage at each site for this CDS/target/mRNA is recorded
        """
        if codonList != []:
          #print (localDic[mRNA][1] - localDic[mRNA][0])
          #print (localDic[mRNA][1] - localDic[mRNA][0]) % 3
          #print codonList[0]
          #print codonList[-1]
          if ((localDic[mRNA][1] - localDic[mRNA][0]) % 3 == 0) and (codonList[0] == "atg"): #check if the transcript is divisible by 3 and starts with an AUG. 
            #print mRNA 3 != 0):# and (codonList[0] == "atg"): #check if the transcript is divisible by 3 and starts with an AUG. 
            if (codonList[-1] == 'taa') or (codonList[-1] == 'tag') or (codonList[-1] == 'tga'):#only those that also end at a stop codon
              #set(codonList) #unique elements
              codonLists[mRNA] = codonList
  
  # Nb of transcripts, Nb which have a sound sequence
  if messages : 
    print
    print "There are ", len(transcript), "transcripts, among"
    print "  which",len(codonLists), " have a correct annotation of start and stop codons."
    print ""
       
  #################################
  # Count reads starting at local start +/- 20
  
  #dictionary containing coverage
  coverStart = defaultdict(lambda:defaultdict(int))
  
  #for each aligned read:
  # - get reference sequence name (refid)
  # - if this reference sequence corresponds to reference genes (localDic):
  #     - determine read's start relatively to the CDS's start
  #     - if this start is in [-20;20], then increment coverStart @ read length
  for alignread in samfile:
    if alignread.tid >= 1:
      refid = samfile.getrname(alignread.tid)

      # keep accession only, not the further information
      refid = refid.split('|')[0]
      
      # keep or remove version number
      if versionStrip:
        #transcript version nb to rm in case transcript names in .bed do not contain it.
        refid = refid.split('.')[0]
      
      # offset length, filter 
      if refid in codonLists:
        cdspos = localDic[refid][0] #start position
        #cdsend = localDic[refid][1] #end position
        if (alignread.pos != None) and (24 <= alignread.alen <= 34):#alignread.alen gives the length of the reads in the BAM file. 
          if 1: #fully unnecessary: if not alignread.is_reverse: ###only for the forward strand. Have to do it also for thereverse strand ###
            startpos =  alignread.pos - cdspos # position relative to the start
            #endpos = cdsend - alignread.pos #position relative to the end
            if -20 <= startpos <= 20: 
              coverStart[len(alignread.seq)][startpos] += 1 #a dictionary containing the sum of the start of the reads relative to the start position (-20 to +20) of the genes and stratified for each length. 
              
              #coverEnd[len(alignread.seq)][]
  ###############
  # Table output
  table=np.array([
                 range(-20,21),
                 [coverStart[25][x] for x in range(-20,21)],
                 [coverStart[26][x] for x in range(-20,21)],
                 [coverStart[27][x] for x in range(-20,21)],
                 [coverStart[28][x] for x in range(-20,21)],
                 [coverStart[29][x] for x in range(-20,21)],
                 [coverStart[30][x] for x in range(-20,21)],
                 [coverStart[31][x] for x in range(-20,21)],
                 [coverStart[32][x] for x in range(-20,21)]
                 ])
  
  if messages : 
    print "Coverage around AUG is written in : "
    print outFile
    print ""
  
  with open(outFile, 'w') as fout :
    fout.write("pos\len\t25\t26\t27\t28\t29\t30\t31\t32\n")
    for row in table.T: 
      for elt in row :
        fout.write( str(elt)+'\t' )
      
      fout.write( '\n' )

  
    

