#!/usr/bin/python

# Read gtf file, derive CDS annotation, write to table file

import gzip

def gtf2table(gtfFile,
              outFile, 
              endinside1 = 1, 
              stopoutside1 = 1, 
              stopminusoutside1 = 1,
              verbose=False):
  """
  Read GTF annotation file in input, for example :
  1	havana	transcript	65419	71585	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	tag	basic						
  1	havana	exon	65419	65433	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	exon_number	1	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	exon_id	ENSE00003812156	exon_version	1	tag	basic
  1	havana	exon	65520	65573	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	exon_number	2	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	exon_id	ENSE00003813641	exon_version	1	tag	basic
  1	havana	CDS	65565	65573	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	exon_number	2	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	protein_id	ENSP00000493376	protein_version	2	tag	basic
  1	havana	start_codon	65565	65567	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	exon_number	2	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	tag	basic				
  1	havana	exon	69037	71585	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	exon_number	3	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	exon_id	ENSE00003813949	exon_version	1	tag	basic
  1	havana	CDS	69037	70005	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	exon_number	3	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	protein_id	ENSP00000493376	protein_version	2	tag	basic
  1	havana	stop_codon	70006	70008	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	exon_number	3	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	tag	basic				
  1	havana	five_prime_utr	65419	65433	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	tag	basic						
  1	havana	five_prime_utr	65520	65564	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	tag	basic						
  1	havana	three_prime_utr	70009	71585	0	+	0	gene_id	ENSG00000186092	gene_version	6	transcript_id	ENST00000641515	transcript_version	2	gene_name	OR4F5	gene_source	ensembl_havana	gene_biotype	protein_coding	transcript_name	OR4F5-202	transcript_source	havana	transcript_biotype	protein_coding	tag	basic						

  To table annotation for CDS start/stop and length, for example :
  ID	                localStart	localEnd	strandDir Exonlen
  ENST00000641515     59          1034      +                       # h. sapiens
  SPAC110.03	        54	        633	      +                       # s. pombe
  YAL008W	            0	          597	      +	        597           # s. cerevisiae

  
  'global' means here: 'relative to genomic coordinates'
  
  'local'              'relative to transcript coordinates' (ie position relative to transcript's first nt 
                        (whose position is 0, whereas last nt in transcript has position transcript length-1).
  
  'localStart'          : Coding sequence start, in transcript coordinates, STARTING AT 0.
                           (START codon should correspond to range [localStart:localStart+3])
                           localStart = "start" field from "start_codon" feature
                           
  'localEnd'            : Coding sequence   end, in transcript coordinates, 1st base AFTER STOP codon.
                           (STOP  codon should correspond to range [localEnd-3:localEnd])
                           localEnd = 1+"stop" field from "stop_codon" feature
                           or
                           localEnd = "start" field from "three_prime_utr" feature
  
  Parameters :
  ----------
                       
  endinside1            : 0 if field end (5th) is first base outside of feature, in GTF file,
                          1                        last base in feature (exon, CDS, transcript, etc.).
                          
  stopoutside1          : 0 if stop is included in CDS, in GTF file, in '+' strand case.
                          1            outside the CDS as annotated in GTF file, in '+' strand case.
                          
  stopminusoutside1     : 0 if stop is included in CDS, in GTF file, in '-' strand case.
                          1            outside the CDS as annotated in GTF file, in '-' strand case.
                          
  """
  gtfDic = {}
  localDic = {}
  
  """ Parse GTF into dictionary """
  
  # detect gzipped file
  if gtfFile[-3:]==".gz" or gtfFile[-5:]==".gzip" :
    gtfSource = gzip.open(gtfFile,'rb')
  else :
    gtfSource = open(gtfFile,'rb')
  
  for row in gtfSource:    
    
    # Consider row content only if not a comment
    if row[0] != "#" :
      
      if verbose : print row
      
      fields    = row.split('\t')
      feature   = fields[2]
      attribute = fields[8]
      
      # protein coding ?
      IsProteinCoding = "protein_coding" in fields[8]
      
      if IsProteinCoding and feature in ["transcript", "exon","CDS"] :
        if verbose : print "if IsProteinCoding and feature in ['transcript', 'exon','CDS']"
        # Attributes as a dictionary
        attributes = attribute.split(';')
        # remove terminal newline if any
        if attributes[-1] == "\n" :
          attributes = attributes[:-1]
        
        if verbose : print "Starting parsing of row, since protein-coding and applicable feature."
        # Parsing attribute subfields (not very robust)
        #    Splitting each field according to blanks
        at2 = [(at.split(" ")) for at in attributes]
        #    Remove non-informative subfields
        at3 = [filter(None, x) for x in at2]
        #    Clean elements inside of subfields
        at4 = [[x.replace(' ','').replace('"','') for x in l] for l in at3 ]
        #    Remove elements besides 2 (should not change anything - if it does, 
        #                       it means something went wrong with the parsing) 
        at5 = [x[:2] for x in at4]
        
        if verbose : print "Update of attributesDict with parsed row = at5 =  ", at5
        attributesDict = dict(at5)
        
        # Transcript ID
        if verbose : print "transcriptID"
        transcriptId = attributesDict['transcript_id'] #fields[8].split('transcript_id "')[1]
        # Initialize if not yet done
        if transcriptId not in gtfDic :
          gtfDic[transcriptId] = {}
        
        if feature == "transcript" :
          if verbose : print "  if feature=transcript"
          # transcript's start and stop
          trStart = fields[3]
          trStop  = fields[4]
          gtfDic[transcriptId]['transcriptStartStop'] = [trStart, trStop]
          # stranddir
          strDir = fields[6]
          gtfDic[transcriptId]['strandDir'] = strDir
        
        if feature == "CDS" :
          # transcript's start and stop
          cdsStarti = fields[3]
          cdsStopi  = fields[4]
          if 'CDS' not in gtfDic[transcriptId] :
            gtfDic[transcriptId]['CDS'] = [cdsStarti, cdsStopi]
          else :
            cdsStart = str(min( [int(gtfDic[transcriptId]['CDS'][0]), int(cdsStarti)]))
            cdsStop  = str(max( [int(gtfDic[transcriptId]['CDS'][1]), int(cdsStopi)] ))  #gtfDic[transcriptId]['CDS'].extend([cdsStarti, cdsStopi]
            gtfDic[transcriptId]['CDS'] = [cdsStart, cdsStop]
        
        if feature == "exon" :
          
          # exon's start and stop
          exonStarti = fields[3]
          exonStopi  = fields[4]
          
          if 'exonStarts' not in gtfDic[transcriptId] :
            gtfDic[transcriptId]['exonStarts'] = [exonStarti]
            gtfDic[transcriptId]['exonStops']  = [exonStopi]
          else :
            gtfDic[transcriptId]['exonStarts'].append( exonStarti )
            gtfDic[transcriptId]['exonStops'].append(  exonStopi )




  
  """ Use gtf-dictionary to derive local end / start and length on transcript and write to file,
          using local dictionary containing gene ID as key with local start and end sites, strand direction as tuple.
          Note :  1st base=0, and STOP is at [localEnd-3:localEnd]"""
  #with open(refGenes,'rb') as source:  
  with open(outFile, 'w') as fichier :
    fichier.write("ID\tlocalStart\tlocalEnd\tstrandDir\texonLength\n")
    #print "transcriptID strandDir CDSstart CDSlength"
    for transcript in gtfDic :        #for rows in source:
      #elements = rows.split('\t')
      transcriptID = transcript    #str(elements[0])  # NCBI RefSeq ID of each canonical transcript. 
      strandDir = gtfDic[transcriptID]['strandDir'] #str(elements[2]) # forward or reverse strand
      startPos = int(gtfDic[transcriptID]['transcriptStartStop'][0]) #int(elements[3])  # global transcript start position
      if verbose: print "startPos: ",startPos
      endPos = int(gtfDic[transcriptID]['transcriptStartStop'][1]) #int(elements[4])   # global transcript end position
      startCDS = int(gtfDic[transcriptID]['CDS'][0]) #int(elements[5]) # global start position of the coding sequence
      if verbose: print "startCDS: ",startCDS
      endCDS = int(gtfDic[transcriptID]['CDS'][1]) #int(elements[6]) #global end position of the coding sequence
      
      #new coordinates of the exon positions:
      exonStartPositions = gtfDic[transcriptID]['exonStarts'] #elements[7].split(',') # list of global exon start positions
      if verbose: print "exonStartPositions: ",exonStartPositions
      #exonStartPositions = exonStartPositions[:-1] # getting rid of the extra space at the end. So that I can turn the elements of the list into integers.
      exonStartPositions = [int(i) for i in exonStartPositions]  # turning strings in a list into numbers.
      
      exonEndPositions = gtfDic[transcriptID]['exonStops'] #elements[8].split(',') # list of global exon end positions
      #exonEndPositions = exonEndPositions[:-1]
      exonEndPositions = [int(i) for i in exonEndPositions]
      
      localStart = 0
      localEnd = 0
      
      exonLength = 0 #adding up all the lengths of introns in the global transcript coordinates (to be used for getting local start and end sites of CDS)
      
      for start, end in zip(exonStartPositions, exonEndPositions): #start and end here indicate pairs of start/end sites for the exons. 
        # endinside1        : 0 if field end (5th) is first base outside of the feature, in GTF file
        #                       -> in this case, exon length = end-start and NOT end-start+1.
        #                     1                        last base of exon, CDS, transcript, etc.                
        #                       -> in this case, exon length = end-start+1.
        # stopoutside1      : 0 if stop is included in CDS, in GTF file
        #                       -> in this case, STOP is at positions [endCDS-3:endCDS] (+1, i.e. [endCDS-2:endCDS+1], if endinside1 = 1 )
        #                     1            outside the CDS as annotated in GTF file
        #                       -> in this case, STOP is at positions [endCDS:endCDS+3] (+1, i.e. [endCDS+1:endCDS+4], if endinside1 = 1 )
        # stopminusoutside1 : 0 if stop is included in CDS, in GTF file, in '-' strand case.
        #                     1            outside the CDS as annotated in GTF file, in '-' strand case.
        #                     -> same as stopoutside1 for the minus strand instead of the plus strand 
        
        exonLength += (end - start + endinside1) #exon length is the cumulative sum of the exons for each transcript
        #                            
        if strandDir == '+':
           # Calculate local start site :
           if (start <= startCDS) and (end < startCDS):   #exon is in the 5' UTR (start and end before CDS start site)
             localStart += (end - start + endinside1)       #add exon length to localStart position
             
           elif (start <= startCDS) and (end > startCDS): #CDS start site is within the current exon. 
             localStart += (startCDS - start)               #add {in CDS}-exon length #CL nota bene: 
             # startCDS-start and NOT startCDS-start+1, because python indices start at 0.
             
           else:
             pass
                
           # Calculate local end site:
           if (end <= endCDS):                            #exon is before the end of the CDS. 
             localEnd += (end - start + endinside1)         #add exon length to localEnd position
              
           elif (end >= endCDS) and (start < endCDS):     #if the exon ends after CDS end-site but starts within the CDS. 
             localEnd += (endCDS + 3*stopoutside1 - start + endinside1) 
             
           else:
             pass
           
        elif strandDir == '-':
           
           #NotaBene: fasta sequences are in start >> stop order, which is reversed wrt exon coordinates on the gene. 
           
           #calculate local start-site:
           if (start <= (endCDS)) and (end >= endCDS):      #exon contains endCDS
             # which might mean either : - start site A is at endCDS (endinside1 = 1)
             # or                      : - start site A is in previous exon (endinside1 = 0)
             # => procedure works correctly in both cases
             localStart += (end - endCDS)  # CL nota bene: endCDS is the start, and NOT start+1, 
             #                               because python indexing starts at 0.   
             
           elif (start >= endCDS):                          #exon is in the 5'UTR
             localStart += (end - start + endinside1)
             
           else:
             pass
           
           #calculate local end-site:
           if (start >= startCDS):                          #exon is in 5' UTR or CDS
             localEnd += (end - start + endinside1)
             
           elif (start < startCDS) and (end >= startCDS):   #exon contains startCDS
             localEnd += (end - startCDS + endinside1 + 3*stopminusoutside1)
             
           else:
             pass
         
      if verbose: print "localStart: ",localStart
      if verbose: print "localEnd: ",localEnd
      if verbose: print "strandDir: ",strandDir
      localDic[transcriptID] = (localStart, localEnd, strandDir, exonLength) #There are 28661 transcripts in the dictionary of which, all of the 28661 of them have the length of their CDS (localEnd - localStart) divisible by 3. 
      #print transcriptID,strandDir,(localStart+1),(localEnd-localStart)
      #fichier.write(str(transcriptID)+'\t'+str(strandDir)+'\t'+str((localStart+1))+'\t'+str((localEnd-localStart))+'\n')
      fichier.write(str(transcriptID)+'\t'+str(localStart)+'\t'+str(localEnd)+'\t'+str(strandDir)+'\t'+str(exonLength)+'\n')
      #fichier.write('\t'.join(map(str, positions))+'\n')

### END of gtf2table

             
             
             
             
             
