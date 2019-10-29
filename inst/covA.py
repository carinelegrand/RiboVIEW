#!/usr/bin/python
from collections import defaultdict

def covA(inSCOfull, codRes, outSCObin):
  """
  Coverage by mRNA and codon position, binned according to codRes 
  
    in input :                              
      inSCOfull type=str  Input file "single codon occupancy" (SCO) for one sample
      codRes  type=int  codon resolution for SCO
    
    in output :      
      outSCObin type=str  Output file coverage by mRNA and binned position
  """ 
  
  ### Parameters
  verbose=False

  ### Read coverage by mRNA and codon position
  SCOcodres = defaultdict(lambda:defaultdict(int))
  #
  with open(inSCOfull, mode='r') as f :
    commnt=f.readline()
    header=f.readline()
    #Further content of inSCOfull :  (mRNA, codonPosition, codon, coverage, SingleCodonOccupancy, pcCDScov)
    
    for line in f:
      elements=line.strip('\n').split('\t')
      if verbose : print elements 
      ## Bin acc. to codRes ; Note : AUG has codon position 0, 1st codon after AUG has position 1, etc.
      #Note : /codRes is an integer division, while *codRes is a normal multiplication
      codonPositionBin = int(elements[1]) /int(codRes) *codRes
      #Read into dict.   :  SCOcodres[mRNA][codonPositionBin] =  coverage
      if codonPositionBin not in SCOcodres[elements[0]] :
        SCOcodres[elements[0]][codonPositionBin] = elements[3]
      else :
        SCOcodres[elements[0]][codonPositionBin] += elements[3]
  
  ### Write to file
  with open(outSCObin, 'w') as fileCovA :
    fileCovA.write("mRNA\tcodonPosition\tcoverage\n") 
    #
    for mRNA in SCOcodres :
      for CPB in SCOcodres[mRNA] :
        fileCovA.write( mRNA                     +'\t'+
                        str(CPB)                 +'\t'+
                        str(SCOcodres[mRNA][CPB])+'\n')
  
  ### End

