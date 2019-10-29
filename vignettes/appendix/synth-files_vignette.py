#!/usr/bin/python
#
# Purpose : generate fasta, tsv file and bam files for package-included 
#   synthetic files.
#
# Setting : Codon Proline CCC enriched (slower) in condition 2 : 
#             sampling probability *2 for this codon. 
#
#
import pysam 
import os
import random

## Parameters
NbReads = 20000
filedir  = "./synth-files-vignette/"
# Create input directory
os.system("mkdir -p "+filedir)
tsvFile = filedir+"synth.tsv" 
faFile  = filedir+"synth.fasta"
seed    = 7234
# Set random seed
random.seed(seed)
Nmrna  = 100
# t and not u because samtools considers only t.
bases  = ['a','c','g','t']
# Codons to sample from
codons = ['aaa','aac','aag','aat','aca','acc','acg','act','aga','agc',
          'agg','agt','ata','atc','atg','att','caa','cac','cag','cat',
          'cca','ccc','ccg','cct','cga','cgc','cgg','cgt','cta','ctc',
          'ctg','ctt','gaa','gac','gag','gat','gca','gcc','gcg','gct',
          'gga','ggc','ggg','ggt','gta','gtc','gtg','gtt','tac','tat',
          'tca','tcc','tcg','tct','tgc','tgg','tgt','tta','ttc','ttg',
          'ttt']
stops  = ['taa','tag','tga']    
NbCond = 2
NbRep  = 3

## Initialize synthetic mRNAs in virtual table name_i,length_i,sequence_i
#  Note : no UTRs and here lenNt = lenCDS (+3nt of stop codon)
mRNAname     = [None] * Nmrna
CDSlenC      = [None] * Nmrna
CDSlenNt     = [None] * Nmrna
mRNAnt5UTR   = [None] * Nmrna
mRNAnt3UTR   = [None] * Nmrna
mRNAlenNt    = [None] * Nmrna
mRNAseq      = [None] * Nmrna
Pos_CCC      = [None] * Nmrna
outFA        = open(faFile, 'w')
outTSV       = open(tsvFile, 'w')
outTSV.write("ID\tlocalStart\tlocalEnd\tstrandDir\texonLength\n")

## Generate synthetic mRNAs
for i in range(Nmrna) :
    
    mRNAname[i]    = "mRNA"+str(i)
    CDSlenC[i]     = random.randint(50,500)
    CDSlenNt[i]    = 3*CDSlenC[i]
    mRNAnt5UTR[i]  = random.randint(20,200)
    mRNAnt3UTR[i]  = random.randint(0,100)
    mRNAlenNt[i]   = mRNAnt5UTR[i] + CDSlenNt[i] + mRNAnt3UTR[i]
    
    # Elements for sequence construction
    mRNA_Pre5      = ''.join([random.choice(bases) for ii in range(mRNAnt5UTR[i])])
    mRNA_PreC0     = [random.choice(codons) for ii in range(CDSlenC[i]-2)]
    mRNA_PreC      = ''.join(mRNA_PreC0)
    mRNA_Stop      = random.choice(stops)
    mRNA_Pre3      = ''.join([random.choice(bases) for ii in range(mRNAnt3UTR[i])])
    
    # Construct a sequence
    mRNAseq[i]     = mRNA_Pre5 + 'atg' + mRNA_PreC + mRNA_Stop + mRNA_Pre3
    # ii+1 since AUG is not yet in mRNA_PreC
    # (ii+1)<(CDSlenC[i]-15) to yield valid footprints later on
    Pos_CCC[i]     = [ii+1 for ii,x in enumerate(mRNA_PreC) if x=="ccc" if (ii+1)<(CDSlenC[i]-15)]
    
    # Write entry to fasta
    outFA.write(">"+mRNAname[i]+"\n")
    for j in range(mRNAlenNt[i]/100) :
        outFA.write(mRNAseq[i][(j*100):(j*100+100)]+"\n")
    
    j=j+1
    outFA.write(mRNAseq[i][(j*100):]+"\n")
    
    # Write entry to TSV
    outTSV.write(mRNAname[i]+"\t"+str(mRNAnt5UTR[i])+"\t"+str(mRNAnt5UTR[i]+CDSlenNt[i])+"\t+\t"+str(mRNAlenNt[i])+"\n")

outFA.close()
outTSV.close()


## Generate SAM/BAM file for each condition and replicate
for Cond in range(1,1+NbCond) :
    for Rep in range(1,1+NbRep) :
        
        samFile = filedir+"Cond"+str(Cond)+"-Rep"+str(Rep)+".sam"
        bamFile = filedir+"Cond"+str(Cond)+"-Rep"+str(Rep)+".bam"
        
        ## Open SAM output file for writing
        with open(samFile, 'w') as fout :
            
            ## Create sam file header
            fout.write("@HD\tVN:1.0\tSO:unsorted\n")
            for i in range(Nmrna) :
                fout.write("@SQ	SN:"+mRNAname[i]+"	LN:"+str(mRNAlenNt[i])+"\n")
            
            fout.write("@PG	ID:Synthetic\n")
            
            ## Create sam file entries
            for i in range(NbReads) :
                imRNA = random.randint(0,Nmrna-1)
                
                Qname = "C"+str(Cond)+"-R"+str(Rep)+"."+str(i)
                Flag  = "0"
                Rname = mRNAname[imRNA]
                
                # Sample positions that would yield valid footprints : 
                #   pos cannot be in the last 35nt of CDS,
                #    - with periodicity ("1+3*random.randint(...)"),
                #    - locate in P-site if AUG, else in A (AP_Offset),
                #    - mostly in-frame ("random.choice([0]*50+...").
                
                # More weight (pause) on Pro-CCC, in Cond2
                k = 2
                indices = (range(CDSlenC[imRNA] - 15) + (k-1)*Pos_CCC[imRNA])
                Pos_codonInCDS_inequal_pre = random.choice(indices)
                Pos_codonInCDS_inequal  = 1 + 3 * Pos_codonInCDS_inequal_pre  #0,(CDSlenC[imRNA]-15))
                
                # Sample equally along CDS
                Pos_codonInCDS_equally = 1 + 3 * random.randint(0,(CDSlenC[imRNA]-15)) 
                
                # Attribute footprint reads with equal or unequal sampling depending on experimental condition
                if Cond==2 :
                    Pos_codonInCDS = Pos_codonInCDS_inequal
                else :
                    Pos_codonInCDS = Pos_codonInCDS_equally
                
                AP_Offset = -15
                if Pos_codonInCDS == 1 :
                    AP_Offset = -12
                
                # Elements for SAM/BAM file entry
                Pos = mRNAnt5UTR[imRNA] + Pos_codonInCDS + AP_Offset + random.choice([0]*50+[1]*10+[2]*1)
                Mapq  = "255"
                Rnext = "*"
                Pnext = "0"
                Tlen  = random.randint(25,32) #     Tlen  = "30"
                Cigar = str(Tlen)+"M"
                Seq   = mRNAseq[imRNA][(Pos-1):(Pos-1+Tlen)]
                Qual  = "".join(Tlen*["Z"])
                
                # Write to SAM
                fout.write(Qname+"\t"+Flag+"\t"+Rname+"\t"+str(Pos)+"\t"+Mapq+"\t"+
                           Cigar+"\t"+Rnext+"\t"+Pnext+"\t"+str(Tlen)+"\t"+Seq+"\t"+Qual+"\n")
        
        
        ## Convert sam to bam
        os.system("samtools view -hb "+samFile+" > "+bamFile)

