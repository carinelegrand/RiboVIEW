\name{repl.correl.codon}
\alias{repl.correl.codon}
\title{Correlation of footprint coverage at codon resolution}
\usage{
repl.correl.codon(list.bam, refCDS, refFASTA, mini, maxi, 
                    XP.names, XP.conditions, pathout)
}
\arguments{
    \item{list.bam }{List of bam files containing aligned reads for 
                                   each sample (same order as in XP.names)}
    \item{refCDS}{Address of file containing coding sequence annotation. This file
    should contain tab-separated values for mRNA name, nt position of start codon,
    nt position of first codon after stop codon, for example :  ID	localStart	localEnd
                                        NM_001276351	145	1348
                                        NM_001276352	145	799
                                        NM_000299	251	2495}
    \item{refFASTA}{Address of reference sequences for mRNA of the 
                                   studied organism}
    \item{mini}{Minimum footprint length to consider (as selected by user)}
    \item{maxi}{Maximum footprint length to consider (as selected by user)}
    \item{XP.names}{Vector of names for each sample}
    \item{XP.conditions}{Vector of experimental conditions for each sample}
    \item{pathout}{Address where output files will be written}
}
\value{
This function returns a list containing the following :
\describe{
    \item{cor}{Correlation between replicates}
    \item{plot}{Address of plot file in png format}
    \item{value}{Median of correlation across different conditions}
    \item{color}{Color white/orange/red corresponding to good/warning/poor level of quality}
    \item{recommendation}{Description and recommendation based on value}
  }
}
\description{
\code{repl.correl.codon} This function takes a list of samples and experimental conditions in input,
 calculates correlations at codon level between replicates.
}
\examples{
\donttest{

# Sequenced reads aligned to mRNA (and containing no rRNA, depleted previously),
#   in bam format
readsBAM.1.1  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                              "/extdata/Cond1-Rep1.bam",sep="")
readsBAM.1.2  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                              "/extdata/Cond1-Rep2.bam",sep="")
readsBAM.1.3  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                              "/extdata/Cond1-Rep3.bam",sep="")
readsBAM.2.1  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                              "/extdata/Cond2-Rep1.bam",sep="")
readsBAM.2.2  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                              "/extdata/Cond2-Rep2.bam",sep="")
readsBAM.2.3  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                              "/extdata/Cond2-Rep3.bam",sep="")

list.bam <- list(readsBAM.1.1, readsBAM.1.2, readsBAM.1.3, 
                 readsBAM.2.1, readsBAM.2.2, readsBAM.2.3)


#
## Experimental conditions, in text and as indicators :
#    0 for control
#    1 for a condition, treatment, case, etc...
#    2, 3, etc. for further conditions

XP.conditions   <- c("cond1","cond1","cond1","cond2", "cond2","cond2")
XP.conditions.i <- c( 1,1,1,2,2,2)
XP.names        <- c("C1.R1", "C1.R2", "C1.R3", 
                     "C2.R1", "C2.R2", "C2.R3")

#
## Reference annotation for mRNAs' CDS.
#

refCDS <- paste(system.file(package="RiboVIEW", mustWork = TRUE), "/extdata/synth.tsv", sep="")
# Note : CDS annotation can be obtained from a GTF file, 
#        using gtf2table(my-gtf-file, outfile = my-cds-file)
#        (for example GTF file as provided by Ensembl.org work well with gtf2table)

#
## Reference sequences for mRNAs.
#

refFASTA <- paste(system.file(package="RiboVIEW", mustWork = TRUE), "/extdata/synth.fasta", sep="")

#
## Work and output folder.
#

pathout  <-  paste(tempdir(),"/", sep="")
  ## !! This is a temporary directory, which will be erased when you leave R !!
  ##   For your own analyses you would probably prefer to point to a permanent repository :
  #      pathout <- /home/me/address-to-my-output-repository/ # Define address, 
  #                                                   #including a final slash.
  #      system(paste('mkdir',pathout)) # Create folder at said address.
  #      setwd(pathout)  # Go to this directory. This is useful if you want to 
  #                                         #save additional tables or figures.

# 
## A-site coverage periodicity by length
#

periodicity(list.bam, refCDS, refFASTA, pathout, XP.names, versionStrip = FALSE)

# 
## Select footprint length with sufficient periodicity
#

attach(listminmax <- select.FPlen(list.bam, pathout, XP.names))

#
## Codon occupancy, codon enrichment.
# 

enrichmentNoccupancy(list.bam, refCDS, refFASTA, mini, maxi, XP.names,  
                       pathout, versionStrip = FALSE)

#
## Replicates.
#

repl.correl.counts.Venn.res <- repl.correl.counts.Venn(XP.conditions, XP.names, 
                                                       pathout)
repl.correl.counts.Venn.res

repl.correl.gene.res <- repl.correl.gene(XP.conditions, XP.names, pathout)
repl.correl.gene.res

repl.correl.codon.res <- repl.correl.codon(list.bam, refCDS, refFASTA, 
                                           mini, maxi, 
                                           XP.names, XP.conditions, pathout)
repl.correl.codon.res

repl.correl.heatmap.res <- repl.correl.heatmap(XP.conditions.i, XP.names, pathout)
repl.correl.heatmap.res

}
\dontshow{

  # Sequenced reads aligned to mRNA (and containing no rRNA, depleted previously),
  #   in bam format
  readsBAM.1.1  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                                "/extdata/Cond1-Rep1.bam",sep="")
  readsBAM.1.2  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                                "/extdata/Cond1-Rep2.bam",sep="")
  readsBAM.1.3  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                                "/extdata/Cond1-Rep3.bam",sep="")
  readsBAM.2.1  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                                "/extdata/Cond2-Rep1.bam",sep="")
  readsBAM.2.2  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                                "/extdata/Cond2-Rep2.bam",sep="")
  readsBAM.2.3  <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                                                "/extdata/Cond2-Rep3.bam",sep="")

  list.bam <- list(readsBAM.1.1, readsBAM.1.2, readsBAM.1.3, 
                   readsBAM.2.1, readsBAM.2.2, readsBAM.2.3)


  #
  ## Experimental conditions, in text and as indicators :
  #    0 for control
  #    1 for a condition, treatment, case, etc...
  #    2, 3, etc. for further conditions

  XP.conditions   <- c("cond1","cond1","cond1","cond2", "cond2","cond2")
  XP.conditions.i <- c( 1,1,1,2,2,2)
  XP.names        <- c("C1.R1", "C1.R2", "C1.R3", 
                       "C2.R1", "C2.R2", "C2.R3")

  #
  ## Reference annotation for mRNAs' CDS.
  #

  refCDS <- paste(system.file(package="RiboVIEW", mustWork = TRUE), "/extdata/synth.tsv", sep="")
  # Note : CDS annotation can be obtained from a GTF file, 
  #        using gtf2table(my-gtf-file, outfile = my-cds-file)
  #        (for example GTF file as provided by Ensembl.org work well with gtf2table)

  #
  ## Reference sequences for mRNAs.
  #

  refFASTA <- paste(system.file(package="RiboVIEW", mustWork = TRUE), "/extdata/synth.fasta", sep="")

  #
  ## Work and output folder.
  #

  pathout  <-  paste(tempdir(),"/", sep="")
    ## !! This is a temporary directory, which will be erased when you leave R !!
    ##   For your own analyses you would probably prefer to point to a permanent repository :
    #      pathout <- /home/me/address-to-my-output-repository/ # Define address, 
    #                                                   #including a final slash.
    #      system(paste('mkdir',pathout)) # Create folder at said address.
    #      setwd(pathout)  # Go to this directory. This is useful if you want to 
    #                                         #save additional tables or figures.

  # 
  ## A-site coverage periodicity by length
  #

  suppressMessages(periodicity(list.bam, refCDS, refFASTA, pathout, XP.names, 
                               versionStrip = FALSE,
                               python.messages=FALSE))

  # 
  ## Select footprint length with sufficient periodicity
  #

  attach(listminmax <- select.FPlen(list.bam, pathout, XP.names))

  #
  ## Codon occupancy, codon enrichment.
  # 

  enrichmentNoccupancy(list.bam, refCDS, refFASTA, mini, maxi, XP.names,  
                         pathout, versionStrip = FALSE,
                         r.messages=FALSE,
                         python.messages=FALSE)
   
  #
  ## Replicates.
  #

  repl.correl.counts.Venn.res <- repl.correl.counts.Venn(XP.conditions, XP.names, 
                                                         pathout, r.messages=FALSE)
 
  repl.correl.gene.res <- repl.correl.gene(XP.conditions, XP.names, pathout)
 
  repl.correl.codon.res <- repl.correl.codon(list.bam, refCDS, refFASTA, 
                                             mini, maxi, 
                                             XP.names, XP.conditions, pathout)

  repl.correl.heatmap.res <- repl.correl.heatmap(XP.conditions.i, XP.names, pathout)

}
}

