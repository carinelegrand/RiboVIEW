#enrichment and and BCO
#
# Contents :
#           enrichmentNoccupancy
#           occupancyNwaveParallel_one
#           

# enrichmentNoccupancy : Calculates unbiased codon enrichment, Hussmann's codon
#                          enrichment and bulk codon occupancy.
#
#
#
# This function takes a list of samples and their conditions as input and 
#   visualizes enrichment around AUG +/- 90 codons, where possible artefacts due
#   to drugs used in the experiment should be visible.
#
# In:
#        list.bam                List of bam files containing aligned reads for 
#                                  each sample (same order as in XP.names)
#        refCDS                  Address of file containing coding sequence
#                                  annotation. This file should contain tab-separated
#                                  values for mRNA name, nt position of start 
#                                  codon, nt position of first codon after stop 
#                                  codon, for example : 
#                                       ID	localStart	localEnd
#                                       NM_001276351	145	1348
#                                       NM_001276352	145	799
#                                       NM_000299	251	2495
#        refFASTA                Address of reference sequences for mRNA of the 
#                                  studied organism
#        mini                    Minimum footprint length to consider (as selected by user)
#        maxi                    Maximum footprint length to consider (as selected by user)
#        pathout                 Address where output files will be written.
#        versionStrip            Indicates if version number should be trimmed, 
#                                  for example NM_001276351.1 should be trimmed 
#                                  to NM_001276351 (defaults to FALSE)
#
# Out:
#       Output files ar written directly to pathout. These files are :
#          Sample-*_Codon-enrichment-unbiased_mean  Mean codon enrichment (unbiased)
#          Sample-*_Codon-enrichment-unbiased_sd    Standard deviation of codon enrichment
#          Sample-*_Hussmann_mean       Mean codon enrichment (Hussmann)
#          Sample-*_Hussmann_sd         Standard deviation of Hussmann's codon enrichment
#          Sample-*_BCO                   Bulk Codon Occupancy
#          Sample-*_limCod                Codon counts at footprint limits
#          Sample-*_limNt                 Nucleotide counts at footprint limits
#          Sample-*_nbreads               Nb of reads per mRNA
#          Sample-*_SCO                   Coverage for A site per mRNA and position, and single codon occupancy
#          Sample-*_metagene              Coverage for A site per mRNA and metagene position
#                                            Note that this contains also UTRs and 
#                                            regions close to AUG or STOP codons,  
#                                            contrary to other output files.
#          Sample-*_paused                List of paused sites with their sequence context
#          Sample-*_RPKM                  Number of reads per kilo base per 
#                                            million mapped reads, for each mRNA.
#

enrichmentNoccupancy <- function(list.bam, refCDS, refFASTA, mini, maxi, XP.names, 
                      pathout, versionStrip = FALSE, r.messages=TRUE, 
                      python.messages=TRUE) {
  
  # Determine number of cores for parallel computations
  nbCores.rec <- max( (parallel::detectCores() - 1) , 1)
  nbCores.pardefaut <- min(2,nbCores.rec)
  if (r.messages) {
    nbCores <- as.integer(readline(
                       prompt = paste( "Please enter the number of cores you would like to use for parallel computations :\n",
                        " (recommended value for your computer : ",nbCores.rec,"), \n",
                        " then press [enter] to continue.\n",sep="")))
    # Use nbCores.pardefaut as default value
    if (is.na(nbCores)) { cat("Number of cores was not set. Using default value (",nbCores.pardefaut,") :\n\n")
                          nbCores <- nbCores.pardefaut }
   } else {
     nbCores <- nbCores.pardefaut
   }
                  
  if (r.messages) {
    cat(paste("Starting codon occupancy and codon enrichment calculations. \n",
            "  This should take no more than ",ceiling(length(XP.names)/nbCores),"h",
            "  (~=1h per ~60Mo bam file on 1 core).\n",sep=""))
    cat("NOTE : the tables of values of enrichment, etc., will be available in ",
      pathout,
      " under the following names, for sample 1 :\n", 
      "   - Sample-", XP.names[1],"_BCO       Codon occupancy for A-site, P-site, \n",
      "                                         E-site and 3 codons upstream and downstream,\n",
      "   - Sample-", XP.names[1],"_Codon-enrichment-unbiased_mean     \n",
      "                                       Unbiased codon enrichment for A-site,\n",
      "                                         P-site, E-site and 90 codons upstream \n",
      "                                         and downstream, average over mRNAs  \n",
      "   - Sample-", XP.names[1],"_Codon-enrichment-unbiased_sd      \n",
      "                                       Unbiased codon enrichment for A-site,\n",
      "                                         P-site, E-site and 90 codons upstream \n",
      "                                         and downstream, standard deviation over mRNAs  \n",
      "   - Sample-", XP.names[1],"_Hussmann_mean  \n",
      "                                       Codon enrichment following Hussmann\n",
      "                                         et al. 2015 (doi : 10.1371/journal.pgen.1005732), average\n",
      "   - Sample-", XP.names[1],"_Hussmann_sd  \n",
      "                                       Codon enrichment following Hussmann\n",
      "                                         et al. 2015 (doi : 10.1371/journal.pgen.1005732), standard deviation\n",
      "   - Sample-", XP.names[1],".limCod    Codon count at footprint start and footprint end, \n",
      "   - Sample-", XP.names[1],".limNt     Nucleotide count at footprint start and footprint end, \n",
      "   - Sample-", XP.names[1],".metagene  Coverage of codons in A-site at standardized positions along CDS, to construct metagene plot, \n",
      "   - Sample-", XP.names[1],".nbreads   Count of footprints per mRNA,\n",
      "   - Sample-", XP.names[1],".paused    Repertoire of paused sites (in A site, codons 3 times more paused than codons in the vicinity), \n",
      "   - Sample-", XP.names[1],".RPKM      Count of reads per kilobase exon per million mapped,\n",
      "   - Sample-", XP.names[1],".SCO       Single-Codon Occupancy = single-codon coverage in A-site / average coverage in an mRNA. \n\n",
      " Calculations now ongoing... \n\n",
     sep="")
  }          
  # Initialize cluster of parallel cores
  cl <- parallel::makeCluster(nbCores)
  
  # Make libraries available inside cluster
  parallel::clusterEvalQ(cl, { library(rPython) })

  ## Determine number of samples
  n <- length(list.bam)
  #
  parallel::parLapply(cl, 1:n, occupancyNwaveParallel_one , list.bam, refCDS, 
      refFASTA, mini, maxi, XP.names, pathout, versionStrip, python.messages)
  #
  parallel::stopCluster(cl)
  
}

occupancyNwaveParallel_one <- function(i, list.bam, refCDS, refFASTA, mini, maxi, 
  XP.names, pathout, versionStrip = FALSE, python.messages = TRUE) {

  readsBAM <- list.bam[[i]]
  outBCO <- paste(pathout, "Sample-", XP.names[i], sep="")
  
  pyEnrichment_filename <- system.file("enrichment.py", package="RiboVIEW", mustWork = TRUE)
  PythonInR::pyExecfile(pyEnrichment_filename)
  pyEnrichment <- PythonInR::pyFunction("enrichment")
  pyEnrichment(readsBAM, refCDS, refFASTA, mini, maxi, outBCO, 
                versionStrip=FALSE,
                messages=python.messages)
  #rPython::python.load( system.file("enrichment.py", package="RiboVIEW", mustWork = TRUE))
  #rPython::python.call("enrichment", readsBAM, refCDS, refFASTA, mini, maxi, outBCO, versionStrip, ex3, ex5)

}


