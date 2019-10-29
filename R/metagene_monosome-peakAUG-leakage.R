# Monosome selection, START and STOP leakage, metagene
#
#
#
# This function generates plots and indicative values for monosome selection, 
#   inflation around AUG codon and around STOP codon.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#        res1           Resolution of metagene
#        res2           Resolution of metagene around AUG and STOP
#
# Out:
#   List of lists :
#     Monosome selection (UTR reads) :
#          - plot                Address of plot file in png format
#          - value               Fraction of mRNAs where footprint cover only CDS and no UTRs
#                                  (median and inter-quartile range)
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#     Inflation of reads around AUG :
#          - plot                Address of plot file in png format
#          - value               Median ratio of coverage in 50 first codons in CDS / all CDS
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#     Leakage at STOP codon :
#          - plot                Address of plot file in png format
#          - value               list : - value$start  if any positive slope of metagene after AUG, p-value for this slope
#                                       - value$stop   Ratio of coverage in [1 ; 1.3] just after STOP, as compared to [0.9 ; 1], just before STOP codon
#                                                        (median and standard deviation)
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#
metagene.all <- function(XP.conditions, XP.names, pathout, res1=0.1, res2=0.01) {
  # Inputs
  # - XP.conditions list of positive integers   Experimental condition by sample
  # - XP.names      list of strings             Sample names
  # - res1          real positive number        resolution of metagene
  # - res2          real positive number        resolution of metagene around AUG and STOP

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))

  # Parameters
  colSamples.pre1 <- RColorBrewer::brewer.pal(9, "Set1")
  colSamples.pre2 <- colSamples.pre1[-6]

  # Number of arguments, reconstruct listSamples
  n <- length(XP.conditions)
  colSamples = rep(colSamples.pre2, ceiling(n/9))

  # Initialisations
  covUTR                <- c()
  metaG.m.all           <- c()
  metaG.iqr.all         <- c()
  metagene.start.all    <- c()
  leakage.start.all     <- c()
  leakage.stop.all      <- c()
  infl.all              <- c()
  fx.0covUTR.all        <- c()
  fxCovUTR_med          <- c()
  fxCovUTR_IQR          <- c()
  pcCovUTR              <- c()

  # Per sample, Nb footprints in UTRs and CDS, as well as in small intervals along a standard transcripts
  for (i in 1:n) {

    inFile     <- paste(pathout, "Sample-",XP.names[i],".metagene", sep="")
    metagene.i <- metageneExpl(inFile, res1, res2)

    # gather values over several samples
    covUTRi                   <- metagene.i$covUTR
    fx.0covUTR.all[[i]]       <- covUTRi[1]
    fxCovUTR_med[[i]]         <- covUTRi[2]
    fxCovUTR_IQR[[i]]         <- covUTRi[3]
    pcCovUTR[[i]]             <- covUTRi[4]
    infl.all[[i]]             <- metagene.i$infl
    metaG.m.all[[i]]          <- metagene.i$metaG_m
    metaG.iqr.all[[i]]        <- metagene.i$metaG_iqr
    metagene.start.all[[i]]   <- metagene.i$metaG_start
    leakage.start.all[[i]]    <- metagene.i$leakStart
    leakage.stop.all[[i]]     <- metagene.i$leakStop
  }

  #
  ## monosome selection : amount of UTR reads
  #

  # plot - overview metagene
  metagene.monosome.plot.png <- paste(pathout,"Metagene-monosome.png",sep="") #tempfile()
  metagene.monosome.plot.eps <- paste(pathout,"Metagene-monosome.eps",sep="") #tempfile()
  #
  for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(metagene.monosome.plot.png, width = 800 ,height = 800)
      } else {
        grDevices::cairo_ps(metagene.monosome.plot.eps,width = 10,height = 10)
      }
    #grDevices::png(metagene.monosome.plot.png, width = 800 ,height = 800)
      breaks <- seq(-1,2,res1)
      mids   <- seq((-1+res1/2),(2-res1/2),res1)
      #
      graphics::par(mar=c(3,4,2,2),mfrow=c(1,1))
      ymax <- 1.2 * (max(sapply(metaG.m.all, max)) + 0.0 * max(sapply(metaG.iqr.all, max)))
      graphics::plot(NA, xlab="", ylab="Probability", main=' ',
           type='l', col='dodgerblue4', lwd=2, axes=FALSE,
           xlim=c(-1,2),
           ylim=c(0, ymax))   #lwd=2
      graphics::lines(x = c(-1,2), y = rep(0,2), col='grey', lty=3)
      graphics::lines(x = c(0,0), y = c(0,ymax), col='grey', lty=3)
      graphics::lines(x = c(1,1), y = c(0,ymax), col='grey', lty=3)
      for (i in 1:n) {
        graphics::lines(mids, metaG.m.all[[i]], col=colSamples[i], lwd=2, lty=1)
      }
      graphics::par(las=1)
      graphics::axis(2)
      offset=1
      graphics::axis(1, at = c(-1,0,1,2), labels=rep("",4),tcl=-0.75, line=-0.8+offset)
      graphics::axis(1, at=c(-0.5,0.5,1.5), labels = c("5'UTR", "CDS", "3'UTR"), tick = FALSE, line=-1.5+offset)
      graphics::axis(1, at = c(-1,0,1,2), labels=rep("",4), tcl=0.75, line=0.7+offset)
      #
      graphics::par(xpd=TRUE) #allow legend outside plot window
      #
      # index of start and stop of CDS
      mi.i <- 1 + round(1/res1)
      ma.i <- round(2/res1)
      leg.text <- c()
      for (i in 1:n) {
        metaG.iqr.i <- round(stats::median(metaG.iqr.all[[i]][mi.i:ma.i]),3)
        leg.text    <- c(leg.text, paste(XP.names[i], ", IQR=",metaG.iqr.i,sep=""))
      }
      graphics::legend("topright",col=colSamples[seq(1,n)],
             lty=1, lwd=2,
             cex = 0.8, #bty='n',
             horiz = FALSE,
             x.intersp = 0.5,
             legend=leg.text)
  grDevices::dev.off()
  } # end loop EPS or PNG

  # Fraction of mRNAs where footprint cover only the CDS and no locus at all in UTRs
  #metagene.monosome.Value.pre <- median(fx.0covUTR.all)
  metagene.monosome.Value <- c(stats::median(pcCovUTR), stats::IQR(pcCovUTR))   #100*(1 - metagene.monosome.Value.pre)

  # Color
  metagene.monosome.Color <- ifelse(metagene.monosome.Value[1] > 10,
                                    yes = "red",
                                    no = ifelse(metagene.monosome.Value[1] > 1,
                                                yes="orange",
                                                no="green"))

  # Recommendation ; use cat(Recommendation)
  metagene.monosome.Recommendation <- ifelse(metagene.monosome.Value[1] > 10,
                           yes = paste("Monosome selection : More than 10% of mRNAs have footprints not only in the CDS,\n",
                                       "   but also in the UTRs. This points to impaired monosome \n",
                                       "   selection due to incomplete digestion.\n",
                                       "Drop-off in case unwanted aminoacid starvation : see Johnson and Li 2018.\n",
                                       "Further info on monitoring footprints using metagene : Sin et al 2016.\n",
                                     sep=""),
                           no = ifelse(metagene.monosome.Value[1] > 1,
                                       yes=paste("More than 1% of mRNAs have footprint not only in the CDS,\n",
                                       "   but also in the UTRs. This may point to impaired monosome \n",
                                       "   selection due to incomplete digestion.\n",
                                     sep=""),
                                       no=paste("About ", round(metagene.monosome.Value) ,"% of mRNAs have footprints not only\n",
                                       "   in the CDS, but also in the UTRs. Though this is not high enough to suggest\n",
                                       "   impaired monosome selection, you might want to compare this value \n",
                                       "   between replicates or conditions.\n",
                                     sep="")))

  #
  ##
  ### Inflation of CDS-start codons coverage
  ##
  #

  # plot metagene zoom at AUG
  metagene.inflation.plot.png <- paste(pathout,"Metagene-inflation.png",sep="") #tempfile()
  metagene.inflation.plot.eps <- paste(pathout,"Metagene-inflation.eps",sep="") #tempfile()
  #
  for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(metagene.inflation.plot.png, width = 800 ,height = 800)
      } else {
        grDevices::cairo_ps(metagene.inflation.plot.eps,width = 10,height = 10)
      }
    #grDevices::png(metagene.inflation.plot.png, width = 800 ,height = 800)
      graphics::par(mar=c(4,4,2,2))
      #
      maxc <- max(sapply(metagene.start.all,max))
      dmax <- log(maxc, base = 10)
      # manip to have useful but also round upper bound
      ymax <- 10^(floor(dmax)-1) * ceiling(1.1 * maxc / (10^(floor(dmax)-1)) )
      # void annotated plot
      graphics::plot(NA,xlab="codon position / START", ylab="Frequency",
           xlim=c(-20,50), ylim=c(0,ymax),
           axes=FALSE, main="")
      graphics::lines(x = c(0,0), y=c(0,ymax), col="grey10", lty=2)
      graphics::axis(1)
      graphics::par(las=1)
      graphics::axis(2, at = seq(0, ymax, round(ymax/4)))
      for (i in 1:n) { graphics::lines(x = seq(-20,50,1), metagene.start.all[[i]], col=colSamples[i]) }
      #
      graphics::legend("topleft",col=colSamples[seq(1,n)],
             lty=1, lwd=2,
             cex = 0.8, #bty='n',
             horiz = FALSE,
             x.intersp = 0.5,
             legend=XP.names)
  grDevices::dev.off()
  } # end loop EPS or PNG

  # Value
  # Ratio CDS 50 first codons in CDS / overall
  metagene.inflation.Value <- list(median=stats::median(infl.all), sd=stats::sd(infl.all))

  # Color
  metagene.inflation.Color <- ifelse(metagene.inflation.Value$median > 0.01, yes = "red", no = ifelse(metagene.inflation.Value$median > 0.005, yes="orange", no="green"))

  # Recommendation ; use cat(Recommendation)
  # compare CHX / non-CHX datasets
  metagene.inflation.Recommendation <- ifelse(metagene.inflation.Value$median > 0.01,
                           yes = paste("More than 1% of footprint coverage is located near CDS Start.\n",
                                       "  This strongly suggests continued initiation happened during the experiment.\n",
                                       "  Results will most likely be biased, unless you exclude a large amount of codons near AUG.\n",
                                     sep=""),
                           no = ifelse(metagene.inflation.Value$median > 0.005,
                                       yes=paste("More than 0.5% of footprint coverage is located near CDS Start.\n",
                                       "  This possibly indicates continued initiation during the experiment.\n",
                                       "  Results might be biased, check the number of codons to exlude int he vicinity of AUG.\n",
                                     sep=""),
                                       no=paste("About ", round(100*metagene.inflation.Value$median,2) ,
                                                          "% (sd=",round(100*metagene.inflation.Value$sd,4),
                                                          "%) of footprint coverage is located near CDS start.\n",
                                                "  This is in line with most RIBO-seq experiments.\n",
                                     sep="")))



  #
  ##
  ### Leakage at start and stop codons
  ##
  #

  # plot metagene zoom around AUG and Stop
  # AUG+-20 profile and TAA,TGA.TAG +/-20 profile
  metagene.leakage.plot.png <- paste(pathout,"Metagene-leakage.png",sep="") #tempfile()
  metagene.leakage.plot.eps <- paste(pathout,"Metagene-leakage.eps",sep="") #tempfile()
  #
  for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(metagene.leakage.plot.png, width = 800 ,height = 800)
      } else {
        grDevices::cairo_ps(metagene.leakage.plot.eps,width = 10,height = 10)
      }
    #grDevices::png(metagene.leakage.plot.png, width = 800 ,height = 800)
      graphics::par(mar=c(4,4,2,2), mfrow=c(2,1))
      #
      # top panel, START
      maxc <- max(sapply(leakage.start.all,max))
      dmax <- log(maxc, base = 10)
      # manip to have useful but also round upper bound
      ymax <- 10^(floor(dmax)-1) * ceiling(1.1 * maxc / (10^(floor(dmax)-1)) )
      # void annotated plot
      graphics::plot(NA,xlab="standardized position / START", ylab="Frequency",
           xlim=c(-0.1,0.3+res2), ylim=c(0,ymax),
           axes=FALSE, main="")
      graphics::lines(x = c(0,0), y=c(0,ymax), col="grey10", lty=2)
      graphics::axis(1)
      graphics::par(las=1)
      graphics::axis(2, at = seq(0, ymax, round(ymax/4)), cex.axis = 0.8)
      for (i in 1:n) {
        graphics::lines(x = seq(-0.1+res2/2 , 0.3-res2/2 , res2),
              y = leakage.start.all[[i]],
              col = colSamples[i])
      }
      #
      # bottom panel, STOP
      maxc <- max(sapply(leakage.stop.all,max))
      dmax <- log(maxc, base = 10)
      # manip to have useful but also round upper bound
      ymax <- 10^(floor(dmax)-1) * ceiling(1.1 * maxc / (10^(floor(dmax)-1)) )
      # void annotated plot
      graphics::plot(NA,xlab="standardized position / STOP", ylab="Frequency",
           xlim=1+c(-0.1,0.3+res2), ylim=c(0,ymax),
           axes=FALSE, main="")
      graphics::lines(x = c(0,0), y=c(0,ymax), col="grey10", lty=2)
      graphics::axis(1)
      graphics::par(las=1)
      graphics::axis(2, at = seq(0, ymax, round(ymax/4)), cex.axis = 0.8)
      for (i in 1:n) { graphics::lines(x = 1+seq(-0.1+res2/2 , 0.3-res2/2 , res2), leakage.stop.all[[i]], col=colSamples[i]) }
      #
      graphics::legend("topright",col=colSamples[seq(1,n )],
             lty=1, lwd=2,
             cex = 0.8, #bty='n',
             horiz = FALSE,
             x.intersp = 0.5,
             legend=XP.names)
  grDevices::dev.off()
  } # end loop EPS or PNG

  # Value
  # Fisher test and robust-lm-slope with pvalue
  leak.start <- 1.
  leak.stop <- c()
	absc.meta <- seq(-0.1+res2/2 , 0.3-res2/2 , res2)
	# Iterate over samples
  for (i in 1:n) {
    # Slope after AUG from robust linear fit ; retain if positive with significant p-value
    fit.Start <- MASS::rlm(leakage.start.all[[i]] ~ absc.meta , subset = which(absc.meta >0))
    # Retry with more iterations if non convergence
    if (!fit.Start$converged) {
      print("Retrying rlm with >20 steps (in function metagene.all).")
      fit.Start <- MASS::rlm(leakage.start.all[[i]] ~ absc.meta , subset = which(absc.meta >0),
                             maxit = 40)
      }
    fit.Start.expl <- exploitfitrob(fit.Start)
    # leak.start : if any positive slope after AUG, p-value corresponding to this slope.
    if (fit.Start.expl$Coeff > 0 & fit.Start.expl$p < 0.05) {
      leak.start <- c(leak.start ,
                      fit.Start.expl$p)
    }
    # leak.stop : Ratio of coverage in [1 ; 1.3] just after STOP, as compared to [0.9 ; 1], just before STOP codon
    leak.stop <- c(leak.stop,
                   stats::median(leakage.stop.all[[i]][ (1+absc.meta) > 1 ]) / stats::median(leakage.stop.all[[i]][ (1+absc.meta) <= 1 ]))
  }
  #
  metagene.leakage.Value <- list(start = list(p = min(leak.start)),
                                 stop  = list(median = stats::median(leak.stop), sd = stats::sd(leak.stop)/sqrt(n)) )

  # Color and recommendation ; use cat(Recommendation)
  # Note : message is consituted of one recommendation about leakage after the start codon,
  #                                 and a second for leakage after the stop codon.
  metagene.leakage.Recommendation <- ""
  if ( metagene.leakage.Value$start$p < 0.05 | metagene.leakage.Value$stop$median > 5/100 ) {
      metagene.leakage.Color          <- "red"
      if ( metagene.leakage.Value$start$p < 0.05 ) {
          metagene.leakage.Recommendation <- paste(metagene.leakage.Recommendation,
                                             "Footprint frequency increases significantly (p=",metagene.leakage.Value$start$p,") after AUG.\n",
                                             "  This strongly suggests continued elongation after initiation has been blocked, possibly ribosome run-off for short transcripts.\n",
                                             "  In some cases this might be an organism- or sample-dependent characteristic. Unless this is the case,\n",
                                             "  you might remove the affected codons before any other calculation,\n",
                                             "  or to adapt the experimental protocol.\n",
                                           sep="") }
      if ( metagene.leakage.Value$stop$median > 5/100 ) {
          metagene.leakage.Recommendation <- paste(metagene.leakage.Recommendation,
                                             "Footprint frequency after STOP is higher than 5% of its value before STOP (median=",as.character(round(100*metagene.leakage.Value$stop$median,2)),", sd=",
                                                                                              as.character(round(100*metagene.leakage.Value$stop$sd,4)),
                                                                                              ").\n",
                                             "  This strongly suggests leakage of the stop codon. This might be linked to incomplete annotation a\n",
                                             "  or to incomplete ribosome arrest at UAA/UAG or UGA codons.\n",
                                             sep="") }
  } else if ( metagene.leakage.Value$start$p < 0.1 | metagene.leakage.Value$stop$median > 1/100) {
      metagene.leakage.Color          <- "orange"
      if ( metagene.leakage.Value$start$p < 0.1 ) {
          metagene.leakage.Recommendation <- paste(metagene.leakage.Recommendation,
                                             "Footprint frequency increases (p=",metagene.leakage.Value$start$p,") after AUG.\n",
                                             "  This possibly suggests continued elongation after initiation has been blocked.\n",
                                             "  In some cases this might be an organism- or sample-dependent characteristic. Unless this is the case,\n",
                                             "  we recommend removing as many codons as necessary after AUG before further calculation,\n",
                                             "  or to adapt the experimental protocol.\n",
                                           sep="") }
      if ( metagene.leakage.Value$stop$median > 1/100 ) {
          metagene.leakage.Recommendation <- paste(metagene.leakage.Recommendation,
                                             "Footprint frequency after STOP is higher than 1% of its value before STOP  (median=",
                                                                                              as.character(round(100*metagene.leakage.Value$stop$median,2)),
                                                                                              ", sd=",
                                                                                              as.character(round(100*metagene.leakage.Value$stop$sd,4)),
                                                                                              ").\n",
                                             "  This possibly suggests leakage of the stop codon. This might be linked to incomplete annotation a\n",
                                             "  or to incomplete ribosome arrest at UAA/UAG or UGA codons.\n",
                                             sep="") }
  } else {
      metagene.leakage.Color          <- "green"
      metagene.leakage.Recommendation <- paste(
                                         "Footprint frequency does not increase (p=",metagene.leakage.Value$start$p,") after AUG.\n",
                                         "  This is in agreement with most RIBO-Seq experiments, where frequency steeply decreases after a peak at AUG.\n",
                                         "Footprint frequency  after STOP is less than 1% of its value before STOP (median=",
                                                                                          as.character(round(100*metagene.leakage.Value$stop$median,2)),
                                                                                          ", sd=",
                                                                                          as.character(round(100*metagene.leakage.Value$stop$sd,4)),
                                                                                          ").\n",
                                         "  This corresponds to no or minimal leakage of the stop codon, due to effective ribosome arrest at UAA/UAG or UGA codons.\n",
                                         sep="")
  }

  # Result
  metagene.res <- list(list(plot = metagene.monosome.plot.png,    value = metagene.monosome.Value,
                            color = metagene.monosome.Color,  recommendation = metagene.monosome.Recommendation),
                       list(plot = metagene.inflation.plot.png,   value = metagene.inflation.Value,
                            color = metagene.inflation.Color, recommendation = metagene.inflation.Recommendation),
                       list(plot = metagene.leakage.plot.png,     value = metagene.leakage.Value,
                            color = metagene.leakage.Color,   recommendation = metagene.leakage.Recommendation)
                     )
  # Save values for later reference, and return
  save(metagene.res, 
       file=paste(pathout, "valNrec_", "metagene.res",".RData", sep=""))
  return(metagene.res)

}



