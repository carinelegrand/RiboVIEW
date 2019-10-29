#!/bin/R


# enricht.aroundA : Plot enrichment of all codons in large [-90 ; +90] or closed in ()
#                     windows around A (= position 0)
#
#
# This function takes a list of samples and experimental conditions in input and
#   plots nerichment around A site.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#
# Out:
#       This function returns a list containing the following : 
#         - plot        Address of png file containing track plot
#
enricht.aroundA <- function(XP.conditions, XP.names, pathout) {

  nsamples <- length(XP.names)
  couleurs3 <- "dodgerblue3"
  # Codon reordering for this plot
  codgen.order <- c("uuu","uuc","uua","uug", "cuu","cuc","cua","cug",
                    "auu","auc","aua","aug", "guu","guc","gua","gug",
                    "ucu","ucc","uca","ucg", "ccu","ccc","cca","ccg",
                    "acu","acc","aca","acg", "gcu","gcc","gca","gcg",
                    "uau","uac","uaa","uag", "cau","cac","caa","cag",
                    "aau","aac","aaa","aag", "gau","gac","gaa","gag",
                    "ugu","ugc","uga","ugg", "cgu","cgc","cga","cgg",
                    "agu","agc","aga","agg", "ggu","ggc","gga","ggg")

  ### loop on samples
  for (isample in 1:nsamples) {

    ### Read enrichment
    ce.m <- read.table(paste("Sample-",XP.names[isample],"_Codon-enrichment-unbiased_mean",sep=""), 
                        header=TRUE, row.names = 1)
    
    ### Filter out STOP codons and reorder
    ce.m.2 <- ce.m[ !(row.names(ce.m) %in% c("uaa","uag","uga")), ]
    ce.m.3 <- ce.m.2[codgen.order,]
    
      
    ### Plot all, separately, in+/-90 window around A
    ymin <- min( 0.5 , min(ce.m.3, na.rm=TRUE) )
    ymax <- max( 1.8 , max(ce.m.3, na.rm=TRUE) )
    png.enrichment.all.png <- paste(pathout, "enrichment-all_",XP.names[isample],".png", sep="")
    png.enrichment.all.eps <- paste(pathout, "enrichment-all_",XP.names[isample],".eps", sep="")
    for (plotformat in c("png","eps")) {
        if (plotformat=="png") {
          grDevices::png(file = png.enrichment.all.png,  width = 400, height = 1600)
        } else {
          grDevices::cairo_ps(png.enrichment.all.eps,width = 5,height = 20)
        }
        # Divide in 64 boxes, +1 line and +1 column for axes
        par(mfrow=c(17,5), las=1, mar=0.2*c(1,1,1,1))
        #par(mfrow=c(16,4), las=1, graphics::par(mar=0.1*c(1,1,1,1)))
        for (i in 1:64) { 
          # Plot axis if first box in the row of plots
          if (i%%4 == 1) {
            plot(-90:90, rep(NA,181), col=NULL, 
               ylim=c(ymin,ymax), xlim=c(-90,90), axes=FALSE, main='', xlab='', 
               ylab='', frame=FALSE)
            axis(side=2, line=-10)
            par(las=3) ; mtext(text='Enrichment (-)', side=2, line=-7.5, cex=0.5) ; par(las=1)
          }        
          # Plot of enrichment
          if (codgen.order[i] != "uaa" & codgen.order[i] != "uag" & codgen.order[i] != "uga") {
          plot(-90:90, ce.m.3[i,], #/median(as.numeric(ce.m.3[i,])), 
               col=couleurs3, type='l', ylim=c(ymin,ymax),
               xlim=c(-90,90), axes=FALSE, main='', xlab='', ylab='', frame=TRUE)
          legend("topleft",legend=codgen.order[i], bty="n", adj=c(0.8,0))
          abline(h=1, lwd=0.5)
          abline(v=0, lty=2, col='grey25', lwd=0.5)
          lines(-90:90, ce.m.3[i,], col=couleurs3)  
          } else {plot(0, 0, col="white", type='l', ylim=c(0,3),
               xlim=c(-15,15), axes=FALSE, main='', xlab='', ylab='', frame=TRUE)
               text(x=0, y=1.5, labels="STOP codon", cex=0.7)
          }
        }
        # One blank box
        plot(-90:90,rep(NA,181), col=NULL, ylim=c(ymin,ymax),
             xlim=c(-90,90), axes=FALSE, main='', xlab='', ylab='', frame=FALSE)
        
        # Plot axes below each column of plot boxes
        for (ii in 1:4) {
          plot(-90:90,rep(NA,181), col=NULL, ylim=c(ymin,ymax),
               xlim=c(-90,90), axes=FALSE, main='', 
               xlab='', 
               ylab='', frame=FALSE)
          axis(side=1, line=-4.5, at=c(-90,-45,0,45,90))
          mtext(text="Codon position relative to A-site", 
                  side=1, line=-2.5, cex=0.5)
          mtext(text="(3' to 5')", 
                  side=1, line=-2, cex=0.5)
        }
      dev.off()
    }
    
    ### Plot all, separately, closed in on ribosome [-15 ; +15]
    ymin <- min( 0.5 , min(ce.m.3[91+c(-15:15)], na.rm=TRUE) )
    ymax <- max( 1.8 , max(ce.m.3[91+c(-15:15)], na.rm=TRUE) )
    png.enrichment.all.zoom.png <- paste(pathout, "enrichment-all-zoom_",XP.names[isample],".png", sep="")
    png.enrichment.all.zoom.eps <- paste(pathout, "enrichment-all-zoom_",XP.names[isample],".eps", sep="")
    for (plotformat in c("png","eps")) {
        if (plotformat=="png") {
          grDevices::png(file = png.enrichment.all.zoom.png,  width = 400, height = 1600)
        } else {
          grDevices::cairo_ps(png.enrichment.all.zoom.eps,width = 5,height = 20)
        }
        par(mfrow=c(16,4), las=1, graphics::par(mar=0.1*c(1,1,1,1)))
        for (i in 1:64) { 
          if (i != 35 & i != 36 & i != 51) {
          plot(-90:90, ce.m.3[i,], 
               col=couleurs3, type='l', ylim=c(ymin,ymax),
               xlim=c(-15,15), axes=FALSE, main='', xlab='', ylab='', frame=TRUE)
          legend("topleft",legend=codgen.order[i], bty="n", adj=c(0.8,0))
          abline(h=1, lwd=0.5)
          abline(v=0, lty=2, col='grey25', lwd=0.8)
          abline(v=c(-4,5), col="lightgrey", lty=2, lwd=0.8)
          lines(-90:90, ce.m.3[i,], col=couleurs3)  
          } else {plot(0, 0, col="white", type='l', ylim=c(0,3),
               xlim=c(-15,15), axes=FALSE, main='', xlab='', ylab='', frame=TRUE)
               text(x=0, y=1.5, labels="STOP codon", cex=0.7)
          }
        }
      dev.off()
    }
    
  } # End of loop on samples

  # Save values for later reference, and return
  enricht.aroundA.res <- list(plot=png.enrichment.all.png)
  save(enricht.aroundA.res, 
       file=paste(pathout, "valNrec_", "enricht.aroundA.res",".RData", sep=""))
  return(enricht.aroundA.res)

}
