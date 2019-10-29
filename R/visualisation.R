# Visualisation of enrichment and tracks
#
# Contents
#          col2trspcol
#          read.pytable
#          generate.m.s
#          visu.m.s.enrichmnt
#          visu.tracks

col2trspcol <- function(couleur, alpha=0.8) {grDevices::rgb(t(grDevices::col2rgb(couleur)/255), alpha=alpha)}

# An example of color-coding by codon
# Nice resources to select colors from : 
#          - Color palettes in R package "RColorBrewer", 
#          - A wealth of color shades in https://www.toutes-les-couleurs.com/en/html-color-code.php
#                 => Use html color codes, for instance "#007FFF" for "azur", 
#                 since html color names are not all built-in in R.
codon.colors=rbind( c('aaa','darkgrey'),  c('aac','orange'),    c('aag','darkgrey'),  c('aau','orange'),
                    c('aca','darkgrey'),  c('acc','darkgrey'),  c('acg','darkgrey'),  c('acu','darkgrey'),
                    c('aga','darkgrey'),  c('agc','darkgrey'),  c('agg','darkgrey'),  c('agu','darkgrey'),
                    c('aua','darkgrey'),  c('auc','darkgrey'),  c('aug','darkgrey'),  c('auu','darkgrey'),
                    c('caa','darkgrey'),  c('cac','orange'),    c('cag','darkgrey'),  c('cau','orange'),
                    c('cca','darkgrey'),  c('ccc','darkgrey'),  c('ccg','darkgrey'),  c('ccu','darkgrey'),
                    c('cga','darkgrey'),  c('cgc','darkgrey'),  c('cgg','darkgrey'),  c('cgu','darkgrey'),
                    c('cua','darkgrey'),  c('cuc','darkgrey'),  c('cug','darkgrey'),  c('cuu','darkgrey'),
                    c('gaa','darkgrey'),  c('gac','orange'),    c('gag','darkgrey'),  c('gau','orange'),
                    c('gca','darkgrey'),  c('gcc','darkgrey'),  c('gcg','darkgrey'),  c('gcu','darkgrey'),
                    c('gga','darkgrey'),  c('ggc','darkgrey'),  c('ggg','darkgrey'),  c('ggu','darkgrey'),
                    c('gua','darkgrey'),  c('guc','darkgrey'),  c('gug','darkgrey'),  c('guu','darkgrey'),
                    c('uaa','red'),       c('uac','orange'),    c('uag','red'),       c('uau','orange'),
                    c('uca','darkgrey'),  c('ucc','darkgrey'),  c('ucg','darkgrey'),  c('ucu','darkgrey'),
                    c('uga','red'),       c('ugc','darkgrey'),  c('ugg','darkgrey'),  c('ugu','darkgrey'),
                    c('uua','darkgrey'),  c('uuc','darkgrey'),  c('uug','darkgrey'),  c('uuu','darkgrey') )

#### Read from py-generated tab-separated file

read.pytable <- function(filename) {
  e.pre           <- utils::read.table(filename, sep="\t", header=TRUE, row.names = 1, check.names=FALSE, fill=FALSE)
  colnames.e      <- colnames(e.pre)[colnames(e.pre) != ""]
  e.fin           <- e.pre[,1:length(colnames.e)]
  colnames(e.fin) <- colnames.e
  return(e.fin)
}

#### Mean and SD of enrichment

# generate.m.s : Estimates codon enrichment mean, standard deviation and standard 
#                  error per experimental condition and between conditions.
#
#
#
# This function takes a list of samples and experimental conditions in input and
#   calculates mean, standard deviation and standard error for codon enrichment  
#   of all replicates of the same condition, and for comparison between conditions.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#        B              Number of bootstrap resampling occurences
#
# Out:
#       This function writes to the following files, under the address pathout :
#          - Enrichment-per-condition__weighted-mean for the mean
#          - Enrichment-per-condition__weighted-sdev for the standard deviation
#          - Enrichment-per-condition__weighted-serr for the standard error
#          - Condition_*_relative-to_*_mean-and-stderr for each comparison ratio's 
#                                                       mean and standard error
#
generate.m.s <- function(XP.conditions, XP.names, pathout, B=1000) {

  # Nb conditions and Nb replicates per condition
  conditions <- unique(XP.conditions)
  nCond      <- length(conditions)

  # Mean and SD per condition
  c.m <- c()
  allc.m <- c()
  c.s <- c()
  c.se <- c()
  compt <- 0
  for (ci in conditions) {
    compt <- compt + 1
    # Indices of replicates which match condition ci
    idces <- which(XP.conditions==ci)
    nb.idces <- length(idces)
    # list  of mean (sd, se) per sample
    l.m <- c()
    l.s <- c()
    l.se <- c()
    for (i.idces in 1:nb.idces) {
      ii <- idces[i.idces]
      # mean
      fi.m <- paste(pathout, "Sample-", XP.names[ii], "_Codon-enrichment-unbiased_mean", sep="")
      e.m <- read.pytable(fi.m)
      e.m.A <- e.m[!(row.names(e.m) %in% c("uaa","uga","uag")) ,"0", drop=FALSE]
      # sd over mRNAs
      fi.s <- paste(pathout, "Sample-", XP.names[ii], "_Codon-enrichment-unbiased_sd", sep="")
      e.s <- read.pytable(fi.s)
      e.s.A <- e.s[!(row.names(e.s) %in% c("uaa","uga","uag")) ,"0", drop=FALSE]
      # Estimate of standard error from positions outside [-10;+10]
      excl.pos <- c("-10","-9","-8","-7","-6","-5","-4","-3","-2","-1","0","1","2","3","4","5","6","7","8","9","10")
      e.n.pre <- e.m[!(row.names(e.s) %in% c("uaa","uga","uag")) ,!(colnames(e.s) %in% excl.pos)]
      e.n <- apply(e.n.pre, 1, stats::sd)
      # se simplified (1/racine(n))
      #TBD
      # se from bootstrap
      #TBD
      # Store in lists
      ni <- XP.names[ii]
      l.m[[ni]] <- e.m.A
      l.s[[ni]] <- e.s.A
      l.se[[ni]] <- e.n
    }
    
    # Overall mean, s, se
    m.m <-  as.matrix(as.data.frame(l.m))
    allc.m <- cbind(allc.m, m.m)
    m.s <-  as.matrix(as.data.frame(l.s))
    m.se <-  as.matrix(as.data.frame(l.se))
    #   weighted mean, weights=1/se**2
    enrich.m <- rowSums( (m.m / m.se**2) / apply(m.se, 1, function(x) {sum(1/x**2)} ) )
    #   weighted sd 
    enrich.s   <- ( rowSums( (m.m**2 / m.se**2) / apply(m.se, 1, function(x) {sum(1/x**2)} )) - enrich.m**2 )**0.5
    #   standard error taken as sd/sqrt(#replicates)
    enrich.se    <- enrich.s / (nb.idces)**0.5
    # Mean, standard deviation and standard error over replicates,
    #   by condition ci
    c.m[[paste(ci)]]  <- enrich.m
    c.s[[paste(ci)]]  <- enrich.s
    c.se[[paste(ci)]] <- enrich.se
  }
  
  ## Write to file
  utils::write.table(x = c.m, file=paste(pathout, "Enrichmnt-per-condition_weighted-mean.tsv", sep=""),
             sep="\t", quote=FALSE, col.names=conditions)
  utils::write.table(x = c.s, file=paste(pathout, "Enrichmnt-per-condition_weighted-sdev.tsv", sep=""),
             sep="\t", quote=FALSE, col.names=conditions)
  utils::write.table(x = c.se, file=paste(pathout, "Enrichmnt-per-condition_weighted-serr.tsv", sep=""),
             sep="\t", quote=FALSE, col.names=conditions)
             
  ## Quotient between conditions, mean and std error ; write to file
  if (nCond>1) {
    for (i in 1:(nCond-1)) {
      for (j in (i+1):nCond) {
        #
        ci1  <- conditions[i]
        idces1 <- which(XP.conditions==ci1)
        #
        ci2 <- conditions[j]
        idces2 <- which(XP.conditions==ci2)
        #
        nums <- sample(idces1, size=B, replace=TRUE)
        denoms <- sample(idces2, size=B, replace=TRUE)
        quotients <- allc.m[,nums]/allc.m[,denoms]
        # mean, SE
        q.m <- apply(quotients, 1, mean)
        q.se <- apply(quotients, 1, stats::sd)
        # Write to file
        utils::write.table(x = cbind(q.m, q.se), 
                    file=paste(pathout, "Condition_",ci1,"_relative-to_",ci2,"_mean-and-stderr.tsv", sep=""),
                    sep="\t", quote=FALSE, col.names=c("mean","stderr")) #, row.names=row.names(q.m))
      }
    }
  }
             
}

# visu.m.s.enrichmnt : Plot enrichment by condition and comparisons
#
#
#
# This function takes a list of samples and experimental conditions in input and
#   generates plots for codon enrichment of all replicates of the same condition, 
#   and for each comparison between conditions.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#
# Out:
#       This function returns a list containing the following : 
#         - plot.enrich  Address of the png plot file for codon enrichment 
#                          for each condition across replicates
#         - plot.compar  Address of one of the png plot file for codon enrichment 
#                          comparison between two different experimental conditions
#
visu.m.s.enrichmnt <- function(XP.conditions, XP.names, pathout) {

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))

  ## Nb conditions
  # Nb conditions and Nb replicates per condition
  conditions <- unique(XP.conditions)
  nCond      <- length(conditions)

  fi.m <- paste(pathout, "Enrichmnt-per-condition_weighted-mean.tsv", sep="")
  fi.s <- paste(pathout, "Enrichmnt-per-condition_weighted-sdev.tsv", sep="")
  e.M    <- utils::read.table(fi.m, sep="\t", header=TRUE, check.names=FALSE, row.names=1)
  e.SD    <- utils::read.table(fi.s, sep="\t", header=TRUE, check.names=FALSE, row.names=1)
  
  ## Plot per condition and then compare between conditions
  couleurs.pre <- c('black','orange','darkslateblue','red', 'forestgreen',
                      'purple','gold','brown','magenta','midnightblue')
  couleurs <- rep(couleurs.pre,7)

  repl.m.s.enrichmnt.png <- paste(pathout, "Visu-enrichmnt-mean-sdev-by-condition.png",sep="")
  repl.m.s.enrichmnt.eps <- paste(pathout, "Visu-enrichmnt-mean-sdev-by-condition.eps",sep="")
  for (plotformat in c("png","eps")) {
    if (plotformat=="png") {
      grDevices::png(filename = repl.m.s.enrichmnt.png, width = 800, height = 800/2)
    } else {
      grDevices::cairo_ps(repl.m.s.enrichmnt.eps,width = 10,height = 10/2)
    }
    # Points with error bars (SD) per condition
    par(mar=c(5.1,4.1,0.5,0.5))
    graphics::plot( 0, col="white", axes=FALSE, xlab="codon", ylab="Codon Enrichment (-)",
         xlim=c(0,nrow(e.M)+1),
         ylim=c( min(0,min(e.M-e.SD)), (max(e.M+e.SD)+stats::median(unlist(e.SD)))  ) )
    graphics::par(las=2) ; graphics::axis(1, at=nCond/2*0.1 + 1:nrow(e.M), labels=row.names(e.M), family="Courier New")
    graphics::par(las=1) ; graphics::axis(2)
    graphics::abline(v=nCond/2*0.1/2 + 1:nrow(e.M) , col="lightgrey")
    graphics::lines(x=c(1,nrow(e.M)), y=c(1,1) , col="lightgrey")
    for (i in 1:nCond) {
      # Actual plot of values and Std-error bars
      xi <- (i-1)*0.1 + 1:nrow(e.M)
      graphics::points(xi, e.M[,i], col=col2trspcol(couleurs[i]), pch=20)
      graphics::arrows(x0= xi, x1= xi, y0= e.M[,i] -e.SD[,i], y1= e.M[,i] +e.SD[,i], code=0, col=col2trspcol(couleurs[i]), lwd=2.5)
    }
    #par(xpd=TRUE) # for plotting legend in margins, not on plot.
    graphics::legend("bottom", legend=conditions, pch=20, col=couleurs[1:nCond], bty="y", bg="white", horiz=TRUE, cex=0.8)
    grDevices::dev.off()
    
  }
  
  # Sorted codon plot per comparison
  iComp <- 0
  if (nCond>1) {
    for (i in 1:(nCond-1)) {
      for (j in (i+1):nCond) {
        ci1  <- conditions[i]
        ci2  <- conditions[j]
        iComp <- iComp + 1
        #
        fci.1.2 <- paste(pathout, "Condition_",ci1,"_relative-to_",ci2,"_mean-and-stderr.tsv", sep="")
        q.m.se <- utils::read.table(fci.1.2, sep="\t", header=TRUE, check.names=FALSE, row.names=1)
        q.m.se.sort <- q.m.se[order(-q.m.se$mean),]
        # Plot
        compar.codon.enrichmnt.png <- paste(pathout, "Visu-enrichmnt_compare-condition_",ci1,"_relative-to_",ci2,".png",sep="")
        compar.codon.enrichmnt.eps <- paste(pathout, "Visu-enrichmnt_compare-condition_",ci1,"_relative-to_",ci2,".eps",sep="")
        for (plotformat in c("png","eps")) {
          if (plotformat=="png") {
            grDevices::png(filename = compar.codon.enrichmnt.png, width = 800, height = 800/2)
          } else {
            grDevices::cairo_ps(compar.codon.enrichmnt.eps,width = 10,height = 10/2)
          }

          graphics::plot( 0, col="white", axes=FALSE, xlab="codon", ylab="Relative Codon Enrichment (-)",
               main=paste("Condition",ci1,"relative to ",ci2),
               xlim=c(0,nrow(q.m.se.sort)+1),
               ylim=c( min(q.m.se.sort$mean-q.m.se.sort$stderr), 
                       max(q.m.se.sort$mean+q.m.se.sort$stderr) ) )
          graphics::par(las=2) ; graphics::axis(1, at=1:nrow(q.m.se.sort), labels=row.names(q.m.se.sort), family="Courier New")
          graphics::par(las=1) ; graphics::axis(2)
          graphics::abline(h=1 , col="lightgrey")
          xi <- 1:nrow(e.M)
          #
          graphics::points(xi, q.m.se.sort$mean, pch=20, col=couleurs[nCond+iComp])
          graphics::arrows(x0= xi, x1= xi, y0= q.m.se.sort$mean-q.m.se.sort$stderr, y1= q.m.se.sort$mean+q.m.se.sort$stderr, code=0, col=couleurs[nCond+iComp], lwd=2.5)
          grDevices::dev.off()
        }
      }
    }
    # Return only 1 comparison
    ci1.ret <- conditions[1]
    ci2.ret <- conditions[2]
    compar.codon.enrichmnt.png.ret <- paste(pathout, "Visu-enrichmnt_compare-condition_",ci1.ret,"_relative-to_",ci2.ret,".png",sep="")
  } else {
    graphics::legend("center", legend="No comparison since 1 condition only.", bty="n")
  }
  # Save values for later reference, and return
  visu.m.s.enrichmnt.res <- list(plot.enrich=repl.m.s.enrichmnt.png, 
                                 plot.compar=compar.codon.enrichmnt.png.ret)
  save(visu.m.s.enrichmnt.res, 
       file=paste(pathout, "valNrec_", "visu.m.s.enrichmnt.res",".RData", sep=""))
  return(visu.m.s.enrichmnt.res)
}




# visu.tracks : Plot tracks of 1 mRNA over several samples
#
#
#
# This function takes a list of samples and experimental conditions in input and
#   plots coverage in A site.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#        refCDS         Address of file containing coding sequence
#                                  annotation. This file should contain tab-separated
#                                  values for mRNA name, nt position of start 
#                                  codon, nt position of first codon after stop 
#                                  codon, for example : 
#                                       ID	localStart	localEnd
#                                       NM_001276351	145	1348
#                                       NM_001276352	145	799
#                                       NM_000299	251	2495
#        mRNA           Name of mRNA for which to plot tracks, defaults to
#                       "random"
#        codon.labels   Indicates if codon identity should be added on the plot,
#                         defaults to FALSE
#        codon.col      Color by codon identity, defaults to "darkslateblue"
#
# Out:
#       This function returns a list containing the following : 
#         - plot        Address of png file containing track plot
#
visu.tracks <- function(XP.conditions, XP.names, pathout, refCDS, mRNA="random", 
                          codon.labels=FALSE, codon.col="darkslateblue") {

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))
  
  ## Check presence of mRNA, or pick a random mRNA from the upper quartile 
  nbr.pre <- utils::read.table(paste(pathout, "Sample-",XP.names[1],".nbreads", sep=""), 
                    sep="\t", header=TRUE, row.names=1)
  if (mRNA == "random") {
    nbr <- nbr.pre[nbr.pre$read.count > stats::quantile(nbr.pre$read.count, 0.75),,drop=FALSE]                      
    mRNA <- sample(row.names(nbr),1)
  } else {
    nbr.mRNA <- sum(row.names(nbr.pre)==mRNA)
    if (nbr.mRNA == 0) print("There is no footprint corresponding to this mRNA, at least in sample 1.")
  }
  
  ## CDS length
  cds <- utils::read.table(refCDS, sep="\t", header=TRUE, row.names=1)
  cdslen <- ( cds[mRNA,"localEnd"] - cds[mRNA,"localStart"] )/3
  
  ## Plot tracks for each sample
  n <- length(XP.conditions)
  track <- c()
  for (i in 1:n) {
    track.i.pre <- utils::read.table(paste(pathout, "Sample-",XP.names[i],".SCO", sep=""), sep="\t", header=TRUE)
    track[[i]] <- track.i.pre[track.i.pre$mRNA==mRNA,c("mRNA","codonPosition", "codon", "coverage")]
  }
  
  # Plot in png and eps
  visu.tracks.png <- paste(pathout, "Visu-tracks_",mRNA,".png",sep="")
  visu.tracks.eps <- paste(pathout, "Visu-tracks_",mRNA,".eps",sep="")
  for (plotformat in c("png","eps")) {
    if (plotformat=="png") {
      grDevices::png(filename = visu.tracks.png, width = 800, height = 400)
    } else {
      grDevices::cairo_ps(visu.tracks.eps,width = 10,height = 10)
    }
  
    # Actual plot
    graphics::par(las=1,mar=c(4,4,0.5,4), xpd=TRUE) #2,4)) #1,4,0,1))
    #graphics::par(mfrow=c(n,1))
    ymax <- 1.25 * max(sapply(track, function(L) {max(L$coverage)})) 
    y.gap <- 1.5 * max(sapply(track, function(L) {max(L$coverage)})) 
    y.max <- n*y.gap
      
    # Fake first plot window, just in order to retrieve proper y-axis annotation
    graphics::plot(c(0,cdslen), c(0,0), type='l', col=NULL, 
                   xlim=c(0,cdslen), ylim=c(0,ymax), xlab="", ylab="", axes=FALSE)
    axis.y <- axis(side=2, lwd=0,tick=FALSE, labels=FALSE, cex.axis=0.8)
    
    # Plot window
    graphics::plot(c(0,cdslen*1.), c(0,0), type='l', col=NULL,
                   xlim=c(0,cdslen*1.), ylim=c(0,y.max), 
                   xlab=paste("mRNA : ", mRNA, sep=""), 
                   ylab="Coverage in A site", axes=FALSE)
    # One track per sample
    for (i in 1:n) {
      # Prepare
      #graphics::plot(c(0,cdslen), c(0,0), type='l', xlim=c(0,cdslen), ylim=c(0,ymax), xlab=paste("mRNA : ", mRNA, sep=""), ylab="Coverage in A site", axes=TRUE)
      y.offset <- (n-(i))*y.gap
      graphics::lines(c(0,cdslen), y.offset+c(0,0), type='l')
      all.x <- track[[i]]$codonPosition
      all.y <- track[[i]]$coverage
      all.c <- track[[i]]$codon

      # Map colors according to codons, for current sample
      if (all(codon.col == "darkslateblue")) {
        codoncol = rep("darkslateblue", length(all.x))
      } else {
        codoncol <- codon.col[match(all.c,codon.col[,1]),2]
      }
      
      # Plot elements, for current sample
      for (ii in 1:length(all.x)) { 
        graphics::polygon(x = all.x[ii]+c(0,0,1,1, 0), #*2.5 width 1nt because line width is otherwise too large2,2, 0),
                y = y.offset + (all.y[ii]*c(0,1,1,0, 0)), 
                border = NA, #rep(C5lq[i], 5), 
                col=codoncol[ii])
      }
      #graphics::legend("topright", XP.names[i], bty='n')
      graphics::text(cdslen*1.01, y.offset,labels=XP.names[i], adj=c(0,0), cex=0.8)
      
      if (codon.labels) {
        graphics::text(all.x, y.offset + all.y+0.2, labels=track[[i]]$codon, 
                         srt=90, col="lightgrey", adj=0)
      }
      
      # Plot y-axis for current sample
      axis(side=2, at= (y.offset + axis.y), labels=axis.y, cex.axis=0.8 )
      
    } # End loop on n
    
    # Plot X-legend under last plot
    axis.x <- axis(side=1, lwd=0,tick=FALSE, col.ticks=NULL, col=NULL, col.axis=NULL,labels=FALSE)
    if (cdslen > axis.x[length(axis.x)]) {
      axis(side=1, at=c(axis.x,cdslen), cex.axis=0.8)
    } else {
      axis(side=1, at=c(axis.x[1:(length(axis.x)-1)],cdslen), cex.axis=0.8)
    }
    #axis(side=1, at=c(0,cdslen),labels=c("","CDSlength"),cex.axis=0.8,line=1,lwd=0)
    axis(side=1, at=c(cdslen),labels=c("CDSlength"),cex.axis=0.8,line=1,lwd=0)

    grDevices::dev.off() 
    
  } #End loop on plotformat
  
  # Save values for later reference, and return
  visu.tracks.res <- list(plot=visu.tracks.png)
  save(visu.tracks.res, 
       file=paste(pathout, "valNrec_", "visu.tracks.res",".RData", sep=""))
  return(visu.tracks.res)
}




# enricht.aroundA : Plot enrichment of all codons in large [-90 ; +90] or closed in [-15 ; +15]
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

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))

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
    ce.m <- read.table(paste(pathout, "Sample-",XP.names[isample],"_Codon-enrichment-unbiased_mean",sep=""), 
                        header=TRUE, row.names = 1)
    
    ### Filter out STOP codons and reorder
    ce.m.2 <- ce.m[ !(row.names(ce.m) %in% c("uaa","uag","uga")), ]
    ce.m.3 <- ce.m.2[codgen.order,]
    
      
    ### Plot all, separately, in+/-90 window around A
    # Prepare bounds for Y-axis
    ymin <- min( 0.5 , min(ce.m.3, na.rm=TRUE) )
    ymax <- max( 1.8 , max(ce.m.3, na.rm=TRUE) )
    # Filenames
    enrichment.all.png <-         paste(pathout, "enrichment-all_",XP.names[isample],".png", sep="")
    enrichment.all.extract.png <- paste(pathout, "enrichment-all-extract_",XP.names[isample],".png", sep="")
    enrichment.all.eps <-         paste(pathout, "enrichment-all_",XP.names[isample],".eps", sep="")
    # Graphical devices
    for (plotformat in c("png","shortpng","eps")) {
        if (plotformat=="png") {
          grDevices::png(file = enrichment.all.png,  width = 400, height = 1600)
          nb.row   <- 17
          nb.codon <- 64
        } else if (plotformat=="shortpng") {
          # Preview of enrichment.all, for inclusion in html report, with 4*4 codons instead of all.
          grDevices::png(file = enrichment.all.extract.png,  width = 400, height = 400)
          nb.row   <- 5
          nb.codon <- 16
        } else {
          grDevices::cairo_ps(enrichment.all.eps,width = 5,height = 20)
          nb.row   <- 17
          nb.codon <- 64
        }
        # Divide plot window in 64 boxes, +1 line and +1 column for axes 
        #   (except for 'shortpng' : 16 boxes)
        par(mfrow=c(nb.row,5), las=1, mar=0.2*c(1,1,1,1))
        # Loop over codons : 64 (except for shortpng : 16)
        for (i in 1:nb.codon) { 
          # Y-axis if first box in the row of plots
          if (i%%4 == 1) {
            plot(-90:90, rep(NA,181), col=NULL, 
               ylim=c(ymin,ymax), xlim=c(-90,90), axes=FALSE, main='', xlab='', 
               ylab='', frame=FALSE)
            axis(side=2, line=-7, cex.axis=0.6)
            par(las=3) ; mtext(text='Enrichment (-)', side=2, line=-4.9, cex=0.5) ; par(las=1)
          }        
          # Enrichment line plot
          if (codgen.order[i] != "uaa" & codgen.order[i] != "uag" & codgen.order[i] != "uga") {
          plot(-90:90, ce.m.3[i,], #/median(as.numeric(ce.m.3[i,])), 
               col=couleurs3, type='l', ylim=c(ymin,ymax),
               xlim=c(-90,90), axes=FALSE, main='', xlab='', ylab='', frame=TRUE)
          legend("topleft",legend=codgen.order[i], bty="n", adj=c(0.8,0))
          abline(h=1, lwd=0.5)
          # A-site line
          abline(v=0, lty=2, col='grey25', lwd=0.8) 
          # Enrichment again, on top of lines plotted inbetween
          lines(-90:90, ce.m.3[i,], col=couleurs3)  
          # 3' and 5' annotation
          #text(0,ymin, adj=0.2, labels="A", cex=0.5)
          text(-90,ymin,adj=0,labels="3'", cex=0.5) ; text(+90,ymin,adj=1,labels="5'", cex=0.5)
          } else {plot(0, 0, col="white", type='l', ylim=c(0,3),
               xlim=c(-15,15), axes=FALSE, main='', xlab='', ylab='', frame=TRUE)
               text(x=0, y=1.5, labels="STOP codon", cex=0.7)
          }
        }
        # One blank box bottom left
        plot(-90:90,rep(NA,181), col=NULL, ylim=c(ymin,ymax),
             xlim=c(-90,90), axes=FALSE, main='', xlab='', ylab='', frame=FALSE)
        
        # X-axis below each column of plot boxes
        for (ii in 1:4) {
          plot(-90:90,rep(NA,181), col=NULL, ylim=c(ymin,ymax),
               xlim=c(-90,90), axes=FALSE, main='', 
               xlab='', 
               ylab='', frame=FALSE)
          axis(side=1, line=-8.3, at=seq(-90,90, by=45), cex.axis=0.6)
          # A-site annotation ; note : hadj 0.4, not 0.5, for correct centering (could be platform or version-dependent, though)
          axis(side=1, line=-8.9, at=c(0), labels=c("A"), tick=FALSE, hadj=0.4, 
                 cex.axis=0.6)
          mtext(text="Codon position", 
                  side=1, line=-6.3, cex=0.5)
        }
      dev.off()
    }
    
    ### Plot enrichment of all codons, separately, closed in on ribosome [-15 ; +15]
    # Prepare bounds of Y-axis
    ymin <- min( 0.5 , min(ce.m.3[91+c(-15:15)], na.rm=TRUE) )
    ymax <- max( 1.8 , max(ce.m.3[91+c(-15:15)], na.rm=TRUE) )
    # Filenames
    enrichment.all.zoom.png         <- paste(pathout, "enrichment-all-zoom_",XP.names[isample],".png", sep="")
    enrichment.all.zoom.extract.png <- paste(pathout, "enrichment-all-zoom-extract_",XP.names[isample],".png", sep="")
    enrichment.all.zoom.eps         <- paste(pathout, "enrichment-all-zoom_",XP.names[isample],".eps", sep="")
    # Graphical devices
    for (plotformat in c("png","shortpng","eps")) {
        if (plotformat=="png") {
          grDevices::png(file = enrichment.all.zoom.png,  width = 400, height = 1600)
          nb.row   <- 17
          nb.codon <- 64
        } else if (plotformat=="shortpng") {
          # Preview with less codon, for inclusion in html report
          grDevices::png(file = enrichment.all.zoom.extract.png,  width = 400, height = 400)
          nb.row   <- 5
          nb.codon <- 16
        } else {
          grDevices::cairo_ps(enrichment.all.zoom.eps,width = 5,height = 20)
          nb.row   <- 17
          nb.codon <- 64
        }
        # Divide in 64 boxes, +1 line and +1 column for axes
        #   except for shortpng : 16 codons, +1 line, +1 column
        par(mfrow=c(nb.row,5), las=1, mar=0.2*c(1,1,1,1))
        # Loop over the 64 codons (except shortpng : 16 codons
        for (i in 1:nb.codon) { 
          # Plot axis if first box in the row of plots
          if (i%%4 == 1) {
            plot(-90:90, rep(NA,181), col=NULL, 
               ylim=c(ymin,ymax), xlim=c(-90,90), axes=FALSE, main='', xlab='', 
               ylab='', frame=FALSE)
            axis(side=2, line=-7, cex.axis=0.6)
            par(las=3) ; mtext(text='Enrichment (-)', side=2, line=-4.9, cex=0.5) ; par(las=1)
          }        
          # Plot of enrichment
          if (codgen.order[i] != "uaa" & codgen.order[i] != "uag" & codgen.order[i] != "uga") {
          plot(-90:90, ce.m.3[i,], 
               col=couleurs3, type='l', ylim=c(ymin,ymax),
               xlim=c(-15,15), axes=FALSE, main='', xlab='', ylab='', frame=TRUE)
          legend("topleft",legend=codgen.order[i], bty="n", adj=c(0.8,0))
          abline(h=1, lwd=0.5)
          # A and P-sites lines
          abline(v=0, lty=2, col='grey25', lwd=0.8) 
          abline(v=c(1), col="lightgrey", lty=2, lwd=0.8) 
          # Enrichment again, on top of lines
          lines(-90:90, ce.m.3[i,], col=couleurs3)  
          # A and P-sites, 3' and 5' annotation
          text(0,ymin, adj=1, labels="A", cex=0.5) ; text(1,ymin, adj=0, labels="P", cex=0.5)
          text(-15,ymin,adj=0,labels="3'", cex=0.5) ; text(+15,ymin,adj=1,labels="5'", cex=0.5)
          } else {plot(0, 0, col="white", type='l', ylim=c(0,3),
               xlim=c(-15,15), axes=FALSE, main='', xlab='', ylab='', frame=TRUE)
               text(x=0, y=1.5, labels="STOP codon", cex=0.7)
          }
        } # end of loop on codons
        # One blank box
        plot(-90:90,rep(NA,181), col=NULL, ylim=c(ymin,ymax),
             xlim=c(-90,90), axes=FALSE, main='', xlab='', ylab='', frame=FALSE)
        
        # Axes below each column of plot boxes
        for (ii in 1:4) {
          plot(-90:90,rep(NA,181), col=NULL, ylim=c(ymin,ymax),
               xlim=c(-15,15), axes=FALSE, main='', 
               xlab='', 
               ylab='', frame=FALSE)
          axis(side=1, line=-8.3, at=c(-15,-10,-5,0,5,10,15), 
                 labels=c("-15","","",0,"","",15),cex.axis=0.6)
          mtext(text="Codon position", 
                  side=1, line=-6.3, cex=0.5)
        }
      dev.off()
    } #end of loop on graphical devices
    
  } # End of loop on samples

  # Save values for later reference, and return
  enricht.aroundA.res <- list(plot=enrichment.all.png, 
                              plot.extract=enrichment.all.extract.png, 
                              plot.zoom=enrichment.all.zoom.png,
                              plot.zoom.extract=enrichment.all.zoom.extract.png)
  save(enricht.aroundA.res, 
       file=paste(pathout, "valNrec_", "enricht.aroundA.res",".RData", sep=""))
  return(enricht.aroundA.res)

}

