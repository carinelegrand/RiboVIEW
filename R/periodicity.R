# Functions related to periodicity around AUG (periodicity and recurrence plot)
#
# Contents : 
#            periodicity 
#            periodicity_one
#            periodicity_plot
#            select.FPlen
#

# periodicity : Generate coverage counts around AUG for each footprint length.
#
# This function takes a list of bam files as input and annotations in input,
#   generate coverage counts around AUG and produces recurrence and periodicity 
#   plots.
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
#        pathout                 Address where output files will be written.
#        XP.names                Vector of names for each sample
#        versionStrip            Indicates if version number should be trimmed, 
#                                  for example NM_001276351.1 should be trimmed 
#                                  to NM_001276351 (defaults to FALSE)
#
# Out: 
#        The following files are written to pathout :
#           - OL-*.txt           Coverage around AUG for each sample
#           - OL-*_recurr.png    Periodicity recurrence plot for each sample
#           - OL-*_period.png    Periodicity barplot for each sample
#
periodicity <- function(list.bam, refCDS, refFASTA, pathout, XP.names, 
                        versionStrip = FALSE, python.messages = TRUE) {

  ## Read arguments and determine their number
  argl  <- list.bam
  n     <- length(argl)

  for (i in 1:n) {
    readsBAM <- argl[[i]]
    outfile <- paste(pathout, "Sample-", XP.names[i], "_OL.txt" , sep="")
    periodicity_one(readsBAM, refCDS, refFASTA, outfile, versionStrip, 
                    python.messages)

    # Plots :
    # Recurrence plot and Periodicity plot
    periodicity_plot(outfile, i, pathout, XP.names)

  }
}

periodicity_one <- function(readsBAM, refCDS, refFASTA, outfile, 
                            versionStrip = FALSE,
                            python.messages = TRUE) {
  
  pyPeriodicity_filename <- system.file("periodicity.py", package="RiboVIEW", mustWork = TRUE)
  PythonInR::pyExecfile(pyPeriodicity_filename)
  pyPeriodicity <- PythonInR::pyFunction("periodicity")
  pyPeriodicity(readsBAM, refCDS, refFASTA, outfile, 
                versionStrip=FALSE, 
                messages=python.messages)
  #rPython::python.load(system.file("periodicity.py", package="RiboVIEW", mustWork = TRUE))
  #rPython::python.call("periodicity", readsBAM, refCDS, refFASTA, outfile, versionStrip)

  }

##### PLOTS ######
#
#!/bin/R

periodicity_plot <- function(outfile, i, pathout, XP.names) {
  
  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))

  # inputs
  nom <- outfile
  
  OL <- utils::read.table(nom, header = TRUE, row.names = 1, as.is = TRUE)
  colnames(OL) <- seq(25,32)


  #
  ## recurrence plot
  #
  for (pos in seq(25,32)) {

    OL.pos.ts <- stats::ts(data=as.numeric(as.vector(OL[,colnames(OL)==pos])), start=-20, end=20, frequency=1)         #deltat=1/3)

    plot.rec.tmp.png <- paste(pathout, "Sample-", XP.names[i], "_OL-recurr_", pos, ".png" , sep="")
    plot.rec.tmp.eps <- paste(pathout, "Sample-", XP.names[i], "_OL-recurr_", pos, ".eps" , sep="")
    #
    for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = plot.rec.tmp.png, width = 500, height = 500)  
      } else {
        grDevices::cairo_ps(plot.rec.tmp.eps,width = 6.25,height = 6.25)   
      }
      graphics::par(mar=c(4,4,2,0))
      CLrecurr(OL.pos.ts, m=4, d=0,nlevels=40, plot.axes = { graphics::axis(1, seq(-18,18,by=3))
                                                             graphics::axis(2, seq(-18,18,by=3)) },
                maincolor="black", invcol=TRUE, pos.name = pos)
      # Close last graphic window, checking if still open
      grDevices::dev.off()

    }
  }
  plot.rec.png <- paste(pathout, "Sample-", XP.names[i], "_OL-recurr.png" , sep="")
  plot.rec.eps <- paste(pathout, "Sample-", XP.names[i], "_OL-recurr.eps" , sep="")
  system(paste("montage -mode concatenate -tile 4x2 ",
             pathout, "Sample-", XP.names[i], "_OL-recurr_*.png ",
             plot.rec.png,
             sep=""))
  system(paste("montage -mode concatenate -tile 4x2 ",
             pathout, "Sample-", XP.names[i], "_OL-recurr_*.eps ",
             plot.rec.eps,
             sep=""))

  #
  ## Bar plot for periodicity
  #
  plot.per.png <- paste(pathout, "Sample-", XP.names[i], "_OL-period.png" , sep="")
  plot.per.eps <- paste(pathout, "Sample-", XP.names[i], "_OL-period.eps" , sep="")
  for (plotformat in c("png","eps")) {
    if (plotformat=="png") {
      grDevices::png(filename = plot.per.png, width = 1900, height = 1000)  
    } else {
      grDevices::cairo_ps(plot.per.eps,width = 23.75,height = 12.50)   
    }
  #grDevices::png(filename = plot.per, width = 1900,height = 1000)
    graphics::par(mar=c(8,4,1,1),las=1)
    maxloc <- max(as.numeric(as.matrix(OL)))
    graphics::plot(row.names(OL), as.vector(OL[,1]),type='h', xlim=c(-20,7*45+41),ylim=c(0,maxloc),
         col='red',axes=FALSE,xlab='',ylab='coverage')
    graphics::axis(2)
    graphics::par(las=3)
    graphics::axis(1,at = seq(-18,18,by=9),cex.axis=0.5)
    graphics::par(las=1)
    graphics::mtext(1, text='position / A-site', line=1.5)
    graphics::axis(1, at= 45*0:7, line = 3.1, labels=colnames(OL))
    graphics::mtext(1,text = 'footprint length', line = 5.1)
    graphics::par(las=3)
    for (i in 2:8) {
      graphics::points((i-1)*45+ as.numeric(row.names(OL)),
             as.vector(OL[,i]),type='h',col='red')
      graphics::axis(1, at = (i-1)*45 + seq(-18,18,by=9),labels = seq(-18,18,by=9),cex.axis=0.5)
    }
  grDevices::dev.off()
  }
  }
  
  
# select.FPlen : Display periodicity plots for user selection of admissible footprint lengths.
#
# This function takes a list of bam files as input and annotations in input,
#   displays periodicity diagnostic plots, and invites the user to select admissible
#   footprint lengths for subsequent calculations and plots. 
#
# In:
#        list.bam                List of bam files containing aligned reads for 
#                                  each sample (same order as in XP.names)
#        pathout                 Address where output files will be written.
#        XP.names                Vector of names for each sample
#
# Out: 
#        A list containing the following :
#           - plot.rec  Periodicity recurrence plot, address of the plot file in png format
#           - plot.cov  Periodicity barplot, address of the plot file in png format
#           - mini      Minimum footprint length to consider, as selected by the user
#           - maxi      Maximum footprint length to consider, as selected by the user
#
select.FPlen <- function(list.bam, pathout, XP.names) {

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))

  # Prompt user to select footprint length and give instructions how to do so
  cat(" Please select footprint lengths with - sufficient data, \n")
  cat("                                      - sufficient periodicity, and\n")
  cat("                                      - peak at -12 nt. \n")
  cat(" (Expected peak at -12 not -15 since transition A-P at initiation is almost immediate.)\n")
  
  name.plot.timestamp <- paste("recurrPlots__",
                               format(Sys.time(), "%Y_%m_%d__%k_%M_%S"),
                               sep="")

  for (i in 1:length(list.bam)) {
    # Fetch image from png file generated earlier, and display in specific graphic window.
    im <- png::readPNG(paste(pathout,"Sample-",XP.names[i],"_OL-recurr.png",sep=""))
    R.devices::devSet(name.plot.timestamp, width=9, height=7)
    graphics::par(mar=c(0,0,0,0))
    graphics::plot(1:2, type='n',axes=FALSE, xlab="", ylab="")
    #rasterImage(im, 1.2, 1.27, 1.8, 1.73)
    graphics::rasterImage(im, 1, 1.25, 2, 1.75)
    graphics::text(1.5,1.8, XP.names[i], cex=2)
    invisible(readline(prompt="Note acceptable footprint lengths, then press [enter] to continue."))  
  }
  
  # Setting mini and maxi, or resorting to default values
  mini <- as.numeric(invisible(readline(
       prompt=paste("Enter minimum acceptable footprint length \n",
                    " (this will be applied to all samples), \n",
                    " then press [enter] to continue :         ",sep=""))))
  maxi <- as.numeric(invisible(readline(
       prompt=paste("Similarly, enter MAXimum acceptable footprint length \n",
                    " then press [enter] to continue :         ",sep=""))))
  # Defaults to 25 and 32 if not set, and print useful information about it :
  if (is.na(mini)) { cat("Minimum was not set. Using default value (25nt) :\n")
                     mini <- 25 }
  if (is.na(maxi)) { cat("Maximum was not set. Using default value (32nt) :\n") 
                     maxi <- 32 }
  cat(paste("The interval [",mini,",",maxi,"] of valid footprint lengths has been selected.\n",sep=""))
  cat("    The bounds of this interval are stored in variables 'mini' and 'maxi'.\n")
    
  # Close graphic window.
  R.devices::devOff(name.plot.timestamp)
  
  # Save
  periodicity.res <- list(plot.rec = paste(pathout,"Sample-",XP.names[1],"_OL-recurr.png",sep=""), 
                          plot.cov = paste(pathout,"Sample-",XP.names[1],"_OL-period.png",sep=""),
                          mini = mini,
                          maxi = maxi)
  save(periodicity.res,file=paste(pathout,"valNrec_","periodicity",".RData",sep=""))
  
  return(list(mini=mini, maxi=maxi))
}

