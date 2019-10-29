# Artefacts visible on codon enrichment, due to cycloheximide or other drug
#
#
#
# This function takes a list of samples and their conditions as input and 
#   visualizes enrichment around AUG +/- 90 codons, where possible artefacts due
#   to drugs used in the experiment should be visible.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#        algo           Algorithm used, either : Unbiased (default)
#                                           or : Hussmann
#
# Out:
#       A list containing:
#          - plot                Address of png plot file
#          - value               Standard deviation of enrichment for each codon
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#
chx.artefacts <- function(XP.conditions, XP.names, pathout, algo="unbiased") {

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))

  ## Read arguments and determine their number
  n     <- length(XP.conditions)

  ## colors
  couleurs.pre <- c('black','red','forestgreen','blue','orange',
                    'purple','gold','brown','magenta','midnightblue')
  couleurs <- rep(couleurs.pre,7)


  ### Arginine codons - Hussmann 3b
  # postscript(file = paste(pathout, "chx-artefacts-Arg_", j, ".eps", sep=""),
  #              onefile = TRUE,
  #              title = sample,
  #              width = 3.937, height = 3.937,
  #              paper = "special",
  #              horizontal=FALSE)
  maxArg <- 2
  png.chx.arg.png <- paste(pathout, "chx-artefacts-Arg.png", sep="")
  png.chx.arg.eps <- paste(pathout, "chx-artefacts-Arg.eps", sep="")

  for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(file = png.chx.arg.png,  width = 400, height = 400)
      } else {
        grDevices::cairo_ps(png.chx.arg.eps,width = 10,height = 10)
    }
    #grDevices::png(file = png.chx.arg,  width = 800, height = 800)
    graphics::par(mfrow=c(2,2), mar=c(4,4,2,2), las=1) #default mar = 5.1 4.1 4.1 2.1
    # cga
    graphics::plot(0, type='l', col=NULL, xlim=c(-90,90), ylim =c(0,maxArg),
                      xlab="3' to 5' codon position relative to A-site",
                      ylab="Mean codon enrichment",
                      main="cga") ; graphics::abline(h = 1)
    for (j in 1:n) {
      if (algo=="unbiased") {
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Codon-enrichment-unbiased_mean",sep=""), header=TRUE, row.names = 1)
        graphics::lines(-90:90, ce.m['cga',], type='l', col=couleurs[j]) 
      } else if (algo=="Hussmann") {
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Hussmann_mean",sep=""), header=TRUE, row.names = 1)
        mednorm <- stats::median(as.numeric(ce.m['cga',]))
        if (mednorm > 0) { graphics::lines(-90:90, ce.m['cga',]/mednorm, type='l', col=couleurs[j]) }
      } else {
        print("Option 'algo' should be either 'unbiased' or 'Hussmann'. ")
        print("Using default option.")
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Codon-enrichment-unbiased_mean",sep=""), header=TRUE, row.names = 1)
        graphics::lines(-90:90, ce.m['cga',], type='l', col=couleurs[j]) 
      }
    }
    #cgg
    graphics::plot(0, type='l', col=NULL, xlim=c(-90,90), ylim =c(0,maxArg),
                      xlab="3' to 5' codon position relative to A-site",
                      ylab="Mean codon enrichment",
                      main="cgg") ; graphics::abline(h = 1)
    for (j in 1:n) {
      if (algo=="unbiased") {
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Codon-enrichment-unbiased_mean",sep=""), header=TRUE, row.names = 1)
        graphics::lines(-90:90, ce.m['cgg',], type='l', col=couleurs[j]) 
      } else if (algo=="Hussmann") {
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Hussmann_mean",sep=""), header=TRUE, row.names = 1)
        mednorm <- stats::median(as.numeric(ce.m['cgg',]))
        if (mednorm > 0) { graphics::lines(-90:90, ce.m['cgg',]/mednorm, type='l', col=couleurs[j]) }
      } else {
        print("Option 'algo' should be either 'unbiased' or 'Hussmann'. ")
        print("Using default option.")
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Codon-enrichment-unbiased_mean",sep=""), header=TRUE, row.names = 1)
        graphics::lines(-90:90, ce.m['cgg',], type='l', col=couleurs[j]) 
      }
    }
    # aga
    graphics::plot(0, type='l', col=NULL, xlim=c(-90,90), ylim =c(0,maxArg),
                      xlab="3' to 5' codon position relative to A-site",
                      ylab="Mean codon enrichment",
                      main="aga") ; graphics::abline(h = 1)
    for (j in 1:n) {
      if (algo=="unbiased") {
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Codon-enrichment-unbiased_mean",sep=""), header=TRUE, row.names = 1)
        graphics::lines(-90:90, ce.m['aga',], type='l', col=couleurs[j]) 
      } else if (algo=="Hussmann") {
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Hussmann_mean",sep=""), header=TRUE, row.names = 1)
        mednorm <- stats::median(as.numeric(ce.m['aga',]))
        if (mednorm > 0) { graphics::lines(-90:90, ce.m['aga',]/mednorm, type='l', col=couleurs[j]) }
      } else {
        print("Option 'algo' should be either 'unbiased' or 'Hussmann'. ")
        print("Using default option.")
        ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Codon-enrichment-unbiased_mean",sep=""), header=TRUE, row.names = 1)
        graphics::lines(-90:90, ce.m['aga',], type='l', col=couleurs[j]) 
      }
    }
    graphics::par(mar=c(4,0,1,2))
    graphics::plot(0, type='l', col=NULL, xlim=c(0,1), ylim =c(0,1), xlab="",ylab="", axes=FALSE)
    graphics::legend("top", legend = XP.names, col = couleurs[1:n], lty = 1, bty='n', bg = NULL)
    graphics::legend('bottom', legend = "Samples with many missing \n points are not displayed.", bty='n', bg=NULL)
  grDevices::dev.off()
  } #end loop PNG and EPS

  sd.enrich <- c()
    for (j in 1:n) {
      ce.m <- utils::read.table(paste(pathout,"Sample-",XP.names[j],"_Hussmann_mean",sep=""), header=TRUE, row.names = 1)
      sd.enrich[[j]] <- apply(ce.m,1,function(x) {
                                          mednorm <- stats::median(as.numeric(x))
                                          if (mednorm > 0) {
                                            sd.normed.BCO <- stats::sd( as.numeric(x) / mednorm )
                                          } else {
                                            sd.normed.BCO <- NA
                                            } })
    }

  sd.enrich.mat <- matrix(unlist(sd.enrich), nrow=64, ncol=n)
  colnames(sd.enrich.mat) <- XP.names
  rownames(sd.enrich.mat) <- rownames(ce.m)

  max.sd.normed <- max(sd.enrich.mat, na.rm=TRUE)

  if (max.sd.normed > 1) {
    couleur <- "red"
    recommendation <- paste("Some codons have standard deviation >1 after normalisation.\n",
                            " This figure is plotted using ",algo," enrichment (plots for each\n",
                            " sample and each codon are available in the output directory :\n",
                            " ",pathout, ").", sep="") 
  } else if (max.sd.normed > 0.5) {
    couleur <- "orange"
    recommendation <- paste("Some codons have standard deviation >0.5 after normalisation.\n",
                            " This figure is plotted using ",algo," enrichment (plots for each\n",
                            " sample and each codon are available in the output directory :\n",
                            " ",pathout, ").", sep="") 
  } else {
    couleur <- "green"
    recommendation <- paste("All codons have standard deviation <0.5 after normalisation.\n",
                            " This figure is plotted using ",algo," enrichment (plots for each\n",
                            " sample and each codon are available in the output directory :\n",
                            " ",pathout, ").", sep="") 
  }

  # Save values for later reference, and return
  chx.artefacts.res <- list(plot=png.chx.arg.png,
                            value=sd.enrich.mat,
                            color=couleur,
                            recommendation=recommendation)
  save(chx.artefacts.res, 
       file=paste(pathout, "valNrec_", "chx.artefacts.res",".RData", sep=""))
  return(chx.artefacts.res)
}


