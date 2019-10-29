# Replicates consistency
#
# Contents
#           repl.correl.counts.Venn
#           repl.correl.gene
#           covA
#           repl.correl.codon.singleres
#           repl.correl.codon
#           BCO.norm
#           repl.correl.heatmap

# repl.correl.counts.Venn : Venn Diagram for mRNA per sample
#
#
#
# This function takes a list of samples and experimental conditions in input,
#   calculates counts of mRNA present in several samples, and presents these in a
#   Venn diagram.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#
# Out:
#       This function returns a list containing the following
#          - plot   Address of the Venn diagram plot in png format
#
repl.correl.counts.Venn <- function(XP.conditions, XP.names, pathout, r.messages=TRUE) {

  couleurs2.pre <- c('red', 'gold', 'orange', 'blue', 'magenta', 'green', 'purple')
  couleurs2 <- rep(couleurs2.pre,7)

  # Nb conditions and Nb replicates per condition
  conditions <- unique(XP.conditions)
  nCond      <- length(conditions)
  
  # One Venn Diagram per condition
  for (i in conditions) {
    idces <- which(XP.conditions==i)
    nb.idces <- length(idces)
    if (nb.idces > 5) {
      if (r.messages) {
      cat(paste("Plotting of a Venn diagram over more than 5 replicates is not supported \n",
                "  (limitation in VennDiagram package).\n",
                "  -> only the first 5 replicates are shown.\n"))
      }
    }
    # list of list
    # i.e. list by replicate of list of mRNAs in this replicate
    lli <- c()
    for (i.idces in 1:min(5,nb.idces)) {
      ii <- idces[i.idces]
      fi <- paste(pathout, "Sample-",XP.names[ii],".nbreads", sep="")
      li <- utils::read.table(fi, sep="\t", header=TRUE)$mRNA
      ni <- XP.names[ii]
      lli[[ni]] <- li
    }
    
    # Store venn diagram for later plot
    temp <- CLvenn.diagram(lli,
                  fill = couleurs2[1:min(5,nb.idces)], alpha = rep(0.5, min(5,nb.idces)), 
                  cex = 1, cat.fontface = 4,
                  lty =0, fontfamily =3, filename = NULL)
    
    # Generate venn diagram (per condition i) in png or eps format
    for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = paste(pathout,"Replicates-venn-",i,".png",sep=""),
            width = 800/2, height = 800/2)
      } else {
        grDevices::cairo_ps(paste(pathout,"Replicates-venn-",i,".eps",sep=""),
            width = 10/2,height = 10/2)
      }
      # Actual plot
      grid::grid.draw(temp)
      #
      grDevices::dev.off()
    }
  }
  repl.correl.counts.Venn.res <- list(plot=paste(pathout,"Replicates-venn-",i,".png",sep=""))
  # Save values for later reference, and return
  save(repl.correl.counts.Venn.res, 
       file=paste(pathout, "valNrec_", "repl.correl.counts.Venn",".RData", sep=""))

  # return(list(plot=paste(pathout,"Replicates-venn-",i,".png",sep="")))
  return(repl.correl.counts.Venn.res)
}

# repl.correl.counts.Venn.all : Venn Diagram for mRNA per sample, all samples all conditions
#
#
#
# This function takes a list of samples and experimental conditions in input,
#   calculates counts of mRNA present in all samples, and presents these in a
#   Venn diagram.
#
# In:
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#
# Out:
#       This function returns a list containing the following
#          - plot   Address of the Venn diagram plot in png format
#
Venn.all <- function(XP.names, pathout, r.messages=TRUE) {

  couleurs2.pre <- c('red', 'gold', 'orange', 'blue', 'magenta', 'green', 'purple')
  couleurs2 <- rep(couleurs2.pre,7)

  # Venn Diagram all samples
    nb.idces <- length(XP.names)
    if (nb.idces > 5) {
      if (r.messages) {
      cat(paste("Plotting of a Venn diagram over more than 5 replicates is not supported \n",
                "  (limitation in VennDiagram package).\n",
                "  -> only the first 5 samples are shown.\n"))
      }
    }
    # list of list
    # i.e. list by replicate of list of mRNAs in this replicate
    lli <- c()
    for (i.idces in 1:min(5,nb.idces)) {
      ii <- i.idces
      fi <- paste(pathout, "Sample-",XP.names[ii],".nbreads", sep="")
      li <- utils::read.table(fi, sep="\t", header=TRUE)$mRNA
      ni <- XP.names[ii]
      lli[[ni]] <- li
    }
    
    # Store venn diagram for later plot
    temp <- CLvenn.diagram(lli,
                  fill = couleurs2[1:min(5,nb.idces)], alpha = rep(0.5, min(5,nb.idces)), 
                  cex = 1, cat.fontface = 4,
                  lty =0, fontfamily =3, filename = NULL)
    
    # Generate venn diagram (per condition i) in png or eps format
    for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = paste(pathout,"Venn-all.png",sep=""),
            width = 800/2, height = 800/2)
      } else {
        grDevices::cairo_ps(paste(pathout,"Venn-all.eps",sep=""),
            width = 10/2,height = 10/2)
      }
      # Actual plot
      grid::grid.draw(temp)
      #
      grDevices::dev.off()
    }
  Venn.all.res <- list(plot=paste(pathout,"Venn-all.png",sep=""))
  # Save values for later reference, and return
  save(Venn.all.res, 
       file=paste(pathout, "valNrec_", "Venn.all",".RData", sep=""))

  # return(list(plot=paste(pathout,"Replicates-venn-",i,".png",sep="")))
  return(Venn.all.res)
}


# repl.correl.gene : Correlation at gene resolution
#
#
#
# This function takes a list of samples and experimental conditions in input,
#   calculates correlations at gene (mRNA level) between replicates.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#        use0           Calculates correlation taking 0 values into account,
#                         defaults to FALSE
#
# Out:
#       A list containing:
#          - cor                 Correlation between replicates
#          - plot                Address of plot file in png format
#          - value               Median of correlation across different conditions
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#
repl.correl.gene <- function(XP.conditions, XP.names, pathout, use0=FALSE) {

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))
  
  
  #### additional data
  correl.gene.DiamentTuller2016 <- c(0.945,0.94,0.93,0.945,0.96,0.96,0.97,0.855,0.9,
                                     0.995,0.965,0.965,0.99,0.96,0.975,0.95,0.98,0.985,
                                     0.975,0.975,0.98,0.99,0.995,0.98,0.99)  
  ## Read arguments and determine their number
  n     <- length(XP.conditions)

  list.rpkm <- c()
  for (i in 1:n) {
    list.rpkm.files.i <- paste(pathout, "Sample-",XP.names[i],".RPKM", sep="")
    list.rpkm[[i]] <- utils::read.table(list.rpkm.files.i, header=TRUE,row.names = 1)
  }

  argl <- list.rpkm

  ID.union <- row.names(argl[[1]])
  for (i in 2:n) {ID.union <- union(ID.union, row.names(argl[[i]]))}

  i<-1
  ri  <- as.data.frame(argl[[i]])
  rii <- data.frame(rep1 = ri[ID.union,], row.names = ID.union)
  rpkm.all <- rii
  #
  for (i in 2:n) {ri  <- as.data.frame(argl[[i]])
                  rii <- data.frame(rep1 = ri[ID.union,], row.names = ID.union)
                  names(rii) = XP.names[i]
                  #names(rii) = paste("rep", as.character(i), sep="")
                  rpkm.all <- data.frame(rpkm.all, rii, row.names = ID.union)
  }
  colnames(rpkm.all)[1] <- XP.names[1]
  
  # correlation, gene level
  rpkm.all.emptyasNA <- rpkm.all
  rpkm.cor.onlycomplete <- stats::cor(rpkm.all.emptyasNA, method = "spearman", use = "complete.obs")

  #assume NA means no detected copies of a specific mRNA
  rpkm.all[is.na(rpkm.all)] <- 0

  rpkm.cor.incl0 <- stats::cor(rpkm.all, method = "spearman")

  rpkm.cor.onlycomplete
  rpkm.cor.incl0

  # Nb conditions and Nb replicates per condition
  conditions <- unique(XP.conditions)
  nCond      <- length(conditions)
  repV       <- apply(cbind(conditions), 1,
                 function(k) {sum(XP.conditions==k)})
  nRepM <- max(repV)

  # plot to file, gene level
  repl.correl.gene.png <- paste(pathout,"Replicates-gene.png",sep="")
  repl.correl.gene.eps <- paste(pathout,"Replicates-gene.eps",sep="")

  for (plotformat in c("png","eps")) {
    if (plotformat=="png") {
      grDevices::png(filename = repl.correl.gene.png, width = 800*(nRepM-1)/nCond, height = 800)
    } else {
      grDevices::cairo_ps(repl.correl.gene.eps,width = 10*(nRepM-1)/nCond,height = 10)
    }

  #grDevices::png(filename = repl.correl.gene.png, width = 800*(nRepM-1)/nCond, height = 800)
    resu <- c()
    if (!use0) {

        graphics::par(mfrow=(c(nCond, (nRepM-1))), las=1, mar=c(4,4.5,2,0.5))
        for (i in conditions) {
          idces <- which(XP.conditions==i)
          for (i.idces in 1:(length(idces)-1)) {
            ii <- idces[i.idces]
            jj <- idces[(i.idces + 1)]
            graphics::plot(log(rpkm.all.emptyasNA[,c(ii,jj)], base=10), pch=20, cex=0.8, col=grDevices::rgb(0.4, 0.4, 0.4, 0.5), cex.axis = 1.5, cex.lab = 1.5)
               rhotext <- paste("$\\rho_S = ",as.character(round(rpkm.cor.onlycomplete[ii,jj], digits = 4)), " $", sep="")
               resu <- c(resu, rpkm.cor.onlycomplete[ii,jj])
               graphics::text(x = log(min(rpkm.all.emptyasNA[,ii], na.rm = TRUE),base=10),
                    y = log(max(rpkm.all.emptyasNA[,jj], na.rm = TRUE),base=10),
                    labels=latex2exp::TeX(rhotext),   #TeX('$\\rho_S = $'),
                    #labels = paste(expression(rho ["S"]),
                    #               as.character(round(rpkm.cor.onlycomplete[i,j], digits = 4), sep="=")),
                    adj = c(0,1), cex=1.5)
          # fill if nb replicates in current condition than < nRepMax
          if (length(idces)<nRepM) {
            for (i.idces in length(idces):(nRepM-1)) {
              graphics::plot(0, col=grDevices::rgb(0,0,0,0), axes=FALSE, xlab="", ylab="")
            }
          }
          }
        }

        # graphics::par(mfrow=(c((n-1),(n-1))), las=1, mar=c(4,4,2,1))
        # for (i in 1:(n-1)) {
        #     for (j in 1:(i-1)) {
        #       if(i!=1) graphics::plot(0, col=grDevices::rgb(0,0,0,0), axes=FALSE, xlab="", ylab="")
        #       }
        #     for (j in (i+1):n) {
        #        graphics::plot(log(rpkm.all.emptyasNA[,c(i,j)], base=10), pch=20, cex=0.8, col=grDevices::rgb(0.4, 0.4, 0.4, 0.5))
        #        rhotext <- paste("$\\rho_S = ",as.character(round(rpkm.cor.onlycomplete[i,j], digits = 4)), " $", sep="")
        #        if (XP.conditions[i] == XP.conditions[j]) {
        #           resu <- c(resu, rpkm.cor.onlycomplete[i,j])
        #        }
        #        text(x = log(min(rpkm.all.emptyasNA[,i], na.rm = TRUE),base=10),
        #             y = log(max(rpkm.all.emptyasNA[,j], na.rm = TRUE),base=10),
        #             labels=TeX(rhotext),   #TeX('$\\rho_S = $'),
        #             #labels = paste(expression(rho ["S"]),
        #             #               as.character(round(rpkm.cor.onlycomplete[i,j], digits = 4), sep="=")),
        #             adj = c(0,1))
        #     }
        # }
    # same, including 0

    } else {

        graphics::par(mfrow=(c(nCond, (nRepM-1))), las=1, mar=c(4,4.5,2,0.5))
        #graphics::par(mfrow=(c((n-1),(n-1))), las=1, mar=c(4,4,2,1))
        for (i in conditions) {
          idces <- which(XP.conditions==i)
          for (i.idces in 1:(length(idces)-1)) {
            ii <- idces[i.idces]
            jj <- idces[(i.idces + 1)]
            # since 0 are included, log(x) is replaced by log(1+x)
            graphics::plot(log(1+rpkm.all.emptyasNA[,c(ii,jj)], base=10), pch=20, cex=0.8, col=grDevices::rgb(0.4, 0.4, 0.4, 0.5), cex.axis = 1.5, cex.lab = 1.5)
               rhotext <- paste("$\\rho_S = ",as.character(round(rpkm.cor.incl0[ii,jj], digits = 4)), " $", sep="")
               resu <- c(resu, rpkm.cor.incl0[ii,jj])
               graphics::text(x = log(1+min(rpkm.all.emptyasNA[,ii], na.rm = TRUE),base=10),
                    y = log(1+max(rpkm.all.emptyasNA[,jj], na.rm = TRUE),base=10),
                    labels=latex2exp::TeX(rhotext),   #TeX('$\\rho_S = $'),
                    #labels = paste(expression(rho ["S"]),
                    #               as.character(round(rpkm.cor.onlycomplete[i,j], digits = 4), sep="=")),
                    adj = c(0,1), cex=1.5)
          }
          # fill if nb replicates in current condition than < nRepMax
          if (length(idces)<nRepM) {
            for (i.idces in length(idces):(nRepM-1)) {
              graphics::plot(0, col=grDevices::rgb(0,0,0,0), axes=FALSE, xlab="", ylab="")
            }
          }
        }
    }
  grDevices::dev.off()
  } #end loop on PNG and EPS

  repl.correl.gene.test <- stats::wilcox.test(resu, correl.gene.DiamentTuller2016, alternative = "less")
  if (repl.correl.gene.test$p.value>0.1) {
    couleur <- "green"
    recommendation <- paste("Correlation between replicates is equal or higher than those observed,\n",
                                       "   in a variety of datasets (Diament Tuller 2016).\n",
                                     sep="")
  } else if (repl.correl.gene.test$p.value<0.1) {

      if (repl.correl.gene.test$p.value>0.05) {
          couleur <- "orange"
          recommendation <- paste("Correlation between replicates is a bit lower than those observed,\n",
                                             "   in a variety of datasets (Diament Tuller 2016).\n",
                                           sep="")
      } else if (repl.correl.gene.test$p.value<0.05) {
          couleur <- "red"
          recommendation <- paste("Correlation between replicates is notably lower than those observed,\n",
                                             "   in a variety of datasets (Diament Tuller 2016).\n",
                                           sep="")
      } else {
          message("something went wrong in function repl.correl.gene")
      }

  } else {
      message("something went wrong in function repl.correl.gene")
  }
  
  # Save values for later reference, and return
  repl.correl.gene.res <- list( cor=resu,
                                plot=repl.correl.gene.png,
                                value=round(stats::median(resu),digits = 4),
                                color=couleur,
                                recommendation=recommendation)
  save(repl.correl.gene.res,file=paste(pathout,"valNrec_","repl.correl.gene.res",".RData",sep=""))
  #  return(list(cor=resu,
  #              plot=repl.correl.gene.png,
  #              value=round(stats::median(resu),digits = 4),
  #              color=couleur,
  #              recommendation=recommendation))
  return(repl.correl.gene.res)              
}


#### correlation at {1;3;10;30;100;300}-codon resolution

covA <- function(i, outcov, codon.res, XP.names, pathout) {

  # Call Python to compute correlation at single- or grouped-codon resolution
  PythonInR::pyExecfile(system.file("covA.py", package="RiboVIEW", mustWork = TRUE))
  PycovA <- PythonInR::pyFunction("covA")
  #rPython::python.load( system.file("covA.py", package="RiboVIEW", mustWork = TRUE))
  inSCO <- paste(pathout, "Sample-",XP.names[i],".SCO", sep="")
  PycovA(inSCO, codon.res, outcov)
  #rPython::python.call("covA", inSCO, codon.res, outcov)
}

repl.correl.codon.singleres <- function(list.bam, refCDS, refFASTA, mini, maxi, 
  XP.names, XP.conditions, pathout, codon.res=1, versionStrip=FALSE) {

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))

  ## Read arguments and determine their number
  argl <- list.bam   #[1:3]
  n     <- length(argl)

  list.cov <- c()
  for (i in 1:n) {
    outcov <- tempfile()
    #covA(argl[[i]], refCDS, refFASTA, mini, maxi, outcov, codon.res, versionStrip)
    covA(i, outcov, codon.res, XP.names, pathout)
    list.cov <- c(list.cov, outcov)
  }

  # Prepare matrix containing all (mRNA;pos) pairs as ID, in rownames
  #   and each sample as a column containing coverage, or 0 if the
  #   corresponding (mRNA,pos) pair is absent from a specific sample.

  all.cov <- tempfile()

  unionMat(listInputFilesR = list.cov,
           outFile         = all.cov)

  codcov.all <- utils::read.table(all.cov, header=TRUE, sep='\t', row.names = 1)
  colnames(codcov.all) <- XP.names
  #head(codcov.all)

  # correlation of coverage between replicates
  codcov.cor.incl0 <- stats::cor(codcov.all, method = "spearman")
  codcov.cor.incl0

  # Nb conditions and Nb replicates per condition
  conditions <- unique(XP.conditions)
  nCond      <- length(conditions)
  repV       <- apply(cbind(conditions), 1,
                 function(k) {sum(XP.conditions==k)})
  nRepM <- max(repV)

  # plot, k*codon resolution

  repl.correl.codon.png <- paste(pathout,"Replicates-codon-res-",codon.res,".png",sep="") #tempfile()
  repl.correl.codon.eps <- paste(pathout,"Replicates-codon-res-",codon.res,".eps",sep="") 

  for (plotformat in c("png","eps")) {
    if (plotformat=="png") {
      grDevices::png(filename = repl.correl.codon.png, width = 800*(nRepM-1)/nCond, height = 800)
    } else {
      grDevices::cairo_ps(repl.correl.codon.eps,width = 10*(nRepM-1)/nCond,height = 10)
    }
  #grDevices::png(filename = repl.correl.codon.png, width = 800*(nRepM-1)/nCond, height = 800)

    graphics::par(mfrow=c(nCond, (nRepM-1)), las=1, mar=c(4,4.5,2,0.5))
    resu <- c()
    for (i in conditions) {
      idces <- which(XP.conditions==i)
      for (i.idces in 1:(length(idces)-1)) {
        ii <- idces[i.idces]
        jj <- idces[(i.idces + 1)]
        graphics::plot(log(1+codcov.all[,c(ii,jj)], base=10), pch=20, cex=0.8, col=grDevices::rgb(0.4, 0.4, 0.4, 0.5), cex.axis = 1.5, cex.lab = 1.5)
        rhotext <- paste("$\\rho_S = ",as.character(round(codcov.cor.incl0[ii,jj], digits = 4)), " $", sep="")
        resu <- c(resu, codcov.cor.incl0[ii,jj])
        #
        graphics::text(x = log(min(1+codcov.all[,ii], na.rm = TRUE),base=10),
              y = log(max(1+codcov.all[,jj], na.rm = TRUE),base=10),
              labels=latex2exp::TeX(rhotext),   #TeX('$\\rho_S = $'),
              #labels = paste(expression(rho ["S"]),
              #               as.character(round(rpkm.cor.onlycomplete[i,j], digits = 4), sep="=")),
              adj = c(0,1), cex=1.5 )
      }
      # fill if nb replicates in current condition than < nRepMax
      if (length(idces)<nRepM) {
        for (i.idces in length(idces):(nRepM-1)) {
          graphics::plot(0, col=grDevices::rgb(0,0,0,0), axes=FALSE, xlab="", ylab="")
        }
      }

    }

  # n.fen <- n - 1
  # graphics::par(mfrow=c(n.fen, n.fen), las=1, mar=c(4,4,1,1))
  # resu <- c()
  # for (i in 1:(n-1)) {
  #     for (j in 1:(i-1)) {
  #       if(i!=1) graphics::plot(0, col=grDevices::rgb(0,0,0,0), axes=FALSE, xlab="", ylab="")
  #       }
  #     for (j in (i+1):n) {
  #        graphics::plot(log(1+codcov.all[,c(i,j)], base=10), pch=20, cex=0.8, col=grDevices::rgb(0.4, 0.4, 0.4, 0.5))
  #        rhotext <- paste("$\\rho_S = ",as.character(round(codcov.cor.incl0[i,j], digits = 4)), " $", sep="")
  #        #
  #        if (XP.conditions[i] == XP.conditions[j]) {
  #           resu <- c(resu, codcov.cor.incl0[i,j])
  #        }
  #        #
  #        text(x = log(min(1+codcov.all[,i], na.rm = TRUE),base=10),
  #             y = log(max(1+codcov.all[,j], na.rm = TRUE),base=10),
  #             labels=TeX(rhotext),   #TeX('$\\rho_S = $'),
  #             #labels = paste(expression(rho ["S"]),
  #             #               as.character(round(rpkm.cor.onlycomplete[i,j], digits = 4), sep="=")),
  #             adj = c(0,1) )
  #     }
  # }

  grDevices::dev.off()
  } #end loop on PNG and EPS

  if (min(resu)>=0.60) {
    couleur <- "green"
  } else if (min(resu)<0.60) {

      if (min(resu)>=0.40) {
          couleur <- "orange"
      } else if (min(resu)<0.40) {
          couleur <- "red"
      } else {
          message("something went wrong in function repl.correl.codon")
      }

  } else {
      message("something went wrong in function repl.correl.codon")
  }

  
  # Save values for later reference, and return
  repl.correl.codon.singleres <- list(cor=resu,
                                      plot=repl.correl.codon.png,
                                      value=round(stats::median(resu),digits = 4),
                                      color=couleur)
  save(repl.correl.codon.singleres, 
       file=paste(pathout, "valNrec_", "repl.correl.codon.singleres-", codon.res,".RData", sep=""))
  #return(list(cor=resu,# à mettre à jour
  #            plot=repl.correl.codon.png,
  #            value=round(stats::median(resu),digits = 4),
  #            color=couleur))  return(repl.correl.codon.res)     
  return(repl.correl.codon.singleres)         
}


# repl.correl.codon : Correlation of footprint coverage at codon resolution
#
#
#
# This function takes a list of samples and experimental conditions in input,
#   calculates correlations at codon level between replicates.
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
#        XP.names       Vector of names for each sample
#        XP.conditions  Vector of experimental conditions for each sample
#        pathout        Address where output files will be written
#
# Out:
#       A list containing:
#          - cor                 Correlation between replicates
#          - plot                Address of plot file in png format
#          - value               Median of correlation across different conditions
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#
repl.correl.codon <- function(list.bam, refCDS, refFASTA, mini, maxi, XP.names, XP.conditions, pathout) {
  
  # init
  codon.res.orange <- ">100"
  codon.res.vert   <- ">100"
  
  #for (codon.res in c(5,10,20)) {
  for (codon.res in c(3,10,30,100)) {
    
    repl.correl.codon.res <- repl.correl.codon.singleres(list.bam, refCDS, refFASTA, mini, maxi, XP.names, XP.conditions, pathout, codon.res)
    
    # Update indicators (corr Spearman > cutoff) if some codon resolution satisfies them
    if (repl.correl.codon.res$color == "orange" & codon.res.orange == ">100") {
        codon.res.orange <- codon.res
    } else if (repl.correl.codon.res$color == "green"  & codon.res.vert == ">100") {
        codon.res.vert <- codon.res
        if (codon.res.orange == ">100") {codon.res.orange <- codon.res}
        break
    }
  }
  
  repl.correl.codon.res$recommendation <- paste("Spearman correlation >0.60 (resp. 0.40) between replicates requires averaging over ",
                                                as.character(codon.res.vert),
                                                " (resp. ",
                                                as.character(codon.res.orange),
                                                ") codons. Tested codon resolutions are : 3,10,30,100 (More info : fig3c,4b, Diament and Tuller 2016).\n",
                                                sep="")
  # Save values for later reference, and return
  save(repl.correl.codon.res, 
       file=paste(pathout, "valNrec_", "repl.correl.codon.res",".RData", sep=""))
  return(repl.correl.codon.res)

}

#### Similarity analysis between replicates, data
# data : codon occupancy, normalized
BCO.norm <- function(fileBCO) {
    BCO <- utils::read.table(fileBCO)
    p1 <- BCO/matrix(rep(colSums(BCO),64),ncol=ncol(BCO),byrow=TRUE)
    p1 <- p1[!(rownames(p1) == "uaa" | rownames(p1) == "uag" | rownames(p1) == "uga"),]
    n1 <- 3 * p1$A / (p1$up3 + p1$up2 + p1$up1)
    names(n1) <- rownames(p1)
    return(n1)
}

# repl.correl.heatmap : Adequation of replicates using a heatmap with hierarchical clustering
#
#
#
# This function takes a list of samples and experimental conditions in input,
#   generates a heatmap with hierarchical clustering, and compares this clustering
#   with the groups of replicates.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#
# Out:
#       A list containing:
#          - plot                Address of plot file in png format
#          - value               Spearman correlation (in absolute value) between clusters and groups of replicates
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#
repl.correl.heatmap <- function(XP.conditions, XP.names, pathout) {

    argl.pre <- XP.conditions   #[1:3]
    n     <- length(argl.pre)

    # Colors
    col.base <- c('bisque', 'darkolivegreen1', 'lightpink1', 'cadetblue1',
                  'azure1','gold1','indianred1','darkslategray1',
                  'gainsboro','darkseagreen1','plum1','lavender')
    # set of conditions and nb of conditions
    set.cond     <- unique(XP.conditions)
    nb.cond     <- length(set.cond)
    col.repl <- rep(NA, length(XP.conditions))
    for (l in 1:nb.cond) {
      cond <- set.cond[l]
      col.repl[XP.conditions==cond] <- col.base[l]
    }
    col.codons  <- rep('grey', 61)

    # Construction matrix containing codons in rows, samples in columns
    i <- 1
    outBCOi <- paste(pathout, "Sample-",XP.names[i],"_BCO", sep="")
    mrep <- cbind(BCO.norm(outBCOi))
    for (i in 2:n) {
        outBCOi <- paste(pathout, "Sample-",XP.names[i],"_BCO", sep="")
        mrep <- cbind(mrep, BCO.norm(outBCOi))
    }

    colnames(mrep) <- XP.names

    # Plot
    repl.correl.heatmap.png <- paste(pathout,"Replicates-heatmap.png",sep="") #tempfile()
    repl.correl.heatmap.eps <- paste(pathout,"Replicates-heatmap.eps",sep="") #tempfile()

    # plot PNG and EPS
    for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = repl.correl.heatmap.png, width = 800, height = 800)
      } else {
        grDevices::cairo_ps(repl.correl.heatmap.eps,width = 10,height = 10)
      }
      #grDevices::png(filename = repl.correl.heatmap.png, width = 800, height = 800)
          rep.heat <- gplots::heatmap.2(as.matrix(mrep),
                      RowSideColors = col.codons,
                      ColSideColors = col.repl)
          # exploit hierarchical clustering of samples
          rd <- rep.heat$colDendrogram #str(rd)
          nb.exp <- length(unique(XP.conditions))
          clusters <- dendextend::cutree(rd, k = nb.exp) #from dendextend package
          # correlation between expected and observed
          rs <- abs(stats::cor(x = XP.conditions, y = clusters, method="spearman"))

      grDevices::dev.off() #for either PNG or EPS
    }

    # List of outputs
    repl.correl.heatmap.res <- c()
    repl.correl.heatmap.res$plot <- repl.correl.heatmap.png
    repl.correl.heatmap.res$value <- rs
    repl.correl.heatmap.res$color <- ifelse(rs>0.9, yes = "green", no = ifelse(rs>0.7, yes = "orange", no = "red"))
    repl.correl.heatmap.res$recommendation <- ifelse(rs>=0.8,
                                                 yes = "Replicates cluster well together (clusters members correlate with defined replicates - spearman rho > 0.8)",
                                                 no = ifelse(rs>=0.5,
                                                             yes = "Replicates cluster somewhat together (clusters members correlate moderately with defined replicates - spearman rho>0.5)",
                                                             no = "Replicates cluster poorly together (clusters members differ from replicates order - spearman rho < 0.5)"))
    # Save values for later reference, and return
    save(repl.correl.heatmap.res, 
         file=paste(pathout, "valNrec_", "repl.correl.heatmap.res",".RData", sep=""))
    return(repl.correl.heatmap.res)
}

