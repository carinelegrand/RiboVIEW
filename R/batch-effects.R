# Batch effects
#
#
# Contents :
#              batch.effects.lm
#              batch.effects.pca
#              

# batch.effects.lm.e : Apply a linear fit to codon enrichment with a, c, g or u 
#                      as explanatory variable.
#
#
# This function takes a list of samples and their conditions as input and applies
#   a linear fit to each of these samples' codon occupancies. Explanatory variable
#   is successively a, c, g and u. A plot of coefficients, p-value and a 
#   description are returned
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#        nvert          Number of legend items which can be displayed corectly, 12 by default
#
# Out:
#       A list containing:
#          - plot                Address of png plot file
#          - value               Minimum p-value obtained from the fits
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#         -  recommendation.mat  Matrix of coefficients for each fit
#

batch.effects.lm.e <- function(XP.conditions, XP.names, pathout, nvert=12) {

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))
  
  # parameter nvert, number of legend items that can be displayed corectly :
  #              -> 12 works well on linux, to adapt wherever necessary

  # Number of arguments, reconstruct listSamples
  n <- length(XP.conditions)

  listSamples <- c()
  for (i in 1:n) {
    listSamples[[i]] <- paste(pathout, "Sample-",XP.names[i],"_Codon-enrichment-unbiased_mean", sep="")
  }
  argl <- listSamples

  # Parameter file
  #NOK : acgu <- system.file("acgu-content.tsv", package="RiboVIEW", mustWork = TRUE)
  acgu <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                "/extdata/acgu-content.tsv", sep="")
  acgu.content <- utils::read.table(acgu) #"inst/extdata/acgu-content.tsv"

  # Dependence on nucleobases : coeff, standard error and raw p-value per sample and nucleobase
  for (i in 1:n) {
    Enrich.pre1 <- utils::read.table(argl[[i]])
    Enrich.pre2 <- Enrich.pre1[,'X0']
    names(Enrich.pre2) <- rownames(Enrich.pre1)
    #BCOcurr.pre <- BCO.norm(argl[[i]])
    #
    #if (all(names(Enrich.pre2) == row.names(acgu.content))) {  Enrich <- Enrich.pre2
    #} else                                                  {  
    Enrich <- Enrich.pre2[row.names(acgu.content)] #}
    #
    marg.fit <- unlist(exploitfitrob(MASS::rlm(Enrich ~ a , data=data.frame(Enrich, acgu.content))))
    marg.fit <- rbind( marg.fit , unlist(exploitfitrob(MASS::rlm(Enrich ~ c , data=data.frame(Enrich, acgu.content)))) )
    marg.fit <- rbind( marg.fit , unlist(exploitfitrob(MASS::rlm(Enrich ~ g , data=data.frame(Enrich, acgu.content)))) )
    marg.fit <- rbind( marg.fit , unlist(exploitfitrob(MASS::rlm(Enrich ~ u , data=data.frame(Enrich, acgu.content)))) )
    row.names(marg.fit) <- c("a","c","g","u")
    #
    if (i==1) {marg.fit.mat <- marg.fit
    } else    {marg.fit.mat <- cbind(marg.fit.mat, marg.fit) }
    ech <- XP.names[i]
    colnames(marg.fit.mat)[(1+(i-1)*3)] <- paste(ech, ", Coeff", sep="")
  }
  #marg.fit.mat

  # plot preparation : filename and elements
  hauteurs <- t(marg.fit.mat[,seq(1,n*3,3)])
  hauteurs.se <- t(marg.fit.mat[,seq(2,n*3,3)])
  ylim <- c(-1,1)*1.2*round(max(abs(hauteurs)+abs(hauteurs.se)),1)

  # Plot preparation : colors, recycled if necessary
  couleurs.pre <- c("lightblue", "mistyrose", "lavender", "lightcyan")
  if (n<5) {couleurs <- couleurs.pre[1:n]
  } else   {couleurs <- rep(couleurs.pre, ceiling(n/4))[1:n] }

  # plot preparation : color for border and error bars
  col.border <- "gray20"

  # plot
  batcheffects.plot.png <- paste(pathout,"Batch-effects-lm-e.png",sep="") #tempfile()
  batcheffects.plot.eps <- paste(pathout,"Batch-effects-lm-e.eps",sep="")

  for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = batcheffects.plot.png, width = 800, height = 800)
      } else {
        grDevices::cairo_ps(batcheffects.plot.eps,width = 10,height = 10)
      }
    #grDevices::png(filename = batcheffects.plot.png, width = 800, height = 800)

    # Split plotting device in 2 ; store index of created device
    #dev.new()
    graphics::layout(cbind(c(1,2)), heights = c(6,3.5))
    #devlist <- dev.list()
    #layout.dev <- devlist[length(devlist)]

    # Plot margins
    graphics::par(las=1, mar=c(3, 4.1, 2.1, 2.1))
    # 1st panel : barplot
    mp <- graphics::barplot(hauteurs, beside = TRUE,
            col = couleurs, border = col.border, ylab = "Coefficient (-)",
            #legend= ..., args.legend = list(bty = "n", x = "top"),
            ylim = ylim,
            main = "")
    graphics::segments(mp, hauteurs - hauteurs.se, mp, hauteurs + hauteurs.se, col = col.border, lwd = 1.5)
    graphics::axis(1, labels = NA, at = colMeans(mp))
    # 2nd panel : legend
    graphics::par(mar=c(0,0,0,0))
    graphics::plot(NA, axes=FALSE, xlab="", ylab="", main="", xlim=c(0,1), ylim=c(0,1))
    if (n <= nvert ) {
      graphics::legend(x = "center",legend = apply(cbind(XP.names), 1,
                           function(x) {nom <- strsplit(x, split = ",")[[1]][1] } ),
             pch=22, pt.cex = 1.5, pt.bg=couleurs, col=col.border, bty='n')
    } else {
      legend.text <- apply(cbind(XP.names), 1,
                            function(x) {nom <- strsplit(x, split = ",")[[1]][1] } )
      legend.col  <- couleurs
      legend.x    <- seq(0+1/ceiling(n/nvert)/2, 1, 1/ceiling(n/nvert))
      for (ii in 1:ceiling(n/nvert)) {
        mi <- 1 + (ii-1)*nvert
        ma <- min(ii*nvert, n)
        graphics::legend(x = legend.x[ii], y = 0.5, xjust = 0.5, yjust = 0.5,
               text.col = "grey30",
               legend = legend.text[mi:ma],
               pch=22, pt.cex = 1.5, pt.bg=legend.col[mi:ma], col=col.border, bty='n')
      }
    }
  grDevices::dev.off()
  } #end loop PNG or EPS

  # Close dev opened for layout
  #grDevices::dev.off(as.numeric(layout.dev))

  ## Restore margins as initial
  #graphics::par(mar=c(5.1, 4.1, 4.1, 2.1))

  # Value
  pvalues <- marg.fit.mat[,seq(3,n*3,3)]
  padj <- rbind(stats::p.adjust(pvalues["a",], method = "BH"),
            stats::p.adjust(pvalues["c",], method = "BH"),
            stats::p.adjust(pvalues["g",], method = "BH"),
            stats::p.adjust(pvalues["u",], method = "BH"))
  row.names(padj) <- row.names(pvalues)
  Value = min(padj)

  # Color
  Color <- ifelse(Value < 0.05, yes = "red", no = ifelse(Value<0.1, yes="orange", no="green"))

  # Recommendation ; use cat(Recommendation)
  Recommendation <- ifelse(Value < 0.05,
                           yes = paste("Codon occupancy depends significantly (cutoff 0.05) on at least 1 nucleobase.\n",
                                       "-> When dependence is similar over replicates this suggests an\n",
                                       "   actual preference of specific nucleobases in the corresponding\n",
                                       "   experimental condition.\n",
                                       "-> If dependence happens in some replicates only,\n",
                                       "   then this points to a batch effect.\n",
                                       "   We recommend adjusting the affected replicates.\n",
                                       "   Coefficients, standard error and raw p-value are as follows : ",
                                     sep=""),
                           no = ifelse(Value<0.1,
                                       yes=paste("Codon occupancy might depend on at least 1 nucleobase.\n",
                                       "-> When dependence is similar over replicates this suggests an\n",
                                       "   actual result triggered by experimental conditions.\n",
                                       "-> If dependence happens in some replicates and not others,\n",
                                       "   then this points to a batch effect.\n",
                                       "   We recommend adjusting the affected replicates.\n",
                                       "   Coefficients, standard error and raw p-value are as follows : ",
                                     sep=""),
                                       no=paste("Codon occupancy seems independent of individual nucleobases.",
                                     sep="")))

  # Result
  batch.effects.lm.e.res <- list(plot = batcheffects.plot.png, 
                               value = Value, 
                               color = Color, 
                               recommendation = Recommendation, 
                               recommendation.mat = marg.fit.mat)
  # Save values for later reference, and return
  save(batch.effects.lm.e.res, 
       file=paste(pathout, "valNrec_", "batch.effects.lm.e.res",".RData", sep=""))
  return(batch.effects.lm.e.res)
}

# batch.effects.pca : Apply a PCA and tSNE to visualize batch effects or groups 
#                       of samples by condition.
#
#
# This function takes a list of samples and their conditions as input and applies
#   a PCA dimension reduction and a tSNE visualisation. A plot of principal axes, 
#   visualisation axes, a p-value for PCA principal axes, corresponding color and 
#   text are returned.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#
# Out:
#       A list containing:
#          - plot                Address of png plot file combining PCA and tSNE
#          - plot.pca            Address of png plot file for PCA principal axes
#          - plot.tsne           Address of png plot file for tSNE
#          - value               p-value for PCA principal axes
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#
batch.effects.pca <- function(XP.conditions, XP.names, pathout) {

  # Save and restore graphical parameters if unexpected exit 
  par.prev <- suppressGraphics(par(no.readonly = TRUE))
  on.exit(suppressGraphics(par(par.prev, new=FALSE)))

  # Number of arguments, reconstruct listSamples
  n <- length(XP.conditions)
  
  # Approx nb of replicates
  ncond <- length(unique(XP.conditions))
  nrep <- n / ncond

  listSamples <- c()
  for (i in 1:n) {
    listSamples[[i]] <- paste(pathout, "Sample-",XP.names[i],"_BCO", sep="")
  }
  argl <- listSamples

  # Parameters
  #NOK : colf <- system.file("codonColor.tsv"; package="RiboVIEW", mustWork = TRUE)
  colf <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                "/extdata/codonColor.tsv", sep="")
  codonColors <- utils::read.table(colf) #"inst/extdata/codonColor.tsv"
  colSamples <- RColorBrewer::brewer.pal(9, "Set1")
  colSamples[6] <- "gold"


  # Parameter file
  acgu <- paste(system.file(package="RiboVIEW", mustWork = TRUE), 
                "/extdata/acgu-content.tsv", sep="")
  acgu.content <- utils::read.table(acgu) #"inst/extdata/acgu-content.tsv"

  # Reformat
  for (i in 1:n) {
    BCOcurr.pre <- BCO.norm(argl[[i]])
    #
    if (all(names(BCOcurr.pre) == row.names(codonColors))) {  BCOcurr <- BCOcurr.pre
    } else                                                      {  BCOcurr <- BCOcurr.pre[row.names(acgu.content)] }
    #
    if (i==1) {pourPca <- cbind(BCOcurr)
    } else    {pourPca <- cbind(pourPca, cbind(BCOcurr)) }
    str.spl   <- strsplit(argl[[i]], split = "/")[[1]]
    str.spl2  <- str.spl[length(str.spl)]
    ech       <- strsplit(str.spl2, split = ".", fixed=TRUE)[[1]][1]
    colnames(pourPca)[i] <- ech
  }
  #head(pourPca)

  # Pca plot
  # Pca plot : preparation
  # Here, observations (occupancy on 61 codons, excluding stop codons, already normalised)
  # are supposed directly comparable, which is why no further scaling is performed.
  pca.bco.pre <- svd(stats::cov(pourPca))  #svd(cor(pourPca))
  pca.bco <- pca.bco.pre$u
  pca.bco.var <- sqrt(pca.bco.pre$d)
  pca.bco.pc.var <- 100 * pca.bco.var / sum(pca.bco.var)
  # Pca plot : add fake points (unplotted) to control xlim and ylim and pass variance to pairs
  mini <- apply(pca.bco, 2, min)
  maxi <- apply(pca.bco, 2, max)
  #
  fake.min      <- mini - abs(maxi-mini)*0.1
  fake.max      <- maxi + abs(maxi-mini)*0.1
  hidden.pc.var <- fake.min + pca.bco.pc.var * (fake.max - fake.min)/100 #fake.min + 0.001*pca.bco.pre$d
  #
  pca.bco.augm <- rbind(pca.bco,fake.min, fake.max, hidden.pc.var ) #fake.min + 0.001*pca.bco.pre$d)
  # Pca plot : upper panel and diagonal panel custom functions
  panel.up <- function(x, y, col = graphics::par("col"), bg = NA, pch = graphics::par("pch"), cex = 1) {
      # remove fake min and max, only there to control xlim and ylim
      x2 <- x[1:(length(x)-3)]
      y2 <- y[1:(length(y)-3)]
      graphics::text(x2, y2, labels = 1:ncol(pourPca), col = col, cex=cex, font = 2)
  }
  #
  panel.var <- function(x, panel.var.sum, ...) {
    usr <- graphics::par("usr"); on.exit(graphics::par(usr))
    graphics::par(usr = c(usr[1:2], 0, 1.5) )
    fake.min.int       <- x[(length(x)-2)]
    fake.max.int       <- x[(length(x)-1)]
    hidden.pc.var.int  <- x[length(x)]
    pca.bco.pc.var.int <- 100* (hidden.pc.var.int - fake.min.int) / (fake.max.int - fake.min.int)
    leg.text <- sprintf("%.1f",pca.bco.pc.var.int)
    #print(leg.text)
    graphics::legend("center", legend = paste("%var=",leg.text,sep=""), bty='n', cex=1)
  }
  # Pca plot : plot and write to tmp file
  batcheffects.pca.plot.png <-  paste(pathout,"Batch-effects-prelim-pca.png",sep="") #tempfile()
  batcheffects.pca.plot.eps <-  paste(pathout,"Batch-effects-prelim-pca.eps",sep="") #tempfile()

  for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = batcheffects.pca.plot.png, width = 800, height = 800)
      } else {
        grDevices::cairo_ps(batcheffects.pca.plot.eps,width = 10,height = 10)
      }
    #grDevices::png(filename = batcheffects.pca.plot.png, width = 800, height = 800)
     graphics::pairs(pca.bco.augm, lower.panel = NULL, col=c(colSamples, 2.88863), pch=20, cex=1.,
        upper.panel = panel.up, diag.panel = panel.var, main="PCA axes 1 to n",
        labels = apply(cbind(1:n), 1, function(x) {paste("PC",x,sep="")}) )
  grDevices::dev.off()
  } # end loop PNG or EPS

  # tSNE plot
  # tSNE plot : preparation
  # Nota : perplexity is the expected nb of points per group, therefore the nb of
  #   replicates is used
  suppressMessages(tsne.bco <- tsne::tsne(t(pourPca), perplexity=nrep))
  tmi <- apply(tsne.bco, 2, min)
  tma <- apply(tsne.bco, 2, max)
  tmi2 <- tmi-0.1*(tma-tmi)
  tma2 <- tma+0.1*(tma-tmi)
  # tSNE plot : plot and write to file
  batcheffects.tsne.plot.png <-  paste(pathout,"Batch-effects-prelim-tsne.png",sep="") #tempfile()
  batcheffects.tsne.plot.eps <-  paste(pathout,"Batch-effects-prelim-tsne.eps",sep="") #tempfile()

  for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = batcheffects.tsne.plot.png, width = 260, height = 260)
      } else {
        grDevices::cairo_ps(batcheffects.tsne.plot.eps,width = 3.25,height = 3.25)
      }
    #grDevices::png(filename = batcheffects.tsne.plot.png, width = 260, height = 260)
    graphics::par(mar=c(3,3,1,0.5), las=1)
    graphics::plot(tsne.bco, col=colSamples, pch=20, type='n',
         xlab="", ylab="", #, xlab="tSNE_1", ylab="tSNE_2",
         cex.axis =0.8, cex.lab = 0.8,
         xlim=c(tmi2[1],tma2[1]), ylim=c(tmi2[2],tma2[2]),
         main="tSNE", cex.main=0.9)
    graphics::text(tsne.bco, col=colSamples, labels = 1:n, font=2)
    graphics::mtext(1, text = "tSNE-1", line = 2, cex=0.8)
    graphics::par(las=3)
    graphics::mtext(2, text = "tSNE-2", line = 2, cex=0.8)
    graphics::par(las=1)
  grDevices::dev.off()
  } # end loop EPS or PNG

  pca.plot.png <- paste(pathout,"Batch-effects-pca.png",sep="") #tempfile()
  #pca.plot.eps <- paste(pathout,"Batch-effects-pca.eps",sep="") #tempfile()

  # Assemble PCA and tSNE plot in one png image.
  #composite -gravity southwest -compose atop tsne.plot pca.plot batcheffects.pca.plot.png
  system(paste("composite -gravity southwest -compose atop  ",
               batcheffects.tsne.plot.png, " ", batcheffects.pca.plot.png, " ", pca.plot.png,
               sep=""))
  ### p-value from Resampling / Bootstrap
  permut.res <- permut.pca.mat(pourPca, pca.bco.var)
  Value = permut.res[3,]
  
  # Color
  Color <- ifelse(min(Value) < 0.01, yes = "red", no = ifelse(min(Value) < 0.05, yes="orange", no="green"))
  # old : for TW : Color <- ifelse(Value > 2.0236, yes = "red", no = ifelse(Value>0.9794, yes="orange", no="green"))

  # Recommendation ; use cat(Recommendation)
  Recommendation <- ifelse(min(Value) < 0.01,
                           yes = paste("At least one eigenvalue is significant at 0.01 level. This corresponds either\n",
                                       "   to batch effects or to biological variations.\n",
                                       "   If batches (for example samples sequenced on the same lane)\n",
                                       "   seem to cluster together on PCA or tSNE plots, you should\n",
                                       "   probably correct for batch effects using principal components\n",
                                       "   Full list of eigenvalues from PC1 to PCn (in % of variance explained) : \n",
                                       "      ",Value,".\n",
                                       "Further information on batch effects in ribosome profiling experiments :\n",
                                       "    Flis et al. 2018, Liu et al. 2019.\n",
                                     sep=""),
                           no = ifelse(min(Value) < 0.05,
                                       yes=paste("At least one eigenvalue is significant at 0.05 level. This might\n",
                                       "   suggest either batch effects or biological variation.\n",
                                       "   If batches (for example samples sequenced on the same lane)\n",
                                       "   seem to cluster together on PCA or tSNE plots, you should\n",
                                       "   probably correct for batch effects using principal components\n",
                                       "   Full list of eigenvalues from PC1 to PCn (in % of variance explained) : \n",
                                       "      ",Value,".\n",
                                       "Further information on batch effects in ribosome profiling experiments :\n",
                                       "    Flis et al. 2018, Liu et al. 2019.\n",
                                     sep=""),
                                       no=paste("No batch effect or other structure could be detected (at 0.05\n",
                                                "  significance level).\n",
                                     sep="")))
  # old : tracy widom
  #  Recommendation <- ifelse(Value > 2.0236,
  #                           yes = paste("First eigenvalue is significant at 0.01 level. This corresponds either\n",
  #                                       "   to batch effects or to biological variations.\n",
  #                                       "   If batches (for example samples sequenced on the same lane)\n",
  #                                       "   seem to cluster together on PCA or tSNE plots, you should\n",
  #                                       "   probably correct for batch effects using principal components\n",
  #                                       "   Note : tSNE is useful for visualization but isn't appropriate\n",
  #                                       "          for adjusting against a batch effect.\n",
  #                                     sep=""),
  #                           no = ifelse(Value>0.9794,
  #                                       yes=paste("First eigenvalue is significant at 0.05 level. This might\n",
  #                                       "   suggest either batch effects or biological variation.\n",
  #                                       "   If batches (for example samples sequenced on the same lane)\n",
  #                                       "   seem to cluster together on PCA or tSNE plots, you should\n",
  #                                       "   probably correct for batch effects using principal components\n",
  #                                       "   Note : tSNE is useful for visualization but isn't appropriate\n",
  #                                       "   for adjusting against a batch effect.\n",
  #                                     sep=""),
  #                                       no=paste("No batch effect or other structure could be detected (at 0.05\n",
  #                                                "  significance level).\n",
  #                                     sep="")))

  # Result
  batch.effects.pca.res <- list(plot = pca.plot.png, 
                                plot.pca = batcheffects.pca.plot.png,
                                plot.tsne = batcheffects.tsne.plot.png, 
                                value = Value, 
                                color = Color, 
                                recommendation = Recommendation)
  # Save values for later reference, and return
  save(batch.effects.pca.res, 
       file=paste(pathout, "valNrec_", "batch.effects.pca.res",".RData", sep=""))
  return(batch.effects.pca.res)
  }

