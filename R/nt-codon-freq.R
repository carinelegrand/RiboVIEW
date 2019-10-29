# Nucleotide or codon frequencies at footprint boundaries
#
# Content :
#           ntcodon.freq.nt
#           ntcodon.freq.ntsingle
#           map.triNt.2.Codon
#           ntcodon.freq.cod
#           ntcodon.freq.cod.single
#

# ntcodon.freq.nt : Calculate nt counts at footprint boundaries and generates logoplot
#
#
#
# This function takes a list of samples and experimental conditions in input and
#   generates nucleotide counts at footprint boundaries and a logoplot for each sample.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#
# Out:
#       A list containing:
#          - plot                Address of png logoplot file
#          - value               Bit value correspondign to the logoplot
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#
ntcodon.freq.nt <- function(XP.conditions, XP.names, pathout) {

    n <- length(XP.conditions)

    ntcodon.freq.nt.res.i <- c()

    for (i in 1:n) {
        BCOi <- paste(pathout, "Sample-",XP.names[i], sep="")
        ntcodon.freq.nt.res.i[[i]] <- ntcodon.freq.nt.single(BCOi, i, XP.names, pathout)
    }

    ### Summarize over samples
    ntcodon.freq.nt.res <- c()
    ##  Plot
    ntcodon.freq.nt.Plot.summary.pre <- paste(pathout,"Nt-Freq-pre.png",sep="")
    ntcodon.freq.nt.Plot.summary <- paste(pathout,"Nt-Freq.png",sep="")
    # Assemble all samples
    system(paste("montage -mode concatenate -tile 1x",n," ",
             pathout, "Nt-Freq-*.png ",
             ntcodon.freq.nt.Plot.summary.pre,
             sep=""))
    # Fit image in 800x800 limits
    system(paste("convert ",
                 ntcodon.freq.nt.Plot.summary.pre,
                 " -resize 800x800 ",
                 ntcodon.freq.nt.Plot.summary, sep=""))
    # Remove obsolete files
    system(paste("rm ", pathout, "Nt-Freq-*.png", sep=""))
    #
    ntcodon.freq.nt.res$plot <- ntcodon.freq.nt.Plot.summary
    #
    ## Value
    # Each sample's bit value
    val.summary <- sapply(ntcodon.freq.nt.res.i, function(x) x$value)
    # Return all
    ntcodon.freq.nt.res$value <- val.summary
    #
    ## Color
    # Each sample's color
    col.summary <- sapply(ntcodon.freq.nt.res.i, function(x) x$color)
    # Reduce to highest color level as soon as it is found
    if (sum(col.summary=="red") >= 1)  {col.summary = "red"}
    if (sum(col.summary=="orange") >= 1) {col.summary = "orange"}
    if (sum(col.summary=="green") >= 1)   {col.summary = "green"}
    #
    ntcodon.freq.nt.res$color <- col.summary
    #
    ## Recommendation
    # Each sample's recommendation
    rec.summary.0 <- sapply(ntcodon.freq.nt.res.i, function(x) x$recommendation)
    rec.summary.1 <- paste(XP.names, rec.summary.0, sep = " : ")
    rec.summary.2 <- paste(rec.summary.1, collapse=" \n")
    #
    ntcodon.freq.nt.res$recommendation <- rec.summary.2
    #
    ## Save values for later reference, and return
    save(ntcodon.freq.nt.res, 
         file=paste(pathout, "valNrec_", "ntcodon.freq.nt.res",".RData", sep=""))
    return(ntcodon.freq.nt.res)
}

ntcodon.freq.nt.single <- function(BCO, i, XP.names, pathout) {

    limNt <- utils::read.table(paste(BCO,".limNt",sep=""), header=TRUE, row.names = 1)

    ntcodon.freq.nt.Plot.png <- paste(pathout,"Nt-Freq-",i,".png",sep="") #tempfile()
    ntcodon.freq.nt.Plot.eps <- paste(pathout,"Nt-Freq-",i,".eps",sep="") #tempfile()

    nr.Nt <- nrow(limNt)
    nc.Nt <- ncol(limNt)

    pwm.pre <- as.matrix(limNt) / matrix(colSums(limNt), nrow = nr.Nt, ncol=nc.Nt, byrow = TRUE)
    pwm.5 <- pwm.pre[,1:(nc.Nt/2)]
    pwm.3.rev <- pwm.pre[,rev((nc.Nt/2 + 1):nc.Nt)] #pwm.3 <- pwm.pre[,(nc.Nt/2 + 1):nc.Nt]

    N <- 4
    nseqs <- min(colSums(limNt))
    seq_type <- "RNA"

    p.5 <- CLggseqlogo_pfm( pwm.5, N, nseqs, seq_type, method = 'bits', 
                            x_lab = paste("nt from 5', ", XP.names[i], sep="") )
    infoContent.5 <- CLcomputeBits(pwm.5, N=N, Nseqs=nseqs)
    p.3 <- CLggseqlogo_pfm( pwm.3.rev, N, nseqs, seq_type, method = 'bits', 
                            x_lab = paste("nt from 3', ", XP.names[i], sep=""), )
    infoContent.3 <- CLcomputeBits(pwm.3.rev, N=N, Nseqs=nseqs)

    for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = ntcodon.freq.nt.Plot.png,width = 800,height = 100)
      } else {
        grDevices::cairo_ps(ntcodon.freq.nt.Plot.eps,width = 10,height = 1.25)
      }
    #grDevices::png(filename = ntcodon.freq.nt.Plot,width = 800,height = 100)
      gridExtra::grid.arrange(p.5, p.3, nrow=1)
    grDevices::dev.off()
    } #end loop PNG or EPS

    ntcodon.freq.nt.Value <- max(c(infoContent.5, infoContent.3))
    ntcodon.freq.nt.Color <- ifelse(ntcodon.freq.nt.Value < 0.2, yes = "green",
                                             no = ifelse(ntcodon.freq.nt.Value < 0.4, yes = "orange", no = "red"))
    ntcodon.freq.nt.Recommendation <- ifelse(ntcodon.freq.nt.Value < 0.2,
                                              yes = "No position has any overrepresented nucleotide",
                                              no  = ifelse(ntcodon.freq.nt.Value < 0.4,
                                                           yes = paste("Position(s) ",
                                                              list(
                                                              apply(cbind(which( c(infoContent.5, infoContent.3) >= 0.2 )),1, function(x) {
                                                                xx <- x
                                                                if (x>(nc.Nt/2)) { xx <- (x-nc.Nt-1) }
                                                                return(xx)
                                                                })
                                                               ),
                                                              " might have overrepresented nucleotides.",
                                                              sep=""),
                                                           no =paste("Position(s) ",
                                                              list(
                                                              apply(cbind(which( c(infoContent.5, infoContent.3) >= 0.2 )),1, function(x) {
                                                                xx <- x
                                                                if (x>(nc.Nt/2)) { xx <- (x-nc.Nt-1) }
                                                                return(xx)
                                                                })
                                                               ),
                                                              " have overrepresented nucleotides ; these shouldn't be used for normalisation.\n",
                                                              " Note that stringent adaptor trimming might systematically erase one terminal base.\n",
                                                              sep="")
                                                ))
    #
    ntcodon.freq.nt.singleres <- list(plot=ntcodon.freq.nt.Plot.png,
                                      value=ntcodon.freq.nt.Value,
                                      color=ntcodon.freq.nt.Color,
                                      recommendation=ntcodon.freq.nt.Recommendation)
    ## Save values for later reference, and return
    save(ntcodon.freq.nt.singleres, 
         file=paste(pathout, "valNrec_", "ntcodon.freq.nt.singleres-",i,".RData", sep=""))
    return(ntcodon.freq.nt.singleres)
 }


#
## Map triplets to codons
#

map.triNt.2.Codon <- function(pwm.TriNt) {
  #
  pwm.Cod <- data.frame(A = cbind(unlist(pwm.TriNt["gcu",] + pwm.TriNt["gcc",] + pwm.TriNt["gca",] + pwm.TriNt["gcg",])),
  C = cbind(unlist(pwm.TriNt["ugu",] + pwm.TriNt["ugc",])),
  D = cbind(unlist(pwm.TriNt["gau",] + pwm.TriNt["gac",])),
  E = cbind(unlist(pwm.TriNt["gaa",] + pwm.TriNt["gag",])),
  F = cbind(unlist(pwm.TriNt["uuu",] + pwm.TriNt["uuc",])),
  G = cbind(unlist(pwm.TriNt["ggu",] + pwm.TriNt["ggc",] + pwm.TriNt["gga",] + pwm.TriNt["ggg",])),
  H = cbind(unlist(pwm.TriNt["cau",] + pwm.TriNt["cac",])),
  I = cbind(unlist(pwm.TriNt["auu",] + pwm.TriNt["auc",] + pwm.TriNt["aua",])),
  K = cbind(unlist(pwm.TriNt["aaa",] + pwm.TriNt["aag",])),
  L = cbind(unlist(pwm.TriNt["uua",] + pwm.TriNt["uug",] + pwm.TriNt["cuu",] + pwm.TriNt["cuc",] + pwm.TriNt["cua",] + pwm.TriNt["cug",])),
  M = cbind(unlist(pwm.TriNt["aug",])),
  N = cbind(unlist(pwm.TriNt["aau",] + pwm.TriNt["aac",])),
  P = cbind(unlist(pwm.TriNt["ccu",] + pwm.TriNt["ccc",] + pwm.TriNt["cca",] + pwm.TriNt["ccg",])),
  Q = cbind(unlist(pwm.TriNt["caa",] + pwm.TriNt["cag",])),
  R = cbind(unlist(pwm.TriNt["cgu",] + pwm.TriNt["cgc",] + pwm.TriNt["cga",] + pwm.TriNt["cgg",] + pwm.TriNt["aga",] + pwm.TriNt["agg",])),
  S = cbind(unlist(pwm.TriNt["ucu",] + pwm.TriNt["ucc",] + pwm.TriNt["uca",] + pwm.TriNt["ucg",] + pwm.TriNt["agu",] + pwm.TriNt["agc",])),
  T = cbind(unlist(pwm.TriNt["acu",] + pwm.TriNt["acc",] + pwm.TriNt["aca",] + pwm.TriNt["acg",])),
  V = cbind(unlist(pwm.TriNt["guu",] + pwm.TriNt["guc",] + pwm.TriNt["gua",] + pwm.TriNt["gug",])),
  W = cbind(unlist(pwm.TriNt["ugg",])),
  Y = cbind(unlist(pwm.TriNt["uau",] + pwm.TriNt["uac",])))
  # check : pwm.TriNt <- limCod.1a ; colSums(limCod.1a) - rowSums(pwm.Cod)
  return(t(pwm.Cod))
  }



# ntcodon.freq.cod : Calculate codon counts at footprint boundaries and generates logoplot
#
#
#
# This function takes a list of samples and experimental conditions in input and
#   generates codon counts at footprint boundaries and a logoplot for each sample.
#
# In:
#        XP.conditions  Vector of experimental conditions for each sample
#        XP.names       Vector of names for each sample
#        pathout        Address where output files will be written
#
# Out:
#       A list containing:
#          - plot                Address of png logoplot file
#          - value               Bit value correspondign to the logoplot
#          - color               Color white/orange/red corresponding to good/warning/poor level of quality
#          - recommendation      Description and recommendation based on value
#
ntcodon.freq.cod <- function(XP.conditions, XP.names, pathout) {

    n <- length(XP.conditions)

    ntcodon.freq.cod.res.i <- c()

    for (i in 1:n) {
        BCOi <- paste(pathout, "Sample-",XP.names[i], sep="")
        ntcodon.freq.cod.res.i[[i]] <- ntcodon.freq.cod.single(BCOi, i, XP.names, pathout)
    }

    ### Summarize over samples
    ntcodon.freq.cod.res <- c()
    ##  Plot
    ntcodon.freq.cod.Plot.summary.pre <- paste(pathout,"Codon-Freq_pre.png",sep="")
    ntcodon.freq.cod.Plot.summary <- paste(pathout,"Codon-Freq.png",sep="")
    # Assemble all samples and legend -density 400 file1.ps file2.ps +append -resize 25%
    #system(paste("montage -mode concatenate -tile 1x",n," ",
    system(paste("montage -mode concatenate -density 800 -tile 1x",n," -resize 800x800 ",
             pathout, "Codon-Freq-*.png ",
             ntcodon.freq.cod.Plot.summary.pre,
             sep=""))
    # Fit image in 800x800 limits
    system(paste("convert ",
                 ntcodon.freq.cod.Plot.summary.pre,
                 " -resize 800x800 ",
                 ntcodon.freq.cod.Plot.summary, sep=""))
    # Remove obsolete files
    system(paste("rm ", pathout, "Codon-Freq-*.png", sep=""))
    #
    ntcodon.freq.cod.res$plot <- ntcodon.freq.cod.Plot.summary
    #
    ## Value
    # Each sample's bit value
    val.summary <- sapply(ntcodon.freq.cod.res.i, function(x) x$value)
    # Return all
    ntcodon.freq.cod.res$value <- val.summary
    #
    ## Color
    # Each sample's color
    col.summary <- sapply(ntcodon.freq.cod.res.i, function(x) x$color)
    # Reduce to highest color level as soon as it is found
    if (sum(col.summary=="red") >= 1)  {col.summary = "red"}
    if (sum(col.summary=="orange") >= 1) {col.summary = "orange"}
    if (sum(col.summary=="green") >= 1)   {col.summary = "green"}
    #
    ntcodon.freq.cod.res$color <- col.summary
    #
    ## Recommendation
    # Each sample's recommendation
    rec.summary.0 <- sapply(ntcodon.freq.cod.res.i, function(x) x$recommendation)
    rec.summary.1 <- paste(XP.names, rec.summary.0, sep = " : ")
    rec.summary.2 <- paste(rec.summary.1, collapse=" \n")
    #
    ntcodon.freq.cod.res$recommendation <- rec.summary.2
    #
    ## Save values for later reference, and return
    save(ntcodon.freq.cod.res, 
         file=paste(pathout, "valNrec_", "ntcodon.freq.cod.res",".RData", sep=""))
    return(ntcodon.freq.cod.res)
}

ntcodon.freq.cod.single <- function(BCO, i, XP.names, pathout) {

    # Input
    limCod.1a <- utils::read.table(paste(BCO,".limCod", sep=""), header=TRUE, row.names = 1)
    limAa.1a <- map.triNt.2.Codon(limCod.1a)

    # Initialize output
    ntcodon.freq.cod.singleres <- c()

    # Prepare
    nr.Aa <- nrow(limAa.1a)
    nc.Aa <- ncol(limAa.1a)

    pwm.pre <- as.matrix(limAa.1a) / matrix(colSums(limAa.1a), nrow = nr.Aa, ncol=nc.Aa, byrow = TRUE)
    pwm.5 <- pwm.pre[,1:(nc.Aa/2)]
    pwm.3.rev <- pwm.pre[,rev((nc.Aa/2 + 1):nc.Aa)] #pwm.3 <- pwm.pre[,(nc.Aa/2 + 1):nc.Aa]

    N<-20
    nseqs <- min(colSums(limAa.1a))
    seq_type <- "other" #"AA"

    p.5 <- CLggseqlogo_pfm(data = pwm.5, N=N, nseqs=nseqs, seq_type=seq_type, 
             method = 'bits', x_lab = paste("codon from 5', ", XP.names[i], sep="") )
    infoContent.5 <- CLcomputeBits(pwm.5, N=N, Nseqs=nseqs)
    p.3 <- CLggseqlogo_pfm( pwm.3.rev, N, nseqs, seq_type, 
             method = 'bits', x_lab = paste("codon from 3', ", XP.names[i], sep="") )
    infoContent.3 <- CLcomputeBits(pwm.3.rev, N=N, Nseqs=nseqs)

    # Plot
    ntcodon.freq.cod.singleres.png <- paste(pathout,"Codon-Freq-",i,".png",sep="") #tempfile()
    ntcodon.freq.cod.singleres.eps <- paste(pathout,"Codon-Freq-",i,".eps",sep="") #tempfile()
    ntcodon.freq.cod.singleres$plot <- ntcodon.freq.cod.singleres.png

    for (plotformat in c("png","eps")) {
      if (plotformat=="png") {
        grDevices::png(filename = ntcodon.freq.cod.singleres.png,width = 600,height = 150)
      } else {
        grDevices::cairo_ps(ntcodon.freq.cod.singleres.eps,width = 10,height = 2)
      }
    #grDevices::png(filename = ntcodon.freq.cod.singleres$plot,width = 1000,height = 250)
      gridExtra::grid.arrange(p.5, p.3, nrow=1)
    grDevices::dev.off()
    } # end loop PNG or EPS

    # Value, Color and Recommendation
    ntcodon.freq.cod.singleres$value <- max(c(infoContent.5, infoContent.3))
    ntcodon.freq.cod.singleres$color <- ifelse(ntcodon.freq.cod.singleres$value < 0.2, yes = "green",
                                             no = ifelse(ntcodon.freq.cod.singleres$value < 0.4, yes = "orange", no = "red"))
    ntcodon.freq.cod.singleres$recommendation <- ifelse(ntcodon.freq.cod.singleres$value < 0.2,
                                              yes = "No position has any overrepresented codons",
                                              no  = ifelse(ntcodon.freq.cod.singleres$value < 0.4,
                                                           yes = paste("Position(s) ",
                                                              list(
                                                              apply(cbind(which( c(infoContent.5, infoContent.3) >= 0.2 )),1, function(x) {
                                                                xx <- x
                                                                if (x>(nc.Aa/2)) { xx <- (x-nc.Aa-1) }
                                                                return(xx)
                                                                })
                                                               ),
                                                              " might have overrepresented codons (note : minus sign means from 3prime end).",
                                                              sep=""),
                                                           no =paste("Position(s) ",
                                                              list(
                                                              apply(cbind(which( c(infoContent.5, infoContent.3) >= 0.2 )),1, function(x) {
                                                                xx <- x
                                                                if (x>(nc.Aa/2)) { xx <- (x-nc.Aa-1) }
                                                                return(xx)
                                                                })
                                                               ),
                                                              " have overrepresented codons  (note : minus sign means from 3prime end) \n",
                                                              "; these codons shouldn't be used for normalisation.",
                                                              sep="")
                                                ))


    ## Save values for later reference, and return
    save(ntcodon.freq.cod.singleres, 
         file=paste(pathout, "valNrec_", "ntcodon.freq.cod.singleres-",i,".RData", sep=""))
    return(ntcodon.freq.cod.singleres)
}

