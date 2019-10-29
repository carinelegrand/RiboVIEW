# Web output to visualize quality control diagnostics in RiboVIEW
#
#
#
# This function takes a list of sample conditions as input and 
#   generates a web page for visualisation of quality control diagnostics.
#
# In:
#        pathout        Address where output files will be written
#        XP.conditions  Vector of experimental conditions for each sample
#
# Out:
#       The file "Results-Qc.html", readable in a suitable internet browser,
#         is written to pathout. Compatible browsers are up-to-date version of
#         firefox, safari, chrome and IE.
#
outputQc <- function(pathout, XP.conditions) {
  
  outFile   <- paste(pathout, "Results-Qc.html", sep="")
  
  # Bind to python script for text image generation
  PythonInR::pyExecfile(system.file("base64png.py", package="RiboVIEW", mustWork = TRUE))
  Pybase64png <- PythonInR::pyFunction("base64png")
  #rPython::python.load( system.file("base64png.py", package="RiboVIEW", mustWork = TRUE))
  
  # Initialization of variables metagene.res, etc. (later erased by RData import)
  metagene.res <- c()
  periodicity.res <- c()
  repl.correl.heatmap.res <- c()
  repl.correl.codon.res <- c()
  repl.correl.gene.res <- c()
  ntcodon.freq.nt.res <- c()
  ntcodon.freq.cod.res <- c()
  chx.artefacts.res <- c()
  
  # Load png, value and text by category from RData
  load(paste(pathout, "valNrec_","periodicity",".RData",sep=""))
  load(paste(pathout, "valNrec_", "batch.effects.lm.e.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "batch.effects.pca.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "chx.artefacts.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "metagene.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "ntcodon.freq.nt.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "ntcodon.freq.nt.singleres-",1,".RData", sep=""))
  load(paste(pathout, "valNrec_", "ntcodon.freq.cod.res",".RData", sep=""))
  #load(paste(pathout, "valNrec_", "ntcodon.freq.cod.singleres-",i,".RData", sep=""))
  load(paste(pathout, "valNrec_", "repl.correl.counts.Venn",".RData", sep=""))
  load(paste(pathout, "valNrec_","repl.correl.gene.res",".RData",sep=""))
  #load(paste(pathout, "valNrec_", "repl.correl.codon.singleres-", codon.res,".RData", sep=""))
  load(paste(pathout, "valNrec_", "repl.correl.codon.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "repl.correl.heatmap.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "visu.m.s.enrichmnt.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "visu.tracks.res",".RData", sep=""))
  metagene.monosome.res <- metagene.res[[1]]
  metagene.inflation.res <- metagene.res[[2]]
  metagene.leakage.res <- metagene.res[[3]]

  # Panels, tabs, plots and text
  Panels <- c("Periodicity", "Replicates", "Footprints", "Drugs")
  #
  Tabs <- c()
  Tabs[["Periodicity"]] <- c("Recurrence", "Coverage")
  Tabs[["Replicates"]]  <- c("Heatmap and clustering", "Correlation between codons", "Correlation between genes")
  Tabs[["Footprints"]]  <- c("Ligation bias (nt logo)", "Ligation bias (codon logo)", "Metagene (Monosome selection, drop-off)")
  Tabs[["Drugs"]]       <- c("Cycloheximide", "Inflation at start codon", "Leakage")
  #
  Plots <- c()
  Plots[["Recurrence"]] <- periodicity.res$plot.rec
  Plots[["Coverage"]] <- periodicity.res$plot.cov
  Plots[["Heatmap and clustering"]] <- repl.correl.heatmap.res$plot
  Plots[["Correlation between codons"]] <- repl.correl.codon.res$plot
  Plots[["Correlation between genes"]] <- repl.correl.gene.res$plot
  Plots[["Ligation bias (nt logo)"]] <- ntcodon.freq.nt.res$plot
  Plots[["Ligation bias (codon logo)"]] <- ntcodon.freq.cod.res$plot
  Plots[["Metagene (Monosome selection, drop-off)"]] <- metagene.monosome.res$plot
  Plots[["Cycloheximide"]] <- chx.artefacts.res$plot
  Plots[["Inflation at start codon"]] <- metagene.inflation.res$plot
  Plots[["Leakage"]] <- metagene.leakage.res$plot
  #
  Texts <- c()
  Texts[["Recurrence"]] <- paste("You selected the following footprint lengths : from ",
                                 periodicity.res$mini, "nt to ", periodicity.res$maxi,"nt.", sep="")
  Texts[["Coverage"]] <- Texts[["Recurrence"]]
  Texts[["Heatmap and clustering"]] <- repl.correl.heatmap.res$recommendation
  Texts[["Correlation between codons"]] <- repl.correl.codon.res$recommendation
  Texts[["Correlation between genes"]] <- repl.correl.gene.res$recommendation
  Texts[["Ligation bias (nt logo)"]] <- ntcodon.freq.nt.res$recommendation
  Texts[["Ligation bias (codon logo)"]] <- ntcodon.freq.cod.res$recommendation
  Texts[["Metagene (Monosome selection, drop-off)"]] <- metagene.monosome.res$recommendation
  Texts[["Cycloheximide"]] <- chx.artefacts.res$recommendation
  Texts[["Inflation at start codon"]] <- metagene.inflation.res$recommendation
  Texts[["Leakage"]] <- metagene.leakage.res$recommendation
  #
  Vals <- c()
  Vals[["Heatmap and clustering"]] <- repl.correl.heatmap.res$value
  Vals[["Correlation between codons"]] <- repl.correl.codon.res$value
  Vals[["Correlation between genes"]] <- repl.correl.gene.res$value


  # Metadata, title
  write('<!DOCTYPE html>', outFile, append=FALSE)
  write('<html>', outFile, append=TRUE)
  write('<head>', outFile, append=TRUE)
  write('<meta charset="UTF-8" />', outFile, append=TRUE)
  write('<title>RiboVIEW</title>', outFile, append=TRUE)
  # Styles
  #NOK : styles <- system.file("output-style.css", package="RiboVIEW", mustWork = TRUE)
  styles <- paste( system.file(package="RiboVIEW", mustWork = TRUE), "/extdata/output-style.css", sep="")
  file.append( outFile, styles)
  write('</head>', outFile, append=TRUE)
  # Body
  write('<body bgcolor="#FBFCFA">', outFile, append=TRUE)
  write('<div id="header">', outFile, append=TRUE)
  write('<br><h0 class=riboview >RiboQC</h0><br>', outFile, append=TRUE)
  write('<nav class="header" id="my_centered_buttons">', outFile, append=TRUE)
  
  
  # Buttons for panels
  for (i in 1:length(Panels)) {
     write(paste('<a href="#A',i,'"  class="StyleButton">', sep=''), outFile, 
           append=TRUE)
     write('<button class="button btn5">', outFile, append=TRUE)
     write(Panels[i], outFile, append=TRUE)
     write('</button>', outFile, append=TRUE)
     write('</a>', outFile, append=TRUE)
  }
  write('</nav>', outFile, append=TRUE)
  write('</div>', outFile, append=TRUE)
  
  # Panels
  for (i in 1:length(Panels)) {
     write(paste('<!-- ',Panels[i],' -->', sep=''), outFile, append=TRUE)
     write(paste('<h4 id="A',i,'"> <br> <br> <br> <br> <br> <br> </h4>', sep=''), 
           outFile, append=TRUE)
     write('<div align="center"><form>', outFile, append=TRUE)
     # Tabs
     for (j in 1:length(Tabs[[i]])) {
        if (j==1) {taille='small' ;  checkradio=' checked="checked"'
        } else if (j==2) {taille='medium' ; checkradio=''
        } else if (j==3) {taille='large' ; checkradio=''
        } else if (j==4) {taille='xlarge' ; checkradio=''
        } else {warnings("More categories than expected in output-page.r")
        } 
        write(paste('<div><input type="radio" name="size',i,'" id="',taille,i,
                    '" value="',taille,'"',checkradio,' />', sep=''), outFile, 
              append=TRUE)
        write(paste('<label for="',taille,i,'"><div class="rotate">',Tabs[[i]][j],
                    '</div></label>', sep=''), outFile, append=TRUE)
        write(paste('<article><p>',Tabs[[i]][j],'</p>', sep=''), outFile, append=TRUE)
        write('<p> <img class="floatC" src="data:image/png;base64,', outFile, 
              append=TRUE)
        # base64 image from Python
        img   <- Plots[[Tabs[[i]][j]]]
        img64 <- Pybase64png(img)
        #img64 <- rPython::python.call("base64png", img)
        write(img64, outFile, append=TRUE)
        #
        write(paste('" alt="Missing image for : ',Tabs[[i]][j],'" /></p>',sep=''), 
              outFile, append=TRUE)
        write(paste('<p class="StyleD">', Texts[[Tabs[[i]][j]]], '</p>',sep=''), 
              outFile, append=TRUE)
        write('</article></div>', outFile, append=TRUE)
     }
     write('</form></div>', outFile, append=TRUE)
  }     
  write('<br><br><br><br><br><br>', outFile, append=TRUE)
  write('<div id="footer">', outFile, append=TRUE)
  write('RiboVIEW version 1.0', outFile, append=TRUE)
  write('</div>', outFile, append=TRUE)
  write('</body>', outFile, append=TRUE)
  write('</html>', outFile, append=TRUE)
 
}

