# Web output to visualize mining results in RiboVIEW
#
#
#
# This function takes a list of sample conditions as input and 
#   generates a web page for visualisation of mining results.
#
# In:
#        pathout        Address where output files will be written
#        XP.conditions  Vector of experimental conditions for each sample
#
# Out:
#       The file "Results-Mine.html", readable in a suitable internet browser,
#         is written to pathout. Compatible browsers are up-to-date version of
#         firefox, safari, chrome and IE.
#


outputMine <- function(pathout, XP.conditions) {
  
  outFile   <- paste(pathout, "Results-Mine.html", sep="")
  
  # Bind to python script for text image generation
  PythonInR::pyExecfile(system.file("base64png.py", package="RiboVIEW", mustWork = TRUE))
  Pybase64png <- PythonInR::pyFunction("base64png")
  #rPython::python.load( system.file("base64png.py", package="RiboVIEW", mustWork = TRUE))

  # Initialization of variables metagene.res, etc. (later erased by RData import)
  metagene.res <- c()
  metagene.res[[1]] <- c()
  metagene.res[[2]] <- c()
  metagene.res[[3]] <- c()
  visu.m.s.enrichmnt.res <- c()
  visu.tracks.res <- c()
  Venn.all.res <- c()
  batch.effects.lm.e.res <- c()
  batch.effects.pca.res <- c()
  enricht.aroundA.res <- c()
  
  # Load png, value and text by category from RData
  load(paste(pathout, "valNrec_","periodicity",".RData",sep=""))
  load(paste(pathout, "valNrec_", "batch.effects.lm.e.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "batch.effects.pca.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "chx.artefacts.res",".RData", sep=""))
  load(paste(pathout, "valNrec_", "enricht.aroundA.res",".RData", sep=""))
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
  load(paste(pathout, "valNrec_", "Venn.all",".RData", sep=""))
  metagene.monosome.res <- metagene.res[[1]]
  metagene.inflation.res <- metagene.res[[2]]
  metagene.leakage.res <- metagene.res[[3]]

  # Panels, tabs, plots and text
  Panels <- c("Within conditions", "Between conditions", "Codon enrichment")
  #
  Tabs <- c()
  Tabs[["Within conditions"]]  <- c("Enrichment by sample", "Nucleobases", "Tracks")
  Tabs[["Between conditions"]] <- c("Comparisons between conditions", "Venn diagram", "PCA", "tSNE")
  Tabs[["Codon enrichment"]]   <- c("A-site +/-90", "A-site +/-15")
  #
  Plots <- c()
  Plots[["Enrichment by sample"]]           <- visu.m.s.enrichmnt.res$plot.enrich
  Plots[["Comparisons between conditions"]] <- visu.m.s.enrichmnt.res$plot.compar
  Plots[["Tracks"]]                         <- visu.tracks.res$plot
  Plots[["Venn diagram"]]                   <- Venn.all.res$plot #repl.correl.counts.Venn.res$plot
  Plots[["Nucleobases"]]                    <- batch.effects.lm.e.res$plot
  Plots[["PCA"]]                            <- batch.effects.pca.res$plot.pca
  Plots[["tSNE"]]                           <- batch.effects.pca.res$plot.tsne
  Plots[["A-site +/-90"]]                 <- enricht.aroundA.res$plot.extract
  Plots[["A-site +/-15"]]                 <- enricht.aroundA.res$plot.zoom.extract
  #
  Texts <- c()
  Texts[["Enrichment by sample"]] <- c("")
  Texts[["Comparisons between conditions"]] <- c("")
  Texts[["Tracks"]] <- paste(c("Obtain this plot for any mRNA using visu.tracks"),
    c("(XP.conditions, XP.names, pathout, refCDS, mRNA=<my-preferred-mRNA>, codon.labels=FALSE)"), sep="")
  Texts[["Venn diagram"]] <- c("")
  Texts[["Nucleobases"]] <- batch.effects.lm.e.res$recommendation
  Texts[["PCA"]] <- batch.effects.pca.res$recommendation
  Texts[["tSNE"]] <- c("")
  Texts[["A-site +/-90"]] <- paste("Codon enrichment around A +/- 90 codons, for NUN codons. ",
                            "Full plots ",
                            "can be found in the output folder, with file names : ",
                            "enrichment-all_*.eps (or enrichment-all_*.png)",sep="")
  Texts[["A-site +/-15"]] <- paste("Codon enrichment around A +/- 15 codons, for NUN codons. ",
                            "Full plots ",
                            "can be found in the output folder, with file names : ",
                            "enrichment-all-zoom_*.eps (or enrichment-all-zoom_*.png)",sep="")

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
  write('<br><h0 class=riboview >RiboMINE</h0><br>', outFile, append=TRUE)
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
  

