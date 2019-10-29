### R code from vignette source 'Intro_RiboVIEW.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Intro_RiboVIEW.Rnw:143-144 (eval = FALSE)
###################################################
## install.packages("RiboVIEW")


###################################################
### code chunk number 2: Intro_RiboVIEW.Rnw:149-150
###################################################
library(RiboVIEW)


###################################################
### code chunk number 3: Intro_RiboVIEW.Rnw:157-158 (eval = FALSE)
###################################################
## install.packages("PATH/TO/FILE", repos = NULL, type="source")


###################################################
### code chunk number 4: Intro_RiboVIEW.Rnw:222-224 (eval = FALSE)
###################################################
## idir <- paste(system.file(package="RiboVIEW"), "/extdata/",sep="")
## idir


###################################################
### code chunk number 5: Intro_RiboVIEW.Rnw:231-232 (eval = FALSE)
###################################################
## pathout <- paste(tempdir(),"/output-files-vignette-RiboVIEW/",sep="")


###################################################
### code chunk number 6: Intro_RiboVIEW.Rnw:239-240 (eval = FALSE)
###################################################
## pathout <- "/Path/to/my/RiboVIEW/output/directory/"


###################################################
### code chunk number 7: Intro_RiboVIEW.Rnw:246-248 (eval = FALSE)
###################################################
## system(paste("mkdir -p ",pathout,sep=""))
## setwd(pathout)


###################################################
### code chunk number 8: Intro_RiboVIEW.Rnw:256-265 (eval = FALSE)
###################################################
## readsBAM.1.1  <- paste(idir,"Cond1-Rep1.bam",sep="")
## readsBAM.1.2  <- paste(idir,"Cond1-Rep2.bam",sep="")
## readsBAM.1.3  <- paste(idir,"Cond1-Rep3.bam",sep="")
## readsBAM.2.1  <- paste(idir,"Cond2-Rep1.bam",sep="")
## readsBAM.2.2  <- paste(idir,"Cond2-Rep2.bam",sep="")
## readsBAM.2.3  <- paste(idir,"Cond2-Rep3.bam",sep="")
## 
## list.bam <- list(readsBAM.1.1, readsBAM.1.2, readsBAM.1.3, 
##                  readsBAM.2.1, readsBAM.2.2, readsBAM.2.3)


###################################################
### code chunk number 9: Intro_RiboVIEW.Rnw:272-277 (eval = FALSE)
###################################################
## # Reference sequences for mRNAs
## refFASTA=paste(idir,"synth.fasta",sep="")
## 
## # Reference annotation for mRNAs' CDS
## refCDS=paste(idir,"synth.tsv",sep="")


###################################################
### code chunk number 10: Intro_RiboVIEW.Rnw:315-319 (eval = FALSE)
###################################################
## XP.conditions   <- c("cond1","cond1","cond1","cond2", "cond2","cond2")
## XP.conditions.i <- c( 1,1,1,2,2,2)
## XP.names        <- c("C1.R1", "C1.R2", "C1.R3", 
##                      "C2.R1", "C2.R2", "C2.R3")


###################################################
### code chunk number 11: Intro_RiboVIEW.Rnw:332-334 (eval = FALSE)
###################################################
## periodicity(list.bam, refCDS, refFASTA, pathout, XP.names, 
##   versionStrip = FALSE)


###################################################
### code chunk number 12: Intro_RiboVIEW.Rnw:366-367 (eval = FALSE)
###################################################
## attach(listminmax <- select.FPlen(list.bam, pathout, XP.names))


###################################################
### code chunk number 13: Intro_RiboVIEW.Rnw:507-508 (eval = FALSE)
###################################################
## generate.m.s(XP.conditions, XP.names, pathout, B=1000)


###################################################
### code chunk number 14: Intro_RiboVIEW.Rnw:513-528 (eval = FALSE)
###################################################
## visu.m.s.enrichmnt.res <- visu.m.s.enrichmnt(XP.conditions, 
##   XP.names, pathout)
## visu.m.s.enrichmnt.res
## 
## visu.tracks.res <- visu.tracks(XP.conditions, XP.names, pathout,
##   refCDS, mRNA="random", codon.labels=FALSE,
##   codon.col="darkslateblue")
## visu.tracks.res
## 
## Venn.all.res <- Venn.all(XP.names, pathout)
## Venn.all.res
## 
## enricht.aroundA.res <- enricht.aroundA(XP.conditions, 
##   XP.names, pathout)
## enricht.aroundA.res


###################################################
### code chunk number 15: Intro_RiboVIEW.Rnw:576-591 (eval = FALSE)
###################################################
## repl.correl.counts.Venn.res <- repl.correl.counts.Venn(XP.conditions, 
##   XP.names, pathout)
## repl.correl.counts.Venn.res
## 
## repl.correl.gene.res <- repl.correl.gene(XP.conditions, XP.names, 
##   pathout)
## repl.correl.gene.res
## 
## repl.correl.codon.res <- repl.correl.codon(list.bam, refCDS, 
##   refFASTA, mini, maxi, XP.names, XP.conditions, pathout)
## repl.correl.codon.res
## 
## repl.correl.heatmap.res <- repl.correl.heatmap(XP.conditions.i, 
##   XP.names, pathout)
## repl.correl.heatmap.res


###################################################
### code chunk number 16: Intro_RiboVIEW.Rnw:620-631 (eval = FALSE)
###################################################
## chx.artefacts.res <- chx.artefacts(XP.conditions, XP.names, 
##   pathout)
## chx.artefacts.res
## 
## ntcodon.freq.nt.res <- ntcodon.freq.nt(XP.conditions, XP.names, 
##   pathout)
## ntcodon.freq.nt.res
## 
## ntcodon.freq.cod.res <- ntcodon.freq.cod(XP.conditions, 
##   XP.names, pathout)
## ntcodon.freq.cod.res


###################################################
### code chunk number 17: Intro_RiboVIEW.Rnw:650-655 (eval = FALSE)
###################################################
## batch.effects.lm.e.res <- batch.effects.lm.e(XP.conditions, XP.names, pathout)
## batch.effects.lm.e.res
## 
## batch.effects.pca.res <- batch.effects.pca(XP.conditions, XP.names, pathout)
## batch.effects.pca.res


###################################################
### code chunk number 18: Intro_RiboVIEW.Rnw:680-682 (eval = FALSE)
###################################################
## library(gridExtra)
## grid.table(round(batch.effects.lm.res$recommendation.mat, 4)[,1:3])  


###################################################
### code chunk number 19: Intro_RiboVIEW.Rnw:692-693 (eval = FALSE)
###################################################
## batch.effects.lm.res$value


###################################################
### code chunk number 20: Intro_RiboVIEW.Rnw:718-728 (eval = FALSE)
###################################################
## metagene.res <- metagene.all(XP.conditions, XP.names, pathout) 
## 
## metagene.monosome.res <- metagene.res[[1]]
## metagene.monosome.res
## 
## metagene.inflation.res <- metagene.res[[2]]
## metagene.inflation.res
## 
## metagene.leakage.res <- metagene.res[[3]]
## metagene.leakage.res


###################################################
### code chunk number 21: Intro_RiboVIEW.Rnw:741-742 (eval = FALSE)
###################################################
## metagene.leakage.res$value


###################################################
### code chunk number 22: Intro_RiboVIEW.Rnw:769-770 (eval = FALSE)
###################################################
## outputQc(pathout, XP.conditions)    


###################################################
### code chunk number 23: Intro_RiboVIEW.Rnw:790-791 (eval = FALSE)
###################################################
## outputMine(pathout, XP.conditions)  


###################################################
### code chunk number 24: Intro_RiboVIEW.Rnw:842-843 (eval = FALSE)
###################################################
## wdir <- "/Path/to/my/working/directory/"


###################################################
### code chunk number 25: Intro_RiboVIEW.Rnw:866-869 (eval = FALSE)
###################################################
## listminmax <- select.FPlen(list.bam, pathout, XP.names)
## mini <- listminmax[[1]]
## maxi <- listminmax[[2]]


