
metageneExpl <- function(inFile, res, res2) {

  PythonInR::pyExecfile(system.file("metageneExpl.py", package="RiboVIEW", mustWork = TRUE))
  PymetageneExpl <- PythonInR::pyFunction("metageneExpl")
  liste <- PymetageneExpl(inFile, res, res2)
  #rPython::python.load( system.file("metageneExpl.py", package="RiboVIEW", mustWork = TRUE))
  #liste <- rPython::python.call("metageneExpl", inFile, res, res2)

  metaGi <- c()
  metaGi$covUTR      <- liste[[1]]
  metaGi$infl        <- liste[[2]]
  metaGi$metaG_m     <- liste[[3]]
  metaGi$metaG_iqr   <- liste[[4]]
  metaGi$metaG_start <- liste[[5]]
  metaGi$leakStart   <- liste[[6]]
  metaGi$leakStop    <- liste[[7]]

  return(metaGi)
  }

