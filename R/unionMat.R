# Merges input matrices by union of the rows and concatenation of columns
#
#
# This function takes a list of matrices in input and merges them by union of 
#     row names and concatenation of columns.
#
# In:
#        listInputFilesR   list of files containing each an input matrix
#        outFile           name of output file for merged matrix
#
# Out:
#        none - output has to be aread from file outFile  
#
unionMat <- function(listInputFilesR, outFile) {

  listInputFiles = paste(listInputFilesR, collapse=";")

  PythonInR::pyExecfile(system.file("unionMat.py", package="RiboVIEW", mustWork = TRUE))
  PyunionMat <- PythonInR::pyFunction("unionMat")
  PyunionMat(listInputFiles, outFile)
  #rPython::python.load( system.file("unionMat.py", package="RiboVIEW", mustWork = TRUE))
  #rPython::python.call("unionMat", listInputFiles, outFile)

  }
