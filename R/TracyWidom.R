#Distribution for eigenvalues

TracyWidom <- function(variances, nv, no) {
# Inputs for R function TracyWidom :
# - Variances adding up to 1, corresponding to eigenvalues from a PCA or SVD
# - number of variables (equal to number of variances)
# - number of observations in matrix from which pca and variances were computed.

  # Load interface to python function
  PythonInR::pyExecfile(system.file("TracyWidom.py", package="RiboVIEW", mustWork = TRUE))
  TracyWidom <- PythonInR::pyFunction("TracyWidom")
  #rPython::python.load( system.file("TracyWidom.py", package="RiboVIEW", mustWork = TRUE))
  
  # Inputs for Python function TracyWidom :
  # - Variances adding up to 1, corresponding to eigenvalues from a PCA or SVD
  # - number of variables (equal to number of variances)
  # - number of observations in matrix from which pca and variances were computed.
  # - name of output file
  # - optionally, tol = tolerance for sum of variances
  outTW <- tempfile()
  TracyWidom(variances, nv, no, outTW)
  #rPython::python.call("TracyWidom", variances, nv, no, outTW)

  # Read output file into a table and return
  TW.res <- utils::read.table(outTW, header=TRUE)
  unlink(outTW)
  return(TW.res)
  }

# exampleTW <- TracyWidom(variances, nv, no)
