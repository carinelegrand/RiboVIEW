# Convert an annotation file in GTF format in tabular format as required by RiboVIEW
#
#
#
# This function converts coding sequence annotation from GTF format to tabular format
#   as required for RiboVIEW.
#
# In:
#        refGTF             Coding sequence annotation in GTF format
#        outfile            File name for output file in tabular format
#        endinside1         0 if field end (5th) is first base outside of feature, in GTF file,
#         (optional)        1                        last base in feature (exon, CDS, transcript, etc.).
#       
#        stopoutside1       0 if stop is included in CDS, in GTF file, in '+' strand case.
#         (optional)        1            outside the CDS as annotated in GTF file, in '+' strand case.
#       
#        stopminusoutside1  0 if stop is included in CDS, in GTF file, in '-' strand case.
#         (optional)        1            outside the CDS as annotated in GTF file, in '-' strand case.
#       
#
# Out:
#       Output file with CDS annotation in tabular format is written to the 
#         file name "outfile".
#
gtf2table <- function(refGTF, 
                      outfile,  
                      endinside1 = 1, 
                      stopoutside1 = 1, 
                      stopminusoutside1 = 1,
                      verbose=FALSE) {
  filepy <- system.file("gtf2table.py", package="RiboVIEW", mustWork = TRUE)
  
  PythonInR::pyExecfile(filepy)
  Pygtf2table <- PythonInR::pyFunction("gtf2table")
  Pygtf2table(refGTF, 
              outfile, 
              endinside1, 
              stopoutside1, 
              stopminusoutside1,
              verbose=verbose)
  #rPython::python.load(filepy)
  #rPython::python.call("gtf2table", refGTF, outfile, endinside1, stopoutside1, stopminusoutside1)
  
  message(paste("A table of CDS annotation was written to file : \n\t",
                outfile,
                "\n from input GTF file : \n\t",
                refGTF,".",sep=""))
}


