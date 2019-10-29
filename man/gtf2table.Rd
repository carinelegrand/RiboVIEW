\name{gtf2table}
\alias{gtf2table}
\title{Convert an annotation file in GTF format in tabular format as required by RiboVIEW}
\usage{
gtf2table(refGTF, outfile,  
          endinside1 = 1, stopoutside1 = 1, stopminusoutside1 = 1,
          verbose=FALSE)
}
\arguments{
\item{refGTF}{Coding sequence annotation in GTF format}
\item{outfile}{File name for output file in tabular format}
\item{endinside1}{Optional indicator, 0 if field end (5' side) is first base outside of feature, in GTF file,
        1 if field end (5' side) is  last base in feature (exon, CDS, transcript, etc.).}
\item{stopoutside1}{Optional indicator, 0 if stop is included in CDS, in GTF file, in '+' strand case.
        1 if stop is outside the CDS as annotated in GTF file, in '+' strand case.}
\item{stopminusoutside1}{Optional indicator, 0 if stop is included in CDS, in GTF file, in '-' strand case.
        1 if stop is outside the CDS as annotated in GTF file, in '-' strand case.}
\item{verbose}{Print or not detailed messages (defaults to FALSE).}
}
\value{
Output file with CDS annotation in tabular format is written to the 
 file name "outfile".
  }
\description{
\code{gtf2table} This function converts coding sequence annotation from GTF format to tabular format
 as required for RiboVIEW. In case GTF file format differs slightly from the expected, options
 endindside1, etc. allow to adjust this script. This should be facilitated using 
 verbose=TRUE. Also, this home-made GTF parser is typically not very robust to unexpected
 characteristics in the attributes field (9th field, containing several semicolon-separated 
 attributes). In this case, setting verbose=TRUE will help settling the problem.
}
\details{
This function relies on Python script gtf2table.py.
}
\examples{

refGTF <- paste(system.file(package="RiboVIEW", mustWork = TRUE), "/extdata/example.gtf.gz", sep="")
refCDS <- tempfile()
gtf2table(refGTF, outfile = refCDS)

cds <- read.table(refCDS)
head(cds)

}

