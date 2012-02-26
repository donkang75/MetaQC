\name{cleanup}
\alias{cleanup}

\title{
Cleaning up resources.
}

\description{
It is to shutdown the workers used for parallel processing and release resources. It is only necessary in windows. (Deprecated)
}

\usage{
cleanup(QC)
}

\arguments{
  \item{QC}{
A proto R object which obtained by MetaQC function.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
NA
}
\references{
Dongwan D. Kang, Etienne Sibille, Naftali Kaminski, and George C. Tseng. (Nucleic Acids Res. 2012) Meta-QC: Objective Quality Control and Inclusion/Exclusion Criteria for Genomic Meta-Analysis. 
}
\author{
Don Kang (donkang75@gmail.com) and George Tseng (ctseng@pitt.edu)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MetaQC}}
%\code{\link{plot.catdes}},\code{\link{condes}}
}
\examples{
\dontrun{
    requireAll(c("proto", "foreach"))

   ## Toy Example
    data(brain) #already hugely filtered
    #Two default gmt files are automatically downloaded, 
	#otherwise it is required to locate it correctly.
    #Refer to http://www.broadinstitute.org/gsea/downloads.jsp
    ## For parallel computation with only 2 cores
	## R >= 2.11.0 in windows to use parallel computing
    brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", 
			filterGenes=FALSE, verbose=TRUE, isParallel=TRUE, nCores=2)
    #B is recommended to be >= 1e4 in real application
	runQC(brainQC, B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt")
    plot(brainQC)
    cleanup(brainQC) #neccessary for windows after using parallel processing

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QualityControl }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
