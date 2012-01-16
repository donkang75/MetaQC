\name{print.proto}
\alias{print.proto}

\title{
Print MetaQC results.
}

\description{
It prints out the results of all QC measures and standardized mean rank of each study. CQCg and AQCg are combined to be CAQCg, and CQCp and AQCp are combined to be CAQCp to reduce the dominance of CQC and AQC due to their greater correlation.
}

\usage{
print.proto(x, ...)
}

\arguments{
  \item{x}{
A proto R object which obtained by MetaQC function.
}
  \item{...}{
Further arguments to print function.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
NA
}
\references{
Dongwan D. Kang, Etienne Sibille, Naftali Kaminski, and George C. Tseng. (Nucleic Acids Res. 2012) MetaQC: Objective Quality Control and Inclusion/Exclusion Criteria for Genomic Meta-Analysis. 
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
    brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", 
						filterGenes=FALSE, verbose=TRUE)
	#B is recommended to be >= 1e4 in real application					
    runQC(brainQC, B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt")
    brainQC

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QualityControl }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
