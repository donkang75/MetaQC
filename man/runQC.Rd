\name{runQC}
\alias{runQC}

\title{
Command to execute quality control procedures.
}

\description{
It it a utility function to RunQC method in MetaQC object.
}

\usage{
runQC(QC, nPath=NULL, B=1e4, pvalCut=.05, 
	pvalAdjust=FALSE, fileForCQCp="c2.all.v3.0.symbols.gmt")
}

\arguments{
  \item{QC}{
A proto R object which obtained by MetaQC function.
}
  \item{nPath}{
The number of top pathways which would be used for EQC calculation. The top pathways are automatically determined by their mean rank of over significance among given studies. It is important that gene sets used for EQC are expected to have higher correlation than background. For better performance, this should be set as a reasonably small number.
}
  \item{B}{
The number of permutation tests used for EQC calculation. More than 1e4 is recommended.
}
  \item{pvalCut}{
P-value threshold used for AQC calculation.
}
  \item{pvalAdjust}{
Whether to apply p-value adjustment due to multiple testing (B-H procedure is used).
}
  \item{fileForCQCp}{
Gene set used for CQCp calculation. Usually larger gene set is used than EQC calculation. 
}

}
\value{
A data frame showing a summary of each quality control score.
}
\references{
Dongwan D. Kang, Etienne Sibille, Naftali Kaminski, and George C. Tseng. (Nucleic Acids Res. 2012) MetaQC: Objective Quality Control and Inclusion/Exclusion Criteria for Genomic Meta-Analysis. 
}
\author{
Don Kang (donkang75@gmail.com) and George Tseng (ctseng@pitt.edu)
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
    plot(brainQC)

    ## For parallel computation with only 2 cores
	## R >= 2.11.0 in windows to use parallel computing
    brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", 
			filterGenes=FALSE, verbose=TRUE, isParallel=TRUE, nCores=2)
    #B is recommended to be >= 1e4 in real application
    runQC(brainQC, B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt")
    plot(brainQC)

    ## For parallel computation with all cores
	## In windows, only 2 cores are used if not specified explicitly
    brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", 
			filterGenes=FALSE, verbose=TRUE, isParallel=TRUE)
	#B is recommended to be >= 1e4 in real application					
    runQC(brainQC, B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt")
    plot(brainQC)

	## Real Example which is used in the paper
	#download the brainFull file 
	#from https://github.com/downloads/donkang75/MetaQC/brainFull.rda
	load("brainFull.rda")
    brainQC <- MetaQC(brainFull, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=TRUE, 
			verbose=TRUE, isParallel=TRUE)
    runQC(brainQC, B=1e4, fileForCQCp="c2.all.v3.0.symbols.gmt") #B was 1e5 in the paper
    plot(brainQC)

	## Survival Data Example
	#download Breast data 
	#from https://github.com/downloads/donkang75/MetaQC/Breast.rda
	load("Breast.rda")
    breastQC <- MetaQC(Breast, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=FALSE, 
			verbose=TRUE, isParallel=TRUE, resp.type="Survival")
    runQC(breastQC, B=1e4, fileForCQCp="c2.all.v3.0.symbols.gmt") 
    breastQC
    plot(breastQC)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QualityControl }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
