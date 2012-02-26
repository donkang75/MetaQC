\name{MetaQC}
\alias{MetaQC}

\title{
MetaQC: Objective Quality Control and Inclusion/Exclusion Criteria for Genomic Meta-Analysis
}

\description{
MetaQC implements our proposed quantitative quality control measures: (1) internal homogeneity of co-expression structure among studies (internal quality control; IQC); (2) external consistency of co-expression structure correlating with pathway database (external quality control; EQC); (3) accuracy of differentially expressed gene detection (accuracy quality control; AQCg) or pathway identification (AQCp); (4) consistency of differential expression ranking in genes (consistency quality control; CQCg) or pathways (CQCp). (See the reference for detailed explanation.)
For each quality control index, the p-values from statistical hypothesis testing are minus log transformed and PCA biplots were applied to assist visualization and decision. Results generate systematic suggestions to exclude problematic studies in microarray meta-analysis and potentially can be extended to GWAS or other types of genomic meta-analysis. The identified problematic studies can be scrutinized to identify technical and biological causes (e.g. sample size, platform, tissue collection, preprocessing etc) of their bad quality or irreproducibility for final inclusion/exclusion decision.
}

\usage{
MetaQC(DList, GList, isParallel = FALSE, nCores = NULL, 
	useCache = TRUE, filterGenes = TRUE, 
	maxNApctAllowed=.3, cutRatioByMean=.4, cutRatioByVar=.4, minNumGenes=5,
	verbose = FALSE)
}

\arguments{
  \item{DList}{
A list of all data matrices; Each data name should be set as the name of each list element. Each data should be a numeric matrix that has genes in the rows and samples in the columns. Row names should be official gene symbols and column names be class labels (currently, only two classes are supported).
}
  \item{GList}{
The location of a file which has sets of gene symbol lists such as gmt files. By default, the gmt file will be converted to list object and saved with the same name with ".rda". Alternatively, a list of gene sets is allowed; the name of each element of the list should be set as a unique pathway name, and each pathway should have a character vector of gene symbols. 
}
  \item{isParallel}{
Whether to use multiple cores in parallel for fast computing. By default, it is false.
}
  \item{nCores}{
When isParallel is true, the number of cores can be set. By default, all cores in the machine are used in the unix-like machine, and 2 cores are used in windows.
}
  \item{useCache}{
Whether imported gmt file should be saved for the next use. By default, it is true.
}
  \item{filterGenes}{
Whether to use gene filtering (recommended).
}
  \item{maxNApctAllowed}{
Filtering out genes which have missing values more than specified ratio (Default .3). Applied if filterGenes is TRUE.
}
  \item{cutRatioByMean}{
Filtering out specified ratio of genes which have least expression value (Default .4). Applied if filterGenes is TRUE.
}
  \item{cutRatioByVar}{
Filtering out specified ratio of genes which have least sample wise expression variance (Default .4). Applied if filterGenes is TRUE. 
}
  \item{minNumGenes}{
Mininum number of genes in a pathway. A pathway which has members smaller than the specified value will be removed.
}
  \item{verbose}{
Whether to print out logs.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A proto R object.
Use RunQC function to run QC procedure.
Use Plot function to plot PCA figure.
Use Print function to view various information.
See examples below.
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
\code{\link{runQC}}
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

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QualityControl }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
