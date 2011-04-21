\name{MetaQC}
\alias{MetaQC}

\title{
MetaQC: Quantitative Quality Assessment for Inclusion/Exclusion Criteria of Genomic Meta-Analysis
}

\description{
MetaQC implements our proposed quantitative quality control measures: (1) internal homogeneity of co-expression structure among studies (internal quality control; IQC); (2) external consistency of co-expression structure correlating with pathway database (external quality control; EQC); (3) accuracy of differentially expressed gene detection (accuracy quality control; AQCg) or pathway identification (AQCp); (4) consistency of differential expression ranking in genes (consistency quality control; CQCg) or pathways (CQCp). (See the reference for detailed explanation.)
For each quality control index, the p-values from statistical hypothesis testing are minus log transformed and PCA biplots were applied to assist visualization and decision. Results generate systematic suggestions to exclude problematic studies in microarray meta-analysis and potentially can be extended to GWAS or other types of genomic meta-analysis. The identified problematic studies can be scrutinized to identify technical and biological causes (e.g. sample size, platform, tissue collection, preprocessing etc) of their bad quality or irreproducibility for final inclustion/exclusion decision.
}

\usage{
MetaQC(DList, GList, isParallel = FALSE, nCores = NULL, 
	useCache = TRUE, filterGenes = TRUE, verbose = FALSE)
}

\arguments{
  \item{DList}{
A list of all data matrices; Each data name should be set as the name of each list element. Each data should be a numeric matrix that has genes in the rows and samples in the columns. Row names should be official gene symbols and column names be class labels (currently, only two classes are supported).
}
  \item{GList}{
The location of a file which has sets of gene symbol lists such as gmt files. By default, the gmt file will be converted to list object and saved with the same name with ``.rda''. Alternatively, a list of gene sets is allowed; the name of each element of the list should be set as a unique pathway name, and each pathway should have a character vector of gene symbols. 
}
  \item{isParallel}{
Whether to use multiple cores in parallel for fast computing. By default, it is false.
}
  \item{nCores}{
When isParallel is true, the number of cores can be set. By default, all cores in the machine are used in the unix-like machine, and 3 cores are used in windows.
}
  \item{useCache}{
Whether imported gmt file should be saved for the next use. By default, it is true.
}
  \item{filterGenes}{
Whether to use gene filtering (recommended).
}
  \item{verbose}{
Whether to print out logs.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
proto R object.
Use RunQC function to run QC procedure.
Use Plot function to plot PCA figure.
Use Print function to view various information.
See examples below.
}
\references{
Dongwan D. Kang, Etienne Sibille, Naftali Kaminski, and George C. Tseng. (2011) Meta-QC: Quantitative Quality Assessment for Inclusion/Exclusion Criteria of Genomic Meta-Analysis. 
}
\author{
Don Kang (donkang75@gmail.com)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
    library(MetaQC)
    libs <- c("proto", "foreach", ifelse(.Platform$OS.type == "unix", "doMC", "doSMP")) #doSMP is available to R > 2.11.0
    requireAll(libs)

    ## Toy Example
    data(brain) #already hugely filtered
    #Two default gmt files are automatically downloaded, otherwise it is required to locate it correctly.
    #Refer to http://www.broadinstitute.org/gsea/downloads.jsp
    brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=FALSE, verbose=TRUE)
    brainQC$RunQC(B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt") #B is recommended to be >= 1e4 in real application
    brainQC$Plot()
    brainQC$Print()

    ## For parallel computation with all cores
    brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=FALSE, verbose=TRUE, isParallel=TRUE)
    brainQC$RunQC(B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt") #B is recommended to be >= 1e4 in real application
    brainQC$Plot()
    brainQC$Cleanup() #neccessary for windows after using parallel processing

    ## For parallel computation with only 2 cores
    brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=FALSE, verbose=TRUE, isParallel=TRUE, nCores=2)
    brainQC$RunQC(B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt") #B is recommended to be >= 1e4 in real application
    brainQC$Plot()
    brainQC$Cleanup() #neccessary for windows after using parallel processing

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QualityControl }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
