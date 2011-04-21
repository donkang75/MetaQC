\name{MetaQC-package}
\alias{MetaQC-package}
\docType{package}
\title{
MetaQC: Quantitative Quality Assessment for Inclusion/Exclusion Criteria of Genomic Meta-Analysis
}
\description{
MetaQC implements our proposed quantitative quality control measures: (1) internal homogeneity of co-expression structure among studies (internal quality control; IQC); (2) external consistency of co-expression structure correlating with pathway database (external quality control; EQC); (3) accuracy of differentially expressed gene detection (accuracy quality control; AQCg) or pathway identification (AQCp); (4) consistency of differential expression ranking in genes (consistency quality control; CQCg) or pathways (CQCp). (See the reference for detailed explanation.)
For each quality control index, the p-values from statistical hypothesis testing are minus log transformed and PCA biplots were applied to assist visualization and decision. Results generate systematic suggestions to exclude problematic studies in microarray meta-analysis and potentially can be extended to GWAS or other types of genomic meta-analysis. The identified problematic studies can be scrutinized to identify technical and biological causes (e.g. sample size, platform, tissue collection, preprocessing etc) of their bad quality or irreproducibility for final inclustion/exclusion decision.
}
\details{
\tabular{ll}{
Package: \tab MetaQC\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2011-04-21\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Don Kang <donkang75@gmail.com>
}
\references{
Dongwan D. Kang, Etienne Sibille, Naftali Kaminski, and George C. Tseng. (2011) Meta-QC: Quantitative Quality Assessment for Inclusion/Exclusion Criteria of Genomic Meta-Analysis. 
}
\keyword{ QualityControl }
\keyword{ MetaAnalysis }% __ONLY ONE__ keyword per line
\keyword{ Microarray }
\seealso{
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
