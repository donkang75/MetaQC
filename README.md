MetaQC: Quantitative Quality Assessment for Inclusion/Exclusion Criteria of Genomic Meta-Analysis
=================================================================================================

Introduction
------------
__Genomic meta-analyses__ for combining microarray studies have been widely applied to increase statistical power and to validate results from individual studies. Currently, the inclusion/exclusion criteria mostly depend on ad-hoc expert opinion or naïve decision by sample size or array platform. No objective quality assessment is available. 

__MetaQC__ implements our proposed quantitative quality control measures: 

* (1) internal homogeneity of co-expression structure among studies (internal quality control; __IQC__)
* (2) external consistency of co-expression structure correlating with pathway database (external quality control; __EQC__)
* (3) accuracy of differentially expressed gene detection (accuracy quality control; __AQCg__) or pathway identification (__AQCp__)
* (4) consistency of differential expression ranking in genes (consistency quality control; __CQCg__) or pathways (__CQCp__). 

For each quality control index, the p-values from statistical hypothesis testing are minus log transformed and PCA biplots were applied to assist visualization and decision. Results generate systematic suggestions to exclude problematic studies in microarray meta-analysis and potentially can be extended to GWAS or other types of genomic meta-analysis. The identified problematic studies can be scrutinized to identify technical and biological causes (e.g. sample size, platform, tissue collection, preprocessing etc) of their bad quality or irreproducibility for final inclusion/exclusion decision.

Installation
--------
To install this package, download a proper file for the target OS. Then run:

### Windows            
`> install.packages("MetaQC_0.1.1.zip", repos=NULL, type="win.binary")`

### Mac OS X            
`> install.packages("MetaQC_0.1.1.tgz", repos=NULL, type="mac.binary")`

### Linux            
`> install.packages("MetaQC_0.1.1.tar.gz", repos=NULL, type="source")`

Examples
--------
library(MetaQC)
libs <- c("proto", "foreach", ifelse(.Platform$OS.type == "unix", "doMC", "doSMP"))
requireAll(libs)

## Toy Example
data(brain) #already hugely filtered
brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", 
					filterGenes=FALSE, verbose=TRUE)
brainQC$RunQC(B=1e2) #need to be more than this (at least 1e4)
brainQC$Plot()
brainQC$Print()

## For parallel computation with all cores
brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", 
					filterGenes=FALSE, verbose=TRUE, isParallel=TRUE)
brainQC$Cleanup() #neccessary for windows after using parallel processing

## For parallel computation with only 2 cores
brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", 
					filterGenes=FALSE, verbose=TRUE, isParallel=TRUE, nCores=2)
}
brainQC$Cleanup() #neccessary for windows after using parallel processing

References
----------
Dongwan D. Kang, Etienne Sibille, Naftali Kaminski, and George C. Tseng. (2011) Meta-QC: Quantitative Quality Assessment for Inclusion/Exclusion Criteria of Genomic Meta-Analysis. 
