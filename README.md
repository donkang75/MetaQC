MetaQC: Objective Quality Control and Inclusion/Exclusion Criteria for Genomic Meta-Analysis
============================================================================

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
--------------
To install this package, run:

        source("http://bioconductor.org/biocLite.R")
        biocLite('MetaQC')

Examples
-------------
	library(MetaQC)
	
	## Toy Example
	data(brain) #already hugely filtered
	
	#Two default gmt files are automatically downloaded, 
	#otherwise it is required to locate it correctly.
	#Refer to http://www.broadinstitute.org/gsea/downloads.jsp
	brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=FALSE, verbose=TRUE)
	#B is recommended to be >= 1e4 in real application					
	runQC(brainQC, B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt") 
	brainQC
	plot(brainQC)
	
	## For parallel computation with only 2 cores
	brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=FALSE, verbose=TRUE, isParallel=TRUE, nCores=2)
	#B is recommended to be >= 1e4 in real application
	runQC(brainQC, B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt") 
	plot(brainQC)
	
	## For parallel computation with all cores
	## In windows, only 2 cores are used if not specified explicitly
	brainQC <- MetaQC(brain, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=FALSE, verbose=TRUE, isParallel=TRUE)
	#B is recommended to be >= 1e4 in real application					
	runQC(brainQC, B=1e2, fileForCQCp="c2.all.v3.0.symbols.gmt") 
	plot(brainQC)
	
	## Real Example which is used in the paper
	#download the brainFull file from https://github.com/downloads/donkang75/MetaQC/brainFull.rda
	load("brainFull.rda")
	brainQC <- MetaQC(brainFull, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=TRUE, verbose=TRUE, isParallel=TRUE)
	runQC(brainQC, B=1e4, fileForCQCp="c2.all.v3.0.symbols.gmt") #B was 1e5 in the paper 
	plot(brainQC)
	
	## Survival Data Example
	#download Breast data 
	#from https://github.com/downloads/donkang75/MetaQC/Breast.rda
	load("Breast.rda")
	breastQC <- MetaQC(Breast, "c2.cp.biocarta.v3.0.symbols.gmt", filterGenes=FALSE, verbose=TRUE, isParallel=TRUE, resp.type="Survival")
	runQC(breastQC, B=1e4, fileForCQCp="c2.all.v3.0.symbols.gmt") 
	breastQC
	plot(breastQC)


References
----------
Dongwan D. Kang, Etienne Sibille, Naftali Kaminski, and George C. Tseng. (Nucleic Acids Res. 2012) MetaQC: Objective Quality Control and Inclusion/Exclusion Criteria for Genomic Meta-Analysis. ([doi: 10.1093/nar/gkr1071](http://nar.oxfordjournals.org/cgi/content/abstract/gkr1071))
