cleanup <- function(QC) {
	if(!is.null(QC$.workers))
		stopWorkers(QC$.workers)
}
