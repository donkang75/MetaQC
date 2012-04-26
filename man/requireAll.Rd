\name{requireAll}
\alias{requireAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MetaQC: Quantitative Quality Assessment for Inclusion/Exclusion Criteria of Genomic Meta-Analysis
}
\description{
requireAll description
}
\usage{
requireAll(packages)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{packages}{
A character vector of required packages. Unavailable packages are going to be installed.
}
}
\value{
None
}
\author{
Don Kang (donkang75@gmail.com) and George Tseng (ctseng@pitt.edu)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
	\dontrun{
	    libs <- c("proto", "foreach", ifelse(.Platform$OS.type == "unix", 
					"doMC", "doSNOW")) 
	    requireAll(libs)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ require }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
