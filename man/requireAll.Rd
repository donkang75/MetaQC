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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
None
}
\references{
%% ~put references to the literature/web site here ~
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
	    libs <- c("proto", "foreach", ifelse(.Platform$OS.type == "unix", 
					"doMC", "doSMP")) #doSMP is available to R > 2.11.0
	    requireAll(libs)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ require }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
