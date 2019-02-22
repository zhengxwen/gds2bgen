\name{seqBGEN_Info}
\alias{seqBGEN_Info}
\title{BGEN File Information}
\description{
    Gets the information of a bgen file.
}

\usage{
seqBGEN_Info(bgen.fn, verbose=TRUE)
}

\arguments{
    \item{bgen.fn}{the file name of a bgen file}
    \item{verbose}{if \code{TRUE}, show information}
}

\value{
    Return a list with the following components:
    \code{num.sample}: the total number of samples;
    \code{num.variant}: the total number of variants;
    \code{compression}: the compression algorithm, e.g., zlib, zstd;
    \code{layout}: the version of layout;
    \code{sample.id}: NULL (anonymized) or a vector of sample IDs.
}

\references{
    \url{https://bitbucket.org/gavinband/bgen},
    \url{http://bgenformat.org},
    \url{https://github.com/zhengxwen/gds2bgen}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqBGEN2GDS}}
}

\examples{
bgen.fn <- system.file("extdata", "example.8bits.bgen", package="gds2bgen")
seqBGEN_Info(bgen.fn)
}

\keyword{GDS}
\keyword{BGEN}
\keyword{utilities}
