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
    An object of class \code{\link{gdsn.class}} of the new node.
}

\references{
    \url{https://bitbucket.org/gavinband/bgen},
    \url{http://bgenformat.org}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqBGEN2GDS}}
}

\examples{

# delete the temporary file
unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{BGEN}
\keyword{utilities}
