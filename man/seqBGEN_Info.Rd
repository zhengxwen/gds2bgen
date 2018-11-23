\name{seqBGEN_Info}
\alias{seqBGEN_Info}
\title{Add a new GDS node}
\description{
    Add a new GDS node to the GDS file.
}

\usage{
seqBGEN_Info(bgen.fn, verbose=TRUE)
}

\arguments{
    \item{bgen.fn}{the file name of a bgen file}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    XX
}

\value{
    An object of class \code{\link{gdsn.class}} of the new node.
}

\references{
    \url{http://github.com/zhengxwen/gdsfmt}
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
