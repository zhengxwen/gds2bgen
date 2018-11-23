\name{seqBGEN2GDS}
\alias{seqBGEN2GDS}
\title{Add a new GDS node}
\description{
    Add a new GDS node to the GDS file.
}

\usage{
seqBGEN2GDS(bgen.fn, out.fn, storage.option="LZMA_RA",
    float.type=c("packed", "packed2", "single", "double"),
    dosage=FALSE, prob=TRUE, start=1L, count=-1L,
    optimize=TRUE, digest=TRUE, parallel=FALSE, verbose=TRUE)
}

\arguments{
    \item{bgen.fn}{the file name of a bgen file}
    \item{verbose}{if TRUE, show information}
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
