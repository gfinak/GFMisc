\name{asinh_trans}
\alias{asinh_trans}
\title{Arcsinh Transformation for ggplot}
\usage{
asinh_trans(c)
}
\arguments{
  \item{x}{\code{numeric} cofactor}
}
\value{
a transformation to be used with \code{coord_trans}
}
\description{
Transform using the arcsinh with a cofactor.
}
\examples{
\dontrun{
 qplot(foo)+coord_trans(ylim=asinh_trans(1))
}
}
\seealso{
\link{coord_trans}
}

