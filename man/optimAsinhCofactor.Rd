\name{optimAsinhCofactor}
\alias{optimAsinhCofactor}
\title{Optimize over the parameter of the arcsine transformation
Inverse of hyperbolic arcsine}
\usage{
optimAsinhCofactor(data, upper = 5e+05, lower = 1e-04)
}
\arguments{
  \item{data}{\code{numeric} vector of the data}

  \item{upper}{\code{numeric} boundary of the interval}

  \item{lower}{\code{numeric} boundary of the interval}
}
\value{
Optimal cofactor
}
\description{
Optimize over the parameter of the arcsine transformation
Inverse of hyperbolic arcsine
}
\examples{
K<-100
sigma<-2
mu<-0
data<-sinh(rnorm(10000,mean=mu,sd=sigma))*1/K
K.hat<-optimAsinhCofactor(data=data)
dtrans<-asinh(data*K.hat)
mu.hat<-mean(dtrans)
sigma.hat<-sd(dtrans)
sprintf("K=\%s, K.hat=\%s",K,signif(K.hat,4))
sprintf("mu=\%s, mu.hat=\%s",mu,signif(mu.hat,4))
sprintf("sigma=\%s, sigma.hat=\%s",sigma,signif(sigma.hat,4))
boxplot(data.frame(data,dtrans),names=c("Raw","Transformed"))
}

