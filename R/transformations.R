#'Optimize over the parameter of the arcsine transformation
#'Inverse of hyperbolic arcsine
#'@importFrom stats optimize
#'@param data \code{numeric} vector of the data 
#'@param upper \code{numeric} boundary of the interval
#'@param lower \code{numeric} boundary of the interval
#'@return Optimal cofactor
#'@examples
#'
#'K<-100
#'sigma<-2
#'mu<-0
#'data<-sinh(rnorm(10000,mean=mu,sd=sigma))*1/K
#'K.hat<-optimAsinhCofactor(data=data)
#'dtrans<-asinh(data*K.hat)
#'mu.hat<-mean(dtrans)
#'sigma.hat<-sd(dtrans)
#'sprintf("K=%s, K.hat=%s",K,signif(K.hat,4))
#'sprintf("mu=%s, mu.hat=%s",mu,signif(mu.hat,4))
#'sprintf("sigma=%s, sigma.hat=%s",sigma,signif(sigma.hat,4))
#'boxplot(data.frame(data,dtrans),names=c("Raw","Transformed"))
optimAsinhCofactor<-function(data,upper=500000,lower=0.0001){
  .asinh_c <- function(x=data, c){  # Inverse IHS transformation
    asinh(x*c)
  }
  
  .asinh.ll<-function(c,x=data){
    xtrans<-.asinh_c(x=x,c=c)
    sum(-log(dnorm((xtrans-mean(xtrans)),mean=0,sd=sd(xtrans))))-sum(log(abs((c/sqrt(1+x^2*c^2)))))
  }
  
  optimize(f=.asinh.ll, lower=lower, upper=upper, x=data)$minimum 
}



#' Arcsinh Transformation for ggplot
#' 
#' Transform using the arcsinh with a cofactor.
#' @param x \code{numeric} cofactor
#' @return a transformation to be used with \code{coord_trans}
#' @seealso \link{coord_trans}
#' @examples
#'\dontrun{
#'  qplot(foo)+coord_trans(ylim=asinh_trans(1))
#'}
#' @export
asinh_trans <- function(c){
  trans_new(name = 'asinh', transform = function(x) asinh(x*c), 
            inverse = function(x) sinh(x)/c)
}
