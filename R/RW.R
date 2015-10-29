#' @name RW
#' 
#' @title 
#' The Reflected Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the reflected weibull distribution with
#' parameters \code{alpha} and \code{theta}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha         ##############BUSCAR############
#' @param theta         ##############BUSCAR############   
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#'
#' @details 
#' The reflected weibull distribution with parameters \code{alpha} and
#' \code{theta} has density given by
#' 
#' f(x) = alpha*theta*(-x)^(theta-1)*exp(-alpha*(-x)^alpha)
#' 
#' for - inf < x < 0.  ### También x < 0.
#' 
#' @return 
#' \code{dRW} gives the density, \code{pRW} gives the distribution 
#' function, \code{qRW} gives the quantile function, and \code{rRW}
#' generates random deviates.
#' 
#' @export
#' @examples  
#' ## Curve 
#' curve(dRW(x,1,1), from=-5, to=0, ylim=c(0,1), col="red",ylab="Density")
#' 
#' ## Curve
#' curve(pRW(x,alpha=1,theta=1), from=-5, to=0, ylim=c(0,1), col="blue",ylab="Density")
#' 
#' ## Curve
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qRW(p=p,alpha=0.5,theta=2), y=p)
#' 
#' ## Curve
#' hist(rRW(10000,alpha=0.5,theta=2), freq=F,xlab="x", main= "Histogram of rRW")
#' curve(dRW(x,alpha=0.5,theta=2),  from=-5, to=0, add=T)

dRW<-function(x,alpha,theta, log=FALSE){
  if (any(x>0)) 
    stop(paste("x must be negative", "\n", ""))
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  loglik<- log(alpha) + log(theta) + (theta-1)*log(-x) -
    alpha*((-x)^theta)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname RW
pRW <- function(q,alpha,theta, lower.tail=TRUE, log.p = FALSE){
  # if (any(q<0)) 
  #  stop(paste("q must be positive", "\n", ""))
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  cdf <- exp(-alpha*(-q)^theta)
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname RW
qRW <- function(p,alpha,theta, lower.tail = TRUE, log.p = FALSE){
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- -{((-1/alpha)*log(p))^(1/theta)}
  q
}

#' @export
#' @rdname RW
rRW <- function(n,alpha,theta){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qRW(p,alpha,theta)
  r
}

