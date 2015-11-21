#' @name OW
#' 
#' @title 
#' The Odd Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the odd weibull distribution with
#' parameters \code{alpha}, \code{theta} and \code{lambda}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter.    
#' @param theta parameter.
#' @param lambda parameter.        
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The generalized power weibull with parameters \code{alpha}, \code{theta}
#' and \code{lambda} has density given by
#' 
#' f(x) = alpha*theta*lambda*x^(theta-1)*exp(alpha*(x^theta))*(exp(alpha*(x^theta))*-1)^(lambda-1)*(1+(exp(alpha*(x^theta))-1)^lambda)^-2
#' ### revisar
#' for x > 0.
#' 
#' @return 
#' \code{dOW} gives the density, \code{pOW} gives the distribution 
#' function, \code{qOW} gives the quantile function, and \code{rOW}
#' generates random deviates.
#' 
#' @export
#' @examples  
#' ## Curve 
#' curve(dOW(x,2,3,0.2), from=0, to=4, ylim=c(0,2), col="blue",ylab="Density")
#' 
#' ## Curve
#' curve(pOW(x, alpha=2,theta=3,lambda=0.2), from=0, to=4, col="blue")
#' 
#' ## Curve
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qOW(p,alpha=2,theta=3,lambda=0.2), y=p)
#' 
#' ## Curve
#' hist(rOW(10000,alpha=2,theta=3,lambda=0.2),freq=F,xlab="x", main= "Histogram of rOW")
#' curve(dOW(x,alpha=2,theta=3,lambda=0.2),  from=0, add=T)
#' 
#' ##
#' curve(hOW(x,alpha=5,theta=3,lambda=0.7), from=0, to=2.5, ylim=c(0,3), col="red",ylab="The hazard function")


dOW<-function(x,alpha,theta,lambda, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  loglik<- log(alpha) +log(theta) + log(lambda) + (theta-1)*log(x) +
    alpha*(x^theta) + (lambda-1)*log(exp(alpha*(x^theta))-1) -
    2*log(1+(exp(alpha*(x^theta))-1)^lambda)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname OW
pOW <- function(q,alpha,theta,lambda, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  cdf <- 1 - (1 + (exp(alpha*(q^theta))-1)^lambda )^(-1)
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname OW
qOW <- function(p,alpha,theta,lambda, lower.tail = TRUE, log.p = FALSE){
  
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- {(1/alpha)*log(1+(-1 + (1-p)^(-1))^(1/lambda))}^(1/theta)
  q
}

#' @export
#' @rdname OW
rOW <- function(n,alpha,theta,lambda){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qOW(p,alpha,theta,lambda)
  r
}
#' @export
#' @rdname OW
hOW<-function(x,alpha,theta,lambda){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  h <- dOW(x,alpha,theta,lambda, log = FALSE)/pOW(q=x,alpha,theta,lambda, lower.tail=FALSE, log.p = FALSE)
  h
}

