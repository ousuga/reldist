#' @name GPW
#' 
#' @title 
#' The Generalized Power Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the generalized power weibull distribution with
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
#' f(x) = alpha*theta*lambda^(-1)*x^(theta-1)*(1+alpha*(x^theta))^(1/lambda-1)*exp(1-(1+alpha*(x^theta))^(1/lambda))
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dGPW} gives the density, \code{pGPW} gives the distribution 
#' function, \code{qGPW} gives the quantile function, and \code{rGPW}
#' generates random deviates.
#' 
#' @export
#' @examples  
#' ## Curve 
#' curve(dGPW(x,0.5,0.5,0.25), from=0, to=2.5, ylim=c(0,3), col="red",ylab="Density")
#' 
#' ## Curve
#' curve(pGPW(x,alpha=0.5,theta=0.5,lambda=0.25), from=0, to=2.5,  col="red")
#' 
#' ## Curve
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qGPW(p,alpha=0.5,theta=0.5,lambda=0.25), y=p)
#' 
#' ## Curve
#' hist(rGPW(10000,alpha=0.5,theta=0.5,lambda=0.25),freq=F,xlab="x", main= "Histogram of rGEP")
#' curve(dGPW(x,alpha=0.5,theta=0.5,lambda=0.25),  from=0, add=T)
#' 
#' ##
#' curve(hGPW(x,alpha=1,theta=2,lambda=1.25), from=0, to=3.2, ylim=c(0,6), col="red",ylab="The hazard function")


dGPW<-function(x,alpha,theta,lambda,log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  loglik <- log(alpha) + log(theta) - log(lambda) + (theta-1)*log(x) +
    ((1/lambda)-1)*log(1+alpha*(x^theta)) + 
    (1-(1+alpha*(x^theta))^(1/lambda))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname GPW
pGPW <- function(q,alpha,theta,lambda, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  cdf <- 1- exp(1-(1 + alpha*(q^theta))^(1/lambda))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf 
}

#' @export
#' @rdname GPW
qGPW <- function(p,alpha,theta,lambda, lower.tail = TRUE, log.p = FALSE){
  
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
  
  term <- 1-log(1-p)
  q <- (((term^lambda)-1) /alpha)^(1/theta)
  q
}

#' @export
#' @rdname GPW
rGPW <- function(n,alpha,theta,lambda){
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
  r <- qGPW(p, alpha,theta,lambda)
  r
}

#' @export
#' @rdname GPW
# Hazard function
hGPW<-function(x,alpha,theta,lambda){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  h <- dGPW(x,alpha,theta,lambda, log = FALSE)/pGPW(q=x,alpha,theta,lambda, lower.tail=FALSE, log.p = FALSE)
  h
}
