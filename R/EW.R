#' @name EW
#' 
#' @title 
#' The Exponentiated Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the exponentiated weibull  distribution with
#' parameters \code{alpha}, \code{theta} and \code{lambda}.
#' 
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha scale parameter.
#' @param theta,lambda shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Exponentiated Weibull Distribution with parameters \code{alpha}, 
#' \code{theta} and \code{lambda} has density given by
#' 
#' f(x)=lambda*alpha*theta*x^(theta-1)*exp(-alpha*(x^theta))*(1-exp(-alpha*(x^theta)))^(lambda-1)
#' 
#' 
#' @return Sabes una cosa...
#' @export
#' @examples  
#' ## Curve 
#' curve(dEW(x,2,1.5,0.5), from=0, to=2, col="red",ylab="Density") 
#' 
#' ## Curve
#' curve(pEW(x,alpha=2,theta=1.5,lambda=0.5), from=0, to=2,  col="red")
#' 
#' ## Curve
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qEW(p,alpha=2,theta=1.5,lambda=0.5), y=p)
#' 
#' ## Curve
#' hist(rEW(10000,alpha=2,theta=1.5,lambda=0.5),freq=F,xlab="x", main="")
dEW<-function(x,alpha,theta,lambda, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  loglik<- log(lambda) +log(alpha) + log(theta) + (theta-1)*log(x)-
    alpha*(x^theta) + (lambda-1)*log(1-exp(-alpha*(x^theta)))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname EW
pEW <- function(q,alpha,theta,lambda, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  cdf <- (1-exp(-alpha*(q^theta)))^lambda
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname EW
qEW <- function(p,alpha,theta,lambda, lower.tail = TRUE, log.p = FALSE){
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
  
  q <- ((-1/alpha)*log(1-p^(1/lambda)))^(1/theta)
  q
}

#' @export
#' @rdname EW
rEW <- function(n,alpha,theta, lambda){
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
  r <- qEW(p, alpha,theta,lambda)
  r
}
