#' @name IW
#' 
#' @title 
#' The Inverse Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the inverse weibull distribution with
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
#' The inverse weibull distribution with parameters \code{alpha} and
#' \code{theta} has density given by
#' 
#' f(x) = alpha*theta*x^(-theta-1)*exp(-alpha*(x^-theta))
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dIW} gives the density, \code{pIW} gives the distribution 
#' function, \code{qIW} gives the quantile function, and \code{rIW}
#' generates random deviates.
#' 
#' @export
#' @examples  
#' ## Curve 
#' curve(dIW(x,alpha=2,   theta=1), from=0, to=10, ylim=c(0,0.6), col="red",ylab="Density")
#' 
#' ## Curve
#' curve(pIW(x,alpha=2,theta=1), from=0, to=10, ylim=c(0,1), col="blue",ylab="Density")
#' 
#' ## Curve
#' p <- seq(0,0.998, length.out=100)
#' plot(x=qIW(p,alpha=5,theta=2.5), y=p)
#' 
#' ## Curve
#' hist(rIW(10000,alpha=5,theta=2.5),freq=F,xlab="x", main= "Histogram of rIW")
#' curve(dIW(x,alpha=5,theta=2.5),  from=0, add=T)

dIW<-function(x,alpha,theta, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha <= 0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  loglik<- log(alpha) + log(theta) - (theta+1)*log(x) - 
    alpha*(x^-theta)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)  
}

#' @export
#' @rdname IW
pIW <- function(q,alpha,theta, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  cdf <- exp((-alpha)*(q^(-theta)))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname IW
qIW <- function(p,alpha,theta, lower.tail = TRUE, log.p = FALSE){
  if (any(alpha<=0 )) 
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
  
  q <- ((-1/alpha)*log(p))^(-1/theta)
  q
}

#' @export
#' @rdname IW
rIW <- function(n,alpha,theta){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qIW(p, alpha,theta)
  r
}

