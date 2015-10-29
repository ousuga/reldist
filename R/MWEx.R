#' @name MWEx
#' 
#' @title 
#' The Modified Weibull Extension Distribution 
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the modified weibull extension distribution with
#' parameters \code{alpha}, \code{beta} and \code{lambda}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter.    
#' @param beta parameter.
#' @param lambda parameter.        
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The modified weibull extension distribution with parameters \code{alpha}, \code{beta}
#' and \code{lambda} has density given by
#' 
#' f(x) = 
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dMWEx} gives the density, \code{pMWEx} gives the distribution 
#' function, \code{qMWEx} gives the quantile function, and \code{rMWEx}
#' generates random deviates.
#' 
#' @export
#' @examples  
#' ## Curve 
#' curve(dMWEx(x,alpha=1/0.5,beta=3,lambda=2),from=0, to=2.5, ylim=c(0,2.5), col="red", ylab="Density")
#' 
#' ## Curve
#' curve(pMWEx(x,alpha=1/0.5,beta=3,lambda=2),from=0, to=2.5, col="red")
#' 
#' ## Curve
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qMWEx(p,alpha=1/0.7,beta=2,lambda=0.2), y=p)
#' 
#' ## Curve
#' hist(rMWEx(10000,alpha=1/0.7,beta=2,lambda=0.2),freq=F,xlab="x", main= "Histogram of rGEP")
#' curve(dMWEx(x,alpha=1/0.7,beta=2,lambda=0.2),  from=0, add=T)

dMWEx<-function(x,alpha,beta,lambda, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  
  loglik<-log(lambda*beta) + (beta-1)*log(x/alpha) +
    (x/alpha)^beta + lambda*alpha*(1-exp((x/alpha)^beta))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname MWEx

pMWEx<- function(q,alpha,beta,lambda, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  cdf <- 1- exp(lambda*alpha*(1-exp((q/alpha)^beta)))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname MWEx
qMWEx <- function(p,alpha,beta,lambda, lower.tail = TRUE, log.p = FALSE){
  
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
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
  
  q <- alpha*(log(1-((1/(lambda*alpha))*log(1-p))))^(1/beta)
  q
}


#' @export
#' @rdname MWEx
rMWEx <- function(n,alpha,beta,lambda){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qMWEx(p, alpha,beta,lambda)
  r
}

