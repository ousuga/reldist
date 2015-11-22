#' @name KW
#' 
#' @title 
#' The Kumaraswamy Weibull Distribution 
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the kumaraswamy weibull distribution with
#' parameters \code{alpha}, \code{beta}, \code{a} and \code{b}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter.    
#' @param beta parameter.
#' @param a parameter.
#' @param b parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The kumaraswamy weibull with parameters \code{alpha}, \code{beta}, \code{a} and \code{b}
#' has density given by
#' 
#' f(x)=a*b*alpha*theta*x^(theta-1)*exp(-alpha*x^theta)*(1-exp(-alpha*x^theta))^(a-1)*
#' (1-(1-exp(-aplha*x^theta))^a)^(b-1)
#' 
#' for x>0.
#' 
#' @return 
#' \code{dKW} gives the density, \code{pKW} gives the distribution 
#' function, \code{qKW} gives the quantile function, and \code{rKW}
#' generates random deviates.
#' 
#' @export
#' @examples 
#' ## Curve 
#' curve(dKW(x,alpha=3,beta=0.8,a=2,b=1.5), from=0, to=3, ylim=c(0,2.5), col="red",ylab="Density")
#' 
#' ## Curve 
#' curve(dKW(x,alpha=3,theta=0.8,a=2,b=1.5), from=0, to=3, ylim=c(0,2.5), col="red",ylab="Density")
#' 
#' ## Curve 
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qKW(p,alpha=3,theta=0.8,a=2,b=1.5), y=p)
#' 
#' ## Curve
#' hist(rKW(10000,alpha=1.2,theta=1.5,a=2,b=1.2),freq=F,xlab="x", main= "Histogram of rKW")
#' curve(dKW(x,alpha=1.2,theta=1.5,a=2,b=1.2),  from=0, add=T)
#' 
#' ## Curve
#' curve(hKW(x,alpha=3,theta=0.5,a=3,b=3), from=0, to=2.1, ylim=c(0,10), col="red",ylab="The hazard function")

dKW<-function(x,alpha,theta,a,b, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  term <- -alpha*(x^theta)
  loglik<- log(a*b*alpha*theta) + (theta-1)*log(x) + term + 
    (a-1)*log(1-exp(term)) + (b-1)*log(1-(1-exp(term))^a)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density) 
}
#' @export
#' @rdname KW
pKW <- function(q,alpha,theta,a,b, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  cdf <- 1-(1-(1-exp(-alpha*(q^theta)))^a)^b
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
  
}
#' @export
#' @rdname KW
qKW <- function(p,alpha,theta,a,b, lower.tail = TRUE, log.p = FALSE){
  
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- ((-1/alpha)*(log(1-(1-(1-p)^(1/b))^(1/a))))^(1/theta)
  q
}
#' @export
#' @rdname KW
rKW <- function(n,alpha,theta,a,b){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qKW(p,alpha,theta,a,b)
  r
}
#' @export
#' @rdname KW
hKW<-function(x,alpha,theta,a,b){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  h <- dKW(x,alpha,theta,a,b, log = FALSE)/pKW(q=x,alpha,theta,a,b, lower.tail=FALSE, log.p = FALSE)
  h
}


