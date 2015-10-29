#' @name PL
#' 
#' @title 
#' The Power Lindley Distribution
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the power Lindley distribution with
#' parameters \code{alpha} and \code{beta}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha shape parameter.
#' @param beta scale parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The power Lindley distribution with parameters \code{alpha} and
#' \code{beta} has density given by
#' 
#' f(x) = (((alpha*beta)^2)/(beta+1))*(1+x^(alpha))*x^(alpha-1)*exp(-beta*(x^alpha))
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dPL} gives the density, \code{pPL} gives the distribution 
#' function, \code{qPL} gives the quantile function, and \code{rPL}
#' generates random deviates.
#' 
#' @export
#' @examples  
#' ## Curve 
#' curve(dPL(x,1,0.5), from=0, to=15, ylim=c(0,0.6), col="red",ylab="Density")
#' 
#' ## Curve
#' curve(pPL(x,alpha=1,beta=0.5), from=0, to=15,  ylim=c(0,1), col="blue")
#' 
#' ## Curve
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qPL(p=p,alpha=0.8,beta=0.3), y=p)
#' 
#' ## Curve
#' hist(rPL(1000,alpha=2,beta=0.5),freq=F,xlab="x", main= "Histogram of rGEP")
#' curve(dPL(x,alpha=2,beta=0.5),  from=0, add=T)

dPL<-function(x,alpha,beta, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))  
  
  loglik<- log(alpha) + 2*log(beta) - log(beta+1) +
    log(1+(x^alpha)) + (alpha-1)*log(x) - beta*(x^alpha)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname PL
pPL <- function(q,alpha,beta, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))  
  
  cdf <- 1-(1+((beta/(beta+1))*q^alpha))*exp(-beta*(q^alpha))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname PL
qPL <- function(p,alpha,beta, lower.tail = TRUE, log.p = FALSE){
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))  
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,alpha, beta){
    1-(1+((beta/(beta+1))*x^alpha))*exp(-beta*(x^alpha))
  }
  
  fda1 <- function(x, alpha, beta, p) {fda(x, alpha, beta) - p}
  
  r_de_la_funcion <- function(alpha, beta, p) {
    uniroot(fda1, interval=c(0,1e+06), alpha, beta, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(alpha, beta, p)
  q
}


#' @export
#' @rdname PL
rPL <- function(n,alpha,beta){
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(beta<=0)) 
    stop(paste("beta must be positive", "\n", ""))  
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qPL(p, alpha,beta)
  r
}