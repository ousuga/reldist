#' @name GMW
#' 
#' @title 
#' The Generalized Modified Weibull Distribution 
#' 
#' @description 
#' Density, distribution function, quantile function and
#' random generation for the generalized modified weibull distribution with
#' parameters \code{beta}, \code{theta}, \code{gamma} and \code{lambda}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param beta parameter.
#' @param theta parameter.
#' @param gamma parameter.
#' @param lambda parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The generalized modified weibull with parameters \code{beta}, \code{theta}, 
#' \code{gamma} and \code{lambda} has density given by
#' 
#' f(x)= beta*theta*x^(gamma-1)*(gamma+lambda*x)*exp(lambda*x-beta*x^(gamma)*exp(lambda*x))*
#' (1-exp(beta*x^(gamma)*exp(lambda*x)))^(theta-1)
#' 
#' for x>0.
#' 
#' #' @return 
#' \code{dGMW} gives the density, \code{pGMW} gives the distribution 
#' function, \code{qGMW} gives the quantile function, \code{rGMW}
#' generates random deviatesand and \code{hGMW} gives the hazard function.
#' 
#' @export
#' @examples  
#'## The probability density function
#' curve(dGMW(x,beta=2,theta=0.5,gamma=2,lambda=1.5), from=0, to=0.8, ylim=c(0,6), col="red", las=1, ylab="The probability density function") 
#' 
#' ## The cumulative distribution and the Reliability function
#' curve(pGMW(x,beta=2,theta=0.5,gamma=2,lambda=1.5), from=0, to=1.2, col="red", las=1, ylab="The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(0,0.99999, length.out=100)
#' plot(x=qGMW(p=p,beta=2,theta=0.5,gamma=2,lambda=0.3), y=p, xlab="Quantile", las=1, ylab="Probability")
#' curve(pGMW(x,beta=2,theta=0.5,gamma=2,lambda=0.3),  from=0, add=T)
#' 
#' ## The random function
#' hist(rGMW(1000,beta=2,theta=0.5,gamma=2,lambda=0.3),freq=F,xlab="x", main=""las=1)
#' curve(dGMW(x,beta=2,theta=0.5,gamma=2,lambda=0.3),  from=0, add=T, col="red")
#' 
#' ## The Hazard function
#'curve(hGMW(x,beta=2,theta=1.5,gamma=2,lambda=0.8), from=0, to=1, ylim=c(0,16), col="red",ylab="The hazard function", las=1)
#'

dGMW<-function(y,beta,theta,gamma,lambda, log = FALSE){
  if (any(y<0)) 
    stop(paste("y must be positive", "\n", ""))
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  log_fy<-log(beta*theta) + (gamma-1)*log(y) + log(gamma +lambda*y) +
<<<<<<< HEAD
    lambda*y - beta*(y^gamma)*exp(lambda*y)  +
    (theta-1)*log(1-exp(-beta*(y^gamma)*exp(lambda*y) ))
  
=======
          lambda*y - beta*(y^gamma)*exp(lambda*y)  +
         (theta-1)*log(1-exp(-beta*(y^gamma)*exp(lambda*y) ))
   
>>>>>>> 828ff1df7811bb20fe2f5902a14aa15ad3c483ed
  if (log == FALSE) 
    density<- exp(log_fy)
  else 
    density <- log_fy
  return(density)
}

#' @export
#' @rdname GMW
pGMW <- function(q,beta,theta,gamma,lambda, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  cdf <- (1-exp(-beta*(q^gamma)*exp(lambda*q)))^theta
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}


#' @export
#' @rdname GMW
qGMW <- function(p, beta,theta,gamma,lambda, lower.tail = TRUE, log.p = FALSE) {
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
<<<<<<< HEAD
  
=======

>>>>>>> 828ff1df7811bb20fe2f5902a14aa15ad3c483ed
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,beta,theta,gamma,lambda){
    (1-exp((-beta*(x^gamma))*exp(lambda*x)))^theta
  }
  
  fda1 <- function(x, beta,theta,gamma,lambda, p) {fda(x, beta,theta,gamma,lambda) - p}
  
  r_de_la_funcion <- function(beta,theta,gamma,lambda, p) {
    uniroot(fda1, interval=c(0,1e+06), beta,theta,gamma,lambda, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(beta,theta,gamma,lambda, p)
  q
  
}


#' @export
#' @rdname GMW
rGMW <- function(n,beta,theta,gamma,lambda){
<<<<<<< HEAD
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
=======
if (any(beta<=0 )) 
  stop(paste("beta must be positive", "\n", ""))
if (any(theta<=0)) 
  stop(paste("theta must be positive", "\n", ""))
if (any(gamma<0)) 
  stop(paste("gamma must be positive", "\n", ""))
if (any(lambda<0)) 
  stop(paste("lambda must be positive", "\n", ""))
>>>>>>> 828ff1df7811bb20fe2f5902a14aa15ad3c483ed
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qGMW(p, beta,theta,gamma,lambda)
  r
}


#' @export
#' @rdname GMW
hGMW<-function(x,beta,theta,gamma,lambda, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(beta<=0 )) 
    stop(paste("beta must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(gamma<0)) 
    stop(paste("gamma must be positive", "\n", ""))
  if (any(lambda<0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  h <- dGMW(x,beta,theta,gamma,lambda, log = FALSE)/pGMW(q=x,beta,theta,gamma,lambda, lower.tail=FALSE, log.p = FALSE)
  h
}



