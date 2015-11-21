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
#' function, \code{qGMW} gives the quantile function, and \code{rGMW}
#' generates random deviates.
#' 
#' @export
#' @examples  
#' 
#' 