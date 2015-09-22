#' Predict the value of the function at the given x-values
#' @param x Vector of x-values
#' @param spline An object produced by fiducialSpline
E.spline <- function(x, spline){
  degree <- length(spline$alpha) - length(spline$knots) - 1
  X <- tp(x, knots=spline$knots, degree=degree)
  Ey <- X %*% spline$alpha
  return(Ey)
}


#' Sample from fiducial b-spline distribution
#' 
#' @param x Vector of x-values
#' @param y Vector of y-values
#' @param degree Degree of the b-spline to fit to the data
#' @param num.knots A scalar or vector of the number of knot points to be considered. e.g. 3:5 or just 3.
#' @param sigma.jitter MCMC tuning parameter: the size of jumps for the knot point locations.
#' @param sigma.split MCMC tuning parameter: how far apart knots end up when we split one.
#' @param chain.length How long a MCMC to create.
#' @param start.point A list of starting points 
#' @param burn.in The length of chain to be considered 'burn in'.
#' @param num.chains How many chains to create.
#' @param num.jacobians How many randomly selected jacobians to average for each step of the MCMC.
#' @examples 
#' require(SemiPar)
#' data(lidar)
#' x <- lidar$range
#' y <- lidar$logratio
#' plot(x, y)

#' foo <- fiducial.spline(x, y, 3, 2)
#' plot(foo)  # plot the three chains...
#' confint(foo)
#' @export
fiducial.spline <- function(x, y, degree, num.knots, sigma.jitter=NULL, sigma.split=NULL,
							chain.length=500, start.point=NULL,
							burn.in=100, num.chains=3, num.jacobians=100 ){

	possible.num.knots <- seq( min(num.knots), max(num.knots) )
	
	if(is.null(start.point)){
		k <- floor(median(possible.num.knots))
		knots <- quantile( x, seq(0,1, length=k+2) )
		names(knots) <- NULL
		knots <- knots[-c(1,length(knots))]
		temp <- calc.start.point(x, y, degree, knots)
		start.point <- NULL
		for(i in 1:num.chains){
			start.point[[i]] <- temp                                    
		}		
	}
	if( is.null(sigma.jitter) ){
		x.range <- range(x)
		sigma.jitter <- (x.range[2] - x.range[1]) / 10
	}
	if( is.null(sigma.split) ){
		x.range <- range(x)
		sigma.split <- (x.range[2] - x.range[1]) / 5
	}
	
	out <- NULL
	out$x <- x
	out$y <- y
	out$degree <- degree
	out$num.knots <- possible.num.knots
	out$sigma.jitter <- sigma.jitter
	out$sigma.split  <- sigma.split
	out$num.chains <- num.chains
	out$burn.in <- 0
	out$num.jacobians <- num.jacobians

	out$chains <- NULL
  out$chain.num.knots <- NULL
	# Initial chains
	for(i in 1:num.chains){
		out$chains[[i]] <- calculate.MCMC.chain( x=x, y=y, 
					possible.num.knots=possible.num.knots, 
          sigma.jitter=sigma.jitter, sigma.split=sigma.split, 
					n.iter = burn.in + chain.length, 
					start.point=start.point[[i]], num.jacobians=num.jacobians )
    out$chain.num.knots[[i]] <- sapply( out$chains[[1]], function(x){length(x$knots)})
	}
	class(out) <- 'fiducialSpline'
	return(out)
}


#' Plot the fiducial spline
#' @param obj A fiducialSpline object created using fiducialSpline 
#' @param grid.length The number of points along the x-axis to predict at. 
#' @param xlab The label of the x-axis
#' @export
plot.fiducialSpline <- function(obj, grid.length=401, xlab=''){
	x.grid <- seq(min(obj$x), max(obj$x), length=grid.length)
	n <- length(obj$chains[[1]])
	n.chains <- length(obj$chains)
	old.mfrow <- par()$mfrow 
	par(mfrow=c(n.chains, 1))
	for(i in 1:n.chains){
		plot(obj$x, obj$y, xlab=xlab)
		for(j in 2:n){
		  lines(x.grid, E.spline(x.grid, obj$chains[[i]][[j]]), col=1)
		}
    points(obj$x, obj$y, col='red')    
	}
	par(mfrow=old.mfrow)		
}

#' Calculate the fiducial confidence intervals (credible intervals???)
#' @param obj A fiducialSpline object created using fiducialSpline 
#' @param level The level of the interval
#' @export
confint.fiducialSpline <- function(obj, level=0.95){
	probs <- c( (1-level)/2, 1-(1-level)/2 )
	out <- list()
  n <- length(obj$chains[[1]])
  num.knots <- obj$num.knots
	for( i in 1:length(num.knots) ){
    temp <- NULL
    # splice the chains together
    for( j in 1:length(obj$chains) ){
      temp <- c(temp, obj$chains[[j]][which( sapply(obj$chains[[j]], function(x){length(x$knots)}) == num.knots[i] ) ] ) 
    }
    # remove the z element and convert from list of lists to a matrix
    temp <- t(as.matrix( sapply(temp, function(x){out <- x; out$z <- NULL; return(unlist(out))}) ))  
    
    if( dim(temp)[1] > 50 ){
			out[[i]] <- t(apply(temp, 2, quantile, probs=probs))
		}else{
			out[[i]] <- NA
		}
	}
  if( length(num.knots) == 1 ){
    out <- out[[1]]
  }
	return(out)
}


#' Calculate the MCMC acceptance rate
#' @param obj A fiducialSpline object created using fiducialSpline 
#' @export
acceptance.rate <- function(obj){
  n <- length(obj$chains[[1]])
  num.chains <- obj$num.chains
  count <- 0
  for( i in 1:num.chains ){
    for( j in 2:n ){
      if( !identical( obj$chains[[i]][[j]], obj$chains[[i]][[j-1]] ) ){
        count <- count + 1
      }
    }
  }
  return( count / (num.chains * n) )  
}

#' Trim the beginning of the MCMC chains
#' @param obj A fiducialSpline object created using fiducialSpline 
#' @param amount Number of initial steps to remove
#' @export
trim <- function(obj, amount){
  out <- obj
  n <- length(out$chains[[1]])
  if( n > amount){
    for( i in 1:length(out$chains) ){
      out$chains[[i]] <- out$chains[[i]][-c(1:amount)]
    }
  }else{
    error('Amount to trim is more than the length of the chain')
  }
  return(out)
}

#' Thin the MCMC chains
#' @param obj A fiducialSpline object created using fiducialSpline 
#' @param by The thinning interval
#' @export
thin <- function(obj, by=10){
  out <- obj
  n <- length(out$chains[[1]])
  for( i in 1:length(out$chains) ){
    out$chains[[i]] <- out$chains[[i]][ seq(1,n,by=by) ]
  }
  return(out)
}

#' Convert a fiducialSpline object to a mcmc.list
#' @param obj A fiducialSpline object created using fiducialSpline 
#' @export
as.mcmc.fiducialSpline <- function(obj){
  if( length(obj$num.knots) > 1){
    return( as.mcmc.fiducialSpline.modelSelection(obj) )
  }else{
    return( as.mcmc.fiducialSpline.knownNumberKnots(obj) )
  }
}

#' Convert a fiducialSpline object to a mcmc.list
#' @param obj A fiducialSpline object created using fiducialSpline 
#' @export
as.mcmc.fiducialSpline.modelSelection <- function(obj){
  out <- list()
  n <- length(obj$chains[[1]])
  num.knots <- obj$num.knots
  
  out <- list()
  for( i in 1:length(num.knots) ){
    temp1 <- NULL
    for( j in 1:length(obj$chains) ){
      temp2 <- obj$chains[[j]][which( sapply(obj$chains[[j]], 
                                            function(x){length(x$knots)}) == num.knots[i] ) ]  
      if( length(temp2) > 0 ){
        # remove the z element and convert from list of lists to a matrix
        temp2 <- t(as.matrix( sapply(temp2, function(x){out <- x; out$z <- NULL; return(unlist(out))}) ))  
        temp1 <- rbind(temp1, temp2)
      }
    }
    if( !is.null(temp1) ){
      out[[i]] <- as.mcmc(temp1)
    }else{
      out[[i]] <- NULL
    }
  }
  
  if(length(out) == 1){
    out <- out[[1]]
  }
  return(out)
}

#' Convert a fiducialSpline object to a mcmc.list
#' @param obj A fiducialSpline object created using fiducialSpline 
#' @export
as.mcmc.fiducialSpline.knownNumberKnots <- function(obj){
  out <- list()
  n <- length(obj$chains[[1]])
  num.knots <- obj$num.knots
  for( i in 1:length(num.knots) ){
    out[[i]] <- list()
    # splice the chains together
    for( j in 1:length(obj$chains) ){
      # remove the z element and convert from list of lists to a matrix
      temp <- t(as.matrix( sapply(obj$chains[[j]], function(x){out <- x; out$z <- NULL; return(unlist(out))}) ))  
      out[[i]][[j]] <- as.mcmc(temp)
    }
    out[[i]] <- as.mcmc.list(out[[i]])
  }
  if(length(out) == 1){
    out <- out[[1]]
  }
  return(out)
}



#' Calculate a reasonable starting point for the fiducial MCMC free-knot b-spline
#' @param x A vector of x-values
#' @param y A vector of y-values
#' @param degree The degree of b-spline
#' @param knots A vector of initial knot points
#' @export
calc.start.point <- function(x, y, degree, knots){
  temp <- NULL
  temp$knots <- knots
  
  X <- tp(x, knots=temp$knots, degree=degree)
  XtXinv <- chol2inv(chol(crossprod(X))) 
  beta.hat <- XtXinv %*% crossprod(X, y)
  sqrtXtXinv <- sqrt.matrix(XtXinv)
  df <- length(x) - length(knots) - length(beta.hat)
  e <- X %*% beta.hat - y
  s2 <- sum(e^2)/df
  temp$var <- s2
  temp$z <- rnorm(length(beta.hat))
  temp$alpha <- as.vector(beta.hat + sqrtXtXinv %*% temp$z)
  
  return(temp)
}	
