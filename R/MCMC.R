##########################################################
calculate.MCMC.chain <- function(x, y, possible.num.knots, sigma.jitter, sigma.split, n.iter, 
										start.point, num.jacobians=100){
	temp <- range(x)
	data <- list(x=x, y=y, min.x=temp[1], max.x=temp[2], n=length(x) )
	
	order <- length(start.point$alpha) - length(start.point$knots)
	degree <- order - 1

  chain <- list()
  chain[[1]] <- start.point
		
	curr <- start.point
	
	fid.from.density  <- unscaled.fiducial.density(data, curr$alpha, curr$knots, curr$var, 
												degree, num.jacobians=num.jacobians, log=TRUE);

	for( i in 2:n.iter ){
		jump <- try( rJump(curr, degree, data, sigma.jitter, sigma.split, possible.num.knots), silent=TRUE )
# 		jump <- rJump(curr, degree, data, sigma.jitter, sigma.split, possible.num.knots)
		if( class(jump) != 'try-error' ){
			proposed <- jump$to
			jump.to.density   <- jump$ln.density.forward       # These are log densities!
			jump.from.density <- jump$ln.density.backward
												          
			fid.to.density    <- try(unscaled.fiducial.density(data, proposed$alpha, proposed$knots, 
											proposed$var, degree, num.jacobians=num.jacobians, log=TRUE), silent=TRUE)
			if(class(fid.to.density) == 'try-error'){
				fid.to.density <- 0
				r <- 0
			}else{											
        r <- exp( fid.to.density - fid.from.density + jump.from.density - jump.to.density);
			}
		}else{
			r <- 0
		}
		
		if(is.na(r)){					      # the to.density and from.density are both 0.  
			#print('In the tails!')   # ie we are too far out in the tail
			r <- 0							      # There isn't much we can do but to propose a new value.  
		}			

		if( r > runif(1)){ 
      chain[[i]] <- proposed
			fid.from.density <- fid.to.density
			curr <- proposed
		}else{
      chain[[i]] <- curr
		}
	}

	return(chain);   	
}
	
# models[[ m ]][[ ]]  -> models[[m]][]
unlist.models <- function(obj, k, order){ 
	if( length(obj) >= 1 ){
		out <- t(sapply(obj, function(x){unlist(x[ c('var', 'alpha', 'knots') ])}))
	}else{
		out <- array(NA, dim=c(0, 1+k+(order+k)))
	}	
	colnames(out) <- c('var', paste('alpha', 1:(order+k), sep=''), paste('knot',1:k,sep=''))
	return(out)
}


	
rJump <- function(from, degree, data, sigma.jitter, sigma.split, n.knots){
	jump <- 1:5
	# jitter, split, create, merge, delete,  
	jump.weights     <- c(1/2, 1/8, 1/8, 1/8, 1/8)
# 	jump.weights     <- c(1, 0, 0, 0, 0)
	
	rev.jump.weights <- c(1/2, 1/8, 1/8, 1/8, 1/8)
# 	rev.jump.weights <- c(1, 0, 0, 0, 0)
	
	Green.jacobian <- c(1, 2, 2, 1/2, 1/2)
	k <- length(from$knots)
	if(k == max(n.knots)){
		jump.weights[c(2,3)] <- 0
	}
	if(k == min(n.knots)){
		jump.weights[c(4,5)] <- 0
	}	
				
	which.jump <- sample(jump, size=1, prob=jump.weights)
	if(which.jump == 1){
		out <- rJump.SameNumKnots(from, degree, data, sigma.jitter)
	}else if(which.jump == 2){
		out <- rJump.splitKnot(from, degree, data, sigma.split)
	}else if(which.jump == 3){
		out <- rJump.createKnot(from, degree, data, sigma.split)
	}else if(which.jump == 4){
		out <- rJump.mergeKnot(from, degree, data, sigma.split)
	}else{
		out <- rJump.deleteKnot(from, degree, data, sigma.split)
	}
	out$ln.density.forward  <- out$ln.density.forward  + 
                              log(    jump.weights[which.jump]) + 
                              log(Green.jacobian[which.jump])
	out$ln.density.backward <- out$ln.density.backward + log(rev.jump.weights[which.jump])
	return(out)
}

rJump.deleteKnot <- function(from, degree, data, sd.knot){
	to <- from
	i <- sample(1:length(to$knots), size=1)
	knot <- to$knots[i]
	u <- knot-data$min.x / (data$max.x - data$min.x)
	to$knots <- to$knots[-i]
	z <- to$z[1]
	to$z <- to$z[-1]

  # 	X <- bs(data$x, knots=to$knots, degree=degree, intercept=TRUE)
	X <- tp(data$x, knots=to$knots, degree=degree)
	XtXinv <- chol2inv(chol(crossprod(X))) 
	alpha.hat <- XtXinv %*% crossprod(X, data$y)
	df <- length(data$x) - length(to$knots) - length(alpha.hat)
	y.hat <- X %*% alpha.hat
	e2 <- (data$y - y.hat)^2
	s2 <- sum(e2)/df 
	sqrtXtXinv <- sqrt.matrix(XtXinv)
	to$alpha <- as.vector( alpha.hat + sqrt(s2) * sqrtXtXinv %*% to$z )

	ln.density.bigger  <- dunif(u, log=TRUE) + dnorm(z, log=TRUE) 
	ln.density.smaller <- log( 1/length(from$knots) )
	out <- list(to=to, ln.density.forward=ln.density.smaller, ln.density.backward=ln.density.bigger)
	return(out)
}

rJump.createKnot <- function(from, degree, data, sd.knot){
	u <- runif(1)
	z <- rnorm(1)
	to <- from
	to$knots <- sort(c(to$knots, data$min.x + u*(data$max.x - data$min.x)))
	to$z <- c(to$z, z)
# 	X <- bs(data$x, knots=to$knots, degree=degree, intercept=TRUE)
	X <- tp(data$x, knots=to$knots, degree=degree)
	XtXinv <- chol2inv(chol(crossprod(X))) 
	alpha.hat <- XtXinv %*% crossprod(X, data$y)
	df <- length(data$x) - length(to$knots) - length(alpha.hat)
	y.hat <- X %*% alpha.hat
	e2 <- (data$y - y.hat)^2
	s2 <- sum(e2)/df 
	sqrtXtXinv <- sqrt.matrix(XtXinv)
	to$alpha <- as.vector( alpha.hat + sqrt(s2) * sqrtXtXinv %*% to$z )

	ln.density.bigger  <- dunif(u, log=TRUE) + dnorm(z, log=TRUE) 
	ln.density.smaller <- log( 1/length(to$knots) )
	out <- list(to=to, ln.density.forward=ln.density.bigger, ln.density.backward=ln.density.smaller)
	return(out)
}
	
	
rJump.mergeKnot <- function(from, degree, data, sd.knot){
	index <- sample(1:(length(from$knots)-1), size=1)
	knot <- from$knots[index]
	knot.range <- min(data$max.x - knot, knot-data$min.x) * .9
	temp <- rJump.mergeKnot.determanistic(from, degree, data, index)
	to <- temp$to	
	ln.density.forward <- dTruncatedNorm(temp$u, mean=0, sd=sd.knot, 
									     min=-knot.range, max=knot.range, log=TRUE) +
									 dnorm(temp$z, log=TRUE)
	out <- list(to=to, ln.density.forward=ln.density.forward, ln.density.backward=0)
	return(out)
}

rJump.mergeKnot.determanistic <- function(from, degree, data, index){
	knot <- mean(from$knots[index:(index+1)])
	u <- knot - from$knots[index]
	to <- list()
	to$knots <- sort(c(from$knots[min(1,index-1):(index-1)], knot,
								 from$knots[(index+2):max(index+2, length(from$knots))]))
# 	X <- bs(data$x, knots=to$knots, degree=degree, intercept=TRUE)
	X <- tp(data$x, knots=to$knots, degree=degree)
	XtXinv <- chol2inv(chol(crossprod(X))) 
	alpha.hat <- XtXinv %*% crossprod(X, data$y)

	df <- length(data$x) - length(to$knots) - length(alpha.hat)
	y.hat <- X %*% alpha.hat
	e2 <- (data$y - y.hat)^2
	s2 <- sum(e2)/df 
	sqrtXtXinv <- sqrt.matrix(XtXinv)
	i <- degree+index
	z.small <- from$z[-i]
	to$alpha <- as.vector( alpha.hat + sqrt(s2) * sqrtXtXinv %*% z.small )
	to$var <- from$var
	to$z <- z.small
	return(list(to=to, u=u, z=from$z[i]))
}


rJump.splitKnot <- function(from, degree, data, sd.knot){
	index <- sample(1:length(from$knots), size=1)
	knot <- from$knots[index]
	knot.range <- min(data$max.x - knot, knot-data$min.x) * .9
	u <- rTruncatedNorm(1, mean=0, sd=sd.knot, min=-knot.range, max=knot.range)
	z <- rnorm(1, mean=0, sd=1)
	to <- rJump.splitKnot.determanistic(from, degree, data, index, u, z)
	ln.density.backward <- dTruncatedNorm(u, mean=0, sd=sd.knot, 
									     min=-knot.range, max=knot.range, log=TRUE) +
									 dnorm(z, log=TRUE)
	out <- list(to=to, ln.density.forward=0, ln.density.backward=ln.density.backward)
	return(out)		
}

	
rJump.splitKnot.determanistic <- function(from, degree, data, index, u, z){
	to <- list()
	knot <- from$knots[index]
	to$knots <- sort(c(from$knots[-index], knot-u, knot+u))
# 	X <- bs(data$x, knots=to$knots, degree=degree, intercept=TRUE)
	X <- tp(data$x, knots=to$knots, degree=degree)
	XtXinv <- chol2inv(chol(crossprod(X))) 
	alpha.hat <- XtXinv %*% crossprod(X, data$y)
	sqrtXtXinv <- sqrt.matrix(XtXinv)
	df <- length(data$x) - length(to$knots) - length(alpha.hat)
	y.hat <- X %*% alpha.hat
	e2 <- (data$y - y.hat)^2
	s2 <- sum(e2)/df 
	i <- degree+index
	z <- c(from$z[min(1,i-1):(i-1)], z, from$z[i:length(from$z)])
	if(any(is.na(z))){ 
		z <- z[-which(is.na(z))]
	}
	to$alpha <- as.vector( alpha.hat + sqrt(s2) * sqrtXtXinv %*% z )
	to$var <- from$var
	to$z <- z
	return(to)
}


rJump.SameNumKnots <- function(from, degree, data, sd.knot){
  to <- list()
  to$knots <- from$knots
  ln.density.forward  <- 0
  ln.density.backward <- 0
  
  ######################################
  ## Jitter the knots
  ######################################
  # how many knots to jitter
  k <- rbinom(1, length(to$knots) - 1, prob=.2) + 1
  k <- 1
  jitter.knots <- sort( sample(1:length(to$knots), k) )
  for(i in jitter.knots){
    to$knots[i] <- rTruncatedNorm(1, mean=from$knots[i], sd=sd.knot, 
                                  min=data$min.x, max=data$max.x)
    ln.density.forward <- ln.density.forward + 
      dTruncatedNorm(to$knots[i], mean=from$knots[i], sd=sd.knot, 
                     min=data$min.x, max=data$max.x, log=TRUE)
    ln.density.backward <- ln.density.backward + 
      dTruncatedNorm(from$knots[i], mean=to$knots[i], sd=sd.knot, 
                     min=data$min.x, max=data$max.x, log=TRUE)
  }
  to$knots <- sort(to$knots)
  
  df <- data$n - length(from$knots) - length(from$alpha)
  
  # W is the design matrix going from -> to
  # X is the design matrix going to -> from
  #	W <- create.design.matrix(data$x, knots.star, degree)
  W <- tp(data$x, knots=to$knots, degree=degree)
  WtWinv <- try( chol2inv(chol(crossprod(W))), silent=TRUE )
  if(class(WtWinv) == 'try-error'){
    # Got too close to the edge; try again.
    return(rJump.SameNumKnots(from, degree, data, sd.knot))
  }		
  alpha.hat <- WtWinv %*% crossprod(W, data$y)
  e <- W %*% alpha.hat - data$y
  s2 <- sum(e^2)/df
  to$var <- rscaleinvchisq(1, df, s2)    
  ln.density.forward <- ln.density.forward + dscaleinvchisq(to$var, df, s2, log=TRUE) 
  
  sqrtWtWinv <- sqrt.matrix(WtWinv)
  to$z <- rnorm(length(alpha.hat))
  to$alpha <- as.vector(alpha.hat + sqrt(s2) * sqrtWtWinv %*% to$z)
  
  ln.density.forward <- ln.density.forward + sum(dnorm(to$z, log=TRUE))
  
  X <- tp(data$x, knots=from$knots, degree=degree)
  XtXinv <- try( chol2inv(chol(crossprod(X))) )
  alpha.hat <- XtXinv %*% crossprod(X, data$y)
  e <- X %*% alpha.hat - data$y
  s2 <- sum(e^2)/df
  
  ln.density.backward <- ln.density.backward + dscaleinvchisq(from$var, df, s2, log=TRUE) 
  ln.density.backward <- ln.density.backward + sum(dnorm(from$z, log=TRUE))
  
  return(list(to=to, 
              ln.density.backward=ln.density.backward, 
              ln.density.forward=ln.density.forward))
  
}


sqrt.matrix <- function(x){
	temp <- eigen(x, symmetric=TRUE)
	e.vec <- temp$vectors
	e.val <- temp$values
	out <- e.vec %*% diag(sqrt(e.val)) %*% t(e.vec)
	return(out)
}




###################################################
##  Scaled Inverse Chi-sq distribution           ##
###################################################
dscaleinvchisq <- function(x, df, s2, log=FALSE){
	out <- df/2*log(df/2) - lgamma(df/2) + df*log(s2)/2 - (df/2+1)*log(x) - df*s2 / (2*x);
	if(log==FALSE){
		out <- exp(out)
	}
	return(out)
}

rscaleinvchisq <- function(n, df, s2){
	X <- rchisq(n, df)
	out <- s2*df/X
	return(out)
}



