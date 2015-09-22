I <- function( x, step.value=0 ){
  return( (x > step.value) * 1 )
}

trunc <- function( x, trunc.value=0 ){
  out <- x
  out[ which(x < trunc.value)] <- 0
  return(out)
}

# Create the design matrix for the truncated polynomial basis
tp <- function( x, knots, degree ){
  out <- NULL
  for( i in 0:degree ){
    out <- cbind(out, x^i)
  }
  for( i in 1:length(knots) ){
    out <- cbind( out, trunc(x-knots[i])^degree )
  }
  return(out)
}


# q: degree of polynomial
unscaled.fiducial.density <- function(data, alpha, knots, sigma2, degree, 
								num.jacobians=1000, log=FALSE){
	x <- data$x;
	y <- data$y;
  n <- length(x)
	X <- tp(x, knots=knots, degree=degree)
	yhat <- X %*% alpha;
	resid <- y-yhat;
	logL <- sum(dnorm(resid, mean=0, sd=sqrt(sigma2), log=TRUE));
  mdl <- (length(knots) + (length(alpha)+1)/2) * log(n,2)     
		
	if( is.function(num.jacobians) ){
		prior <- do.call(num.jacobians, list(sigma2, knots, alpha))
    out <- logL + log(prior)
	}else{
		J <- create.jacobian(x, y, alpha, knots, sigma2, degree, num.jacobians);
    out <- logL + log(J) - mdl
	}

	if(log == TRUE){
		return(out);
	}else{
		return(exp(out));
	}	
}


	
create.jacobian <- function(x, y, alpha, knots, sigma2, degree, num.jacobians){
	n <- length(x);
	kappa <- length(knots);
	m <- degree+1
	
	start.ptr <- rep(0, kappa+1);                     # start.ptr and tail.ptr are vectors
	tail.ptr <- rep(0, kappa+1);                      # that tell us what observations fall
	start.ptr[1] <- 1;                                # into the region before the first knot
	for(i in 1:(kappa)){                              # the region between the first and second knots...
		tail.ptr[i]  <- which.max( x[x <= knots[i]] )
	}
	start.ptr[2:(kappa+1)] <- tail.ptr[1:kappa]+1;
	tail.ptr[kappa+1] <- n;

	p <- 1 + length(alpha) + length(knots);

	# If the parameters are such that we can't solve for the parameters
	# the Jacobian should be zero and therefore the density will be 0.
	if( any( (tail.ptr-start.ptr+1) < 2 ) ){
		return(0);
	}

	# Figure out the Jacobian matrix
  
	B <- tp(x, knots=knots, degree=degree)
	J <- B
	for(k in 1:kappa){
    J <- cbind(J, alpha[degree+1+k]*trunc(x-knots[k])^(degree-1) * I(x-knots[k]))
	}
#	J <- cbind(J, ( y - B %*% alpha ) / 2 )
  J <- cbind(J, y/2 )
  J <- (1/sigma2) * degree^kappa * J
					
	# randomly select rows and calculate the determinant.  Then take the average of those.
	jacobians <- rep(NA, num.jacobians)
	for(i in 1:num.jacobians){
		index <- generate.random.indices(start.ptr, tail.ptr, p, n);
		jacobians[i] <- abs(det(J[index,])) 
	}
	jacobians <- jacobians / sigma2 
	out <- median(jacobians)	

	return(out);
}

		
generate.random.indices <- function(start.ptr, tail.ptr, p, n){
	index <- NULL
	num.segments <- length(tail.ptr)
	for(i in 1:num.segments){
		index <- append(index, sample(start.ptr[i]:tail.ptr[i], 2));
	}
	index <- append(index, sample( (1:n)[-index], p-2*num.segments ));
	return(index);
}
	

MDL.penalty <- function(n, p, k){
	mdl <- ((p+1)+2*k)/2 * log(n,2) 
	return(mdl)
}
