############################################################
## Truncated_Normal distribution
############################################################
rTruncatedNorm <- function(n=1, mean=0, sd=1, min=-1, max=1){
	u.min <- pnorm( min, mean=mean, sd=sd);
	u.max <- pnorm( max, mean=mean, sd=sd);
	u <- runif(n, min=u.min, max=u.max);
	z <- qnorm(u, mean=mean, sd=sd);
	return(z);	
}

dTruncatedNorm <- function(x, mean=0, sd=1, min=-1, max=1, log=FALSE){
	if( is.vector(x)){
		out <- sapply(x, dTruncatedNorm.simple, mean, sd, min, max, log);
	} else{
		out <- dTruncatedNorm.simple(x, mean, sd, min, max, log);
	}
	return(out);
}
dTruncatedNorm.simple <- function(x, mean=0, sd=1, min=-1, max=1, log=FALSE){
#	print( paste('x=',x,'  min=',min,'  max=', max, sep='' ))
	if( x >= min & x <= max){
		out <- dnorm(x, mean=mean, sd=sd) / 
	            (  pnorm(max, mean=mean, sd=sd) - pnorm(min, mean=mean, sd=sd) );
	}else{
		out <- 0;
	}
	if(log == TRUE){
		out <- log(out);
	}
	return(out);			
}

