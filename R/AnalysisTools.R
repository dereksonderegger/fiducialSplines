custom.labels <- function(variable, value){
  if( variable == 'n' || variable == 'sigma'){
    out <- paste(variable, value, sep=' = ')
  }else if( variable == 'method' ){
    out <- paste(value)
  }
  return(out)
}



acceptance.rate <- function(chain){
	N <- length(chain)
	current <- chain[1]
	count <- 0
	for(i in 2:N){
		if( !identical( chain[i], current ) ){
			count <- count + 1
			current <- chain[i]
		}
	}
	return( count / (N-1) )
}

######################################################################
##  Creates a CI plot                                               ##
##                                                                  ##
## levels is a vector of CI levels where one element corresponds    ##
## to one sample and is the smallest alpha-level whereby a CI       ##
## still contains the truth                                         ##
######################################################################
calc.coverage.plot <- function(levels, response=NULL, by=NULL, file='cstar.RData'){
  if(is.vector(levels)){
    N <- length(levels)
  }else if(is.data.frame(levels) & is.null(by)){
    N <- dim(levels)[1]
  }else{
    N <- min( ddply(levels, .variables=by, .fun=function(x){dim(x)[1]})$V1 )
  }
  
  ## Calculate the error bars for the QQ plots 
	## the order statistics are Beta(j, N-j+1)
	## which has mu = j/(N+1)
	## and var = j(N-j+1) / ( (N+1)^2 (N+2) )
	x <- seq(0,1,length=N);
	j <- seq(1,N);
	mu <- j / (N+1);
	std.dev <- sqrt( j*(N-j+1) / ( (N+1)^2 * (N+2) ) );
	c.star <- calc.c.star(N, file=file);

  CR <- data.frame(Nominal=x, 
    ymin=sapply(mu-c.star*std.dev, max, 0), 
    ymax=sapply(mu+c.star*std.dev, min, 1))

  nominal <- seq(0,1, length=N+2)[c(-1, -(N+2))]
  if(is.vector(levels)){
    out <- data.frame(Nominal=nominal, Observed=ecdf(levels)(nominal))
  }else if(is.data.frame(levels) & is.null(by)){
      out <- data.frame(Nominal=nominal, Observed = ecdf(levels[,which(colnames(levels)==response)])(nominal))
  }else{
    out <- ddply(levels, by, 
          .fun=function(x, nominal){
              data.frame(Nominal=nominal, Observed=ecdf(x[,which(colnames(x)==response)])(nominal))}, 
          nominal )
  }
  
  return(list(CR=CR, ObsVsNom=out))  
}





############################################
##  If we have already calculated the c.star
##  value for n and it is in the file 'cstar.R'
##  then we'll just use that number, otherwise
##  we'll calculate the value and store it in
##  the 'cstar.R' file.
##  n - number of simulations in the experiment
#############################################
calc.c.star <- function(n, N=100000, file='cstar.RData'){
	if(file.exists(file)){
		c.star.matrix <- read.table(file);
	}else{
		c.star.matrix <- data.frame( n=10, c.star=2.961925 );
		write.table(c.star.matrix, file=file)
	}
	index <- match(n, c.star.matrix[,1]);
	if(is.na(index)){
		# define the mean and standard deviation vectors of the order stats
		mu <- seq(1,n) / (n+1);
		std.dev <- sqrt( (seq(1,n) * seq(n,1)) / ( (n+1)^2 * (n+2) ) );
		c <- rep(0,N);
		for( i in 1:N ){
			c[i] <- max( abs(sort(runif(n)) - mu) / std.dev );
		}
		c.star <- quantile(c, probs=c(.95));
		c.star.matrix[ dim(c.star.matrix)[1] + 1 , ] <- c(n, c.star);
		write.table(c.star.matrix, file=file);
		out <- c.star;
	}else{
		out <- c.star.matrix[index, 2];
	}
	return(out);
}
 
