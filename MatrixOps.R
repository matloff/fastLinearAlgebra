
dotProd <- function(x,y) {
   if (length(x) != length(y)) stop('x,y must be of the same length')
   sum(x*y)
}

mTimesv <- function(M,v)
{
   if (ncol(M) != length(v)) stop('incompatible sizes')
   mvProduct <- vector(length=length(v))
   for (i in 1:nrow(M))
      mvProduct[i] <- dotProd(M[i,],v)
   mvProduct
}

mTimess <- function() 
{

}

