cauchydist <- function(p,x){
  H <- vector(mode = "numeric",length = length(x))
  for(i in 1:length(x)){
    H[i] = (1+(x[i] - p)^2)^(-1)
  }
  return(prod(H))
}
dens <- function(p) vapply(p,cauchydist,0,x=x)
# d <- integrate(dens, -Inf,Inf)
# c <-  d$value^(-1)