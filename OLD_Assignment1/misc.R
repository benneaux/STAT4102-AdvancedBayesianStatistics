test <- tbl_df(waldcover(p))
test <- vapply(p,waldcover,0)

coverage <- 0

n <- 50
a <- 0.05
p <- seq(0.0001,0.9999,1/1000)
z <- abs(qnorm(.5*a,0,1))

for(i in 1:length(p)){
  p = p[i]
  fpx <- dbinom(x, n, p)
  for(j in 1:n){
    phat  = j / n
    low   = phat - z * sqrt(phat * (1 - phat) / n)
    hig   = phat + z * sqrt(phat * (1 - phat) / n)
    abs(hig-low)
  }
  
  
  hpd(qbeta,shape1 = 2, shape2 = 10)
  
  f=function(p){
    
    # find the symmetric
    g=function(x){return(x-p*((1-p)/(1-x))^(115.5/117.5))}
    return(uniroot(g,c(.504,.99))$root)}
  
  ff=function(alpha){
    
    # find the coverage
    g=function(x){return(x-p*((1-p)/(1-x))^(115.5/117.5))}
    return(uniroot(g,c(.011,.49))$root)}
  
  vapply(p,ff(0.95),0)
 ##################################################################### 
  HDIofICDF = function(ICDF,
                       credMass=0.95,
                       tol=1e-8,
                       ...) {

      incredMass    = 1.0 - credMass
      intervalWidth = function(lowTailPr ,
                               ICDF ,
                               credMass,
                               ...) 
      {
        ICDF(credMass + lowTailPr, ...) - ICDF(lowTailPr, ...)
      }
      optInfo        = optimize(intervalWidth, 
                                c(0, incredMass) ,
                                ICDF = ICDF ,
                                credMass = credMass , 
                                tol = tol,
                                ...)
      HDIlowTailPr    = optInfo$minimum
      
      return(c(ICDF(HDIlowTailPr, ...),
               ICDF(credMass + HDIlowTailPr, ...) ) )
      }
  
  HDIofICDF(qbeta, shape1 = 30, shape2 = 12)
  # Arguments:
  # ICDFname is R’s name for the inverse cumulative density function
  # of the distribution.
  # credMass is the desired mass of the HDI region.
  # tol is passed to R’s optimize function.
  # Return value:
  # Highest density iterval (HDI) limits in a vector.
  # Example of use: For determining HDI of a beta(30,12) distribution, type
  # HDIofICDF( qbeta , shape1 = 30 , shape2 = 12 )
  # Notice that the parameters of the ICDFname must be explicitly named;
  # e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
  # Adapted and corrected from Greg Snow’s TeachingDemos package.