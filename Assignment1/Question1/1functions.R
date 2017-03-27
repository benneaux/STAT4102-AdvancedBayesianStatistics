
solve.HPD.beta <-  function(h, y, n, a, b, p, plot=T, ...){
  post_a = y + a
  post_b = n - y + b
  if (post_a < 1 | post_b < 1)
    warning("code assumes mode is not at 0 or 1")
  
  post_mode = (y + a - 1)/(n + a + b - 2)
  
  lt = uniroot(
    f = function(x){
        dbeta(x,post_a, post_b) - h
        },
    lower=0,
    upper=post_mode)$root
  
  ut = uniroot(
    f = function(x){
      dbeta(x,post_a, post_b) - h
      },
    lower = post_mode,
    upper = 1)$root
  
  coverage = pbeta(ut, post_a, post_b) - pbeta(lt, post_a, post_b)
  hpdval = abs(p-coverage) 
  if (plot) {
    th = seq(0, 1, length=5000)
    plot(th, dbeta(th, post_a, post_b),
         t="l", lty=1,xlab=expression(theta),
         ylab="posterior Density", ...)
    abline(h=h)
    segments(ut,0,ut,dbeta(ut,post_a,post_b))
    segments(lt,0,lt,dbeta(lt,post_a,post_b))
    title(bquote(paste("p(", .(round(lt, 5))," < ", theta, " < ",
                       .(round(ut,5)), " | " , y, ") = ",
                       .(round(coverage, 5)), ")")))
  }
  results <<- data.frame(lt,ut,coverage,h)
  return(hpdval)
}

# Coverage Intervals

waldcover <- function(p) {
  x <- 0:n
  fpx <- dbinom(x, n, p)
  phat <- x / n
  low <- phat - z * sqrt(phat * (1 - phat) / n)
  hig <- phat + z * sqrt(phat * (1 - phat) / n)
  inies <- as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
adjwaldcover <- function(p) {
  x <- 0:n
  fpx <- dbinom(x, n, p)
  nadj <- n + (z^2)
  phat <- (1/nadj)*(x + ((z^2)/2))
  low <- phat - z*sqrt(phat * (1 - phat)/nadj)
  hig <- phat + z* sqrt(phat * (1 - phat)/nadj)
  inies <- as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
scorecover <- function(p) {
  x <- 0:n
  fpx <- dbinom(x, n, p)
  foo <- function(x) suppressWarnings(prop.test(x, n = n,
                                                p = p, conf.level = 1-a, correct = FALSE))
  bar <- lapply(x, foo)
  low <- sapply(bar, function(x) x$conf.int[1])
  hig <- sapply(bar, function(x) x$conf.int[2])
  inies <- as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
exactcover <- function(p) {
  x <- 0:n
  fpx <- dbinom(x, n, p)
  foo <- lapply(x, binom.test, n = n, p = p,
                conf.level = 1-a)
  low <- sapply(foo, function(x) x$conf.int[1])
  hig <- sapply(foo, function(x) x$conf.int[2])
  inies <- as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
jeffreyscover <- function(p) {
  x <- 0:n
  fpx <- dbinom(x,n,p)
  low <- qbeta(a/2, x + 0.5, n - x + 0.5)
  hig <- qbeta(1-(a/2), x + 0.5, n - x + 0.5)
  inies <- as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}

blcover <- function(p) {
  x <- 0:n
  fpx <- dbinom(x,n,p)
  low <- qbeta(a/2, x + 1, n - x + 1)
  hig <- qbeta(1-(a/2), x + 1, n - x + 1)
  inies <- as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
