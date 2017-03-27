prior = function(x) dcauchy(x, 1.5, 0.4)
like = function(x) dnorm(x,6.1,.4)

# Posterior
propto = function(x) prior(x)*like(x)
d = integrate(propto, -Inf, Inf)
post = function(x) propto(x)/d$value
# Plot
par(mar=c(0,0,0,0)+.1, lwd=2)
curve(like, 0, 8, col="red", axes=F, frame=T)
curve(prior, add=TRUE, col="blue")
curve(post, add=TRUE, col="seagreen")
legend("bottomleft", c("Prior","Likelihood","Posterior"), col=c("blue","red","seagreen"), lty=1, bg="white")
