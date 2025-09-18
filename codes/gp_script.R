
######
# simulation of phylogenetic partial pooling in occupancy model

library(rethinking)
library(cmdstanr)
library(ape)

#####################
# sim binary occupancy data structure

log_inv_logit <- function(x) log(inv_logit(x))
#log_inv_logit <- function(x) -log( 1 + exp(-x) )

n_spp <- 30
D <- 3 # dimension of random effects
# sim tree and compute distance matrix
a_tree <- rtree(n_spp)
M <- cophenetic(a_tree)
M <- M / max(M)
Omega <- rlkjcorr(1,D,1) # correlation structure of outcome vector

if (FALSE) {
    blank2(w=2)
    par(mfrow=c(1,2))
    plot(a_tree)
    image(M)

}

# build covariance kernel across spp
K <- matrix(NA,n_spp,n_spp)
rho <- 1

# L1 distance kernel
max_cov <- 3
for ( i in 1:n_spp ) for ( j in 1:n_spp )
    K[i,j] <- max_cov*exp( -rho*M[i,j] )
diag(K) <- max_cov + 0.001
LK <- t(chol(K))
LO <- t(chol(Omega))

# we want f = LK * eta * LO'
k <- c( 1.5 , 1 , 2 ) # intercept
f <- matrix(NA,n_spp,D)
eta <- rmvnorm2(n_spp,sigma=rep(1,D),Rho=diag(D))
f <- LK %*% eta %*% t(LO)
for ( i in 1:n_spp ) {
    f[i,] <- f[i,] + k
    if ( exp(f[i,1]) > 35 ) f[i,1] <- log(35)
}

# plot f scales
if (FALSE) {
    par(mfrow=c(1,2))
    plot( exp(f[sid,1]) , exp(f[sid,2]) )
    plot( exp(f[sid,1]) , inv_logit(f[sid,3]) )

}

n_each <- 30
N <- n_spp * n_each # number of surveys
sid <- rep(1:n_spp,each=n_each)
sst <- rep(NA,N)
y <- rep(NA,N)
p <- 0.8 # detect prob when occupied

for ( i in 1:N ) {
    sst[i] <- runif(1,5,30) # temperature at site
    # log prob occupy site
    log_q <- (-0.5)*(( sst[i] - exp(f[sid[i],1]) ) / exp(f[sid[i],2]))^2 + log_inv_logit(f[sid[i],3])
    # sim true occupancy
    Y <- rbern(1,exp(log_q))
    # sim observation
    y[i] <- rbern(1,Y*p)
}

# mean(y)

dat2 <- list(
    N = N,
    D = D,
    N_spp = n_spp,
    sid = sid,
    sst = sst,
    y = y,
    M = M
)

fit2 <- cstan(file="model8_gp_phylo.stan",data=dat2,chains=4,cores=4,cpp_fast=FALSE)
precis(fit2,2)

post <- extract.samples(fit2)
f_est <- apply(post$f,2:3,mean)

blank2(w=2.6)
par(mfrow=c(1,3))

plot( f[,1] , f_est[,1] , xlab="simulated" , ylab="estimated" )
abline(a=0,b=1,lty=2)

plot( f[,2] , f_est[,2] , xlab="simulated" , ylab="estimated" )
abline(a=0,b=1,lty=2)

plot( f[,3] , f_est[,3] , xlab="simulated" , ylab="estimated" )
abline(a=0,b=1,lty=2)

# compare post mean to corr matrix
apply(post$Omega,2:3,mean)
Omega
