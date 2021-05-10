# Generalized Co-Sparse Factor Regression

In multivariate analysis, we routinely observe outcomes which are either continuous, binary, count or mixed types. Often outcomes are interrelated and predictors correlated. In a large dimension, many predictors are irrelevant. Generalized co-sparse factor regression assumes that the outcomes are from the exponential family and models the underlying dependency through a low-rank and sparse coefficient matrix. We express the required coefficient matrix as the sum of unit-rank components with co-sparse singular vectors. 



We propose a greedy procedure which sequentially recovers the sparse unit-rank components of the coefficient matrix, and refers it as **generalized sequential extraction via constrained unit-rank estimation (GOFAR(S))**. In order to be computationally efficient, we also propose an exclusive extraction procedure where the sparse unit-rank components are estimated in parallel and refer it as **generalized exclusive extraction via constrained unit-rank estimation (GOFAR(P))**.



# Example

Here we present examples of using R functions corresponding to the two estimation approach, i.e., **GOFAR(S)** and **GOFAR(P)**. 

### Data simulation

```
rm(list = ls())
# load library
library(gofar)
#
## Model specification:
SD <- 123
set.seed(SD)
n <- 200
p <- 100
pz <- 0
# Model I in the paper
# n <- 200; p <- 300; pz <- 0 ;           # Model II in the paper
# q1 <- 0; q2 <- 30; q3 <- 0               # Similar response cases
q1 <- 15
q2 <- 15
q3 <- 0                   # mixed response cases
nrank <- 3                # true rank
rank.est <- 4             # estimated rank
nlam <- 40                # number of tuning parameter
s <- 1                    # multiplying factor to singular value
snr <- 0.25               # SNR for variance Gaussian error
#
q <- q1 + q2 + q3
respFamily <- c("gaussian", "binomial", "poisson")
family <- list(gaussian(), binomial(), poisson())
familygroup <- c(rep(1, q1), rep(2, q2), rep(3, q3))
cfamily <- unique(familygroup)
nfamily <- length(cfamily)
#
control <- gofar_control()
#
#
## Generate data
D <- rep(0, nrank)
V <- matrix(0, ncol = nrank, nrow = q)
U <- matrix(0, ncol = nrank, nrow = p)
#
U[, 1] <- c(sample(c(1, -1), 8, replace = TRUE), rep(0, p - 8))
U[, 2] <- c(rep(0, 5), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 14))
U[, 3] <- c(rep(0, 11), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 20))
#
if (nfamily == 1) {
  # for similar type response type setting
  V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8, replace =
                                  TRUE) * runif(8, 0.3, 1), rep(0, q - 16))
  V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8, replace =
                                   TRUE) * runif(8, 0.3, 1), rep(0, q - 28))
  V[, 3] <- c(
    sample(c(1, -1), 5, replace = TRUE) * runif(5, 0.3, 1), rep(0, 23),
    sample(c(1, -1), 2, replace = TRUE) * runif(2, 0.3, 1), rep(0, q - 30)
  )
} else {
  # for mixed type response setting
  # V is generated such that joint learning can be emphasised
  V1 <- matrix(0, ncol = nrank, nrow = q / 2)
  V1[, 1] <- c(sample(c(1, -1), 5, replace = TRUE), rep(0, q / 2 - 5))
  V1[, 2] <- c(rep(0, 3), V1[4, 1], -1 * V1[5, 1], 
               sample(c(1, -1), 3, replace = TRUE), rep(0, q / 2 - 8))
  V1[, 3] <- c(V1[1, 1], -1 * V1[2, 1], rep(0, 4), 
               V1[7, 2], -1 * V1[8, 2], sample(c(1, -1), 2, replace = TRUE), 
               rep(0, q / 2 - 10))
  #
  V2 <- matrix(0, ncol = nrank, nrow = q / 2)
  V2[, 1] <- c(sample(c(1, -1), 5, replace = TRUE), rep(0, q / 2 - 5))
  V2[, 2] <- c(rep(0, 3), V2[4, 1], -1 * V2[5, 1], 
               sample(c(1, -1), 3, replace = TRUE), rep(0, q / 2 - 8))
  V2[, 3] <- c(V2[1, 1], -1 * V2[2, 1], rep(0, 4), 
               V2[7, 2], -1 * V2[8, 2], 
               sample(c(1, -1), 2, replace = TRUE), rep(0, q / 2 - 10))
  #
  V <- rbind(V1, V2)
}
U[, 1:3] <- apply(U[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
V[, 1:3] <- apply(V[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#
D <- s * c(4, 6, 5) # signal strength varries as per the value of s
or <- order(D, decreasing = T)
U <- U[, or]
V <- V[, or]
D <- D[or]
C <- U %*% (D * t(V)) # simulated coefficient matrix
intercept <- rep(0.5, q) # specifying intercept to the model:
C0 <- rbind(intercept, C)
#
Xsigma <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
# Simulated data
sim.sample <- gofar_sim(U, D, V, n, Xsigma, C0, familygroup, snr) 
# Dispersion parameter
pHI <- c(rep(sim.sample$sigmaG, q1), rep(1, q2), rep(1, q3)) 
X <- sim.sample$X[1:n, ]
Y <- sim.sample$Y[1:n, ]
# Test data
sim.sample <- gofar_sim(U, D, V, 1000, Xsigma, C0, familygroup, snr)  
Xt <- sim.sample$X
Yt <- sim.sample$Y
#
crossprod(X %*% U)
apply(X, 2, norm, "2")
X0 <- cbind(1, X)
#
# Simulate data with 20% missing entries
miss <- 0.20          # Proportion of entries missing 
t.ind <- sample.int(n * q, size = miss * n * q)
y <- as.vector(Y)
y[t.ind] <- NA
Ym <- matrix(y, n, q)
naind <- (!is.na(Ym)) + 0 # matrix(1,n,q)
misind <- any(naind == 0) + 0
#
```

### Sequential approach: GOFAR(S)
An R user will need the following package to connect to the database: a) “DBI”, for database interface; b) “odbc” for connecting to the database using **DBI** interface.

In addition, user may require some additional R package for downloading, processing and visualizing the data. Run the following commands to install some of the essential packages.


```
# Model fit:
control$epsilon <- 1e-7
control$spU <- 50 / p # estimated percentage of non zero entries in U
control$spV <- 25 / q # estimated percentage of non zero entries in V
control$maxit <- 1000

# Model fitting: GOFAR(S) (full data)
set.seed(SD)
rank.est <- 5
fit.seq <- gofar_s(Y, X, nrank = rank.est, family = family, 
                      nlambda = nlam, familygroup = familygroup, 
                      control = control, nfold = 5)


# Model fitting: GOFAR(S) (missing data)
set.seed(SD)
rank.est <- 5
fit.seq.m <- gofar_s(Ym, X, nrank = rank.est, family = family, 
                      nlambda = nlam, familygroup = familygroup, 
                      control = control, nfold = 5)


```
### Parallel approach: GOFAR(P)


```
# Model fit:
control$epsilon <- 1e-7
control$spU <- 50 / p
control$spV <- 25 / q
control$maxit <- 1000

# Model fitting: GOFAR(P) (full data)
set.seed(SD)
rank.est <- 5
fit.eea <- gofar_p(Y, X, nrank = rank.est, nlambda = nlam,
                       family = family, familygroup = familygroup, 
                       control = control, nfold = 5)

# Model fitting: GOFAR(P) (missing data)
set.seed(SD)
rank.est <- 5
fit.eea.m <- gofar_p(Ym, X, nrank = rank.est, nlambda = nlam,
                      family = family, familygroup = familygroup, 
                      control = control, nfold = 5)
```




## Queries
Please contact authors and creators for any queries related to using the analysis 


-   Aditya Mishra: [mailto](mailto:amishra@flatironinstitute.org)
