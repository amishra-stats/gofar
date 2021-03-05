#' Default control parameters for Generalized co-sparse factor regresion
#'
#' Control parameters for GOFAR(S) and GOFAR(P)
#'
#' @param lamMaxFac a multiplier of calculated lambda_max
#' @param lamMinFac a multiplier of determing lambda_min as a fraction of lambda_max
#' @param gamma0 power parameter in the adaptive weights
#' @param elnetAlpha elastic net penalty parameter
#' @param spU maximum proportion of nonzero elements in each column of U
#' @param spV maximum proportion of nonzero elements in each column of V
#' @param maxit maximum iteration for each sequential steps
#' @param epsilon tolerence value set for convergene of gcure
#' @param se1 apply 1se sule for the model;
#' @param equalphi dispersion parameter for all gaussian outcome equal or not 0/1
#' @param objI 1 or 0 convergence on the basis of objective function or not
#' @param initmaxit maximum iteration for initialization problem
#' @param initepsilon tolerence value for convergene in the initialization problem
#' @param alp scaling factor corresponding to poisson outcomes
#' @return a list of controling parameter.
#' @export
#' @examples
#' \donttest{
#' # control variable for GOFAR(S) and GOFAR(P)
#' control <- gofar_control()
#' }
#' @useDynLib gofar
gofar_control <- function(maxit = 5000, epsilon = 1e-6,
                          elnetAlpha = 0.95,
                          gamma0 = 1, se1 = 1,
                          spU = 0.5, spV = 0.5,
                          lamMaxFac = 1, lamMinFac = 1e-6,
                          initmaxit = 2000, initepsilon = 1e-6,
                          equalphi = 1, objI = 1, alp = 60) {
  list(
    lamMaxFac = lamMaxFac,
    lamMinFac = lamMinFac,
    gamma0 = gamma0, elnetAlpha = elnetAlpha, spU = spU, spV = spV,
    maxit = maxit, epsilon = epsilon,
    objI = objI, se1 = se1,
    initmaxit = initmaxit, initepsilon = initepsilon,
    equalphi = equalphi, alp = alp
  )
}



#' Simulate data
#'
#' Genertate random samples from a generalize sparse factor regression model
#'
#' @param U specified value of U
#' @param V specified value of V
#' @param D specified value of D
#' @param n sample size
#' @param snr signal to noise ratio specified for gaussian type outcomes
#' @param Xsigma covariance matrix for generating sample of X
#' @param C0 Specified coefficient matrix with first row being intercept
#' @param familygroup parameter defining correlated error
#' @return
#'   \item{Y}{Generated response matrix}
#'   \item{X}{Generated predictor matrix}
#'   \item{sigmaG}{standard deviation for gaussian error}
#' @export
#' @useDynLib gofar
#' @import magrittr
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rbinom rpois family gaussian binomial poisson
#' @importFrom stats coef sd glm.control dbinom dnorm dpois glm
#' @examples
#' \donttest{
#' ## Model specification:
#' SD <- 123
#' set.seed(SD)
#' n <- 200
#' p <- 100
#' pz <- 0
#' # Model I in the paper
#' # n <- 200; p <- 300; pz <- 0 ;           # Model II in the paper
#' # q1 <- 0; q2 <- 30; q3 <- 0               # Similar response cases
#' q1 <- 15
#' q2 <- 15
#' q3 <- 0 # mixed response cases
#' nrank <- 3 # true rank
#' rank.est <- 4 # estimated rank
#' nlam <- 40 # number of tuning parameter
#' s <- 1 # multiplying factor to singular value
#' snr <- 0.25 # SNR for variance Gaussian error
#' #
#' q <- q1 + q2 + q3
#' respFamily <- c("gaussian", "binomial", "poisson")
#' family <- list(gaussian(), binomial(), poisson())
#' familygroup <- c(rep(1, q1), rep(2, q2), rep(3, q3))
#' cfamily <- unique(familygroup)
#' nfamily <- length(cfamily)
#' #
#' control <- gofar_control()
#' #
#' #
#' ## Generate data
#' D <- rep(0, nrank)
#' V <- matrix(0, ncol = nrank, nrow = q)
#' U <- matrix(0, ncol = nrank, nrow = p)
#' #
#' U[, 1] <- c(sample(c(1, -1), 8, replace = TRUE), rep(0, p - 8))
#' U[, 2] <- c(rep(0, 5), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 14))
#' U[, 3] <- c(rep(0, 11), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 20))
#' #
#' if (nfamily == 1) {
#'   # for similar type response type setting
#'   V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 16))
#'   V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 28))
#'   V[, 3] <- c(
#'     sample(c(1, -1), 5, replace = TRUE) * runif(5, 0.3, 1), rep(0, 23),
#'     sample(c(1, -1), 2, replace = TRUE) * runif(2, 0.3, 1), rep(0, q - 30)
#'   )
#' } else {
#'   # for mixed type response setting
#'   # V is generated such that joint learning can be emphasised
#'   V1 <- matrix(0, ncol = nrank, nrow = q / 2)
#'   V1[, 1] <- c(sample(c(1, -1), 5, replace = TRUE), rep(0, q / 2 - 5))
#'   V1[, 2] <- c(
#'     rep(0, 3), V1[4, 1], -1 * V1[5, 1],
#'     sample(c(1, -1), 3, replace = TRUE), rep(0, q / 2 - 8)
#'   )
#'   V1[, 3] <- c(
#'     V1[1, 1], -1 * V1[2, 1], rep(0, 4),
#'     V1[7, 2], -1 * V1[8, 2], sample(c(1, -1), 2, replace = TRUE),
#'     rep(0, q / 2 - 10)
#'   )
#'   #
#'   V2 <- matrix(0, ncol = nrank, nrow = q / 2)
#'   V2[, 1] <- c(sample(c(1, -1), 5, replace = TRUE), rep(0, q / 2 - 5))
#'   V2[, 2] <- c(
#'     rep(0, 3), V2[4, 1], -1 * V2[5, 1],
#'     sample(c(1, -1), 3, replace = TRUE), rep(0, q / 2 - 8)
#'   )
#'   V2[, 3] <- c(
#'     V2[1, 1], -1 * V2[2, 1], rep(0, 4),
#'     V2[7, 2], -1 * V2[8, 2],
#'     sample(c(1, -1), 2, replace = TRUE), rep(0, q / 2 - 10)
#'   )
#'   #
#'   V <- rbind(V1, V2)
#' }
#' U[, 1:3] <- apply(U[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' V[, 1:3] <- apply(V[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' #
#' D <- s * c(4, 6, 5) # signal strength varries as per the value of s
#' or <- order(D, decreasing = T)
#' U <- U[, or]
#' V <- V[, or]
#' D <- D[or]
#' C <- U %*% (D * t(V)) # simulated coefficient matrix
#' intercept <- rep(0.5, q) # specifying intercept to the model:
#' C0 <- rbind(intercept, C)
#' #
#' Xsigma <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
#' # Simulated data
#' sim.sample <- gofar_sim(U, D, V, n, Xsigma, C0, familygroup, snr)
#' # Dispersion parameter
#' pHI <- c(rep(sim.sample$sigmaG, q1), rep(1, q2), rep(1, q3))
#' X <- sim.sample$X[1:n, ]
#' Y <- sim.sample$Y[1:n, ]
#' }
gofar_sim <- function(U, D, V, n, Xsigma, C0, familygroup, snr) {
  ## finding basis along more number of columns of data vector
  basis.vec <- function(x) {
    # require(Matrix)
    if (diff(dim(x)) < 0) x <- t(x)
    qd <- qr(x)
    k <- qr.Q(qd) %*% qr.R(qd)[, 1:qd$rank]
    k[abs(k) < 1e-6] <- 0
    b.ind <- vector()
    for (i in 1:qd$rank) {
      b.ind[i] <- which(apply(x, 2, function(x, y) sum(abs(x - y)), k[, i]) < 1e-6)[1]
    }
    return(list(ind = b.ind, vec = x[, b.ind]))
  }

  p <- nrow(U)
  q <- nrow(V)
  nrank <- ncol(U)

  U.t <- diag(max(dim(U)))
  U.t <- U.t[, -basis.vec(U)$ind]
  P <- cbind(U, U.t)
  UtXsUt <- t(U.t) %*% Xsigma %*% U.t
  UtXsU <- t(U.t) %*% Xsigma %*% U
  UXsU <- t(U) %*% Xsigma %*% U
  UXsUinv <- solve(UXsU)
  ## sigma.X2 <- t(U.t)%*%Xsigma%*%U.t - t(U.t)%*%Xsigma%*%U%*%solve(t(U)%*%Xsigma%*%U)%*%t(U)%*%Xsigma%*%U.t
  sigma.X2 <- UtXsUt - UtXsU %*% UXsUinv %*% t(UtXsU)
  sigma.X2 <- (sigma.X2 + t(sigma.X2)) / 2


  X1 <- matrix(nrow = nrank, ncol = n, rnorm(n * nrank))
  ## X1 <- t(mvrnorm(n,rep(0,ncol(U)),diag(ncol(U)) ))
  mean.X2 <- UtXsU %*% UXsUinv %*% X1
  ## mean.X2 <- t(U.t)%*%Xsigma%*%U%*%solve(t(U)%*%Xsigma%*%U)%*%X1
  X2 <- mean.X2 + t(MASS::mvrnorm(ncol(mean.X2), rep(0, nrow(mean.X2)), sigma.X2))
  X <- t(solve(t(P)) %*% rbind(X1, X2)) # /sqrt(n)
  # crossprod(X%*%U)

  X0 <- cbind(1, X)
  MU <- X0 %*% C0
  sigmaG <- 0
  ## Generation  of Y and ty.lam to find lambda max in the selection process
  Y <- matrix(nrow = n, ncol = q, 0)
  # sigma <- rho
  cfamily <- unique(familygroup)
  family <- list(gaussian(), binomial(), poisson())
  for (i in cfamily) {
    qq <- familygroup == i
    if (i == 1) {
      rho <- 0 # snr <- 0.5
      qt <- sum(qq)
      svdrr <- eigen(rho^abs(outer(1:qt, 1:qt, FUN = "-")))
      svdrrinv <- svdrr$vectors %*% diag(svdrr$values^0.5, nrow = qt) %*%
        t(svdrr$vectors)
      UU <- matrix(nrow = n, ncol = qt, rnorm(n * qt, 0, 1)) %*%
        svdrrinv
      Y3 <- X %*% U[, nrank] %*% t(V[1:qt, nrank]) * D[nrank]
      sigma <- sqrt(sum(Y3^2) / sum(UU^2)) / snr
      sigmaG <- sigma
      UU <- UU * sigma
      # Y <- X %*% C + UU
      Y[, qq] <- MU[, qq] + UU
    } else if (i == 2) {
      prob <- as.matrix(family[[2]]$linkinv(MU[, qq]))
      Y[, qq] <- apply(prob, 2, function(a) rbinom(n = n, size = 1, a))
    } else if (i == 3) {
      prob <- as.matrix(family[[3]]$linkinv(MU[, qq]))
      Y[, qq] <- apply(prob, 2, function(a) rpois(n = n, lambda = a))
    }
  }
  return(list(Y = Y, X = X, sigmaG = sigmaG^2))
}


















#' Generalize Exclusive factor extraction via co-sparse unit-rank estimation (GOFAR(P)) using k-fold crossvalidation
#'
#' Divide and conquer approach for low-rank and sparse coefficent matrix estimation: Exclusive extraction
#'
#' @param Yt response matrix
#' @param X covariate matrix; when X = NULL, the fucntion performs unsupervised learning
#' @param nrank an integer specifying the desired rank/number of factors
#' @param nlambda number of lambda values to be used along each path
#' @param family set of family [gaussian, bernoulli, possion]
#' @param familygroup index set of the type of multivariate outcomes: "1" for Gaussian, "2" for Bernoulli, "3" for Poisson outcomes
#' @param cIndex control index, specifying index of control variable in the design matrix X
#' @param ofset offset matrix specified
#' @param control a list of internal parameters controlling the model fitting
#' @param nfold number of fold for cross-validation
#' @param PATH TRUE/FALSE for generating solution path of sequential estimate after cross-validation step
#' @return
#'   \item{C}{estimated coefficient matrix; based on GIC}
#'   \item{Z}{estimated control variable coefficient matrix}
#'   \item{Phi}{estimted dispersion parameters}
#'   \item{U}{estimated U matrix (generalize latent factor weights)}
#'   \item{D}{estimated singular values}
#'   \item{V}{estimated V matrix (factor loadings)}
#'   \item{lam}{selected lambda values based on the chosen information criterion}
#'   \item{lampath}{sequences of lambda values used in model fitting. In each sequential unit-rank estimation step,
#'   a sequence of length nlambda is first generated between (lamMax*lamMaxFac, lamMax*lamMaxFac*lamMinFac) equally
#'   spaced on the log scale, in which lamMax is estimated and the other parameters are specified in gofar_control.
#'   The model fitting starts from the largest lambda and stops when the maximum proportion of nonzero elements is reached in
#'   either u or v, as specified by spU and spV in gofar_control.}
#'   \item{IC}{values of information criteria}
#'   \item{Upath}{solution path of U}
#'   \item{Dpath}{solution path of D}
#'   \item{Vpath}{solution path of D}
#'   \item{ObjDec}{boolian type matrix outcome showing if objective function is monotone decreasing or not.}
#'   \item{familygroup}{spcified familygroup of outcome variables.}
#' @export
#' @import magrittr
#' @useDynLib gofar
#' @examples
#' \donttest{
#' ## Model specification:
#' SD <- 123
#' set.seed(SD)
#' n <- 200
#' p <- 100
#' pz <- 0
#' # Model I in the paper
#' # n <- 200; p <- 300; pz <- 0 ;           # Model II in the paper
#' # q1 <- 0; q2 <- 30; q3 <- 0               # Similar response cases
#' q1 <- 15
#' q2 <- 15
#' q3 <- 0 # mixed response cases
#' nrank <- 3 # true rank
#' rank.est <- 4 # estimated rank
#' nlam <- 40 # number of tuning parameter
#' s <- 1 # multiplying factor to singular value
#' snr <- 0.25 # SNR for variance Gaussian error
#' #
#' q <- q1 + q2 + q3
#' respFamily <- c("gaussian", "binomial", "poisson")
#' family <- list(gaussian(), binomial(), poisson())
#' familygroup <- c(rep(1, q1), rep(2, q2), rep(3, q3))
#' cfamily <- unique(familygroup)
#' nfamily <- length(cfamily)
#' #
#' control <- gofar_control()
#' #
#' #
#' ## Generate data
#' D <- rep(0, nrank)
#' V <- matrix(0, ncol = nrank, nrow = q)
#' U <- matrix(0, ncol = nrank, nrow = p)
#' #
#' U[, 1] <- c(sample(c(1, -1), 8, replace = TRUE), rep(0, p - 8))
#' U[, 2] <- c(rep(0, 5), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 14))
#' U[, 3] <- c(rep(0, 11), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 20))
#' #
#' if (nfamily == 1) {
#'   # for similar type response type setting
#'   V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 16))
#'   V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 28))
#'   V[, 3] <- c(
#'     sample(c(1, -1), 5, replace = TRUE) * runif(5, 0.3, 1), rep(0, 23),
#'     sample(c(1, -1), 2, replace = TRUE) * runif(2, 0.3, 1), rep(0, q - 30)
#'   )
#' } else {
#'   # for mixed type response setting
#'   # V is generated such that joint learning can be emphasised
#'   V1 <- matrix(0, ncol = nrank, nrow = q / 2)
#'   V1[, 1] <- c(sample(c(1, -1), 5, replace = TRUE), rep(0, q / 2 - 5))
#'   V1[, 2] <- c(
#'     rep(0, 3), V1[4, 1], -1 * V1[5, 1],
#'     sample(c(1, -1), 3, replace = TRUE), rep(0, q / 2 - 8)
#'   )
#'   V1[, 3] <- c(
#'     V1[1, 1], -1 * V1[2, 1], rep(0, 4),
#'     V1[7, 2], -1 * V1[8, 2], sample(c(1, -1), 2, replace = TRUE),
#'     rep(0, q / 2 - 10)
#'   )
#'   #
#'   V2 <- matrix(0, ncol = nrank, nrow = q / 2)
#'   V2[, 1] <- c(sample(c(1, -1), 5, replace = TRUE), rep(0, q / 2 - 5))
#'   V2[, 2] <- c(
#'     rep(0, 3), V2[4, 1], -1 * V2[5, 1],
#'     sample(c(1, -1), 3, replace = TRUE), rep(0, q / 2 - 8)
#'   )
#'   V2[, 3] <- c(
#'     V2[1, 1], -1 * V2[2, 1], rep(0, 4),
#'     V2[7, 2], -1 * V2[8, 2],
#'     sample(c(1, -1), 2, replace = TRUE), rep(0, q / 2 - 10)
#'   )
#'   #
#'   V <- rbind(V1, V2)
#' }
#' U[, 1:3] <- apply(U[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' V[, 1:3] <- apply(V[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' #
#' D <- s * c(4, 6, 5) # signal strength varries as per the value of s
#' or <- order(D, decreasing = T)
#' U <- U[, or]
#' V <- V[, or]
#' D <- D[or]
#' C <- U %*% (D * t(V)) # simulated coefficient matrix
#' intercept <- rep(0.5, q) # specifying intercept to the model:
#' C0 <- rbind(intercept, C)
#' #
#' Xsigma <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
#' # Simulated data
#' sim.sample <- gofar_sim(U, D, V, n, Xsigma, C0, familygroup, snr)
#' # Dispersion parameter
#' pHI <- c(rep(sim.sample$sigmaG, q1), rep(1, q2), rep(1, q3))
#' X <- sim.sample$X[1:n, ]
#' Y <- sim.sample$Y[1:n, ]
#'
#'
#'
#' # Test data
#' sim.sample <- gofar_sim(U, D, V, 1000, Xsigma, C0, familygroup, snr)
#' Xt <- sim.sample$X
#' Yt <- sim.sample$Y
#' #
#' crossprod(X %*% U)
#' apply(X, 2, norm, "2")
#' X0 <- cbind(1, X)
#' #
#' # Simulate data with 20% missing entries
#' miss <- 0.20 # Proportion of entries missing
#' t.ind <- sample.int(n * q, size = miss * n * q)
#' y <- as.vector(Y)
#' y[t.ind] <- NA
#' Ym <- matrix(y, n, q)
#' naind <- (!is.na(Ym)) + 0 # matrix(1,n,q)
#' misind <- any(naind == 0) + 0
#' #
#' # Model fitting begins:
#' control$epsilon <- 1e-7
#' control$spU <- 50 / p
#' control$spV <- 25 / q
#' control$maxit <- 1000
#' # Model fitting: GOFAR(P) (full data)
#' set.seed(SD)
#' rank.est <- 5
#' fit.eea <- gofar_p(Y, X,
#'   nrank = rank.est, nlambda = nlam,
#'   family = family, familygroup = familygroup,
#'   control = control, nfold = 5
#' )
#'
#' # Model fitting: GOFAR(P) (missing data)
#' set.seed(SD)
#' rank.est <- 5
#' fit.eea.m <- gofar_p(Ym, X,
#'   nrank = rank.est, nlambda = nlam,
#'   family = family, familygroup = familygroup,
#'   control = control, nfold = 5
#' )
#' }
#' @references
#' Mishra, A., Dey, D., Chen, K. (2019) \emph{ Generalized co-sparse factor regression, In prepration}
gofar_p <- function(Yt, X, nrank = 3, nlambda = 40, family,
                       familygroup = NULL, cIndex = NULL, ofset = NULL,
                       control = list(), nfold = 5, PATH = FALSE) {
  # Yt = Y;nrank =rank.est; nlambda = nlam; cIndex = NULL; ofset=NULL;  nfold = 5
  cat("Initializing...", "\n")
  n <- nrow(Yt)
  p <- ncol(X)
  q <- ncol(Yt)

  if (is.null(cIndex)) {
    cIndex <- 1
    pz <- 0
    p <- ncol(X)
  } else {
    pz <- length(cIndex)
    p <- ncol(X) - pz
    cIndex <- c(1, cIndex + 1)
  }
  if (is.null(familygroup)) familygroup <- rep(1, q)

  U <- matrix(0, nrow = p, ncol = nrank)
  V <- matrix(0, nrow = q, ncol = nrank)
  D <- rep(0, nrank)
  Z <- matrix(0, nrow = pz + 1, ncol = q)
  PHI <- rep(0, q)
  lamSel <- rep(0, nrank)
  X0 <- cbind(1, X)

  xtem <- getScaleGaussian(Yt, X0, familygroup)
  mx <- xtem$ko
  Y <- xtem$Y


  if (is.null(ofset)) ofset <- matrix(0, nrow = n, ncol = q)
  aft <- glmCol(Y, X0, ofset, family, familygroup, q, cIndex)
  Z0 <- aft$Z
  PHI0 <- aft$PHI

  Yin <- Y
  naind <- (!is.na(Y)) + 0 # matrix(1,n,q)
  misind <- any(naind == 0) + 0
  if (misind == 1) Yin[is.na(Y)] <- 0


  ## Initialization
  kappaco <- getKappaC0(X0, familygroup, control$alp)
  ndev <- getNullDev(Yin, ofset, familygroup, naind)
  xx <- gcure_cpp_init2(Yin, X0, nrank,
    cindex = cIndex,
    ofset = ofset, familygroup,
    Zini = Z0, PhiIni = PHI0, kappaC0 = kappaco,
    control, misind, naind, ndev
  )
  XC <- X0[, -cIndex] %*% xx$C[-cIndex, ]
  Z0 <- matrix(xx$C[cIndex, ], ncol = q)
  Phi0 <- xx$PHI
  svdxc <- svd(XC)
  nkran <- sum(svdxc$d > 1e-07)
  V0 <- as.matrix(svdxc$v[, 1:nkran])
  D0 <- svdxc$d[1:nkran] / sqrt(n)
  U0 <- xx$C[-cIndex, ] %*% V0 %*% diag(1 / D0, nrow = nkran, ncol = nkran)
  nrank <- nkran

  # save(list = ls(), file = 'aditya.rda')
  # stop('aditya')

  U <- U0
  V <- V0
  D <- D0
  Z <- Z0
  PHI <- Phi0
  lamSel <- rep(0, nrank)
  X0 <- cbind(1, X)
  # totTime <- 0


  ## Begin exclusive extraction procedure
  # if(is.null(ofset)) ofset <- matrix(0,nrow=n,ncol=q)
  # C0 <- tcrossprod(tcrossprod(U0,diag(D0,nrow = nrank)),V0)

  N.nna <- sum(!is.na(Y))
  naind <- (!is.na(Y)) + 0 # matrix(1,n,q)
  ind.nna <- which(!is.na(Y))
  fit.nlayer <- fit.nfold <- vector("list", nrank)

  for (k in 1:nrank) { # desired rank extraction
    initW <- list(
      wu = abs(U0[, k])^-control$gamma0,
      wd = abs(D0[k])^-control$gamma0,
      wv = abs(V0[, k])^-control$gamma0
    )
    ## define ofset
    C0 <- tcrossprod(tcrossprod(U0[, -k], diag(D0[-k], nrow = nrank - 1)),
                     V0[, -k])
    ofset1 <- tcrossprod(X0[, -cIndex], t(C0))
    lambda.max <- get.lam.max2(Y, X, familygroup, ofset1)
    kappaco <- getKappaC0(X0, familygroup, control$alp)
    ndev <- getNullDev(Y, ofset1, familygroup, naind)
    ## store the deviance of the test data
    dev <- matrix(NA, nfold, nlambda)
    sdcal <- matrix(NA, nfold, nlambda)
    tec <- rep(0, nfold)

    ID <- rep(1:nfold, len = N.nna)
    ID <- sample(ID, N.nna, replace = FALSE)

    fitT <- vector("list", nfold)
    for (ifold in 1:nfold) { # ifold=1
      cat('Fold ', ifold, '\n')
      ind.test <- ind.nna[which(ID == ifold)]
      Yte <- Y
      Yte[-ind.test] <- NA
      Ytr <- Y
      Ytr[ind.test] <- NA

      zerosol <- glmCol(Ytr, X0, ofset1, family, familygroup, q, cIndex)

      naind2 <- (!is.na(Ytr)) + 0 # matrix(1,n,q)
      misind <- any(naind2 == 0) + 0
      Ytr[is.na(Ytr)] <- 0
      ndev <- getNullDev(Ytr, ofset1, familygroup, naind2)

      fitT[[ifold]] <- gcure_cpp_miss(Ytr, X0,
        nlam = nlambda,
        cindex = cIndex,
        ofset = ofset1, familygroup, initw = initW,
        Dini = D0[k], Zini = Z, PhiIni = PHI,
        Uini = as.matrix(U0[, k]),
        Vini = as.matrix(V0[, k]),
        kappaC0 = kappaco,
        lmax = lambda.max, control, misind,
        naind2, ndev, 0, zerosol,
        control$maxit,
        control$epsilon
      )
      fitF <- fitT[[ifold]]
      if (fitF$nkpath == 1) {
        im <- fitF$nkpath
        lam <- fitF$lamKpath[im, 1]
        insel <- which(fitF$lamseq == lam)
        mu.test <- fitF$mukpath[, , im]
        mu.test[-ind.test] <- NA
        # dev[ifold, ] <- objFun2(Yte, mu.test, fitF$phipath[,im],familygroup)
        # dev[ifold, ] <- dev[ifold, ]/sum((!is.na(Yte)))
        tttval <- objFun5(Yte, mu.test, fitF$phipath[, im], familygroup)
        dev[ifold, ] <- tttval[1] #+ tval
        sdcal[ifold, ] <- tttval[2]
        tec[ifold] <- tttval[3]
      } else {
        for (im in 1:fitF$nkpath) {
          lam <- fitF$lamKpath[im, 1]
          insel <- which(fitF$lamseq == lam)
          mu.test <- fitF$mukpath[, , im]
          mu.test[-ind.test] <- NA
          tval <- control$elnetAlpha * lam * fitF$dkpath[im, 1] * initW$wd * sum(
            initW$wu * abs(fitF$ukpath[, im])
          ) * sum(initW$wv * abs(fitF$vkpath[, im]))
          abc <- sum((fitF$vkpath[, im])^2) * sum((fitF$ukpath[, im])^2) * (
            fitF$dkpath[im, 1])^2
          tval <- tval + lam * (1 - control$elnetAlpha) * abc

          # dev[ifold, insel] <- dev[ifold, insel]/sum((!is.na(Yte)))
          tttval <- objFun5(Yte, mu.test, fitF$phipath[, im], familygroup)
          dev[ifold, insel] <- tttval[1] #+ tval
          sdcal[ifold, insel] <- tttval[2]
          tec[ifold] <- tttval[3]
        }
      }
    }
    fit.nfold[[k]] <- fitT



    # select lambda: 1se rule
    dev.mean <- colMeans(dev, na.rm = FALSE)
    sderr <- sqrt(apply(sdcal, 2, sum, na.rm = FALSE) / (nfold * sum(tec - 1)))
    l.mean <- which.min(dev.mean)
    if (is.na(sderr[l.mean])) {
      l.mean <- min(which(dev.mean <= (dev.mean[l.mean])))
    } else {
      l.mean <- min(which(dev.mean <= (dev.mean[l.mean] +
        control$se1 * sderr[l.mean])))
    }
    lamS <- fitF$lamseq[l.mean]


    Yf <- Y
    zerosol <- glmCol(Yf, X0, ofset1, family, familygroup, q, cIndex)

    naind1 <- (!is.na(Yf)) + 0 # matrix(1,n,q)
    misind <- any(naind1 == 0) + 0
    Yf[is.na(Yf)] <- 0
    ndev <- getNullDev(Yf, ofset1, familygroup, naind1)


    if ( PATH ) {
      fit.nlayer[[k]] <- fit.layer <- gcure_cpp_miss(Yf,X0,
                                                     nlam = nlambda,
                                                     cindex=cIndex,
                                                     ofset = ofset1,
                                                     familygroup,
                                                     initw=initW,
                                                     Dini = D0[k],
                                                     Zini = Z, PhiIni=PHI,
                                                     Uini = as.matrix(U0[,k]),
                                                     Vini = as.matrix(V0[,k]),
                                                     kappaC0 = kappaco,
                                                     lmax = lambda.max,
                                                     control,misind,
                                                     naind1,ndev,0,zerosol,
                                                     control$maxit,
                                                     control$epsilon)
      l.mean <- which(fit.layer$lamKpath==lamS)
      if(length(l.mean)==0){
        l.mean <- 1;
      }
    } else {
      control1 <- control
      control1$lamMinFac <- 1
      fit.nlayer[[k]] <- fit.layer <- gcure_cpp_miss(Yf, X0,
                                                     nlam = 1,
                                                     cindex = cIndex,
                                                     ofset = ofset1, familygroup,
                                                     initw = initW,
                                                     Dini = D0[k],
                                                     Zini = Z, PhiIni = PHI,
                                                     Uini = as.matrix(U0[, k]),
                                                     Vini = as.matrix(V0[, k]),
                                                     kappaC0 = kappaco,
                                                     # lmax = lambda.max,
                                                     lmax = lamS,
                                                     control1, misind,
                                                     naind1, ndev, 0, zerosol,
                                                     control$maxit,
                                                     control$epsilon)
      l.mean <- 1
    }

    fit.nlayer[[k]]$initw <- initW
    fit.nlayer[[k]]$lamseq <- fitF$lamseq   ## save sequence of lambda path
    fit.nlayer[[k]]$dev <- dev
    fit.nlayer[[k]]$lamS <- lamS
    fit.nlayer[[k]]$sderr <- sderr

    U[, k] <- fit.layer$ukpath[, l.mean]
    V[, k] <- fit.layer$vkpath[, l.mean]
    D[k] <- fit.layer$dkpath[l.mean, 1]
    cat(D[k], '\n')
    lamSel[k] <- fit.layer$lamKpath[l.mean, 1]

    Ck <- D[k] * tcrossprod(U[, k], V[, k])
    PHI <- fit.layer$phipath[, l.mean]
    if (pz == 0) {
      Z[1, ] <- drop(fit.layer$zpath[, ,l.mean])
    } else {
      Z <- fit.layer$zpath[, , l.mean]
    }
  }











  ind <- D != 0
  cat("Estimated rank =", sum(ind), "\n")
  U <- matrix(U[, ind], ncol = sum(ind))
  V <- matrix(V[, ind], ncol = sum(ind))
  D <- D[ind]
  fit.nlayer <- fit.nlayer[ind]
  if (sum(ind) == 0) {
    U <- matrix(rep(0, p), ncol = 1)
    V <- matrix(rep(0, q), ncol = 1)
    D <- 0
  }
  ind <- order(D, decreasing = T)
  ft1 <- list(
    fit = fit.nlayer, C = U %*% (D * t(V)), Z = Z, Phi = PHI,
    U = matrix(U[, ind], ncol = length(ind)),
    V = matrix(V[, ind], ncol = length(ind)), D = D[ind],
    lam = lamSel, familygroup = familygroup
  )
  if (all(unique(familygroup) == 1:2)) {
    ft1 <- updateFitObject(ft1, mx)
  }
  return(ft1)
}












#' Generalize Sequential factor extraction via co-sparse unit-rank estimation (GOFAR(S)) using k-fold crossvalidation
#'
#' Divide and conquer approach for low-rank and sparse coefficent matrix estimation: Sequential
#'
#' @param Yt response matrix
#' @param X covariate matrix; when X = NULL, the fucntion performs unsupervised learning
#' @param nrank an integer specifying the desired rank/number of factors
#' @param nlambda number of lambda values to be used along each path
#' @param familygroup index set of the type of multivariate outcomes: "1" for Gaussian, "2" for Bernoulli, "3" for Poisson outcomes
#' @param cIndex control index, specifying index of control variable in the design matrix X
#' @param ofset offset matrix specified
#' @param family set of family [gaussian, bernoulli, possion]
#' @param control a list of internal parameters controlling the model fitting
#' @param nfold number of folds in k-fold crossvalidation
#' @param PATH TRUE/FALSE for generating solution path of sequential estimate after cross-validation step
#' @param weightU vector of inputs weights which will be combined with weights in the algorithm
#' @return
#'   \item{C}{estimated coefficient matrix; based on GIC}
#'   \item{Z}{estimated control variable coefficient matrix}
#'   \item{Phi}{estimted dispersion parameters}
#'   \item{U}{estimated U matrix (generalize latent factor weights)}
#'   \item{D}{estimated singular values}
#'   \item{V}{estimated V matrix (factor loadings)}
#'   \item{lam}{selected lambda values based on the chosen information criterion}
#'   \item{familygroup}{spcified familygroup of outcome variables.}
#'   \item{fitCV}{output from crossvalidation step, for each sequential step}
#' @export
#' @useDynLib gofar
#' @examples
#' \donttest{
#' ## Model specification:
#' SD <- 123
#' set.seed(SD)
#' n <- 200
#' p <- 100
#' pz <- 0
#' # Model I in the paper
#' # n <- 200; p <- 300; pz <- 0 ;           # Model II in the paper
#' # q1 <- 0; q2 <- 30; q3 <- 0               # Similar response cases
#' q1 <- 15
#' q2 <- 15
#' q3 <- 0 # mixed response cases
#' nrank <- 3 # true rank
#' rank.est <- 4 # estimated rank
#' nlam <- 40 # number of tuning parameter
#' s <- 1 # multiplying factor to singular value
#' snr <- 0.25 # SNR for variance Gaussian error
#' #
#' q <- q1 + q2 + q3
#' respFamily <- c("gaussian", "binomial", "poisson")
#' family <- list(gaussian(), binomial(), poisson())
#' familygroup <- c(rep(1, q1), rep(2, q2), rep(3, q3))
#' cfamily <- unique(familygroup)
#' nfamily <- length(cfamily)
#' #
#' control <- gofar_control()
#' #
#' #
#' ## Generate data
#' D <- rep(0, nrank)
#' V <- matrix(0, ncol = nrank, nrow = q)
#' U <- matrix(0, ncol = nrank, nrow = p)
#' #
#' U[, 1] <- c(sample(c(1, -1), 8, replace = TRUE), rep(0, p - 8))
#' U[, 2] <- c(rep(0, 5), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 14))
#' U[, 3] <- c(rep(0, 11), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 20))
#' #
#' if (nfamily == 1) {
#'   # for similar type response type setting
#'   V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 16))
#'   V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 28))
#'   V[, 3] <- c(
#'     sample(c(1, -1), 5, replace = TRUE) * runif(5, 0.3, 1), rep(0, 23),
#'     sample(c(1, -1), 2, replace = TRUE) * runif(2, 0.3, 1), rep(0, q - 30)
#'   )
#' } else {
#'   # for mixed type response setting
#'   # V is generated such that joint learning can be emphasised
#'   V1 <- matrix(0, ncol = nrank, nrow = q / 2)
#'   V1[, 1] <- c(sample(c(1, -1), 5, replace = TRUE), rep(0, q / 2 - 5))
#'   V1[, 2] <- c(
#'     rep(0, 3), V1[4, 1], -1 * V1[5, 1],
#'     sample(c(1, -1), 3, replace = TRUE), rep(0, q / 2 - 8)
#'   )
#'   V1[, 3] <- c(
#'     V1[1, 1], -1 * V1[2, 1], rep(0, 4),
#'     V1[7, 2], -1 * V1[8, 2], sample(c(1, -1), 2, replace = TRUE),
#'     rep(0, q / 2 - 10)
#'   )
#'   #
#'   V2 <- matrix(0, ncol = nrank, nrow = q / 2)
#'   V2[, 1] <- c(sample(c(1, -1), 5, replace = TRUE), rep(0, q / 2 - 5))
#'   V2[, 2] <- c(
#'     rep(0, 3), V2[4, 1], -1 * V2[5, 1],
#'     sample(c(1, -1), 3, replace = TRUE), rep(0, q / 2 - 8)
#'   )
#'   V2[, 3] <- c(
#'     V2[1, 1], -1 * V2[2, 1], rep(0, 4),
#'     V2[7, 2], -1 * V2[8, 2],
#'     sample(c(1, -1), 2, replace = TRUE), rep(0, q / 2 - 10)
#'   )
#'   #
#'   V <- rbind(V1, V2)
#' }
#' U[, 1:3] <- apply(U[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' V[, 1:3] <- apply(V[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' #
#' D <- s * c(4, 6, 5) # signal strength varries as per the value of s
#' or <- order(D, decreasing = T)
#' U <- U[, or]
#' V <- V[, or]
#' D <- D[or]
#' C <- U %*% (D * t(V)) # simulated coefficient matrix
#' intercept <- rep(0.5, q) # specifying intercept to the model:
#' C0 <- rbind(intercept, C)
#' #
#' Xsigma <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
#' # Simulated data
#' sim.sample <- gofar_sim(U, D, V, n, Xsigma, C0, familygroup, snr)
#' # Dispersion parameter
#' pHI <- c(rep(sim.sample$sigmaG, q1), rep(1, q2), rep(1, q3))
#' X <- sim.sample$X[1:n, ]
#' Y <- sim.sample$Y[1:n, ]
#'
#'
#'
#' # Test data
#' sim.sample <- gofar_sim(U, D, V, 1000, Xsigma, C0, familygroup, snr)
#' Xt <- sim.sample$X
#' Yt <- sim.sample$Y
#' #
#' crossprod(X %*% U)
#' apply(X, 2, norm, "2")
#' X0 <- cbind(1, X)
#' #
#' # Simulate data with 20% missing entries
#' miss <- 0.20 # Proportion of entries missing
#' t.ind <- sample.int(n * q, size = miss * n * q)
#' y <- as.vector(Y)
#' y[t.ind] <- NA
#' Ym <- matrix(y, n, q)
#' naind <- (!is.na(Ym)) + 0 # matrix(1,n,q)
#' misind <- any(naind == 0) + 0
#' #
#' # Model fitting begins:
#' control$epsilon <- 1e-7
#' control$spU <- 50 / p
#' control$spV <- 25 / q
#' control$maxit <- 1000
#'
#'
#'
#' # Model fitting: GOFAR(S) (full data)
#' set.seed(SD)
#' rank.est <- 5
#' fit.seq <- gofar_s(Y, X,
#'   nrank = rank.est, family = family,
#'   nlambda = nlam, familygroup = familygroup,
#'   control = control, nfold = 5
#' )
#'
#'
#' # Model fitting: GOFAR(S) (missing data)
#' set.seed(SD)
#' rank.est <- 5
#' fit.seq.m <- gofar_s(Ym, X,
#'   nrank = rank.est, family = family,
#'   nlambda = nlam, familygroup = familygroup,
#'   control = control, nfold = 5
#' )
#' }
#' @references
#'  Mishra, A., Dey, D., Chen, K. (2019) \emph{ Generalized co-sparse factor regression, In prepration}
gofar_s <- function(Yt, X, nrank = 3, nlambda = 40, family,
                    familygroup = NULL, cIndex = NULL, ofset = NULL,
                    control = list(), nfold = 5, PATH = FALSE, weightU = NULL) {
  cat("Initializing...", "\n")
  n <- nrow(Yt)
  p <- ncol(X)
  q <- ncol(Yt)

  if (is.null(cIndex)) {
    cIndex <- 1
    pz <- 0
    p <- ncol(X)
  } else {
    pz <- length(cIndex)
    p <- ncol(X) - pz
    cIndex <- c(1, cIndex + 1)
  }
  if (is.null(familygroup)) familygroup <- rep(1, q)

  U <- matrix(0, nrow = p, ncol = nrank)
  V <- matrix(0, nrow = q, ncol = nrank)
  D <- rep(0, nrank)
  Z <- matrix(0, nrow = pz + 1, ncol = q)
  PHI <- rep(0, q)
  lamSel <- rep(0, nrank)
  X0 <- cbind(1, X)

  xtem <- getScaleGaussian(Yt, X0, familygroup)
  mx <- xtem$ko
  Y <- xtem$Y


  # ObjDecpath <- matrix(nrow = nlambda + 1, ncol = nrank, NA)
  # totTime <- 0

  if (is.null(ofset)) ofset <- matrix(0, nrow = n, ncol = q)
  aft <- glmCol(Y, X0, ofset, family, familygroup, q, cIndex)
  Z <- aft$Z
  PHI <- aft$PHI


  ## handiling missing data
  # naind <- (!is.na(Y)) + 0 #matrix(1,n,q)
  # misind <- any(naind==0)+0
  # if(misind==1) Y[is.na(Y)] <- 0

  N.nna <- sum(!is.na(Y))
  ind.nna <- which(!is.na(Y))

  Yf2 <- Y
  naind22 <- (!is.na(Yf2)) + 0 # matrix(1,n,q)
  misind22 <- any(naind22 == 0) + 0
  Yf2[is.na(Yf2)] <- 0
  fit.nlayer <- fit.nfold <- vector("list", nrank)




  # Implementation of k-fold cross validation:
  for (k in 1:nrank) { # desired rank extraction
    cat("Initializing unit-rank unit", k, "\n")
    ## extract unpelized unit rank and estimate and use it for weight construction:
    kappaco <- getKappaC0(X0, familygroup,control$alp)
    ndev <- getNullDev(Yf2, ofset, familygroup, naind22)



    xx <- gcure_cpp_init2(Yf2, X0, 1,
                          cindex = cIndex,
                          ofset = ofset, familygroup,
                          Zini = Z, PhiIni = PHI, kappaC0 = kappaco,
                          control, misind22, naind22, ndev
    )
    # print(c(xx$maxit,xx$converge,ndev))
    XC <- X0[, -cIndex] %*% xx$C[-cIndex, ]
    Z0 <- matrix(xx$C[cIndex, ], ncol = q)
    Phi0 <- xx$PHI
    svdxc <- svd(XC)
    nkran <- sum(svdxc$d > 1e-07)
    xx$V <- as.matrix(svdxc$v[, 1:nkran])
    xx$D <- svdxc$d[1:nkran] / sqrt(n)
    xx$U <- xx$C[-cIndex, ] %*% xx$V %*% diag(1 / xx$D, nrow = nkran, ncol = nkran)
    
    
    if(is.null(weightU)){
      initW <- list(
        wu = abs(xx$U)^-control$gamma0,
        wd = abs(xx$D)^-control$gamma0,
        wv = abs(xx$V)^-control$gamma0
      )
    }
    else{
    initW <- list(
      wu = abs(xx$U)^-weightU,
      wd = abs(xx$D)^-control$gamma0,
      wv = abs(xx$V)^-control$gamma0
    )
    }
    
    lambda.max <- get.lam.max2(Y, X, familygroup, ofset)
    cat(xx$D, '\n')
    cat("Cross validation:", k, "\n")
    # save(list=ls(),file = 'aditya1.rda')

    ## store the deviance of the test data
    dev <- matrix(NA, nfold, nlambda)
    sdcal <- matrix(NA, nfold, nlambda)
    tec <- rep(0, nfold)
    ID <- rep(1:nfold, len = N.nna)
    ID <- sample(ID, N.nna, replace = FALSE)

    fitT <- vector("list", nfold)
    for (ifold in 1:nfold) { # ifold=1
      cat('Fold ', ifold, '\n')
      ind.test <- ind.nna[which(ID == ifold)]
      Yte <- Y
      Yte[-ind.test] <- NA
      Ytr <- Y
      Ytr[ind.test] <- NA
      zerosol <- glmCol(Ytr, X0, ofset, family, familygroup, q, cIndex)

      naind <- (!is.na(Ytr)) + 0 # matrix(1,n,q)
      misind <- any(naind == 0) + 0
      Ytr[is.na(Ytr)] <- 0
      ndev <- getNullDev(Ytr, ofset, familygroup, naind)

      fitT[[ifold]] <- gcure_cpp_miss(Ytr, X0,
                                      nlam = nlambda,
                                      cindex = cIndex,
                                      ofset = ofset, familygroup, initw = initW,
                                      Dini = xx$D, Zini = Z0, PhiIni = Phi0,
                                      Uini = xx$U, Vini = xx$V,
                                      kappaC0 = kappaco, lmax = lambda.max,
                                      control, misind,
                                      naind, ndev, 0, zerosol,
                                      control$maxit, control$epsilon
      )



      fitF <- fitT[[ifold]]
      # naind2 <- (!is.na(Yte)) + 0 #matrix(1,n,q)
      # Yte[is.na(Yte)] <- 0
      if (fitF$nkpath == 1) {
        im <- fitF$nkpath
        lam <- fitF$lamKpath[im, 1]
        insel <- which(fitF$lamseq == lam)
        mu.test <- fitF$mukpath[, , im]
        mu.test[-ind.test] <- NA
        tttval <- objFun5(Yte, mu.test, fitF$phipath[, im], familygroup)
        dev[ifold, ] <- tttval[1] #+ tval
        sdcal[ifold, ] <- tttval[2]
        tec[ifold] <- tttval[3]
        # dev[ifold, ] <- dev[ifold, ]/sum((!is.na(Yte)))
      } else {
        for (im in 1:fitF$nkpath) {
          lam <- fitF$lamKpath[im, 1]
          insel <- which(fitF$lamseq == lam)
          mu.test <- fitF$mukpath[, , im]
          mu.test[-ind.test] <- NA
          tval <- control$elnetAlpha * lam * fitF$dkpath[im, 1] * initW$wd * sum(initW$wu * abs(fitF$ukpath[, im])) * sum(initW$wv * abs(fitF$vkpath[, im]))
          abc <- sum((fitF$vkpath[, im])^2) * sum((fitF$ukpath[, im])^2) * (fitF$dkpath[im, 1])^2
          tval <- tval + lam * (1 - control$elnetAlpha) * abc

          tttval <- objFun5(Yte, mu.test, fitF$phipath[, im], familygroup)
          dev[ifold, insel] <- tttval[1] #+ tval
          sdcal[ifold, insel] <- tttval[2]
          tec[ifold] <- tttval[3]
        }
      }
    }

    #save(list = ls(), file = "aditya2.rda")


    fit.nfold[[k]] <- fitT
    dev.mean <- colMeans(dev, na.rm = FALSE)
    sderr <- sqrt(apply(sdcal, 2, sum, na.rm = FALSE) / (nfold * sum(tec - 1)))
    l.mean <- which.min(dev.mean)
    if (is.na(sderr[l.mean])) {
      l.mean <- min(which(dev.mean <= (dev.mean[l.mean])))
    } else {
      l.mean <- min(which(dev.mean <= (dev.mean[l.mean] +
                                         control$se1 * sderr[l.mean])))
    }
    lamS <- fitF$lamseq[l.mean]


    Yf <- Y
    zerosol <- glmCol(Yf, X0, ofset, family, familygroup, q, cIndex)

    naind <- (!is.na(Yf)) + 0 # matrix(1,n,q)
    misind <- any(naind == 0) + 0
    Yf[is.na(Yf)] <- 0
    ndev <- getNullDev(Yf, ofset, familygroup, naind)

    if ( PATH ) {
      fit.nlayer[[k]] <- fit.layer <- gcure_cpp_miss(Yf,X0,
                                                     nlam = nlambda,
                                                     cindex = cIndex,
                                                     ofset = ofset,
                                                     familygroup,
                                                     initw = initW,
                                                     Dini = xx$D,
                                                     Zini = Z0 ,
                                                     PhiIni = Phi0,
                                                     Uini = xx$U,
                                                     Vini = xx$V,
                                                     kappaC0 = kappaco,
                                                     lmax = lambda.max,
                                                     control,misind,
                                                     naind,ndev,0,zerosol,
                                                     control$maxit,
                                                     control$epsilon)
      l.mean <- which(fit.layer$lamKpath == lamS)
      if ( length(l.mean) == 0) {
        l.mean <- 1;
      }
    } else {
      control1 <- control
      control1$lamMinFac <- 1
      fit.nlayer[[k]] <- fit.layer <- gcure_cpp_miss(Yf, X0,
                                                     nlam = 1,
                                                     cindex = cIndex,
                                                     ofset = ofset,
                                                     familygroup,
                                                     initw = initW,
                                                     Dini = xx$D,
                                                     Zini = Z0,
                                                     PhiIni = Phi0,
                                                     Uini = xx$U,
                                                     Vini = xx$V,
                                                     kappaC0 = kappaco,
                                                     # lmax = lambda.max,
                                                     lmax = lamS,
                                                     control1, misind,
                                                     naind, ndev, 0,
                                                     zerosol,
                                                     control$maxit,
                                                     control$epsilon)
      l.mean <- 1
    }

    ## Same as before
    fit.nlayer[[k]]$initw <- initW
    fit.nlayer[[k]]$lamseq <- fitF$lamseq   ## save sequence of lambda path
    fit.nlayer[[k]]$dev <- dev
    fit.nlayer[[k]]$lamS <- lamS
    fit.nlayer[[k]]$sderr <- sderr
    U[, k] <- fit.layer$ukpath[, l.mean]
    V[, k] <- fit.layer$vkpath[, l.mean]
    D[k] <- fit.layer$dkpath[l.mean, 1]
    cat(D[k], '\n')
    lamSel[k] <- fit.layer$lamKpath[l.mean, 1]

    if (D[k] == 0) {
      U <- matrix(U[, 1:(k - 1)], ncol = k - 1)
      V <- matrix(V[, 1:(k - 1)], ncol = k - 1)
      D <- D[1:(k - 1)]
      lamSel <- lamSel[1:(k - 1)]
      break
    }

    Ck <- D[k] * tcrossprod(U[, k], V[, k])
    PHI <- fit.layer$phipath[, l.mean]
    if (pz == 0) {
      Z[1, ] <-  drop(fit.layer$zpath[, ,l.mean])
    } else {
      Z <- fit.layer$zpath[, , l.mean]
    }
    ofset <- ofset + crossprod(t(X0[, -cIndex]), Ck) # crossprod(t(X), Ck)
  }

  cat("Estimated rank =", sum(D != 0), "\n")
  if (sum(D != 0) == nrank) {
    cat("Increase nrank value!")
  }
  ft1 <- list(fit = fit.nlayer, C = U %*% (D * t(V)), Z = Z, Phi = PHI,
    U = U, V = V, D = D,
    lam = lamSel, familygroup = familygroup
  )
  if (all(unique(familygroup) == 1:2)) {
    ft1 <- updateFitObject(ft1, mx)
  }
  return(ft1)
}




