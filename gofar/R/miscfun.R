# miscelaneous function



#' Obtain scaling constant for monotone decreasing G-CURE problem
#'
#' @param X0 Design matrix
#' @param familygroup indicator {1,2,3} for {gaussian, binary, poisson}
#' @useDynLib gofar
#' @import magrittr
#' @importFrom stats family gaussian binomial poisson
getKappaC0 = function(X0,familygroup){
  family <- list(gaussian(),binomial(),poisson())
  temp <- rep(0,length(family))
  svdX0d1 <- svd(X0)$d[1]
  cfamily <- unique(familygroup)
  for (j in cfamily)
    temp[j] <- switch(family[[j]]$family,'gaussian' = svdX0d1,
                      'binomial' = svdX0d1/2,'poisson' = svdX0d1*10)
  1*max(temp)
}




#' Obtain lambda_max for for generating sequence of tuning parameter in G-CURE problem
#'
#' @param X Design matrix
#' @param Y Multivariate response matrix
#' @param familygroup indicator {1,2,3} for {gaussian, binary, poisson}
#' @param offset offset term matrix
#' @useDynLib gofar
#' @import magrittr
#' @importFrom stats family gaussian binomial poisson
get.lam.max2 = function(Y, X,familygroup,offset){
  Y[is.na(Y)] <- 0
  family <- list(gaussian(), binomial(),poisson())
  cfamily <- unique(familygroup)
  n <- nrow(X);q <- ncol(Y)
  ofset.ini <- offset # matrix(0,nrow=n,ncol=q)
  MU0 <- matrix(0,nrow=n,ncol=q)
  for(i in cfamily)
    MU0[,familygroup == i] <- family[[i]]$linkinv(ofset.ini[,familygroup == i])
  ttt <- abs(crossprod(X,Y-MU0));maxtt <- max(ttt)
  lambda.max <- 2*maxtt
  return(lambda.max)
}







#' Obtain null deviance to utilize in convergence
#'
#' @param Y Multivariate response matrix
#' @param familygroup indicator {1,2,3} for {gaussian, binary, poisson}
#' @param ofset offset matrix
#' @param naind index of missing entries
#' @useDynLib gofar
#' @import magrittr
#' @importFrom stats family gaussian binomial poisson
#' @importFrom stats coef sd glm.control dbinom dnorm dpois glm
getNullDev <- function(Y, ofset, familygroup, naind) {
  dev <- 0
  t1 <- which(familygroup == 1)
  t2 <- which(familygroup == 2)
  t3 <- which(familygroup == 3)
  if (length(t1) > 0) {
    abc <- scale(as.matrix(Y[, t1]) - as.matrix(ofset[, t1]), center = T, scale = F)
    # class(abc)
    dev <- dev + sum((abc * as.matrix(naind[, t1]))^2)
    # print(c(dev,sum(ofset)))
  }
  if (length(t2) > 0) {
    m1 <- as.matrix(Y[, t2])
    m2 <- as.matrix(ofset[, t2])
    m3 <- as.matrix(naind[, t2])

    for (i in 1:length(t2)) {
      qq <- m3[, i] == 1
      glm.fitxx <- glm(m1[qq, i] ~ offset(m2[qq, i]), family = binomial)
      abc <- as.numeric(coef(glm.fitxx))
      p2 <- exp(m2[qq, i] + abc) / (1 + exp(m2[qq, i] + abc))
      dev <- dev - 2 * sum(m1[qq, i] * log(p2) + (1 - m1[qq, i]) * log(1 - p2))
      # print(dev)
    }
  }
  if (length(t3) > 0) {
    qq <- naind[, t3] == 1
    abc <- as.matrix(exp(ofset[, t3]))
    abc[!qq] <- 0
    li <- abc %*% diag(colSums(as.matrix(Y[, t3])) / colSums(abc), length(t3))
    abc <- log(as.matrix(Y[, t3])) - log(li)
    abc[!is.finite(abc)] <- 0
    dev <- dev + sum(2 * (as.matrix(Y[, t3]) * abc + li - as.matrix(Y[, t3])))
    # print(dev)
  }
  dev # /sum(naind)
}






#' Fit glm Columnwise on the control variable
#'
#' @param Y Multivariate response matrix
#' @param X0 design matrix
#' @param familygroup indicator {1,2,3} for {gaussian, binary, poisson}
#' @param ofset offset  matrix
#' @param family {gaussian, Bernouli, poisson}
#' @param q number of outcomes
#' @param cIndex index of control variable in X0
#' @useDynLib gofar
#' @import magrittr
#' @importFrom stats family gaussian binomial poisson df.residual residuals
#' @importFrom stats coef sd glm.control dbinom dnorm dpois glm glm.fit
glmCol <- function(Y, X0, ofset, family, familygroup, q, cIndex) {
  PHI <- rep(1, q)
  Z <- matrix(0, nrow = length(cIndex), ncol = q)
  for (i in 1:q) {
    qqq <- !is.na(Y[, i])
    ft <- glm.fit(X0[qqq, cIndex], Y[qqq, i],
                  family = family[[familygroup[i]]],
                  offset = ofset[qqq, i], intercept = F,
                  control = glm.control(maxit = 5000)
    )
    Z[, i] <- ft$coefficients
    PHI[i] <- ifelse(familygroup[i] == 1,
                     sum(residuals(ft)^2) / df.residual(ft), 1
    )
    if (PHI[i] == 0) PHI[i] <- 1
  }
  return(list(Z = Z, PHI = PHI))
}


#' Suitably scale gaussian response for unit variance
#'
#' @param Y0 Multivariate response matrix
#' @param X0 design matrix
#' @param familygroup indicator {1,2,3} for {gaussian, binary, poisson}
#' @useDynLib gofar
#' @import magrittr
#' @importFrom stats family gaussian binomial poisson
#' @importFrom stats coef sd glm.control dbinom dnorm dpois glm glm.fit
#' @importFrom rrpack rrr
getScaleGaussian1 = function(Y0, X0, familygroup) {
  if (all(unique(familygroup) == 1:2)) {
    Y <- Y0
    Y[is.na(Y)] <- 0
    ft <- rrpack::rrr(Y[, familygroup == 1], X0, "rank", ic.type = "BIC", maxrank = 10)
    res <- Y0[, familygroup == 1] - X0 %*% ft$coef
    mx <- min(apply(res, 2, sd, na.rm = T))
    Y2 <- Y0
    Y2[, familygroup == 1] <- Y2[, familygroup == 1] / mx
    return(list(ko = mx, Y = Y2))
  } else {
    return(list(ko = 1, Y = Y0))
  }
}



#' Suitably scale gaussian response for unit variance
#'
#' @param Y0 Multivariate response matrix
#' @param X0 design matrix
#' @param familygroup indicator {1,2,3} for {gaussian, binary, poisson}
#' @useDynLib gofar
#' @import magrittr
#' @importFrom stats family gaussian binomial poisson
#' @importFrom stats coef sd glm.control dbinom dnorm dpois glm glm.fit
#' @importFrom glmnet cv.glmnet
getScaleGaussian = function(Yt, X0, familygroup) {
  if (all(unique(familygroup) == 1:2)) {
    Ys = Yt[, familygroup == 1]
    mx <- rep(0,ncol(Ys))
    for (i in 1:ncol(Ys)) {
      qqq <- !is.na(Ys[, i])
      cv <- glmnet::cv.glmnet(X0[qqq,-1],Ys[qqq,i],family="gaussian",
                              standardize = F,intercept=T,
                              nfold=10,nlambda = 20,
                              alpha=1, maxit = 10000)
      mx[i] <- sd(Ys[,i] - X0%*%as.vector(coef(cv)))
    }
    mx = min(mx)
    Yt[, familygroup == 1] <- Yt[, familygroup == 1] / mx
    return(list(ko = mx, Y = Yt))
  } else {
    return(list(ko = 1, Y = Yt))
  }
}


#' Rescale gaussian response
#'
#' @param fit fitted object from gsecure and geecure
#' @param mx scaling value
#' @useDynLib gofar
#' @import magrittr
#' @importFrom stats coef sd glm.control dbinom dnorm dpois glm
updateFitObject <- function(fit, mx) {
  # fit=fit.seq2
  vx <- fit$V
  vx[fit$familygroup == 1, ] <- vx[fit$familygroup == 1, ] * mx
  fit$V <- t(t(vx) / sqrt(colSums((vx)^2)))
  fit$D <- fit$D * sqrt(colSums((vx)^2))
  fit$Phi[fit$familygroup == 1] <- (mx^2) * fit$Phi[fit$familygroup == 1]
  fit$Z[, fit$familygroup == 1] <- mx * fit$Z[, fit$familygroup == 1]
  fit$C <- fit$U %*% (fit$D * t(fit$V))
  return(fit)
}



#' loglikelihood of the observation
#'
#' @param Y outcome variables
#' @param MU natural parameter matrix
#' @param Sigma dispersion parameter for gacussian
#' @param family gaussian binomial poisson
#' @importFrom stats coef sd glm.control dbinom dnorm dpois glm
logLikehood <- function(Y, MU, Sigma = 1, family) {
  switch(family,
         "gaussian" = dnorm(Y, MU, Sigma, log = TRUE),
         "poisson" = dpois(Y, MU, log = TRUE),
         "binomial" = dbinom(Y, 1, MU, log = TRUE)
  )
}





#' Evaluate of objective function
#'
#' @param Y outcome variables
#' @param mu natural parameter matrix
#' @param Phi dispersion parameter
#' @param familygroup gaussian binomial poisson
#' @importFrom stats family gaussian binomial poisson
#' @importFrom stats coef sd glm.control dbinom dnorm dpois glm
objFun5 <- function(Y, mu, Phi, familygroup) {
  n <- nrow(Y)
  family <- list(gaussian(), binomial(), poisson())
  cfamily <- unique(familygroup)
  negLogl <- 0 * Y
  for (j in cfamily) {
    # print(Phi[familygroup == j])
    # sum(familygroup==j)
    Sig <- if (j == 1) t(matrix(Phi[familygroup == j], sum(familygroup == j), n)) else 1
    negLogl[, familygroup == j] <- -1 * logLikehood(
      Y[, familygroup == j],
      mu[, familygroup == j],
      sqrt(Sig), family[[j]]$family
    )
  }
  kind <- !is.na(Y)
  teS <- sum(kind)
  mn <- sum(negLogl[kind]) / teS
  sdp <- (sum((negLogl[kind] - mn)^2)) / teS
  c(mn, sdp, teS)
}

