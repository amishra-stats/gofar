#' Plot heatmap of GOFAR output V
#' 
#' @param fit.seq gofar fit object
#' @param Y qxn multivariate output data matrix in GOFAR with annotated rows and columns
#' @return heatmap of output vectors v_1,..., v_r
#' 
#' library(pheatmap)
#' library(RColorBrewer)


plot.gofar.V <- function(fit.seq, Y){
  
  
  D <- fit.seq$D
  order_D <- order(D, decreasing = T)
  D <- D[order_D]
  
  V <- fit.seq$V
  colnames(V) <- paste0("v", 1:ncol(V))
  rownames(V) <- rownames(Y)
  Vplt <- V
  Vplt[which(Vplt==0)] <- NA
  pV <- pheatmap(Vplt#[, order_D]
                 , cluster_rows = F, cluster_cols = F,
                 breaks = seq(from = -max(abs(V)), to = max(abs(V)), length.out = 100),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                 cellwidth = 20, na_col = "grey")
  
}  

# library(pheatmap)
# library(RColorBrewer)
# colnames(Y) <- paste0("Y", 1:ncol(Y))
# plot.gofar.V(fit.seq = fit.seq, Y = t(Y))


#' Plot heatmap of GOFAR output U
#' 
#' @param fit.seq gofar fit object (e.g. GOFAR with additional U weights)
#' @param fit.seq.2 optional second gofar fit object (e.g. GOFAR without additional U weights)
#' @param X pxn input data matrix in GOFAR with annotated rows and columns
#' @param complex optional px1 vector of additional information on variables on the U side (e.g. complex mebmership of a protein)
#' @param own.weights optional px1 vector of input weights in GOFAR
#' @return (complex) heatmap of output vectors u_1,..., u_r
#' 
#' library(randomcoloR)
plot.gofar.U <- function(fit.seq, fit.seq.2 = NULL, X, complex = NULL,
                         own.weights = NULL){
  
  D <- fit.seq$D
  order_D <- order(D, decreasing = T)
  
  U <- fit.seq$U
  rownames(U) <- rownames(X)
  ind <- c()
  for(j in 1:ncol(U)){
    Ui <- U[, j]
    ind <- c(ind, which(Ui != 0))
  }
  ind1 <- unique(ind)
  
  if(!is.null(fit.seq.2)){
    U2 <- fit.seq.2$U
    rownames(U2) <- rownames(X)
    ind2 <- c()
    for(j in 1:ncol(U2)){
      Ui2 <- U2[, j]
      ind2 <- c(ind2, which(Ui2 != 0))
    }
    ind2 <- unique(ind2)
    
    ind <- sort(unique(c(ind1, ind2)))
    
    overlap <- as.numeric(!(ind %in% ind1))
  }
  
  else{ind <- ind1}
  
  tabU <- U[ind, ]
  colnames(tabU) <- paste0("U", 1:ncol(U))
  tabU <- round(tabU, 4)
  
  Uplt <- tabU
  Uplt[which(Uplt==0)] <- NA
  
  
  # create complex heatmap with additional information if provided in input  
  annotdf <- data.frame(row.names = rownames(Uplt))
  
  mycolors <- list()
  
  if(!is.null(complex)){ 
    annotdf$complex <- complex[ind]
    set.seed(123)
    n <- length(unique(complex))
    mycolors1 <- distinctColorPalette(n)
    
    names(mycolors1) <- unique(complex)
    mycolors1 <- mycolors1[unique(annotdf$complex)]
    mycolors$complex <- mycolors1
  }
  
  if(!is.null(own.weights)){ 
    annotdf$encode <- own.weights$encode_count[ind]
    annotdf$encode[which(annotdf$encode == 0)] <- NA
    mycolors2 <- colorRampPalette((brewer.pal(n = 7, name = "Reds")))(
      length(unique(annotdf$encode)))
    names(mycolors2) <- unique(annotdf$encode)
    mycolors$encode <- mycolors2
  }
  
  if(!is.null(fit.seq.2)){ 
    annotdf$overlap <- overlap
    mycolors3 <- c("white", "black")
    names(mycolors3) <- unique(annotdf$overlap)
    mycolors$overlap <- mycolors3
  }
  
  
  for(i in 1:ncol(Uplt)){Uplt[, i] <- Uplt[, i] * D[i]}
  
  colnames(Uplt) <- paste0(colnames(Uplt), "d", 1:ncol(Uplt))
  
  if(length(mycolors) == 0){
    pU <- pheatmap(t(Uplt), cluster_rows = F, cluster_cols = F, fontsize_row = 5,
                   fontsize_col = 5, cellwidth = 5, cellheight = 108, 
                   breaks = seq(from = -max(abs(tabU)), to = max(abs(tabU)), length.out = 100),
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                   fontsize = 4, na_col = "grey"
    )
  }
  else{
    pU <- pheatmap(t(Uplt), cluster_rows = F, cluster_cols = F, fontsize_row = 5,
                   fontsize_col = 5, cellwidth = 5, cellheight = 108, 
                   breaks = seq(from = -max(abs(tabU)), to = max(abs(tabU)), length.out = 100),
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                   annotation_col = annotdf, annotation_colors = mycolors,
                   fontsize = 4, na_col = "grey"
    )
    
  }
  
}

# library(randomcoloR)
# colnames(X) <- paste0("X", sprintf("%02d", 1:ncol(X)))
# plot.gofar.U(fit.seq = fit.seq, X = t(X))
# p <- ncol(X)
# complex.X <- sample(1:10, p, replace = T)
# plot.gofar.U(fit.seq = fit.seq, X = t(X), complex = complex.X)




#' Plot heatmap of GOFAR output coefficient matrix C (C_1,...,C_r)
#' 
#' @param fit.seq gofar fit object (e.g. GOFAR with additional U weights)
#' @param X pxn input data matrix in GOFAR with annotated rows and columns
#' @param Y qxn multivariate output data matrix in GOFAR with annotated rows and columns
#' @param complex optional px1 vector of additional information on variables on the U side (e.g. complex mebmership of a protein)
#' @return (complex) heatmap of estimated coefficient matrices C_1,..., C_r without rows and columns equal to zero

plot.gofar.C <- function(fit.seq, X, Y, complex, cellwidth = 7, cellheight = 16){
  
  D <- fit.seq$D
  V <- fit.seq$V
  U <- fit.seq$U
  rownames(U) <- rownames(X)
  
  for(i in 1:ncol(U)){
    
    Ui <- U[, i]
    ind <- which(Ui != 0)
    Ui <- Ui[ind]
    # reconstruction:
    rec <- V[, i] %*% t(Ui) * D[i]
    rownames(rec) <- rownames(Y)
    colnames(rec) <- rownames(X)[ind]
    
    
    if(!is.null(complex)){
    annotdf <- data.frame(
      row.names = rownames(U)[ind],
      complex = complex[ind]
    )
    
    newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotdf$complex))))
    mycolors <- newCols(length(unique(annotdf$complex)))
    names(mycolors) <- unique(annotdf$complex)
    mycolors <- list(complex = mycolors)
    
    
    recplt <- rec
    recplt[which(recplt == 0)] <- NA
    pheatmap(recplt, cluster_cols = F, cluster_rows = F,
             cellwidth = cellwidth, cellheight = cellheight,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(
               100), fontsize_col = 4, fontsize_row = 6,
             breaks = seq(from = -max(abs(rec)), to = max(abs(rec)), 
                          length.out = 100),
             annotation_col = annotdf, annotation_colors = mycolors,
             fontsize = 5, na_col = "grey")
    }
    else{
      recplt <- rec
      recplt[which(recplt == 0)] <- NA
      pheatmap(recplt, cluster_cols = F, cluster_rows = F,
               cellwidth = cellwidth, cellheight = cellheight,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(
                 100), fontsize_col = 4, fontsize_row = 6,
               breaks = seq(from = -max(abs(rec)), to = max(abs(rec)), 
                            length.out = 100),
               fontsize = 5, na_col = "grey")
    }
  }
}

# plot.gofar.C(fit.seq = fit.seq, X = t(X), Y = t(Y), complex = complex.X,
#              cellheight = 9)



#' This function only works for binary output right now!!!

library(gtools)

plot.gofar.reconstruction <- function(fit.seq, X, Y, linkfnc = inv.logit
                                        ){
  
  C <- fit.seq$C
  XC <- t(X) %*% C # nxq
  Zbeta <- rep(1, nrow(XC)) %*% fit.seq$Z # nxq
  
  rec <- t(Zbeta + XC)
  rownames(rec) <- rownames(Y)
  colnames(rec) <- colnames(Y) 
  Rec <- linkfnc(rec) 

  pheatmap(Rec, cluster_rows = F, cluster_cols = F,
           color = colorRampPalette((brewer.pal(n = 7, name = "Greys")))(100),
           main = expression(paste("Reconstruction:", "logit"^-1,
                                   (Z * beta + X * C)  )))
  
  
  
}

#' This function only works for binary output right now!!!

plot.gofar.residuals <- function(fit.seq, X, Y){
  
  C <- fit.seq$C
  XC <- t(X) %*% C # nxq
  Zbeta <- rep(1, nrow(XC)) %*% fit.seq$Z # nxq
  
  rec <- t(Zbeta + XC)
  rownames(rec) <- rownames(Y)
  colnames(rec) <- colnames(Y) 
  Rec <- inv.logit(rec) 
  Y_flipped <- (Y == 0)
  Resid <- -2 * log(abs(Rec - Y_flipped))
  
  pheatmap(Resid, cluster_rows = F, cluster_cols = F,
           color = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100),
           main = expression(paste("-2log", (abs(hat(y) - y^c)) )))
  
  ## The (squared) deviance of each data point is equal to (-2 times) the 
  ## logarithm of the difference between its predicted probability logit−1(𝑋𝛽)
  ## and the complement of its actual value (1 for a control; a 0 for a case) 
  ## in absolute terms. A perfect fit of a point (which never occurs) gives a 
  ## deviance of zero as log(1) is zero. A poorly fitting point has a large 
  ## residual deviance as -2 times the log of a very small value is a large number.
}
