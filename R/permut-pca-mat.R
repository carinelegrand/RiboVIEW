permut.pca.mat <- function(pourPca, pca.bco.pc.var, B = 10000, seed.b=300, 
                           permut.byrow = TRUE, bootstrap = TRUE) {

  # In :
  #   pourPca              Matrix for resampling and PCA
  #   pca.bco.pc.var       Percentage of variance explained by PC, in original matrix
  #   B                    Number of resampling or bootstrap iterations
  #   seed.b               Seed for reproducible random number generation
  #   permut.byrow         If TRUE, each row is resampled, otherwise the full matrix
  #   bootstrap = TRUE     If TRUE, resampling is done with replacement
  #
  # Out : 
  #   pca.bco.pc.var       Input / Output, for tidy matrix of results ; 
  #                        Percentage of variance explained by PC, in original matrix
  #   pca.bco.pc.var.p
  #   pca.bco.pc.var.p.c
  
  set.seed(seed.b)
  eig <- c()
  nr <- nrow(pourPca)
  nc <- ncol(pourPca)

  for (k in 1:B) {
      # Resample each row (preserves part of potential batch effects),
      #   or resample full matrix
      # If bootstrap=TRUE, resampling is done with replacement.
      if (permut.byrow) {
          pourPca.k <- pourPca[1,]
          for (ir in 2:nr) {
              pourPca.k <- rbind(pourPca.k,
                                  sample(pourPca[ir,], size = nc, replace = bootstrap))
          }
      } else {
          pourPca.k <- matrix(data = sample(pourPca, size=nc*nr, replace=bootstrap),
                               nrow = nr,
                               ncol = nc)
      }
      # Singular value decomposition of kth resampled matrix
      pca.bco.pre.k <- svd(stats::cov(pourPca.k)) #cor : seems not to make any difference on eig, hist(eig)
      pca.bco.k <- pca.bco.pre.k$u
      pca.bco.var.k <- sqrt(pca.bco.pre.k$d)
      pca.bco.pc.var.k <- 100 * pca.bco.var.k / sum(pca.bco.var.k)
      # Eigenvalues of kth resampled matrix stored for later evaluation
      eig[[k]] <- pca.bco.pc.var.k
  }
  # List of eigenvalues, 
  #   one row per resampling, 
  #   columns correspond to principal component 1, principal component 2, etc.
  eig.mat <- matrix(unlist(eig), nrow=B, ncol=nc, byrow = TRUE)
  
  # Init
  pca.bco.pc.var.p <- c()
  pca.bco.pc.var.p.c <- c()
  
  # Loop on principal component 1 to nc
  for (i in 1:nc) {
    # Proportion (in [0;1]) of eigenvalues of PCi larger than the original 
    #   eigenvalue of PCi, aka bootstrap p-value.
    pB.pre <- round(mean(eig.mat[,i] > pca.bco.pc.var[i]),4)
    # Based on B=10000 by default, avoid claiming wrongly that a p-value is 0 or 1.
    if (pB.pre < 0.0001) {
      pB <- "< 1e-4"
    } else if (pB.pre > 0.9999) {
      pB <- "> 0.9999"
    } else {
    pB <- paste("= ",as.character(pB.pre), sep="")

    }
    # p-values
    pca.bco.pc.var.p <- c(pca.bco.pc.var.p, pB.pre)
    # Curated p-values
    pca.bco.pc.var.p.c <- c(pca.bco.pc.var.p.c, pB)
  }
  
  # Percentage of variance explained and associated p.value and curated p-value
  rbind(pca.bco.pc.var,
        pca.bco.pc.var.p,
        pca.bco.pc.var.p.c)
}  
