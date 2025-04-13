
GraphicalModel <- function(Adj, a1 = -0.9, a2 = -0.1, b1 = 0.1, b2 = 0.9, scaleP = 1.5) {
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  library(MASS)
  ## Check adjacency matrix
  if (nrow(Adj) != ncol(Adj)) {
    stop("Error: The adjacency matrix must be squared matrix!")
  } else {
    if (!all.equal(Adj, t(Adj))) {
      stop("Error: The adjacency matrix must be symmetric matrix!")
    } else if (unique(diag(Adj)) != 0) {
      stop("Error: The diagnal elements of adjacency matrix must be all 0!")
    } else if (!all.equal(unique(as.vector(Adj)), c(0, 1))) {
      stop("Error: The elements of adjacency matrix must be 0 or 1!")
    }
  }
  
  ## Calculate Sigma
  n <- nrow(Adj) * ncol(Adj)
  rho <- runif(n, 0, a2 - a1 + b2 - b1)
  rho.1 <- (a1 + rho) * (rho < (a2 - a1)) + (b1 + rho - (a2 - a1)) * (rho >= (a2 - a1))
  Adj <- Adj * matrix(rho.1, nrow(Adj))
  # - Rescale matrix
  rSum <- rowSums(abs(Adj))
  Adj.rescale <- Adj / (scaleP * sapply(1:ncol(Adj), function(t) {
    return(rSum)
  }))
  # - Ensure symmetry
  A <- 0.5 * (Adj.rescale + t(Adj.rescale))
  diag(A) <- 1
  # - Calculate Sigma
  A.inv <- ginv(A)
  A.diag <- as.matrix(diag(A.inv))
  Sigma <- A.inv / sqrt(A.diag %*% t(A.diag))
  return(Sigma)
}




simulateX <- function(y, effect, gve) {
  #y <- matrix(Y[1, 5:ncol(Y)], nrow = 1)
  n.module <- sample(5:15, 1)
  n <- sapply(1:n.module, function(x) {return(sample(10:50, 1))})
  a <- lapply(1:n.module, function(t) {
    temp <- barabasi.game(n[t], directed = FALSE)
    return(data.matrix(as_adjacency_matrix(temp, type = "both")))
  })
  Adj <- data.matrix(bdiag(a)) # need modify
  
  Sigma <- tryCatch(GraphicalModel(Adj, scaleP = 1.5),error = function(e) e, warning = function(w) w)
  temp <- 1.5
  while (is(Sigma,"warning")) {
    print('Retry')
    temp <- temp + 0.05
    Sigma <- tryCatch(GraphicalModel(Adj, scaleP = temp),error = function(e) e, warning = function(w) w)
  }
  
  n_degree <- rowSums(Adj)
  pos <- sample(n.module, round(runif(1, 1, 5)), replace = FALSE) # how many blocks have true effect
  temp <- c(0, cumsum(n))
  true_causal <- sapply(pos, function(t){return(sample(c((temp[t]+1):temp[t+1]), 1))}) # every block with true effect only has 1 true causal variant
  # true_causal <- sapply(pos, function(t){return(sample((temp[t]+1):temp[t+1]), 5)})
  beta <- rep(0, nrow(Adj))
  #s <- 10 # set s that ensure complete y and x has a power of 0.8
  # b <- lapply(true_causal, function(t){return(sqrt(n_degree[t]) * rnorm(1, 0, s))})
  b <- lapply(true_causal, function(t){return(effect)}) # effect size
  beta[unlist(true_causal)] <- unlist(b)
  
  beta <- as.matrix(beta)
  N_sample <- length(y)
  Z <- rmvnorm(N_sample, mean = rep(0, nrow(Adj)), sigma = Sigma)
  #phi <- gve
  y <- unlist(y)
  var_y <- var(y)
  sigma_error <- var_y*effect^2/gve - var_y*effect^2#mean(sapply(as.vector(true_causal), function(x){((1-gve)*var_y - (gve / beta[x]^2))/(gve/beta[x]^2)})) 
  #sigma_error <- length(as.vector(true_causal))*((1-phi)*beta[x]^2*var_y + Sigma[x, x])
  error <- rmvnorm(N_sample, mean = rep(0, nrow(Adj)), sigma = diag(rep(sigma_error, nrow(Adj))))
  #gve_lst1 <- c(gve_lst1, var_y*effect^2/(var_y*effect^2 + mean(diag(cov(error)))))
  #X <- y %*% t(beta) + error
  X <- as.matrix(y) %*% as.matrix(t(beta)) + Z + error
  return(list('X' = X, 'true_causal' = true_causal))
}






