library(softImpute)

soft_imputation <- function(x) {
  # softImputation
  X_mis_C = as(as.matrix(x),"Incomplete")
  ###uses "svd" algorithm
  fit1 = softImpute(X_mis_C,rank=50,lambda=30,type="svd")
  X_imp = complete(as.matrix(x),fit1)
  return(X_imp)
}
