library(impute)

knn_imputation = function(x, k = 20) {
  x_imp = impute.knn(data = as.matrix(x), k=k, rowmax = 1)
  return(x_imp$data)
}