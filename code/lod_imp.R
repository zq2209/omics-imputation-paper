
lod_imputation <- function(x) {
  n <- nrow(x)
  xmis <- x
  for (t.row in 1:n) {
    x[t.row, is.na(x[t.row,])] <- min(x[t.row, ], na.rm = TRUE)
  }
  return(xmis)
}
