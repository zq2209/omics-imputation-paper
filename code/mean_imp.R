mean_imputation <- function(x) {
  n <- nrow(x)
  xmis <- x
  for (t.row in 1:n) {
    xmis[t.row, is.na(x[t.row,])] <- rowMeans(x, na.rm = TRUE)[t.row] #mean(xmis[,t.co], na.rm = TRUE)
    #if (is.numeric(xmis[[t.co]])) {
    # varType[t.co] <- 'numeric'
    #ximp[is.na(xmis[,t.co]),t.co] <- mean(xmis[,t.co], na.rm = TRUE)
    #next()
    #}
  }
  return(xmis)
}