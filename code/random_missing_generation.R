# This function creates an matrix with additional missingness. The input
# should be a matrix in which rows are features and columns are
# samples.

add_missing_values <- function(df, percent_missing = 0.5, missing) {
  n_rows <- nrow(df)
  n_cols <- ncol(df)

  # Calculate the number of missing values to add
  num_missing <- round(n_cols * percent_missing)
  num_add <- round(num_missing - missing*n_cols)

  # Loop through each row and randomly assign missing values
  for (i in 1:n_rows) {
    # Sample indices of columns to set as missing
    missing_cols <- sample(1:n_cols, num_add[i], replace = FALSE)

    # Set missing values in the selected columns for the current row
    df[i, missing_cols] <- NA
  }

  return(df)
}
