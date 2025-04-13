# Calculate three comparison metrics for imputation methods

calculate_metrics <- function(df_to_imp, imput_df, group_vars = "method") {
  # Get SD and range from original data
  phenotype <- df_to_imp
  pheno <- phenotype
  samples <- colnames(pheno)[-(1:4)]
  pheno[, pct.na := rowSums(is.na(.SD))/length(samples), .SDcols = samples]
  pheno <- pheno[pct.na < 0.5]

  pheno[, sd := apply(.SD, 1, sd, na.rm = TRUE), .SDcols = samples]
  pheno[, range := apply(.SD, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)), .SDcols = samples]
  rownames(pheno) <- pheno$ID

  # Calculate errors
  imput_df[, error := imputation - true]

  # Global metrics
  imput_df[, mae := mean(abs(error)), by = group_vars]
  imput_df[, mse := mean(error^2), by = group_vars]
  imput_df[, nmae := mae/(max(true) - min(true)), by = group_vars]
  imput_df[, nmse := mse/var(true), by = group_vars]
  imput_df[, r2 := cor(imputation, true)^2, by = group_vars]

  # Per-feature metrics
  imput_df[, mae_feature := mean(abs(error)), by = c('ID', group_vars)]
  imput_df[, mse_feature := mean(error^2), by = c('ID', group_vars)]
  imput_df[, r2_feature := cor(imputation, true)^2, by = c('ID', group_vars)]
  imput_df[, nmae_feature := mean(abs(error))/(max(true) - min(true)), by = c('ID', group_vars)]
  imput_df[, nmse_feature := mean(error^2)/var(true), by = c('ID', group_vars)]

  # Merge with original data stats
  res_features <- unique(imput_df, by = c('ID', group_vars))
  res_features <- merge(res_features, pheno[, .(ID, sd, range)], by = 'ID')

  # Calculate normalized metrics
  res_features[, NRMSE := sqrt(mse_feature)/sd]
  res_features[, NMAE := mae_feature/range]

  return(res_features)
}

