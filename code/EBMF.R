
deterK<-function(dat,niter=3){
  require('softImpute')
  require('irlba')
  pvelist <- numeric(niter)
  
  message('Automatic determination of K based on sample permutation')
  if(any(is.na(dat))){
    imputed_dat <- softImpute::complete(dat, softImpute::softImpute(dat))
  }else{
    imputed_dat <- dat
    
  }
  
  imputed_dat <- imputed_dat / apply(imputed_dat, 1, sd) # Scale feature-wise.
  total_var <- sum(imputed_dat^2)
  # No additional factors were added:
  
  for (i in 1:niter) {
    message("Iteration ", i,"/",niter)
    set.seed(i)
    permdat <- t(apply(imputed_dat, 1, sample))
    svd_res <- irlba::irlba(permdat, nv = 2)
    pvelist[i] <- svd_res$d[2]^2 / total_var
  }
  
  svd_res <- eigen(crossprod(imputed_dat), symmetric = TRUE, only.values = TRUE)
  K <- sum(svd_res$values / total_var > mean(pvelist))
  
  message('found K=',K)
  return(K)
}

EBMFimp<-function(x,rank=50,lambda=5,Kmax = 'auto',maxiter=300,return.flash=FALSE){
  require('softImpute')
  require('irlba')
  require('flashier')
  
  if(!'matrix'%in%class(x)){
    x<-as.matrix(x)
  }
  
  ###uses "svd" algorithm
  fit1 <- softImpute(as(x, "Incomplete"),rank = rank,lambda = lambda,type = "svd")
  pheno_soft <- complete(x,fit1)
  

  
  if(Kmax=='auto'){
    K<-deterK(pheno_soft)
    Kmax=K+7
    message('Kmax set to K+7=',Kmax)
    
  }
  min_sd <- min(apply(pheno_soft, 1, sd))
 
  pca_res <- irlba::irlba(pheno_soft, nv = Kmax)
  
  pca_res <- list(d=pca_res$d, u=pca_res$u, v=pca_res$v)
  
  fl_pca <- flash_init(x, S =  min_sd, var_type = 1) |>
  flash_factors_init(pca_res, ebnm_fn = ebnm::ebnm_point_laplace) |>
  flash_backfit(maxiter = maxiter)
  pheno.imp <- ifelse(is.na(x), fitted(fl_pca), x)
  
  
  if(return.flash){
    return(list(fl=fl_pca,Ximp=pheno.imp))
    
  }else{
      return(pheno.imp)
      
    }
}
  
gEBMF<-function(x,
                groups,
                nCores,
                probit=TRUE,
                save_flash=NULL,
                Kmax = 'auto',
                backfit_imput_iter=4,null_check=FALSE){
  #devtools::install_github('willwerscheid/flashier', ref = 'parallel-groups')
  
  require(flashier)
  require(snow)
  message('running optimized version of flashier...')
  if(probit)
    x<-qnorm(x)
  
  if(Kmax=='auto'){
    K<-deterK(x)
    Kmax=K+7
    message('Kmax set to K+7=',Kmax)
    
  }
  
  
  cl <- flash_init_cluster_for_grouped_data(x, groups, nCores , K = Kmax)
  
  for (i in seq_len(backfit_imput_iter)) {
    flash_backfit_grouped_data(cl, maxiter = 10)
    flash_impute_grouped_data(cl)
    
  }
  fl <- flash_recover_fl_object_from_cluster(cl)
  snow::stopCluster(cl)
  rm(cl)
  
  if(null_check){
    fl <- flash_nullcheck(fl)
    
  }
  
  if(!is.null(save_flash)){
    saveRDS(fl,save_flash)
  }
  
  
  x_imp <- ifelse(is.na(x), fitted(fl), x)
  if(probit)
    x_imp<-pnorm(x_imp)
  
  message('done.')
  
  return(x_imp)
}
