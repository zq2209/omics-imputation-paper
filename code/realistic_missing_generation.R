
LearnMissingPatterns<-function(mat,sample_cluster=10,feat_cluster=100,nHiMissThr=10){
  #learn missing pattern: determine the probability of the features being missing in group of features-samples (p_kg, P(missing)) 
  #AND the probability of the sample features being missing knowing their value (likelihood/evidence, i.e P(x_gj|missing)/P(x_gj))
  #For each methylation value x_ij from feature i of group g and sample j of group k,
  #we will draw missing using bayesian posterior probability P(missing|x_gj) = P(x_gj|missing)*P(missing)/P(x_gj)
  #require(biclust)
  require(pheatmap)
  require(data.table)
  require(stringr)
  
  #NA PATTERN EXTRACTION
  #1) cluster the data ~missingness 
  #into G features group, K samples clusters
  mat<-as.matrix(mat)
  na_mat<-is.na(mat)
  
  message('Extracting clusters with missing patterns...')
  clusters<-pheatmap(matrix(as.numeric(na_mat),ncol = ncol(na_mat),dimnames = list(rownames(na_mat),colnames(na_mat))),
                     show_rownames = F,show_colnames = F,
                     cluster_rows =   T,
                     clustering_distance_cols = 'euclidean',
                     kmeans_k = feat_cluster,
                     cutree_cols = sample_cluster
  )
  
  #get the samples clustering 
  clusters$tree_col$merge
  clusters_samples<-cutree(clusters$tree_col,k = sample_cluster)
  clusters_samples_dt<-data.table(cluster=paste0('k',clusters_samples),
                                  sample_id_corr = names(clusters_samples))
  
  #get the Features clustering 
  clusters_feats<-clusters$kmeans$cluster
  
  clusters_feats_dt<-data.table(cluster=paste0('g',clusters_feats),
                                feat_id= names(clusters_feats))
  
  #filter for cluster without enough samples/features
  #at least 5 features/samples
  clusters_samples_dt<-clusters_samples_dt[,n.cluster:=.N,by='cluster'][n.cluster>=5]
  clusters_feats_dt<-clusters_feats_dt[,n.cluster:=.N,by='cluster'][n.cluster>=5]
  
  #2)get the NA pattern prior P(missing)
  #extract the proportion of missingness per Samples-Features Clusters, 
  #it is the probability of the feature being missing in the sample-Features cluster 
  message('calculate missing rate in each clusters ...')
  
  missing_clusters<-sapply(unique(clusters_samples_dt$cluster), function(s_c){
    samples<-clusters_samples_dt[cluster==s_c]$sample_id_corr
    
    miss_feats<-sapply(unique(clusters_feats_dt$cluster), function(f_c){
      feats<-clusters_feats_dt[cluster==f_c]$feat_id
      miss<-sum(na_mat[feats,samples])/length(na_mat[feats,samples])
      return(miss)
    })
    
    return(miss_feats)
  })
  missing_clusters<-data.table(missing_clusters,keep.rownames = 'feat_cluster')
  missing_clusters<-melt(missing_clusters,id.vars = 'feat_cluster',variable.name = 'sample_cluster',value.name = 'prop_missing')
  message('missing rates: ')
  #3) get  the density function of the data distribution
  message('get the density functions for each clusters ...')
  plot(density(na.omit(mat)))
  dens_f_all <- approxfun(density(na.omit(mat)))
  
  
  #4) get the density function of highly missing Feature in each cluster of Samples and in each cluster of Features
  #We define highly missing Features being Feature with %NA above the median for features having missing value 
  missingness_FeatMissing<-rowMeans(na_mat[rowMeans(na_mat)>0,])
  med_FeatMissing<-median(missingness_FeatMissing)
  message('missingrate among missed featurew, feature median =',round(med_FeatMissing,digits = 3))
  
  #for Sample Clusters
  
  plot(density(na.omit(mat)))
  dens_f_HiMiss_cluster_samples<-lapply(unique(clusters_samples_dt$cluster), function(s_c){
    samples<-clusters_samples_dt[cluster==s_c]$sample_id_corr
    missingness_cluster<-rowSums(na_mat[,samples])/length(samples)
    
    hiMiss<-names(missingness_cluster)[which(missingness_cluster>med_FeatMissing)]
    
    hiMissMeth<-mat[hiMiss,samples]
    hiMissMeth<-hiMissMeth[!is.na(hiMissMeth)]
    
    
    if(length(hiMissMeth)>nHiMissThr){
      lines(density(hiMissMeth),col=str_extract(s_c,'[0-9]+'))
      
      dens_f_hiMiss<-approxfun(density(hiMissMeth))
      return(dens_f_hiMiss)
      
    }else{
      return(NA)
      
    }
    
  })
  names(dens_f_HiMiss_cluster_samples)<-unique(clusters_samples_dt$cluster)
  message(sum(!is.na(dens_f_HiMiss_cluster_samples)),'samples clusters with enough hiMissing observation to estimate a density distribution')
  
  
  #for features clusters
  missingness_SampleMissing<-colMeans(na_mat[rowMeans(na_mat)>0,])
  med_SampleMissing<-median(missingness_SampleMissing)
  message('missingrate among missed featurew, sample median =',round(med_SampleMissing,digits = 3))
  
  plot(density(na.omit(mat)))
  
  dens_f_HiMiss_cluster_feats<-lapply(unique(clusters_feats_dt$cluster), function(f_c){
    feats<-clusters_feats_dt[cluster==f_c]$feat_id
    missingness_cluster<-colSums(na_mat[feats,])/length(feats)
    
    hiMiss<-names(missingness_cluster)[which(missingness_cluster>med_SampleMissing)]
    
    hiMissMeth<-mat[feats,hiMiss]
    hiMissMeth<-hiMissMeth[!is.na(hiMissMeth)]
    
    if(length(hiMissMeth)>nHiMissThr){
      lines(density(hiMissMeth),col=as.numeric(str_extract(f_c,'[0-9]+')))
      dens_f_hiMiss<-approxfun(density(hiMissMeth))
      return(dens_f_hiMiss)
      
    }else{
      return(NA)
      
    }
    
  })
  names(dens_f_HiMiss_cluster_feats)<-unique(clusters_feats_dt$cluster)
  message(sum(!is.na(dens_f_HiMiss_cluster_feats)),'feature clusters with enough hiMissing observation to estimate a density distribution')
  #Finally get proportions of each clusters 
  propsK<-table(clusters_samples_dt$cluster)/nrow(clusters_samples_dt)
  propsG<-table(clusters_feats_dt$cluster)/nrow(clusters_feats_dt)
  
  #show NA patterns for QCs
  #pattern NA heatmap
  feats<-rownames(na_mat)
  sub_samples<-sample(colnames(na_mat),ifelse(ncol(na_mat)>100,100,ncol(na_mat)))
  sub_feats<-sample(feats,ifelse(length(feats)>10000,10000,length(feats)))
  
  ShowNAPattern(na_mat[sub_feats,sub_samples])
  message('Finished ! ')
  
  
  
  
  return(list(priors=missing_clusters,
              props=list(K=propsK,
                         G=propsG),
              clusters=list(K=clusters_samples_dt,
                            G=clusters_feats_dt),
              data_distrib=dens_f_all,
              missingData_distrib=list(K=dens_f_HiMiss_cluster_samples,
                                       G=dens_f_HiMiss_cluster_feats)))
}


ApplyMissingPattern<-function(mat,patterns,n_rounds=NULL,missrate.threshold=NULL,n_threads=1){
  require(data.table)
  require(pbapply)
  require(parallel)
  require(pheatmap)
  message('Applying miising patterns in the matrix... ')
  
 
  p_x<-function(x,dens_function,domain){
    # function to calculate likelihood of observed data according to a certain density function
    #if within definition domain calculate normally, otherwise reduce to domain
    x<-x[!is.na(x)]
    if(suppressWarnings(is.na(dens_function)))return(NA)
    
    if(length(x)<3)return(NA)
    x_med<-median(x)
    x_mad<-mad(x)
    lower<-x_med-x_mad
    upper<-x_med+x_mad
    if((lower<domain[1]&upper<domain[1])|(lower>domain[2]&upper>domain[2])){
      
      return(NA)
    }else{
      return(integrate(dens_function,
                       lower = max(c(lower,domain[1])),
                       upper =min(c(upper,domain[2])),
                       rel.tol =.Machine$double.eps^.1,
                       subdivisions = 500,
                       stop.on.error	=F)$value)
    }
    
    
  }
  
  
  Posterior<-function(likelihood,prior,evidence){
    # function to calculate posterior probability 
 
   
    return(likelihood*prior/evidence)
  }
  
  
  prodNA <- function(x, pNA ) {
    n <- length(x)
    NAloc <- unique(sample(1:n, floor(n * pNA),replace=T))
    x[NAloc] <- NA
    return(x)
  }
  
  #remove NA density function 
  #and the respective cluster
  k_with_dens<-names(patterns$missingData_distrib$K)[!is.na(patterns$missingData_distrib$K)]
  
  
  g_with_dens<-names(patterns$missingData_distrib$G)[!is.na(patterns$missingData_distrib$G)]
  
  #Get domain of definition of each density functions
  
  pdf(file = NULL)
  global_dom<-range(plot(patterns$data_distrib)$x[!is.na(plot(patterns$data_distrib)$y)])
  clustersK_dom<-lapply(patterns$missingData_distrib$K[k_with_dens],function(dens_func)
    range(plot(dens_func)$x[!is.na(plot(dens_func)$y)],na.rm = T))

  clustersG_dom<-lapply(patterns$missingData_distrib$G[g_with_dens],function(dens_func)
    range(plot(dens_func)$x[!is.na(plot(dens_func)$y)],na.rm = T))

  dev.off()

  # MISSING GENERATION#####
  mat_sim<-mat
  low_miss_rate<-TRUE
  i<-1
  #if is null nround and is null missrate = nround=1, missrate=1
  #if is null nround and  missrate = nround=Inf
  #if is  nround and is null missrate = missrate=1
  if(is.null(n_rounds)&is.null(missrate.threshold)){
    n_rounds<-1
    missrate.threshold<-1
  }else  if(is.null(n_rounds)){
    n_rounds<-Inf
    message('Missing Generation until reach missing rate ',round(missrate.threshold,digits = 3))
  }else if(is.null(missrate.threshold)){
    missrate.threshold<-1
    
  }
  

  while(i<=n_rounds&low_miss_rate){
    
    if(n_rounds>1)message('####### Round ',i,' ########')
    #2Sides missing generation: by samples cluster and then the transposed way
    i<-i+1
    
  #1) By sample Missing Gen####
  #split the matrix X randomly into K' groups of samples 
  #respecting the K sample clusters proportions
  message('Missing Generation by Sample')
  
  samples_clusters<-sample(names(patterns$props$K),size =ncol(mat_sim), replace = T,prob =patterns$props$K )
    
  samples_split<-split(colnames(mat_sim),samples_clusters)
  mat_split<-lapply(samples_split, function(samples)mat_sim[,samples])
  
  #then generate missing for each sample of a group and each feature
  #For each value x_ij from feature i and sample j of group k,
  #=> draw missing using  P(missing|x_gj) 
  #first get clusters of Features present in the matrix
  feats_clusters<-patterns$clusters$G[feat_id%in%rownames(mat_sim)]
  

  #then simulate
  mat_sim_list<-mclapply(mat_split[k_with_dens],mc.cores = n_threads,FUN= function(matf){
    samples<-colnames(matf)
    s_c<-names(samples_split)[sapply(samples_split,function(ss)all(samples%in%ss))]
    message('cluster ',s_c)
    
    
    #for each samples, and feature cluster, draw missing
    mat_sim_cluster<-pbsapply(samples,cl = n_threads, function(s){
      
      mat_sim_sample<-Reduce(`c`,lapply(unique(feats_clusters$cluster), FUN=function(f_c){
        #  message(f_c,' feature cluster')
        
        feats<-feats_clusters[cluster==f_c]$feat_id
        
        x_to_miss=matf[which(rownames(matf)%in%feats),s]
        
        pData<-p_x(x_to_miss,
                   dens_function = patterns$data_distrib,
                   domain = global_dom)
        
        pDataKnowingMissing<-p_x(x_to_miss,
                                 dens_function = patterns$missingData_distrib$K[[s_c]],
                                 domain = clustersK_dom[[s_c]])
        
        pMissing<-patterns$priors[feat_cluster==f_c&sample_cluster==s_c]$prop_missing
        
        pMissingKnowingData<-Posterior(likelihood = pDataKnowingMissing,
                                       prior = pMissing,
                                       evidence = pData)
        
        #  message(p,' probability of being missing')
        
        #draw missingness if probaMissing have been successfully calculated 
        if(length(pMissingKnowingData)>0){
          
          if(!is.na(pMissingKnowingData)){
           # if(pMissingKnowingData>1)pMissingKnowingData<-1
            x_to_miss<-prodNA(x_to_miss,pNA = pMissingKnowingData)
            
          }
        }
        
        return(x_to_miss)
      }))
      
      return(mat_sim_sample)
    })
    
    
    return(mat_sim_cluster)
    
  })
  
  #add clusters without missing
   mat_sim_list<-c(mat_sim_list,mat_split[setdiff(names(mat_split),k_with_dens)])
   #rm non matrix
   mat_sim_list<-mat_sim_list[sapply(mat_sim_list, function(x)'matrix'%in%class(x))]
   
     #get common features
  comm_features<-Reduce(intersect,lapply(mat_sim_list, function(x)rownames(x)))
  mat_sim_sample<-Reduce(cbind,lapply(mat_sim_list,function(mat)mat[match(comm_features,rownames(mat)),]))
  
  #add lacking rows
  row_to_add<-setdiff(rownames(mat),rownames(mat_sim_sample))
  if(length(row_to_add)>0){
    message(length(row_to_add),' rows have not been consider for sample simulation (features not present in the features clusters)')
    comm_samples<-intersect(colnames(mat_sim_sample),colnames(mat))
    mat_sim_sample<-rbind(mat_sim_sample[,comm_samples],mat[row_to_add,comm_samples,drop=F])
    
  }
  
  
  
  #2)BY feature missing gen####
  #split the features randomly into G' groups of features
  #respecting the G features clusters proportions
    message('Missing Generation by Feature')
  
  feats_clusters<-sample(names(patterns$props$G),size =nrow(mat_sim), replace = T,prob =patterns$props$G )
  

  feats_split<-split(rownames(mat_sim),feats_clusters)
  
  mat_split<-lapply(feats_split, function(feats)mat_sim[feats,])
  

  
  #then generate missing for each sample of a group and each feature
  #For each value x_ij from feature i and sample j of group k,
  #=> draw missing using  P(missing|x_gj)
  #first get clusters of Samples present in the matrice
  samples_clusters<-patterns$clusters$K[sample_id_corr%in%colnames(mat_sim)]
  
  mat_sim_list<-mclapply(mat_split[g_with_dens],mc.cores = n_threads,FUN= function(matf){
    feats<-rownames(matf)
    
    f_c<-names(feats_split)[sapply(feats_split,function(ff)all(feats%in%ff))]
    message('cluster ',f_c)
    
    #for each feature, and sample cluster, draw missing
    mat_sim_cluster<-pbsapply(feats, FUN=function(f){
      
      mat_sim_feat<-Reduce(`c`,lapply(unique(samples_clusters$cluster), FUN= function(s_c){
        #  message(f_c,' feature cluster')
        
        
        samples<-samples_clusters[cluster==s_c]$sample_id_corr
        # print(samples)
        x_to_miss=matf[rownames(matf)==f,samples]
        
        pData<-p_x(x_to_miss,
                   dens_function = patterns$data_distrib,
                   domain = global_dom)
        
        pDataKnowingMissing<-p_x(x_to_miss,
                                 dens_function = patterns$missingData_distrib$G[[f_c]],
                                 domain = clustersG_dom[[f_c]])
        
        pMissing<-patterns$priors[feat_cluster==f_c&sample_cluster==s_c]$prop_missing
        
        
        pMissingKnowingData<-Posterior(likelihood = pDataKnowingMissing,
                                       prior = pMissing,
                                       evidence = pData)
        
        
        
        #draw missingness if probaMissing have been successfully calculated 
        if(length(pMissingKnowingData)>0){
          if(!is.na(pMissingKnowingData)){
           # if(pMissingKnowingData>1)pMissingKnowingData<-1
            
            x_to_miss<-prodNA(x_to_miss,pNA =pMissingKnowingData )
            
          }
        }
        
        return(x_to_miss)
      }))
      
      return(mat_sim_feat)
    })
    
    return(t(mat_sim_cluster))
    
  })
  

  #add clusters without missing
  mat_sim_list<-c(mat_sim_list,mat_split[setdiff(names(mat_split),g_with_dens)])
  #rm non matrix
  mat_sim_list<-mat_sim_list[sapply(mat_sim_list, function(x)'matrix'%in%class(x))]
  
  #get common samples
  comm_samples<-Reduce(intersect,lapply(mat_sim_list, function(x)colnames(x)))

  mat_sim_feature<-Reduce(rbind,lapply(mat_sim_list,function(x)x[,comm_samples]))
  
  #add lacking columns
  col_to_add<-setdiff(colnames(mat),colnames(mat_sim_feature))
  
  if(length(col_to_add)>0){
    message(length(col_to_add),' columns have not been consider for sample simulation (samples not present in the sample clusters)')
    comm_feats<-intersect(rownames(mat_sim_feature),rownames(mat))
    mat_sim_feature<-cbind(mat_sim_feature[comm_feats,],mat[comm_feats,col_to_add,drop=F])
  }
  
  
  #merge####
  message('Merging by Sample and By Feature simulated matrices..')
  comm_feats<-intersect(rownames(mat_sim_sample),rownames(mat_sim_feature))
  comm_sample<-intersect(colnames(mat_sim_sample),colnames(mat_sim_feature))
  
  mat_sim<-mat_sim_feature[comm_feats,comm_sample]
  mat_sim[is.na(mat_sim_sample[comm_feats,comm_sample])]<-NA
  

  #add lacking row/columns
  row_to_add<-setdiff(rownames(mat),rownames(mat_sim))
  
  if(length(row_to_add)>0){
    message(length(row_to_add),' rows havent been consider for simulation')
    mat_sim<-rbind(mat_sim,mat[row_to_add,colnames(mat_sim),drop=F])
    
  }
  
  
  cols_to_add<-setdiff(colnames(mat),colnames(mat_sim))
  if(length(cols_to_add)>0){
    message(length(cols_to_add),' columns havent been consider for simulation')
    mat_sim<-cbind(mat_sim,mat[rownames(mat_sim),cols_to_add,drop=F])
    
    
  }

  #order matrix
  mat_sim<-mat_sim[rownames(mat),colnames(mat)]
  
  #show some stats
  NewNas<-is.na(mat)!=is.na(mat_sim)
  miss_rate<-sum(is.na(mat_sim))/length(mat_sim)
  
  new_miss_rate<-sum(NewNas)/length(NewNas)
  
  message(round(new_miss_rate*100,digits = 2),'% missing generated','(+',round(sum(NewNas)/sum(is.na(mat))*100),'%) ','\n total %NA =',round(miss_rate*100,digits = 2),'%' )
  
  low_miss_rate<-miss_rate<missrate.threshold
  
  }
  
  if(missrate.threshold<1){
    message('Missing rate threshold reached after ', i-1,' iteration')
  }
  
  #display matrix
  message('first 10 row and columns matrix:')
  print(head(mat_sim[,1:10]),10)
  
  return(mat_sim)
  
  
}

ApplyMissingPattern2<-function(query,reference,n_rounds=NULL,max.rounds=50,
                               missrate.threshold=NULL,min.new.missrate=0.0001,
                               n_threads=1,
                               masking=TRUE,
                               feat_by_cluster = 500,nfeat.cluster.max=100,nfeat.cluster.min=10,
                               samples_by_cluster =20,nsamples.cluster.max=20,nsamples.cluster.min=4,
                               min.n.missing=10,nHiMissThr=15){
  #Instead of using pMissingKnowingData to generate missingness randomly, 
  #the pMissingKnowingData is used to choose which sample/feature to draw in the reference to apply the missing 
  require(data.table)
  require(pbapply)
  require(parallel)
  require(pheatmap)
  
  #determine number of clusters to create based on wanted average number of samples / features per cluster
  n_samples<-ncol(reference)
  n_features<-nrow(reference)
  N_ref<-length(reference)
  
  n_sample_clusters<-round(n_samples/samples_by_cluster)
  if(n_sample_clusters>nsamples.cluster.max)n_sample_clusters=nsamples.cluster.max
  if(n_sample_clusters<nsamples.cluster.min)n_sample_clusters=nsamples.cluster.min
  
  message('will generate ',n_sample_clusters,' sample clusters')
  
  n_feat_clusters<-round(n_features/feat_by_cluster)
  if(n_feat_clusters>nfeat.cluster.max)n_feat_clusters=nfeat.cluster.max
  if(n_feat_clusters<nfeat.cluster.min)n_feat_clusters=nfeat.cluster.min
  message('will generate ',n_feat_clusters,' feature clusters')
  
  message('Learning missing patterns in the reference matrix... ')
  patterns<-LearnMissingPatterns(reference,sample_cluster =n_sample_clusters,
                                 feat_cluster = n_feat_clusters,nHiMissThr = nHiMissThr )
  
  message('Applying missing patterns in the query matrix... ')
  
  p_x<-function(x,dens_function,domain){
    # function to calculate likelihood of observed data according to a certain density function
    #if within definition domain calculate normally, if outside p = 0
    x<-x[!is.na(x)]
    if(suppressWarnings(is.na(dens_function)))return(NA)
    
    if(length(x)<3)return(NA)
    x_med<-median(x)
    x_mad<-mad(x)
    lower<-x_med-x_mad
    upper<-x_med+x_mad
    
    if((lower<domain[1]&upper<domain[1])|(lower>domain[2]&upper>domain[2])){
      
      return(0)
    }else{
      return(integrate(dens_function,
                       lower = max(c(lower,domain[1])),
                       upper =min(c(upper,domain[2])),
                       rel.tol =.Machine$double.eps^.1,
                       subdivisions = 500,
                       stop.on.error	=F)$value)
    }
    
  }
  
  
  Posterior<-function(likelihood,prior,evidence){
    # function to calculate posterior probability 
    p<-likelihood*prior/evidence
    return(p)
  }
  
  
  ApplyNARow <- function(x, pNA , ref_na) {
    if(pNA>1)pNA=1
    generate.miss<-sample(c(T,F),size = 1,prob = c(pNA,1-pNA))
    if(generate.miss){
      
      miss_rates<-rowMeans(ref_na)
      #for each pattern, the difference between the posterio probablity and the missing rate of this pattern is used as probablity to pick the pattern
      
      patternsP<-1-(pNA-miss_rates)
      patternsP<-ifelse(patternsP<0,0,patternsP)
      #use this Probabilities to pick the pattern
      i<-sample(1:nrow(ref_na),size=1,prob = patternsP)
      
      mask<-ref_na[i,][1:length(x)]
      x[mask] <- NA
        
      
      
    }
    
    return(x)
  }
  
  ApplyNACol <- function(x, pNA , ref_na) {
    if(pNA>1)pNA=1
    
    generate.miss<-sample(c(T,F),size = 1,prob = c(pNA,1-pNA))
    if(generate.miss){
      
      miss_rates<-colMeans(ref_na)
      #for each pattern, the difference between the posterio probablity and the missing rate of this pattern is used as probablity to pick the pattern
      
      patternsP<-1-(pNA-miss_rates)
      patternsP<-ifelse(patternsP<0,0,patternsP)
      #use this Probabilities to pick the pattern
      j<-sample(ncol(ref_na),size=1,prob = patternsP)
      
      mask<-ref_na[,j][1:length(x)]
      x[mask] <- NA
        
    }

    return(x)
  }
  
  #alternative version
  prodNA <- function(x, pNA ) {
    n <- length(x)
    NAloc <- unique(sample(1:n, floor(n * pNA),replace=T))
    x[NAloc] <- NA
    return(x)
  }
  
  #remove NA density function 
  #and the respective cluster
  k_with_dens<-names(patterns$missingData_distrib$K)[!is.na(patterns$missingData_distrib$K)]
  
  
  g_with_dens<-names(patterns$missingData_distrib$G)[!is.na(patterns$missingData_distrib$G)]
  

  #Get domain of definition of each density functions
  
  pdf(file = NULL)
  global_dom<-range(plot(patterns$data_distrib)$x[!is.na(plot(patterns$data_distrib)$y)])
  clustersK_dom<-lapply(patterns$missingData_distrib$K[k_with_dens],function(dens_func)
    range(plot(dens_func)$x[!is.na(plot(dens_func)$y)],na.rm = T))
  
  clustersG_dom<-lapply(patterns$missingData_distrib$G[g_with_dens],function(dens_func)
    range(plot(dens_func)$x[!is.na(plot(dens_func)$y)],na.rm = T))
  
  dev.off()
  
  # MISSING GENERATION#####
  mat_sim<-query
  ref_na<-is.na(reference)
  low_miss_rate<-TRUE
  i<-1
  rate_is_not_increasing<-c()
  old_miss_rate<-is.na(mat_sim)/length(mat_sim)
  
  #if is null nround and is null missrate = nround=1, missrate=1
  #if is null nround and  missrate = nround=Inf
  #if is  nround and is null missrate = missrate=1
  if(is.null(n_rounds)&is.null(missrate.threshold)){
    n_rounds<-1
    missrate.threshold<-1
  }else  if(is.null(n_rounds)){
    n_rounds<-Inf
    message('Missing Generation until reach missing rate ',round(missrate.threshold,digits = 3))
  }else if(is.null(missrate.threshold)){
    missrate.threshold<-1
    
  }
  
  
  while(i<=n_rounds&i<=max.rounds&low_miss_rate&sum(rate_is_not_increasing)<3){
    
    if(n_rounds>1)message('####### Round ',i,' ########')
    #2Sides missing generation: by samples cluster and then the transposed way
    i<-i+1
    
    #1) By sample Missing Gen####
    #split the matrix X randomly into K' groups of samples 
    #respecting the K sample clusters proportions
    message('Missing Generation by Sample')
    
    samples_clusters<-sample(names(patterns$props$K),size =ncol(mat_sim), replace = T,prob =patterns$props$K )
    
    samples_split<-split(colnames(mat_sim),samples_clusters)
    mat_split<-lapply(samples_split, function(samples)query[,samples])
    
    #then generate missing for each sample of a group and each feature
    #For each value x_ij from feature i and sample j of group k,
    #=> draw missing using  P(missing|x_gj) 
    #first get clusters of Features present in the matrix
    feats_clusters<-patterns$clusters$G[feat_id%in%rownames(mat_sim)]
    
    
    #then simulate
    mat_sim_list<-mclapply(mat_split[k_with_dens],mc.cores = n_threads,FUN= function(matf){
      samples<-colnames(matf)
      s_c<-names(samples_split)[sapply(samples_split,function(ss)all(samples%in%ss))]
      message('cluster ',s_c)
      
      
      #for each samples, and feature cluster, draw missing
      mat_sim_cluster<-pbsapply(samples,cl = n_threads, function(s){
        
        mat_sim_sample<-Reduce(`c`,lapply(unique(feats_clusters$cluster), FUN=function(f_c){
          #  message(f_c,' feature cluster')
          feats<-feats_clusters[cluster==f_c]$feat_id
          
          x_to_miss=matf[which(rownames(matf)%in%feats),s]
          
          pMissing<-patterns$priors[feat_cluster==f_c&sample_cluster==s_c]$prop_missing
          #do pattern simulation only if n missing in cluster > x
          
          nMissing<-pMissing*(N_ref*patterns$props$G[f_c]*patterns$props$K[s_c])
          if(nMissing>min.n.missing){
            pData<-p_x(x_to_miss,
                       dens_function = patterns$data_distrib,
                       domain = global_dom)
            
            pDataKnowingMissing<-p_x(x_to_miss,
                                     dens_function = patterns$missingData_distrib$K[[s_c]],
                                     domain = clustersK_dom[[s_c]])
            
            
            pMissingKnowingData<-Posterior(likelihood = pDataKnowingMissing,
                                           prior = pMissing,
                                           evidence = pData)
            
            #  message(p,' probability of being missing')
            
            #draw missingness if probaMissing have been successfully calculated 
            if(length(pMissingKnowingData)>0){
              
              if(!is.na(pMissingKnowingData)){
                # if(pMissingKnowingData>1)pMissingKnowingData<-1
                if(masking){
                  x_to_miss<-ApplyNACol(x_to_miss,pMissingKnowingData,
                                        ref_na = ref_na[patterns$clusters$G[cluster==f_c]$feat_id,patterns$clusters$K[cluster==s_c]$sample_id_corr])
                  
                }else{
                  x_to_miss<-prodNA(x_to_miss,pNA = pMissingKnowingData)
                  
                }
                
              }
            }
          }
          
          return(x_to_miss)
        }))
        
        return(mat_sim_sample)
      })
      
      
      return(mat_sim_cluster)
      
    })
    
    #add clusters without missing
    mat_sim_list<-c(mat_sim_list,mat_split[setdiff(names(mat_split),k_with_dens)])
    #rm non matrix
    mat_sim_list<-mat_sim_list[sapply(mat_sim_list, function(x)'matrix'%in%class(x))]
    #get common features
    comm_features<-Reduce(intersect,lapply(mat_sim_list, function(x)rownames(x)))
    mat_sim_sample<-Reduce(cbind,lapply(mat_sim_list,function(mat)mat[match(comm_features,rownames(mat)),]))
    
    
    #add lacking rows
    row_to_add<-setdiff(rownames(query),rownames(mat_sim_sample))
    if(length(row_to_add)>0){
      message(length(row_to_add),' rows have not been consider for sample simulation (features not present in the features clusters)')
      comm_samples<-intersect(colnames(mat_sim_sample),colnames(query))
      mat_sim_sample<-rbind(mat_sim_sample[,comm_samples],query[row_to_add,comm_samples,drop=F])
      
    }
    
    
    
    #2)BY feature missing gen####
    #split the features randomly into G' groups of features
    #respecting the G features clusters proportions
    message('Missing Generation by Feature')
    
    feats_clusters<-sample(names(patterns$props$G),size =nrow(mat_sim), replace = T,prob =patterns$props$G )
    
    
    feats_split<-split(rownames(mat_sim),feats_clusters)
    
    mat_split<-lapply(feats_split, function(feats)query[feats,])
    
    
    
    #then generate missing for each sample of a group and each feature
    #For each value x_ij from feature i and sample j of group k,
    #=> draw missing using  P(missing|x_gj)
    #first get clusters of Samples present in the matrice
    samples_clusters<-patterns$clusters$K[sample_id_corr%in%colnames(mat_sim)]
    
    mat_sim_list<-mclapply(mat_split[g_with_dens],mc.cores = n_threads,FUN= function(matf){
      feats<-rownames(matf)
      
      f_c<-names(feats_split)[sapply(feats_split,function(ff)all(feats%in%ff))]
      message('cluster ',f_c)
      
      #for each feature, and sample cluster, draw missing
      mat_sim_cluster<-pbsapply(feats, FUN=function(f){
        
        mat_sim_feat<-Reduce(`c`,lapply(unique(samples_clusters$cluster), FUN= function(s_c){
          #  message(f_c,' feature cluster')
          
          samples<-samples_clusters[cluster==s_c]$sample_id_corr
          # print(samples)
          x_to_miss=matf[rownames(matf)==f,samples]
          
          pMissing<-patterns$priors[feat_cluster==f_c&sample_cluster==s_c]$prop_missing
          #do pattern simulation only if prior missing rate > x%
          nMissing<-pMissing*(N_ref*patterns$props$G[f_c]*patterns$props$K[s_c])
          
          if(nMissing>min.n.missing){
            
            
          pData<-p_x(x_to_miss,
                     dens_function = patterns$data_distrib,
                     domain = global_dom)
          
          pDataKnowingMissing<-p_x(x_to_miss,
                                   dens_function = patterns$missingData_distrib$G[[f_c]],
                                   domain = clustersG_dom[[f_c]])
          print(pDataKnowingMissing)

          pMissingKnowingData<-Posterior(likelihood = pDataKnowingMissing,
                                         prior = pMissing,
                                         evidence = pData)
          
          
          
          #draw missingness if probaMissing have been successfully calculated 
          if(length(pMissingKnowingData)>0){
            if(!is.na(pMissingKnowingData)){
              # if(pMissingKnowingData>1)pMissingKnowingData<-1
              if(masking){
              x_to_miss<-ApplyNARow(x_to_miss,pMissingKnowingData,
                                    ref_na = ref_na[patterns$clusters$G[cluster==f_c]$feat_id,patterns$clusters$K[cluster==s_c]$sample_id_corr])
              
              }else{
                x_to_miss<-prodNA(x_to_miss,pNA = pMissingKnowingData)
                
              }
            }
          }
          
          }
          
          
          return(x_to_miss)
        }))
        
        return(mat_sim_feat)
      })
      
      return(t(mat_sim_cluster))
      
    })
    
    #add clusters without missing
    mat_sim_list<-c(mat_sim_list,mat_split[setdiff(names(mat_split),g_with_dens)])
    #rm non matrix
    mat_sim_list<-mat_sim_list[sapply(mat_sim_list, function(x)'matrix'%in%class(x))]
    
    #get common samples
    comm_samples<-Reduce(intersect,lapply(mat_sim_list, function(x)colnames(x)))
    mat_sim_feature<-Reduce(rbind,lapply(mat_sim_list,function(x)x[,comm_samples]))
    
    #add lacking columns
    col_to_add<-setdiff(colnames(query),colnames(mat_sim_feature))
    
    if(length(col_to_add)>0){
      message(length(col_to_add),' columns have not been consider for sample simulation (samples not present in the sample clusters)')
      comm_feats<-intersect(rownames(mat_sim_feature),rownames(query))
      mat_sim_feature<-cbind(mat_sim_feature[comm_feats,],query[comm_feats,col_to_add,drop=F])
    }
    
    
    #merge####
    message('Merging by Sample and By Feature simulated matrices..')
    comm_feats<-intersect(rownames(mat_sim_sample),rownames(mat_sim_feature))
    comm_sample<-intersect(colnames(mat_sim_sample),colnames(mat_sim_feature))
    
    mat_sim_new<-mat_sim_feature[comm_feats,comm_sample]
    mat_sim_new[is.na(mat_sim_sample[comm_feats,comm_sample])]<-NA
    
    
    
    #add lacking row/columns
    row_to_add<-setdiff(rownames(query),rownames(mat_sim_new))
    
    if(length(row_to_add)>0){
      message(length(row_to_add),' rows havent been consider for simulation')
      mat_sim_new<-rbind(mat_sim_new,query[row_to_add,colnames(mat_sim_new),drop=F])
      
    }
    
    
    cols_to_add<-setdiff(colnames(query),colnames(mat_sim_new))
    if(length(cols_to_add)>0){
      message(length(cols_to_add),' columns havent been consider for simulation')
      mat_sim_new<-cbind(mat_sim_new,query[rownames(mat_sim_new),cols_to_add,drop=F])
      
      
    }
    
    #order matrix
    mat_sim_new<-mat_sim_new[rownames(query),colnames(query)]
    
    #merge with NA generated in previous round
    mat_sim[is.na(mat_sim_new)]<-NA
 
    #save and show some missing rates stats
    NewNas<-is.na(query)!=is.na(mat_sim)
    miss_rate<-sum(is.na(mat_sim))/length(mat_sim)
    
    new_miss_rate<-sum(NewNas)/length(NewNas)
    
    message(round(new_miss_rate*100,digits = 2),'% missing generated','(+',round(sum(NewNas)/sum(is.na(query))*100),'%) ','\n total %NA =',round(miss_rate*100,digits = 2),'%' )
    
    low_miss_rate<-miss_rate<missrate.threshold
    
    rate_increase<-miss_rate-old_miss_rate
    
    rate_is_not_increasing<-c(rate_is_not_increasing,rate_increase<min.new.missrate)
    old_miss_rate<-miss_rate
  }
  
  if(missrate.threshold<1){
    message('Missing rate threshold reached after ', i-1,' iteration')
  }
  
  #display matrix
  message('first 10 row and columns matrix:')
  print(head(mat_sim[,1:10]),10)
  
  return(mat_sim)
  
  
}



ShowNAPattern<-function(mat,main='NA Patterns',cluster_rows=TRUE,cluster_cols=TRUE){
  require(pheatmap)

  if(!is.logical(mat)){
    mat<-is.na(mat)
    
  }
  #pattern NA heatmap
  matf<-mat[rowSums(mat)>0,]
  samples<-colnames(matf)
  feats<-rownames(matf)
  
  
  if(ncol(matf)>100){
    samples<-sample(samples,100)
    
  }
  
  if(nrow(matf)>10000){
    feats<-sample(feats,10000)
    
  }
  
  
  breaks <- c(0,0.5,1)
  colors <- c("#4575B4",'#D73027')

  
  message('NAs  pattern:')
  clusters<-pheatmap(matrix(as.numeric(matf[feats,samples]),ncol = length(samples),
                  dimnames = list(feats,samples)),
           clustering_distance_rows = 'euclidean',
           clustering_distance_col = 'euclidean',
           cluster_rows=cluster_rows,
           cluster_cols=cluster_cols,
           main=main,
           breaks = breaks,
           color = colors,
           show_rownames = F,show_colnames = F, legend = FALSE)
  
  invisible(list(clusters=clusters,features_shown=feats,samples_shown=samples))
  
}


