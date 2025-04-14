#ALL PCA BASED ANALYSIS UTILITIES

require('ggplot2')

RunPca<-function(norm_mat,scale=TRUE,imputation=NULL,center=TRUE){
  #return  the PCA of the transposed matrix : if the samples are the columns, factors will be estimated by samples
  n_features_init<-nrow(norm_mat)
  if(is.null(imputation)){
    norm_mat<-norm_mat[rowSums(is.na(norm_mat))==0,] #remove feature with NA
    norm_mat<-norm_mat[which(apply(norm_mat, 1, var)!=0),] #remove feature without variance
    n_features_qc<-nrow(norm_mat)
    message('removing ',n_features_init-n_features_qc,' features (',n_features_qc, 'remaining) because missing value or without variance')
    
  }else if(stringr::str_to_lower(imputation)%in%c('soft','softimpute')){
    norm_mat<-softImpute::complete(norm_mat, softImpute::softImpute(norm_mat))
    
  }
  
  
    return(prcomp(t(norm_mat),scale.=scale,center=center))
    
}

pctPC<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)
  names(pct.varPCs)<-paste0('PC',rngPCs)
  return( pct.varPCs)
}


PcaPlot<-function(pca,mtd,group.by, pc_x='PC1', pc_y='PC2',
                  sample_col='sample_id',scale=TRUE,return_pcs_mtd=FALSE,label=FALSE,ncol=NULL){
  require('ggrepel')
  require('patchwork')
  
  
  pca_dt<-data.table(pca$x,keep.rownames = sample_col)[mtd,on=sample_col]
  pctpcs<-pctPC(pca)
  
  pca_dt_toplot<-unique(pca_dt,by=sample_col)
  pca_dt_toplot[,sample_id:=.SD,.SDcols=sample_col]
  if(all(is.character(label))){
    pca_dt_toplot[,sample_to_label:=sample_id%in%label]
    label<-TRUE
  }else{
    pca_dt_toplot[,sample_to_label:=label]
    label<-TRUE
    
  }
  ps<-lapply(group.by, function(c){
    p<-ggplot(pca_dt_toplot,aes_string(x=pc_x,y=pc_y))+geom_point(aes_string(col=c))+
      labs(x=paste0(pc_x,' (',round(pctpcs[pc_x]*100),'%)'),
           y=paste0(pc_y,' (',round(pctpcs[pc_y]*100),'%)'))+
      theme_bw()
    return(p+ggtitle(c))
  })
  
  if(label){
    ps<-lapply(ps, function(p){
      p<-p+geom_text_repel(aes(label=ifelse(sample_to_label,sample_id,'')),
                                      max.overlaps =1000)
      return(p)
    })
  }
  
  
  if(return_pcs_mtd){
    print(wrap_plots(ps,ncol = ncol))
    return(pca_dt)
  }
  else{
    return(wrap_plots(ps,ncol = ncol))
  }
}



CorrelCovarPCs<-function(pca,mtd,
                         sample_col='sample_id',vars_num=NULL,
                         vars_cat=NULL,rngPCs=1:10){
  require(data.table)
  
  if(!sample_col%in%colnames(mtd)){
    stop(sample_col,' is not present in the metadata, please specify correct sample id column name')
  }
  
  if(!'data.table'%in%class(mtd)){
    mtd<-data.table(mtd)
  }
  pcs<-data.frame(pca$x)
  names(rngPCs)<-paste0('PC',rngPCs)
  
  if(is.null(vars_num))
    vars_num=colnames(mtd)[sapply(mtd, is.numeric)]
  if(is.null(vars_cat))
    vars_cat=colnames(mtd)[!sapply(mtd, is.numeric)&colnames(mtd)!=sample_col]

  mtd[,(vars_cat):=lapply(.SD,as.factor),.SDcols=vars_cat]
    
  if(length(vars_num)>0){
    
    res_num<-rbindlist(lapply(vars_num,function(f){
      
      res<-rbindlist(lapply(rngPCs,function(i){
        mod<-lm(pcs[,i]~as.numeric(unlist(mtd[rownames(pcs),..f,on=sample_col])))
        
        summstats<-summary(mod)
        data.table(PC=paste0('PC',i),
                   p=summstats$coefficients[2,4],
                   beta=summstats$coefficients[2,1],
                   R2=summstats$adj.r.squared)
        
      }))
      
      return(res[,factor:=f])
      
    }))
    
   
  }else{
    
    res_num<-data.table()
  }
  
  if(length(vars_cat)>0){
    
    res_fac<-rbindlist(lapply(vars_cat,function(f){
      
      res<-rbindlist(lapply(rngPCs,function(i){
        mod<-lm(pcs[,i]~as.factor(unlist(mtd[rownames(pcs),..f,on=sample_col])))
        
        summstats<-summary(mod)
        data.table(PC=paste0('PC',i),
                   p=anova(mod)$Pr[1],
                   R2=summstats$adj.r.squared)
        
      }))
      
      return(res[,factor:=f])
      
    }))
    
  }else{
    res_fac<-data.table()
    
    
  }
  res_all<-rbind(res_num,res_fac,fill=TRUE)
  res_all[,padj:=p.adjust(p),by='PC']
  return(res_all)
  
}


plotPvalsHeatMap<-function(x,main='-log10(Pvalue)',
                           p_col='p',p.thr=0.1,col_breaks=c(20,10:1, 0.5,0.1),
                           legend_breaks=NA,cluster_rows = F,cluster_cols = F,
                           labels_PC=NULL,
                           labels_Cov=NULL,
                           fontsize=10,
                           fontsize_number = 0.8*fontsize){
  require(pheatmap)
  
  if(p_col%in%colnames(x)){
    pvals_mat<-as.matrix(data.frame(dcast(x,factor~PC,value.var = p_col),row.names = 'factor'))
    pvals_mat<-pvals_mat[,paste0('PC',1:ncol(pvals_mat))]
  }else{
    pvals_mat<-x
  }
  pvals_mat[which(pvals_mat>p.thr)]<-1 #### put them to 1 if less than 0.1
  pvals_mat<--log10(pvals_mat)
  
  vars<-rownames(pvals_mat)

  pheatmap(pvals_mat,main = main,
           cluster_rows = cluster_rows,cluster_cols = cluster_cols,
           display_numbers = T,
           fontsize = fontsize,
           fontsize_number = fontsize_number,
           labels_col=labels_PC,
           labels_row=labels_Cov,
           color = colorRampPalette(c("white", "red"))(length(col_breaks)-1),
           breaks = sort(col_breaks),legend_breaks = legend_breaks)
  
}


correl<-function(x,y=NULL,ret="pval",verbose=F){
  if(is.null(y)){
    y=unlist(x[,2],use.names = F)
    x=unlist(x[,1],use.names = F)
  }
  
  if(is.numeric(x)){
    
    if(verbose){
      message("linear modeling ")
    }
    
    res<-lm(x~y)
    if(verbose){
      print(summary(res))
    }
    if(ret=="r2"){
      
      return(summary(res)$adj.r.squared)
    }else if(ret=="pval"){
      return(anova(res)$Pr[1])
      
    }else{
      return(
        list(p=anova(res)$Pr[1],r2=summary(res)$adj.r.squared)
      )
    }
    
  }else if(all(sapply(list(x,y),is.factor))){
    if(verbose){
      message("Chi-squared independance test")
    }
    
    x<-factor(x,levels = unique(x))
    y<-factor(y,levels = unique(y))
    tableF<-table(x,y)
    test<-chisq.test(tableF)
    if(verbose){
      print(test)
    }
    
    return(test$p.value)
  }
  
}


#PCsignCorr: allign the sign of the PC coordinate with the expression matrix. ie. a sample with a high expression will have the higher PC coord
#x: pc coords with sample names
#mat= scaled normalized expression matrix used for PC. if not for correct the PC1, should be filtered to have only the feature contributing tho the PC to be corrected
PCsignCorr<-function(x,mat,pc_col='PC1',sample_col='sample_id'){
  if(colSums(mat)[which.max(x)]-colSums(mat)[which.min(x)]<0){
    message('top and bottom sample have been inverted, correcting')
    return(-x)
  }else{
    message('PC and expression matrix well aligned')
    return(x)
  }
}