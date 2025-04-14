#R utils packages and functions for every project

#constantly used packages####
library("data.table")
setDTthreads(threads = 0)
message(getDTthreads(),' threads available for data.table')

library("stringr")
library("ggplot2")
library("ggrepel")
library("patchwork")
#library("here")
#library('qs')

#basic utils functions####

fp<-function(...)file.path(...)

ps<-function(...,sep="",collapse = NULL)paste(...,sep=sep,collapse = collapse)

source<-function(file,chdir=TRUE)base::source(file,chdir = chdir)

qs<-function()system('qstat -u adpelle1')


#detach_package
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}
#usage: detach_package(vegan)
#or detach_package("vegan", TRUE)

# Stats Related Functions ####
#zscore to pval
getPval<-function(zscore){
  return(2 * (1 - pnorm(abs(zscore))))
}


#CategoricalsToDummy
#inputs: covariates in data.table format
#outputs: dummy matrix : all categorical factors are transformed to binary (0/1) outcomes

CategoricalsToDummy<-function(covs){
  covs2<-copy(covs)
  id_col<-colnames(covs2)[as.vector(unlist(covs2[,lapply(.SD, function(x)length(unique(x))==.N)]))]
  num_vars<-colnames(covs2[,.SD,.SDcols=is.numeric])
  cat_vars<-setdiff(colnames(covs2)[unlist(covs2[,lapply(.SD, function(x)all(table(x)>1))])],num_vars)
  others_vars<-setdiff(colnames(covs2),c(num_vars,cat_vars,id_col))
  message(paste(others_vars,collapse = ', '),' are categorical variables with only 1 sample falling in one category, so wont be Dummyfied')
  dum_vars<-c()
  
  for(cat in cat_vars){
    t<-table(as.factor(unlist(covs2[,..cat])))
    lvls<-names(t)[-1]
    cols<-paste0(cat,lvls)
    covs2[,(cols):=lapply(lvls, function(l)as.numeric(.SD==l)),.SDcols=cat]
    
    dum_vars<-c(dum_vars,cols) 
  }
  
  return(covs2[,.SD,.SDcols=c(id_col,num_vars,others_vars,dum_vars)])
}



GetVarPCs<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)
  names(pct.varPCs)<-rngPCs
  return( pct.varPCs)
}


#Reformatting classical R data/results ####
UniqueClean<-function(x,key_cols='sample_id',pattern_to_exclude=NULL){
  #keep only unconstant variables
  nums_to_keep<-names(which(x[,sapply(.SD,function(y)(!all(is.na(y)))&var(y,na.rm = T)!=0)&length(unique(y))!=.N,.SDcols=is.numeric]))
  
  cats_to_keep<-names(which(x[,sapply(.SD,function(y)(!all(is.na(y)))&length(unique(y))!=.N&length(unique(y))!=1),.SDcols=!is.numeric]))
  
  if(!is.null(pattern_to_exclude)){
    nums_to_keep<-nums_to_keep[!str_detect(nums_to_keep,pattern_to_exclude)]
    cats_to_keep<-cats_to_keep[!str_detect(cats_to_keep,pattern_to_exclude)]
    
  }
  cols_to_keep<-c(key_cols,nums_to_keep,cats_to_keep)
  
  return(unique(x,by=c(key_cols))[,.SD,.SDcols=cols_to_keep])
}

RemoveUselessColumns<-function(x,key_cols='sample_id',pattern_to_exclude=NULL){
  #keep only unconstant variables
  nums_to_keep<-names(which(x[,sapply(.SD,function(x)(!all(is.na(x)))&var(x,na.rm = T)!=0),.SDcols=is.numeric]))
  
  cats_to_keep<-names(which(x[,sapply(.SD,function(x)(!all(is.na(x)))&length(unique(x))!=.N&length(unique(x))!=1),.SDcols=!is.numeric]))
  
  if(!is.null(pattern_to_exclude)){
    nums_to_keep<-nums_to_keep[!str_detect(nums_to_keep,pattern_to_exclude)]
    cats_to_keep<-cats_to_keep[!str_detect(cats_to_keep,pattern_to_exclude)]
    
  }
  cols_to_keep<-c(key_cols,nums_to_keep,cats_to_keep)
  
  return(x[,.SD,.SDcols=cols_to_keep])
}

DetectDEResFormat<-function(res_de){
  
  if(all(c('log2FoldChange','padj')%in%colnames(res_de)))return('DESEQ2')
  else if(any(c('avg_log2FC','avg_logFC')%in%colnames(res_de))&'p_val_adj'%in%colnames(res_de))return('SEURAT')
  else return('unknown')
  
}

FormatDEResToSeurat<-function(res_de){
  res_forms<-list(DESEQ2=c(padj='padj',
                           pval='pvalue',
                           stat='stat',
                           FC='log2FoldChange'),
                  SEURAT=c(padj='p_val_adj',
                           pval='p_val',
                           stat='stat',
                           FC='avg_log2FC'))
  
  if(DetectDEResFormat(res_de)=='DESEQ2'){
    setnames(res_de,old = res_forms[['DESEQ2']],
             new = res_forms[['SEURAT']])
    return(res_de)
  }else{
    stop('unknown format')
  }
  return(res_de)
}

ReFormatDERes<-function(res_de,to='SEURAT'){
  if(to=='SEURAT') ifelse(DetectDEResFormat(res_de)==to,return(res_de),return(FormatDEResToSeurat(res_de)))
  else{
    stop('unsupported formatting')
  }
}


#GGPLOT####

bar_bw<-function()scale_fill_manual(values=c('black','grey'))
bar_rb<-function(invert=FALSE)ifelse(invert,return(scale_fill_manual(values=c('royalblue3','orangered3'))),return(scale_fill_manual(values=c('orangered3','royalblue3'))))

#BIOMART####


GetMartGenes<-function(mirror=NULL)biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror=mirror)
GetMartMouseGenes<-function(mirror=NULL)biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl",mirror=mirror)

GetMartReg<-function()biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")
GetMartMotif<-function()biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_motif_feature")



GetBiomartAttrs<-function(mart)data.table::data.table(biomaRt::listAttributes(mart))
GetBiomartFilter<-function()data.table::data.table(biomaRt::listFilter(mart))

GetBiomartAttrs_reg<-function()data.table::data.table(listAttributes(biomaRt::useEnsembl(biomart = "regulation",dataset = "hsapiens_regulatory_feature")))
GetBiomartFilter_reg<-function()data.table::data.table(biomaRt::listFilter(biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")))

TransNMtoSymbol<-function(refseq_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(getBM(attributes = c('refseq_mrna', 'hgnc_symbol',add.attributes),
                                      filters = 'refseq_mrna', 
                                      values = refseq_ids, 
                                      mart = GetMartGenes())))
}

TransTranscriptToGene<-function(transcript_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id',add.attributes),
                                               filters = 'ensembl_transcript_id', 
                                               values = transcript_ids, 
                                               mart = GetMartGenes())))
}

TransSymboltoNM<-function(hgnc_symbols,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('refseq_mrna', 'hgnc_symbol',add.attributes),
                                               filters = 'hgnc_symbol', 
                                               values = hgnc_symbols, 
                                               mart = GetMartGenes())))
}

TransEnsembltoSymbol<-function(ensembl_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol',add.attributes),
                                               filters = 'ensembl_gene_id', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

TransSymboltoEnsembl<-function(hgnc_symbols,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol',add.attributes),
                                               filters = 'hgnc_symbol', 
                                               values = hgnc_symbols, 
                                               mart = GetMartGenes())))
}

TransEnsemblVerstoSymbol<-function(ensembl_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol',add.attributes),
                                               filters = 'ensembl_gene_id_version', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

TransSymboltoEnsemblVers<-function(hgnc_symbols,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol',add.attributes),
                                               filters = 'hgnc_symbol', 
                                               values = hgnc_symbols, 
                                               mart = GetMartGenes())))
}



TransEnsembltoNM<-function(ensembl_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'refseq_mrna',add.attributes),
                                               filters = 'ensembl_gene_id', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

TransNMtoEnsembl<-function(refseq_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'refseq_mrna',add.attributes),
                                               filters = 'refseq_mrna', 
                                               values = refseq_ids, 
                                               mart = GetMartGenes())))
}

TransEnsemblVerstoNM<-function(ensembl_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'refseq_mrna',add.attributes),
                                               filters = 'ensembl_gene_id_version', 
                                               values = ensembl_ids, 
                                               mart = GetMartGenes())))
}

TransNMtoEnsemblVers<-function(refseq_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'refseq_mrna',add.attributes),
                                               filters = 'refseq_mrna', 
                                               values = refseq_ids, 
                                               mart = GetMartGenes())))
}


TransTranscriptToGene<-function(transcript_ids,add.attributes=NULL){
  require(biomaRt)
  require(data.table)
  return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id',
                                                              'ensembl_transcript_id',add.attributes),
                                               filters = 'ensembl_transcript_id', 
                                               values = transcript_ids, 
                                               mart = GetMartGenes())))
}



tr<-function(ids_sepBySlash,retourne="all",sep="/",tradEntrezInSymbol=FALSE,uniqu=TRUE){
  IDs<-as.vector(strsplit(ids_sepBySlash,sep)[[1]])
  if(retourne=="all"){
    ret<-1:length(IDs)
  }else{
    ret<-retourne
  }
  if(tradEntrezInSymbol){
    require(clusterProfiler)
    library(org.Hs.eg.db)
    if(retourne=="all"){
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL)
    }
    else{
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL[ret])
    }
    
  }else{
    if(uniqu){
      return(unique(IDs[ret]))
      return(IDs[ret])
    }
    
    
  }
  
}







# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x,return_dt=T,mirror=NULL){
  
  require("biomaRt")
  human =GetMartGenes(mirror=mirror)
  mouse = GetMartMouseGenes(mirror=mirror)
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  if(return_dt)return(data.table(genesV2))
  
  mousex <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

convertMouseGeneList <- function(x,return_dt=T,mirror=NULL){
  
  require("biomaRt")
  human =GetMartGenes(mirror=mirror)
  mouse = GetMartMouseGenes(mirror=mirror)
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  if(return_dt)return(data.table(genesV2))
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}



#Over repre / GSEA####

OR<-function(set1,set2,size_universe){
  if(any(duplicated(set1))){
    set1<-unique(set1)
  }
  if(any(duplicated(set2))){
    set2<-unique(set2)
  }
  phyper(q=sum(set1%in%set2)-1, 
         #number of white balls drawn without replacement from an urn which contains both black and white balls. 
         #<=> number of tcf4 target in DEGs. "-1" because normally give P[X > x] but we want P[X >= x])
         m=length(set2), #the number of white balls in the urn. <=> number of DEGs 
         n=size_universe-length(set2), #the number of black balls in the urn. <=> number of genes tested - number of DEGs
         k=length(set1), #the number of balls drawn from the urn <=> numbers of tcf4 targets 
         lower.tail=FALSE)  # if TRUE (default), probabilities are P[X ≤ x] (under-representation), otherwise, P[X > x] (over-representation).
  
}



OR2<-function(querys,terms_list,size_universe,min.term.size=0,max.term.size=Inf,
              overlap_column=TRUE,verbose=TRUE){
  if(is.list(querys)){
    return(Reduce(rbind,lapply(names(querys),
                               function(q)OR2(querys = querys[[q]],
                                              terms_list = terms_list,
                                              size_universe = size_universe,
                                              min.term.size = min.term.size,
                                              max.term.size = max.term.size,
                                              overlap_column = overlap_column,
                                              verbose = verbose )[,query:=q])))
  }else{
    res_or<-data.table(term=names(terms_list),term.size=sapply(terms_list,length))
    res_or<-res_or[term.size<=max.term.size]
    n_terms<-nrow(res_or)
    if(verbose)message(length(terms_list)-n_terms, " terms were filtered due to term.size above the limit of ",max.term.size," genes")
    res_or<-res_or[term.size>=min.term.size]
    if(verbose)message(n_terms-nrow(res_or), " terms were filtered due to term.size below the limit of ",min.term.size," genes")
    
    res_or[,n.query:=length(querys)]
    res_or[,n.overlap:=sum(querys%in%terms_list[[term]]),by="term"]
    if(overlap_column==TRUE){
      res_or[,genes.overlap:=paste(querys[querys%in%terms_list[[term]]],collapse="|"),by="term"]
    }
    res_or[,pct.query.overlap:=n.overlap/n.query]
    res_or[,precision:=pct.query.overlap]
    
    res_or[,pct.term.overlap:=n.overlap/term.size]
    
    res_or[,background_size:=size_universe]
    
    res_or[,pct.background:=term.size/size_universe] #TO IMPROVE (Here size universe can be the intersection between 2 "univers" while n.query is the n.query is the query universe. 
    
    
    res_or[,pval:=phyper(q=n.overlap-1, 
                         m=n.query, 
                         n=size_universe-n.query, 
                         k=term.size, 
                         lower.tail=FALSE),
           by="term"]
    res_or[,padj:=p.adjust(pval,method = 'BH')]
    if(verbose)message(nrow(res_or[padj<0.05])," terms enriched in your genes of interest with padj<0.05")
    return(res_or)
  }
  
}

OR3<-function(querys,terms_list,background,
              min.term.size=0,
              max.term.size=Inf,
              min.query.size=5,
              overlap_column=TRUE,verbose=FALSE,underrepr=FALSE){
  if(is.list(querys)){
    querys_sizes=sapply(querys,length)
    dt<-Reduce(rbind,lapply(names(querys)[querys_sizes>=min.query.size],
                            function(q)OR3(querys = querys[[q]],
                                           terms_list = terms_list,
                                           background = background,
                                           min.term.size = min.term.size,
                                           max.term.size = max.term.size,
                                           overlap_column = overlap_column,
                                           verbose = verbose )[,query:=q]))
    
    return(dt[,query.:=query][,.SD,.SDcols=c(ncol(dt),1:(ncol(dt)-1))])
  }else{
    queryf<-intersect(querys,background)
    terms_listf<-lapply(terms_list,function(x)intersect(x,background))
    res_or<-data.table(term=names(terms_listf),term.size=sapply(terms_listf,length))
    res_or<-res_or[term.size<=max.term.size]
    n_terms<-nrow(res_or)
    if(verbose)message(length(terms_listf)-n_terms, " terms were filtered due to term.size above the limit of ",max.term.size," genes")
    res_or<-res_or[term.size>=min.term.size]
    n_terms<-nrow(res_or)
    if(verbose)message(n_terms-nrow(res_or), " terms were filtered due to term.size below the limit of ",min.term.size," genes")
    
    res_or[,n.query:=length(queryf)]
    res_or[,n.overlap:=sum(queryf%in%terms_listf[[term]]),by="term"]
    
    res_or[,pct.query.overlap:=n.overlap/n.query]
    res_or[,precision:=pct.query.overlap]
    
    res_or[,pct.term.overlap:=n.overlap/term.size]
    
    res_or[,background_size:=length(background)]
    
    res_or[,pct.term.background:=term.size/background_size] 
    
    
    res_or[,pval:=phyper(q=ifelse(underrepr,n.overlap,n.overlap-1), 
                         m=term.size, 
                         n=background_size-term.size, 
                         k=n.query, 
                         lower.tail=underrepr),
           by="term"]
    res_or[,padj:=p.adjust(pval,method = 'BH')]
    res_or[,fold.enrichment:=pct.query.overlap/pct.term.background]
    if(overlap_column==TRUE){
      res_or[,genes.overlap:=paste(queryf[queryf%in%terms_listf[[term]]],collapse="|"),by="term"]
    }
    if(verbose)message(nrow(res_or[padj<0.05])," terms enriched in your genes of interest with padj<0.05")
    return(res_or[order(padj)])
  }
  
}

#Liftover / Coordinates changes
LiftOver<-function(bed,chainfile,
                   col.names=c('chr.new','start.new','end.new','id'),
                   select=NULL){
  require(data.table)
#pip install CrossMap
  #if 'ImportError: libhts.so.3: cannot open shared object file: No such file or directory'
  #htslib folder should be declare in ~/.Renviron : i.e. 'LD_LIBRARY_PATH=/path/to/htslib-x.xx.x:$LD_LIBRARY_PATH'
  if(is.null(select)){
    select<-1:length(col.names)
  }else{
    col.names<-col.names[select]
  }
  
  in_file<-"temp_in.bed"
  out_file<-"temp_out.bed"
  
  fwrite(bed,in_file,col.names = F,sep="\t")
  
  system(paste("micromamba run -n utils CrossMap bed",chainfile,in_file,out_file))
  
  bed_trans<-fread(out_file,select = select,col.names = col.names)
  file.remove(c(in_file,out_file))
  
  return(bed_trans)
}

hg19to38<-function(bed){
  chainfile="/projectnb/tcwlab/RefData/liftover/hg19ToHg38.over.chain.gz"
  bed_trans<-LiftOver(bed,chainfile =chainfile ,col.names=c('chr','start.hg38','end.hg38','id'),
)
  return(bed_trans)
}

hg38to19<-function(bed){
  chainfile="/projectnb/tcwlab/RefData/liftover/hg38ToHg19.over.chain.gz"
  bed_trans<-LiftOver(bed,chainfile =chainfile,col.names=c('chr','start.hg19','end.hg19','id'))
  return(bed_trans)
}


FindGOGenes<-function(terms_or_ids){
  require("biomaRt")
  require("stringr")
  
  if(!all(str_detect(terms_or_ids,"^GO:")))terms_or_ids=FindGO_ID(term_description=terms_or_ids)
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
  #gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
  gene.data <- getBM(attributes=c('hgnc_symbol','go_id'),
                     filters = 'go', values = terms_or_ids, mart = ensembl,uniqueRows = T)
  return(data.table(gene.data)[go_id%in%terms_or_ids])
}



FindGO_ID<-function(term_descriptions){
  require("GO.db")
  terms<-Term(GOTERM)
  ids<-names(terms)[match(term_descriptions,unlist(terms))]
  return(ids)
}




#genomics coordinates manipulation ####
start<-function(x,start_pos=2)sapply(x,function(x)as.numeric(strsplit(x,"\\.|-|:|_|,|\\[|\\]")[[1]][start_pos]))
end<-function(x,end_pos=3)sapply(x,function(x)as.numeric(strsplit(x,"\\.|-|:|_|,|\\[|\\]")[[1]][end_pos]))

seqid<-function(x,only_num=FALSE){
  if(only_num){
    str_extract(x,'[0-9]+')|>as.numeric()
  }else{
    sapply(x,function(x)strsplit(x,"\\.|-|:|_|,|\\[|\\]")[[1]][1])
    
  }
}

pos<-function(x)sapply(x,function(x)as.numeric(strsplit(x,"\\.|-|:|_|,|\\[|\\]")[[1]][2]))


ref<-function(x)sapply(x,function(x){
  vec=strsplit(x,"-|:|_|,|\\[|\\]")[[1]]
  return(vec[length(vec)-1])
})

alt<-function(x)sapply(x,function(x){
  vec=strsplit(x,"-|:|_|,|\\[|\\]")[[1]]
  return(vec[length(vec)])
})


#FlipGeno
#flip e.g. 0/1 in 1/0 in a vector. works also with the separator '|'
#warning: do not support non-biallelic allele
FlipGeno<-function(x){
  require(stringr)
  sapply(x, function(y){
    sep=str_extract(y,'[/|]')
    info_vec=strsplit(y,':')[[1]]
    geno_vec<-info_vec[1]|>strsplit(paste0("\\",sep))|>unlist()
    geno_vec<-ifelse(geno_vec=='0','1',ifelse(geno_vec=='1','0',geno_vec))
    
    if(length(info_vec)>1){
      geno_info=paste(paste(geno_vec,collapse  = sep),paste(info_vec[-1],collapse = ':'),sep=':')
    }else{
      geno_info=paste(geno_vec,collapse  = sep)
    }
    
    return(geno_info)
  })
}


#AlleleFlipping
#flip alleles of vcf-like table (modifying `ID` (optional), `REF`, `ALT` and `$genotype` columns) according to a reference of variants with the format #CHROM:POS:REF:ALT
#arguments:
#vcf: the vcf-like table in data.table::data.table() format 
#ref_snps: the SNPs to be flipped accordingly
#genotype_columns=names of the genotypes columns.
#edit_ID: if need to edit the SNPs IDs, default TRUE
#all.x= if need to return the full vcf, even SNPs not found in the reference. default FALSE

#Note: make sure chromosome ID have been formated the same way (e.g. both with the suffixe 'chr')
#Value: Return a data.table of the vcf like table reformatted according to the SNPs reference
AlleleFlipping<-function(vcf,ref_snps,genotype_columns='genotype',edit_ID=TRUE,
                         all.x=FALSE,verbose=NULL,nThreads=NULL){
  require(data.table)
  require(parallel)
  
  if(!'data.table'%in%class(vcf)){
    message('converting to data.table')
    vcf<-data.table(vcf)
  }
  vcf[,`#CHROM`:=as.character(`#CHROM`)]
  
  
  if(is.null(nThreads))nThreads=1
  ref_snps<-unique(as.vector(unlist(ref_snps)))
  #filter for SNPs pos present in the vcf
  
  ref_snps_list<-strsplit(ref_snps,"-|:|_|,|\\[|\\]")
  
  ref_snps<-data.table(`#CHROM`=sapply(ref_snps_list,function(x)x[1]),
                       POS=sapply(ref_snps_list,function(x)as.numeric(x[2])),
                       ID=ref_snps)
  
  ref_snpsf<-merge(ref_snps,unique(vcf[,.(`#CHROM`,POS)]),by=c('#CHROM','POS'))
  message(nrow(ref_snpsf),' ref SNPs found in the vcf (',round(nrow(ref_snpsf)/nrow(ref_snps)*100,digits = 1),'% of the refs, ',
          round(nrow(ref_snpsf)/nrow(vcf)*100,digits = 1),'% of the vcf)')
  
  ref_snpsf[,REF:=ref(ID)]
  ref_snpsf[,ALT:=alt(ID)]
  
  ref_snps_matching<-merge(ref_snpsf,
                           unique(vcf[,.(`#CHROM`,POS,REF,ALT)]),by=c('#CHROM','POS','REF','ALT'))
  
  vcf2<-vcf[ref_snps_matching[,.(`#CHROM`,POS,REF,ALT)],on=c('#CHROM','POS','REF','ALT')]
 
  vcf2[,flipped_allele:=F]
  
  if(edit_ID){
    #rename ID for matching variants
    message('renaming IDs for matching variants based on reference')
    cols_order=colnames(vcf2)
    vcf2<-merge(vcf2[,-'ID'],ref_snps_matching[,.(`#CHROM`,POS,REF,ALT,ID)],by=c('#CHROM','POS','REF','ALT'))
    setcolorder(vcf2,cols_order)
  }
  
  ref_snpsf_toflip<-ref_snpsf[!ID%in%ref_snps_matching$ID]
  
  message(nrow(ref_snps_matching),' SNPs matching alleles.')
  
  if(nrow(ref_snpsf_toflip)>0){
    message('try flipping for ',nrow(ref_snpsf_toflip),' variants')
    
    if(nrow(ref_snpsf_toflip)>1000&is.null(verbose))verbose=FALSE
    else verbose=TRUE
    
    vcf_flipped<-rbindlist(mclapply(ref_snpsf_toflip$ID,function(snp){
      
      chr=ref_snpsf[ID==snp]$`#CHROM`
      p=ref_snpsf[ID==snp]$POS
      r=ref(snp)
      a=alt(snp)
      
      vcff=vcf[`#CHROM`%in%chr&POS==p]
      
      if (all(vcff[,REF==a&ALT==r])){
        if(verbose)message('flipping ', snp)
        vcff[,flipped_allele:=T]
        vcff[,REF:=r]
        vcff[,ALT:=a]
        
        if(edit_ID&all(vcff[,ID!=snp])){
          vcff[,ID:=snp]
        }
        
        vcff[,(genotype_columns):=lapply(.SD,function(x)as.vector(FlipGeno(x))),.SDcols=genotype_columns]
        
      }else{
        vcff[,flipped_allele:=NA]
        message(snp, ' uncorrect alt/ref, returning NA (query = ',unique(vcff$REF),'/',unique(vcff$ALT),')')    
        
      }
      
      return(vcff)
      
    },mc.cores=nThreads))
    
    message(nrow(vcf_flipped),'/',  nrow(vcf),' alleles flipped')
    
    vcf2<-rbind(vcf_flipped,vcf2)
  }
  
  if(all.x){
    tested<-vcf[unique(vcf2[,.(`#CHROM`,POS)]),on=c('#CHROM','POS')]$ID
    vcf2<-rbind(vcf2,vcf[!ID%in%tested],fill=T)
    #reorder cols
    first_cols=c('#CHROM','POS',"ID","REF","ALT")
    setcolorder(vcf2,c(first_cols,setdiff(colnames(vcf2),first_cols)))
  }

  return(vcf2[order(`#CHROM`,POS)])
}




MethChangeReg<-function(res_meth,region){
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  res_meth_reg<-res_meth[chr==chromosome&pos>start.pos&pos<end.pos]
  res_meth_reg[,start:=pos][,end:=pos+1]
  return(res_meth_reg)
}

MethChangePlot<-function(res_meth,region,limits=NULL,breaks=waiver()){
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  res_meth_reg<-res_meth[chr==chromosome&pos>start.pos&pos<end.pos]
  res_meth_reg[,start:=pos][,end:=pos+1]
  p<-ggplot(data = res_meth_reg) + geom_segment(aes(x = start, y = 0, 
                                                    xend = end, yend = logFC,col=-log10(P.Value)), size = 2, data = res_meth_reg)+
    scale_color_gradient(low = "white",high = "black",limits=limits,breaks=breaks)
  
  p<-p+ theme_classic() + ylab(label = "Methylation change") + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))
  
  return(p)
}


#DATA.TABLE ADD IN
freadvcf<-function(file){
  cmd<-paste('zcat',file,'| grep -v "^##"')
  fread(cmd)
}




#QSUB FILES CREATION####


CreateJobFile<-function(cmd_list,file,proj_name='tcwlab',modules=NULL,
                        loadBashrc=FALSE,conda_env=NULL,
                        micromamba_env=NULL,
                        cwd='.',
                        nThreads=NULL,memPerCore=NULL,maxHours=24,
                        parallelize=FALSE,maxChildJobs=60){
  
  template_header='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_header.txt'
  template_tail='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_tail.txt'
  if(!str_detect(file,'\\.qsub$')){
    file=paste0(file,'.qsub')
  }
  filename<-basename(file)
  projdir<-dirname(file)
  
  
  while(str_detect(projdir,'scripts')){
    projdir<-dirname(projdir)
  }
  
  dir.create(file.path(projdir,'logs'),showWarnings = F)
  file_path<-file
  log_file=file.path(projdir,'logs',paste0(str_remove(filename,'.qsub$'),'.log'))
  
  #create the qsub file
  cat('#!/bin/bash -l\n',file = file_path)
  message('qsub file created at ',file_path)
  
  #set the general job parameters
  proj_name_opt<-paste('-P ',proj_name)
  CombStdOutErr_opt<-'-j y'
  maxHours_opt<-paste0('-l h_rt=',maxHours,':00:00')
  qlog<-paste('-o ',log_file)
  if(!is.null(nThreads)){
    nThreads_opt=paste('-pe omp',nThreads)
  }else{
    nThreads_opt=NULL
    
  }
  
  if(!is.null(memPerCore)){
    memPerCore_opt=paste0('-l mem_per_core=',memPerCore)
  }else{
    memPerCore_opt=NULL
    
    }
  #create child qsub files if parallelize mode
  if(parallelize&length(cmd_list)>1){
    
    script_id=str_extract(filename,'^[0-9A-Za-z]+')
    scripts_dir<-file.path(projdir,'scripts')
    childs_dir<-file.path(scripts_dir,paste0(script_id,'childs'))
    
   if(!dir.exists(childs_dir))dir.create(childs_dir,recursive = T)

    #add the main job parameters
    cat( '#Parameters of the Jobs :',file = file_path,append = T)
    cat( c('\n',proj_name_opt,CombStdOutErr_opt,maxHours_opt,qlog),
         file = file_path,append = T,sep = '\n#$')
    cat( '\n',file = file_path,append = T,sep = '\n')
    
    #add the header
    system(paste('cat',template_header,'>>',file_path))
    
    
      if(!is.null(names(cmd_list))){
        child_jobnames<-make.names(names(cmd_list))
      }else{
        child_jobnames<-paste0(script_id,'-child',1:length(cmd_list))
      }
    
    child_jobfiles<-file.path(childs_dir,paste0(child_jobnames,'.qsub'))
   
     #add command to run the childs qsub
    n_jobs<-ifelse(length(cmd_list)<=maxChildJobs,length(cmd_list),maxChildJobs)
    
   
    cmd_list_jobs<-split(cmd_list,1:n_jobs)
    
    
    cmds<-sapply(1:n_jobs,function(i){
      
      cmds_child<-cmd_list_jobs[[i]]
      CreateJobFile(cmd_list =cmds_child,file = child_jobfiles[i],
                    cwd = cwd,
                    nThreads = nThreads ,proj_name = proj_name ,
                    loadBashrc = loadBashrc,modules = modules,conda_env = conda_env,micromamba_env = micromamba_env,
                    maxHours =  maxHours,memPerCore = memPerCore,
                    parallelize = FALSE)
      
      # Submit the job using qsub and capture the job ID
      cmd_job<-paste(paste0('job_id',i,'=$(qsub'), '-N',paste0('j',child_jobnames[[i]]),child_jobfiles[[i]],proj_name,' | grep -Ewo [0-9]+)')
      # cmd<-RunQsub(child_jobfiles[[i]],job_name = paste0('j',script_id,child_jobnames[[i]]),dryrun = T)
      return(cmd_job)
      
    })
    
    cat( cmds,file = file_path,append = T,sep = '\n')
    
    # Wait until the jobs are completed
    for(i in 1:n_jobs){
      cmd_wait<-c(paste('while qstat | grep -w',paste0('"',paste0('$job_id',i),'"'),  '> /dev/null; do'),
                  paste('echo', '"waiting jobs to complete.."'),
                  paste('sleep', '30'), # Adjust sleep time as needed
                  'done',
                  paste('echo',paste0('"Job',i),'have completed."'))
      cat( c('\n',cmd_wait),file = file_path,append = T,sep = '\n')
      
    }
    
    #add the tail
    system(paste('cat',template_tail,'>>',file_path))
    
   
    
  }else{
    #add the core job parameters
  cat( '#Parameters of the Jobs :',file = file_path,append = T)
  cat( c('\n',proj_name_opt,CombStdOutErr_opt,maxHours_opt,qlog,nThreads_opt,memPerCore_opt),file = file_path,append = T,sep = '\n#$')
  cat( '\n',file = file_path,append = T,sep = '\n')
  
  #add the module to loads
  if(!is.null(modules)){
    modules<-ifelse(modules=='R','R/4.2.1',modules)
    
    if('gatk'%in%modules){
      
      modules<-c(modules[modules!='gatk'],
                 'miniconda/23.1.0',
            'java/17.0.8',
            'gatk/4.4.0.0')
      
      conda_env<-union(conda_env,'/share/pkg.8/gatk/4.4.0.0/install/gatk-4.4.0.0')
      
    }
    
    cat( '#Modules to load:',file = file_path,append = T)
    cat( c('\n',modules),file = file_path,append = T,sep = '\nmodule load ')
    cat( '\n',file = file_path,append = T,sep = '\n')
    
  }
  #load .bashrc
  if(loadBashrc){
    cat( '# loading of bashrc profile:',file = file_path,append = T)
    cat('\n source $HOME/.bashrc',file = file_path,append = T)
    cat( '\n',file = file_path,append = T,sep = '\n')
  }
  
  #activate conda environment 
  if(!is.null(conda_env)){
    cat( '#Conda environment activation:',file = file_path,append = T)
    cat( c('\n',conda_env),file = file_path,append = T,sep = '\nconda activate ')
    cat( '\n',file = file_path,append = T,sep = '\n')
    
  }
  
  #activate micromamba environment 
  if(!is.null(micromamba_env)){
    cat( '#Micromamba environment activation:',file = file_path,append = T)
    
    micromamba_env<-ifelse(micromamba_env=='fungen','pisces-rabbit',micromamba_env)
    
    cat( c('\n',micromamba_env),file = file_path,append = T,sep = '\nmicromamba activate ')
    
    if(any(c('pisces-rabbit','fungen')%in%micromamba_env)){
      cat( c('\n# singularity specifics configuration: ',
             'export SINGULARITY_BIND="$TMP,/restricted/projectnb/tcwlab-adsp/,/projectnb/tcwlab-adsp/,/projectnb/tcwlab/"',
             'touch .Rprofile'),
           file = file_path,append = T,sep = '\n')

    }
    cat( '\n',file = file_path,append = T,sep = '\n')
    
  }
  
  if('simpleaf'%in%micromamba_env){
    cat( c('\n# simpleaf specifics configuration: ',
           'export ALEVIN_FRY_HOME="$PWD"', #define env variable
           'simpleaf set-paths', #simpleaf configuration
           'ulimit -n 2048'), #update the maximum number of file descriptors a program can create
         file = file_path,append = T,sep = '\n')
    
    
  }
  
  
  #add the header
  system(paste('cat',template_header,'>>',file_path))
  
  #go to cwd
  if(cwd!='.'){
    cat( paste('mkdir -p',cwd),file = file_path,append = T,sep = '\n')
    cat( paste('cd',cwd),file = file_path,append = T,sep = '\n')
    
    cat( '\n',file = file_path,append = T,sep = '\n')
    
  }
  
  #add the commandes to exec

  if(!is.null(names(cmd_list))){
    cmds<-unlist(lapply(names(cmd_list),
                        function(s)return(c(paste('echo',paste0('"----- Processing of ',s,' -----"')),cmd_list[[s]]))))
    
  }else{
    cmds<-unlist(cmd_list)
  }
  
  cat( cmds,file = file_path,append = T,sep = '\n')
  
  
  
  
  if(any(c('pisces-rabbit','fungen')%in%micromamba_env)){
    #go back to cwd
    if(cwd!='.'){
      cwd<-getwd()
      cat( paste('cd',cwd),file = file_path,append = T,sep = '\n')
      
      cat( '\n',file = file_path,append = T,sep = '\n')
      
    }
    
    cat( c('\n# Singularity specifics clean up: ',
           'rm .Rprofile'),file = file_path,append = T,sep = '\n')
  }
  cat( '\n',file = file_path,append = T,sep = '\n')
  
  #add the tail
  system(paste('cat',template_tail,'>>',file_path))
  
 
  #show the file header
  message('header:')
  system(paste('head -n 15',file_path))
  }
  
  #show the first commands
  message('5 first commands:')
  cat( head(cmds,5),sep = '\n')
  
  return(invisible(file_path))
}



CreateJobForPyFile<-function(python_file,proj_name='tcwlab',modules=NULL,
                        loadBashrc=FALSE,conda_env=NULL,
                        micromamba_env=NULL,
                        nThreads=4,memPerCore=NULL,maxHours=24){
  template_header='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_header.txt'
  template_tail='/projectnb/tcwlab/LabMember/adpelle1/utils/template/qsub_file_tail.txt'
  
  filename<-basename(python_file)
  projdir<-ifelse(str_detect(python_file,'scripts/'),dirname(dirname(python_file)),dirname(python_file))
  
  
  #create qsub file
  qsub_file<-str_replace(python_file,'\\.py$','.qsub')
  cat('#!/bin/bash -l\n',file = qsub_file)
  
  #create log file
  dir.create(file.path(projdir,'logs'),showWarnings = F)
  log_file=file.path(projdir,'logs',str_replace(filename,'\\.py$','.log'))
  
  
  #add the job parameters
  proj_name<-paste('-P ',proj_name)
  CombStdOutErr<-'-j y'
  maxHours<-paste0('-l h_rt=',maxHours,':00:00')
  qlog<-paste('-o ',log_file)
  if(!is.null(nThreads))nThreads=paste('-pe omp',nThreads)
  if(!is.null(memPerCore))memPerCore=paste0('-l mem_per_core=',memPerCore)
  
  cat( '#Parameters of the Jobs :',file = qsub_file,append = T)
  cat( c('\n',proj_name,CombStdOutErr,maxHours,qlog,nThreads,memPerCore),file = qsub_file,append = T,sep = '\n#$')
  cat( '\n',file = qsub_file,append = T,sep = '\n')
  
  #add the module to loads
  if(!is.null(modules)){
    modules<-ifelse(modules=='R','R/4.2.1',modules)
    
    if('gatk'%in%modules){
      
      modules<-c(modules[modules!='gatk'],
                 'miniconda/23.1.0',
                 'java/17.0.8',
                 'gatk/4.4.0.0')
      
      conda_env<-union(conda_env,'/share/pkg.8/gatk/4.4.0.0/install/gatk-4.4.0.0')
      
    }
    
    cat( '#Modules to load:',file = qsub_file,append = T)
    cat( c('\n',modules),file = qsub_file,append = T,sep = '\nmodule load ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  #load .bashrc
  if(loadBashrc){
    cat( '# loading of bashrc profile:',file = qsub_file,append = T)
    cat('\n source $HOME/.bashrc',file = qsub_file,append = T)
    cat( '\n',file = qsub_file,append = T,sep = '\n')
  }
  
  #activate conda environment 
  if(!is.null(conda_env)){
    cat( '#Conda environment activation:',file = qsub_file,append = T)
    cat( c('\n',conda_env),file = qsub_file,append = T,sep = '\nconda activate ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  
  #activate micromamba environment 
  if(!is.null(micromamba_env)){
    cat( '#Micromamba environment activation:',file = qsub_file,append = T)
    cat( c('\n',micromamba_env),file = qsub_file,append = T,sep = '\nmicromamba activate ')
    cat( '\n',file = qsub_file,append = T,sep = '\n')
    
  }
  

  #add the header
  system(paste('cat',template_header,'>>',qsub_file))
  
  
  #add the Python script to exec
  cmd<-paste('python',python_file,'>>',log_file)
  
  cat( cmd,file = qsub_file,append = T,sep = '\n')
  cat( '\n',file = qsub_file,append = T,sep = '\n')
  
  #add the tail
  system(paste('cat',template_tail,'>>',qsub_file))
  
  #show the file header
  message('qsub file created at ',qsub_file)
  message('header:')
  system(paste('head -n 15',qsub_file))
  
  #show the  bash command
  message('the bash command to execute:')
  cat( cmd,sep = '\n')
  
  #show the first R lines
  message('the 15 firsts lines of Python scripts to execute:')
  system(paste('head -n 15',python_file))
  
  
  
}

#CreateJobForRfile~~~~
#create job to execute an R file using Rscript rfile.R [args]
#Arguments:
#args: should be a list of argument to pass to Rscript command. If the first argument is a vector, will create one Rscript command by element of this vector
# eg Rscript rfile.R element1
# Rscript rfile.R element2


CreateJobForRfile<-function(r_file,qsub_file=NULL,args=NULL,
                            parallelize=NULL,maxChildJobs=60,
                            proj_name='tcwlab',
                            modules='R',
                            loadBashrc=FALSE,conda_env=NULL,micromamba_env=NULL,
                            nThreads=4,memPerCore=NULL,maxHours=24){
  
  if(is.null(parallelize)){
    parallelize=FALSE
    if(!is.null(args)){
      if(length(args[[1]])>1){
        parallelize=TRUE
      }
    }
  }
  
  if(is.null(qsub_file)){
    qsub_file<-str_replace(r_file,'\\.R$','.qsub')
  }
  
  filename<-basename(r_file)
  projdir<-dirname(r_file)
  while(str_detect(projdir,'scripts')){
    projdir<-dirname(projdir)
  }
  
  #create log file
  dir.create(file.path(projdir,'logs'),showWarnings = F)
  log_file=file.path(projdir,'logs',str_replace(filename,'\\.R$','.log'))
  
  
  #create the Rscript cmds to exec
  if(is.null(args)){
    cmds<-paste('Rscript',r_file,'>>',log_file)
    
  }else{
  cmds<-lapply(as.character(args[[1]]), function(arg1){
    log_file=paste0(str_replace(log_file,'.log$','_'),make.names(basename(arg1)),'.log')
    cmd<-paste('Rscript',r_file,arg1,paste(unlist(args[-1]),
                                           collapse = ' '),'>>',log_file)
    return(cmd)
  })
  
  }
  #show the first R lines
  message('the 15 firsts R lines to executes')
  system(paste('head -n 15',r_file))
  
  #pass the argument to createJobFile
  CreateJobFile(cmds,file = qsub_file,proj_name = proj_name ,
                modules = modules ,micromamba_env = micromamba_env ,conda_env = conda_env ,
                nThreads=nThreads,memPerCore=memPerCore,maxHours=maxHours,
                loadBashrc = loadBashrc,
                parallelize = parallelize,
                maxChildJobs = maxChildJobs)
  
}

RunQsub<-function(qsub_file=NULL,job_name=NULL,proj_name=NULL,wait_for=NULL,dryrun=FALSE){
  if(is.null(qsub_file)){
    qsub_file=.Last.value
  }
  if(!str_detect(qsub_file,'.qsub$'))qsub_file=paste0(tools::file_path_sans_ext(qsub_file),'.qsub')
  
  if(!file.exists(qsub_file))stop(qsub_file,' does not exist')
  
  if(is.null(job_name)){
    job_name<-str_sub(str_remove(basename(qsub_file),'[0-9A-Za-z]+-'),end = 15)
  }
  message('running job file ', qsub_file,' with name ', job_name)

  if(is.null(wait_for)){
    cmd<-paste('qsub','-N',job_name,qsub_file,proj_name)

  }else{
    cmd<-paste('qsub',
               '-N',job_name,
               '-hold_jid',wait_for,
               qsub_file,proj_name)
  }
  if(!dryrun){
    message<-system(cmd,intern = TRUE)
    message(paste(message,collapse = '\n'))

    jobid<-ifelse(str_detect(message,'has been submitted'),str_extract(message,'[0-9]+'),NA)
    return(jobid[!is.na(jobid)])
  }else{
    return(cmd)
  }

}



WaitQsub<-function(qsub_file,jobid=NULL,max_hours=24){
  if(!str_detect(qsub_file,'\\.qsub$'))qsub_file=paste0(tools::file_path_sans_ext(qsub_file),'.qsub')
  
  if(!file.exists(qsub_file))stop(qsub_file,' file does not exist')
  
  filename<-basename(qsub_file)
  projdir<-dirname(qsub_file)
  t0 <- Sys.time()
  
  while(str_detect(projdir,'scripts')){
    projdir<-dirname(projdir)
  }
  log_file=file.path(projdir,'logs',paste0(str_remove(filename,'.qsub$'),'.log'))
  
  max_sec=max_hours*60*60
  
  
  while(!file.exists(log_file)){
    Sys.sleep(60)
    t1 <- Sys.time()
    d<-t1-t0
    if(d>max_sec)stop('time reached max limits')
    message('waiting job ',jobid,' to be executed..')
    
    
  }
  
  jobid_inlog=str_extract(tail(grep('job ID :',
                                    readLines(con = log_file,skipNul = TRUE),
                                    value = T),n = 1),'[0-9]+$')
  if(is.null(jobid)){
    message('jobid not provided, using last jobid found in log file. ')
    jobid=jobid_inlog
  }

  while(jobid!=jobid_inlog){
    message('waiting job ',jobid,' to be executed..')
    Sys.sleep(60)
    t1 <- Sys.time()
    d<-t1-t0
    if(d>max_sec)stop('time reached max limits')
    jobid_inlog=str_extract(tail(grep('job ID :',
                                      readLines(con = log_file,skipNul = TRUE),
                                      value = T),n = 1),'[0-9]+$')
    
    
  }
  
  pattern=paste('Finished Analysis for job',jobid)
  lastlines<-tail(readLines(con = log_file,skipNul = TRUE))

  while(!any(str_detect(lastlines,pattern)) ){
    message('waiting job ',jobid,' to be completed..')
    Sys.sleep(120)
    lastlines<-tail(readLines(con = log_file,skipNul = TRUE))
    t1 <- Sys.time()
    d<-t1-t0
    if(d>max_sec)stop('time reached max limits')
  }
  t1 <- Sys.time()
  d<-t1-t0
  message(pattern,' after ', round(d/360,2),' hours')
  
  
}


#bedtools/vcf tools wrapper####
bed_inter<- function(a, b, opt1="-wa",
                     opt2="-wb",out_dir=".",
                     select=NULL, col.names=NULL){
  require(data.table)
  l<-list(a,b)
  #order chr and pos
  l<-lapply(l, function(x){
    chr_cols<-colnames(x)[1:2]
    setorderv(x,chr_cols)
  })
  
  files_to_rm<-c(FALSE,FALSE)
  file_paths<-sapply(1:2, function(i){
    x<-l[[i]]
    if(is.data.frame(x)){
      files_to_rm[i]<-TRUE
      file_path<-fp(out_dir,paste0("temp",i,".bed"))
      fwrite(x,file_path,sep="\t",col.names = FALSE,scipen = 999)
      
    }else{
      file_path<-x
    }
    return(file_path)
  })
  
  out_file<-fp(out_dir,"temp_inter.bed")
  cmd<-paste("bedtools intersect -a",file_paths[1],"-b",file_paths[2], opt1, opt2,">",out_file)
  message("run in shell : ",cmd)
  system(cmd)
  message("done.")
  
  if(!is.null(col.names)){
    dt<-fread(out_file,select = select,col.names = col.names)
    file.remove(out_file)
    file.remove(file_paths[files_to_rm])
    return(dt)
  }else{
    dt<-fread(out_file,select = select)
    file.remove(out_file)
    file.remove(file_paths[files_to_rm])
    return(dt)
  }
}


#CountBEDOverlap
#inputs : 1) bed files of the reads and 2) bed files of the genomics regions to count overlap
#for each genomic region, count how many read overlap using bedtools intersect -c 

#outputs : bed file <bed_file_name>.<genomics_regions_name>.overlap.count.bed.gz of the number of read falling in each genomic region
CountBEDOverlap<-function(bed_files,genomic_regions_file,
                          out_dir=NULL,job_file='X-count_bed_overlap.qsub',
                          job_name = 'countbedoverlap',
                          write_genomic_regions=TRUE, write_bed_file_intervals=FALSE,
                          nThreads=NULL,parallelize=F,
                          maxChildJobs=40,wait_job=TRUE,
                          wait_for=NULL,dryrun=FALSE){
  bed_dir=unique(dirname(bed_files))
  if(is.null(out_dir)){
    out_dir=bed_dir
    if(length(out_dir)>1){
      stop('specify outputs directory')
    }
  }
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  
  opt1=ifelse(write_genomic_regions,"-wa",'')
  opt2=ifelse(write_bed_file_intervals,"-wb",'')
  
  
  cmds<-lapply(bed_files, function(b){
    paste('bedtools intersect -c',opt1,opt2,'-a',genomic_regions_file,
          '-b',b,'|gzip -c >',
          fp(out_dir,ps(tools::file_path_sans_ext(basename(b)),'.',tools::file_path_sans_ext(basename(genomic_regions_file)),
                        '.overlap.count.bed.gz')))
  })
  
  if(is.null(job_file)){
    
    #job_file<-fp('scripts',ps('counts_overlap_',tools::file_path_sans_ext(bed_files[1],'_andCo_with_',tools::file_path_sans_ext(basename(genomic_regions_file)),'.qsub')))
    
    for(cmd in cmds){
      message('running ',cmd)
      
      if(!dryrun){
        system(cmd)
      }      
      
    }
    
  }else{
    
    CreateJobFile(cmds,file = job_file,
                  loadBashrc = T,modules = c('bedtools'),
                  nThreads = nThreads,parallelize =parallelize,
                  maxChildJobs=maxChildJobs)
    
    if(!dryrun){
      jobid<-RunQsub(job_file,job_name = job_name,wait_for = wait_for)
      
      if(wait_job)WaitQsub(job_file,jobid =jobid )
      
    }
    
  }
  
  count_files<- fp(out_dir,ps(tools::file_path_sans_ext(basename(bed_files)),'.',tools::file_path_sans_ext(basename(genomic_regions_file)),
                              '.overlap.count.bed.gz'))
  
  return(count_files)
  
}

#Example:
# introns_counts<-CountBEDOverlap(bed_files =bed_files,
#                                 genomic_regions_file = 'examples_data/some_introns.bed.gz')
# 
# exons_counts<-CountBEDOverlap(bed_files =bed_files,
#                                 genomic_regions_file = 'examples_data/some_exons.bed.gz')
# 



RunGatk<-function(cmd){
  message('run in a terminal:')
  
  cat(paste('module load miniconda/23.1.0',
            'module load java/17.0.8',
            'module load gatk/4.4.0.0',
            'conda activate /share/pkg.8/gatk/4.4.0.0/install/gatk-4.4.0.0',
            cmd,sep = '\n'))
}

#plink wrapper####
#TransToVCF
#inputs bedfile to convert in vcf
#outputs : vcf file (compressed with bgzipped if bgzip=T), at the same location than the bed file
#note : return in R the path to the new vcf file. 
#if dry_run=TRUE, do not execute the code but raterh return the codes in addition to the path to the future vcf (in a list format where $cmds contain the commands, and $vcf contain the path to the vcf)
TransToVCF<-function(bed_file,bgzip=T,dry_run=FALSE,threads=8){
  file_pref<-tools::file_path_sans_ext(bed_file)
  vcf_file=paste0(file_pref,'.vcf')  
  vcfgz_file=paste0(vcf_file,'.gz')
  
  
  tovcf<-paste('plink',
               '--bfile',file_pref,
               paste('--threads',threads),
               '--recode vcf',
               '--out',file_pref)
  zip=paste('bgzip -f',vcf_file)
  tabix=paste('tabix',vcfgz_file)
  clean=paste('rm -f',vcf_file)
  
  if(!dry_run){
    system(tovcf)
    if(bgzip){
      system(zip)
      system(tabix)
      system(clean)
      return(vcfgz_file)
    }else{
      return(vcf_file)
      
    }
    
  }else{
    if(bgzip){
      return(list(cmds=c(tovcf,zip,tabix,clean),vcf=vcfgz_file))
      
    }else{
      return(list(cmds=c(tovcf),vcf=vcf_file))
      
    }
  }
  
}

#genomics annotations access####

GetGTF<-function(gtf='/projectnb/tcwlab/RawData/Genome2.7.9a/hg38/gencode.v26.annotation.gtf',features=NULL){
  if(!is.null(features)){
    if(str_detect(gtf,'.gz$')){
      cmd=paste('zcat',gtf,
                "|grep -wE",paste0("'",paste(features,collapse = '|'),"'"))
      
      
    }else{
      cmd=paste("grep -wE",paste0("'",paste(features,collapse = '|'),"'"),gtf)
      
    }

    return(fread(cmd = cmd,col.names = c('chr','source','feature_type','start','end','score','strand','frame','anno')))
    
  }else{
    return(fread(gtf,col.names = c('chr','source','feature_type','start','end','score','strand','frame','anno')))
    
  }
  
}
GetGenesGTF<-function(genes,genome='hg38',gtf=NULL,anno_sources=c('ENSEMBL','HAVANA'),feature_types='gene'){
  #feature_type : gene, transcript, exon, CDS, UTR, start_codon, stop_codon
  if(genome=='hg38'&is.null(gtf)){
    gtf='/projectnb/tcwlab/RefData/gencode/hg38/gencode.v45.annotation.gtf'
  }
  if(genome=='hg19'&is.null(gtf)){
    gtf='/projectnb/tcwlab/RefData/gencode/hg19/gencode.v19.annotation.gtf.gz'
  }
  
  gtf<-GetGTF(features=genes,gtf = gtf)
  gtff<-gtf[source%in%anno_sources&feature_type%in%feature_types]
  gtff[,gene_name:=sapply(anno,function(x){
    annos<-strsplit(x,'; ')[[1]]
    annos<-annos[str_detect(annos,'gene_name')]
    anno=str_remove_all(annos,'"| |gene_name')
    return(anno)
  })]
  gtff[,gene_id:=sapply(anno,function(x){
    annos<-strsplit(x,'; ')[[1]]
    annos<-annos[str_detect(annos,'gene_id')]
    anno=str_remove_all(annos,'"| |gene_id')
    return(anno)
  })]
  if('transcript'%in%feature_types){
    
    gtff[feature_type=='transcript',transcript_name:=sapply(anno,function(x){
      annos<-strsplit(x,'; ')[[1]]
      annos<-annos[str_detect(annos,'transcript_name')]
      anno=str_remove_all(annos,'"| |transcript_name')
      return(anno)
    })]
    gtff[feature_type=='transcript',transcript_id:=sapply(anno,function(x){
      annos<-strsplit(x,'; ')[[1]]
      annos<-annos[str_detect(annos,'transcript_id')]
      anno=str_remove_all(annos,'"| |transcript_id')
      return(anno)
      
    })]
    gtff[feature_type=='transcript',transcript_type:=sapply(anno,function(x){
      annos<-strsplit(x,'; ')[[1]]
      annos<-annos[str_detect(annos,'transcript_type')]
      anno=str_remove_all(annos,'"| |transcript_type')
      return(anno)
    })]
  }
  genes_not_found<-setdiff(genes,unique(gtff$gene_name))
  if(length(genes_not_found)>0){
    warning(paste(genes_not_found,collapse=','),' coordinates not found in ',anno_sources, 'with feature',features ,'. Try with other source/feature')
    warning('available sources: ',paste(unique(gtf$source),collapse = ','))
    warning('available gene features: ',paste(unique(gtf$feature_type),collapse = ','))
    
  }
  return(gtff)
}

#GetSNPsInfos
#extract SNPs info from SNP database
#if flip_alleles=TRUE, the allele in strand '-' are converted in alleles in strand '+'
GetSNPsInfos<-function(rsids,snp_db='/projectnb/tcwlab/RefData/SNPs/snp141.txt.gz',flip_alleles=TRUE){
  
  cmd=paste('zcat',snp_db,
            "|grep -wE",paste0("'",paste(rsids,collapse = '|'),"'"))
  #or : zcat snp137.txt.gz | grep -Fwf myRsIDs.txt > mySnps.txt
  # fwrite(data.table(rsids),tempfile,col.names =FALSE)
  # cmd=paste('zcat',snp_db,
  #           "|grep -Fwf",tempfile)
  snps_infos<-fread(cmd = cmd,select=c(2:5,7,8,9,10,12,14,15,16,18),
        col.names =c('chr','start','end','rsid','strand','REF','ref_ucsc','alleles','type','maf','maf_se','function','align_qual'))
  snps_infos<-snps_infos[chr%in%paste0('chr',c(1:22,'X','Y'))]
  
  if(flip_alleles){
    message('converting alleles strands if strand==-')
    
    snps_infos[strand=='-',alleles:=strandflip(alleles)]
    snps_infos[strand=='-',note:='strand flipped alleles']
    snps_infos[,strand:='+']
    
  }
  
  snps_infos[,ALT:=paste(setdiff(strsplit(alleles,'/')[[1]],REF),collapse='/'),by='rsid']
  
  return(snps_infos)
  #col.names = c('chr','source','feature_type','start','end','score','strand','frame','anno')
  
}

strandflip<-function(nt_or_allele){
 seps=str_extract_all(nt_or_allele,'\\/|\\||\\:|\\-|_')[[1]]
 sep=seps[length(seps)]
vecs=strsplit(nt_or_allele,sep)
flips=lapply(vecs,strflip)

return(sapply(flips, function(x)paste(x,collapse = sep)))

}

strflip<-function(x){
  conv=list('A'='T',
             'T'='A',
             'C'='G',
             'G'='C')
  return(unlist(conv[x]))
}
#Interaction with bash####
#bgzip
#compressed in bed.gz like table
#input: x: data.table with 1: chr, 2: start, 3: end
#bgz_file=  bed.gz filename
#output; bed.gz file compressed with bgzip and indexed with tabix
bgzip<-function(x,bgz_file,sort_coord=FALSE,col.names=TRUE,add_header_of=NULL){
  require('data.table')
  bed_file<-tools::file_path_sans_ext(bgz_file)
  if(sort_coord)
    setorderv(x,cols = colnames(x)[1:2],order = 1)
  
  if(!is.null(add_header_of)){
    append=TRUE
    col.names=FALSE
    #extract header from vcf file
    cmd=paste('bcftools head', add_header_of,'>',bed_file)
    system(cmd)
  }
  fwrite(x,bed_file,sep='\t',col.names = col.names,append = append)
  system(paste('bgzip -f',bed_file))
  system(paste('tabix',bgz_file))
  system(paste('rm -f',bed_file))
  message(bgz_file,' created')
  return(bgz_file)
  
}

