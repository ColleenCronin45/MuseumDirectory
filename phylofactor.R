## phylfactor museum samples
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(ape)
library(caper)
library(phylofactor)
library(data.table)
library(treeio)
library(ggtree)

## load data
setwd("~/Desktop/MuseumDirectory")
data1=read.csv('SpeciesFrequency.csv')
data2=read.csv('SpeciesFTBfrequency3.4.22.csv')

## combinee
data1$order=NULL
data=merge(data1,data2,by="species")
rm(data1,data2)

## load taxonomy and phy
setwd("~/Desktop/BeckerLabOU/phylos")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
taxonomy=read.csv("taxonomy_mamPhy_5911species.csv")

## subset to orders
unique(data$order)
taxonomy=taxonomy[taxonomy$ord%in%c("CHIROPTERA","RODENTIA"),]

## fix names
taxonomy$tip=sapply(strsplit(taxonomy$tiplabel,'_'),function(x) paste(x[1],x[2],sep='_'))
data$tip=gsub(" ","_",data$species)

## mark 0/1 for sampled
tdata=taxonomy
tdata$sampled=ifelse(tdata$tip%in%data$tip,1,0)

## subset data
data$tree=ifelse(data$tip%in%taxonomy$tip,1,0)
table(data$tree)
data=data[data$tree==1,]

## make simple names
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))

## trim phylo
rtree=tree
tree=keep.tip(rtree,data$tip)
atree=keep.tip(rtree,taxonomy$tip)

## merge
data=merge(data,taxonomy,by="tip")

## match
bdata=data[match(tree$tip.label,data$tip),]
adata=tdata[match(atree$tip.label,tdata$tip),]

## save
bdata$label=bdata$tip
bdata$Species=bdata$tip
adata$label=adata$tip
adata$Species=adata$tip

## merge
cdata=comparative.data(phy=tree,data=bdata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)
adata=comparative.data(phy=rtree,data=adata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)

## new tip
cdata$data$tip=cdata$data$label
adata$data$tip=adata$data$label

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
pfsum=function(pf){
  
  ## get formula
  chars=as.character(pf$frmla.phylo)[-1]
  
  ## response
  resp=chars[1]
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## save model
  model=chars[2]
  
  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$Species%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}

## taxonomy
cdata$data$taxonomy=paste(cdata$data$ord,cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')
adata$data$taxonomy=paste(adata$data$ord,adata$data$fam,adata$data$gen,adata$data$Species,sep='; ')

## set taxonomy
taxonomy=data.frame(adata$data$taxonomy)
names(taxonomy)="taxonomy"
taxonomy$Species=rownames(adata$data)
taxonomy=taxonomy[c("Species","taxonomy")]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## binary
set.seed(1)
bpf=gpf(Data=adata$data,tree=adata$phy,
        frmla.phylo=sampled~phylo,
        family=binomial,algorithm='phylo',nfactors=10,min.group.size=5)

## summarize
bpf_results=pfsum(bpf)$results

## negbin model fcn
model.fcn2 <- function(formula,data,...){
  fit <- tryCatch(MASS::glm.nb(formula,data,...),
                  error=function(e) NA)
  #fit <- do.call
  return(fit)
}

## negbin objective function
obj.fcn2 <- function(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...){
  #if (!'negbin' %in% class(fit) & !'glm' %in% class(fit) & !'lm' %in% class(fit))
  if (!'negbin' %in% class(fit))
  {
    return(0)
  }
  else 
  {
    #fit2 <- MASS::glm.nb(Z.poisson~1,data = fit$model)
    fit$null.deviance-fit$deviance %>% return()
    #fit$twologlik %>% return()
  }
}

## for tissue
set.seed(1)
tpf=gpf(Data=cdata$data,tree=cdata$phy,
        frmla.phylo=tissUbuff~phylo,
        model.fcn = model.fcn2,objective.fcn = obj.fcn2,
        cluster.depends='library(MASS)',
        algorithm='phylo',nfactors=3,min.group.size=5)

## as binomial
cdata$data$pos=cdata$data$tissUbuff
cdata$data$neg=cdata$data$n-cdata$data$pos
bat=cdata[cdata$data$ord=="CHIROPTERA",]
set.seed(1)
tpf=gpf(Data=bat$data,tree=bat$phy,
        frmla.phylo=cbind(pos,neg)~phylo,
        family=binomial,algorithm='phylo',nfactors=3,min.group.size=5)

