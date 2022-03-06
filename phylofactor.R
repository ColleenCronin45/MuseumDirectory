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
tree=keep.tip(tree,data$tip)

## merge
data=merge(data,taxonomy,by="tip")

## match
bdata=data[match(tree$tip.label,data$tip),]

## save
bdata$label=bdata$tip
bdata$Species=bdata$tip

## merge
cdata=comparative.data(phy=tree,data=bdata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)

## new tip
cdata$data$tip=cdata$data$label
