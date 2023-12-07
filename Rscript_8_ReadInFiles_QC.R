
#=== LOAD LIBRARIES AND DEFAULTS

library(tcR)
library(ggplot2)
library(Hmisc)
library(survival)
library(ggfortify)
library(pander)
library(ggpubr)
library(survminer)
library(formattable)
library(ggExtra)
library(corrplot)
library(pheatmap)
library(Kendall)
library(knitr)
library(kableExtra)
library(factoextra)
library(dplyr)
library(cowplot)
library(tidyverse)
library(caret)
library(ModelGood)
library(pROC)
library(mlbench)
library(randomForest)
library(data.table)

myplot_parameters <-function() { font("title", size = 30) +
    font("xlab", size = 35) +
    font("ylab", size = 35) +
    font("x.text",size = 30) +
    font("y.text",size = 30) +
    font("legend.title", size = 25) + 
    font("legend.text", size = 25) }

set.seed(41)

#=== SET YOUR WORKING DIRECTORY 

setwd("Set_WD")

#=== READ TCR REPERTOIRE MIXCR FILES
#=== READ PHENOTYPE FILE

dir.path="../EBV_Rcode/data/"
pheno.file="../EBV_Rcode/data/phenotype_hla.txt"

suffix.name="_clones_full.txt"
repertoire.tool="mixcr"

combined.list<-list()
list.filenames<-list.files(path=dir.path, pattern = paste0("*",suffix.name,"$",sep=""), recursive = FALSE, full.names = TRUE)
if(repertoire.tool=="mixcr"){
  combined.list = lapply(list.filenames, parse.mixcr)
}
list.filenames<-sub(".*/","",list.filenames)
names(combined.list)<-sub(suffix.name,"",list.filenames)

pheno.orig<-read.table(pheno.file,header=T,sep="\t",stringsAsFactors = F,na.strings=c(""," ","NA"))
pheno<-pheno.orig[,c("ID","response")]
phenoR<-pheno[pheno$response=="R" & !is.na(pheno$response),]
phenoNR<-pheno[pheno$response=="NR" & !is.na(pheno$response),]


combined.list.norm<-lapply(combined.list, function(x) {
  total.reads<-sum(x$Read.count); 
  rc1<-x$Read.count; 
  newreads<-(x$Read.count/total.reads)*10^6; 
  x$Read.count<-newreads; 
  total.reads<-sum(x$Read.count); 
  x$Read.proportion<-x$Read.count/total.reads; 
  return(x)})

for(i in 1:length(combined.list.norm)){
  combined.list.norm[[i]]<-combined.list.norm[[i]][order(combined.list.norm[[i]]$Read.count,decreasing = T),]
}


#=== QC INPUT FILE

combined.list.R<- get.inframes(combined.list.norm[phenoR$ID])
combined.list.NR<- get.inframes(combined.list.norm[phenoNR$ID])

#=== EXTRACT TOP 10% CLONES
#=== USE ONLY CLONES THAT START WITH C AND WITH F

perc=(10/100)

if(!is.na(perc)) {
  combined.list.Rtemp<-lapply(combined.list.R,function(x){ 
    new_rows<-round(nrow(x)*perc); 
    newx<-x[1:new_rows,];
    newx<-newx[grep("^C[A-Z]*F$",newx$CDR3.amino.acid.sequence),]; 
    newx[grep("[*, ~]", newx$CDR3.amino.acid.sequence, invert = T),]; 
    return(newx) } )
  combined.list.NRtemp<-lapply(combined.list.NR,function(x){ 
    new_rows<-round(nrow(x)*perc); 
    newx<-x[1:new_rows,]; 
    newx<-newx[grep("^C[A-Z]*F$",newx$CDR3.amino.acid.sequence),]; 
    newx[grep("[*, ~]", newx$CDR3.amino.acid.sequence, invert = T),]; 
    return(newx) } )
}


combined.list.R<-combined.list.Rtemp
combined.list.NR<-combined.list.NRtemp

