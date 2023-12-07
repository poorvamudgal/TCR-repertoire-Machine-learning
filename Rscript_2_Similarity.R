

#====== CREATE SIMILARITY HLA TABLE FOR PLOTS
#=== USE 0.1 PERC, 1 PERC , 10 PERC AND ALL CLONES FOR CALCULATING SIMILARITY
#=== AMINO ACID CLONES JACCARD SIMILARITY
#=== FOR CDR3 AND FULL TCR (CDR3 + VGENE)

#=== ADD DUMMY ROW

similarity.table.byPerc=data.frame("CTL000","CTL000","3","12","R","0.1","nucT",stringsAsFactors = F)
colnames(similarity.table.byPerc)<-c("ID2","ID","Similarity","OS","Resp","Perc.clones","method")

#=== ONLY INFRAME CLONES

combined.list.in <- get.inframes(combined.list.norm)

for (perc in c(0.001, 0.01, 0.1, 1)) {
  
  combined.list_per<-lapply(combined.list.in,function(x){ new_rows<-nrow(x)*perc; newx<-x[1:new_rows,]; return(newx) } )
  
  #=== NON RESPONSE SAMPLES FROM PHENOTYPE FILE AND TCR REPERTOIRE
  nresp<-pheno[pheno$response=="NR" & !is.na(pheno$response),]
  combined.list.nresp_per<-combined.list_per[nresp$ID]
  
  
  #=== RESPONSE SAMPLES FROM PHENOTYPE FILE AND TCR REPERTOIRE
  resp<-pheno[pheno$response=="R" & !is.na(pheno$response),]
  combined.list.resp_per<-combined.list_per[resp$ID]
  
  method<-c("jaccard")
  seq<-c("aa")
  vgene.val<-c("F","T")
  
  for (i in method)
  {
    for (j in seq)
    {
      for (k in vgene.val)
      {
        
        
        output.nresp<-repOverlap(combined.list.nresp_per, i, j , .vgene = k, .verbose = F, .norm=F, .do.unique = T)
        
        output.resp<-repOverlap(combined.list.resp_per, i, j , .vgene = k, .verbose = F, .norm=F, .do.unique = T)
        
        k2<- ifelse(k=="T",".With.Vgene",".Without.Vgene")
        
        simi.nresp<-as.data.frame(as.table(output.nresp))
        colnames(simi.nresp)<-c("Var1","ID","Freq")
        simi.nresp$OS<-rep(nresp$OS, times=1, each=nrow(nresp))
        simi.nresp$Resp<-"NR"
        simi.nresp$perc<-perc
        simi.nresp$method<-paste0(j,k2)
        
        simi.resp<-as.data.frame(as.table(output.resp))
        colnames(simi.resp)<-c("Var1","ID","Freq")
        simi.resp$OS<-rep(resp$OS, times=1, each=nrow(resp))
        simi.resp$Resp<-"R"
        simi.resp$perc<-perc
        simi.resp$method<-paste0(j,k2)
        
        colnames(simi.nresp)<-c("ID2","ID","Similarity","OS","Resp","Perc.clones","method")
        colnames(simi.resp)<-c("ID2","ID","Similarity","OS","Resp","Perc.clones","method")
        similarity.table.byPerc<-rbind(similarity.table.byPerc,simi.nresp)
        similarity.table.byPerc<-rbind(similarity.table.byPerc,simi.resp)
        
      }
    }
  }
  
  cat("\n\n")
}
#=== REMOVE DUMMY ROW

similarity.table.byPerc<-similarity.table.byPerc[-1,]
similarity.table.byPerc$Similarity<-as.numeric(similarity.table.byPerc$Similarity)

#====== HLA COUNTING

#=== CREATE HLA A ALLELE TABLE

newpheno<-pheno[!is.na(pheno$response),]

del1<-select(newpheno,-c("HLA.A2"))
del2<-select(newpheno,-c("HLA.A1"))
colnames(del2)<-colnames(del1)
newpheno1<-rbind(del1,del2)

newpheno1<-newpheno1[!is.na(newpheno1$HLA.A1),]

hla.a1.table<-table(newpheno1$ID,newpheno1$HLA.A1)

#=== CREATE HLA B ALLELE TABLE

del1<-select(newpheno,-c("HLA.B2"))
del2<-select(newpheno,-c("HLA.B1"))
colnames(del2)<-colnames(del1)

newpheno1<-rbind(del1,del2)
newpheno1<-newpheno1[!is.na(newpheno1$HLA.B1),]

hla.b1.table<-table(newpheno1$ID,newpheno1$HLA.B1)

#=== CREATE HLA C ALLELE TABLE

del1<-select(newpheno,-c("HLA.C2"))
del2<-select(newpheno,-c("HLA.C1"))
colnames(del2)<-colnames(del1)

newpheno1<-rbind(del1,del2)
newpheno1<-newpheno1[!is.na(newpheno1$HLA.C1),]

hla.c1.table<-table(newpheno1$ID,newpheno1$HLA.C1)

#=== MERGE ALL HLA A B AND C TABLES

hla.table<-cbind(hla.a1.table,hla.b1.table,hla.c1.table)

#=== ADD DUMMY ROW

hla.table.count=data.frame("CTL000","CTL000","3",stringsAsFactors = F)
colnames(hla.table.count)<-c("ID","ID2","HLA.count")

#=== COUNT SHARED HLA BETWEEN SAMPLES
#=== SAMPLES WITH COUNT 2 (e.g. A02;A02 FOR BOTH IDS) ARE COUNTED AS 1 SHARED ALLELE

for(i in 1:nrow(hla.table)) {
  for(j in 1:nrow(hla.table)) {
    count=0
    for(k in 1:ncol(hla.table)) {
      if(hla.table[i,k]==1 & hla.table[j,k]==1) { count=count+1}
      if(hla.table[i,k]==1 & hla.table[j,k]==2) { count=count+1}
      if(hla.table[i,k]==2 & hla.table[j,k]==1) { count=count+1}
      if(hla.table[i,k]==2 & hla.table[j,k]==2) { count=count+1}
    }
    if(rownames(hla.table)[i]==rownames(hla.table)[j]){count=0} 
    
    temp<-data.frame(rownames(hla.table)[i],rownames(hla.table)[j],count,stringsAsFactors = F)
    colnames(temp)<-c("ID","ID2","HLA.count")
    hla.table.count<-rbind(hla.table.count,temp)
    
  }
}
#=== REMOVE DUMMY ROW

hla.table.count<-hla.table.count[-1,]

similarity.hla<-merge(similarity.table.byPerc,hla.table.count,by.x=c("ID","ID2"))

#=== SAVE SIMILARITY HLA TABLE AS RDATA
save(similarity.hla,file = "similarity.hla.RData")

