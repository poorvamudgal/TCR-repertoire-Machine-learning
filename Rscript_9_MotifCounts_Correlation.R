
#====== CREATE THE R TO NR RATIO OF MOTIFS
#=== KMERS LENGTH OF 2 TO 4
#=== REMOVE 3 LEADING AND 3 TRAILING AMINO ACIDS


#=== NON RESPONSE SAMPLES KMER COUNTS

tableNR<-c()
for(i in 2:4) { 
  tableNR<-rbind(tableNR,get.kmers(combined.list.NR, .head = -1, .k = i, .clean = F, .meat = F, .verbose = F, .left.shift = 3, .right.shift = 3))
}
tableNR[tableNR==0]<-NA

#=== RESPONSE SAMPLES KMER COUNTS

tableR<-c()
for(i in 2:4) {
  tableR<-rbind(tableR,get.kmers(combined.list.R, .head = -1, .k = i, .clean = F, .meat = F, .verbose = F, .left.shift = 3, .right.shift = 3))
}
tableR[tableR==0]<-NA


temp<-sapply(combined.list.NR, '[[', "CDR3.amino.acid.sequence")
temp2<-unique(sort(unlist(temp,use.names = FALSE)))
cdr3.list.NR<-str_sub(temp2,start=4,end=-4)

temp<-sapply(combined.list.R, '[[', "CDR3.amino.acid.sequence")
temp2<-unique(sort(unlist(temp,use.names = FALSE)))
cdr3.list.R<-str_sub(temp2,start=4,end=-4)

#=== MERGE TABLES TO CREATE ONE LIST OF KMERS

table<-merge(tableR,tableNR,by="Kmers",all=T)
table<-table[,-c(3:ncol(table))]
table$motifR<-NA
table$motifNR<-NA
table$motifR<-sapply(table$Kmers,function(x) {y<-length(grep(x,cdr3.list.R,value = TRUE)) } )
table$motifNR<-sapply(table$Kmers,function(x) {y<-length(grep(x,cdr3.list.NR,value = TRUE)) } )

#=== GET THE RATIO OF KMER COUNTS R/NR

table$ratio<-table$motifR/table$motifNR
table[,2]<-NULL
table[table=="Inf"]<-0

#=== NON RESPONSE SAMPLES KMER ACCUMULATIVE CDR3

tableNR.cloneNum<-c()
for(i in 2:4) {
  tableNR.cloneNum<-rbind(tableNR.cloneNum,get.kmers(combined.list.NR, .head = -1 ,.k = i, .clean = F, .meat = T, .verbose = F, .left.shift = 3, .right.shift = 3))
}
tableNR.cloneNum[tableNR.cloneNum==0]<-NA
rownames(tableNR.cloneNum)<-tableNR.cloneNum$Kmers
tableNR.cloneNum$Kmers<-NULL
tableNR.cloneNum<-log(tableNR.cloneNum)

#=== RESPONSE SAMPLES KMER ACCUMULATIVE CDR3

tableR.cloneNum<-c()
for(i in 2:4) {
  tableR.cloneNum<-rbind(tableR.cloneNum,get.kmers(combined.list.R, .head = -1 ,.k = i, .clean = F, .meat = T, .verbose = F, .left.shift = 3, .right.shift = 3))
}

tableR.cloneNum[tableR.cloneNum==0]<-NA
rownames(tableR.cloneNum)<-tableR.cloneNum$Kmers
tableR.cloneNum$Kmers<-NULL
tableR.cloneNum<-log(tableR.cloneNum)


#=== CALCULATE CORRELATION OF RESPONSE SAMPLES WITH OVERALL SURVIVAL FOR EACH KMER
#=== SPEARMAN CORRELATION

tableR.cloneNum.T<-data.frame(t(tableR.cloneNum))
tableR.cloneNum.T$Sample<-rownames(tableR.cloneNum.T)

pheno<-pheno.orig[!is.na(pheno.orig$response),c("ID","OS")]
colnames(pheno)<-c("Sample","OS")

tableR.cloneNum.pheno<-merge(tableR.cloneNum.T,pheno,by="Sample")
rownames(tableR.cloneNum.pheno)<-tableR.cloneNum.pheno$Sample
tableR.cloneNum.pheno$Sample<-NULL
names(tableR.cloneNum.pheno)[names(tableR.cloneNum.pheno) == 'NA.'] <- 'NA'

#=== ADD DUMMY ROW

cor.dt=data.frame("ABC",0.87,0.001,stringsAsFactors = F)
colnames(cor.dt)<-c("Kmers","Cor","P")


n=1
for( var in rownames(tableR.cloneNum)) {
  temp<-tableR.cloneNum.pheno[,names(tableR.cloneNum.pheno) %in% c(var,"OS")]
  if(nrow(temp[!is.na(temp[,var]),])>2) {
    res <- cor.test(temp[,var], temp[,"OS"], method = "spearman")
    cor.dt[n,]<-c(var,res$estimate,res$p.value)
  }
  else{
    cor.dt[n,]<-c(var,NA,NA)
  }
  n=n+1
}

cor.dt[,2:3] <- sapply(cor.dt[,2:3], as.numeric)


#=== MERGE R/NR RATIO AND CORRELATION TABLES
#=== SAVE RATIO TABLE

table.cor<-merge(table,cor.dt)
save(table.cor,file = "RtoNRratio_Rcorrelation.RData")

#=== SAVE MOTIF CDR3 COUNTS TABLE

table.num<-merge(tableR,tableNR,by="Kmers",all=T)
save(table.num,file = "motif_cdr3_counts.RData")

#=== #save MOTIF ACCUMULATIVE CDR3 TABLE
tableR.cloneNum$Kmers<-rownames(tableR.cloneNum)
tableNR.cloneNum$Kmers<-rownames(tableNR.cloneNum)

table.cloneNum<-merge(tableR.cloneNum,tableNR.cloneNum,by="Kmers",all=T)
save(table.cloneNum,file = "motif_accu_cdr3_counts.RData")

