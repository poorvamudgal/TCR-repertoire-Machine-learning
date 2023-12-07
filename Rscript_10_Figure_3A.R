

#====== MACHINE LEARNING, RANDOM FOREST AND FIGURE 3A

#=== INPUT FILE (IF REQUIRED)
load("../Set_WD/RtoNRratio_Rcorrelation.RData")
load("../Set_WD/motif_cdr3_counts.RData")
pheno.file="../EBV_Rcode/data/phenotype_hla.txt"
pheno.orig<-read.table(pheno.file,header=T,sep="\t",stringsAsFactors = F,na.strings=c(""," ","NA"))

#=== GET THE R/NR RATIO RANGE FOR SELECTING MOTIFS

cutoff.range.orig<-unique(sort(round(table.cor$ratio)))
cutoff.range<-cutoff.range.orig[round(length(cutoff.range.orig)/2):length(cutoff.range.orig)]

pheno<-pheno.orig[!is.na(pheno.orig$response),c("ID","response")]
colnames(pheno)<-c("Sample","response")

#=== ADD DUMMY ROW

accuracy<-data.frame(1,2,3,stringsAsFactors = F)
colnames(accuracy)<-c("R/NR cutoff","accuracy","number of motif")

#=== RUN RANDOM FOREST OF DIFFERENT SET OF MOTIFS
#=== MOTIFS HAVING RATIO GREATER THAN CUTOFF VALUES AND POSITIVE CORRELATION OF RESPONSE WITH OVERALL SURVIVAL
#=== REMOVE MOTIFS THAT ARE HIGHLY CORRELATED WITH EACH OTHER

n=1
for(cutoff in cutoff.range){
  pos_cor_motif<-table.cor[table.cor$ratio > cutoff & !is.na(table.cor$ratio) & table.cor$Cor > 0 ,"Kmers"]
  if(length(pos_cor_motif)>1) {
    df<-table.num[table.num$Kmers %in% c(pos_cor_motif),]
    rownames(df)<-df$Kmers
    df$Kmers<-NULL
    df2<-as.data.frame.matrix(t(apply(df,2,function(x){ x/sum(x,na.rm = T) })))
    
    df2$Sample<-rownames(df2)
    df3<-merge(df2,pheno,by="Sample",all.x=T)
    #rownames(df3)<-paste0(df3$Sample,df3$response)
    rownames(df3)<-df3$ID
    df3[is.na(df3)]<-runif(1,1e-20,1e-5)
    
    correlationMatrix <- cor(df3[,-c(1,ncol(df3))], method = "spearman")
    highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7)
    
    if(length(highlyCorrelated)>0){
      df3<-df3[,-c(highlyCorrelated+1)]
    }
    
    df3<-df3[,!names(df3) %in% c("Sample")]
    
    model <- train(
      response ~., data = df3, method = "rf",
      trControl = trainControl("repeatedcv", number = 10, repeats = 10),
      importance = TRUE
    )
    temp<-model$results
    accuracy[n,]<-c(cutoff,temp[nrow(temp),"Accuracy"],length(pos_cor_motif))
    
    n=n+1
  }
}
accuracy$`R/NR cutoff`<-as.numeric(accuracy$`R/NR cutoff`)
accuracy$accuracy<-as.numeric(accuracy$accuracy)

#=== FIGURE 3A

ggplot(accuracy) + 
  geom_point(aes(`R/NR cutoff`,accuracy),shape=19,size=2) + 
  ylim(0.5,1) +
  scale_x_continuous(breaks=c(cutoff.range)) +
  labs(x="CDR3 ratio of R/NR",y="Model accuracy") +
  font("x.text",size = 20) +
  theme_pubr(border = FALSE) +
  myplot_parameters() 

ggsave("3A.Accuracy.vs.RtoNR.cdr3.cutoff.svg",device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.8)
