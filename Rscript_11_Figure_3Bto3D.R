
#====== FIGURE 3B 3C AND 3D
#=== SELECT THE R/NR RATIO CUTOFF THAT GAVE HIGHEST ACCURACY IN PREDICTING R AND NR GROUPS

#=== INPUT FILE (IF REQUIRED)
load("../Set_WD/RtoNRratio_Rcorrelation.RData")
load("../Set_WD/motif_cdr3_counts.RData")
pheno.file="../EBV_Rcode/data/phenotype_hla.txt"
pheno.orig<-read.table(pheno.file,header=T,sep="\t",stringsAsFactors = F,na.strings=c(""," ","NA"))

pheno<-pheno.orig[!is.na(pheno.orig$response),c("ID","response")]
colnames(pheno)<-c("Sample","response")

cutoff<-data.frame(accuracy[order(accuracy$accuracy,decreasing = T),])

#=== FIGURE 3B
#=== GET THE LIST OF MOTIFS THAT HAVE R/NR RATIO GREATER THAN THE CUTOFF VALUE

final_motifs<-table.cor[table.cor$ratio > cutoff[1,c("R.NR.cutoff")] & !is.na(table.cor$ratio) & table.cor$Cor > 0 ,"Kmers"]

table.cor$highlight<-ifelse(table.cor$Kmers %in% final_motifs,"Possible motif","ALL")


ggplot(table.cor, aes(Cor, ratio))+
  stat_binhex (data = table.cor[table.cor$highlight=="ALL",],color="lightgray",fill="lightgray") +
  geom_point (data = table.cor[table.cor$highlight=="Possible motif",], aes(Cor, ratio),size=2,color="black") +
  ylim(0,16) + xlim(-1,1) +
  labs(y = "CDR3 ratio of R/NR", x = "Correlation") +
  theme_pubr(border = FALSE) +
  myplot_parameters()

ggsave("3B.Correlation.vs.RtoNR.cdr3.ratio.svg",device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.8)


df<-table.num[table.num$Kmers %in% c(final_motifs),]
rownames(df)<-df$Kmers
df$Kmers<-NULL

df2<-as.data.frame.matrix(t(apply(df,2,function(x){ x/sum(x,na.rm=T) })))
df2$Sample<-rownames(df2)
df3<-merge(df2,pheno,by="Sample",all.x=T)
#rownames(df3)<-paste0(df3$Sample,df3$response)
rownames(df3)<-df3$Sample

#=== REMOVE MOTIFS THAT ARE HIGHLY CORRELATED WITH OTHER MOTIFS

correlationMatrix <- cor(df3[,-c(1,ncol(df3))], method = "spearman")
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7)

if(length(highlyCorrelated)>0){
  df3<-df3[,-c(highlyCorrelated+1)]
}
df3<-df3[order(df3$response),]

df3[is.na(df3)]<-runif(1,1e-20,1e-5)

df3<-df3[,!names(df3) %in% c("Sample")]

#=== RUN RANDOM FOREST TO GET IMPORTANCE VALUES

model <- train(
  response ~., data = df3, method = "rf",
  trControl = trainControl("repeatedcv", number = 10, repeats = 10),
  importance = TRUE
)
importance.dt<-as.data.frame(importance(model$finalModel))
importance.dt$Motif<-rownames(importance.dt)

#=== FIGURE 3D
#=== MOTIF PLOT USING ACCURACY AND GINI VALUES

ggdotchart(importance.dt,x="Motif",y="MeanDecreaseAccuracy",sorting="descending",rotate = TRUE, dot.size = 2, color = rep(c("lightgray","black"),times=c(9,2)), ylab="Mean Decrease Accuracy")+
  myplot_parameters()+
  font("x.text",size = 20) +
  font("y.text",size = 20) +
  font("xlab",size = 20) 

ggsave("3D.uniqueCDR3.meanDecreaseAccuracy.svg",device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.6)

ggdotchart(importance.dt,x="Motif",y="MeanDecreaseGini",sorting="descending",rotate = TRUE, dot.size = 2, color = rep(c("lightgray","black"),times=c(9,2)), ylab="Mean Decrease Gini" ) +
  myplot_parameters() +
  font("x.text",size = 20) +
  font("y.text",size = 20) +
  font("xlab",size = 20) 

ggsave("3D.uniqueCDR3.meanDecreaseGini.svg",device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.6)

#=== FIGURE 3C
#=== UNSUPERVISED CLUSTERING
#=== MANHATTAN DISTANCE CALCULATION AND WARD D2 CLUSTERING

df3<-df3[,!names(df3) %in% c("response")]
res.dist <- dist(df3, method = "manhattan")
res.hc <- hclust(d = res.dist, method = "ward.D2")

final_motifs<-names(df3)

color_ids_vt<-rep(c("tomato","slategrey","tomato","slategrey","tomato","slategrey"),times=c(15,6,2,1,1,6))
fviz_dend(res.hc, cex = 1.8, k=2, type = "circular",label_cols = color_ids_vt, k_colors = c("tomato","slategrey"),lwd = 1) 

ggsave("3C.dendrogram.svg",device = "svg",width = 12,height = 12,dpi = "retina", scale = 0.6)

