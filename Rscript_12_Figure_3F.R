

#=== Figure 3F
#=== ROC CURVE FOR SPFS AND SPDQ
#=== USING CDR3 COUNTS TABLE

#=== INPUT FILE (IF REQUIRED)
load("../Set_WD/motif_cdr3_counts.RData")
pheno.file="../EBV_Rcode/data/phenotype_hla.txt"
pheno.orig<-read.table(pheno.file,header=T,sep="\t",stringsAsFactors = F,na.strings=c(""," ","NA"))

pheno<-pheno.orig[!is.na(pheno.orig$response),c("ID","response")]
colnames(pheno)<-c("Sample","response")

rownames(table.num)<-table.num$Kmers
table.num$Kmers<-NULL
temp<-as.data.frame(t(table.num))

temp<-temp[,names(temp) %in% c("SPDQ","SPFS","DGAG","EES","EVAG","GSRS","KTGE","LADN","LYL","SERR","TTNT")]
temp2<-as.data.frame.matrix(t(apply(temp,1,function(x){x/sum(x,na.rm=T)})))
temp2$Sample<-rownames(temp2)
temp3<-merge(temp2,pheno,by="Sample")
temp3[is.na(temp3)]<-runif(1,1e-20,1e-5)

for(var in c("SPFS","SPDQ")) {
  
  forPlots<-temp3[,names(temp3) %in% c(var,"response")]
  fmla1 <- as.formula(paste0("response ~ ", var))
  forPlots$response<-as.factor(forPlots$response)
  
  #=== RUN GENERALIZED LINEAR MODEL ON HIGHLY SIGNIFICANT MOTIFS
  
  model <- glm(fmla1,data=forPlots,family=binomial(link="logit"))
  
  #=== GET ROC CURVE BETWEEN TO ACTUAL AND PREDICTED RESPONSE VALUES
  
  res.roc <- roc(forPlots$response,model$fitted.values)
  roc.dt<-data.frame(res.roc$sensitivities,res.roc$specificities,stringsAsFactors = F)
  roc.dt<-roc.dt[order(roc.dt$res.roc.sensitivities),]
  SE<-as.data.frame(ci.se(res.roc,specificities = roc.dt$res.roc.specificities))
  
  print(ggplot(roc.dt) + 
          geom_ribbon(aes(ymin = SE$`2.5%`, ymax = SE$`97.5%`,x=1-res.roc.specificities), fill = "gray95") + 
          geom_line(aes(x=1-res.roc.specificities,y=res.roc.sensitivities),size=2) + 
          geom_text(x=0.8,y=0.5,label=paste0("AUC :",round(res.roc$auc,digits = 3)),size=7,family = "Arial") +
          geom_abline(intercept = 0, slope = 1) +
          ggtitle(var) +
          labs(x="1-Specificity",y="Sensitivity") +
          theme_pubr(border = FALSE) +
          myplot_parameters() +
          font("title",size = 13) +
          font("x.text",size = 18))
  
  ggsave(paste0("3F.",var,".ROC.unique.cdr3.counts.svg"),device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.8)
  
}