
#=== FIGURE 3E 3G AND 3H
#=== BOXPLOT OF LOG ACCUMULATIVE CDR3 IN R AND NR SAMPLES
#=== SCATTERPLOT OF LOG ACCUMULATIVE CDR3 COUNTS AND OVERALL SURVIVAL
#=== BOXPLOT OF LOG ACCUMULATIVE CDR3 VS HLA ALLELES IN R AND NR SAMPLES

#=== INPUT FILE (IF REQUIRED)
load("../Set_WD/motif_accu_cdr3_counts.RData")
pheno.file="../EBV_Rcode/data/phenotype_hla.txt"
pheno.orig<-read.table(pheno.file,header=T,sep="\t",stringsAsFactors = F,na.strings=c(""," ","NA"))

myplot_parameters <-function() { font("title", size = 15) +
    font("xlab", size = 20) +
    font("ylab", size = 20) +
    font("x.text",size = 15) +
    font("y.text",size = 20) +
    font("legend.title", size = 15) + 
    font("legend.text", size = 15) }

pheno<-pheno.orig[!is.na(pheno.orig$response),c("ID","response","OS","HLA.A1","HLA.A2","HLA.B1","HLA.B2")]
colnames(pheno)<-c("Sample","response","OS","HLA.A1","HLA.A2","HLA.B1","HLA.B2")

rownames(table.cloneNum)<-table.cloneNum$Kmers
table.cloneNum$Kmers<-NULL

temp<-as.data.frame(t(table.cloneNum))
#temp[is.na(temp)]<-0
temp$Sample<-rownames(temp)
temp2<-merge(temp,pheno,by="Sample")


for(var in c("SPDQ","SPFS")) {
  
  forPlots<-temp2[,names(temp2) %in% c(var,"OS","response")]
  
  #=== FIGURE 3E
  
  print(ggboxplot(forPlots,x = "response", y = var, 
                  color = "response", palette = c("slategrey","tomato"), 
                  add = "jitter", size=1, order = c("NR","R"),
                  ylab = "log Accumulative CDR3", xlab = "Response",title = var)+ 
          stat_compare_means(label="p.format", hide.ns = TRUE, size = 10,label.y = max(forPlots[,var],na.rm=T)+1, label.x = 1.25) + 
          theme(legend.position="none") +
          expand_limits(x = 0, y = 0) +
          myplot_parameters() +
          scale_y_continuous(limit = c(0, max(forPlots[,var],na.rm=T)+2)) )
  
  ggsave(paste0("3E.",var,".resp.vs.logclonecount.svg"),device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.6)
  
  #=== FIGURE 3G
  
  print(ggscatter(forPlots,x = "OS", y = var,
                  palette = c("slategrey","tomato"), facet.by = "response", color="response",
                  ylab = "log Accumulative CDR3", xlab = "Overall Survival (weeks)",
                  title = var,
                  add = "reg.line", cor.method = "spearman",
                  conf.int=TRUE,conf.int.level = 0.95,
                  cor.coef = TRUE, cor.coef.size = 5,
                  space = "free", scale="free_x",
                  add.params = list(linetype="dashed",fill = "gray95")) +
          scale_y_continuous(limit = c(0, max(forPlots[,var],na.rm=T)+1)) +
          theme(strip.text.x = element_text(size = 20), legend.position="none") +
          myplot_parameters() 
  )
  
  ggsave(paste0("3G.",var,".OS.vs.logclonecount.svg"),device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.6)
  
  #=== FIGURE 3H
  
  forHla<-temp2[,c("Sample","response","HLA.A1","HLA.A2",var)]
  del1<-select(forHla,-c("HLA.A2"))
  del2<-select(forHla,-c("HLA.A1"))
  colnames(del2)<-colnames(del1)
  
  forHla.v2<-rbind(del1,del2)
  
  forHla.v3<-unique(forHla.v2[!is.na(forHla.v2[,var]),])
  forHla.v3.A<-unique(forHla.v3[!is.na(forHla.v3[,"HLA.A1"]),])
  
  
  forHla<-temp2[,c("Sample","response","HLA.B1","HLA.B2",var)]
  del1<-select(forHla,-c("HLA.B2"))
  del2<-select(forHla,-c("HLA.B1"))
  colnames(del2)<-colnames(del1)
  
  forHla.v2<-rbind(del1,del2)
  
  forHla.v3<-unique(forHla.v2[!is.na(forHla.v2[,var]),])
  forHla.v3.B<-unique(forHla.v3[!is.na(forHla.v3[,"HLA.B1"]),])
  colnames(forHla.v3.B)<-c("Sample","response","HLA.A1",var)
  
  hlaA.hlaB<-rbind(forHla.v3.A,forHla.v3.B)
  hlaA.hlaB<-hlaA.hlaB[hlaA.hlaB$HLA.A1 %in% c("A*02","A*11","B*46"),]
  
  print(ggboxplot(hlaA.hlaB,x = "HLA.A1", y = var, 
                  color = "response", palette = c("slategrey","tomato"), 
                  add = "jitter",
                  ylab = "log Accumulative CDR3", xlab = "HLA allele",title = var)+ 
          stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, size = 15,label.y = max(forHla.v3[,var])+1) + 
          expand_limits(x = 0, y = 0) +
          myplot_parameters() +
          theme(legend.title=element_blank()) +
          scale_y_continuous(limit = c(0, max(forHla.v3[,var],na.rm=T)+1.5)) )
  
  ggsave(paste0("3H.",var,".hla.alleles.vs.logclonecount.svg"),device = "svg",width = 8,height = 10,dpi = "retina", scale = 0.8)
  
}