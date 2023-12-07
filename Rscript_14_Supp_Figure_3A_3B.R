
#=== SUPPLEMENTARY FIGURE 3A AND 3B
#=== BOXPLOT OF LOG ACCUMULATIVE CDR3 IN R AND NR SAMPLES
#=== SCATTERPLOT OF LOG ACCUMULATIVE CDR3 COUNTS AND OVERALL SURVIVAL

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

pheno<-pheno.orig[!is.na(pheno.orig$response),c("ID","response","OS")]
colnames(pheno)<-c("Sample","response","OS")

rownames(table.cloneNum)<-table.cloneNum$Kmers
table.cloneNum$Kmers<-NULL

temp<-as.data.frame(t(table.cloneNum))
#temp[is.na(temp)]<-0
temp$Sample<-rownames(temp)
temp2<-merge(temp,pheno,by="Sample")

for(var in c(c("DGAG","EES","EVAG","GSRS","KTGE","LADN","LYL","SERR","TTNT"))) {

  forPlots<-temp2[,names(temp2) %in% c(var,"OS","response")]
  
  #=== SUPPLEMENTARY FIGURE 3A
  
  print(ggboxplot(forPlots,x = "response", y = var, 
                  color = "response", palette = c("slategrey","tomato"), 
                  add = "jitter", size=1, order = c("NR","R"),
                  ylab = "log Accumulative CDR3", xlab = "Response",title = var)+ 
          stat_compare_means(label="p.format", hide.ns = TRUE, size = 10,label.y = max(forPlots[,var],na.rm=T)+1, label.x = 1.25) + 
          theme(legend.position="none") +
          expand_limits(x = 0, y = 0) +
          myplot_parameters() +
          scale_y_continuous(limit = c(0, max(forPlots[,var],na.rm=T)+2)) )
  
  ggsave(paste0("Supp.3A.",var,".resp.vs.logclonecount.svg"),device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.6)
  
  #=== SUPPLEMENTARY FIGURE 3B
  
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
  
  ggsave(paste0("Supp.3B.",var,".OS.vs.logclonecount.svg"),device = "svg",width = 8,height = 8,dpi = "retina", scale = 0.6)
  
}
