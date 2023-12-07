
#====== FIGURES 2A TO 2E

#=== INPUT FILE (IF REQUIRED)
#=== SIMILARITY HLA TABLE CREATED IN RSCRIPT 2
load("../Set_WD/similarity.hla.RData")


#=== FIGURE 2A

similarity.hla.subset<-similarity.hla[similarity.hla$method=="aa.Without.Vgene" & similarity.hla$Perc.clones==1,]
similarity.hla.subset<-similarity.hla.subset[order(similarity.hla.subset$OS,decreasing = F),]

ggviolin(similarity.hla.subset,x = "OS", y = "Similarity",
         color = "Resp", palette = c("slategrey","tomato"), add = "boxplot", fill = "Resp",
         add.params = list(fill = "white"),
         ylab = "CDR3 Similarity", xlab = "Overall Survival (weeks)")+
  #title = "Overall Survival and Similarity" ) +
  facet_grid(. ~ Resp, scales = "free", space = "free") +
  expand_limits(x = 0, y = 0) +
  theme(strip.text.x = element_text(size = 30), legend.position="none") +
  theme(text=element_text(family="Arial")) +
  myplot_parameters() +
  font("x.text",angle = 90,size = 18) +
  font("y.text",size = 28) 

ggsave("2A.survival.similarity.byResp.full.tcr.svg",device = "svg",width = 15,height = 15,dpi = "retina", scale = 0.8)

#=== FIGURE 2B

similarity.hla.subset<-similarity.hla[similarity.hla$method=="aa.Without.Vgene" & similarity.hla$Perc.clones==1,]
similarity.hla.subset<-similarity.hla.subset[order(similarity.hla.subset$OS,decreasing = F),]

ggboxplot(similarity.hla.subset,x = "Resp", y = "Similarity", 
          color = "Resp", palette = c("slategrey","tomato"), 
          add = "jitter", order = c("NR","R"), size=1,
          ylab = "CDR3 Similarity", xlab = "") +
  #,title = "Overall Survival and Similarity" ) + 
  stat_compare_means(size = 15, label.y = 0.25, label.x = 1.25, label = "p.signif") + 
  expand_limits(x = 0, y = 0) +
  theme(legend.position="none") +
  myplot_parameters() 

ggsave("2B.survival.similarity.byResp.cdr3.svg",device = "svg",width = 12,height = 15,dpi = "retina", scale = 0.8)

#=== FIGURE 2C

ggboxplot(similarity.hla.subset,x = "Resp", y = "Similarity", 
          color = "Resp", palette = c("slategrey","tomato"), 
          add = "jitter", order = c("NR","R"), size=1,
          ylab = "Full TCR Similarity", xlab = "") +
  #,title = "Overall Survival and Similarity" ) + 
  stat_compare_means(size = 15, label.y = 0.25, label.x = 1.25, label = "p.signif") + 
  expand_limits(x = 0, y = 0) +
  theme(legend.position="none") +
  myplot_parameters() 

ggsave("2C.survival.similarity.byResp.full.tcr.svg",device = "svg",width = 12,height = 15,dpi = "retina", scale = 0.8)

#=== FIGURE 2D

similarity.hla.subset<-similarity.hla[similarity.hla$method=="aa.Without.Vgene" & similarity.hla$Perc.clones %in% c(0.001,0.01,0.1),]
similarity.hla.subset<-similarity.hla.subset[order(similarity.hla.subset$Perc.clones,decreasing = T),]

ggboxplot(similarity.hla.subset, x = "Perc.clones", y = "Similarity", 
          color = "Resp", palette = c("slategrey","tomato"), 
          add = "jitter", size=1,
          ylab = "CDR3 Similarity", xlab = "Top% clones")+
  #,title = "Top% clones and Similarity" )  + 
  scale_x_discrete(breaks=c("0.1","0.01","0.001"), labels=c("10%", "1%", "0.1%")) +
  stat_compare_means(aes(group=Resp), size = 15, label.y = 0.6, label = "p.signif") + 
  myplot_parameters()+
  theme(legend.title=element_blank()) +
  font("legend.text", size = 30) 

ggsave("2D.perclones.similarity.byResp.cdr3.svg",device = "svg",width = 15,height = 15,dpi = "retina", scale = 0.8)

#=== FIGURE 2E

similarity.hla.subset<-similarity.hla[similarity.hla$method=="aa.Without.Vgene" & similarity.hla$Perc.clones==1 & !is.na(similarity.hla$Similarity),]
similarity.hla.subset<-similarity.hla.subset[order(similarity.hla.subset$HLA.count,similarity.hla.subset$OS,decreasing = F),]

similarity.hla.subset.v1<-as.data.frame(table(similarity.hla.subset$OS,similarity.hla.subset$HLA.count))
colnames(similarity.hla.subset.v1)<-c("OS","HLA.count","HLA.count.freq")

resp.val<-unique(similarity.hla.subset[,c("OS","Resp")])

similarity.hla.subset.freq<-merge(similarity.hla.subset.v1,resp.val,by.x="OS",all.x=T)
similarity.hla.subset.freq<-similarity.hla.subset.freq[similarity.hla.subset.freq$HLA.count.freq!=0,]

ggscatter(similarity.hla.subset.freq,x = "OS", y = "HLA.count", 
          color = "Resp", palette = c("slategrey","tomato"), size = "HLA.count.freq",
          ylab = "Number of HLA shared", xlab = "Overall Survival (weeks)") +
  #, title = "Overall Survival (weeks) and Number of HLA shared") +
  facet_grid(. ~ Resp, scales = "free", space = "free") +
  myplot_parameters() +
  font("x.text",angle = 90,size=18) +
  theme(strip.text.x = element_text(size = 28), legend.position="none") 


ggsave("2E.survival.hlaShared.byResp.svg",device = "svg",width = 15,height = 10,dpi = "retina", scale = 0.8)


#=== FIGURE 2F

similarity.hla.subset<-similarity.hla[similarity.hla$method=="aa.Without.Vgene" & similarity.hla$Perc.clones==1 & !is.na(similarity.hla$Similarity),]
similarity.hla.subset<-subset(similarity.hla.subset,select = -OS)
similarity.hla.subset.v1<-similarity.hla.subset[!duplicated(t(apply(similarity.hla.subset, 1, sort))), ]

similarity.hla.subset.v1<-similarity.hla.subset.v1[order(similarity.hla.subset.v1$HLA.count,decreasing = F),]

ggboxplot(similarity.hla.subset.v1, x = "HLA.count", y = "Similarity", 
          color = "Resp", palette = c("slategrey","tomato"), 
          add = "jitter", size=1,
          ylab = "CDR3 Similarity", xlab = "Number of HLA shared") +
  #title = "Number of HLA shared and Similarity" )  +
  stat_compare_means(aes(group=Resp), size = 12, label.y = 0.25, label = "p.signif") +
  expand_limits(x = 0, y = 0) +
  theme(legend.position="none")  +
  myplot_parameters()

ggsave("2F.hla.shared.similarity.byResp.all.clones.cdr3.svg", device = "svg",width = 15,height = 15,dpi = "retina", scale = 0.8)
