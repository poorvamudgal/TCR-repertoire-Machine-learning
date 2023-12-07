
#====== FIGURES 1A TO 1C

#=== INPUT FILE (IF REQUIRED)
#=== DIVERSITY HLA TABLE CREATED IN SCRIPT 3
#load("../Set_WD/diversity.hla.RData")

#=== FIGURE 1A

ggboxplot(diversity.hla, x = "response", y = "entropy",
          color = "response", palette = c("slategrey","tomato"), 
          add = "jitter",  order = c("NR","R"), size=1,
          ylab = "CDR3  Diversity", xlab="") +
  #, title = "Diversity") + 
  stat_compare_means(label.y = 13, label.x = 1.25, size = 12, label = "p.signif") +
  theme(legend.position="none") +
  expand_limits(x = 0, y = 0) +
  myplot_parameters()

#ggsave("1A.diversity.cdr3.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.8)

#=== FIGURE 1B

ggboxplot(diversity.hla, x = "response", y = "vgene.entropy",
          color = "response", palette = c("slategrey","tomato"), size=1,
          add = "jitter",  order = c("NR","R"),
          ylab = "Full TCR Diversity", xlab="" )+
  #, title = "Diversity") + 
  stat_compare_means(label.y = 4.5, label.x = 1.25, size = 12, label = "p.signif") +
  theme(legend.position="none") +
  expand_limits(x = 0, y = 0) +
  myplot_parameters()


#ggsave("1B.diversity.full.tcr.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.8)


#=== FIGURE 1C

del1<-select(diversity.hla,-c("HLA.A2"))
del2<-select(diversity.hla,-c("HLA.A1"))
colnames(del2)<-colnames(del1)

diversity.hla.subset<-rbind(del1,del2)
diversity.hla.subset<-unique(diversity.hla.subset[!is.na(diversity.hla.subset$HLA.A1),])
diversity.hla.subset$facetgroup<- ifelse(diversity.hla.subset$HLA.A1 %in% c("A*02","A*11","A*33","A*24"),1,2)

ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==1), 
          x = "HLA.A1", y = "vgene.entropy", 
          color="response", size=1,
          ylab="Full TCR Diversity", xlab = "HLA-A alleles",
          add = "jitter", palette = c("slategrey","tomato"), 
          ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 12) +
  myplot_parameters() +
  theme(legend.title=element_blank()) +
  font("legend.text", size = 30)

#ggsave("1C.diversity.hlaA.alleles.full.tcr.svg",device = "svg",width = 15,height = 12,dpi = "retina", scale = 0.8)

