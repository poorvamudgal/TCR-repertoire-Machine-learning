
#====== SUPPLEMENTARY FIGURES 2A TO 2E

#=== INPUT FILE (IF REQUIRED)
#=== DIVERSITY HLA TABLE CREATED IN SCRIPT 3
load("../Set_WD/diversity.hla.RData")


#=== SETTING PARAMETERS

myplot_parameters <-function() { font("title", size = 30) +
    font("xlab", size = 35) +
    font("ylab", size = 35) +
    font("x.text",size = 22) +
    font("y.text",size = 35) +
    font("legend.title", size = 25) + 
    font("legend.text", size = 30) }

#=== SUPP FIGURE 2A

counts<-table(diversity.hla$HLA.A)
diversity.hla.subset<-diversity.hla[diversity.hla$HLA.A %in% names(counts)[counts>1],]
diversity.hla.subset$facetgroup<- ifelse(diversity.hla.subset$HLA.A %in% c("A*02;A*11","A*02;A*02","A*02;A*24"),1,2)

p1<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==1), x = "HLA.A", y = "vgene.entropy", color="response",
               ylab="Full TCR  Diversity", 
               #title="HLA-A diversity",
               add = "jitter", palette = c("slategrey","tomato"), size = 1,
               ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank()) 

p2<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==2), x = "HLA.A", y = "vgene.entropy", color="response",
               xlab="HLA-A type", ylab="Full TCR  Diversity", size = 1,
               add = "jitter", palette = c("slategrey","tomato"), ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank(), axis.title.y = element_text(vjust=2, hjust=-0.5))

ggarrange(p1 + rremove(c("xylab")),p2,ncol=1,nrow=2,align = "v", common.legend = TRUE) 

ggsave("Supp.2A.diversity.hlaA.full.tcr.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.6)

#=== SUPP FIGURE 2C

counts<-table(diversity.hla$HLA.B)
diversity.hla.subset<-diversity.hla[diversity.hla$HLA.B %in% names(counts)[counts>1],]
table(diversity.hla.subset$HLA.B)
diversity.hla.subset$facetgroup<- ifelse(diversity.hla.subset$HLA.B %in% c("B*15;B*46","B*38;B*46","B*40;B*46"),1,2)

p1<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==1), x = "HLA.B", y = "vgene.entropy", color="response",
               ylab="Full TCR  Diversity", 
               #title="HLA-B diversity",
               add = "jitter", palette = c("slategrey","tomato"), size = 1,
               ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank()) +
  font("x.text",size = 20)

p2<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==2), x = "HLA.B", y = "vgene.entropy", color="response",
               xlab="HLA-B type", ylab="Full TCR  Diversity", size = 1,
               add = "jitter", palette = c("slategrey","tomato"), ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank(), axis.title.y = element_text(vjust=2, hjust=-0.5)) +
  font("x.text",size = 17)

ggarrange(p1 + rremove(c("xylab")),p2,ncol=1,nrow=2,align = "v", common.legend = TRUE)

ggsave("Supp.2C.diversity.hlaB.full.tcr.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.6)


#=== SUPP FIGURE 2E

counts<-table(diversity.hla$HLA.C)
diversity.hla.subset<-diversity.hla[diversity.hla$HLA.C %in% names(counts)[counts>1],]
table(diversity.hla.subset$HLA.C)
diversity.hla.subset$facetgroup<- ifelse(diversity.hla.subset$HLA.C %in% c("C*01;C*07","Cw*01","Cw*01;Cw*03"),1,2)

p1<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==1), x = "HLA.C", y = "vgene.entropy", color="response",
               ylab="Full TCR  Diversity", 
               #title="HLA-C diversity",
               add = "jitter", palette = c("slategrey","tomato"), size = 1,
               ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank()) +
  font("x.text",size = 20)

p2<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==2), x = "HLA.C", y = "vgene.entropy", color="response",
               xlab="HLA-C type", ylab="Full TCR  Diversity", size = 1,
               add = "jitter", palette = c("slategrey","tomato"), ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank(), axis.title.y = element_text(vjust=2, hjust=-0.5)) +
  font("x.text",size = 17)

ggarrange(p1 + rremove(c("xylab")),p2,ncol=1,nrow=2,align = "v", common.legend = TRUE)

ggsave("Supp.2E.diversity.hlaC.full.tcr.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.6)


#=== SUPP FIGURE 2B

del1<-select(diversity.hla,-c("HLA.B2"))
del2<-select(diversity.hla,-c("HLA.B1"))
colnames(del2)<-colnames(del1)

diversity.hla.subset<-rbind(del1,del2)
diversity.hla.subset<-unique(diversity.hla.subset[!is.na(diversity.hla.subset$HLA.B1),])
diversity.hla.subset$facetgroup<- ifelse(diversity.hla.subset$HLA.B1 %in% c("B*46","B*40","B*51","B*38","B*15"),1,2)

p1<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==1), x = "HLA.B1", y = "vgene.entropy", color="response",
               ylab="Full TCR  Diversity", 
               #title="HLA-B alleles diversity",
               add = "jitter", palette = c("slategrey","tomato"),size = 1,
               ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank())

p2<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==2), x = "HLA.B1", y = "vgene.entropy", color="response",
               xlab="HLA-B alleles", ylab="Full TCR  Diversity", size = 1,
               add = "jitter", palette = c("slategrey","tomato"), ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank(), axis.title.y = element_text(vjust=2, hjust=-0.5))

ggarrange(p1 + rremove(c("xylab")),p2,ncol=1,nrow=2,align = "v", common.legend = TRUE)

ggsave("Supp.2B.diversity.hlaB.alleles.full.tcr.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.6)


#=== SUPP FIGURE 2D

del1<-select(diversity.hla,-c("HLA.C2"))
del2<-select(diversity.hla,-c("HLA.C1"))
colnames(del2)<-colnames(del1)

diversity.hla.subset<-rbind(del1,del2)
diversity.hla.subset<-unique(diversity.hla.subset[!is.na(diversity.hla.subset$HLA.C1),])
diversity.hla.subset$facetgroup<- ifelse(diversity.hla.subset$HLA.C1 %in% c("Cw*01","Cw*03","Cw*04","Cw*07","Cw*08"),1,2)

p1<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==1), x = "HLA.C1", y = "vgene.entropy", color="response",
               ylab="Full TCR  Diversity", 
               #title="HLA-C alleles diversity",
               add = "jitter", palette = c("slategrey","tomato"), size = 1,
               ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank()) +
  font("x.text",size = 20)

p2<- ggboxplot(subset(diversity.hla.subset,diversity.hla.subset$facetgroup==2), x = "HLA.C1", y = "vgene.entropy", color="response",
               xlab="HLA-C alleles", ylab="Full TCR  Diversity", size = 1,
               add = "jitter", palette = c("slategrey","tomato"), ylim = c(0,max(diversity.hla.subset$vgene.entropy))) +
  stat_compare_means(aes(group=response),label="p.signif", hide.ns = TRUE, label.y=10, size = 7) +
  myplot_parameters() +
  theme(legend.title=element_blank(), axis.title.y = element_text(vjust=2, hjust=-0.5)) +
  font("x.text",size = 16)

ggarrange(p1 + rremove(c("xylab")),p2,ncol=1,nrow=2,align = "v", common.legend = TRUE)

ggsave("Supp.2D.diversity.hlaC.alleles.full.tcr.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.6)
