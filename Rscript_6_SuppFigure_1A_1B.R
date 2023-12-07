
#====== SUPPLEMENTARY FIGURES 1A AND 1B

#=== INPUT FILE (IF REQUIRED)
#=== DIVERSITY HLA TABLE CREATED IN SCRIPT 3
load("../Set_WD/diversity.hla.RData")

#=== SUPP FIGURE 1A

ggscatter(diversity.hla, x = "OS", y = "entropy", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 10,
          add.params = list(fill = "gray95"),
          xlab = "Survival (weeks)", ylab = "CDR3 Diversity") +
  expand_limits(x = 0, y = 0) +
  scale_y_continuous(limit = c(0, max(diversity.hla[,"entropy"],na.rm=T)+1.5)) +
  myplot_parameters()

ggsave("Supp.1A.survival.diversity.cdr3.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.8)

#=== SUPP FIGURE 1B

ggscatter(diversity.hla, x = "OS", y = "vgene.entropy", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 10,
          add.params = list(fill = "gray95"),
          xlab = "Survival (weeks)", ylab = "Full TCR Diversity") +
  expand_limits(x = 0, y = 0) +
  scale_y_continuous(limit = c(0, max(diversity.hla[,"vgene.entropy"],na.rm=T)+1.5)) +
  myplot_parameters()

ggsave("Supp.1B.survival.diversity.full.tcr.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.8)

#=== SUPP 

ggboxplot(diversity.hla, x = "Sex", y = "entropy",
          color = "Sex", palette = c("slategrey","tomato"), 
          add = "jitter",  order = c("M","F"), size=1,
          ylab = "CDR3  Diversity", xlab="") +
  #, title = "Diversity") + 
  stat_compare_means(label.y = 13, label.x = 1.25, size = 12, label = "p.format") +
  theme(legend.position="none") +
  expand_limits(x = 0, y = 0) +
  myplot_parameters()

ggsave("Supp.gender.diversity.cdr3.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.8)

#=== Supp

ggscatter(diversity.hla, x = "Age", y = "entropy", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 10,
          add.params = list(fill = "gray95"),
          xlab = "Age", ylab = "CDR3 Diversity") +
  expand_limits(x = 0, y = 0) +
  scale_y_continuous(limit = c(0, max(diversity.hla[,"entropy"],na.rm=T)+1.5)) +
  myplot_parameters()

ggsave("Supp.age.diversity.cdr3.svg",device = "svg",width = 10,height = 12,dpi = "retina", scale = 0.8)