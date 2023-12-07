

#====== CREATE DIVERSITY TABLE
#=== SHANNON ENTROPY
#=== FOR CDR3 AND FULL TCR (CDR3 + VGENE)

diversity.hla<-pheno[!is.na(pheno$response),]
combined.list.all<-combined.list.in[diversity.hla$ID]

diversity.hla$entropy<-sapply(combined.list.all, function (x) entropy(x$Read.count))
diversity.hla$vgene.entropy<-entropy.seq.internal(combined.list.all,HUMAN_TRBV,.quant = "read.count")


#=== SAVE DIVERSITY TABLE AS RDATA
save(diversity.hla,file = "diversity.hla.RData")




