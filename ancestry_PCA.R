#!/usr/bin/env Rscript

load("../result/icogs_flashpca.RData")

pheno_icogs=read.table("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T,sep="\t")
colnames(pheno_icogs)[which(colnames(pheno_icogs)=="SG_ID")]="ID"
tmp=as.data.frame(icogspca$projection)
colnames(tmp)=paste0("PC",1:10)
tmp1=unlist(strsplit(rownames(tmp),":"))
tmp1=unlist(strsplit(tmp1,"_"))
tmp$ID=tmp1[seq(2,length(tmp1),3)]

# tmp=read.table("../result/imp_icogs/merged1.eigenvec")
# colnames(tmp)=c("ID",paste0("PC",1:20))
# tmp1=unlist(strsplit(tmp$ID,"_"))
# tmp$ID=tmp1[seq(1,length(tmp1),2)]
pheno_icogs=pheno_icogs[pheno_icogs$ID %in% tmp$ID,]
idx=match(pheno_icogs$ID,tmp$ID)
pheno_icogs=cbind(pheno_icogs,tmp[idx,1:10])

load("../result/onco_flashpca.RData")
pheno_onco=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
colnames(pheno_onco)[which(colnames(pheno_onco)=="Onc_ID")]="ID"
tmp=as.data.frame(oncopca$projection)
colnames(tmp)=paste0("PC",1:10)
tmp1=unlist(strsplit(rownames(tmp),":"))
tmp1=unlist(strsplit(tmp1,"_"))
tmp$ID=tmp1[seq(2,length(tmp1),3)]
pheno_onco=pheno_onco[pheno_onco$ID %in% tmp$ID,]
idx=match(pheno_onco$ID,tmp$ID)
pheno_onco=cbind(pheno_onco,tmp[idx,1:10])

#colnames(pheno_onco)[5:19]=paste0("PC",1:15)
table(pheno_icogs$EthnicityGeno)
# African    Asian European 
# 1758    10822    72653
table(pheno_onco$EthnicityGeno)
# African    Asian European    other 
# 5569    25214   118419     2413 

source("/data/BB_Bioinformatics/Kevin/PRS_EASLC/code/theme_publication.R")
library(ggplot2)
plot_pca=function(dat=pheno_icogs,title="iCOGS")
{
  # pop="European"
  # dat1=dat[dat$EthnicityGeno==pop & is.na(dat$Behaviour1),]
  # p1=ggplot(dat1, aes(x=PC1, y=PC2)) + geom_point() + 
  #   labs(title=pop) + theme_Publication()
  # pop="Asian"
  # dat1=dat[dat$EthnicityGeno==pop & is.na(dat$Behaviour1),]
  # p2=ggplot(dat1, aes(x=PC1, y=PC2)) + geom_point() + 
  #   labs(title=pop) + theme_Publication()
  # pop="African"
  # dat1=dat[dat$EthnicityGeno==pop & is.na(dat$Behaviour1),]
  # p3=ggplot(dat1, aes(x=PC1, y=PC2)) + geom_point() + 
  #   labs(title=pop) + theme_Publication()
  # library(patchwork)
  # p1 + p2 + p3 
  dat1=dat[dat$EthnicityGeno %in% c("European","Asian","African") & is.na(dat$Behaviour1),]
  if (title=="iCOGS")
  {
    idx=which(dat1$EthnicityGeno=="African" & dat1$PC2>-0.2)
  }
  ggplot(dat1, aes(x=PC1, y=PC2,color=EthnicityGeno)) + geom_point() + 
     theme_Publication()+ labs(title=title,color='Ancestry') 
  ggsave(paste0("../result/PCAplot_",title,".pdf"),width = 8, height = 8)
}
