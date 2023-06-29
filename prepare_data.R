#!/usr/bin/env Rscript

library(data.table)
plink="/usr/local/apps/plink/1.9/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

pvar_icogs=fread("../result/imp_icogs/euro.pvar")
pvar_icogs$ID1=paste0(pvar_icogs$`#CHROM`,":",pvar_icogs$POS)
sum(duplicated(pvar_icogs$ID1))
# [1] 4211
pvar_onco=fread("../result/imp_onco/euro.pvar")
pvar_onco$ID1=paste0(pvar_onco$`#CHROM`,":",pvar_onco$POS)
#no need to align snps
table(pvar_icogs$ID %in% pvar_onco$ID)
# FALSE    TRUE 
# 170074 8845712
table(pvar_icogs$ID1 %in% pvar_onco$ID1)
# FALSE    TRUE 
# 169994 8845792 

#prepare a small dataset
options(scipen=99)
snppos=data.frame(chr=22,start=10000000,end=14000000,name="test")
write.table(snppos,file="../result/test.range",row.names = F,col.names = F,sep=" ",quote=F)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro --extract range ../result/test.range --recode A-transpose --out ../result/onco_test")
system(cmd)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro --extract range ../result/test.range --recode A-transpose --out ../result/icogs_test")
system(cmd)
