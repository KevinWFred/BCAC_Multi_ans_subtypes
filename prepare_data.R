#!/usr/bin/env Rscript

library(data.table)
plink="/usr/local/apps/plink/1.9/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

pvar_icogs=fread("../result/imp_icogs/euro/euro.pvar")
tmp=unlist(strsplit(pvar_icogs$ID,":"))
ALT=tmp[seq(3,length(tmp),4)]
all(ALT==pvar_icogs$ALT)
pvar_icogs$ID1=paste0(pvar_icogs$`#CHROM`,":",pvar_icogs$POS)
sum(duplicated(pvar_icogs$ID1))
# [1] 4211
pvar_onco=fread("../result/imp_onco/euro/euro.pvar")
tmp=unlist(strsplit(pvar_onco$ID,":"))
ALT=tmp[seq(3,length(tmp),4)]
all(ALT==pvar_onco$ALT)
pvar_onco$ID1=paste0(pvar_onco$`#CHROM`,":",pvar_onco$POS)
#no need to align snps
table(pvar_icogs$ID %in% pvar_onco$ID)
# FALSE    TRUE 
# 170074 8845712
table(pvar_icogs$ID1 %in% pvar_onco$ID1)
# FALSE    TRUE 
# 169994 8845792 
compare2dat=function(dat1=pvar_icogs,dat2=pvar_onco)
{
  dat1$ID1=paste0(dat1$`#CHROM`,":",dat1$POS)
  dat2$ID1=paste0(dat2$`#CHROM`,":",dat2$POS)
  idx1=dat1$ID1 %in% dat2$ID1
  idx2=dat2$ID1 %in% dat1$ID1
  dat1=dat1[idx1,]
  dat2=dat2[idx2,]
  print(paste0("dat1 in dat2: ",nrow(dat1)))
  print(paste0("dat2 in dat1: ",nrow(dat2)))
  idx1=dat1$ID %in% dat2$ID
  idx2=dat2$ID %in% dat1$ID
  dat11=dat1[!idx1,]
  dat22=dat2[!idx2,]
  tmp=unlist(strsplit(dat11$ID,":"))
  chr=tmp[seq(1,length(tmp),4)]
  pos=tmp[seq(2,length(tmp),4)]
  alt=tmp[seq(3,length(tmp),4)]
  ref=tmp[seq(4,length(tmp),4)]
  tmp1=paste0(chr,":",pos,":",ref,":",alt)
  tmp2=tmp1[tmp1 %in% dat2$ID]
  print(paste0("dat1 mismatch in dat2: ",length(tmp2))) #number of mismatched snps in dat1
  tmp=unlist(strsplit(dat22$ID,":"))
  chr=tmp[seq(1,length(tmp),4)]
  pos=tmp[seq(2,length(tmp),4)]
  alt=tmp[seq(3,length(tmp),4)]
  ref=tmp[seq(4,length(tmp),4)]
  tmp1=paste0(chr,":",pos,":",ref,":",alt)
  tmp2=tmp1[tmp1 %in% dat1$ID]
  print(paste0("dat2 mismatch in dat1: ",length(tmp2)))
}
compare2dat()
# [1] "dat1 in dat2: 10151030"
# [1] "dat2 in dat1: 10151155"
# [1] "dat1 mismatch in dat2: 43"
# [1] "dat2 mismatch in dat1: 56"
pvar_onco_asian=fread("../result/imp_onco/asian/asian.pvar")
pvar_onco_african=fread("../result/imp_onco/african/african.pvar")
pvar_onco_hispanic=fread("../result/imp_onco/hispanic/hispanic.pvar")
pvar_icogs_asian=fread("../result/imp_icogs/asian/asian.pvar")
pvar_icogs_african=fread("../result/imp_icogs/african/african.pvar")
compare2dat(dat1=pvar_onco_african,dat2=pvar_onco)
# [1] "dat1 in dat2: 8317609"
# [1] "dat2 in dat1: 8306629"
# [1] "dat1 mismatch in dat2: 3440"
# [1] "dat2 mismatch in dat1: 1635"
compare2dat(dat1=pvar_onco,dat2=pvar_onco_asian)
# [1] "dat1 in dat2: 7811728"
# [1] "dat2 in dat1: 7811514"
# [1] "dat1 mismatch in dat2: 732"
# [1] "dat2 mismatch in dat1: 715"
compare2dat(dat1=pvar_onco,dat2=pvar_icogs_asian)
# [1] "dat1 in dat2: 8217973"
# [1] "dat2 in dat1: 8218510"
# [1] "dat1 mismatch in dat2: 702"
# [1] "dat2 mismatch in dat1: 845"
compare2dat(dat1=pvar_onco,dat2=pvar_icogs_african)
# [1] "dat1 in dat2: 8280668"
# [1] "dat2 in dat1: 8291335"
# [1] "dat1 mismatch in dat2: 1633"
# [1] "dat2 mismatch in dat1: 3396"
compare2dat(dat1=pvar_onco,dat2=pvar_onco_hispanic)
# [1] "dat1 in dat2: 9174277"
# [1] "dat2 in dat1: 9175788"
# [1] "dat1 mismatch in dat2: 426"
# [1] "dat2 mismatch in dat1: 797"
#prepare a small dataset
options(scipen=99)
snppos=data.frame(chr=22,start=10000000,end=14000000,name="test")
write.table(snppos,file="../result/test.range",row.names = F,col.names = F,sep=" ",quote=F)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro --extract range ../result/test.range --recode A-transpose --out ../result/onco_test")
system(cmd)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro --extract range ../result/test.range --recode A-transpose --out ../result/icogs_test")
system(cmd)
