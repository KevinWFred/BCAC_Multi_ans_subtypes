#!/usr/bin/env Rscript
# to select SNPs for EUR/ASN separately. The SNPs are in icogs and onco; having MAF>0.01 in icogs/onco.
#create genotype data for tuning/validation, snpid use rsid.
#use snps from onco

plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
library(data.table)
library(dplyr)

#QCed SNP lists:
qc_african_onco=as.data.frame(fread("../result/imp_QC/onco/african_info.snplist",header=F))
qc_asian_onco=as.data.frame(fread("../result/imp_QC/onco/asian_info.snplist",header=F))
qc_euro_onco=as.data.frame(fread("../result/imp_QC/onco/euro_info.snplist",header=F))
qc_hispanic_onco=as.data.frame(fread("../result/imp_QC/onco/hispanic_info.snplist",header=F))
qc_african_icogs=as.data.frame(fread("../result/imp_QC/icogs/african_info.snplist",header=F))
qc_asian_icogs=as.data.frame(fread("../result/imp_QC/icogs/asian_info.snplist",header=F))
qc_euro_icogs=as.data.frame(fread("../result/imp_QC/icogs/euro_info.snplist",header=F))

#to process EUR sumstats
#generate genotype data for tunning and validation (EUR and Asian separately)
#Freq (MAF >0.01 in icogs/onco EUR EAS)
eur_tr_icogs=read.table("../result/PRS_icogs_trainingpheno.txt",header=T,sep="\t")
eur_tr_icogs=eur_tr_icogs[eur_tr_icogs$EthnicityGeno=="European",]
tmp=data.frame(FID=eur_tr_icogs$SG_ID,IID=eur_tr_icogs$SG_ID)
write.table(tmp,file="../result/euro_training_icogs.samples",row.names = F,col.names = F,sep="\t",quote=F)

eur_tr_onco=read.table("../result/PRS_onco_trainingpheno.txt",header=T,sep="\t")
eur_tr_onco=eur_tr_onco[eur_tr_onco$EthnicityGeno=="European",]
tmp=data.frame(FID=eur_tr_onco$Onc_ID,IID=eur_tr_onco$Onc_ID)
write.table(tmp,file="../result/euro_training_onco.samples",row.names = F,col.names = F,sep="\t",quote=F)

#/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian.afreq
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro --keep ../result/euro_training_icogs.samples --freq --out /data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_icogs --memory 100000 --threads 8 ")
system(cmd)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro --keep ../result/euro_training_onco.samples --freq --out /data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_onco --memory 100000 --threads 8 ")
system(cmd)


#Missingness
#/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro.vmiss
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro --keep ../result/euro_training_onco.samples --missing variant-only --genotyping-rate dosage --out /data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_onco --memory 100000 --threads 8 ")
system(cmd)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro --keep ../result/euro_training_icogs.samples --missing variant-only --genotyping-rate dosage --out /data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_icogs --memory 100000 --threads 8 ")
system(cmd)


pvar_onco=fread("../result/imp_onco/euro/euro.pvar")
pvar_onco=pvar_onco[pvar_onco$ID %in% qc_euro_onco$V1,]
# pvar_icogs=fread("../result/imp_icogs/euro/euro.pvar")
# pvar_icogs=pvar_icogs[pvar_icogs$ID %in% qc_euro_icogs$V1,]
# eurcomsnps=intersect(pvar_icogs$ID,pvar_onco$ID)
# length(eurcomsnps)
# [1] 7803064
eurmetafile="../result/meta_euro_lr.meta.rsid"
eursumdat=as.data.frame(fread(eurmetafile))
# table(eurcomsnps %in% eursumdat$SNP[eursumdat$N==2])
# table(eursumdat$SNP[eursumdat$N==2] %in% eurcomsnps)
# selsnps=eursumdat$SNP[!is.na(eursumdat$rsid) & eursumdat$N==2]
selsnps=eursumdat$SNP[eursumdat$SNP %in% pvar_onco$ID]
length(selsnps) #9928681
eur_icogs_freq=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_icogs.afreq"))
eur_onco_freq=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_onco.afreq"))
tmp=eur_icogs_freq %>% filter(ALT_FREQS<0.99 & ALT_FREQS>0.01) %>% select(ID)
# tmp2=eur_onco_freq %>% filter(ALT_FREQS<0.99 & ALT_FREQS>0.01) %>% select(ID)
# tmp=unique(c(tmp$ID,tmp2$ID))
selsnps=intersect(selsnps,tmp$ID)
length(selsnps) #8551464
#eur_icogs_missing=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_icogs.vmiss"))
eur_onco_missing=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_onco.vmiss"))
# idx=order(eur_icogs_missing$F_MISS,decreasing = T)
# tmp=eur_icogs_missing$ID[idx][1:10]
# write.table(tmp,file="../result/test.samples",row.names = F,col.names = F,quote=F)
# cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro --keep ../result/euro_training_icogs.samples --extract ../result/test.samples --recode A-transpose --out ../result/test")
# system(cmd)
# tmp=as.data.frame(fread("../result/test.traw"))
# sum(is.na(tmp)) #0 dosage data have no missing
# cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro --keep ../result/euro_training_icogs.samples --extract ../result/test.samples --make-bed --out ../result/test")
# system(cmd)
# cmd=paste0(plink," --bfile ../result/test --recode A-transpose --out ../result/test")
# system(cmd)
# tmp=as.data.frame(fread("../result/test.traw"))
# tmp1=eur_icogs_missing %>% filter(F_MISS<0.05) %>% select(ID)
# tmp2=eur_onco_missing %>% filter(F_MISS<0.05) %>% select(ID)
# tmp=unique(c(tmp1$ID,tmp2$ID))
# selsnps=intersect(selsnps,tmp2$ID)
# length(selsnps)
# tmp1=eur_icogs_missing %>% filter(F_MISS<0.1) %>% select(ID)
#tmp11=eur_icogs_missing %>% filter(F_MISS>0.9) %>% select(ID)
# tmp2=eur_onco_missing %>% filter(F_MISS<0.1) %>% select(ID)
tmp22=eur_onco_missing %>% filter(F_MISS>0.5) %>% select(ID)#snps with high missing rate in tuning
# tmp=unique(c(tmp1$ID,tmp2$ID))
# tmp=tmp[!tmp %in% tmp22]
# selsnps=intersect(selsnps,tmp)
selsnps=selsnps[!selsnps %in% c(tmp22$ID)]
length(selsnps) #8548410
idx=match(selsnps,eursumdat$SNP)
eursumdat1=eursumdat[idx,]
eur_neff=as.integer(4*(1/(1/33505+1/33961)+1/(1/(26699+1/35150)))) #for LDpred2
#most SNPid are chr:pos:alt:ref; some are chr:pos:ref:alt. As long as use rsid and A1, no problem.
formsumstats=function(sumdat=eursumdat1,outprefix="euro_training",neff=eur_neff)
{
  idx=which(colnames(sumdat) %in% c("chromosome","CHR"))
  if (length(idx)>0) colnames(sumdat)[idx]="chr"
  idx=which(colnames(sumdat) %in% c("position","BP"))
  if (length(idx)>0) colnames(sumdat)[idx]="pos"
  idx=which(colnames(sumdat) %in% c("alleleB","Effect.allele","Allele1","A1"))
  if (length(idx)>0) colnames(sumdat)[idx]="a1" #effect
  idx=which(colnames(sumdat) %in% c("alleleA","Reference.allele","Allele2","A2"))
  if (length(idx)>0) colnames(sumdat)[idx]="a0"
  
  sumdat$n_eff=neff #For LDpred2
  sumdat$beta=log(sumdat$OR)
  #sumdat$zscore <- qnorm(1 - sumdat$P/2)
  sumdat$zscore=abs(qnorm(sumdat$P/2))
  sumdat$beta_se=abs(sumdat$beta/sumdat$zscore)
  
  colnames(sumdat)[which(colnames(sumdat)=="P")]="p"
  colnames(sumdat)[which(colnames(sumdat)=="P(R)")]="P_R"
  sumdat=sumdat[!is.na(sumdat$rsid),]
  # LDpred 2 require the header to follow the exact naming. a1 :effect allelle
  #sumstats=data.frame(chr=sumdat$chr,pos=sumdat$pos,rsid=sumdat$rsid,a0=toupper(sumdat$a0),a1=toupper(sumdat$a1),n_eff=sumdat$n_eff,beta_se=sumdat$beta_se,p=sumdat$p,beta=sumdat$beta)
  
  fwrite(sumdat,file=paste0("../result/PRS1/",outprefix,"_sumstats.txt"),row.names = F,sep="\t",quote=F)
  #return(sumstats)
}
formsumstats()

#to process asian sumstats
asn_tr_icogs=read.table("../result/PRS_icogs_trainingpheno.txt",header=T,sep="\t")
asn_tr_icogs=asn_tr_icogs[asn_tr_icogs$EthnicityGeno=="Asian",]
tmp=data.frame(FID=asn_tr_icogs$SG_ID,IID=asn_tr_icogs$SG_ID)
write.table(tmp,file="../result/asian_training_icogs.samples",row.names = F,col.names = F,sep="\t",quote=F)

asn_tr_onco=read.table("../result/PRS_onco_trainingpheno.txt",header=T,sep="\t")
asn_tr_onco=asn_tr_onco[asn_tr_onco$EthnicityGeno=="Asian",]
tmp=data.frame(FID=asn_tr_onco$Onc_ID,IID=asn_tr_onco$Onc_ID)
write.table(tmp,file="../result/asian_training_onco.samples",row.names = F,col.names = F,sep="\t",quote=F)

#/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian.afreq
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian --keep ../result/asian_training_icogs.samples --freq --out /data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_icogs --memory 100000 --threads 8 ")
system(cmd)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian --keep ../result/asian_training_onco.samples --freq --out /data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_onco --memory 100000 --threads 8 ")
system(cmd)

#Missingness
#/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian.vmiss
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian --keep ../result/asian_training_onco.samples --missing variant-only --genotyping-rate dosage --out /data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_onco --memory 100000 --threads 8 ")
system(cmd)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian --keep ../result/asian_training_icogs.samples --missing variant-only --genotyping-rate dosage --out /data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_icogs --memory 100000 --threads 8 ")
system(cmd)

pvar_onco_asian=fread("../result/imp_onco/asian/asian.pvar") #10083629
pvar_onco_asian=pvar_onco_asian[pvar_onco_asian$ID %in% qc_asian_onco$V1,] #8405662
# pvar_icogs=fread("../result/imp_icogs/asian/asian.pvar")
# pvar_icogs=pvar_icogs[pvar_icogs$ID %in% qc_asian_icogs$V1,]
# asncomsnps=intersect(pvar_icogs$ID,pvar_onco$ID)
# length(asncomsnps)
# # [1] 5343414
asnmetafile="../result/meta_asian_lr.meta.rsid"
asnsumdat=as.data.frame(fread(asnmetafile)) #9283968
# table(asncomsnps %in% asnsumdat$SNP[asnsumdat$N==2])
# table(asnsumdat$SNP[asnsumdat$N==2] %in% asncomsnps)
#selsnps=asnsumdat$SNP[!is.na(asnsumdat$rsid) & asnsumdat$N==2]
selsnps=asnsumdat$SNP[asnsumdat$SNP %in% pvar_onco_asian$ID]
length(selsnps) #8405662
#asn_icogs_freq=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_icogs.afreq"))
asn_onco_freq=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_onco.afreq"))
#tmp1=asn_icogs_freq %>% filter(ALT_FREQS<0.99 & ALT_FREQS>0.01) %>% select(ID)
tmp=asn_onco_freq %>% filter(ALT_FREQS<0.99 & ALT_FREQS>0.01) %>% select(ID)
#tmp=unique(c(tmp1$ID,tmp2$ID))
selsnps=intersect(selsnps,tmp$ID)
length(selsnps) #7040049

#asn_icogs_missing=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_icogs.vmiss"))
asn_onco_missing=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_onco.vmiss"))
#tmp1=asn_icogs_missing %>% filter(F_MISS<0.1) %>% select(ID)
#tmp11=asn_icogs_missing %>% filter(F_MISS>0.9) %>% select(ID)
#tmp2=asn_onco_missing %>% filter(F_MISS<0.1) %>% select(ID)
tmp22=asn_onco_missing %>% filter(F_MISS>0.5) %>% select(ID)
# tmp=unique(c(tmp1$ID,tmp2$ID))
# tmp=tmp[!tmp %in% tmp22]
# selsnps=intersect(selsnps,tmp)
selsnps=selsnps[!selsnps %in% c(tmp22$ID)]
length(selsnps) #7008849
idx=match(selsnps,asnsumdat$SNP)
asnsumdat1=asnsumdat[idx,]
asn_neff=as.integer(4*(1/(1/5940+1/4882)+1/(1/(6859+1/7545)))) #for LDpred2
formsumstats(sumdat=asnsumdat1,outprefix="asian_training",neff=asn_neff)

#generate tuning/validation data
eursumdat2=as.data.frame(fread("../result/PRS1/euro_training_sumstats.txt")) #8245578
asnsumdat2=as.data.frame(fread("../result/PRS1/asian_training_sumstats.txt")) #6756271
alltrainsnps=unique(c(eursumdat2$SNP,asnsumdat2$SNP)) 
length(alltrainsnps) #9719954
#to map SNP to rsid
tmp1=eursumdat2[,c("SNP","rsid")]
tmp1$type="euro"
tmp2=asnsumdat2[,c("SNP","rsid")]
tmp2$type="asian"
tmp3=merge(tmp1,tmp2,by="SNP")
table(tmp3$rsid.x==tmp3$rsid.y) #T
tmp2=tmp2[!tmp2$SNP %in% tmp1$SNP,]
tmp=rbind(tmp1,tmp2)
idx=which(duplicated(tmp$rsid))
idx=which(tmp$rsid %in% tmp$rsid[idx])
tmp1=tmp[idx,] #5674 SNPs have duplicated rsid between eur and asn
# idx1=match(tmp1$SNP,eursumdat1$SNP)
# idx2=match(tmp1$SNP,asnsumdat1$SNP)
# intersect(idx1[!is.na(idx1)],idx2[!is.na(idx2)])
# tmp1$p_euro=eursumdat1$p[idx1]
# tmp1$p_asian=asnsumdat1$p[idx2]
# tmp1$p=rowSums(tmp1[,c("p_euro","p_asian")],na.rm=T)
# tmp2=tmp1 %>% group_by(rsid) %>% arrange(p) %>% filter(row_number()>1) #SNPs to be removed
# eur_snp_2rm=tmp2$SNP[tmp2$type=="euro"]
# asn_snp_2rm=tmp2$SNP[tmp2$type=="asian"]
eur_snp_2rm=asn_snp_2rm=tmp1$SNP
eursumdat3=eursumdat2[!eursumdat2$SNP %in% eur_snp_2rm,]
nrow(eursumdat3) #8242741
fwrite(eursumdat3,file="../result/PRS1/euro_training_sumstats.txt",row.names = F,sep="\t",quote=F)
asnsumdat3=asnsumdat2[!asnsumdat2$SNP %in% asn_snp_2rm,]
nrow(asnsumdat3) #6753434
fwrite(asnsumdat3,file="../result/PRS1/asian_training_sumstats.txt",row.names = F,sep="\t",quote=F)

alltrainsnps=unique(c(eursumdat3$SNP,asnsumdat3$SNP))
length(alltrainsnps) #9714280
fwrite(data.frame(alltrainsnps),file="../result/PRS1/trainsnps.txt",row.names = F,col.names = F,quote=F)
# #to change to rsid
tmp1=eursumdat3[,c("SNP","rsid")]
tmp2=asnsumdat3[,c("SNP","rsid")]
tmp2=tmp2[!tmp2$SNP %in% tmp1$SNP,]
tmp=rbind(tmp1,tmp2)
sum(duplicated(tmp$rsid)) #0
fwrite(tmp,file="../result/PRS1/trainsnps_updatename.txt",row.names = F,col.names = F,quote=F,sep="\t")
# cmd=paste0(plink2," --pfile ../result/imputation/merged --update-name ../result/merged_updatename.txt --make-pgen --out ../result/imputation/merged_rsid")
# system(cmd)

formgenotype=function(phenofile="../result/PRS_onco_euro_tuningpheno.txt",
                      prefix="../result/imp_onco/euro/euro",
                      outprefix="../result/PRS1/euro_onco_tuning")
{
  pheno=read.table(phenofile,header=T,sep="\t")
  if (sum(colnames(pheno)=="Onc_ID")>0)
  {
    tmp=data.frame(FID=pheno$Onc_ID,IID=pheno$Onc_ID)
  }else
  {
    tmp=data.frame(FID=pheno$SG_ID,IID=pheno$SG_ID)
  }
  
  samplefile=paste0("../result/PRS1/",basename(outprefix),".samples")
  write.table(tmp,file=samplefile,row.names = F,col.names = F,sep="\t",quote=F)
  cmd=paste0(plink2," --pfile ",prefix," --keep ",samplefile," --extract ../result/PRS1/trainsnps.txt --make-pgen --out ",outprefix," --memory 100000 --threads 8 ")
  system(cmd)
  cmd=paste0(plink2," --pfile ",outprefix, " --update-name ../result/PRS1/trainsnps_updatename.txt --make-pgen --out ",outprefix," --memory 100000 --threads 8 ")
  system(cmd)
  cmd=paste0(plink2," --pfile ",outprefix, " --make-bed --out ",outprefix," --memory 100000 --threads 8 ")
  system(cmd)
}
formgenotype()
formgenotype(phenofile="../result/PRS_onco_asian_tuningpheno.txt",
                      prefix="../result/imp_onco/asian/asian",
                      outprefix="../result/PRS1/asian_onco_tuning")
formgenotype(phenofile="../result/PRS_icogs_african_validationpheno.txt",
             prefix="../result/imp_icogs/african/african",
             outprefix="../result/PRS1/african_icogs_validation")
formgenotype(phenofile="../result/PRS_onco_african_validationpheno.txt",
             prefix="../result/imp_onco/african/african",
             outprefix="../result/PRS1/african_onco_validation")
formgenotype(phenofile="../result/PRS_onco_asian_validationpheno.txt",
             prefix="../result/imp_onco/asian/asian",
             outprefix="../result/PRS1/asian_onco_validation")
formgenotype(phenofile="../result/PRS_onco_euro_validationpheno.txt",
             prefix="../result/imp_onco/euro/euro",
             outprefix="../result/PRS1/euro_onco_validation")
formgenotype(phenofile="../result/PRS_onco_other_validationpheno.txt",
             prefix="../result/imp_onco/hispanic/hispanic",
             outprefix="../result/PRS1/hispanic_onco_validation")

#check if they have the same allele coding in different datasets
compare2dat=function(prefix1="../result/PRS1/euro_onco_tuning",prefix2="../result/PRS1/asian_onco_tuning")
{
  dat1=fread(paste0(prefix1,".bim"))
  dat2=fread(paste0(prefix2,".bim"))
  dat1$ID=paste0(dat1$V1,":",dat1$V4,":",dat1$V5,":",dat1$V6)
  dat2$ID=paste0(dat2$V1,":",dat2$V4,":",dat2$V5,":",dat2$V6)
  dat1$ID1=paste0(dat1$V1,":",dat1$V4)
  dat2$ID1=paste0(dat2$V1,":",dat2$V4)
  idx1=dat1$ID1 %in% dat2$ID1
  idx2=dat2$ID1 %in% dat1$ID1
  dat1=dat1[idx1,]
  dat2=dat2[idx2,]
  print(paste0("dat1 in dat2: ",nrow(dat1)))
  print(paste0("dat2 in dat1: ",nrow(dat2)))
  idx1=dat1$ID %in% dat2$ID
  idx2=dat2$ID %in% dat1$ID
  if (sum(!idx1)>0)
  {
    dat11=dat1[!idx1,]
    dat22=dat2[!idx2,]
    tmp=unlist(strsplit(dat11$ID,":"))
    chr=tmp[seq(1,length(tmp),4)]
    pos=tmp[seq(2,length(tmp),4)]
    alt=tmp[seq(3,length(tmp),4)]#this is not right, can be ref
    ref=tmp[seq(4,length(tmp),4)]#this is not right, can be alt
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
}

compare2dat(prefix1="../result/PRS1/euro_onco_tuning",prefix2="../result/PRS1/asian_onco_tuning")
compare2dat(prefix1="../result/PRS1/euro_onco_tuning",prefix2="../result/PRS1/african_icogs_validation")
compare2dat(prefix1="../result/PRS1/euro_onco_tuning",prefix2="../result/PRS1/african_onco_validation")
compare2dat(prefix1="../result/PRS1/euro_onco_tuning",prefix2="../result/PRS1/asian_onco_validation")
compare2dat(prefix1="../result/PRS1/euro_onco_tuning",prefix2="../result/PRS1/euro_onco_validation")
compare2dat(prefix1="../result/PRS1/euro_onco_tuning",prefix2="../result/PRS1/hispanic_onco_validation")
compare2dat(prefix1="../result/PRS1/euro_onco_tuning",prefix2="../result/PRS1/tuning_validation")
