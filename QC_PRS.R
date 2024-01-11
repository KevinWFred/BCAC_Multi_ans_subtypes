#!/usr/bin/env Rscript
# to select SNPs for EUR/ASN separately. The SNPs are in icogs and onco; having MAF>0.01 in icogs/onco.
#create genotype data for tuning/validation, snpid use rsid.

plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
library(data.table)
library(dplyr)
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
pvar_icogs=fread("../result/imp_icogs/euro/euro.pvar")
eurcomsnps=intersect(pvar_icogs$ID,pvar_onco$ID)
length(eurcomsnps)
# [1] 10150740
eurmetafile="../result/meta_euro_lr.meta.rsid"
eursumdat=as.data.frame(fread(eurmetafile))
table(eurcomsnps %in% eursumdat$SNP[eursumdat$N==2])
table(eursumdat$SNP[eursumdat$N==2] %in% eurcomsnps)
selsnps=eursumdat$SNP[!is.na(eursumdat$rsid) & eursumdat$N==2]
length(selsnps) #9786199
eur_icogs_freq=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_icogs.afreq"))
eur_onco_freq=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_onco.afreq"))
tmp1=eur_icogs_freq %>% filter(ALT_FREQS<0.99 & ALT_FREQS>0.01) %>% select(ID)
tmp2=eur_onco_freq %>% filter(ALT_FREQS<0.99 & ALT_FREQS>0.01) %>% select(ID)
tmp=unique(c(tmp1$ID,tmp2$ID))
selsnps=intersect(selsnps,tmp)
length(selsnps) #8770007
eur_icogs_missing=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/euro_training_icogs.vmiss"))
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
# length(selsnps) #8153686
# tmp1=eur_icogs_missing %>% filter(F_MISS<0.1) %>% select(ID)
tmp11=eur_icogs_missing %>% filter(F_MISS>0.9) %>% select(ID)
# tmp2=eur_onco_missing %>% filter(F_MISS<0.1) %>% select(ID)
tmp22=eur_onco_missing %>% filter(F_MISS>0.5) %>% select(ID)#snps with high missing rate in tuning
# tmp=unique(c(tmp1$ID,tmp2$ID))
# tmp=tmp[!tmp %in% tmp22]
# selsnps=intersect(selsnps,tmp)
selsnps=selsnps[!selsnps %in% c(tmp22$ID,tmp11$ID)]
length(selsnps) #8708422
idx=match(selsnps,eursumdat$SNP)
eursumdat1=eursumdat[idx,]
eur_neff=as.integer(4*(1/(1/33505+1/33961)+1/(1/(26699+1/35150)))) #for LDpred2
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
  sumdat$zscore <- qnorm(1 - sumdat$P/2)
  sumdat$beta_se=abs(sumdat$beta/sumdat$zscore)
  
  colnames(sumdat)[which(colnames(sumdat)=="P")]="p"
  colnames(sumdat)[which(colnames(sumdat)=="P(R)")]="P_R"
  # LDpred 2 require the header to follow the exact naming. a1 :effect allelle
  #sumstats=data.frame(chr=sumdat$chr,pos=sumdat$pos,rsid=sumdat$rsid,a0=toupper(sumdat$a0),a1=toupper(sumdat$a1),n_eff=sumdat$n_eff,beta_se=sumdat$beta_se,p=sumdat$p,beta=sumdat$beta)
  
  fwrite(sumdat,file=paste0("../result/",outprefix,"_sumstats.txt"),row.names = F,sep="\t",quote=F)
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

pvar_onco=fread("../result/imp_onco/asian/asian.pvar")
pvar_icogs=fread("../result/imp_icogs/asian/asian.pvar")
asncomsnps=intersect(pvar_icogs$ID,pvar_onco$ID)
length(asncomsnps)
# [1] 9332432
asnmetafile="../result/meta_asian_lr.meta.rsid"
asnsumdat=as.data.frame(fread(asnmetafile))
table(asncomsnps %in% asnsumdat$SNP[asnsumdat$N==2])
table(asnsumdat$SNP[asnsumdat$N==2] %in% asncomsnps)
selsnps=asnsumdat$SNP[!is.na(asnsumdat$rsid) & asnsumdat$N==2]
length(selsnps) #8996998
asn_icogs_freq=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_icogs.afreq"))
asn_onco_freq=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_onco.afreq"))
tmp1=asn_icogs_freq %>% filter(ALT_FREQS<0.99 & ALT_FREQS>0.01) %>% select(ID)
tmp2=asn_onco_freq %>% filter(ALT_FREQS<0.99 & ALT_FREQS>0.01) %>% select(ID)
tmp=unique(c(tmp1$ID,tmp2$ID))
selsnps=intersect(selsnps,tmp)
length(selsnps) #8423019

asn_icogs_missing=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_icogs.vmiss"))
asn_onco_missing=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/asian_training_onco.vmiss"))
#tmp1=asn_icogs_missing %>% filter(F_MISS<0.1) %>% select(ID)
tmp11=asn_icogs_missing %>% filter(F_MISS>0.9) %>% select(ID)
#tmp2=asn_onco_missing %>% filter(F_MISS<0.1) %>% select(ID)
tmp22=asn_onco_missing %>% filter(F_MISS>0.5) %>% select(ID)
# tmp=unique(c(tmp1$ID,tmp2$ID))
# tmp=tmp[!tmp %in% tmp22]
# selsnps=intersect(selsnps,tmp)
selsnps=selsnps[!selsnps %in% c(tmp11$ID,tmp22$ID)]
length(selsnps) #8037079
idx=match(selsnps,asnsumdat$SNP)
asnsumdat1=asnsumdat[idx,]
asn_neff=as.integer(4*(1/(1/5940+1/4882)+1/(1/(6859+1/7545)))) #for LDpred2
formsumstats(sumdat=asnsumdat1,outprefix="asian_training",neff=asn_neff)

#generate tuning/validation data
eursumdat1=as.data.frame(fread("../result/euro_training_sumstats.txt"))
asnsumdat1=as.data.frame(fread("../result/asian_training_sumstats.txt"))
alltrainsnps=unique(c(eursumdat1$SNP,asnsumdat1$SNP))
length(alltrainsnps) #10209527
#write.table(alltrainsnps,file="../result/PRS/trainsnps.txt",row.names = F,col.names = F,quote=F)
#to map SNP to rsid
tmp1=eursumdat1[,c("SNP","rsid")]
tmp1$type="euro"
tmp2=asnsumdat1[,c("SNP","rsid")]
tmp2$type="asian"
tmp2=tmp2[!tmp2$SNP %in% tmp1$SNP,]
tmp=rbind(tmp1,tmp2)
idx=which(duplicated(tmp$rsid))
idx=which(tmp$rsid %in% tmp$rsid[idx])
tmp1=tmp[idx,] #5868 SNPs have duplicated rsid between eur and asn
idx1=match(tmp1$SNP,eursumdat1$SNP)
idx2=match(tmp1$SNP,asnsumdat1$SNP)
intersect(idx1[!is.na(idx1)],idx2[!is.na(idx2)])
tmp1$p_euro=eursumdat1$p[idx1]
tmp1$p_asian=asnsumdat1$p[idx2]
tmp1$p=rowSums(tmp1[,c("p_euro","p_asian")],na.rm=T)
tmp2=tmp1 %>% group_by(rsid) %>% arrange(p) %>% filter(row_number()>1) #SNPs to be removed
eur_snp_2rm=tmp2$SNP[tmp2$type=="euro"]
asn_snp_2rm=tmp2$SNP[tmp2$type=="asian"]
eursumdat2=eursumdat1[!eursumdat1$SNP %in% eur_snp_2rm,]
nrow(eursumdat2) #8707002
write.table(eursumdat2,file="../result/euro_training_sumstats.txt",row.names = F,sep="\t",quote=F)
asnsumdat2=asnsumdat1[!asnsumdat1$SNP %in% asn_snp_2rm,]
nrow(asnsumdat2) #8035565
write.table(asnsumdat2,file="../result/asian_training_sumstats.txt",row.names = F,sep="\t",quote=F)
rm(eursumdat,eursumdat1,asnsumdat,asnsumdat1)
alltrainsnps=unique(c(eursumdat2$SNP,asnsumdat2$SNP))
length(alltrainsnps) #10206593
write.table(alltrainsnps,file="../result/PRS/trainsnps.txt",row.names = F,col.names = F,quote=F)
#to change to rsid
tmp1=eursumdat2[,c("SNP","rsid")]
tmp2=asnsumdat2[,c("SNP","rsid")]
tmp2=tmp2[!tmp2$SNP %in% tmp1$SNP,]
tmp=rbind(tmp1,tmp2)
sum(duplicated(tmp$rsid)) #0
write.table(tmp,file="../result/PRS/trainsnps_updatename.txt",row.names = F,col.names = F,quote=F,sep="\t")
# cmd=paste0(plink2," --pfile ../result/imputation/merged --update-name ../result/merged_updatename.txt --make-pgen --out ../result/imputation/merged_rsid")
# system(cmd)

formgenotype=function(phenofile="../result/PRS_onco_euro_tuningpheno.txt",
                      prefix="../result/imp_onco/euro/euro",
                      outprefix="../result/PRS/euro_onco_tuning")
{
  pheno=read.table(phenofile,header=T,sep="\t")
  if (sum(colnames(pheno)=="Onc_ID")>0)
  {
    tmp=data.frame(FID=pheno$Onc_ID,IID=pheno$Onc_ID)
  }else
  {
    tmp=data.frame(FID=pheno$SG_ID,IID=pheno$SG_ID)
  }
  
  samplefile=paste0("../result/PRS/",basename(outprefix),".samples")
  write.table(tmp,file=samplefile,row.names = F,col.names = F,sep="\t",quote=F)
  cmd=paste0(plink2," --pfile ",prefix," --keep ",samplefile," --extract ../result/PRS/trainsnps.txt --make-pgen --out ",outprefix," --memory 100000 --threads 8 ")
  system(cmd)
  cmd=paste0(plink2," --pfile ",outprefix, " --update-name ../result/PRS/trainsnps_updatename.txt --make-pgen --out ",outprefix," --memory 100000 --threads 8 ")
  system(cmd)
  cmd=paste0(plink2," --pfile ",outprefix, " --make-bed --out ",outprefix," --memory 100000 --threads 8 ")
  system(cmd)
}
formgenotype()
formgenotype(phenofile="../result/PRS_onco_asian_tuningpheno.txt",
                      prefix="../result/imp_onco/asian/asian",
                      outprefix="../result/PRS/asian_onco_tuning")
formgenotype(phenofile="../result/PRS_icogs_african_validationpheno.txt",
             prefix="../result/imp_icogs/african/african",
             outprefix="../result/PRS/african_icogs_validation")
formgenotype(phenofile="../result/PRS_onco_african_validationpheno.txt",
             prefix="../result/imp_onco/african/african",
             outprefix="../result/PRS/african_onco_validation")
formgenotype(phenofile="../result/PRS_onco_asian_validationpheno.txt",
             prefix="../result/imp_onco/asian/asian",
             outprefix="../result/PRS/asian_onco_validation")
formgenotype(phenofile="../result/PRS_onco_euro_validationpheno.txt",
             prefix="../result/imp_onco/euro/euro",
             outprefix="../result/PRS/euro_onco_validation")
formgenotype(phenofile="../result/PRS_onco_other_validationpheno.txt",
             prefix="../result/imp_onco/hispanic/hispanic",
             outprefix="../result/PRS/hispanic_onco_validation")



