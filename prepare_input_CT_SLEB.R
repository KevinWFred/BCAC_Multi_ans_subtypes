#!/usr/bin/env Rscript
#prapare input fies for CT-SLEB

#transform sumdat to work with CTSLEB,use rsid
#EAStarget use EUR as reference, EURtarget use EAS as reference (allele codings are in agreement),make sumstats consistent with genotpye data
generate_sumdat=function(sumdatfile="../result/PRS1/euro_training_sumstats.txt",prefix="EUR",bimprefix="../result/PRS1/euro_onco_tuning")
{
  sumdat=fread(sumdatfile)
  idx=which(colnames(sumdat) %in% c("pos","position"))
  if (length(idx)>0) colnames(sumdat)[idx]="BP"
  idx=which(colnames(sumdat) %in% c("p","Pvalue_fixed"))
  if (length(idx)>0) colnames(sumdat)[idx]="P"
  idx=which(colnames(sumdat) %in% c("beta_se","StdError_fixed"))
  if (length(idx)>0) colnames(sumdat)[idx]="SE"
  idx=which(colnames(sumdat)=="chr")
  if (length(idx)>0) colnames(sumdat)[idx]="CHR"
  idx=which(colnames(sumdat)=="OR_fixed")
  if (length(idx)>0) colnames(sumdat)[idx]="OR"
  # idx=which(colnames(sumdat) %in% c("rs_number","rsid"))
  # if (length(idx)>0) colnames(sumdat)[idx]="SNP"
  sumdat$SNP=sumdat$rsid
  idx=which(colnames(sumdat) %in% c("effect_allele","a1"))
  if (length(idx)>0) colnames(sumdat)[idx]="A1"
  idx=which(colnames(sumdat) %in% c("reference_allele","a0"))
  if (length(idx)>0) colnames(sumdat)[idx]="A2"
  sumdat$BETA=log(sumdat$OR)
  
  bim=as.data.frame(fread(paste0(bimprefix,".bim")))
  idx=match(sumdat$SNP,bim$V2)
  print(sum(is.na(idx)))
  bim=bim[idx,]
  bim$name=paste0(bim$V1,":",bim$V4,":",bim$V5,":",bim$V6)
  sumdat$name=paste0(sumdat$CHR,":",sumdat$BP,":",sumdat$A1,":",sumdat$A2)
  sumdat$name2=paste0(sumdat$CHR,":",sumdat$BP,":",sumdat$A2,":",sumdat$A1)
  sum(sumdat$name %in% bim$name)
  sum(sumdat$name2 %in% bim$name)
  print(table(sumdat$name %in% bim$name |sumdat$name2 %in% bim$name))
  idx=which(sumdat$name2 %in% bim$name)
  tmp=sumdat$A1[idx]
  sumdat$A1[idx]=sumdat$A2[idx]
  sumdat$A2[idx]=tmp
  sumdat$BETA[idx]=-sumdat$BETA[idx]
  sumdat$OR[idx]=1/sumdat$OR[idx]
  sumdat$name=paste0(sumdat$CHR,":",sumdat$BP,":",sumdat$A1,":",sumdat$A2)
  print(table(sumdat$name==bim$name))
  sumdat1=data.frame(CHR=sumdat$CHR,SNP=sumdat$SNP,BP=sumdat$BP,A1=sumdat$A1,BETA=sumdat$BETA,SE=sumdat$SE,P=sumdat$P,rs_id=sumdat$rsid)
  fwrite(sumdat1,file=paste0("../result/PRS1/ctsleb/",prefix,"/",prefix,"_CTSLEB.sumdat"),row.names = F,sep="\t",quote=F)
  
}

generate_sumdat(sumdatfile="../result/PRS1/euro_training_sumstats.txt",prefix="EUR")
generate_sumdat(sumdatfile="../result/PRS1/asian_training_sumstats.txt",prefix="EAS",bimprefix="../result/PRS1/asian_onco_tuning")

#merge tuning and validation data into one, to be used for CTSLEB, to compute all the PRS together.
prefix_val1="../result/PRS1/african_onco_validation"
prefix_val2="../result/PRS1/african_icogs_validation"
prefix_val3="../result/PRS1/asian_onco_validation"
prefix_val4="../result/PRS1/euro_onco_validation"
prefix_val5="../result/PRS1/hispanic_onco_validation"
validationprefix=c(prefix_val1,prefix_val2,prefix_val3,prefix_val4,prefix_val5)
allprefix=c("../result/PRS1/euro_onco_tuning","../result/PRS1/asian_onco_tuning",
            validationprefix)
write.table(allprefix,file="../result/PRS1/tuning_validation_plink_mergelist.txt",row.names = F,col.names = F,quote=F)
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
cmd=paste0(plink," --keep-allele-order --allow-no-sex --merge-list ../result/PRS1/tuning_validation_plink_mergelist.txt --make-bed --out ../result/PRS1/tuning_validation --memory 100000 --threads 8")
system(cmd)
# #genotype data align with 1kG A1
# bim=fread("../../tools/CTSLEB/EAS/chr_all.bim")
# #use chr:pos:A1:A2 to align sumdat and bim
# bim$name=bim$V2
# idx=which(grepl("^rs",bim$V2))
# bim$name[idx]=paste0(bim$V1[idx],":",bim$V4[idx],":",bim$V5[idx],":",bim$V6[idx])
# bim$SNP=NA
# tmp=unlist(strsplit(bim$V2[idx],":"))
# bim$SNP[idx]=tmp[seq(1,length(tmp),4)]
# pvar=fread("/gpfs/gsfs12/users/BB_Bioinformatics/Kevin/PRS_EASLC/data/FLCCA/merged_filterd.pvar") 
# colnames(pvar)=c("CHR", "BP", "SNP", "A1", "A2")
# pvar$name=paste0(pvar$CHR,":",pvar$BP,":",pvar$A1,":",pvar$A2)
# pvar$name2=paste0(pvar$CHR,":",pvar$BP,":",pvar$A2,":",pvar$A1)
# sum(pvar$name %in% bim$name)
# sum(pvar$name2 %in% bim$name)
# sum(pvar$SNP %in% bim$SNP)
# idx=which(pvar$name2 %in% bim$name)
# tmp=pvar$A1[idx]
# pvar$A1[idx]=pvar$A2[idx]
# pvar$A2[idx]=tmp
# pvar$name[idx]=pvar$name2[idx]
# pvar=pvar[pvar$name %in% bim$name,]
# #remove ambiguous snps
# idx=(pvar$A1=="A" & pvar$A2=="T")|(pvar$A1=="T" & pvar$A2=="A")|(pvar$A1=="G" & pvar$A2=="C")|(pvar$A1=="C" & pvar$A2=="G")
# pvar=pvar[!idx,]
# fwrite(pvar[,c("SNP", "A1")], "../data/FLCCA/merged.bim_pvar_1kg.a1", col.names=F, sep="\t")
# write.table(pvar$SNP,"../data/FLCCA/merged.bim_pvar_1kg.snps",row.names = F,col.names = F,quote=F)
# cmd=paste0("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../data/FLCCA/merged_filterd --make-pgen ",
#            "--out ../data/FLCCA/plink2_1kg ",
#            "--ref-allele ../data/FLCCA/merged.bim_pvar_1kg.a1 ",
#            "--extract ../data/FLCCA/merged.bim_pvar_1kg.snps")
# system(cmd)
# #updata snp name to rsid:bp:a1:a2
# tmp=fread("../data/FLCCA/plink2_1kg.pvar")
# tmp$ID=paste0(tmp$ID,":",tmp$POS,":",tmp$REF,":",tmp$ALT)
# fwrite(tmp,file="../data/FLCCA/plink2_1kg.pvar",row.names = F,sep="\t",quote=F)
# cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../data/FLCCA/plink2_1kg --keep ../data/FLCCA/EAS_never_train_plinksample.txt ",
#           "--make-pgen --out ../result/FLCCA_1kg_train")
# system(cmd)
# 
# cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../data/FLCCA/plink2_1kg --keep ../data/FLCCA/EAS_never_train_plinksample.txt ",
#           "--make-bed --out ../result/FLCCA_1kg_train")
# system(cmd)
# 
# cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../data/FLCCA/plink2_1kg --keep ../data/FLCCA/EAS_never_test_plinksample.txt ",
#           "--make-pgen --out ../result/FLCCA_1kg_test")
# system(cmd)
# 
# cmd=paste("/usr/local/apps/plink/2.3-alpha/plink2 --pfile ../data/FLCCA/plink2_1kg --keep ../data/FLCCA/EAS_never_test_plinksample.txt ",
#           "--make-bed --out ../result/FLCCA_1kg_test")
# system(cmd)
# 
# #the 1kg reference data is too big to clump
# tmp=fread("../result/FLCCA_1kg_train.bim")
# write.table(tmp$V2,file="../result/FLCCA_1kg.snps",row.names = F,col.names=F,quote=F)
# cmd=paste0(plink2," -bfile ../../tools/CTSLEB/EUR/chr_all --extract ../result/FLCCA_1kg.snps --make-bed --out ../../tools/CTSLEB/EUR/chr_all_FLCCA")
# system(cmd)
