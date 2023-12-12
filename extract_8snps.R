#!/usr/bin/env Rscript

.libPaths(c("/data/wangx53",.libPaths()))
library(data.table)
library(readxl)
snpdat=as.data.frame(read_excel("../data/eight_snps_info.xlsx"))
#change coordinates to hg38
library(rtracklayer)
library(GenomicRanges)
chain=import.chain("../../tools/liftover/hg19ToHg38.over.chain")

gr_dat=GRanges(seqnames = paste0("chr",snpdat$Chr.2),ranges=IRanges(start=snpdat$Position,width=1))
tmp=liftOver(gr_dat,chain)
tmp=as.data.frame(tmp)
snpdat$Position=tmp$start
tmp=unlist(strsplit(snpdat$Alleles3,"/"))
snpdat$allele1=tmp[seq(1,length(tmp),2)]
snpdat$allele2=tmp[seq(2,length(tmp),2)]
# snpdat$name1=paste0("chr",snpdat$Chr.2,":",snpdat$Position,":",snpdat$allele1,":",snpdat$allele2)
# snpdat$name2=paste0("chr",snpdat$Chr.2,":",snpdat$Position,":",snpdat$allele2,":",snpdat$allele1)
# snpdat$chrposname=paste0(snpdat$Chr.2,":",snpdat$Position)


plink="/usr/local/apps/plink/1.9/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
infolder="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/"
bim=as.data.frame(fread(paste0(infolder,"euro.pvar")))
sum(snpdat$name1 %in% bim$ID) #4
sum(snpdat$name2 %in% bim$ID) #1
tmp=paste0(bim$`#CHROM`,":",bim$POS)
sum(snpdat$chrposname %in% tmp) #5

infolder="/data/BB_Bioinformatics/ProjectData/BCAC/icogs/"
prefix="zhang_750_euro_icogs_topmed_"
outfolder="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/"
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

## this step is optional!
## here we just simplify the names of the objects, making the code neater
# genome <- BSgenome.Hsapiens.UCSC.hg38
# all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
# seqlevelsStyle(genome) <- "NCBI"

extractsnp=function(infolder,outfolder,prefix)
{
  #gerante range file
  # tmp=data.frame(chr=snpdat$Chr.2,start=snpdat$Position,end=NA,name=snpdat$`Lead variant1`)
  # for (i in 1:nrow(tmp))
  # {
  #   tmp$end[i]=max(nchar(snpdat$allele1[i]),nchar(snpdat$allele2[i]))-1+tmp$start[i]
  # }
  # snpinfo=snpsById(all_snps,snpdat$`Lead variant1`,ifnotfound="drop")
  # snps=snpinfo@elementMetadata@listData$RefSNP_id
  # tmp=snpinfo@seqnames
  # chrs=rep(as.character(tmp@values),times=tmp@lengths)
  # pos=snpinfo@ranges@pos
  # snppos=data.frame(snp=snps,chr=chrs,pos=pos)
  # missingsnps=snpdat$`Lead variant1`[!snpdat$`Lead variant1` %in% snppos$snp]
  library(biomaRt)
  mart2 = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")
  
  results = getBM(attributes = c("refsnp_id","chr_name",
                                 "chrom_start","chrom_end","allele"),
                  filters = "snp_filter", values = snpdat$`Lead variant1`, mart =
                    mart2)
  snppos=data.frame(chr=results$chr_name,start=results$chrom_start,end=results$chrom_end,snp=results$refsnp_id)
  missingsnps=snpdat$`Lead variant1`[!snpdat$`Lead variant1` %in% snppos$snp]
  if (length(missingsnps)>0) warning(paste0(missingsnps," not found"))
  missingpos=data.frame(chr=1,start=220558134,end=220558136,snp=missingsnps) #https://genome.ucsc.edu
  snppos=rbind(snppos,missingpos)
  #write.table(snppos,file="../result/eightsnps.range",row.names = F,col.names = F,sep=" ",quote=F)
  allchrs=unique(snpdat$Chr.2)
  for (i in 1:length(allchrs))
  {
    chr=allchrs[i]
    allfiles=list.files(infolder,paste0(prefix,chr,"_\\w*.bgen"))
    if (length(allfiles)==0) allfiles=list.files(infolder,paste0(prefix,chr,".bgen"))
    for (j in 1:length(allfiles))
    {
      prefix1=gsub(".bgen","",allfiles[j])
      cmd=paste0(plink2," --bgen ",infolder,prefix1,".bgen ref-unknown --extract range ../result/eightsnps.range --make-pgen --out ",outfolder,"eightsnp_chr",chr,"_",j,
                 " --memory 128000 --threads 24")
      system(cmd)
    }
  }
}
infolder="/data/BB_Bioinformatics/ProjectData/BCAC/onco/"
prefix="zhang_750_topmed_"
outfolder="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/"

mergedat=function(outfolder)
{
  allfiles=list.files(outfolder,"eightsnp_chr\\w*.pvar")
  tmp=data.frame(prefix=paste0(outfolder,allfiles))
  tmp$prefix=gsub(".pvar","",tmp$prefix)
  write.table(tmp,file=paste0(outfolder,"eightsnp_merglist.txt"),row.names = F,col.names = F,quote=F)
  cmd=paste0(plink2," --pmerge-list ",outfolder,"eightsnp_merglist.txt --make-pgen --out ",outfolder,"eightsnp")
  system(cmd)
  cmd=paste0(plink2," --pfile ",outfolder,"eightsnp --recode A-transpose --out ",outfolder,"eightsnp")
  system(cmd)
}

getgendat=function(outfolder)
{
  gendat=as.data.frame(fread(paste0(outfolder,"eightsnp.traw"))) #the last allele in SNPID is the counted one
  colnames(gendat)=gsub("^0_","",colnames(gendat))
  gendat1=gendat[,7:ncol(gendat)]
  tmp=unlist(strsplit(colnames(gendat1),"_"))
  colnames(gendat1)=tmp[seq(1,length(tmp),2)]
  rownames(gendat1)=gendat$SNP
  AF=rowSums(gendat1)/2/ncol(gendat1)
  idx=which(AF>0.01)
  AF=AF[idx]
  gendat=gendat[idx,]
  gendat1=gendat1[idx,]
  idx=match(paste0(gendat$CHR,":",gendat$POS),paste0(snpdat$Chr.2,":",snpdat$Position))
  snpid=snpdat$`Lead variant1`[idx]
  snpinfo=data.frame(rsid=snpid,alleles=snpdat$Alleles3[idx],AF=snpdat$MAF4[idx],AF_dat=AF,ID=rownames(gendat1))
  return(list(snpinfo=snpinfo,dat=gendat1))
}

oncodat=getgendat(outfolder = "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/")
icogsdat=getgendat(outfolder = "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/")

library(readxl)
pheno_icogs=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_icogs_pheno_v15_02.txt",header=T,sep="\t")
pheno_onco=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_onco_pheno_v15_02.txt",header=T,sep="\t")

oncodat1=pheno_onco[,c(1,3,5:14,27,28,34,39:41,43,45)]
for(i in 13:ncol(oncodat1))
{
  idx=which(oncodat1[,i]==888)
  oncodat1[idx,i]=NA
}
oncodat1$Behaviour1[is.na(oncodat1$Behaviour1)]=0
oncodat2=t(oncodat$dat)
idx=match(rownames(oncodat2),oncodat1$Onc_ID)
oncodat2=cbind(oncodat1[idx,],oncodat2)
write.csv(oncodat2,file="../result/eightsnp_oncodata.csv",row.names = F,quote=F)

icogsdat1=pheno_icogs[,c(1,3,6:15,28,29,35,40:42,44,46)]
for(i in 13:ncol(icogsdat1))
{
  idx=which(icogsdat1[,i]==888)
  icogsdat1[idx,i]=NA
}
icogsdat1$Behaviour1[is.na(icogsdat1$Behaviour1)]=0
icogsdat2=t(icogsdat$dat)
idx=match(rownames(icogsdat2),icogsdat1$SG_ID)
icogsdat2=cbind(icogsdat1[idx,],icogsdat2)
write.csv(icogsdat2,file="../result/eightsnp_icogsdata.csv",row.names = F,quote=F)

#for this snp
snp="chr9:119662503:T:C"
checkpos=data.frame(chr=9,start=119662502,end=119662504,snp="chr9:119662503:T:C")
write.table(checkpos,file="../result/tmp/checksnp.range",row.names = F,col.names = F,sep=" ",quote=F)
infolder="/data/BB_Bioinformatics/ProjectData/BCAC/onco/"
prefix1="zhang_750_asian_topmed_9"
cmd=paste0(plink2," --bgen ",infolder,prefix1,".bgen ref-unknown --extract range ../result/tmp/checksnp.range --make-pgen --out ../result/tmp/checkpos ",
           " --memory 28000 --threads 10")
system(cmd)
cmd=paste0(plink2," --pfile ../result/tmp/checkpos --recode A-transpose --out ../result/tmp/checkpos")
system(cmd)
tmp=as.data.frame(fread("../result/tmp/checkpos.traw"))
tmp=tmp[,7:ncol(tmp)]
colnames(tmp)=gsub("^0_","",colnames(tmp))
tmp1=unlist(strsplit(colnames(tmp),"_"))
colnames(tmp)=tmp1[seq(1,length(tmp1),2)]
pheno=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
tmp=tmp[,colnames(tmp) %in% pheno$Onc_ID]
idx=match(colnames(tmp),pheno$Onc_ID)
idx1=which(is.na(pheno$Behaviour1[idx]))
#MAF 0.00487
sum(unlist(tmp[,idx1]))/2/length(idx1)  
sum(unlist(tmp[]))/2/ncol(tmp) #0.0044

infolder="/data/BB_Bioinformatics/ProjectData/BCAC/onco/"
prefix1="zhang_750_topmed_9_p3"
cmd=paste0(plink2," --bgen ",infolder,prefix1,".bgen ref-unknown --extract range ../result/tmp/checksnp.range --make-pgen --out ../result/tmp/checkpos1 ",
           " --memory 28000 --threads 10")
system(cmd)
cmd=paste0(plink2," --pfile ../result/tmp/checkpos1 --recode A-transpose --out ../result/tmp/checkpos1")
system(cmd)
tmp=as.data.frame(fread("../result/tmp/checkpos1.traw"))
tmp=tmp[,7:ncol(tmp)]
colnames(tmp)=gsub("^0_","",colnames(tmp))
tmp1=unlist(strsplit(colnames(tmp),"_"))
colnames(tmp)=tmp1[seq(1,length(tmp1),2)]
tmp=tmp[,colnames(tmp) %in% pheno$Onc_ID]
idx=match(colnames(tmp),pheno$Onc_ID)
idx1=which(is.na(pheno$Behaviour1[idx]))
#MAF 0.003637787
sum(unlist(tmp[,idx1]))/2/length(idx1)  
sum(unlist(tmp[]))/2/ncol(tmp) #0.003652263

prefix1="zhang_750_african_topmed_9"
cmd=paste0(plink2," --bgen ",infolder,prefix1,".bgen ref-unknown --extract range ../result/tmp/checksnp.range --make-pgen --out ../result/tmp/checkpos2 ",
           " --memory 28000 --threads 10")
system(cmd)
cmd=paste0(plink2," --pfile ../result/tmp/checkpos2 --recode A-transpose --out ../result/tmp/checkpos2")
system(cmd)
tmp=as.data.frame(fread("../result/tmp/checkpos2.traw"))
tmp=tmp[,7:ncol(tmp)]
colnames(tmp)=gsub("^0_","",colnames(tmp))
tmp1=unlist(strsplit(colnames(tmp),"_"))
colnames(tmp)=tmp1[seq(1,length(tmp1),2)]
tmp=tmp[,colnames(tmp) %in% pheno$Onc_ID]
idx=match(colnames(tmp),pheno$Onc_ID)
idx1=which(is.na(pheno$Behaviour1[idx]))
#MAF 0.00174171
sum(unlist(tmp[,idx1]))/2/length(idx1)  
sum(unlist(tmp[]))/2/ncol(tmp) #0.001401526
