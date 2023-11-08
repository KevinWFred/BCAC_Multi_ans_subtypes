#!/usr/bin/env Rscript
#aggregate logodds and sigma across populations
#get the final global test p-values (ACAT p-values)
.libPaths(c("/data/wangx53",.libPaths()))
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
library(data.table)
library(bc2)
startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
# read_metap=function(pop="euro")
# {
#   allbeta=NULL
#   allsigma=list()
#   for (i in 1:1000)
#   {
#     if (i %% 100 ==0) cat(i,'..')
#     load(paste0("../result/metaresult/",pop,"/res_",i,".RData"))
#     allbeta=rbind(allbeta,allmetares$metares)
#     allsigma=c(allsigma,allmetares$metaallsigma)
#     
#   }
#   res=list(allbeta=allbeta,allsigma=allsigma)
#   save(res,file=paste0("../result/",pop,"_meta_beta_sigma.RData"))
# }
#for euro
pop="euro"
load(paste0("../result/",pop,"_beta_sigma.RData"))
euro_allsigma_icogs=allsigma_icogs
euro_allsigma_onco=allsigma_onco
euro_allbeta_icogs=beta_icogs
euro_allbeta_onco=beta_onco

pop="asian"
load(paste0("../result/",pop,"_beta_sigma.RData"))
asian_allsigma_icogs=allsigma_icogs
asian_allsigma_onco=allsigma_onco
asian_allbeta_icogs=beta_icogs
asian_allbeta_onco=beta_onco
#for african
pop="african"
#save(beta_onco,sigma_onco,allsigma_onco,file=paste0(outprefix,"_beta_sigma.RData"))
load(paste0("../result/",pop,"_beta_sigma.RData"))
african_allsigma_onco=allsigma_onco
african_allbeta_onco=beta_onco

rm(allsigma_icogs,allsigma_onco,beta_icogs,beta_onco,sigma_icogs,sigma_onco)
gc()
#check how much size used
#sum( sapply(ls(),function(x){object.size(get(x))}))/1E6

LogoddsMetaAnalysis <- function(logoddslist,sigmalist){
  if (length(logoddslist)>1)
  {
    sigma.inv.sum=0
    sigma.inv.logodds.sum=0
    for(i in 1:length(sigmalist))
    {
      sigma.inv.sum=sigma.inv.sum+solve(sigmalist[[i]])
      sigma.inv.logodds.sum=sigma.inv.logodds.sum+solve(sigmalist[[i]])%*%unlist(logoddslist[[i]])
    }
    sigma.meta=solve(sigma.inv.sum)
    
    logodds.meta <- sigma.meta%*% sigma.inv.logodds.sum
  }else #==1
  {
    logodds.meta=as.numeric(logoddslist[[1]])
    sigma.meta=sigmalist[[1]]
  }
  
  
  return(list(logodds.meta = logodds.meta,
              sigma.meta = sigma.meta))
  
}



# LogoddsMetaAnalysis2 <- function(logodds1,sigma1,logodds2,sigma2){
#   sigma1.inv <- solve(sigma1)
#   sigma2.inv <- solve(sigma2)
#   sigma.meta <- solve(sigma1.inv+sigma2.inv)
#   logodds.meta <- sigma.meta%*%(sigma1.inv%*%logodds1+sigma2.inv%*%logodds2 )
#   
#   return(list(logodds.meta = logodds.meta,
#               sigma.meta = sigma.meta))
# }

#Global test for association using logodds and sigma
GlobalTestForAssoc <- function(logodds,sigma){
  sigma <- as.matrix(sigma)
  df <- length(logodds)
  GTA.stat <- t(logodds)%*%solve(sigma)%*%logodds
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  if(p.value.GTA!=0)
  {
    places = 3
    power.number <- floor(-log10(p.value.GTA))+places
    ###format the output with three digits in total
    p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)
  }
  
  
  return(p.value.GTA)
}


allvar=unique(c(rownames(euro_allbeta_icogs),rownames(euro_allbeta_onco),
                rownames(asian_allbeta_icogs),rownames(asian_allbeta_onco),
                rownames(african_allbeta_onco)))

args = commandArgs(trailingOnly = T)
print(Sys.time())
print(args)

i1 = as.numeric(args[1]) #block ID, starts with 1
ntask=1000 #total number of jobs
print(i1)
startendidx=startend(length(allvar),ntask,i1)
print(paste0("startendidx: ",startendidx))
idxstart=startendidx[1]
idxend=startendidx[2]
allres=data.frame(matrix(NA,nrow=idxend-idxstart+1,ncol=1))
rownames(allres)=allvar[idxstart:idxend]
colnames(allres)="P" #combined P-value
euro_allbeta_icogs=euro_allbeta_icogs[rownames(euro_allbeta_icogs) %in% rownames(allres),]
euro_allbeta_onco=euro_allbeta_onco[rownames(euro_allbeta_onco) %in% rownames(allres),]
asian_allbeta_icogs=asian_allbeta_icogs[rownames(asian_allbeta_icogs) %in% rownames(allres),]
asian_allbeta_onco=asian_allbeta_onco[rownames(asian_allbeta_onco) %in% rownames(allres),]
african_allbeta_onco=african_allbeta_onco[rownames(african_allbeta_onco) %in% rownames(allres),]
euro_allsigma_icogs=euro_allsigma_icogs[names(euro_allsigma_icogs) %in% rownames(allres)]
euro_allsigma_onco=euro_allsigma_onco[names(euro_allsigma_onco) %in% rownames(allres)]
asian_allsigma_icogs=asian_allsigma_icogs[names(asian_allsigma_icogs) %in% rownames(allres)]
asian_allsigma_onco=asian_allsigma_onco[names(asian_allsigma_onco) %in% rownames(allres)]
african_allsigma_onco=african_allsigma_onco[names(african_allsigma_onco) %in% rownames(allres)]
gc()

for (i in 1:nrow(allres))
{
  if (i %% 1000 ==0) cat(i,'..')
  varname=rownames(allres)[i]
  logoddslist=sigmalist=list()
  if (varname %in% rownames(euro_allbeta_icogs))
  {
    logoddslist=c(logoddslist,list(euro_allbeta_icogs[which(rownames(euro_allbeta_icogs)==varname),]))
    sigmalist=c(sigmalist,euro_allsigma_icogs[which(names(euro_allsigma_icogs)==varname)])
  }
  if (varname %in% rownames(euro_allbeta_onco))
  {
    logoddslist=c(logoddslist,list(euro_allbeta_onco[which(rownames(euro_allbeta_onco)==varname),]))
    sigmalist=c(sigmalist,euro_allsigma_onco[which(names(euro_allsigma_onco)==varname)])
  }

  if (varname %in% rownames(asian_allbeta_icogs))
  {
    logoddslist=c(logoddslist,list(asian_allbeta_icogs[which(rownames(asian_allbeta_icogs)==varname),]))
    sigmalist=c(sigmalist,asian_allsigma_icogs[which(names(asian_allsigma_icogs)==varname)])
  }
  if (varname %in% rownames(asian_allbeta_onco))
  {
    logoddslist=c(logoddslist,list(asian_allbeta_onco[which(rownames(asian_allbeta_onco)==varname),]))
    sigmalist=c(sigmalist,asian_allsigma_onco[which(names(asian_allsigma_onco)==varname)])
  }
  if (varname %in% rownames(african_allbeta_onco))
  {
    logoddslist=c(logoddslist,list(african_allbeta_onco[which(rownames(african_allbeta_onco)==varname),]))
    sigmalist=c(sigmalist,african_allsigma_onco[which(names(african_allsigma_onco)==varname)])
  }
  
  metares=LogoddsMetaAnalysis(logoddslist,sigmalist)
  allres$P[i]=GlobalTestForAssoc(metares$logodds.meta,metares$sigma.meta)
}
save(allres,file=paste0("../result/metaresult/all/","res_",i1,".RData"))
print(Sys.time())
print("done")

# swarmjobs=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/compute_metapvalues.R",1000),
#                      i1=1:1000)
# write.table(swarmjobs,file=paste0("compute_metapvalues.swarm"),row.names=F,col.names=F,sep="\t",quote=F)
#9757248
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/compute_metapvalues.swarm -g 56 --module R/4.3 --partition=quick --time=4:00:00 --gres=lscratch:40 -p 2

#collect the results
allpvalues=NULL
for (i in 1:1000)
{
  if (i %% 100 ==0) cat(i,'..')
  load(paste0("../result/metaresult/all","/res_",i,".RData"))
  allpvalues=rbind(allpvalues,allres)
}
intrinsicpvalues=allpvalues

load("../result/scoretestresult/allpvalues.Rdata")
scoretestallpvaues=data.frame(P=allpvalues)
rownames(scoretestallpvaues)=names(allpvalues)
tmp=intersect(rownames(scoretestallpvaues),rownames(intrinsicpvalues))
library(ACAT)
allpvalues=cbind(scoretestallpvaues[match(tmp,rownames(scoretestallpvaues)),],intrinsicpvalues[match(tmp,rownames(intrinsicpvalues)),])
idx=complete.cases(allpvalues)
allpvalues=allpvalues[idx,]
allpvalues=ACAT(t(allpvalues))
names(allpvalues)=tmp[idx]
#new scoretestallpvaues (MAF>0.006)-->allpvalues
save(intrinsicpvalues,scoretestallpvaues,allpvalues,file="../result/compute_metapvalues_new.RData")
#check novel variants:
library(GenomicRanges)
knownvar=read.csv("../data/218_known_discovery_snp_paper_order_101323.csv")
for(i in 1:nrow(knownvar)) knownvar$position[i]=gsub(",","",knownvar$position[i])
knownvar$position=as.integer(knownvar$position)
#change coordinates to hg38
library(rtracklayer)
chain=import.chain("../../tools/liftover/hg19ToHg38.over.chain")
gr_dat=GRanges(seqnames = paste0("chr",knownvar$CHR),ranges=IRanges(start=knownvar$position,width=1))
tmp=liftOver(gr_dat,chain)
tmp=as.data.frame(tmp)
knownvar$hg38position=tmp$start
table(knownvar$hg38position-5e5>0)
tmp=unlist(strsplit(names(allpvalues),":"))
chr=tmp[seq(1,length(tmp),4)]
chr=as.integer(gsub("chr","",chr))
pos=as.integer(tmp[seq(2,length(tmp),4)])
gr_allpvalues=GRanges(seqnames = chr,ranges = IRanges(start=pos,width = 10))
gr_knownvar=GRanges(seqnames = knownvar$CHR,ranges=IRanges(start=pmax(0,knownvar$hg38position-5e5),width=1e6))
gr_knownvar1=GRanges(seqnames = knownvar$CHR,ranges=IRanges(start=knownvar$hg38position,width=1))
#dist=distance(gr_allpvalues,gr_knownvar)
neardist=distanceToNearest(gr_allpvalues,gr_knownvar)
neardist1=distanceToNearest(gr_allpvalues,gr_knownvar1)
sum(neardist@elementMetadata@listData$distance==0) #1504407
idx=which(neardist@elementMetadata@listData$distance!=0)
novelpvalues=allpvalues[idx]
sum(novelpvalues<5e-8) #627
idx1=which(novelpvalues<5e-8)
tmp=data.frame(snp=names(novelpvalues[idx1]),p=novelpvalues[idx1])
write.table(tmp,file="../result/novelpvalues.txt",row.names = F,sep="\t",quote=F)
tmp=data.frame(snp=names(novelpvalues[idx1]))
write.table(tmp,file="../result/novelsnp.txt",row.names = F,sep="\t",quote=F)

#find independent novel regions
novelpvaules1=novelpvalues[which(novelpvalues<5e-8)]
idx=match(names(novelpvaules1),names(allpvalues))
chr1=chr[idx]
pos1=pos[idx]
gr_novel1=GRanges(seqnames = chr1,ranges=IRanges(start=pos1,width=1))
independentnovelvar=NULL
n=1
while (min(novelpvaules1)<5e-8 & n<length(novelpvaules1))
{
  idx=which.min(novelpvaules1)
  if (novelpvaules1[idx]<5e-8) independentnovelvar=c(independentnovelvar,novelpvaules1[idx])
  novelpvaules1[idx]=1
  gr_snp=GRanges(seqnames = chr1[idx],ranges=IRanges(start=pos1[idx],width=1))
  dist1=distance(gr_novel1,gr_snp)
  idx1=which(dist1<1e6)
  if(length(idx1)>0) novelpvaules1[c(idx1)]=1
  n=n+1
}
length(independentnovelvar) #30
tmp=data.frame(snp=names(independentnovelvar),p=independentnovelvar)

tmp=data.frame(snp=names(independentnovelvar))
write.table(tmp,file="../result/independentnovelsnp.txt",row.names = F,sep="\t",quote=F)

#LD pruning
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
pheno_onco=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_onco_pheno_v15_02_corrected.txt",header=T,sep="\t")
fam=read.table("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1.psam")
tmp=unlist(strsplit(fam$V1,"_"))
fam$V3=tmp[seq(1,length(tmp),2)]
idx=match(fam$V3,pheno_onco$Onc_ID)
table(pheno_onco$EthnicityGeno[idx])
onco_contrl=pheno_onco$Onc_ID[which(is.na(pheno_onco$Behaviour1))]
onco_contrl=intersect(onco_contrl,fam$V3)
onco_contrl1=paste0(onco_contrl,"_",onco_contrl)
write.table(onco_contrl1,file="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/onco_control.txt",row.names = F,col.names = F,quote=F)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1 --keep /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/onco_control.txt --extract ../result/independentnovelsnp.txt --make-bed --out ../result/onco_control_independentnovelsnp --memory 128000 --threads 8")
system(cmd)
cmd=paste0(plink," --bfile ../result/onco_control_independentnovelsnp --recode A-transpose --out ../result/onco_control_independentnovelsnp")
cmd=paste0(plink," --bfile ../result/onco_control_independentnovelsnp --freq --out ../result/onco_control_independentnovelsnp")

tmp=as.data.frame(fread("../result/onco_control_independentnovelsnp.bim"))
onco_pvar=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1.pvar"))
table(names(independentnovelvar) %in% onco_pvar$ID)
# FALSE  TRUE 
# 7    23 
oncomissing=names(independentnovelvar)[!names(independentnovelvar) %in% onco_pvar$ID]

pthr=1
r2thr=0.1
kbpthr=10000

# cmd=paste0(plink," --bfile ../result/onco_control_independentnovelsnp"," --clump ../result/novelpvalues.txt"," --clump-p1 ",
#            pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field snp --clump-field p --out ../result/onco_independentnovel")
cmd=paste0(plink," --bfile ../result/onco_control_independentnovelsnp"," --clump ../result/novelpvalues.txt"," --clump-p1 ",
           pthr," --clump-r2 ",r2thr," --clump-snp-field snp --clump-field p --out ../result/onco_independentnovel")

system(cmd)
#23 --> 23
onco_clumping=read.table("../result/onco_independentnovel.clumped",header=T)


icogs_pvar=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1.pvar"))
table(names(independentnovelvar) %in% icogs_pvar$ID)
# FALSE  TRUE 
# 1    29
table(oncomissing %in% icogs_pvar$ID)
pheno_icogs=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_icogs_pheno_v15_02.txt",header=T,sep="\t")
fam_icogs=read.table("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1.psam")
tmp=unlist(strsplit(fam_icogs$V1,"_"))
fam_icogs$V3=tmp[seq(1,length(tmp),2)]
idx=match(fam_icogs$V3,pheno_icogs$SG_ID)
table(pheno_icogs$EthnicityGeno[idx])
icogs_contrl=pheno_icogs$SG_ID[which(is.na(pheno_icogs$Behaviour1))]
icogs_contrl=intersect(icogs_contrl,fam_icogs$V3)
icogs_contrl1=paste0(icogs_contrl,"_",icogs_contrl)
write.table(icogs_contrl1,file="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/icogs_control.txt",row.names = F,col.names = F,quote=F)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1 --keep /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/icogs_control.txt --extract ../result/independentnovelsnp.txt --make-bed --out ../result/icogs_control_independentnovelsnp --memory 128000 --threads 8")
system(cmd)
cmd=paste0(plink," --bfile ../result/icogs_control_independentnovelsnp --recode A-transpose --out ../result/icogs_control_independentnovelsnp")
cmd=paste0(plink," --bfile ../result/icogs_control_independentnovelsnp --freq --out ../result/icogs_control_independentnovelsnp")

# cmd=paste0(plink," --bfile ../result/icogs_control_independentnovelsnp"," --clump ../result/novelpvalues.txt"," --clump-p1 ",
#            pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field snp --clump-field p --out ../result/icogs_independentnovel")
cmd=paste0(plink," --bfile ../result/icogs_control_independentnovelsnp"," --clump ../result/novelpvalues.txt"," --clump-p1 ",
           pthr," --clump-r2 ",r2thr," --clump-snp-field snp --clump-field p --out ../result/icogs_independentnovel")
system(cmd)
icogs_independentnovel=read.table("../result/icogs_control_independentnovelsnp.bim")
names(independentnovelvar)[!names(independentnovelvar) %in% icogs_independentnovel$V2]
#chr6:52468337:G:A
#30 --> 30
#total is 30 (29+1)
icogs_clumping=read.table("../result/icogs_independentnovel.clumped",header=T)
allnovelsnps=c("chr6:52468337:G:A",icogs_clumping$SNP)
allnovelsnps=data.frame(ID=allnovelsnps,chr=NA,pos=NA,rsid=NA,pvalue=NA)
idx=match(allnovelsnps$ID,names(allpvalues))
allnovelsnps$pvalue=allpvalues[idx]
tmp=unlist(strsplit(allnovelsnps$ID,":"))
allnovelsnps$chr=tmp[seq(1,length(tmp),4)]
allnovelsnps$chr=gsub("chr","",allnovelsnps$chr)
allnovelsnps$pos=as.integer(tmp[seq(2,length(tmp),4)])
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
seqlevelsStyle(genome) <- "NCBI"

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = allnovelsnps$chr, pos = allnovelsnps$pos)

## query the genome with positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)

idx = match(paste0(allnovelsnps$chr,":",allnovelsnps$pos),paste0(my_snps$seqnames,":",my_snps$pos))
allnovelsnps$rsid=my_snps$RefSNP_id[idx]
library(biomaRt)
ensembl <- useEnsembl("snp",dataset = "hsapiens_snp",GRCh = "37")
#get genomic position
tmp <- getBM(attributes=c("refsnp_id",
                           "chr_name",
                           "chrom_start",
                           "chrom_end"),
              filters ="snp_filter", values =allnovelsnps$rsid, mart = ensembl, uniqueRows=TRUE)
tmp=tmp[tmp$chr_name %in% 1:22,]
idx=match(tmp$refsnp_id,allnovelsnps$rsid)
allnovelsnps$hg19pos[idx]=tmp$chrom_start
#get hg19 position
library(rtracklayer)
chain=import.chain("../../tools/liftover/hg38ToHg19.over.chain")
allnovelsnps$hg19pos=NA
for (i in 1:nrow(allnovelsnps))
{
  gr_dat=GRanges(seqnames = paste0("chr",allnovelsnps$chr[i]),ranges=IRanges(start=allnovelsnps$pos[i],width=1))
  tmp=liftOver(gr_dat,chain)
  tmp=as.data.frame(tmp)
  if (nrow(tmp)>0 & sum(grepl("chr",tmp$seqnames))>0)
  {
    allnovelsnps$hg19pos[i]=tmp$start[1]
  }
}
idx=match(allnovelsnps$ID,names(metascoreinfo4$freq))
allnovelsnps$freq=metascoreinfo4$freq[idx]
write.csv(allnovelsnps,file="../result/allnovelsnps.csv",row.names = F,quote=F)
allnovelsnps=read.csv("../result/allnovelsnps.csv")

gr_allnovelsnps=GRanges(seqnames = allnovelsnps$chr,ranges = IRanges(start=allnovelsnps$pos,width=1))
tmp=distanceToNearest(gr_allnovelsnps,gr_knownvar1)
allnovelsnps$dist2nearestknown=tmp@elementMetadata@listData$distance
all(knownvar$hg38position==start(gr_knownvar1))
table(allnovelsnps$dist2nearestknown<2e6)
# FALSE  TRUE 
# 20    10
quantile(allnovelsnps$dist2nearestknown)
# 0%      25%      50%      75%     100% 
# 550583  1403882  2769633  6586420 28452496
#for those within 2MB of known variants, find the knownvar in +-2MB
allnovelsnps1=allnovelsnps[allnovelsnps$dist2nearestknown<2e6,]

allnovelsnps1$knownvar_rsid=allnovelsnps1$knownvar=allnovelsnps1$dist2nearestknown1=NA
allnovelsnps1$knownvar_pvalue=NA
gr_allnovelsnps1=GRanges(seqnames = allnovelsnps1$chr,ranges=IRanges(start=allnovelsnps1$pos,width=1))
# for(i in 1:nrow(allnovelsnps1))
# {
#   tmp1=distance(gr_knownvar,gr_allnovelsnps1[i])
#   idx=which.min(tmp1)
#   tmp2=distance(gr_knownvar[idx],gr_allpvalues)
#   idx1=which(tmp2==0)
#   if (length(idx1)>0)
#   {
#     allnovelsnps1$knownvar_rsid[i]=knownvar$Best.published.SNP[idx]
#     allnovelsnps1$knownvar[i]=names(allpvalues)[idx1]
#     allnovelsnps1$knownvar_pvalue[i]=allpvalues[idx1]
#     tmp3=distance(gr_allnovelsnps1[i],gr_allpvalues[idx1])
#     allnovelsnps1$dist2nearestknown1[i]=tmp3
#   }
# }
gr_knownvar2=intersect(gr_knownvar1,gr_allpvalues)
idx=match(gr_knownvar2,gr_knownvar1)
knownvar2=knownvar[idx,]
for(i in 1:nrow(allnovelsnps1))
{
  tmp1=distance(gr_knownvar2,gr_allnovelsnps1[i])
  idx=which(tmp1<2e6)
  tmp2=distanceToNearest(gr_knownvar2[idx],gr_allpvalues)
  idx1=tmp2@to[which(tmp2@elementMetadata$distance==0)]
  if (length(idx1)>0)
  {
    allnovelsnps1$knownvar_rsid[i]=paste0(knownvar2$Best.published.SNP[idx],collapse = ",")
    allnovelsnps1$knownvar[i]=paste0(names(allpvalues)[idx1],collapse = ",")
    allnovelsnps1$knownvar_pvalue[i]=paste0(allpvalues[idx1],collapse = ",")
    tmp3=distance(gr_allnovelsnps1[i],gr_allpvalues[idx1])
    allnovelsnps1$dist2nearestknown1[i]=paste0(tmp3,collapse = ",")
  }
}
#all the closest knowvar can be found in the result
write.table(allnovelsnps1,file="../result/allnovelsnps1.txt",row.names = F,sep="\t",quote=F)

tmp1=as.data.frame(fread("../result/onco_control_independentnovelsnp.traw"))
rownames(tmp1)=tmp1$SNP
tmp1=tmp1[,-c(1:6)]
tmp2=read.table("../result/onco_control_independentnovelsnp.frq",header = T)
quantile(tmp2$MAF)
# 0%       25%       50%       75%      100% 
# 0.0009262 0.0085740 0.0491500 0.2545750 0.4734000
idx1=which(rownames(tmp1)=="chr1:120953531:C:T")
idx2=which(rownames(tmp1)=="chr1:120425254:C:G")
idx3=which(rownames(tmp1)=="chr1:143209281:C:A")
cor(unlist(tmp1[idx1,]),unlist(tmp1[idx2,]),use="complete") #0.58
cor(unlist(tmp1[idx1,]),unlist(tmp1[idx3,]),use="complete") #0.48
idx1=which(rownames(tmp1)=="chr1:144450288:G:A")
idx2=which(rownames(tmp1)=="chr1:143209281:C:A")
cor(unlist(tmp1[idx1,]),unlist(tmp1[idx2,]),use="complete") #0.23

onco_sample=pheno_onco$Onc_ID[c(which(is.na(pheno_onco$Behaviour1)),which(pheno_onco$Behaviour1==1))]
onco_sample=intersect(onco_sample,fam$V3)
onco_sample1=paste0(onco_sample,"_",onco_sample)
write.table(onco_sample1,file="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/onco_sample.txt",row.names = F,col.names = F,quote=F)
cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1 --keep /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/onco_sample.txt --extract ../result/novelsnp.txt --make-pgen --out ../result/onco_novelsnp --memory 128000 --threads 8")
system(cmd)
cmd=paste0(plink2," --pfile ../result/onco_novelsnp --recode A-transpose --out ../result/onco_novelsnp")
onconovel=as.data.frame(fread("../result/onco_novelsnp.traw"))
rownames(onconovel)=onconovel$SNP
onconovel=onconovel[,-c(1:6)]
fam1=read.table("../result/onco_novelsnp.psam")
tmp=unlist(strsplit(fam1$V1,"_"))
fam1$V3=tmp[seq(1,length(tmp),2)]
idx=match(fam1$V3,pheno_onco$Onc_ID)
y=pheno_onco$Behaviour1[idx]
y[is.na(y)]=0
idx=which(rownames(onconovel)=="chr9:119662503:T:C")
fm=glm(y~unlist(onconovel[idx,]),family = binomial)
summary(fm)$coefficients

boxplot(unlist(onconovel[idx,])~y)
tmp=unlist(onconovel[idx,])
idx=y==0
hist(tmp[idx])
hist(tmp[!idx])
quantile(tmp[idx])
quantile(tmp[!idx])


write.table(independentnovelvar,file="../result/independentnovelvar.txt",row.names = F,quote=F)
tmp=unlist(strsplit(names(independentnovelvar),":"))
chr2=tmp[seq(1,length(tmp),4)]
pos2=tmp[seq(2,length(tmp),4)]


knownvar$foundminp=knownvar$foundpos=knownvar$founddist=knownvar$nearestdist=knownvar$nearesp=NA
for(i in 1:nrow(knownvar))
{
  tmp1=distance(gr_allpvalues,gr_knownvar[i])
  idx=which(tmp1==0)
  if (length(idx)>0)
  {
    idx1=which.min(allpvalues[idx])
    knownvar$foundminp[i]=allpvalues[idx[idx1]]
    knownvar$foundpos[i]=pos[idx[idx1]]
    knownvar$founddist[i]=abs(knownvar$foundpos[i]-knownvar$hg38position[i])
  }
  tmp1=distance(gr_allpvalues,gr_knownvar1[i])
  idx=which.min(tmp1)
  knownvar$nearesp[i]=allpvalues[idx]
  knownvar$nearestdist[i]=tmp1[idx]
}
write.table(knownvar,file ="../result/218_known_discovery_snp_paper_order_101323_overlap.txt",row.names = F,quote=F,sep="\t")
#knownvar=read.table("../result/218_known_discovery_snp_paper_order_101323_overlap.txt",header = T,sep="\t")
table(knownvar$foundminp<5e-8,useNA="ifany")
# FALSE  TRUE   
#   68   150  
table(knownvar$nearesp[knownvar$nearestdist==0]<5e-8,useNA="ifany")
# FALSE  TRUE 
# 137    81
quantile(knownvar$nearesp[knownvar$nearestdist==0])
table(knownvar$nearestdist==0)
# FALSE  TRUE 
# 19   199
quantile(knownvar$foundminp,na.rm=T)
idx=which(neardist@elementMetadata@listData$distance==0)
knownpvalues=allpvalues[idx]
sum(knownpvalues<5e-8) #12961



#position_nudge()
#qqplot,manhattan plot
load("../result/metascoreinfo4.RData")

library(data.table)
library("plotrix") #axis.break
library("RColorBrewer")
library("optparse")
library(readr)
library(dplyr)
library(ggplot2)
library(qqman)

source("/data/BB_Bioinformatics/Kevin/PRS_EASLC/code/theme_publication.R")
#draw QQ plot
qqplotdata <- function(logpvector){
  o = sort(logpvector,decreasing=T)
  e = -log10(ppoints(length(o)))
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]

  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA

  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)

  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}
#optbreak=1, break top p-values, used for very low p-value cases
plotqq=function(data,optbreak=1,title="")
{
  dat = data %>%
    mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>%
    select(rsid,CHR,BP,P,MAF) %>%
    rename(SNP = rsid)
  dat=dat[!is.na(dat$P),]
  x = dat$P
  z = qnorm(x / 2)
  lambda = round(median(z^2) / qchisq(0.5,1), 3)
  N.effect = median(data$N)
  lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)

  yLine <- c(-log10(5E-8))
  colLine <- c("red")
  dat$log10P = -log10(dat$P)
  gwas = as.data.frame(dat)
  # Determine frequency bins and create variable for binned QQ plot

  minMAF <- min(gwas$MAF)

  freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
  gwas$freqbin <- cut(gwas$MAF, freqbins,include.lowest=T)
  freqtable <- table(gwas$freqbin)
  freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
  freqtable <- freqtable[freqtable > 0]

  ## Generate QQ plot data by frequency bin
  fbin <- character(0)
  fN <- integer(0)
  fx <- numeric(0)
  fy <- numeric(0)
  fcol <- character(0)
  legendcol <- character(0)
  conf <- list()
  allcols <- brewer.pal(4,"Set1")
  ycol <- "log10P"
  for(f in 1:length(freqtable)){
    fbin <- c(fbin,names(freqtable)[f])
    fsnps <- which(gwas$freqbin ==names(freqtable)[f])
    plotdata <- qqplotdata(gwas[[ycol]][fsnps])
    fN <- c(fN,freqtable[f])
    fx <- c(fx,plotdata$e)
    fy <- c(fy,plotdata$o)
    fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
    conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                            'y'=c(plotdata$c975,rev(plotdata$c025)))
    legendcol <- c(legendcol,allcols[f])
  }
  legendtext <- paste0("MAF=",fbin,"; #SNPs=",format(fN,big.mark=",",scientific=FALSE))

  opt =  list(break.top = ifelse(optbreak==1,15,ceiling(max(fy))+1),
              top.size = 0.125)


  #png(filename = paste0(outpath,"/QQ_",eth[i1],"_",trait[i2],".png"), width = 8, height = 8, units = "in",res=300)
  xlim <- c(0,max(fx,na.rm=T))
  ylim <- c(0,max(fy,na.rm=T))
  maxY <- max(fy,na.rm=T)

  par(mar=c(5.1,5.1,4.1,1.1))

  lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
  lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
  #top
  lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
  if (any(lab2>max(lab1))) lab2 <- lab2[lab2 > max(lab1)]

  # resulting range of top scale in bottom scale units
  top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
  top.data = max(lab2)-opt$break.top

  # function to rescale the top part
  rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}

  plot(0,0,
       ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
       xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
       ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
       cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
       main=opt$maintitle,pch=19)

  # Plot confidence intervals
  for(p in 1:length(conf)){
    polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
            col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
            border = NA)
  }

  # add points below top
  points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)

  # identify line & add axis break
  lines(xlim,xlim,col="black",lty = 2)
  axis(1,cex.axis=1.5,cex.lab=1.5)
  par(las=1)
  axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
  box()
  if (optbreak==1)
  {
    rescaled.y = rescale(fy[fy>opt$break.top])
    par(las=0)
    points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
    par(las=1)
    axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
    axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
    axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
    lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
  }

  abline(h=ifelse(yLine<opt$break.top,
                  yLine,
                  rescale(yLine)),
         col=colLine,lwd=1.5,lty=2)
  legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
  text(5,1,expression(paste(lambda[1000]," = ")),cex = 1.5)
  text(5.9,1,paste(lambda_1000),cex = 1.5)

  title(title)
  #dev.off()

}

plotmanhattan=function(data,title="",filename="../result/man_ILCCOnever.png")
{
  dat = data %>%
    mutate(MAF = ifelse(FREQ_A1<=0.5,FREQ_A1,1-FREQ_A1)) %>%
    select(rsid,CHR,BP,P,MAF) %>%
    rename(SNP = rsid)
  dat$CHR1=factor(dat$CHR,levels = c(1:22,"X","Y"))
  idx=order(dat$CHR1,dat$BP)
  dat=dat[idx,]
  p.pwas <- 5E-08

  nCHR <- length(unique(dat$CHR))
  dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(dat$CHR)){
    nbp[i] <- max(dat[dat$CHR == i,]$BP)
    dat$BPcum[dat$CHR == i] <- dat$BP[dat$CHR == i] + s
    s <- s + nbp[i]
  }

  axis.set <- dat %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ylim <- abs(floor(log10(min(dat$P)))) + 2
  sig1 <- p.pwas


  sigline <- data.frame(sig=c(-log10(sig1)),val=c(paste0("P=",signif(sig1,2))))
  library(ggplot2)
  manhplot <- ggplot(dat, aes(x = BPcum, y = -log10(P),
                              color = as.factor(CHR), size = -log10(P))) +
    geom_point(alpha = 0.8, size=0.8) +
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    geom_hline(data = sigline, aes(yintercept = sig), color= "red", linetype="dashed") +
    guides(color = FALSE) +
    labs(x = NULL,
         y = "-log10(p)",
         linetype = "",
         title = title)+
    #subtitle = "A2: Critically ill COVID19+ vs. population controls;\nB1: Hospitalized COVID19+ vs non-hospitalized COVID19+;\nB2: Hospitalized COVID19+ vs. population controls;\nC2: Reported SARS-CoV-2 infection vs. population controls") +
    theme_Publication()+
    theme(
      legend.position = "top",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )

  #outpath <-"../result"


  ggsave(filename=filename,
         plot=manhplot, device="png",
         width=9, height=4, units="in", dpi=300)
}
#table(rownames(allpvalues0) %in% names(metascoreinfo4$freq))

myplot=function(myallpvalues=allpvalues,outprefix="waldtest_meta")
{
  idx0=match(names(myallpvalues),names(metascoreinfo4$freq))
  all(names(metascoreinfo4$freq[idx0])==names(myallpvalues)) #T

  tmp=unlist(strsplit(names(myallpvalues),":"))
  data=data.frame(SNP=names(myallpvalues),CHR=tmp[seq(1,length(tmp),4)],BP=tmp[seq(2,length(tmp),4)],EAF=metascoreinfo4$freq[idx0],P=myallpvalues,N=0)
  data$CHR=gsub("chr","",data$CHR)
  data$BP=as.integer(data$BP)
  data$CHR=as.integer(data$CHR)
  #euronco
  idx=which(metascoreinfo4$avail[idx0,1]==1)
  data$N[idx]=data$N[idx]+1/(1/34367+1/52500)
  #euricogs
  idx=which(metascoreinfo4$avail[idx0,2]==1)
  data$N[idx]=data$N[idx]+1/(1/13283+1/33505)
  #asianicogs
  idx=which(metascoreinfo4$avail[idx0,3]==1)
  data$N[idx]=data$N[idx]+1/(1/1682+1/5940)
  idx=which(metascoreinfo4$avail[idx0,4]==1)
  data$N[idx]=data$N[idx]+1/(1/9064+1/12344)
  idx=which(metascoreinfo4$avail[idx0,5]==1)
  data$N[idx]=data$N[idx]+1/(1/1815+1/2088)
  data$EAF[is.na(data$EAF)]=0.01

  data1 = data %>%
    mutate(CHR = as.integer(CHR),
           BP = as.integer(BP),
           FREQ_A1 = as.numeric(EAF),
           P = as.numeric(P),
           rsid = SNP,
           N = as.integer(N)) %>%
    mutate(P = ifelse(P==0,1E-300,P)) %>%
    select(rsid,CHR,BP,FREQ_A1,P,N)
  data1=data1[data1$FREQ_A1>0.01 & data1$FREQ_A1<0.99,]
  png(filename = paste0("../result/QQplot_",outprefix,".png"), width = 8, height = 8, units = "in",res=300)
  plotqq(data=data1,optbreak=1,title="")
  dev.off()
  plotmanhattan(data=data1,title="",filename=paste0("../result/man_",outprefix,".png"))
}

myplot(myallpvalues=allpvalues,outprefix="ACAT_combined")
tmp=intrinsicpvalues$P
names(tmp)=rownames(intrinsicpvalues)
myplot(myallpvalues=tmp,outprefix="Intrinsic_meta")

myplot(myallpvalues=novelpvalues,outprefix="Novel_variants")
myplot(myallpvalues=knownpvalues,outprefix="Known_variants")
