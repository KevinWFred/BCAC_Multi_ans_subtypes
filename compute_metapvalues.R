#!/usr/bin/env Rscript
#get the p-values of global test
#aggregate logodds and sigma across populations
#get the final global test p-values (ACAT p-values)
#work on the QCed SNPs (HWE and info score)
.libPaths(c("/data/wangx53",.libPaths()))
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

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
#QCed SNP lists:
qc_african_onco=as.data.frame(fread("../result/imp_QC/onco/african_info.snplist",header=F))
qc_asian_onco=as.data.frame(fread("../result/imp_QC/onco/asian_info.snplist",header=F))
qc_euro_onco=as.data.frame(fread("../result/imp_QC/onco/euro_info.snplist",header=F))
qc_hispanic_onco=as.data.frame(fread("../result/imp_QC/onco/hispanic_info.snplist",header=F))
qc_african_icogs=as.data.frame(fread("../result/imp_QC/icogs/african_info.snplist",header=F))
qc_asian_icogs=as.data.frame(fread("../result/imp_QC/icogs/asian_info.snplist",header=F))
qc_euro_icogs=as.data.frame(fread("../result/imp_QC/icogs/euro_info.snplist",header=F))

#for euro
pop="euro"
load(paste0("../result/",pop,"_beta_sigma.RData"))
euro_allsigma_icogs=allsigma_icogs
euro_allsigma_icogs=euro_allsigma_icogs[names(euro_allsigma_icogs) %in% qc_euro_icogs[,1]]
euro_allsigma_onco=allsigma_onco
euro_allsigma_onco=euro_allsigma_onco[names(euro_allsigma_onco) %in% qc_euro_onco[,1]]

euro_allbeta_icogs=beta_icogs
euro_allbeta_icogs=euro_allbeta_icogs[rownames(euro_allbeta_icogs) %in% qc_euro_icogs[,1],]
all(rownames(euro_allbeta_icogs)==names(euro_allsigma_icogs))
euro_allbeta_onco=beta_onco
euro_allbeta_onco=euro_allbeta_onco[rownames(euro_allbeta_onco) %in% qc_euro_onco[,1],]
all(rownames(euro_allbeta_onco)==names(euro_allsigma_onco))

pop="asian"
load(paste0("../result/",pop,"_beta_sigma.RData"))
asian_allsigma_icogs=allsigma_icogs
asian_allsigma_icogs=asian_allsigma_icogs[names(asian_allsigma_icogs) %in% qc_asian_icogs[,1]]
asian_allsigma_onco=allsigma_onco
asian_allsigma_onco=asian_allsigma_onco[names(asian_allsigma_onco) %in% qc_asian_onco[,1]]
asian_allbeta_icogs=beta_icogs
asian_allbeta_icogs=asian_allbeta_icogs[rownames(asian_allbeta_icogs) %in% qc_asian_icogs[,1],]
all(rownames(asian_allbeta_icogs)==names(asian_allsigma_icogs))
asian_allbeta_onco=beta_onco
asian_allbeta_onco=asian_allbeta_onco[rownames(asian_allbeta_onco) %in% qc_asian_onco[,1],]
all(rownames(asian_allbeta_onco)==names(asian_allsigma_onco))
#for african
pop="african"
#save(beta_onco,sigma_onco,allsigma_onco,file=paste0(outprefix,"_beta_sigma.RData"))
load(paste0("../result/",pop,"_beta_sigma.RData"))
african_allsigma_onco=allsigma_onco
african_allsigma_onco=african_allsigma_onco[names(african_allsigma_onco) %in% qc_african_onco[,1]]
african_allbeta_onco=beta_onco
african_allbeta_onco=african_allbeta_onco[rownames(african_allbeta_onco) %in% qc_african_onco[,1],]
all(rownames(african_allbeta_onco)==names(african_allsigma_onco))

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
save(allres,file=paste0("../result/metaresult/all/","QCres_",i1,".RData"))
print(Sys.time())
print("done")

# swarmjobs=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/compute_metapvalues_QC.R",1000),
#                      i1=1:1000)
# write.table(swarmjobs,file=paste0("compute_metapvalues_QC.swarm"),row.names=F,col.names=F,sep="\t",quote=F)
#17215969
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/compute_metapvalues_QC.swarm -g 64 --module R/4.3 --partition=quick --time=4:00:00 --gres=lscratch:40 -p 2

#collect the results
allpvalues=NULL
for (i in 1:1000)
{
  if (i %% 100 ==0) cat(i,'..')
  load(paste0("../result/metaresult/all","/QCres_",i,".RData"))
  allpvalues=rbind(allpvalues,allres)
}
intrinsicpvalues=allpvalues

load("../result/scoretestresult/QCallpvalues.Rdata")
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
save(intrinsicpvalues,scoretestallpvaues,allpvalues,file="../result/compute_metapvalues_newQC.RData")
quantile(allpvalues)
# 0%           25%           50%           75%          100%
# 3.304934e-301  2.282329e-01  4.891243e-01  7.520828e-01  1.000000e+00
load("../result/metascoreinfo4_newQC.RData")
table(names(allpvalues) %in% rownames(metascoreinfo4$score)) #T
idx=match(names(allpvalues),rownames(metascoreinfo4$score))
# metascoreinfo=metascoreinfo4
# metascoreinfo$score=metascoreinfo4$score[idx,]
# metascoreinfo$infor=metascoreinfo4$infor[idx,]
# metascoreinfo$avail=metascoreinfo4$avail[idx,]
# metascoreinfo$allfreq=metascoreinfo4$allreq[idx,]
# metascoreinfo$freq=metascoreinfo4$freq[idx]
# metascoreinfo$maxfreq=metascoreinfo4$maxfreq[idx]

sum(allpvalues< 5e-8) #12289
#sum(allpvalues[metascoreinfo$freq>0.01]< 5e-8) #11841
allpvalues1=allpvalues[names(allpvalues) %in% names(metascoreinfo4$freq[metascoreinfo4$freq>0.01])] #16454816

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
gr_allpvalues=GRanges(seqnames = chr,ranges = IRanges(start=pos,width = 1))
gr_knownvar=GRanges(seqnames = knownvar$CHR,ranges=IRanges(start=pmax(0,knownvar$hg38position-5e5),width=1e6))
sum(width(reduce(gr_knownvar))) #196005593
gr_knownvar1=GRanges(seqnames = knownvar$CHR,ranges=IRanges(start=knownvar$hg38position,width=1))
#dist=distance(gr_allpvalues,gr_knownvar)
neardist=distanceToNearest(gr_allpvalues,gr_knownvar)
neardist1=distanceToNearest(gr_allpvalues,gr_knownvar1)
sum(neardist@elementMetadata@listData$distance==0) #1495557 within 500KB of known
sum(allpvalues<5e-8) #12289,11841
idx=which(neardist@elementMetadata@listData$distance!=0)
novelpvalues=allpvalues[idx]
sum(novelpvalues<5e-8) #257
#with freq>0.01
novelpvalues1=novelpvalues[names(novelpvalues) %in% names(allpvalues1)]
sum(novelpvalues1<5e-8) #197

#compute LD between novels and known
ge_knownsnps=function(pvarfile="../result/imp_onco/merged1.pvar",outprefix="../result/onco_218knownvar")
{
  pvar=as.data.frame(fread(pvarfile))
  gr_pvar=GRanges(seqnames=pvar$`#CHROM`,ranges=IRanges(start=pvar$POS,width=1))
  tmp=distanceToNearest(gr_knownvar1,gr_pvar)
  selsnps=pvar$ID[tmp@to]
  write.table(selsnps,file=paste0(outprefix,".snplist"),row.names = F,col.names=F,quote=F)
}

filter_novelsig_ld=function(novelsig,knownsnpfile="../result/onco_218knownvar.snplist",
                            allsnpfile="../result/onco_novel_knkownsnps.bim",ldfile="../result/test.ld")
{
  ldres=as.data.frame(fread(ldfile))
  allsnps=read.table(allsnpfile)
  rownames(ldres)=colnames(ldres)=allsnps$V2
  idx1=which(colnames(ldres) %in% novelsig)
  knownsnps=read.table(knownsnpfile)
  idx2=which(colnames(ldres) %in% knownsnps$V1)
  tmp=ldres[idx1,idx2]
  tmp1=apply(tmp,1,function(x) {
    max(x,na.rm=T)>0.1
  })
  tmp2=apply(tmp,1,function(x) {
    max(x,na.rm=T)
  })
  novelsig2rm=rownames(tmp)[tmp1==T]
  return(novelsig2rm)
}
get_novel_snps=function(mynovelpvalues=novelpvalues)
{
  idx=which(mynovelpvalues<5e-8)
  novelsig=names(mynovelpvalues)[idx]
  print(paste0("number of var with p < 5e-8:",length(novelsig))) #257
  #work on onco
  ge_knownsnps()
  onco_knownsnps=read.table("../result/onco_218knownvar.snplist")
  onco_novel_knownsnps=c(onco_knownsnps$V1,novelsig)
  write.table(onco_novel_knownsnps,file="../result/onco_novel_knkownsnps.snplist",row.names = F,col.names=F,quote=F)
  cmd=paste0(plink2, " --pfile ../result/imp_onco/merged1 --extract ../result/onco_novel_knkownsnps.snplist --make-bed --out ../result/onco_novel_knkownsnps")
  system(cmd)           
  cmd=paste0(plink, " --bfile ../result/onco_novel_knkownsnps --r2 square  --out ../result/onco_novel_knkownsnps")
  system(cmd)
  onco_novelsig2rm=filter_novelsig_ld(novelsig,knownsnpfile="../result/onco_218knownvar.snplist",
                              allsnpfile="../result/onco_novel_knkownsnps.bim",ldfile="../result/onco_novel_knkownsnps.ld")
  #work on icogs
  ge_knownsnps(pvarfile="../result/imp_icogs/merged1.pvar",outprefix="../result/icogs_218knownvar")
  icogs_knownsnps=read.table("../result/icogs_218knownvar.snplist")
  icogs_novel_knownsnps=c(icogs_knownsnps$V1,novelsig)
  write.table(icogs_novel_knownsnps,file="../result/icogs_novel_knkownsnps.snplist",row.names = F,col.names=F,quote=F)
  cmd=paste0(plink2, " --pfile ../result/imp_icogs/merged1 --extract ../result/icogs_novel_knkownsnps.snplist --make-bed --out ../result/icogs_novel_knkownsnps")
  system(cmd)           
  cmd=paste0(plink, " --bfile ../result/icogs_novel_knkownsnps --r2 square  --out ../result/icogs_novel_knkownsnps")
  system(cmd)
  icogs_novelsig2rm=filter_novelsig_ld(novelsig,knownsnpfile="../result/icogs_218knownvar.snplist",
                                      allsnpfile="../result/icogs_novel_knkownsnps.bim",ldfile="../result/icogs_novel_knkownsnps.ld")
  print(table(icogs_novelsig2rm %in% onco_novelsig2rm))
  novelsig2rm=unique(c(icogs_novelsig2rm,onco_novelsig2rm))
  novelsig1=novelsig[!novelsig %in% novelsig2rm]
  print(paste0("number of var with ld<0.1:",length(novelsig1))) #201 
  
  #find independent novel regions
  idx=match(novelsig1,names(mynovelpvalues))
  novelsig1pvalue=mynovelpvalues[idx]
  tmp=unlist(strsplit(novelsig1,":"))
  chr1=tmp[seq(1,length(tmp),4)]
  chr1=gsub("chr","",chr1)
  pos1=as.integer(tmp[seq(2,length(tmp),4)])
  gr_novelsig1=GRanges(seqnames = chr1,ranges=IRanges(start=pos1,width=1))
  novelsig2=NULL
  n=1
  while (min(novelsig1pvalue)<5e-8 & n<length(novelsig1pvalue))
  {
    idx=which.min(novelsig1pvalue)
    if (novelsig1pvalue[idx]<5e-8) novelsig2=c(novelsig2,names(novelsig1pvalue[idx]))
    novelsig1pvalue[idx]=1
    gr_snp=GRanges(seqnames = chr1[idx],ranges=IRanges(start=pos1[idx],width=1))
    dist1=distance(gr_novelsig1,gr_snp)
    idx1=which(dist1<1e6)
    if(length(idx1)>0) novelsig1pvalue[c(idx1)]=1
    n=n+1
  }
  print(paste0("number of var with dist>1MB:",length(novelsig2))) #22
  return(list(novelsig=novelsig,novelsig1=novelsig1,novelsig2=novelsig2))
}

allnovelsig=get_novel_snps(mynovelpvalues=novelpvalues)#257,201,22
#freq>0.01
allnovelsig1=get_novel_snps(mynovelpvalues=novelpvalues1) #197,144,13
table(allnovelsig1$novelsig2 %in% allnovelsig$novelsig2)
# FALSE  TRUE 
# 3    10
table(allnovelsig1$novelsig1 %in% allnovelsig$novelsig1) #TRUE 144

# tmp=data.frame(snp=names(independentnovelvar),p=independentnovelvar)
# write.table(tmp,file="../result/independentnovelvarpvalues.txt",row.names = F,sep="\t",quote=F)
# tmp=data.frame(snp=names(independentnovelvar))
# write.table(tmp,file="../result/independentnovelvarsnp.txt",row.names = F,sep="\t",quote=F)
# 
# #LD pruning
# pheno_onco=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_onco_pheno_v15_02_corrected.txt",header=T,sep="\t")
# fam=read.table("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1.psam")
# tmp=unlist(strsplit(fam$V1,"_"))
# fam$V3=tmp[seq(1,length(tmp),2)]
# idx=match(fam$V3,pheno_onco$Onc_ID)
# table(pheno_onco$EthnicityGeno[idx])
# onco_contrl=pheno_onco$Onc_ID[which(is.na(pheno_onco$Behaviour1))]
# onco_contrl=intersect(onco_contrl,fam$V3)
# onco_contrl1=paste0(onco_contrl,"_",onco_contrl)
# #write.table(onco_contrl1,file="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/onco_control.txt",row.names = F,col.names = F,quote=F)
# cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1 --keep /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/onco_control.txt --extract ../result/independentnovelvarsnp.txt --make-bed --out ../result/QConco_control_independentnovelvar --memory 64000 --threads 8")
# system(cmd)
# 
# #tmp=as.data.frame(fread("../result/QConco_control_independentnovelvar.bim"))
# onco_pvar=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1.pvar"))
# table(names(independentnovelvar) %in% onco_pvar$ID)
# # FALSE  TRUE
# # 6    17
# oncomissing=names(independentnovelvar)[!names(independentnovelvar) %in% onco_pvar$ID]
# 
# pthr=1
# r2thr=0.1
# kbpthr=500
# 
# # cmd=paste0(plink," --bfile ../result/onco_control_independentnovelsnp"," --clump ../result/novelpvalues.txt"," --clump-p1 ",
# #            pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field snp --clump-field p --out ../result/onco_independentnovel")
# cmd=paste0(plink," --bfile ../result/QConco_control_independentnovelvar"," --clump ../result/independentnovelvarpvalues.txt"," --clump-p1 ",
#            pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field snp --clump-field p --out ../result/QConco_control_independentnovelvar")
# 
# system(cmd)
# #17
# onco_clumping=read.table("../result/QConco_control_independentnovelvar.clumped",header=T)
# 
# icogs_pvar=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1.pvar"))
# table(names(independentnovelvar) %in% icogs_pvar$ID)
# # FALSE  TRUE 
# # 1    22
# table(oncomissing %in% icogs_pvar$ID) #T
# pheno_icogs=read.table("/data/BB_Bioinformatics/ProjectData/BCAC/phenotype/concept_750_zhang_icogs_pheno_v15_02.txt",header=T,sep="\t")
# fam_icogs=read.table("/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1.psam")
# tmp=unlist(strsplit(fam_icogs$V1,"_"))
# fam_icogs$V3=tmp[seq(1,length(tmp),2)]
# idx=match(fam_icogs$V3,pheno_icogs$SG_ID)
# table(pheno_icogs$EthnicityGeno[idx])
# icogs_contrl=pheno_icogs$SG_ID[which(is.na(pheno_icogs$Behaviour1))]
# icogs_contrl=intersect(icogs_contrl,fam_icogs$V3)
# icogs_contrl1=paste0(icogs_contrl,"_",icogs_contrl)
# write.table(icogs_contrl1,file="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/icogs_control.txt",row.names = F,col.names = F,quote=F)
# cmd=paste0(plink2," --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1 --keep /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/icogs_control.txt --extract ../result/QCnovelsnp.txt --make-bed --out ../result/QCicogs_control_novelpvalues1 --memory 64000 --threads 8")
# system(cmd)
# # cmd=paste0(plink," --bfile ../result/icogs_control_independentnovelsnp --recode A-transpose --out ../result/icogs_control_independentnovelsnp")
# # cmd=paste0(plink," --bfile ../result/icogs_control_independentnovelsnp --freq --out ../result/icogs_control_independentnovelsnp")
# 
# # cmd=paste0(plink," --bfile ../result/icogs_control_independentnovelsnp"," --clump ../result/novelpvalues.txt"," --clump-p1 ",
# #            pthr," --clump-r2 ",r2thr," --clump-kb ",kbpthr," --clump-snp-field snp --clump-field p --out ../result/icogs_independentnovel")
# cmd=paste0(plink," --bfile ../result/QCicogs_control_novelpvalues1"," --clump ../result/QCnovelpvalues.txt"," --clump-p1 ",
#            pthr," --clump-r2 ",r2thr," --clump-snp-field snp --clump-field p --out ../result/QCicogs_novelpvalues1")
# system(cmd)
# #28
# icogs_clumping=read.table("../result/QCicogs_novelpvalues1.clumped",header=T)
# table(onco_clumping$SNP %in% icogs_clumping$SNP)
# # TRUE 
# # 24
# novelpvalues2=novelpvalues1[names(novelpvalues1) %in% icogs_clumping$SNP]
# # idx=match(names(novelpvalues2),names(allpvalues))
# # chr1=chr[idx]
# # pos1=pos[idx]
# # gr_novel1=GRanges(seqnames = chr1,ranges=IRanges(start=pos1,width=1))
# # independentnovelvar=NULL
# # n=1
# # while (min(novelpvalues2)<5e-8 & n<length(novelpvalues2))
# # {
# #   idx=which.min(novelpvalues2)
# #   if (novelpvalues2[idx]<5e-8) independentnovelvar=c(independentnovelvar,novelpvalues2[idx])
# #   novelpvalues2[idx]=1
# #   gr_snp=GRanges(seqnames = chr1[idx],ranges=IRanges(start=pos1[idx],width=1))
# #   dist1=distance(gr_novel1,gr_snp)
# #   idx1=which(dist1<1e6)
# #   if(length(idx1)>0) novelpvalues2[c(idx1)]=1
# #   n=n+1
# # }
# # length(independentnovelvar) #14
# table(names(independentnovelvar) %in% onco_clumping$SNP)
# # FALSE  TRUE 
# # 1    13 
# tmp=data.frame(snp=names(independentnovelvar),p=independentnovelvar)
# 
# tmp=data.frame(snp=names(independentnovelvar))
# write.table(tmp,file="../result/QCindependentnovelsnp.txt",row.names = F,sep="\t",quote=F)

annotate_novelsig=function(allnovelsnps=allnovelsig$novelsig2,outprefix="QCallnovel")
{
  idx=match(allnovelsnps,names(allpvalues))
  allnovelsnps=data.frame(ID=allnovelsnps,chr=NA,pos=NA,rsid=NA,pvalue=allpvalues[idx])
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
  write.csv(allnovelsnps,file=paste0("../result/",outprefix,"snps.csv"),row.names = F,quote=F)
  gr_allnovelsnps=GRanges(seqnames = allnovelsnps$chr,ranges = IRanges(start=allnovelsnps$pos,width=1))
  tmp=distanceToNearest(gr_allnovelsnps,gr_knownvar1)
  allnovelsnps$dist2nearestknown=tmp@elementMetadata@listData$distance
  all(knownvar$hg38position==start(gr_knownvar1))
  print(table(allnovelsnps$dist2nearestknown<2e6))
  # FALSE  TRUE
  # 13    9
  quantile(allnovelsnps$dist2nearestknown)
  # 0%      25%      50%      75%     100%
  # 550583  1008492  2731904  6674517 28452496
  allnovelsnps0=allnovelsnps[allnovelsnps$dist2nearestknown>2e6,]
  write.table(allnovelsnps0,file=paste0("../result/",outprefix,"snps0.txt"),row.names = F,sep="\t",quote=F)
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
  write.table(allnovelsnps1,file=paste0("../result/",outprefix,"snps1.txt"),row.names = F,sep="\t",quote=F)
  #
  
}
annotate_novelsig(allnovelsnps=allnovelsig$novelsig2,outprefix="QCallnovel")
annotate_novelsig(allnovelsnps=allnovelsig1$novelsig2,outprefix="QCfreq01allnovel")

#
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
write.table(knownvar,file ="../result/QC218_known_discovery_snp_paper_order_101323_overlap.txt",row.names = F,quote=F,sep="\t")
#knownvar=read.table("../result/QC218_known_discovery_snp_paper_order_101323_overlap.txt",header = T,sep="\t")
table(knownvar$foundminp<5e-8,useNA="ifany")
# FALSE  TRUE
#   82   136
table(knownvar$nearesp[knownvar$nearestdist==0]<5e-8,useNA="ifany")
# FALSE  TRUE
# 74    100
# quantile(knownvar$nearesp[knownvar$nearestdist==0])
table(knownvar$nearestdist==0)
# FALSE  TRUE
# 44   174
# quantile(knownvar$foundminp,na.rm=T)
idx=which(neardist@elementMetadata@listData$distance==0)
knownpvalues=allpvalues[idx]
sum(knownpvalues<5e-8) #11644



#position_nudge()
#qqplot,manhattan plot
load("../result/metascoreinfo4_newQC.RData")

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
  ggsave(filename=filename,
         plot=manhplot, device="png",
         width=12, height=4, units="in", dpi=300)
  #outpath <-"../result"
  myknownvar=knownvar
  myallnovelsnps=allnovelsnps
  if (!is.null(myknownvar))
  {
    gr_myknownvar=GRanges(seqnames=myknownvar$CHR,ranges = IRanges(myknownvar$hg38position,width=1))
    gr_dat=GRanges(seqnames = dat$CHR,ranges = IRanges(start=dat$BP,width = 1))
    tmp=distanceToNearest(gr_myknownvar,gr_dat)
    manhplot1=manhplot+
      geom_vline(data = dat[tmp@to,], aes(xintercept = BPcum), color= "green", linetype="solid",size=0.3,alpha=0.3)
    filename1=gsub(".png","",filename)
    filename1=paste0(filename1,"_knownvar.png")
    ggsave(filename=filename1,
           plot=manhplot1, device="png",
           width=12, height=4, units="in", dpi=300)
    gr_myallnovelsnps=GRanges(seqnames = myallnovelsnps$chr,ranges = IRanges(start=myallnovelsnps$pos,width=1))
    tmp=distanceToNearest(gr_myallnovelsnps,gr_dat)
    manhplot2=manhplot+
      geom_vline(data = dat[tmp@to,], aes(xintercept = BPcum), color= "red", linetype="solid",size=0.3,alpha=0.3)
    filename1=gsub(".png","",filename)
    filename2=paste0(filename1,"_novelvar.png")
    ggsave(filename=filename2,
           plot=manhplot2, device="png",
           width=12, height=4, units="in", dpi=300)
  }

  
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

myplot(myallpvalues=allpvalues,outprefix="QC_ACAT_combined")
tmp=intrinsicpvalues$P
names(tmp)=rownames(intrinsicpvalues)
myplot(myallpvalues=tmp,outprefix="QC_Intrinsic_meta")

myplot(myallpvalues=novelpvalues,outprefix="QC_Novel_variants")
myplot(myallpvalues=knownpvalues,outprefix="QC_Known_variants")

#check frequnces of SNPs in each data
check_freq=function(snps=c("chr1:144450288:G:A","chr5:53186485:T:C"))
{
  res=data.frame(matrix(NA,nrow=length(snps),ncol=5))
  rownames(res)=snps
  colnames(res)=c("icogs_euro","icogs_asian","onco_euro","onco_asian","onco_african")
  tmp=as.data.frame(fread("../result/imp_icogs/euro/euro.afreq"))
  idx=match(snps,tmp$ID)
  res$icogs_euro=tmp$ALT_FREQS[idx]
  tmp=as.data.frame(fread("../result/imp_icogs/asian/asian.afreq"))
  idx=match(snps,tmp$ID)
  res$icogs_asian=tmp$ALT_FREQS[idx]
  tmp=as.data.frame(fread("../result/imp_onco/euro/euro.afreq"))
  idx=match(snps,tmp$ID)
  res$onco_euro=tmp$ALT_FREQS[idx]
  tmp=as.data.frame(fread("../result/imp_onco/asian/asian.afreq"))
  idx=match(snps,tmp$ID)
  res$onco_asian=tmp$ALT_FREQS[idx]
  tmp=as.data.frame(fread("../result/imp_onco/african/african.afreq"))
  idx=match(snps,tmp$ID)
  res$onco_african=tmp$ALT_FREQS[idx]
}