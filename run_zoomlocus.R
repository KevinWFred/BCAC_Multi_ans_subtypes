#!/usr/bin/env Rscript
#to run standalone zoomlocus
#source /data/BB_Bioinformatics/Kevin/tools/python2/bin/activate
#export PATH=:/data/BB_Bioinformatics/Kevin/tools/generic-new_fugue/bin:/data/BB_Bioinformatics/Kevin/tools/locuszoom/bin:$PATH
#ml samtools
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")

.libPaths(c("/data/wangx53",.libPaths()))
system("source '/data/BB_Bioinformatics/Kevin/tools/python2/bin/activate'")
tmp=Sys.getenv("PATH")
tmp1=paste0("/data/BB_Bioinformatics/Kevin/tools/generic-new_fugue/bin:/data/BB_Bioinformatics/Kevin/tools/locuszoom/bin:/data/BB_Bioinformatics/Kevin/tools/python2/bin:")
Sys.setenv(PATH=paste0(tmp1,tmp))
system("ml samtools")
print(Sys.getenv("PATH"))

plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

library(GenomicRanges)
library(data.table)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
seqlevelsStyle(genome) <- "NCBI"
chain=import.chain("../../tools/liftover/hg38ToHg19.over.chain")
asian_icogs_pvar=as.data.frame(fread("../result/imp_icogs/asian/asian.pvar"))
asian_onco_pvar=as.data.frame(fread("../result/imp_onco/asian/asian.pvar"))
euro_icogs_pvar=as.data.frame(fread("../result/imp_icogs/euro/euro.pvar"))
euro_onco_pvar=as.data.frame(fread("../result/imp_onco/euro/euro.pvar"))
african_onco_pvar=as.data.frame(fread("../result/imp_onco/african/african.pvar"))

load("../result/compute_metapvalues_population_QC.RData")
#load("../result/metascoreinfo4_newQC.RData")

update_bimhg19pos=function(bimfile="../result/zoomlocus/Snpchr12_81015799_G_T_asian.bim",
                           hg19posfile=hg19posfile)
{
  tmp=read.table(bimfile)
  hg19pos=read.table(hg19posfile)
  idx=match(tmp$V2,hg19pos$V1)
  tmp$V4=hg19pos$V2[idx]
  tmp=tmp[!is.na(tmp$V4),]
  write.table(tmp,file=bimfile,row.names = F,col.names = F,sep=" ",quote=F)
}

#allnovelcandidatesnps=read.csv("../result/QCallcandidatenovelsnps.csv")
get_zoomlocus=function(snpfile="../result/QCfreq01allcandidatenovelsnps.csv",
                       i1=1,outprefix="QCfreq01")
{
  allnovelcandidatesnps=read.csv(snpfile)
  snp=allnovelcandidatesnps$ID[i1] #"chr12:81015799:G:T"
  snpid=gsub(":","_",snp)
  tmp=unlist(strsplit(snp,":"))
  snpchr=tmp[1] #chr12
  rsid=allnovelcandidatesnps$rsid[i1]
  gr_snp=GRanges(seqnames = gsub("chr","",tmp[1]),ranges = IRanges(start=as.integer(tmp[2]),width = 1))
  idx=which.min(allnovelcandidatesnps[i1,c("P_EUR","P_EAS","P_AFR")])
  if (idx==1) pop="euro"
  if (idx==2) pop="asian"
  if (idx==3) pop="african"
  
  outprefix1=paste0("../result/zoomlocus/",outprefix,"_Snp",snpid,"_",pop)
  #myallpvalues=get(paste0("collect_pvalues_",pop))
  #myallpvalues=myallpvalues$acatpvalues
  myallpvalues=allpvalues
  tmp=unlist(strsplit(names(myallpvalues),":"))
  mychr=tmp[seq(1,length(tmp),4)]
  mychr=gsub("chr","",mychr)
  mypos=as.integer(tmp[seq(2,length(tmp),4)])
  gr_myallpvalues=GRanges(seqnames = mychr,ranges = IRanges(start=mypos,width=1))
  mydist=distance(gr_myallpvalues,gr_snp)
  idx=which(mydist<5e6)
  myallpvalues1=myallpvalues[idx]
  tmp=unlist(strsplit(names(myallpvalues1),":"))
  mychr1=tmp[seq(1,length(tmp),4)]
  mypos1=as.integer(tmp[seq(2,length(tmp),4)])
  metadat=data.frame(snpid=names(myallpvalues1),MarkerName=paste0(mychr1,":",mypos1),Pvalue=myallpvalues1,hg38pos=mypos1)
  gr_myallpvalues1=GRanges(seqnames = gsub("chr","",mychr1),ranges = IRanges(start=mypos1,width=1))
  my_snps <- snpsByOverlaps(all_snps, gr_myallpvalues1, genome = genome)
  my_snps = as.data.frame(my_snps)
  idx = match(paste0(gsub("chr","",mychr1),":",mypos1),paste0(my_snps$seqnames,":",my_snps$pos))
  metadat$MarkerName=my_snps$RefSNP_id[idx]
  # idx=which(is.na(metadat$MarkerName))
  # tmp=unlist(strsplit(metadat$snpid[idx],":"))
  # gr_metadat=GRanges(seqnames = tmp[seq(1,length(tmp),4)],ranges = IRanges(as.integer(tmp[seq(2,length(tmp),4)]),width=1))
  # tmp=as.data.frame(liftOver(gr_metadat,chain))
  # tmp1=paste0(snpchr,":",mypos1)
  # tmp1[tmp$group]=paste0(snpchr,":",tmp$start)
  # metadat$MarkerName[idx]=tmp1[idx]
  gr_allmetadat=GRanges(seqnames = mychr1,ranges = IRanges(start=mypos1,width=1))
  tmp=as.data.frame(liftOver(gr_allmetadat,chain))
  tmp1=as.integer(mypos1)
  tmp1[tmp$group]=tmp$start
  metadat$hg19pos=tmp1
  metadat=metadat[grepl("^rs",metadat$MarkerName),]
  metalfile=paste0(outprefix1,"_metalfile.txt")
  write.table(metadat[,c(2,3)],file=metalfile,col.names=T,row.names = F,sep="\t",quote=F)
  snpmapfile=paste0(outprefix1,"_snpmap.txt")
  write.table(metadat[,c(2,1)],file=snpmapfile,row.names = F,col.names = F,sep="\t",quote=F)
  snpfile=paste0(outprefix1,"_snp.txt")
  write.table(metadat$snpid,file=snpfile,row.names = F,col.names = F,quote=F)
  hg19posfile=paste0(outprefix1,"_hg19pos.txt")
  write.table(metadat[,c("MarkerName","hg19pos")],file=hg19posfile,row.names = F,col.names = F,quote=F)
  if (snp %in% get(paste0(pop,"_onco_pvar"))$ID)
  {
    cmd=paste0(plink2," --pfile ../result/imp_onco/",pop,"/",pop," --extract ",snpfile," --make-bed --out ",outprefix1," --memory 64000 --threads 8")
    system(cmd)
  }else
  {
    cmd=paste0(plink2," --pfile ../result/imp_icogs/",pop,"/",pop," --extract ",snpfile," --make-bed --out ",outprefix1," --memory 64000 --threads 8")
    system(cmd)
  }
  
  cmd=paste0(plink," --bfile ",outprefix1," --update-name ",snpmapfile," 1 2 --make-bed --out ",outprefix1," --memory 64000 --threads 8")
  system(cmd)
  update_bimhg19pos(bimfile=paste0(outprefix1,".bim"),
                    hg19posfile=hg19posfile)
  
  #tmp1=read.table(paste0(outprefix1,".bim"))
  cmd=paste0(plink2," --bfile ",outprefix1," --sort-vars "," --make-pgen --out ",outprefix1," --memory 64000 --threads 8")
  system(cmd)
  #tmp1=read.table(paste0(outprefix1,".pvar"))
  cmd=paste0(plink2," --pfile ",outprefix1," --recode vcf --out ",outprefix1," --memory 64000 --threads 8")
  system(cmd)
  tmp=as.data.frame(fread(paste0(outprefix1,".vcf")))
  #tmp=read.table(paste0(outprefix1,".vcf"))
  tmp$FILTER="PASS"
  #tmp[tmp=="./."]="0/0"
  #tmp$V1=paste0("chr",tmp$V1)
  fwrite(tmp,sep="\t",file=paste0(outprefix1,"_pass",".vcf"),row.names = F,col.names = F,quote=F)
  cmd=paste0("head -7 ",outprefix1,".vcf > ", outprefix1,"1.vcf")
  system(cmd)
  cmd=paste0("cat ",outprefix1,"_pass",".vcf >>",outprefix1,"1.vcf")
  system(cmd)
  cmd=paste0("bgzip -c ",outprefix1,"1.vcf > ",outprefix1,"1.vcf.gz")
  system(cmd)
  
  cmd=paste0("tabix -f -p vcf ",outprefix1,"1.vcf.gz")
  system(cmd)
  vcf_gzfile=paste0(outprefix1,"1.vcf.gz")
  #cmd=paste0("ml locuszoom; locuszoom"," --metal ",metalfile," --flank 1000kb"," --build hg19 --pop ASN --source 1000G_March2012 --refsnp rs138955280 --pvalcol Pvalue --plotonly --prefix zoomlocustest --ld-vcf ",vcf_gzfile," ymax=",20," signifLine=",-log10(5*1e-8)," signifLineColor=red ", " axisTextSize=1.4 legendSize=1.0  axisSize=1.4 xlabPos=-2.9")
  cmd=paste0("locuszoom"," --metal ",metalfile," --flank 500kb"," --build hg19 --pop ASN --source 1000G_March2012 --refsnp ",rsid," --pvalcol Pvalue --plotonly --prefix ",outprefix1," --ld-vcf ",vcf_gzfile," ymax=",20," signifLine=",-log10(5*1e-8)," signifLineColor=red ", " axisTextSize=1.4 legendSize=1.0  axisSize=1.4 xlabPos=-2.9")
  system(cmd)
}
args <- commandArgs(trailingOnly=T)
snpfile=as.character(args[1])
i1=as.integer(args[2])
outprefix=as.character(args[3])
print(args)
get_zoomlocus(snpfile=snpfile,
                       i1=i1,outprefix=outprefix)
  
tmp1=read.csv("../result/QCfreq01allcandidatenovelsnps.csv")
tmp=data.frame(code="/data/BB_Bioinformatics/Kevin/BCAC/code/run_zoomlocus.R",snpfile="../result/QCfreq01allcandidatenovelsnps.csv",
               i1=1:nrow(tmp1),outprefix=rep("QCfreq01",nrow(tmp1)))
tmp1=read.csv("../result/QCallcandidatenovelsnps.csv")
tmp2=data.frame(code="/data/BB_Bioinformatics/Kevin/BCAC/code/run_zoomlocus.R",snpfile="../result/QCallcandidatenovelsnps.csv",
                i1=1:nrow(tmp1),outprefix=rep("QC",nrow(tmp1)))
tmp=rbind(tmp,tmp2)
#write.table(tmp,file="run_zoomlocus.swarm",row.names=F,col.names=F,sep="\t",quote=F)
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/run_zoomlocus.swarm -g 72 --module R --time 1-00:00:00 --gres=lscratch:64
#
# for (i in 1:nrow(tmp))
# {
#   print(i)
#   res = tryCatch(
#     expr = {
#       get_zoomlocus(snpfile=tmp$snpfile[i],
#                     i1=tmp$i1[i],outprefix=tmp$outprefix[i])
#     },
#     error = function(e){ 
#       return(NULL)
#     })
#   
#   print(Sys.time())
# }

# snp=allnovelcandidatesnps$ID[1] #"chr12:81015799:G:T"
# snpid=gsub(":","_",snp)
# tmp=unlist(strsplit(snp,":"))
# snpchr=tmp[1] #chr12
# gr_snp=GRanges(seqnames = gsub("chr","",tmp[1]),ranges = IRanges(start=as.integer(tmp[2]),width = 1))
# # generage metafile with columns: MarkerName	Pvalue	ref	alt
# pop="asian"
# myallpvalues=get(paste0("collect_pvalues_",pop))
# myallpvalues=myallpvalues$acatpvalues
# tmp=unlist(strsplit(names(myallpvalues),":"))
# mychr=tmp[seq(1,length(tmp),4)]
# mychr=gsub("chr","",mychr)
# mypos=as.integer(tmp[seq(2,length(tmp),4)])
# gr_myallpvalues=GRanges(seqnames = mychr,ranges = IRanges(start=mypos,width=1))
# mydist=distance(gr_myallpvalues,gr_snp)
# idx=which(mydist<4e6)
# myallpvalues1=myallpvalues[idx]
# tmp=unlist(strsplit(names(myallpvalues1),":"))
# mychr1=tmp[seq(1,length(tmp),4)]
# mypos1=as.integer(tmp[seq(2,length(tmp),4)])
# metadat=data.frame(snpid=names(myallpvalues1),MarkerName=paste0(mychr1,":",mypos1),Pvalue=myallpvalues1)
# gr_myallpvalues1=GRanges(seqnames = gsub("chr","",mychr1),ranges = IRanges(start=mypos1,width=1))
# my_snps <- snpsByOverlaps(all_snps, gr_myallpvalues1, genome = genome)
# my_snps = as.data.frame(my_snps)
# idx = match(paste0(gsub("chr","",mychr1),":",mypos1),paste0(my_snps$seqnames,":",my_snps$pos))
# metadat$MarkerName=my_snps$RefSNP_id[idx]
# idx=which(is.na(metadat$MarkerName))
# tmp=unlist(strsplit(metadat$snpid[idx],":"))
# gr_metadat=GRanges(seqnames = tmp[seq(1,length(tmp),4)],ranges = IRanges(as.integer(tmp[seq(2,length(tmp),4)]),width=1))
# tmp=as.data.frame(liftOver(gr_metadat,chain))
# tmp1=paste0(snpchr,":",mypos1)
# tmp1[tmp$group]=paste0(snpchr,":",tmp$start)
# metadat$MarkerName[idx]=tmp1[idx]
# gr_allmetadat=GRanges(seqnames = mychr1,ranges = IRanges(start=mypos1,width=1))
# tmp=as.data.frame(liftOver(gr_allmetadat,chain))
# tmp1=as.integer(mypos1)
# tmp1[tmp$group]=tmp$start
# metadat$hg19pos=tmp1
# metadat=metadat[grepl("^rs",metadat$MarkerName),]
# metalfile=paste0("../result/zoomlocus/Snp",snpid,"_",pop,"_metalfile.txt")
# write.table(metadat[,c(2,3)],file=metalfile,col.names=T,row.names = F,sep="\t",quote=F)
# snpmapfile=paste0("../result/zoomlocus/Snp",snpid,"_",pop,"_snpmap.txt")
# write.table(metadat[,c(2,1)],file=snpmapfile,row.names = F,col.names = F,sep="\t",quote=F)
# snpfile=paste0("../result/zoomlocus/Snp",snpid,"_",pop,"_snp.txt")
# write.table(metadat$snpid,file=snpfile,row.names = F,col.names = F,quote=F)
# hg19posfile=paste0("../result/zoomlocus/Snp",snpid,"_",pop,"_hg19pos.txt")
# write.table(metadat[,c("MarkerName","hg19pos")],file=hg19posfile,row.names = F,col.names = F,quote=F)
# 
# 
# idx=which(rownames(metascoreinfo4$avail)==snp)
# metascoreinfo4$avail[idx,]
# #work on vcf file
# if (snp %in% asian_onco_pvar$ID)
# {
#   
# }else
# {
#   if (snp %in% asian_icogs_pvar$ID)
#   {
#     cmd=paste0(plink2," --pfile ../result/imp_icogs/asian/asian --extract ",snpfile," --make-bed --out ../result/zoomlocus/Snp",snpid,"_",pop," --memory 64000 --threads 8")
#     system(cmd)
#     cmd=paste0(plink," --bfile ../result/zoomlocus/Snp",snpid,"_",pop," --update-name ",snpmapfile," 1 2 --make-bed --out ../result/zoomlocus/Snp",snpid,"_",pop," --memory 64000 --threads 8")
#     system(cmd)
#     update_bimhg19pos(bimfile=paste0("../result/zoomlocus/Snp",snpid,"_",pop,".bim"),
#                                hg19posfile=hg19posfile)
#     tmp1=read.table(paste0("../result/zoomlocus/Snp",snpid,"_",pop,".bim"))
#     cmd=paste0(plink2," --bfile ../result/zoomlocus/Snp",snpid,"_",pop," --sort-vars "," --make-pgen --out ../result/zoomlocus/Snp",snpid,"_",pop," --memory 64000 --threads 8")
#     system(cmd)
#     tmp2=read.table(paste0("../result/zoomlocus/Snp",snpid,"_",pop,".pvar"))
#     #cmd=paste0(plink2," --bfile ../result/zoomlocus/Snp",snpid,"_",pop," --make-pgen --out ../result/zoomlocus/Snp",snpid,"_",pop," --memory 64000 --threads 8")
#     #system(cmd)
#     cmd=paste0(plink2," --pfile ../result/zoomlocus/Snp",snpid,"_",pop," --recode vcf --out ../result/zoomlocus/Snp",snpid,"_",pop," --memory 64000 --threads 8")
#     system(cmd)
#     tmp=read.table(paste0("../result/zoomlocus/Snp",snpid,"_",pop,".vcf"))
#     tmp$V7="PASS"
#     tmp[tmp=="./."]="0/0"
#     #tmp$V1=paste0("chr",tmp$V1)
#     write.table(tmp,sep="\t",file=paste0("../result/zoomlocus/Snp",snpid,"_",pop,"_pass",".vcf"),row.names = F,col.names = F,quote=F)
#     cmd=paste0("head -7 ","../result/zoomlocus/Snp",snpid,"_",pop,".vcf > ", "../result/zoomlocus/Snp",snpid,"_",pop,"1.vcf")
#     system(cmd)
#     cmd=paste0("cat ","../result/zoomlocus/Snp",snpid,"_",pop,"_pass",".vcf >>","../result/zoomlocus/Snp",snpid,"_",pop,"1.vcf")
#     system(cmd)
#     cmd=paste0("ml samtools;"," bgzip -c ","../result/zoomlocus/Snp",snpid,"_",pop,"1.vcf > ","../result/zoomlocus/Snp",snpid,"_",pop,"1.vcf.gz")
#     system(cmd)
#     
#     cmd=paste0("ml samtools;"," tabix -f -p vcf ","../result/zoomlocus/Snp",snpid,"_",pop,"1.vcf.gz")
#     #cmd=paste0("ml samtools;"," bcftools index ","../result/zoomlocus/Snp",snpid,"_",pop,"1.vcf.gz")
#     
#     system(cmd)
#     vcf_gzfile=paste0("../result/zoomlocus/Snp",snpid,"_",pop,"1.vcf.gz")
#     cmd=paste0("ml locuszoom; locuszoom"," --metal ",metalfile," --flank 1000kb"," --build hg19 --pop ASN --source 1000G_March2012 --refsnp rs138955280 --pvalcol Pvalue --plotonly --prefix zoomlocustest --ld-vcf ",vcf_gzfile," ymax=",20," signifLine=",-log10(5*1e-8)," signifLineColor=red ", " axisTextSize=1.4 legendSize=1.0  axisSize=1.4 xlabPos=-2.9")
#     cmd=paste0("ml locuszoom;ml samtools;ml plink/1.9.0-beta4.4; locuszoom"," --metal ",metalfile," --flank 500kb"," --build hg19 --pop ASN --source 1000G_March2012 --refsnp rs868967 --pvalcol Pvalue --prefix zoomlocustest --ld-vcf ",vcf_gzfile," ymax=",20," signifLine=",-log10(5*1e-8)," signifLineColor=red ", " axisTextSize=1.4 legendSize=1.0  axisSize=1.4 xlabPos=-2.9")
#     system(cmd)
#   }
# }


