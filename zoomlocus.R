#!/usr/bin/env Rscript
#generate the data used for zoomlocus
#https://my.locuszoom.org/gwas/upload/#
# Your summary statistics should be sorted via chromosome and position before uploading. In general, if your data can be tabixed, there is a good chance that our server can accept and parse the contents. For best results, make sure that your file contains the following information for each variant: chromosome, position, reference and alt alleles, and p-value (or -log10 p).

.libPaths(c("/data/wangx53",.libPaths()))
#load global test p-values allpvalues
load("../result/compute_metapvalues_new.RData")
#read novel SNP data
allnovelsnps0=read.table("../result/allnovelsnps0.txt",header=T)
conditionalres=read.csv("../result/conditional_analysis_res.csv")
allnovelsnps1=conditionalres[conditionalres$acatP<1e-6,]
library(GenomicRanges)
tmp=unlist(strsplit(names(allpvalues),":"))
chr=tmp[seq(1,length(tmp),4)]
chr=as.integer(gsub("chr","",chr))
pos=as.integer(tmp[seq(2,length(tmp),4)])
alt=tmp[seq(3,length(tmp),4)]
ref=tmp[seq(4,length(tmp),4)]
gr_allpvalues=GRanges(seqnames = chr,ranges = IRanges(start=pos,width = 1))

outfolder=paste0("../result/zoomlocus/")
if (!dir.exists(outfolder)) dir.create(outfolder)

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
seqlevelsStyle(genome) <- "NCBI"

generate_zoomlocusfile=function(novelsnp=allnovelsnps0$ID[1],cutoff=5e5)
{
  tmp=unlist(strsplit(novelsnp,":"))
  gr_novelsnp=GRanges(seqnames=gsub("chr","",tmp[1]),ranges=IRanges(start=as.integer(tmp[2],width=1)))
  mydist=distance(gr_allpvalues,gr_novelsnp)
  idx=which(mydist<cutoff)
  dat=data.frame(Chrom=chr[idx],Pos=pos[idx],Ref=ref[idx],Alt=alt[idx],rsID=NA,P=allpvalues[idx])
  if (novelsnp=="chr9:119662503:T:C")
  {
    dat$P[which(rownames(dat)=="chr9:119662503:T:C")]=1e-20
  }
  dat=dat[order(dat$Pos),]
  dat=dat[!duplicated(dat$Pos),]
  positions <- GPos(seqnames = dat$Chrom, pos = dat$Pos)
  
  ## query the genome with positions
  my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)
  
  ## this gives us a GPos object
  my_snps = as.data.frame(my_snps)
  
  idx = match(paste0(dat$Chrom,":",dat$Pos),paste0(my_snps$seqnames,":",my_snps$pos))
  dat$rsID=my_snps$RefSNP_id[idx]
  dat$Ref=my_snps$ref_allele[idx]
  dat$Alt=my_snps$alt_alleles[idx]
  dat$Alt=as.character(my_snps$alt_alleles[idx])
  dat$Alt=gsub("c(","",dat$Alt,fixed = T)
  dat$Alt=gsub(")","",dat$Alt,fixed = T)
  dat$Alt=gsub("\"","",dat$Alt,fixed = T)
  dat$Alt=gsub(" ","",dat$Alt,fixed = T)
  dat=dat[!is.na(dat$rsID),]
  write.table(dat,file=paste0(outfolder,"zoomlocus_",gsub(":","_",novelsnp),".txt"),row.names = F,sep="\t",quote=F)
}

for (i in 1:nrow(allnovelsnps0))
{
  generate_zoomlocusfile(novelsnp=allnovelsnps0$ID[i])
}

for (i in 1:nrow(allnovelsnps1))
{
  generate_zoomlocusfile(novelsnp=allnovelsnps1$ID[i])
}