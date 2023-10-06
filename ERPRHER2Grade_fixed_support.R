#!/usr/bin/env Rscript
#Fit the null model for fixed effect two-stage model

.libPaths(c("/data/wangx53",.libPaths()))
library(readr)
library(devtools)
library(data.table)
library(bc2)
library(bcutility)

i1=1

for (pop in c("asian","african"))
{
  print(pop)
  {
    for (dataopt in c("onco","icogs")) 
    {
      print(dataopt)
      #iCOGS first
      if (dataopt=="icogs")
      {
        genofile=paste0("../result/imp_icogs/",pop,"/geno/geno_",i1,".traw.gz")
        pheno=read.table("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T,sep="\t")
        colnames(pheno)[which(colnames(pheno)=="SG_ID")]="ID"
        samplefile="../result/icogs_samples_750.txt"
      }else
      {
        genofile=paste0("../result/imp_onco/",pop,"/geno/geno_",i1,".traw.gz")
        if (!file.exists(genofile)) stop("no variants in this block!")
        pheno=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
        colnames(pheno)[which(colnames(pheno)=="Onc_ID")]="ID"
        samplefile="../result/onco_samples_750.txt"
      }
      
      #samples included in the analysis
      allsamples=read.table(samplefile)$V1
      #read genotype data
      genodat=as.data.frame(fread(genofile))
      rownames(genodat)=genodat$SNP
      genodat=genodat[,7:ncol(genodat)]
      #process sample names
      tmp=unlist(strsplit(colnames(genodat),"_"))
      tmp1=unlist(strsplit(colnames(genodat)[1],"_"))
      tmp2=tmp[seq(1,length(tmp),length(tmp1))]
      if (length(tmp2) !=ncol(genodat)) stop("Sampel names are not correct")
      colnames(genodat)=tmp2
      ##for now, only use EUR samples
      #allsamples=intersect(allsamples,pheno$ID[which(pheno$EthnicityGeno=="European")])
      #table(colnames(genodat) %in% allsamples)
      allsamples=intersect(allsamples,colnames(genodat)) #pick samples from an ancestry
      #make pheno and genodat consistent in sample orders
      # idx=match(allsamples,colnames(genodat))
      # genodat=genodat[,idx]
      pheno=pheno[match(allsamples,pheno$ID),]
      
      data1=pheno
      #behaviour1 to (0,1)
      data1$Behaviour1[which(is.na(data1$Behaviour1))]=0
      #prepare the phenotypes data for iCOGs 
      y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
      #ER, PR, HER2 is binary with negative as 0 and positive as 1
      #Grade is ordinal with 1, 2, 3
      #controls don't have tumor characteristics data (all NA)
      #cases with missing tumor characteristics marked as 888
      colnames(y.pheno.mis1) <- c("Behaviour","ER","PR","HER2","Grade")
      
      if ("pc1" %in% colnames(pheno)) #iCOGS
      {
        allcov=c(paste0("plinkPC",1:10),"age")
        tmp=read.table("../result/imp_icogs/merged1.eigenvec")
        colnames(tmp)=c("ID",paste0("plinkPC",1:20))
        tmp1=unlist(strsplit(tmp$ID,"_"))
        tmp$ID=tmp1[seq(1,length(tmp1),2)]
        idx=match(pheno$ID,tmp$ID)
        pheno=cbind(pheno,tmp[idx,2:11])
        
      }else
      {
        allcov=c(paste0("PC_",1:10),"age")
      }
      
      x.covar.mis1 =pheno[,allcov]
      
      score.test.support.ERPRHER2Grade <- ScoreTestSupport(
        y.pheno.mis1,
        baselineonly = NULL,
        additive = x.covar.mis1,
        pairwise.interaction = NULL,
        saturated = NULL,
        missingTumorIndicator = 888
      )
      
      outfolder=paste0("../result/imp_",dataopt,"/",pop,"/scoretest")
      if (!dir.exists(outfolder)) dir.create(outfolder)
      save(score.test.support.ERPRHER2Grade,file=paste0(outfolder,"/score.test.support.ERPRHER2Grade.Rdata"))
    }
  }
  
}


