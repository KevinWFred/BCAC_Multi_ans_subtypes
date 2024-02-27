#!/usr/bin/env Rscript
#https://github.com/Jingning-Zhang/PROSPER

setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
#where is the package
package='/data/BB_Bioinformatics/Kevin/tools/PROSPER'
#path_example='/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/prosper/'
path_result='/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/prosper/'
path_plink='/usr/local/apps/plink/2.3-alpha/plink2'

# #generate summary files,columns:"rsid"    "chr"     "a1"      "a0"      "beta"    "beta_se" "n_eff"
generate_sum=function(sumstatfile="../result/PRS1/euro_training_sumstats.txt",outprefix="euro_training")
{
  sumstat=as.data.frame(fread(sumstatfile))
  mycolnames=c("rsid","chr","a1","a0","beta","beta_se","n_eff")
  sumstat1=sumstat[,mycolnames]
  sumstat1=sumstat1[sumstat1$beta!=0,]
  outfile=paste0(path_result,outprefix,"_prosper_sumstat.txt")
  fwrite(sumstat1,file=outfile,row.names = F,sep="\t",quote=F)

}
generate_sum(sumstatfile="../result/PRS1/euro_training_sumstats.txt",outprefix="euro_training")
generate_sum(sumstatfile="../result/PRS1/asian_training_sumstats.txt",outprefix="asian_training")

phenoicogs=read.table("../data/concept_750_zhang_icogs_pheno_v15_02_age.txt",header=T,sep="\t")
phenoicogs$ID=phenoicogs$SG_ID
idx=which(colnames(phenoicogs) %in% paste0("pc",1:10))
phenoicogs[,idx]=NULL
tmp=read.table("../result/imp_icogs/merged1.eigenvec")
colnames(tmp)=c("ID",paste0("pc",1:20))
tmp1=unlist(strsplit(tmp$ID,"_"))
tmp$ID=tmp1[seq(1,length(tmp1),2)]
idx=match(phenoicogs$ID,tmp$ID)
phenoicogs=cbind(phenoicogs,tmp[idx,2:11])
phenoicogs$y=NA
phenoicogs$y[which(phenoicogs$Behaviour1==1)]=1
phenoicogs$y[which(is.na(phenoicogs$Behaviour1))]=0
phenoonco=read.table("../data/concept_750_zhang_onco_pheno_v15_02_corrected_age.txt",header=T,sep="\t")
phenoonco$ID=phenoonco$Onc_ID
colnames(phenoonco)=gsub("PC_","pc",colnames(phenoonco))
phenoonco$y=NA
phenoonco$y[which(phenoonco$Behaviour1==1)]=1
phenoonco$y[which(is.na(phenoonco$Behaviour1))]=0

call_subtype=function(pheno=phenoicogs)
{
  pheno$LumA=pheno$LumB=pheno$LumB_HN=pheno$Her2E=pheno$TripN=NA
  idx=which(is.na(pheno$Behaviour1))
  pheno$LumA[idx]=pheno$LumB[idx]=pheno$LumB_HN[idx]=pheno$Her2E[idx]=pheno$TripN[idx]=0
  y.pheno.mis1 <- pheno[,c("ER_status1","PR_status1","HER2_status1","Grade1")]
  
  idx.1 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==0
                 &(y.pheno.mis1[,4]==1|y.pheno.mis1[,4]==2))
  pheno$LumA[idx.1]=1
  #define Luminal-B like
  idx.2 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==1)
  pheno$LumB[idx.2]=1
  #for Luminal B HER2 negative-like
  idx.3 <- which((y.pheno.mis1[,1]==1|y.pheno.mis1[,2]==1)
                 &y.pheno.mis1[,3]==0
                 &y.pheno.mis1[,4]==3)
  pheno$LumB_HN[idx.3]=1
  #for HER2 enriched-like
  idx.4 <- which(y.pheno.mis1[,1]==0&y.pheno.mis1[,2]==0
                 &y.pheno.mis1[,3]==1)
  pheno$Her2E[idx.4]=1
  #for Triple negative
  idx.5 <- which(y.pheno.mis1[,1]==0&y.pheno.mis1[,2]==0
                 &y.pheno.mis1[,3]==0)
  pheno$TripN[idx.5]=1
  return(pheno)
  
}
phenoicogs=call_subtype(pheno=phenoicogs)
phenoonco=call_subtype(pheno=phenoonco)

generate_pheno=function(famprefix="../result/PRS1/euro_onco_tuning",subtype=NULL,outprefix="euro_onco_tuning")
{
  fam=read.table(paste0(famprefix,".fam"))
  idx=match(fam$V1,phenoonco$ID)
  print(sum(is.na(idx)))
  if (is.null(subtype))
  {
    fam$V3=phenoonco$y[idx]
  }
  write.table(fam[,1:3],file=paste0(path_result,outprefix,".pheno"),row.names=F,col.names=F,sep="\t",quote=F)
}
generate_pheno(famprefix="../result/PRS1/euro_onco_tuning",subtype=NULL,outprefix="euro_onco_tuning")
generate_pheno(famprefix="../result/PRS1/asian_onco_tuning",subtype=NULL,outprefix="asian_onco_tuning")
generate_pheno(famprefix="../result/PRS1/euro_onco_validation",subtype=NULL,outprefix="euro_onco_validation")
generate_pheno(famprefix="../result/PRS1/asian_onco_validation",subtype=NULL,outprefix="asian_onco_validation")

#lassosum2
# ${package}/scripts/lassosum2.R \
# --PATH_package ${package} \
# --PATH_out ${path_result}/lassosum2 \
# --PATH_plink ${path_plink} \
# --FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
# --pop EUR,AFR \
# --chrom 1-22 \
# --bfile_tuning ${path_example}/sample_data/EUR/tuning_geno,${path_example}/sample_data/AFR/tuning_geno \
# --pheno_tuning ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
# --bfile_testing ${path_example}/sample_data/EUR/testing_geno,${path_example}/sample_data/AFR/testing_geno \
# --pheno_testing ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
# --testing TRUE \
# --NCORES 8
eursumfile=paste0(path_result,"euro_training_prosper_sumstat.txt")
eassumfile=paste0(path_result,"asian_training_prosper_sumstat.txt")
eurtuning="../result/PRS1/euro_onco_tuning"
eastuning="../result/PRS1/asian_onco_tuning"
eurtuningpheno=paste0(path_result,"euro_onco_tuning.pheno")
eastuningpheno=paste0(path_result,"asian_onco_tuning.pheno")

eurvalidation="../result/PRS1/euro_onco_validation"
easvalidation="../result/PRS1/asian_onco_validation"
eurvalidationpheno=paste0(path_result,"euro_onco_validation.pheno")
easvalidationpheno=paste0(path_result,"asian_onco_validation.pheno")

cmd=paste0("Rscript ",package,"/scripts/lassosum2.R --PATH_package ",package," --PATH_out ",path_result,"lassosum2 --PATH_plink ",path_plink,
           " --FILE_sst ",eursumfile,",",eassumfile," --pop EUR,EAS --chrom 1-22 ",
           " --bfile_tuning ",eurtuning,",",eastuning," --pheno_tuning ",eurtuningpheno,",",eastuningpheno,
           " --bfile_testing ",eurvalidation,",",easvalidation," --pheno_testing ",eurvalidationpheno,",",easvalidationpheno,
           " --testing TRUE  --NCORES 8 --verbose 2")
system(cmd)
print(Sys.time())

#For debug
# opt=list()
# opt$PATH_package=package
# opt$PATH_out=paste0(path_result,"PROSPER")
# opt$FILE_sst=paste0(eursumfile,",",eassumfile)
# opt$pop="EUR,EAS"
# opt$lassosum_param=paste0(path_result,"lassosum2/EUR/optimal_param.txt,",path_result,"lassosum2/EAS/optimal_param.txt")
# opt$chrom="1-22"
# opt$Ll=5
# opt$Lc=5
# opt$NCORES=8
# opt$verbose=2
#PROSPER.R was modified (fixed a bug to compute PRS)
cmd=paste0("Rscript ",package,"/scripts/PROSPER.R --PATH_package ",package," --PATH_out ",path_result,"PROSPER ",
           " --FILE_sst ",eursumfile,",",eassumfile," --pop EUR,EAS --chrom 1-22 ",
           " --lassosum_param ",path_result,"lassosum2/EUR/optimal_param.txt,",path_result,"lassosum2/EAS/optimal_param.txt ",
           " --NCORES 8 ")
system(cmd)
print(Sys.time())


# #used for debug
# opt=list()
# opt$PATH_plink=path_plink
# opt$PATH_package=package
# opt$PATH_out=paste0(path_result,"PROSPER")
# opt$prefix="EAS"
# opt$bfile_tuning=eastuning
# opt$pheno_tuning=eastuningpheno
# opt$bfile_testing=easvalidation
# opt$pheno_testing=easvalidationpheno
# opt$cleanup=F
# opt$testing=T
# opt$SL_library="SL.glmnet,SL.ridge,SL.lm"
# opt$linear_score=T
# opt$covar_tuning=NA
# opt$covar_testing=NA
# opt$NCORES=8
# opt$verbose=2

target_pop='EAS'
#tuning_testing.R was modified (the old code doesn't work if fam file has pheno available)
cmd=paste0(" Rscript ",package,"/scripts/tuning_testing.R --PATH_plink ",path_plink," --PATH_out ",path_result,"/PROSPER",
           " --prefix ",target_pop," --testing TRUE ",
           " --bfile_tuning ",eastuning," --pheno_tuning ",eastuningpheno,
           " --bfile_testing ",easvalidation," --pheno_testing ",easvalidationpheno,
           " --cleanup F --NCORES 8 --verbose 2")
system(cmd)
target_pop='EUR'
cmd=paste0(" Rscript ",package,"/scripts/tuning_testing.R --PATH_plink ",path_plink," --PATH_out ",path_result,"/PROSPER",
           " --prefix ",target_pop," --testing TRUE ",
           " --bfile_tuning ",eurtuning," --pheno_tuning ",eurtuningpheno,
           " --bfile_testing ",eurvalidation," --pheno_testing ",eurvalidationpheno,
           " --cleanup F --NCORES 8 --verbose 2")
system(cmd)
