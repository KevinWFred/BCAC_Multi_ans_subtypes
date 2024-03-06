#!/usr/bin/env Rscript

#create config files for CTSLEB
template=read.table("CT_SLEB_EUR.config")
target_test_plink=paste0("/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/tuning_validation")
pops=c("euro","asian")
subtypes=c("LumA","LumB","LumB_HN","Her2E","TripN")
basedir="/data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/"
for (i in 1:length(pops))
{
  pop=pops[i]
  for (j in 1:length(subtypes))
  {
    subtype=subtypes[j]
    configfile=paste0(basedir,"config/CT_SLEB_",pop,"_",subtype,".config")
    pop1=pops[!pops %in% pop] #reference
    reference_sumstat=paste0("/data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/",pop1,"/",pop1,"_",subtype,"_CTSLEB_sumstats.txt")
    target_sumstat=paste0("/data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/",pop,"/",pop,"_",subtype,"_CTSLEB_sumstats.txt")
    reference_clumping_plink=paste0("/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/asian_onco_tuning")
    target_clumping_plink=paste0("/data/BB_Bioinformatics/Kevin/BCAC/result/PRS1/euro_onco_tuning")
    output_dir=paste0(basedir,pop,"/",subtype,"/")
    if (!dir.exists(output_dir)) dir.create(output_dir)
    outprefix=paste0(pop,"_",subtype)
    template$V2[1]=reference_sumstat
    template$V2[2]=target_sumstat
    template$V2[3]=reference_clumping_plink
    template$V2[4]=target_clumping_plink
    template$V2[5]=target_test_plink
    template$V2[6]=output_dir
    template$V2[7]=outprefix
    write.table(template,file=configfile,row.names = F,col.names = F,sep="\t",quote=F)
  }
}

dir.create("/data/BB_Bioinformatics/Kevin/BCAC/code/swarm/ctsleb")
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code/swarm/ctsleb")
#generate swarmfiles
clumptemplate=read.table("/data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EAS_runclump.swarm")
for (i in 1:length(pops))
{
  pop=pops[i]
  for (j in 1:length(subtypes))
  {
    subtype=subtypes[j]
    swarmfile=paste0(basedir,"swarmfile/CT_SLEB_clump_",pop,"_",subtype,".swarm")
    configfile=paste0(basedir,"config/CT_SLEB_",pop,"_",subtype,".config")
    clumptemplate$V3=configfile
    write.table(clumptemplate,file=swarmfile,row.names = F,col.names = F,sep="\t",quote=F)
    cmd=paste0("swarm  -f  ",swarmfile," --module R -g 72 --time=5-00:00:00 --gres=lscratch:72")
    print(cmd)
    #system(cmd)
    system("sleep 2")
  }
}
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_euro_LumA.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096702
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_euro_LumB.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096724
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_euro_LumB_HN.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096731
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_euro_Her2E.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096739
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_euro_TripN.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096744
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_asian_LumA.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096746
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_asian_LumB.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096748
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_asian_LumB_HN.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096750
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_asian_Her2E.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096751
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_clump_asian_TripN.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20096754

scoretemplate=read.table("/data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EAS_runscore.swarm")
for (i in 1:length(pops))
{
  pop=pops[i]
  for (j in 1:length(subtypes))
  {
    subtype=subtypes[j]
    swarmfile=paste0(basedir,"swarmfile/CT_SLEB_score_",pop,"_",subtype,".swarm")
    configfile=paste0(basedir,"config/CT_SLEB_",pop,"_",subtype,".config")
    scoretemplate$V3=configfile
    write.table(scoretemplate,file=swarmfile,row.names = F,col.names = F,sep="\t",quote=F)
    cmd=paste0("swarm  -f  ",swarmfile," --module R -g 72 --time=5-00:00:00 --gres=lscratch:72")
    print(cmd)
    #system(cmd)
    #system("sleep 2")
  }
}
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_euro_LumA.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265505
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_euro_LumB.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265511
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_euro_LumB_HN.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265517
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_euro_Her2E.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265522
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_euro_TripN.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265528
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_asian_LumA.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265533
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_asian_LumB.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265539
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_asian_LumB_HN.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265544
# [1] "swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_score_asian_Her2E.swarm --module R -g 72 --time=5-00:00:00 --gres=lscratch:72"
# 20265549
# [1] template
# 20265555

ebscoretemplate=read.table("/data/BB_Bioinformatics/Kevin/BCAC/code/CTSLEB_EAS_runEBscore.swarm")
for (i in 1:length(pops))
{
  pop=pops[i]
  for (j in 1:length(subtypes))
  {
    subtype=subtypes[j]
    swarmfile=paste0(basedir,"swarmfile/CT_SLEB_ebscore_",pop,"_",subtype,".swarm")
    configfile=paste0(basedir,"config/CT_SLEB_",pop,"_",subtype,".config")
    ebscoretemplate$V3=configfile
    ebscoretemplate$V4=subtype
    write.table(ebscoretemplate,file=swarmfile,row.names = F,col.names = F,sep="\t",quote=F)
    cmd=paste0("swarm  -f  ",swarmfile," --module R -g 72 --time=5-00:00:00 --gres=lscratch:72")
    print(cmd)
    system(cmd)
    system("sleep 2")
  }
}

#super,target EAS or EUR
tmp=NULL
code="/data/BB_Bioinformatics/Kevin/BCAC/code/run_CT_SLEB_super.R"
for (i in 1:length(pops))
{
  pop=pops[i]
  for (j in 1:length(subtypes))
  {
    subtype=subtypes[j]
    configfile=paste0(basedir,"config/CT_SLEB_",pop,"_",subtype,".config")
    tmp=rbind(tmp,data.frame(code=code,config=configfile,opt_tuning2=0,subtype=subtype))
  }
}
write.table(tmp,file=paste0(basedir,"swarmfile/CT_SLEB_super_subtype.swarm"),row.names = F,col.names = F,sep=" ",quote=F)
cmd="swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_super_subtype.swarm --module R -g 32 --time=3-00:00:00 --gres=lscratch:32"
system(cmd)

#super,target EAS or EUR
tmp=NULL
pop="euro"
code="/data/BB_Bioinformatics/Kevin/BCAC/code/run_CT_SLEB_super.R"
for (j in 1:length(subtypes))
{
  subtype=subtypes[j]
  configfile=paste0(basedir,"config/CT_SLEB_",pop,"_",subtype,".config")
  tmp=rbind(tmp,data.frame(code=code,config=configfile,opt_tuning2=1,subtype=subtype))
}
write.table(tmp,file=paste0(basedir,"swarmfile/CT_SLEB_super_subtype_2tuning.swarm"),row.names = F,col.names = F,sep=" ",quote=F)
cmd="swarm  -f  /data/BB_Bioinformatics/Kevin/BCAC/result/PRS_subtype/ctsleb/swarmfile/CT_SLEB_super_subtype_2tuning.swarm --module R -g 32 --time=3-00:00:00 --gres=lscratch:32"
system(cmd)