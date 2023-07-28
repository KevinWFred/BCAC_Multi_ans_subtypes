#!/usr/bin/env Rscript
library(data.table)
#use .o files to check if a job is running
check_swarm=function(prefix="../result/imp_onco/euro/euro",swarm=c(4476647,4476653,4761917,4764222,4864266,4998830,4998831))
{
  nvar=2000
  pvar=as.data.frame(fread(paste0(prefix,".pvar")))
  infolder=paste0(dirname(prefix),"/geno/")
  allfiles=list.files(infolder,"\\w+.traw.gz")
  allposfiles=list.files(infolder,"pos_\\d+.txt")
  allidxs=gsub("pos_","",allposfiles)
  allidxs=as.numeric(gsub(".txt","",allidxs,fixed = T))
  allidxs=allidxs[order(allidxs)]
  n=as.integer(nrow(pvar)/nvar)+1
  tmp=1:n
  missingblocks=tmp[!tmp %in% allidxs]
  if (length(missingblocks)>0 ) warning(paste0(paste0(missingblocks,collapse = ",")," are missing"))
  swarmfolder="./swarm/"
  #get a swarm job table
  #get runnig jobs
  for (i in 1:length(swarm))
  {
    if (i ==1 )
    {
      cmd=paste0("squeue -u wangx53 | grep ",swarm[i],"_"," > ../result/swarmjobs.txt")
      system(cmd)
    }else
    {
      cmd=paste0("squeue -u wangx53 | grep ",swarm[i],"_"," >> ../result/swarmjobs.txt")
      system(cmd)
    }
  }

  swarmtable=data.frame(jobid=1:n,file=NA,swarmid=NA,running=0,host=NA)

  if (file.size("../result/swarmjobs.txt"))
  {
    runningjobs=read.table("../result/swarmjobs.txt")
    for (i in 1:length(swarm))
    {
      if (i==1)
      {
        allswarmfiles=list.files(swarmfolder,paste0("swarm_",swarm[i],"\\w+.o"))
      }else
      {
        allswarmfiles=c(allswarmfiles,list.files(swarmfolder,paste0("swarm_",swarm[i],"\\w+.o")))
      }
    }

    for (i in 1:length(allswarmfiles))
    {
      tmp=read.table(paste0(swarmfolder,allswarmfiles[i]),fill=T)
      idx=which(grepl("---------------------------",tmp$V1))
      jobids=tmp$V5[2:(idx-1)]
      idx1=match(jobids,swarmtable$jobid)
      swarmtable$file[idx1]=allswarmfiles[i]
      swarmtable$host[idx1]=gsub("host:","",tmp$V2[idx+1])
      tmp=unlist(strsplit(allswarmfiles[i],"_"))
      swarmtable$swarmid[idx1]=paste0(tmp[2:3],collapse = "_")
    }
    swarmtable$running=0
    idx=match(runningjobs$V1,swarmtable$swarmid)
    idx=which(swarmtable$swarmid %in% swarmtable$swarmid[idx])
    swarmtable$running[idx]=1
  }
  currenttime=Sys.time()
  for (i in 1:nrow(swarmtable))
  {
    thetime=file.info(paste0("swarm/",swarmtable$file[i]))$mtime
    timediff=as.numeric(difftime(currenttime,thetime,units = "hours"))
    if (timediff>5)
    {
      swarmtable$running[i]=0
    }
  }


  resfolder=paste0(dirname(prefix),"/res/")
  res=data.frame(jobid=1:n,processed=0,running=0,finished=0,cantload=0)
  for (i in 1:nrow(res))
  {
    if (i %% 500==0) cat(i,'..')
    resfile=paste0(resfolder,"res_",i,".RData")
    if (file.exists(resfile))
    {
      loaddata = tryCatch(
        expr = {
          load(resfile)
        },
        error = function(e){
          return(NULL)
        })
      if (!is.null(loaddata))
      {
        res$processed[i]=nrow(all.log.odds)
      }else
      {
        res$cantload[i]=1
      }

    }
  }


  idx=which(res$processed==nvar)
  if (length(idx)>0)
  {
    res$finished[idx]=1
    res$running[idx]=0
  }
  if (any(!is.na(swarmtable$swarmid)))
  {
    res$running=swarmtable$running
    res=cbind(res,swarmtable[,c(2,3,5)])
  }

  #failed jobs (needs to resubmit?)
  #idx=which(res$finished==0 & res$running==0 & res$processed>0)
  #res[idx,]
  return(res)
}

#use res.RData to check if a job is running 
check_swarm=function(prefix="../result/imp_onco/euro/euro")
{
  nvar=2000
  pvar=as.data.frame(fread(paste0(prefix,".pvar")))
  infolder=paste0(dirname(prefix),"/geno/")
  allfiles=list.files(infolder,"\\w+.traw.gz")
  allposfiles=list.files(infolder,"pos_\\d+.txt")
  allidxs=gsub("pos_","",allposfiles)
  allidxs=as.numeric(gsub(".txt","",allidxs,fixed = T))
  allidxs=allidxs[order(allidxs)]
  n=as.integer(nrow(pvar)/nvar)+1
  tmp=1:n
  missingblocks=tmp[!tmp %in% allidxs]
  if (length(missingblocks)>0 ) warning(paste0(paste0(missingblocks,collapse = ",")," are missing"))
  
  currenttime=Sys.time()
  resfolder=paste0(dirname(prefix),"/res/")
  res=data.frame(jobid=1:n,processed=0,running=0,finished=0,cantload=0)
  currenttime=Sys.time()
  for (i in 1:nrow(res))
  {
    if (i %% 500==0) cat(i,'..')
    resfile=paste0(resfolder,"res_",i,".RData")
    if (file.exists(resfile))
    {
      loaddata = tryCatch(
        expr = {
          load(resfile)
        },
        error = function(e){ 
          return(NULL)
        })
      if (!is.null(loaddata))
      {
        res$processed[i]=nrow(all.log.odds)
        thetime=file.info(resfile)$mtime
        timediff=as.numeric(difftime(currenttime,thetime,units = "hours"))
        if (timediff<3)
        {
          res$running[i]=1
        }
        if (exists("ilast"))
        {
          if (ilast>=nvar)
          {
            res$finished[i]=1
            ilast=1
          }
        }
      }else
      {
        res$cantload[i]=1
      }
    }
  }
  
  
  idx=which(res$processed==nvar)
  if (length(idx)>0) 
  {
    res$finished[idx]=1
    #res$running[idx]=0
  }
  
  #failed jobs (needs to resubmit?)
  #idx=which(res$finished==0 & res$running==0 & res$processed>0)
  #res[idx,]
  return(res)  
}
euro_onco=check_swarm()
quantile(euro_onco$processed/2000,c(0,0.01,0.05,0.1,0.5,1))

#euro_icogs=check_swarm(prefix="../result/imp_icogs/euro/euro",swarm=c(4476416,4476544,4761926,4764226,4864280,4864853,4929402,4930181,4947570,4968668,4969145,4970366,4978833,4979109,4996292,4998714,4999414,5059837))
euro_icogs=check_swarm(prefix="../result/imp_icogs/euro/euro")

quantile(euro_icogs$processed/2000,c(0,0.01,0.05,0.1,0.5,1))

resubmitjobs=function(dataopt="onco",pop="euro",swarmres=euro_onco)
{
  #jobs not finished
  idx=which(swarmres$finished==0 & swarmres$running==0 & swarmres$processed>0 & swarmres$processed %%20 ==0)
  idx=which(swarmres$finished==0 & swarmres$running==0 & swarmres$processed>0 & swarmres$processed<1975 & swarmres$jobid<nrow(swarmres))
  tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/intrinsic_subtypes_genome.R",length(idx)),
                 dataopt=dataopt,pop=pop,i1=swarmres$jobid[idx])
  write.table(tmp,file=paste0(pop,"_",dataopt,"3.swarm"),row.names = F,col.names = F,sep="\t",quote=F)
  #jobs never running (it is the second job, and the first job quit)
  idx1=which(swarmres$finished==0 & swarmres$running==0 & swarmres$processed==0)
  if (length(idx1)>0)
  {
    tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/intrinsic_subtypes_genome.R",length(idx1)),
                   dataopt=dataopt,pop=pop,i1=swarmres$jobid[idx1])
    write.table(tmp,file=paste0(pop,"_",dataopt,"4.swarm"),row.names = F,col.names = F,sep="\t",quote=F)
    
  }
}
resubmitjobs()
resubmitjobs(dataopt="icogs",pop="euro",swarmres=euro_icogs)

# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_onco3.swarm -g 10 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2
# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_onco4.swarm -g 10 --module R/4.3 --time=7-00:00:00 --gres=lscratch:8 -p 2

# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_icogs3.swarm -g 15 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2
# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_icogs4.swarm -g 8 --module R/4.3 --time=5-00:00:00 --gres=lscratch:8 -p 2
