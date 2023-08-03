#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
setwd("/data/BB_Bioinformatics/Kevin/BCAC/code")
library(data.table)

startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
LogoddsMetaAnalysis=function (logodds1, sigma1, logodds2, sigma2) 
{
  sigma1.inv <- solve(sigma1)
  sigma2.inv <- solve(sigma2)
  sigma.meta <- solve(sigma1.inv + sigma2.inv)
  logodds.meta <- sigma.meta %*% (sigma1.inv %*% logodds1 + 
                                    sigma2.inv %*% logodds2)
  return(list(logodds.meta = logodds.meta, sigma.meta = sigma.meta))
}

# readres=function(prefixonco="../result/imp_onco/euro/euro",prefixicogs="../result/imp_icogs/euro/euro",
#                  nvar=2000,outprefix="../result/euro")
# {
#   print(Sys.time())
#   pvar=as.data.frame(fread(paste0(prefixicogs,".pvar")))
#   nicogs=as.integer(nrow(pvar)/nvar)+1
#   pvar=as.data.frame(fread(paste0(prefixonco,".pvar")))
#   nonco=as.integer(nrow(pvar)/nvar)+1
#   beta_icogs=sigma_icogs=data.frame(matrix(NA,nrow=nicogs*nvar,ncol=5))
#   beta_onco=sigma_onco=data.frame(matrix(NA,nrow=nonco*nvar,ncol=5))
#   allsigma_icogs=list()
#   allsigma_onco=list()
#   #beta_icogs=beta_onco=sigma_icogs=sigma_onco=NULL
#   icogsresfolder=paste0(dirname(prefixicogs),"/res/")
#   oncoresfolder=paste0(dirname(prefixonco),"/res/")
#   for (i in 1:nonco)
#   {
#     if (i %% 500 ==0) cat(i,'..')
#     resfile=paste0(oncoresfolder,"res_",i,".RData")
#     if (file.exists(resfile))
#     {
#       loaddata = tryCatch(
#         expr = {
#           load(resfile)
#         },
#         error = function(e){ 
#           return(NULL)
#         })
#       if (!is.null(loaddata))
#       {
#         idx=match(names(all.sigma.log.odds),rownames(all.log.odds))
#         all.log.odds=all.log.odds[idx,]
#         allsigma_onco=c(allsigma_onco,all.sigma.log.odds)
#         startidx=(i-1)*nvar+1
#         endidx=startidx+nrow(all.log.odds)-1
#         beta_onco[startidx:endidx,]=all.log.odds
#         rownames(beta_onco)[startidx:endidx]=rownames(all.log.odds)
#         tmp=do.call("rbind",all.sigma.log.odds)
#         colid=rep(1:5,times=nrow(all.log.odds))
#         rowid=1:nrow(tmp)
#         idx=as.matrix(data.frame(row=rowid,col=colid))
#         tmp1=tmp[idx]
#         tmp2=matrix(tmp1,ncol=5,byrow = T)
#         rownames(tmp2)=rownames(all.log.odds)
#         sigma_onco[startidx:endidx,]=tmp2#rbind(sigma_onco,tmp2)
#         rownames(sigma_onco)[startidx:endidx]=rownames(all.log.odds)
#       }else
#       {
#         warning(paste0(resfile," is missing"))
#       }
#     }
#   }
#   idx=complete.cases(beta_onco)
#   beta_onco=beta_onco[idx,]
#   print(dim(beta_onco))
#   sigma_onco=sigma_onco[idx,]
#   print(Sys.time())
#   for (i in 1:nicogs)
#   {
#     if (i %% 500 ==0) cat(i,'..')
#     resfile=paste0(icogsresfolder,"res_",i,".RData")
#     if (file.exists(resfile))
#     {
#       loaddata = tryCatch(
#         expr = {
#           load(resfile)
#         },
#         error = function(e){ 
#           return(NULL)
#         })
#       if (!is.null(loaddata))
#       {
#         idx=match(names(all.sigma.log.odds),rownames(all.log.odds))
#         all.log.odds=all.log.odds[idx,]
#         allsigma_icogs=c(allsigma_icogs,all.sigma.log.odds)
#         startidx=(i-1)*nvar+1
#         endidx=startidx+nrow(all.log.odds)-1
#         beta_icogs[startidx:endidx,]=all.log.odds
#         rownames(beta_icogs)[startidx:endidx]=rownames(all.log.odds)
#         tmp=do.call("rbind",all.sigma.log.odds)
#         colid=rep(1:5,times=nrow(all.log.odds))
#         rowid=1:nrow(tmp)
#         idx=as.matrix(data.frame(row=rowid,col=colid))
#         tmp1=tmp[idx]
#         tmp2=matrix(tmp1,ncol=5,byrow = T)
#         rownames(tmp2)=rownames(all.log.odds)
#         sigma_icogs[startidx:endidx,]=tmp2
#         rownames(sigma_icogs)[startidx:endidx]=rownames(all.log.odds)
#       }else
#       {
#         warning(paste0(resfile," is missing"))
#       }
#     }
#   }
#   idx=complete.cases(beta_icogs)
#   beta_icogs=beta_icogs[idx,]
#   sigma_icogs=sigma_icogs[idx,]
#   print(dim(beta_icogs))
#   save(beta_icogs,beta_onco,sigma_icogs,sigma_onco,allsigma_icogs,allsigma_onco,file=paste0(outprefix,"_beta_sigma.RData"))
#   print(Sys.time())
# }

readres=function(pop="euro",nvar=2000)
{
  print(Sys.time())
  #prefixonco="../result/imp_onco/euro/euro"
  prefixonco=paste0("../result/imp_onco/",pop,"/",pop)
  #prefixicogs="../result/imp_icogs/euro/euro",
  prefixicogs=paste0("../result/imp_icogs/",pop,"/",pop)
  #outprefix="../result/euro"
  outprefix=paste0("../result/",pop)
  pvar=as.data.frame(fread(paste0(prefixicogs,".pvar")))
  nicogs=as.integer(nrow(pvar)/nvar)+1
  pvar=as.data.frame(fread(paste0(prefixonco,".pvar")))
  nonco=as.integer(nrow(pvar)/nvar)+1
  #beta_icogs=sigma_icogs=data.frame(matrix(NA,nrow=nicogs*nvar,ncol=5))
  #beta_onco=sigma_onco=data.frame(matrix(NA,nrow=nonco*nvar,ncol=5))
  allsigma_icogs=list()
  beta_icogs=sigma_icogs=NULL
  allsigma_onco=list()
  beta_onco=sigma_onco=NULL
  icogsresfolder=paste0(dirname(prefixicogs),"/res/")
  oncoresfolder=paste0(dirname(prefixonco),"/res/")
  for (i in 1:nonco)
  {
    if (i %% 500 ==0) cat(i,'..')
    resfile=paste0(oncoresfolder,"res_",i,".RData")
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
        idx=match(names(all.sigma.log.odds),rownames(all.log.odds))
        all.log.odds=all.log.odds[idx,]
        allsigma_onco=c(allsigma_onco,all.sigma.log.odds)
        beta_onco=rbind(beta_onco,all.log.odds)
        tmp=do.call("rbind",all.sigma.log.odds)
        colid=rep(1:5,times=nrow(all.log.odds))
        rowid=1:nrow(tmp)
        idx=as.matrix(data.frame(row=rowid,col=colid))
        tmp1=tmp[idx]
        tmp2=matrix(tmp1,ncol=5,byrow = T)
        rownames(tmp2)=rownames(all.log.odds)
        sigma_onco=rbind(sigma_onco,tmp2)
      }else
      {
        warning(paste0(resfile," is missing"))
      }
    }
  }
  print(dim(beta_onco))
  print(Sys.time())
  for (i in 1:nicogs)
  {
    if (i %% 500 ==0) cat(i,'..')
    resfile=paste0(icogsresfolder,"res_",i,".RData")
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
        idx=match(names(all.sigma.log.odds),rownames(all.log.odds))
        all.log.odds=all.log.odds[idx,]
        allsigma_icogs=c(allsigma_icogs,all.sigma.log.odds)
        beta_icogs=rbind(beta_icogs,all.log.odds)
        tmp=do.call("rbind",all.sigma.log.odds)
        colid=rep(1:5,times=nrow(all.log.odds))
        rowid=1:nrow(tmp)
        idx=as.matrix(data.frame(row=rowid,col=colid))
        tmp1=tmp[idx]
        tmp2=matrix(tmp1,ncol=5,byrow = T)
        rownames(tmp2)=rownames(all.log.odds)
        sigma_icogs=rbind(sigma_icogs,tmp2)
      }else
      {
        warning(paste0(resfile," is missing"))
      }
    }
  }
  print(dim(beta_icogs))
  save(beta_icogs,beta_onco,sigma_icogs,sigma_onco,allsigma_icogs,allsigma_onco,file=paste0(outprefix,"_beta_sigma.RData"))
  print(Sys.time())
}

#library(bc2)
#compute p-values for icogs/onco array
compute_p=function(pop="euro")
{
  print(Sys.time())
  #outprefix="../result/euro"
  outprefix=paste0("../result/",pop)
  #inputfile="../result/euro_beta_sigma.RData"
  inputfile=paste0("../result/",pop,"_beta_sigma.RData")
  load(inputfile)
  comvar=intersect(rownames(beta_icogs),rownames(beta_onco))
  print(length(comvar))
  allvar=unique(c(rownames(beta_icogs),rownames(beta_onco)))
  print(length(allvar))
  icogsvar=rownames(beta_icogs)
  oncovar=rownames(beta_onco)
  
  allres=data.frame(matrix(NA,nrow=length(allvar),ncol=45))
  rownames(allres)=allvar
  colnames(allres)=c(paste0("icogs_beta",1:5),paste0("icogs_se",1:5),paste0("icogs_p",1:5),
                     paste0("onco_beta",1:5),paste0("onco_se",1:5),paste0("onco_p",1:5),
                     paste0("meta_beta",1:5),paste0("meta_se",1:5),paste0("meta_p",1:5))
  idx=match(icogsvar,allvar)
  allres[idx,1:5]=beta_icogs
  allres[idx,6:10]=sqrt(sigma_icogs)
  z=beta_icogs/sqrt(sigma_icogs)
  allres[idx,11:15] = 2*pnorm(as.matrix(-abs(z)), lower.tail = T)
 
  idx=match(oncovar,allvar)
  allres[idx,16:20]=beta_onco
  allres[idx,21:25]=sqrt(sigma_onco)
  z=beta_onco/sqrt(sigma_onco)
  allres[idx,26:30] = 2*pnorm(as.matrix(-abs(z)), lower.tail = T) 
  
  # #add meta results
  # idx=match(comvar,allvar)
  # idx1=match(comvar,icogsvar)
  # idx2=match(comvar,oncovar)
  # ncomvar=length(comvar)
 #  metaallsigma=list()
 #  # metares=data.frame(matrix(NA,nrow=idxend-idxstart+1,ncol=15))
 #  # rownames(metares)=comvar[idxstart:idxend]
 #  # colnames(metares)=c(paste0("meta_beta",1:5),paste0("meta_se",1:5),paste0("meta_p",1:5))
 #  for (i in 1:ncomvar)
 #  {
 #    if (i %% 100000 ==0) cat(i/100000,'..')
 #    varname=rownames(beta_icogs)[idx1[i]]
 #    meta.result = tryCatch(
 #      expr = {
 #        LogoddsMetaAnalysis(unlist(beta_icogs[idx1[i],]),
 #                            as.matrix(allsigma_icogs[[varname]]),
 #                            unlist(beta_onco[idx2[i],]),
 #                            as.matrix(allsigma_onco[[varname]]))
 #      },
 #      error = function(e){ 
 #        return(NULL)
 #      })
 #    if (!is.null(meta.result))
 #    {
 #      second.stage.logodds.meta <- meta.result[[1]]
 #      second.stage.sigma.meta <- meta.result[[2]]
 #      # metares[i,1:5]=second.stage.logodds.meta
 #      # metares[i,6:10]=sqrt(diag(second.stage.sigma.meta))
 #      allres[idx[i],31:35]=second.stage.logodds.meta
 #      allres[idx[i],36:40]=sqrt(diag(second.stage.sigma.meta))
 #      metaallsigma[[varname]]=second.stage.sigma.meta
 #      
 #      #function to calculate z-statistics
 #      z = second.stage.logodds.meta/sqrt(diag(second.stage.sigma.meta))
 #      #function to calculate P-value
 #      p = 2*pnorm(-abs(z), lower.tail = T)
 #      #metares[i,11:15]=p
 #      allres[idx[i],41:45]=p
 #    }
 #  }
  
  save(allres,file=paste0(outprefix,"_beta_sigma_pvalue.RData"))
  print(Sys.time())
  #return(allres)
}


print(paste0("host:",as.character(Sys.info()["nodename"])))
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
print(Sys.time())
print(args)
i1 = as.numeric(args[1]) #block ID, starts with 1
pop = as.character(args[2]) #euro
ntask=1000 #total number of jobs
#inputfile="../result/euro_beta_sigma.RData"
inputfile=paste0("../result/",pop,"_beta_sigma.RData")
load(inputfile)
#inputfile="../result/euro_beta_sigma_pvalue.RData"
inputfile=paste0("../result/",pop,"_beta_sigma_pvalue.RData")
load(inputfile)
comvar=intersect(rownames(beta_icogs),rownames(beta_onco))
ncomvar=length(comvar)
print(length(comvar))
allvar=unique(c(rownames(beta_icogs),rownames(beta_onco)))
print(length(allvar))
icogsvar=rownames(beta_icogs)
oncovar=rownames(beta_onco)
startendidx=startend(ncomvar,ntask,i1)
print(paste0("startendidx: ",startendidx))
tmp=sum( sapply(ls(),function(x){object.size(get(x))}))/1E6
print(paste0("memory size: ",tmp))

compute_metap=function(pop="euro",idxstart=1,idxend=10)
{
  # outprefix=paste0("../result/",pop)
  # #inputfile="../result/euro_beta_sigma.RData"
  # inputfile=paste0("../result/",pop,"_beta_sigma.RData")
  # load(inputfile)
  # #inputfile="../result/euro_beta_sigma_pvalue.RData"
  # inputfile=paste0("../result/",pop,"_beta_sigma_pvalue.RData")
  # load(inputfile)
  # 
  # comvar=intersect(rownames(beta_icogs),rownames(beta_onco))
  # print(length(comvar))
  # allvar=unique(c(rownames(beta_icogs),rownames(beta_onco)))
  # print(length(allvar))
  # icogsvar=rownames(beta_icogs)
  # oncovar=rownames(beta_onco)
  
  idx=match(comvar,allvar)
  idx1=match(comvar,icogsvar)
  idx2=match(comvar,oncovar)
  ncomvar=length(comvar)
  metares=data.frame(matrix(NA,nrow=idxend-idxstart+1,ncol=15))
  rownames(metares)=comvar[idxstart:idxend]
  colnames(metares)=c(paste0("meta_beta",1:5),paste0("meta_se",1:5),paste0("meta_p",1:5))
  metaallsigma=list()
  
  j=1 #index of result
  for (i in idxstart:idxend)
  {
    if (j %% 1000 ==0) cat(j,'..')
    varname=rownames(beta_icogs)[idx1[i]]
    meta.result = tryCatch(
      expr = {
        LogoddsMetaAnalysis(unlist(beta_icogs[idx1[i],]),
                            as.matrix(allsigma_icogs[[varname]]),
                            unlist(beta_onco[idx2[i],]),
                            as.matrix(allsigma_onco[[varname]]))
      },
      error = function(e){ 
        return(NULL)
      })
    if (!is.null(meta.result))
    {
      second.stage.logodds.meta <- meta.result[[1]]
      second.stage.sigma.meta <- meta.result[[2]]
      metares[j,1:5]=second.stage.logodds.meta
      metares[j,6:10]=sqrt(diag(second.stage.sigma.meta))
      metaallsigma[[varname]]=second.stage.sigma.meta
      
      #function to calculate z-statistics
      z = second.stage.logodds.meta/sqrt(diag(second.stage.sigma.meta))
      #function to calculate P-value
      p = 2*pnorm(-abs(z), lower.tail = T)
      metares[j,11:15]=p
    }
    j=j+1
  }
  allmetares=list(metares=metares,metaallsigma=metaallsigma)
  outfolder=paste0("../result/metaresult/",pop)
  if(!dir.exists(outfolder)) dir.create(outfolder)
  save(allmetares,file=paste0(outfolder,"/res_",i1,".RData"))
}

res=compute_metap(pop="euro",idxstart=startendidx[1],idxend=startendidx[2])

print(Sys.time())
print("done")

createswarm=function(pop="euro")
{
  swarmjobs=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/BCAC/code/compute_pvalues.R",1000),
                       i1=1:1000,pop=pop)
  write.table(swarmjobs,file=paste0("compute_metap_",pop,".swarm"),row.names=F,col.names=F,sep="\t",quote=F)
}
#tmp=sum( sapply(ls(),function(x){object.size(get(x))}))/1E6
#sort( sapply(ls(),function(x){object.size(get(x))})) /1E6
#createswarm()
#cd swarm
#5504190
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/compute_metap_euro.swarm -g 30 --module R/4.3 --partition=quick --time=4:00:00 --gres=lscratch:30 -p 2


