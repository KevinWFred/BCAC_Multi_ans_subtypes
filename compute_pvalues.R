#!/usr/bin/env Rscript

readres=function(prefixonco="../result/imp_onco/euro/euro",prefixicogs="../result/imp_icogs/euro/euro",
                 nvar=2000,outprefix="../result/euro")
{
  print(Sys.time())
  pvar=as.data.frame(fread(paste0(prefixicogs,".pvar")))
  nicogs=as.integer(nrow(pvar)/nvar)+1
  pvar=as.data.frame(fread(paste0(prefixonco,".pvar")))
  nonco=as.integer(nrow(pvar)/nvar)+1
  or_icogs=or_onco=sigma_icogs=sigma_onco=NULL
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
        or_onco=rbind(or_onco,all.log.odds)
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
        or_icogs=rbind(or_icogs,all.log.odds)
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
  save(or_icogs,or_onco,sigma_icogs,sigma_onco,file=paste0(outprefix,"_or_sigma.RData"))
  print(Sys.time())
}