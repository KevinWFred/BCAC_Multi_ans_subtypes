library(data.table)
eursum=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/tools/PROSPER/example//summdata/EUR.txt"))
afrsum=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/tools/PROSPER/example//summdata/AFR.txt"))
tmp=intersect(eursum$rsid,afrsum$rsid)
eursum=eursum[match(tmp,eursum$rsid),]
afrsum=afrsum[match(tmp,afrsum$rsid),]
table(eursum$a1==afrsum$a1) #T

eurbim=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/tools/PROSPER/example//sample_data/EUR/tuning_geno.bim"))
idx=grepl(":",eurbim$V2)
tmp=unlist(strsplit(eurbim$V2[idx],":"))
eurbim$rsid=eurbim$V2
eurbim$rsid[idx]=tmp[seq(1,length(tmp),4)]
table(eursum$rsid %in% eurbim$rsid)
idx=match(eursum$rsid,eurbim$rsid)
eurbim=eurbim[idx,]
table(eursum$a1==eurbim$V5)
# FALSE   TRUE 
# 330962 751194
eurpheno=read.table("/data/BB_Bioinformatics/Kevin/tools/PROSPER/example//sample_data/EUR/pheno.fam")
head(eurpheno,2)
# V1 V2       V3
# 1  1  1  0.10247
# 2  2  2 -1.13322

prsres=as.data.frame(fread("/data/BB_Bioinformatics/Kevin/tools/PROSPER/PROSPER_example_results//PROSPER/after_ensemble_EUR/PROSPER_prs_file.txt"))