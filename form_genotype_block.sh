#!/usr/bin/env bash

#form genotype value matrices

cd /data/BB_Bioinformatics/Kevin/BCAC/code
plink=/usr/local/apps/plink/1.9.0-beta4.4/plink
plink2=/usr/local/apps/plink/2.3-alpha/plink2

#number of var in each block
nvar=2000

i1="$1" #starts with 1
echo $i1
prefix="$2" #"/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro"
outfolder=$(dirname $prefix)/geno
mkdir -p ${outfolder}

allsnp=($(awk '{if (NR>1) print $3}' $prefix.pvar))

s1=$((($i1-1)*($nvar)))
s2=$(($s1+$nvar-1))
if [ "$s1" -ge ${#allsnp[@]} ];then exit 1;fi
if [ "$s2" -gt "${#allsnp[@]}" ];then echo hi;fi
if [ "$s2" -gt "${#allsnp[@]}" ];then s2=${#allsnp[@]};fi
echo ${s1}-${s2}
snp=(${allsnp[@]:$s1:$(($s2-$s1+1))})
printf "%s\n" "${snp[@]}" > $outfolder/pos_$i1.txt

$plink2 --pfile $prefix --extract ${outfolder}/pos_${i1}.txt --recode A-transpose --out ${outfolder}/geno_${i1} --memory 32000 --threads 1
gzip ${outfolder}/geno_${i1}.traw
if [ -f ${outfolder}/geno_${i1}.traw ]; then rm ${outfolder}/geno_${i1}.traw;fi
echo `date +%c`
echo "done"
# #creat swarm file
# pop="euro"
# pop="asian"
pop="african"
# if [ -f ${pop}_onco_block.swarm ]; then rm ${pop}_onco_block.swarm; fi
# totaljob=$((${#allsnp[@]}/$nvar+1))
# echo $totaljob
# prefix=/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/${pop}/${pop}
# for (( i =1 ; i<= $totaljob ; i++ ))
# do
#   echo "/data/BB_Bioinformatics/Kevin/BCAC/code/form_genotype_block.sh $i /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/${pop}/${pop}">>${pop}_onco_block.swarm
# done
# if [ -f ${pop}_icogs_block.swarm ]; then rm ${pop}_icogs_block.swarm; fi
# prefix=/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/${pop}/${pop}
# for (( i =1 ; i<= $totaljob ; i++ ))
# do
#   echo "/data/BB_Bioinformatics/Kevin/BCAC/code/form_genotype_block.sh $i /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/${pop}/${pop}">>${pop}_icogs_block.swarm
# done
# cd swarm
# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/${pop}_onco_block.swarm -g 32   --time 0:20:00 --gres=lscratch:16 -b 100 -p 2
# swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/${pop}_icogs_block.swarm -g 32   --time 0:20:00 --gres=lscratch:16 -b 100 -p 2


  
