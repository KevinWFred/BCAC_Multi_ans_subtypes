#!/usr/bin/env bash

#form genotype value matrices

cd /data/BB_Bioinformatics/Kevin/BCAC/code
plink=/usr/local/apps/plink/1.9.0-beta4.4/plink
plink2=/usr/local/apps/plink/2.3-alpha/plink2

#number of var in each block
nvar=2000

i1="$1" #starts with 0
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
# if [ -f euro_onco_block.swarm ]; then rm euro_onco_block.swarm; fi
# totaljob=$((${#allsnp[@]}/$nvar+1))
# echo $totaljob
# for (( i =1 ; i<= $totaljob ; i++ ))
# do
#   echo "/data/BB_Bioinformatics/Kevin/BCAC/code/form_genotype_block.sh $i /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro">>euro_onco_block.swarm
# done
# if [ -f euro_icogs_block.swarm ]; then rm euro_icogs_block.swarm; fi
# for (( i =1 ; i<= $totaljob ; i++ ))
# do
#   echo "/data/BB_Bioinformatics/Kevin/BCAC/code/form_genotype_block.sh $i /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro">>euro_icogs_block.swarm
# done

#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_onco_block.swarm -g 32   --time 0:20:00 --gres=lscratch:16 -b 100 -p 2
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_icogs_block.swarm -g 32   --time 0:20:00 --gres=lscratch:16 -b 100 -p 2


  
