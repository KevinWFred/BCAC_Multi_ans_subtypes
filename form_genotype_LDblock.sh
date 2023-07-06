#!/usr/bin/env bash

#form genotype value matrices

cd /data/BB_Bioinformatics/Kevin/BCAC/code
plink=/usr/local/apps/plink/1.9.0-beta4.4/plink
plink2=/usr/local/apps/plink/2.3-alpha/plink2

#euro starts with an empty line
n=($(grep -n "chr1" /data/BB_Bioinformatics/ProjectData/ld_block/ld_block_EUR |cut -f1 -d: ))
chr=($(awk -v n=${n[0]} '{if (NR>=n) {gsub("chr","",$1); print $1} }' /data/BB_Bioinformatics/ProjectData/ld_block/ld_block_EUR))
startpos=($(awk -v n=${n[0]} '{if (NR>=n) print $2}' /data/BB_Bioinformatics/ProjectData/ld_block/ld_block_EUR))
endpos=($(awk -v n=${n[0]} '{if (NR>=n) print $3}' /data/BB_Bioinformatics/ProjectData/ld_block/ld_block_EUR))
for (( i = 0 ; i < ${#endpos[@]}-1 ; i++ )) do  (( endpos[$i]=${endpos[$i]} -1 )) ; done

i1="$1" #starts with 0
echo $1
prefix="$2" #"/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro"
outfolder=$(dirname $prefix)/geno
mkdir -p ${outfolder}
echo "${chr[$i1]} ${startpos[$i1]} ${endpos[$i1]}" > ${outfolder}/pos_${i1}.txt 
$plink2 --pfile $prefix --extract range ${outfolder}/pos_${i1}.txt --recode A-transpose --out ${outfolder}/geno_${i1} --memory 32000 --threads 1
gzip ${outfolder}/geno_${i1}.traw
if [ -f ${outfolder}/geno_${i1}.traw ]; then rm ${outfolder}/geno_${i1}.traw;fi
echo `date +%c`
echo "done"
#creat swarm file
# if [ -f euro_onco_geno.swarm ]; then rm euro_onco_geno.swarm; fi
# for (( i =0 ; i< ${#endpos[@]} ; i++ ))
# do
#   echo "/data/BB_Bioinformatics/Kevin/BCAC/code/form_genotype_block.sh $i /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro">>euro_onco_geno.swarm
# done
# if [ -f euro_icogs_geno.swarm ]; then rm euro_icogs_geno.swarm; fi
# for (( i =0 ; i< ${#endpos[@]} ; i++ ))
# do
#   echo "/data/BB_Bioinformatics/Kevin/BCAC/code/form_genotype_block.sh $i /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro">>euro_icogs_geno.swarm
# done

#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_onco_geno.swarm -g 32   --time 0:20:00 --gres=lscratch:16 -b 200
#swarm -f /data/BB_Bioinformatics/Kevin/BCAC/code/euro_icogs_geno.swarm -g 32   --time 0:20:00 --gres=lscratch:16 -b 200


  
