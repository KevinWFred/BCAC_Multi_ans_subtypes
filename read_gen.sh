#!/usr/bin/env bash

cd /data/BB_Bioinformatics/Kevin/BCAC/code
plink=/usr/local/apps/plink/1.9.0-beta4.4/plink
plink2=/usr/local/apps/plink/2.3-alpha/plink2
ml samtools
#https://community.dnanexus.com/s/question/0D5t000004AeqshCAB/utility-of-hardyweinberg-equilibrium-filtering-in-ukb-genomic-data-p1e15-is-not-a-good-cutoff
#plinksample and pheno files were generated in check_data.R
read_seg(){
  local infolder="$1"
  local prefix="$2"
  local outfolder="$3"
  local outprefix="$4"
  rm ${outfolder}${outprefix}_mergelist.txt
  for chr in {1..22}
  do
    echo $chr
    if [ -f ${infolder}${prefix}${chr}_p1.bgen ] #chr has segments
    then
      echo $chr
      segs=()
      for seg in ${infolder}${prefix}${chr}_p*.bgen
      do
        seg=$(basename ${seg})
        seg=${seg#${prefix}${chr}_p}
        seg=${seg%.bgen}
        echo $seg
        segs+=($seg)
      done
      sorted=($(printf '%s\n' "${segs[@]}"|sort -n))
      echo ${sorted[@]}
      for seg in ${sorted[@]}
      do
        echo $seg
        #$plink2 --bgen ${infolder}${prefix}${chr}_p${seg}.bgen ref-unknown --sample ${infolder}${prefix}1_p1.sample --oxford-single-chr ${chr} --maf 0.005 --geno 0.1 --keep ${infolder}${prefix}plinksamples.txt --pheno ${infolder}${prefix}pheno.txt --hwe 1e-50 midp keep-fewhet --make-pgen --out ${outfolder}${outprefix}_chr.${chr}.${seg} --memory 300000 --threads 45
        #there are a few duplicated SNPIDs
        $plink2 --bgen ${infolder}${prefix}${chr}_p${seg}.bgen ref-unknown --rm-dup force-first --sample ${infolder}${prefix}.sample --oxford-single-chr ${chr} --maf 0.005 --keep ${infolder}${prefix}plinksamples.txt --pheno ${infolder}${prefix}pheno.txt --make-pgen --out ${outfolder}${outprefix}_chr.${chr}.${seg} --memory 300000 --threads 45
        echo ${outfolder}${outprefix}_chr.${chr}.${seg} >> ${outfolder}${outprefix}_mergelist.txt
      done
    else #no segments
      #$plink2 --bgen ${infolder}${prefix}${chr}.bgen ref-unknown --sample ${infolder}${prefix}1.sample --oxford-single-chr ${chr} --maf 0.005 --geno 0.1 --keep ${infolder}${prefix}plinksamples.txt --pheno ${infolder}${prefix}pheno.txt --hwe 1e-50 midp keep-fewhet --make-pgen --out ${outfolder}${outprefix}_chr.${chr} --memory 300000 --threads 45
      $plink2 --bgen ${infolder}${prefix}${chr}.bgen ref-unknown --rm-dup force-first --sample ${infolder}${prefix}.sample --oxford-single-chr ${chr} --maf 0.005 --keep ${infolder}${prefix}plinksamples.txt --pheno ${infolder}${prefix}pheno.txt --make-pgen --out ${outfolder}${outprefix}_chr.${chr} --memory 300000 --threads 45
      
      echo ${outfolder}${outprefix}_chr.${chr} >> ${outfolder}${outprefix}_mergelist.txt
    fi
  done
  echo "merge data"
  $plink2 --pmerge-list ${outfolder}${outprefix}_mergelist.txt --make-pgen --out ${outfolder}${outprefix} --memory 300000 --threads 45
}

# read_seg(){
#   local infolder="$1"
#   local prefix="$2"
#   local outfolder="$3"
#   local outprefix="$4"
#   rm ${outfolder}${outprefix}_mergelist.txt
#   for chr in {1..22}
#   do
#     echo $chr
#     if [ -f ${infolder}${prefix}${chr}_p1.bgen ] #chr has segments
#     then
#       echo $chr
#       segs=()
#       for seg in ${infolder}${prefix}${chr}_p*.bgen
#       do
#         seg=$(basename ${seg})
#         seg=${seg#${prefix}${chr}_p}
#         seg=${seg%.bgen}
#         echo $seg
#         segs+=($seg)
#       done
#       sorted=($(printf '%s\n' "${segs[@]}"|sort -n))
#       echo ${sorted[@]}
#       for seg in ${sorted[@]}
#       do
#         echo $seg
#         #$plink2 --bgen ${infolder}${prefix}${chr}_p${seg}.bgen ref-unknown --sample ${infolder}${prefix}1_p1.sample --oxford-single-chr ${chr} --maf 0.005 --geno 0.1 --keep ${infolder}${prefix}plinksamples.txt --pheno ${infolder}${prefix}pheno.txt --hwe 1e-50 midp keep-fewhet --make-pgen --out ${outfolder}${outprefix}_chr.${chr}.${seg} --memory 300000 --threads 45
#         $plink2 --bgen ${infolder}${prefix}${chr}_p${seg}.bgen ref-unknown --sample ${infolder}${prefix}.sample --oxford-single-chr ${chr} --keep ${infolder}${prefix}plinksamples.txt --pheno ${infolder}${prefix}pheno.txt --export vcf vcf-dosage=DS --const-fid --out ${outfolder}${outprefix}_chr.${chr}.${seg} --memory 300000 --threads 45
#         #fix missing dosage:1/1: to 1/1:.
#         sed -i -E 's/:\s/:.\t/g' ${outfolder}${outprefix}_chr.${chr}.${seg}.vcf
#         sed -i -E 's/:$/:./g' ${outfolder}${outprefix}_chr.${chr}.${seg}.vcf
#         bcftools view ${outfolder}${outprefix}_chr.${chr}.${seg}.vcf -o ${outfolder}${outprefix}_chr.${chr}.${seg}.bcf
#         rm ${outfolder}${outprefix}_chr.${chr}.${seg}.vcf
#         echo ${outfolder}${outprefix}_chr.${chr}.${seg} >> ${outfolder}${outprefix}_mergelist.txt
#       done
#     else #no segments
#       #$plink2 --bgen ${infolder}${prefix}${chr}.bgen ref-unknown --sample ${infolder}${prefix}1_p1.sample --oxford-single-chr ${chr} --maf 0.005 --geno 0.1 --keep ${infolder}${prefix}plinksamples.txt --pheno ${infolder}${prefix}pheno.txt --hwe 1e-50 midp keep-fewhet --make-pgen --out ${outfolder}${outprefix}_chr.${chr} --memory 300000 --threads 45
#       $plink2 --bgen ${infolder}${prefix}${chr}.bgen ref-unknown --sample ${infolder}${prefix}.sample --oxford-single-chr ${chr} --maf 0.005 --keep ${infolder}${prefix}plinksamples.txt --pheno ${infolder}${prefix}pheno.txt --export vcf vcf-dosage=DS --const-fid --out ${outfolder}${outprefix}_chr.${chr} --memory 300000 --threads 45
#       #fix missing dosage:1/1: to 1/1:.
#       sed -i -E 's/:\s/:.\t/g' ${outfolder}${outprefix}_chr.${chr}.vcf
#       sed -i -E 's/:$/:./g' ${outfolder}${outprefix}_chr.${chr}.vcf
#       bcftools view ${outfolder}${outprefix}_chr.${chr}.vcf -o ${outfolder}${outprefix}_chr.${chr}.bcf
#       rm ${outfolder}${outprefix}_chr.${chr}.vcf
#       echo ${outfolder}${outprefix}_chr.${chr} >> ${outfolder}${outprefix}_mergelist.txt
#     fi
#   done
#   echo "merge data"
#   #$plink2 --pmerge-list ${outfolder}${outprefix}_mergelist.txt --make-pgen --out ${outfolder}${outprefix} --memory 300000 --threads 45
#   allfiles=($(cat ${outfolder}${outprefix}_mergelist.txt))
#   allfiles=( "${allfiles[@]/%/.bcf}" )
#   bcftools concat ${allfiles[@]} -Oz -o ${outfolder}${outprefix}.vcf.gz
#   $plink2 --vcf ${outfolder}${outprefix}.vcf.gz --make-pgen --out ${outfolder}${outprefix} --memory 300000 --threads 45
# }
# infolder="/data/BB_Bioinformatics/ProjectData/BCAC/icogs/"
# prefix="zhang_750_euro_icogs_topmed_"
# outfolder="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/"
# #outfolder="/data/DCEGLeiSongData/Kevin/BCAC/result/imp_icogs/"
# outprefix="euro"
# # $plink2 --pfile ${outfolder}${outprefix} --chr 22 --export A-transpose --out ${outfolder}${outprefix}_chr22 --memory 300000 --threads 45
# # gzip ${outfolder}${outprefix}_chr22.traw
# 
# infolder="/data/BB_Bioinformatics/ProjectData/BCAC/onco/"
# prefix="zhang_750_topmed_"
# #outfolder="/data/DCEGLeiSongData/Kevin/BCAC/result/imp_onco/"
# outfolder="/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/"
# outprefix="euro"
#EURO
read_seg "/data/BB_Bioinformatics/ProjectData/BCAC/icogs/" "zhang_750_euro_icogs_topmed_" "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/" "euro"
read_seg "/data/BB_Bioinformatics/ProjectData/BCAC/onco/" "zhang_750_topmed_" "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/" "euro"
#Asian
read_seg "/data/BB_Bioinformatics/ProjectData/BCAC/icogs/" "zhang_750_asian_icogs_topmed_" "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/" "asian"
read_seg "/data/BB_Bioinformatics/ProjectData/BCAC/onco/" "zhang_750_asian_topmed_" "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/" "asian"
#African
read_seg "/data/BB_Bioinformatics/ProjectData/BCAC/icogs/" "zhang_750_african_icogs_topmed_" "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/" "african"
#read_seg "/data/BB_Bioinformatics/ProjectData/BCAC/onco/" "zhang_750_african_topmed_" "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/" "african"
#Hispanic
#read_seg "/data/BB_Bioinformatics/ProjectData/BCAC/icogs/" "zhang_750_hispanic_icogs_topmed_" "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/hispanic/" "hispanic"
read_seg "/data/BB_Bioinformatics/ProjectData/BCAC/onco/" "zhang_750_hispanic_topmed_" "/data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/hispanic/" "hispanic"

$plink2 --pfile ${outfolder}${outprefix} --chr 22 --export A-transpose --out ${outfolder}${outprefix}_chr22 --memory 64000 --threads 8
gzip ${outfolder}${outprefix}_chr22.traw

# #for icogs
# $plink2 --pfile /data/DCEGLeiSongData/Kevin/BCAC/result/imp_icogs/euro --keep ../result/icogs_keepsamples.txt --make-pgen --out /data/DCEGLeiSongData/Kevin/BCAC/result/imp_icogs/euro --threads 8 --memory 56000

runPCA(){
  local inputfile=$1
  local outprefix=$2
  outfolder=$(dirname $inputfile)
  echo $outfolder
  #$plink2 --pfile $inputfile --maf 0.01 --geno 0.05  --make-bed --out ${outfolder}/$outprefix --threads 48 --memory 300000

  # First, we need to perform prunning 
  $plink2 \
    --pfile ${outfolder}/$outprefix \
    --indep-pairwise 1000 50 0.05 \
    --exclude range /data/BB_Bioinformatics/Kevin/tools/other/exclusion_regions_hg38.txt\
    --threads 48\
    --memory 300000\
    --out ${outfolder}/$outprefix
  
  $plink2 \
  --pfile ${outfolder}/$outprefix \
  --extract ${outfolder}/$outprefix.prune.in \
  -threads 48\
  --memory 300000\
  --make-bed --out ${outfolder}/${outprefix}_pruned 
  
  # Then we calculate the first 20 PCs
  $plink2 \
    --bfile ${outfolder}/$outprefix \
    --extract ${outfolder}/$outprefix.prune.in \
    --pca 20 \
    --threads 48\
    --memory 300000\
    --out ${outfolder}/$outprefix
}
runPCA /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro euro
runPCA /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro euro

