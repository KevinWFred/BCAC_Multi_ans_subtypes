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


#compute AF
$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro --freq --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro --memory 200000 --threads 10
$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian --freq --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian --memory 100000 --threads 10 &
$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african --freq --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african --memory 50000 --threads 10

#merge icogs samples to run PCA
$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro --export bcf vcf-dosage=DS --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro
tabix /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro.bcf

$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian --export bcf vcf-dosage=DS --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian
tabix /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian.bcf

$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african --export bcf vcf-dosage=DS --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african
tabix /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african.bcf

# sed 's/^>chr/>/' /data/BB_Bioinformatics/Kevin/tools/database/hg38.fa > /data/BB_Bioinformatics/Kevin/tools/database/hg38_nochr.fa
# bcftools norm --check-ref s -f /data/BB_Bioinformatics/Kevin/tools/database/hg38_nochr.fa /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african.bcf -o /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african1.bcf
# bcftools stats /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african1.bcf > /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african1.txt
# tabix /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african1.bcf

#bcftools norm --check-ref s -f /data/BB_Bioinformatics/Kevin/tools/database/human_g1k_v37.fasta /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african.bcf -o /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african2.bcf

#merge based on IDs
bcftools merge -m id /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro.bcf \
 /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian.bcf \
 /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african.bcf --threads 10 -o /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged.bcf

$plink2 -bcf /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged.bcf dosage=DS --make-pgen --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged --memory 200000 --threads 10

#merge onco samples to run pooled analysis
$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro --export bcf vcf-dosage=DS --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro
tabix /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro.bcf

$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian --export bcf vcf-dosage=DS --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian
tabix /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian.bcf

$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african --export bcf vcf-dosage=DS --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african
tabix /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african.bcf

$plink2 --pfile /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/hispanic/hispanic --export bcf vcf-dosage=DS --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/hispanic/hispanic
tabix /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/hispanic/hispanic.bcf

#merge based on IDs
bcftools merge -m id /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro/euro.bcf \
 /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/asian/asian.bcf \
 /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/hispanic/hispanic.bcf \
 /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/african/african.bcf --threads 10 -o /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged.bcf

$plink2 -bcf /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged.bcf dosage=DS --make-pgen --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged --memory 200000 --threads 10

runPCA(){
  local inputfile=$1
  local outprefix=$2
  outfolder=$(dirname $inputfile)
  echo $outfolder
  $plink2 --pfile $inputfile --maf 0.05 --geno 0.05  --make-pgen --out ${outfolder}/${outprefix}_MAF05 --threads 12 --memory 200000

  # First, we need to perform prunning 
  $plink2 \
    --pfile ${outfolder}/${outprefix}_MAF05 \
    --indep-pairwise 1000 50 0.05 \
    --exclude range /data/BB_Bioinformatics/Kevin/tools/other/exclusion_regions_hg38.txt\
    --threads 12\
    --memory 200000\
    --out ${outfolder}/$outprefix
  
  $plink2 \
  --pfile ${outfolder}/$outprefix \
  --extract ${outfolder}/$outprefix.prune.in \
  -threads 12\
  --memory 200000\
  --make-bed --out ${outfolder}/${outprefix}_pruned 
  
  $plink2 \
  --pfile ${outfolder}/$outprefix \
  --extract ${outfolder}/$outprefix.prune.in \
  -threads 12\
  --memory 200000\
  --make-pgen --out ${outfolder}/${outprefix}_pruned 
  
  # Then we calculate the first 10 PCs
  $plink2 \
    --pfile ${outfolder}/$outprefix \
    --extract ${outfolder}/$outprefix.prune.in \
    --pca 10 \
    --threads 12\
    --memory 200000\
    --out ${outfolder}/$outprefix
}
runPCA /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/euro euro
runPCA /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro euro
runPCA /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1 merged1
runPCA /data/BB_Bioinformatics/Kevin/BCAC/result/imp_onco/merged1 merged1

#merge the whole datasets to get PCs
$plink2 --pfile ../result/imp_onco/merged1_MAF05 --make-bed --out /data/DCEGLeiSongData/Kevin/tmp/onco_MAF05 --memory 200000 --threads 10
$plink2 --pfile ../result/imp_icogs/merged1_MAF05 --make-bed --out /data/DCEGLeiSongData/Kevin/tmp/icogs_MAF05 --memory 200000 --threads 10
echo "/data/DCEGLeiSongData/Kevin/tmp/onco_MAF05"> /data/DCEGLeiSongData/Kevin/tmp/BCACmergelist.txt
echo "/data/DCEGLeiSongData/Kevin/tmp/icogs_MAF05">> /data/DCEGLeiSongData/Kevin/tmp/BCACmergelist.txt
$plink --merge-list /data/DCEGLeiSongData/Kevin/tmp/BCACmergelist.txt --geno 0.05 --make-bed --out /data/DCEGLeiSongData/Kevin/tmp/BCAC --memory 20000 --threads 5
$plink2 \
    --bfile /data/DCEGLeiSongData/Kevin/tmp/BCAC \
    --indep-pairwise 1000 50 0.05 \
    --exclude range /data/BB_Bioinformatics/Kevin/tools/other/exclusion_regions_hg38.txt\
    --threads 12\
    --memory 20000\
    --out /data/DCEGLeiSongData/Kevin/tmp/BCAC
$plink2 \
  --bfile /data/DCEGLeiSongData/Kevin/tmp/BCAC \
  --extract /data/DCEGLeiSongData/Kevin/tmp/BCAC.prune.in \
  -threads 12\
  --memory 20000\
  --make-bed --out /data/DCEGLeiSongData/Kevin/tmp/BCAC_pruned 
$plink2 \
    --bfile /data/DCEGLeiSongData/Kevin/tmp/BCAC_pruned \
    --pca 10 \
    --threads 40\
    --memory 400000\
    --out ../result/BCAC