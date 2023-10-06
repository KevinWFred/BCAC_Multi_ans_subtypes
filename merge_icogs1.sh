#!/usr/bin/env bash

cd /data/BB_Bioinformatics/Kevin/BCAC/code
plink=/usr/local/apps/plink/1.9.0-beta4.4/plink
plink2=/usr/local/apps/plink/2.3-alpha/plink2
ml samtools
date
bcftools merge -m id /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/euro/euro.bcf \
/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/asian/asian.bcf \
/data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/african/african.bcf --threads 24 -o /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1.bcf
date
echo "onco bcf done"
$plink2 -bcf /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1.bcf dosage=DS --make-pgen --out /data/BB_Bioinformatics/Kevin/BCAC/result/imp_icogs/merged1 --memory 120000 --threads 24
date
echo "onco plink2 done"
#5806097
#sbatch --mem=128g --cpus-per-task=12 --time=3-00:00:00 --gres=lscratch:64  /data/BB_Bioinformatics/Kevin/BCAC/code/merge_icogs1.sh
