#!/usr/bin/env bash

cd /data/BB_Bioinformatics/Kevin/BCAC/code
plink=/usr/local/apps/plink/1.9.0-beta4.4/plink
plink2=/usr/local/apps/plink/2.3-alpha/plink2
ml samtools

$plink2 \
--bfile /data/DCEGLeiSongData/Kevin/tmp/BCAC_pruned \
--pca 10 \
--threads 31 \
--memory 450000 \
--out ../result/BCAC
#sbatch --mem=480g --cpus-per-task=32 --time=7-00:00:00  /data/BB_Bioinformatics/Kevin/BCAC/code/plink2_BCACPCA.sh
