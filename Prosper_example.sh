#!/usr/bin/bash

#add .libPaths(c("/data/BB_Bioinformatics/Kevin/tools/Rpackages",.libPaths())) in R files.

package='/data/BB_Bioinformatics/Kevin/tools/PROSPER'
path_example='/data/BB_Bioinformatics/Kevin/tools/PROSPER/example/'
path_result='/data/BB_Bioinformatics/Kevin/tools/PROSPER/PROSPER_example_results/'
path_plink='/usr/local/apps/plink/2.3-alpha/plink2'

# head -1 ${path_example}/sample_data/EUR/tuning_geno.fam
# 110001  110001  0       0       0       -9
# grep 110001 ${path_example}/sample_data/EUR/pheno.fam
# 110001 110001 0.202136

mkdir ${path_result}

date
Rscript ${package}/scripts/lassosum2.R \
--PATH_package ${package} \
--PATH_out ${path_result}/lassosum2 \
--PATH_plink ${path_plink} \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_example}/sample_data/EUR/tuning_geno,${path_example}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--bfile_testing ${path_example}/sample_data/EUR/testing_geno,${path_example}/sample_data/AFR/testing_geno \
--pheno_testing ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--testing TRUE \
--NCORES 8
date

Rscript ${package}/scripts/PROSPER.R \
--PATH_package ${package} \
--PATH_out ${path_result}/PROSPER \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--lassosum_param ${path_result}/lassosum2/EUR/optimal_param.txt,${path_result}/lassosum2/AFR/optimal_param.txt \
--chrom 1-22 \
--NCORES 8
date

target_pop='AFR'

Rscript ${package}/scripts/tuning_testing.R \
--PATH_plink ${path_plink} \
--PATH_out ${path_result}/PROSPER \
--prefix ${target_pop} \
--testing TRUE \
--bfile_tuning ${path_example}/sample_data/${target_pop}/tuning_geno \
--pheno_tuning ${path_example}/sample_data/${target_pop}/pheno.fam \
--bfile_testing ${path_example}/sample_data/${target_pop}/testing_geno \
--pheno_testing ${path_example}/sample_data/${target_pop}/pheno.fam \
--cleanup F \
--NCORES 8

target_pop='EUR'
Rscript ${package}/scripts/tuning_testing.R \
--PATH_plink ${path_plink} \
--PATH_out ${path_result}/PROSPER \
--prefix ${target_pop} \
--testing TRUE \
--bfile_tuning ${path_example}/sample_data/${target_pop}/tuning_geno \
--pheno_tuning ${path_example}/sample_data/${target_pop}/pheno.fam \
--bfile_testing ${path_example}/sample_data/${target_pop}/testing_geno \
--pheno_testing ${path_example}/sample_data/${target_pop}/pheno.fam \
--cleanup F \
--NCORES 8

