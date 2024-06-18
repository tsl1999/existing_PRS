#!/bin/bash
#SBATCH --job-name="phenotype"
#SBATCH --mem=5GB
#SBATCH --output=phenotype.out

echo start extracting phenotype

cd /well/emberson/users/hma817/projects/existing_PRS/
echo $PWD

module load R/3.6.2-foss-2019b
source /well/emberson/users/hma817/projects/existing_PRS/env_prs/bin/activate

Rscript /gpfs3/well/emberson/shared/workflows/gwas/topmed-imputed/01.1_extract-phenotype-info.R correction/data SEX,FEMALE,AGE,SMOKE,INCOME,\
BMI,BASE_CHD,BASE_CVD,SBP,DBP,EDU_LEVEL,smokegp,smokegp2,EDU_UNI,COYOACAN,BASE_DIABETES,WAISTC,HIPC,\
BASE_CANCER,BASE_EMPHYSEMA,BASE_CIRR,BASE_PEP,BASE_CKD,BASE_PAD,ALCGP,\
WHRATIO,BASE_HBA1C,DRUG_D1,DRUG_D2,DRUG_D3,DRUG_D4,HDL_C,LDL_C,EPA001,EPA001A,\
EPO001,EPO001A,ICD10_OXFORD,DATE_RECRUITED,DATE_OF_DEATH


echo done
