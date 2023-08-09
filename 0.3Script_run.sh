#!/bin/bash
#SBATCH --job-name="analysis"
#SBATCH --mem=20GB
#SBATCH --output=analysis.out

echo start extracting phenotype

cd /well/emberson/users/hma817/projects/existing_PRS/
  echo $PWD

module load R/4.2.1-foss-2022a
source /well/emberson/users/hma817/projects/existing_PRS/env_prs/bin/activate


Rscript 4.2PrimaryAnalysis.R 
Rscript 4.3PrimaryAnalysis_decile.R 
Rscript 4.4PrimaryAnalysis_int.R 
Rscript 5.1CAD_Mortality.R 
Rscript 5.2CAD_Mortality_decile.R 
Rscript 6.1Stratifieed_by_sex_log_EPA.R
Rscript 6.2Stratifieed_by_sex_log_EPO.R 
Rscript 6.3Baseline_CAD.R 
Rscript 6.4Age_over_75.R 
Rscript 6.4Stratified_by_conf.R 
Rscript 6.5Third_degree_unrelated.R    

echo done