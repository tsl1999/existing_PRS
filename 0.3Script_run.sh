#!/bin/bash
#SBATCH --job-name="analysis"
#SBATCH --mem=20GB
#SBATCH --output=analysis.out

echo start

cd /well/emberson/users/hma817/projects/existing_PRS/analysis_code_correction/80
  echo $PWD

module load Anaconda3/2022.05 R/4.2.1-foss-2022a


#Rscript 3.2Datapreprocessing_07Aug2023.R
Rscript 4.2.PrimaryAnalysis_80_17Apr2024.R
Rscript 4.3.PrimaryAnalysis_decile_death80_17Apr2024.R 
Rscript 5.1.CAD_Mortality_80_17Apr2024.R 
Rscript 5.2.CAD_Mortality_decile_17Apr2024.R 
# Rscript 6.1.Stratified_by_sex_log_EPA_80_17Apr2024.R
# Rscript 6.2.Stratified_by_sex_log_EPO_80_17Apr2024.R 
Rscript 6.3.Baseline_CAD_80_17Apr2024.R 
Rscript 6.4.Stratified_by_conf_80_17Apr2024.R 
Rscript 6.5.Third_degree_unrelated_80_17Apr2024.R    
Rscript 6.6.Mediator_analyses_80.R 
Rscript 6.7.Age_over_80_17Apr2024.R

echo done