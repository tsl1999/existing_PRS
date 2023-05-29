echo start
cd /well/emberson/users/hma817/projects/existing_PRS
echo $PWD

source /well/emberson/users/hma817/projects/existing_PRS/env_prs/bin/activate
module load Anaconda3/2022.05 PLINK/2.00a2.3_x86_64 Python/3.10.4-GCCcore-11.3.0 Java/11.0.2 R/4.2.1-foss-2022a yaml-cpp/0.7.0-GCCcore-11.3.0 
pip install pyyaml
for pgs_name in PGS000011 PGS000018 PGS000013 PGS003446 PGS001780;
 do 
 echo $pgs_name start computation
 mkdir PRS_calculated/$pgs_name
 ./nextflow run pgscatalog/pgsc_calc -profile conda\
 --input sample_sheet_UPDATE.csv\
 --pgs_id $pgs_name\
 --target_build GRCh38\
 --parallel\
 --outdir PRS_calculated/$pgs_name ;
 
 echo finished
 done
  
  echo all PRS in the list computed
  