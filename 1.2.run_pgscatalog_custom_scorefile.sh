
echo start
cd /well/emberson/users/hma817/projects/existing_PRS
echo $PWD

source /well/emberson/users/hma817/projects/existing_PRS/env_prs/bin/activate
module load Anaconda3/2022.05 PLINK/2.00a2.3_x86_64 Python/3.10.4-GCCcore-11.3.0 Java/11.0.2 R/4.2.1-foss-2022a yaml-cpp/0.7.0-GCCcore-11.3.0 
pip install pyyaml
for pgs_name in PGS000337 Oni-Orisan;
 do 
 echo $pgs_name start computation
 ./nextflow run pgscatalog/pgsc_calc -profile conda\
 --input sample_sheet_UPDATE.csv\
 --target_build GRCh38\
 --parallel\
 --outdir PRS_calculated/custom_$pgs_name\
 --scorefile PRS_calculated/custom_$pgs_name/$pgs_name-harHeader-GRCh38.txt;
 
 echo finished
 done
  
  echo all PRS in the list computed
  deactivate

#for some reason if letting the program do the liftover the matching variants will be very few (1%)
#I suspect there is some problem with the liftover function they used, therefore I performed manual liftover
#or if there is already harmonised file available, I used that file and change the build name to GRCh38
#all future PRS should adopt the same naming ruld, ie $pgs_name-harHeader-GRCh38 to indicate
# it is harmonised to GRCh38 genome build