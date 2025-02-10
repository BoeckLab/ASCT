#!/bin/bash

#SBATCH --job-name=imaging       #Relate the name of your job with the plate/experiment ID
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved for every task
#SBATCH --mem-per-cpu=20G             #This is the memory reserved per core for every task (might need to be adapted)
#SBATCH --tmp=20G 
 
#SBATCH --time=05:30:00         #This is the time that your task will run (better to adapt when better known)
#SBATCH --qos=6hours
        

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=%j.o     #These are the STDOUT and STDERR files
#SBATCH --error=%j.e


#This job runs from the current working directory

#load your required modules below
#################################


#export your required environment variables below
#################################################


#add your command lines below
#############################
module purge
module load Fiji/20201104-1356 
module load numpy/1.10.1-goolf-1.7.20-Python-2.7.11 

# This SLURM job script automates image analysis by processing .nd2 files, 
# running Fiji, for background correction, metadata extraction and Ilastik pipelines for morphology, viability classification and background assessment, 
# generating H5 files, and extracting key results. It manages data paths, 
# logs resource usage, and exports outputs for further analysis. 

#DEFINE PATHS
fiji_macroBaSicDef="/scicore/projects/rinfsci/Jovanovic/ILASTIK_FIJI_ALGORITHMS/FIJI_scripts/Fiji_TLCI_moc.poc/ASCT_SINGLETIMELAPSE_BaSicDEFAULT_moc.poc_20220815.groovy" # BaSic with Default settings i.e automatic
fiji_macroBaSicFF9FDF0="/scicore/projects/rinfsci/Jovanovic/ILASTIK_FIJI_ALGORITHMS/FIJI_scripts/Fiji_TLCI_moc.poc/TLKK48_SINGLETIMELAPSE_QC_BaSic_MD_moc.poc_20220915_v0.groovy" # BaSic with Default settings i.e automatic
fiji_macroFLuncorr="/scicore/projects/rinfsci/Jovanovic/ILASTIK_FIJI_ALGORITHMS/FIJI_scripts/Fiji_TLCI_moc.poc/ASCT_SINGLETIMELAPSE_FLuncorr_Only.groovy"


project_PC="/scicore/projects/rinfsci/Jovanovic/ILASTIK_FIJI_ALGORITHMS/ILASTIK-algorithms/MOC.POC/PC_bf_v1_20220502.ilp"

# PI based object classification
project_POC="/scicore/projects/rinfsci/Jovanovic/ILASTIK_FIJI_ALGORITHMS/ILASTIK-algorithms/MOC.POC/POC_hyst_v0_20220608_v1_ASCT_FINAL_VERSION.ilp"

project_MOC="/scicore/projects/rinfsci/Jovanovic/ILASTIK_FIJI_ALGORITHMS/ILASTIK-algorithms/MOC.POC/MOC_20221216_ASCT_NewCAMERA_v0.ilp"



#REPRESENTS BaSic skipped Algorithm::: Bright-field background intensity
project_BFbckInt="/scicore/projects/rinfsci/Jovanovic/ILASTIK_FIJI_ALGORITHMS/ILASTIK-algorithms/BOC/ASCT_BF_BckgInf_v0.ilp"

# Background object classification of the rawfL (uncorrected by BaSic) 
project_BOC="/scicore/projects/rinfsci/Jovanovic/ILASTIK_FIJI_ALGORITHMS/ILASTIK-algorithms/BOC/BOC_v1.ilp"



experimentName=$(head -$SLURM_ARRAY_TASK_ID parameters.txt | tail -1 | awk '{print $1}')
experimentID=$(head -$SLURM_ARRAY_TASK_ID parameters.txt | tail -1 | awk '{print $2}')
experimentDate=$(head -$SLURM_ARRAY_TASK_ID parameters.txt | tail -1 | awk '{print $3}')
loop=$(head -$SLURM_ARRAY_TASK_ID parameters.txt | tail -1 | awk '{print $4}')
well=$(head -$SLURM_ARRAY_TASK_ID parameters.txt | tail -1 | awk '{print $5}')
field=$(head -$SLURM_ARRAY_TASK_ID parameters.txt | tail -1 | awk '{print $6}')


output_name=$experimentName"."$experimentID"."$experimentDate"_"$well"_"$field
output_name2=$experimentName"."$experimentID"."$experimentDate"_"$well"_"$field

#echo $output_name2




path=$(pwd | rev | cut -d "_" -f 2- | cut -d "/" -f 2- | rev)

pathA=$path"/"$experimentName"."$experimentID"."$experimentDate"_data/"$experimentName"."$experimentID"."$experimentDate"_data_LoopAt"


output_directory=$path"/"$experimentName"."$experimentID"."$experimentDate"_directories/"$output_name

#Djava.io.tmpdir=$TMPDIR
mkdir $output_name
cd $TMPDIR
ls *

echo "fiji code: "$fiji_macroBaSicDef
echo "Brightfield h5 data path: "$BF_data_path

echo "Corrected Fluorescence h5 data path: "$Cy3corr_data_path

echo "output_name: "$output_name
echo "output_name2: "$output_name2
echo "path_nameA: "$pathA
echo "output directory: "$output_directory



#MOVE ND2 FILES INTO DIRECTORY
path_fileA=$pathA"/"$experimentName"."$experimentID"."$experimentDate"*_"$well"_"$field"*nd2"

echo "path_pathfileA: "$path_fileA

cp $path_fileA $TMPDIR/


echo "Files copied"
ls *




#GENERATE H5 FILE BaSic default settings
ImageJ-linux64 --ij2 --headless --run $fiji_macroBaSicDef "input_path='',export_path=''"
echo "BF.H5 & FLcorrBaSicdef.h5 files generated"
sstat --format=JobID,AveCPU,AveRSS,MaxRSS -j  $SLURM_JOBID.batch
sacct -o JobID,CPUTime -j $SLURM_JOBID
ls *


#GENERATE H5 FILE BaSic FF9 DF 0
ImageJ-linux64 --ij2 --headless --run $fiji_macroBaSicFF9FDF0 "input_path='',export_path=''"
echo "BF.H5 & FLcorr.H5 files generated"
sstat --format=JobID,AveCPU,AveRSS,MaxRSS -j  $SLURM_JOBID.batch
sacct -o JobID,CPUTime -j $SLURM_JOBID
ls *

#GENERATE H5 FILE for Backgroung assessment

ImageJ-linux64 --ij2 --headless --run $fiji_macroFLuncorr "input_path='',export_path=''"
echo "rawFL.H5 file generated"
sstat --format=JobID,AveCPU,AveRSS,MaxRSS -j  $SLURM_JOBID.batch
sacct -o JobID,CPUTime -j $SLURM_JOBID
ls *


#RENAME H5 FILE
rename "__Channel_x40BF-x40Cy3_BFCy3.h5" "_BFFL.h5" *
rename "__Channel_x40BF-x40Cy3_BF.h5" "_BF.h5" *
rename "__Channel_x40BF-40xSytG_BFCy3.h5" "_BFFL.h5" *
rename "__Channel_x40BF-40xSytG_BF.h5" "_BF.h5" *
rename "__Channel_x40BF-40xSytB_BFCy3.h5" "_BFFL.h5" *
rename "__Channel_x40BF-40xSytB_BF.h5" "_BF.h5" *
rename "__Channel_x40BF-x40Cy3_0001_BFCy3.h5" "_BFFL.h5" *
rename "__Channel_x40BF-x40Cy3_0001_BF.h5" "_BF.h5" *
rename "__Channel_40x_BF-40x_Cy3_BF.h5" "_BF.h5" *
rename "__Channel_40x_BF-40x_Cy3_BFCy3.h5" "_BFFL.h5" *
rename "__Channel_40xBF-40xCy3_BFCy3.h5" "_BFFL.h5" *
rename "__Channel_40xBF-40xCy3_BF.h5" "_BF.h5" *
rename "_BFCy3.h5" "_BFFL.h5" *h5 
rename "_Channel_x40BF,x40Cy3_BF.h5" "_BF.h5" *h5
rename "_Channel_x40BF,x40Cy3_BF.h5" "_BFFL.h5" *h5




BF_data_path=$output_name2"_BF.h5"
POCbSd_data_path=$output_name2"_FLcorrBaSicdef.h5"
POC_data_path=$output_name2"_FLcorr.h5"

MOC_data_path=$output_name2"_BF.h5"
BOC_data_path=$output_name2"_raw_FLcorr.h5"


echo $BF_data_path
#echo $BFFL_data_path
echo $MOC_data_path

#GENERATE PIXEL CLASSIFICATION MAP
~/ilastik-1.3.3post3-Linux/run_ilastik.sh \
--headless \
--project $project_PC \
--raw_data $BF_data_path \
--output_filename_format PC-$output_name.h5 \
--export_source "Probabilities"
echo "Pixel classification finished"
sstat --format=JobID,AveCPU,AveRSS,MaxRSS -j  $SLURM_JOBID.batch
sacct -o JobID,CPUTime -j $SLURM_JOBID
ls $TMPDIR/*


echo "PI BaSic default classification STARTED"

#RUN PI OBJECT CLASSIFICATION
~/ilastik-1.3.3post3-Linux/run_ilastik.sh \
--headless \
--project $project_POC \
--raw_data $POCbSd_data_path \
--prediction_maps PC-$output_name.h5 \
--table_filename POCbSd-$output_name.csv \
--export_source "Object Identities" 
echo "PI (BaSic default)Object classification finished"
cat $TMPDIR/POCbSd-*csv | ~/.local/bin/csvcut -c "timestep","labelimage_oid","Predicted Class","Center of the object_0","Center of the object_1","Mean Intensity" > OCsimple_POCbSd_$output_name.csv
#cat $TMPDIR/OCsimple*.csv | sed 's/cluster_2-5_cells/C2/g' | sed 's/false_cell/FC/g' |sed 's/single_cell/SC/g'| sed 's/v-snapp/VS/g' |sed 's/cluster_5-20_cells/C5/g' |sed 's/cluster_20+_cells/C20/g' | tr '.' ',' |cut -f 1,3,4,6 -d , |  sed 's/Mean Intensity_1/Mean Intensity_1,Size in pixels/g' > OCsimple2_$output_name.csv
#rename "_00" "_p" *csv



echo "PI BaSic FF9DF0 classification STARTED"

#RUN PI OBJECT CLASSIFICATION
~/ilastik-1.3.3post3-Linux/run_ilastik.sh \
--headless \
--project $project_POC \
--raw_data $POC_data_path \
--prediction_maps PC-$output_name.h5 \
--table_filename POC-$output_name.csv \
--export_source "Object Identities" 
echo "PI (BaSic FF9DF0)Object classification finished"
cat $TMPDIR/POC-*csv | ~/.local/bin/csvcut -c "timestep","labelimage_oid","Predicted Class","Center of the object_0","Center of the object_1","Mean Intensity" > OCsimple_POC_$output_name.csv


echo "PI classification is COMPLETE"

echo "Morphology Object classification STARTED"
#RUN MORPHOLOGY OBJECT CLASSIFICATION
~/ilastik-1.3.3post3-Linux/run_ilastik.sh \
--headless \
--project $project_MOC \
--raw_data $MOC_data_path \
--prediction_maps PC-$output_name.h5 \
--table_filename MOC-$output_name.csv \
--export_source "Object Identities" 
echo "Object classification finished"
cat $TMPDIR/MOC-*csv | ~/.local/bin/csvcut -c "timestep","labelimage_oid","Predicted Class","Size in pixels","Center of the object_0","Center of the object_1","Mean Intensity" > OCsimple_MOC_$output_name.csv


echo "Morphology Object classification is COMPLETE"


echo "Bright-field Background Object classification STARTED"
#RUN BACKGROUND OBJECT CLASSIFICATION
~/ilastik-1.3.3post3-Linux/run_ilastik.sh \
--headless \
--project $project_BFbckInt \
--raw_data $MOC_data_path \
--prediction_maps PC-$output_name.h5 \
--table_filename BFbckInt_$output_name.csv \
--export_source "Object Identities" 
echo "Bright-field Background Object classification finished"
cat $TMPDIR/BFbckInt_*csv | ~/.local/bin/csvcut -c "timestep","labelimage_oid","Predicted Class","Center of the object_0","Center of the object_1","Size in pixels","Mean Intensity" > OCsimple_BFbckint_$output_name.csv

echo "Brightfield Background Object classification is COMPLETE"



echo "Background (rawFL) Object classification STARTED"

#RUN BACKGROUND OBJECT CLASSIFICATION
~/ilastik-1.3.3post3-Linux/run_ilastik.sh \
--headless \
--project $project_BOC \
--raw_data $BOC_data_path \
--prediction_maps PC-$output_name.h5 \
--table_filename BOC_$output_name.csv \
--export_source "Object Identities" 
echo "Background Object classification finished"
cat $TMPDIR/BOC_*csv | ~/.local/bin/csvcut -c "timestep","labelimage_oid","Predicted Class","Center of the object_0","Center of the object_1","Mean Intensity" > OCsimple_BOC_$output_name.csv
#cat $TMPDIR/OCsimple*.csv | sed 's/cluster_2-5_cells/C2/g' | sed 's/false_cell/FC/g' |sed 's/single_cell/SC/g'| sed 's/v-snapp/VS/g' |sed 's/cluster_5-20_cells/C5/g' |sed 's/cluster_20+_cells/C20/g' | tr '.' ',' |cut -f 1,3,4,6 -d , |  sed 's/Mean Intensity_1/Mean Intensity_1,Size in pixels/g' > OCsimple2_$output_name.csv
#rename "_00" "_p" *csv

echo "Background (rawFL) Object classification is COMPLETE"


sstat --format=JobID,AveCPU,AveRSS,MaxRSS -j  $SLURM_JOBID.batch
sacct -o JobID,CPUTime -j $SLURM_JOBID
ls *

#GENERATE MEMORY TIME OUTPUT
echo "COMPLETELY FINISHED"
sstat --format=JobID,AveCPU,AveRSS,MaxRSS -j  $SLURM_JOBID.batch
sacct -o JobID,CPUTime -j $SLURM_JOBID

#MOVE LOG FILES
du -k * | sort -nr | cut -f2 | xargs -d '\n' du -sh
mv $path"/"$experimentName"."$experimentID"."$experimentDate"_directories/"$SLURM_JOB_ID.? $output_directory"/."


#EXPORT

rm $output_directory"/*csv"

cp POCbSd-*csv $output_directory"/."
cp OCsimple_POCbSd_*csv $output_directory"/."

cp POC-*csv $output_directory"/."
cp OCsimple_POC_*csv $output_directory"/."


cp *_BF.h5 $output_directory"/." # KEEP ME
#cp *.h5 $output_directory"/."
cp MOC-*Identities.h5 $output_directory"/." 


cp *csv $output_directory"/."
