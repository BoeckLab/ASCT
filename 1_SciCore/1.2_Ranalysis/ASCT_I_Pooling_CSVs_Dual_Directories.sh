#!/bin/bash
#SBATCH --job-name=ASCT_I_PoolingCSVs       
#SBATCH --cpus-per-task=2                 
#SBATCH --mem-per-cpu=3G             
#SBATCH --time=00:20:00     
#SBATCH --qos=30min           
#SBATCH --output=%j.o     
#SBATCH --error=%j.e
#SBATCH --array=1-11340


#Loading R app
#module load R/4.0.3-foss-2018b
module load R



# ==================================================================
# Script for Pooling CSV Files from High-Throughput Image Analysis
# ==================================================================
# This script is designed to automate the organization of image 
# analysis results in a high-performance computing (HPC) cluster.
# It runs as a SLURM job array, processing multiple wells and fields
# in parallel.
#
# Functionality:
# - Loads required computing resources and R environment.
# - Extracts metadata (experiment ID, well, field, etc.) from 
#   a parameter file.
# - Identifies and copies CSV data files from image analysis output
#   to a structured destination folder.
# - If enabled, also processes BaSic Viability Object Classification analysis results.
#
# Usage:
# - Modify the experiment ID (ASCTexpID) before running.
# - Ensure the correct SLURM array size based on the number of wells
#   and fields to be processed.
#
# This script is essential for organizing large-scale image analysis 
# data, making it ready for further analysis in R.


#----- SECTION 1 : Preparing variables

# Please edit Number of SBATCH --array=1-11340 (number of wells x number of fields ) 

# Please edit the input variable ASCTexpID with the experiment ID forexample ASCT.08.20230504

ASCTexpID="ASCT.05.20230413"

ExpN="05"
# Please edit the input variable BasicDefault. If the BaSic default POC analysis was done.
# Options are either "Yes" or "No"
BaSicDefault="Yes"
BaSic_dir="_BaSicDefault"




path_main="/scicore/projects/rinfsci/Jovanovic/ASCT_Experiments/ASCT."$ExpN"/"$ASCTexpID"_directories" # image analysis results (were all the csv files are found)
destination_main="/scicore/projects/rinfsci/Jovanovic/TLKK_R_analysis/ASCT_moc.poc_simplified_multipleBaSic/"$ASCTexpID # Where they will be pooled 

parameters_txt_file=$ASCTexpID"_perField_parameters.txt"


# Copy perField parameters text file in the path_main directory. Thats where all the image analysis results are found



cd $path_main

# experiment Name e.g. TLCI
experimentName=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $1}')

# experiment ID e.g 13
experimentID=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $2}')

# experiment Data TLCI. YYMMDD
experimentDate=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $3}')

# Well coordinate
well=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $4}') 

# Field coordinate
field=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $5}') 

path1=$(pwd)


path=$path1"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_"$field

#path to destination directory
well_dir=$destination_main"/"$experimentName"."$experimentID"."$experimentDate"_"$well




#rm -r $outputName
cp $path/*.csv $well_dir/.



# Check the value of BaSicDefault and execute code accordingly
if [ "$BaSicDefault" == "Yes" ]; then
  path_main_BaSic="$path_main$BaSic_dir"
  path_BaSic=$path_main_BaSic"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_"$field
  # rm -r $outputName
  cp $path_BaSic/*.csv $well_dir/.
fi







