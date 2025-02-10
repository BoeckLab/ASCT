#!/bin/bash
#SBATCH --job-name=ASCT_II_R_1_2_3_4     
#SBATCH --cpus-per-task=2                 
#SBATCH --mem-per-cpu=8G             
#SBATCH --time=00:10:00     
#SBATCH --qos=30min           
#SBATCH --output=%j.o     
#SBATCH --error=%j.e
#SBATCH --array=1-1260


#Loading R app
#module load R/4.0.3-foss-2018b
module load R

# ==================================================================
# Script for Running R Analysis on High-Throughput Image Data
# ==================================================================
# This script automates the processing and analysis of image data 
# using R scripts in a high-performance computing (HPC) cluster.
# It runs as a SLURM job array, executing multiple tasks in parallel.
#
# Functionality:
# - Loads required computing resources and the R environment.
# - Extracts metadata (experiment ID, well, etc.) from a parameter file.
# - Copies R analysis scripts into the relevant directories.
# - Executes multiple R scripts sequentially for:
#   1. Ordering CSV files into structured directories.
#   2. Performing population analysis across multiple parameters.
#   3. Analyzing bright-field background intensity.
#   4. Analyzing fluorescence background intensity.
# - Pools the final results into designated result directories.
#
# Usage:
# - Modify the experiment ID (ASCTexpID) before running.
# - Ensure the correct SLURM array size based on the number of wells
#   to be processed.
#
# This script is essential for efficiently processing and analyzing
# high-throughput image data in a structured and automated manner.

#----- SECTION 1 : Preparing variables
# EDIT 1) Number of SBATCH --array=1-1260 & 2) paths

# Please edit Number of SBATCH --array=1-11340 (number of wells x number of fields ) 
# Please edit the input variable ASCTexpID with the experiment ID forexample ASCT.08.20230504

ASCTexpID="ASCT.05.20230413"
ExpN="05"

# Please edit the input variable BasicDefault. If the BaSic default POC analysis was done.
# Options are either "Yes" or "No"
BaSicDefault="Yes"

# Please define the R scripts which will be used
Rscript1="1_R_Ordering_Script.R" # Ordering all CSV files in their respective sub directories.
Rscript2="2_R_Population.Analysis.R" #	R population analysis scripts, across both BaSic parameters: Artefacts, growth, detachment, contamination etc
Rscript3="3_R_Brightfield.Background.Analysis.R" #	R population analysis of bright-field background intensity
Rscript4="4_R_uncorr.Fluorescence.Background.Analysis.R" #	R population analysis of Fluorescence background intensity



#----- SECTION 2 :Soft coded variables

# Soft coded variables below. However note the general paths should not change.

#CHECK BELOW path 
path_main="/scicore/projects/rinfsci/Jovanovic/TLKK_R_analysis/ASCT_moc.poc_simplified_multipleBaSic/"$ASCTexpID # image analysis results (were all the csv files are found)

#parameters_text_dir="/scicore/projects/rinfsci/Jovanovic/TLKK_R_analysis/ASCT_moc.poc_simplified_multipleBaSic/ASCT_Parameters_txt/"$ASCTexpID"_perWell_parameters.txt"
parameters_txt_file=$ASCTexpID"_perWell_parameters.txt"


cd $path_main

#cp $parameters_text_dir $path_main/.

# experiment Name e.g. TLCI
experimentName=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $1}')

# experiment ID e.g 13
experimentID=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $2}')

# experiment Data TLCI. YYMMDD
experimentDate=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $3}')

# Well coordinate
well=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $4}') 

#outputName=$experimentName"."$experimentID"."$experimentDate"/"$experimentName"."$experimentID"."$experimentDate"."$well


path1=$(pwd)
#path=$path1"/"$experimentName"."$experimentID"."$experimentDate"/"$experimentName"."$experimentID"."$experimentDate"."$well
path=$path1"/"$experimentName"."$experimentID"."$experimentDate"_"$well

# List list of R script that will be run
path_ofR1=$path1"/"$Rscript1
path_ofR2=$path1"/"$Rscript2
path_ofR3=$path1"/"$Rscript3
path_ofR4=$path1"/"$Rscript4



#rm -r $outputName
cp $path_ofR1 $path/
cp $path_ofR2 $path/
cp $path_ofR3 $path/
cp $path_ofR4 $path/


#----- SECTION 3: Applying R analysis script in well directories
cd $path


Rscript "$Rscript1"

echo "Script 1: Ordering CSV completed"

find -size 0 -delete
echo "Zero Bytes files deleted"

Rscript "$Rscript2"

echo "Script 2: Population analysis completed"

Rscript "$Rscript3"
echo "Script 3:  Brightfield backgroung analysis completed"

 Rscript "$Rscript4"
echo "Script 4: Fluorescence background analysis completed"


#----- SECTION 4: Pooling results in results sub directory

echo " CONGRATULATIONS ANALYSIS is complete "

echo "Preparing to copy results into respective results directory"

echo "Pooling Populationg analysis results"

cp $path"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_Results"/*_results.csv $path_main"/PerWell-Pop.Results"/.

echo "The Population analysis results is pooled"


cp $path"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_Results"/ASCT_BFint_results*.csv $path_main"/PerWell-BFInt.Results"/.

echo "The brightfield background analysis results is pooled"

cp $path"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_Results"/ASCT_FLuncorr_results_*.csv $path_main"/PerWell-FLInt.Results"/.

echo "The uncorrected background fluorescence analysis results is pooled"

echo " CONGRATULATIONS the POOLING step is complete "







