#!/bin/bash
#SBATCH --job-name=ASCT_III_R_5.2       
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
# Script for Tracking Data Merging & Time-Kill Curve Analysis
# ==================================================================
# This script automates the processing and analysis of tracking data 
# in high-throughput image analysis experiments using R scripts. 
# It runs as a SLURM job array, executing multiple tasks in parallel.
#
# Functionality:
# - Loads required computing resources and the R environment.
# - Extracts metadata (experiment ID, well, etc.) from a parameter file.
# - Copies tracking result files into the relevant well directory.
# - Executes an R script to:
#   1. Merge tracking data with POC analysis.
#   2. Perform time-kill curve analysis.
# - Organizes and pools results into designated directories for further analysis.
#
# Usage:
# - Modify the experiment ID (ASCTexpID) before running.
# - Ensure the correct SLURM array size based on the number of wells
#   to be processed.
#
# This script is essential for tracking-based cell population analysis,
# deriving time-kill curves, and structuring the output for downstream 
# data interpretation.
# ==================================================================

#----- SECTION 1 : Preparing variables
# EDIT 1) Number of SBATCH --array=1-4932 & 2) paths

# Please edit Number of SBATCH --array=1-11340 (number of wells x number of fields ) 

# Please edit the input variable ASCTexpID with the experiment ID forexample ASCT.08.20230504

ASCTexpID="ASCT.05.20230413"
ExpN="05"
TrackingVersion="Trkv2"

# Please edit the input variable BasicDefault. If the BaSic default POC analysis was done.
# Options are either "Yes" or "No"
BaSicDefault="Yes"

# Please define the R scripts which will be used

Rscript5="5.2_R_Single.Cell.Tracking.Intergration.R" # R tracking analysis; merges the POC data to the tracking information v5.2 is the new tracking output with missing frames and no piDynamics



#----- SECTION 2 :Soft coded variables

path_main="/scicore/projects/rinfsci/Jovanovic/TLKK_R_analysis/ASCT_moc.poc_simplified_multipleBaSic/"$ASCTexpID
trk_dir="/scicore/projects/rinfsci/Ahmad/MATLAB/PI/ASCT."$ExpN"_3"

parameters_text_dir="/scicore/projects/rinfsci/Jovanovic/TLKK_R_analysis/ASCT_moc.poc_simplified_multipleBaSic/ASCT_Parameters_txt/"$ASCTexpID"_perWell_parameters.txt"
parameters_txt_file=$ASCTexpID"_perWell_parameters.txt"

cd $path_main

# experiment Name e.g. TLCI
experimentName=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $1}')

# experiment ID e.g 13
experimentID=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $2}')

# experiment Data TLCI. YYMMDD
experimentDate=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $3}')

# Well coordinate
well=$(head -$SLURM_ARRAY_TASK_ID $parameters_txt_file | tail -1 | awk '{print $4}') 



#----- SECTION 2.1 : Copying tracking results to main R analysis directory

path1=$(pwd)

#path=$path1"/"$experimentName"."$experimentID"."$experimentDate"/"$experimentName"."$experimentID"."$experimentDate"."$well
path=$path1"/"$experimentName"."$experimentID"."$experimentDate"_"$well

# Giving myself rights
chmod u+rwx $trk_dir"/"$experimentName"."$experimentID"."$experimentDate"_"$well"/Dynamic_Table_corrected_"$experimentName"."$experimentID"."$experimentDate"_"$well".csv"


path_well_dir=$trk_dir"/"$experimentName"."$experimentID"."$experimentDate"_"$well"/Dynamic_Table_corrected_"$experimentName"."$experimentID"."$experimentDate"_"$well".csv"

# List list of R script that will be run
path_ofR5=$path1"/"$Rscript5



# Copy necessary files into specific well directory
cp $path_well_dir $path/
cp $path_ofR5 $path/

cd $path

mkdir $experimentName"."$experimentID"."$experimentDate"_"$well"_Tracking_labeloid"
mkdir $experimentName"."$experimentID"."$experimentDate"_"$well"_Tracking_LC"

mv "Dynamic_Table_corrected_"$experimentName"."$experimentID"."$experimentDate"_"$well".csv" "Dynamic_Table_corrected_"$experimentName"."$experimentID"."$experimentDate"_"$well"-"$TrackingVersion".csv"
mv "Dynamic_Table_corrected_"$experimentName"."$experimentID"."$experimentDate"_"$well"-"$TrackingVersion".csv"  $experimentName"."$experimentID"."$experimentDate"_"$well"_Tracking_labeloid"/

echo "Starting part ASCT part II: tracking analysis"
echo $experimentName"."$experimentID"."$experimentDate"_"$well

#----- SECTION 3: Applying R analysis script in well directories

# Merging should fail because Im using data from TLKK77 and merging it with TLKK 83 

#Rscript 5_R_Tracking_merge_LCF_20230609_FINAL_dplyr.R
Rscript "$Rscript5"

echo "Script 5 tracking data merging and time-kill curve derivation is complete"



#----- SECTION 4: Pooling results in results sub directory


# Tracking data pooling
echo "Preparing to copy results into respective results directory"

echo "Pooling Tracking time-kill curves analysis results"

cp $path"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_Tracking_LC"/*_Trkv2_LC.csv $path_main"/PerWell-Trk.Results"/.

echo "Tracking time-kill curves analysis results is pooled"

echo "Pooling tracking based time-kill curves results: cell numbers"

cp $path"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_Tracking_LC"/*_Trkv2_LC.N.csv $path_main"/PerWell-Trk.n.Results"/.

echo "Tracking time-kill curves cell numbers results is pooled"


# Killing features pooling
echo " CONGRATULATIONS the POOLING step is complete "

