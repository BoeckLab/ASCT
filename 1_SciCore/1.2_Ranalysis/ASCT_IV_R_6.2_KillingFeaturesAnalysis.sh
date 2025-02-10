#!/bin/bash
#SBATCH --job-name=ASCT_IV_R_6.2     
#SBATCH --cpus-per-task=2                 
#SBATCH --mem-per-cpu=8G             
#SBATCH --time=00:20:00     
#SBATCH --qos=30min           
#SBATCH --output=%j.o     
#SBATCH --error=%j.e
#SBATCH --array=1-1260

#Loading R app
#module load R/4.0.3-foss-2018b
module load R

#----- SECTION 1 : Preparing variables
# Please edit Number of SBATCH --array=1-11340 (number of wells x number of fields) 

# Please edit the input variable ASCTexpID with the experiment ID forexample ASCT.08.20230504

ASCTexpID="ASCT.05.20230413"
ExpN="05"
maxLCF_file="ASCT.Avg_Population_&_Tracking_maxLCF_20231001.csv"


# Please edit the input variable BasicDefault. If the BaSic default POC analysis was done.
# Options are either "Yes" or "No"
BaSicDefault="Yes"

# Please define the R scripts which will be used

Rscript6="6.2_R_Killing.Features.Extraction.R" # R tracking killing features analysis; here the new tracking algorithm has been used and no pI dynamics is being assessed. Killing Features analysis is applied across the different basic flavours

#----- SECTION 2 : Soft coded variables

path_main="/scicore/projects/rinfsci/Jovanovic/TLKK_R_analysis/ASCT_moc.poc_simplified_multipleBaSic/"$ASCTexpID
cc_file_dir="/scicore/projects/rinfsci/Jovanovic/TLKK_R_analysis/CC_csvFiles"

maxLCF_dir="/scicore/projects/rinfsci/Jovanovic/TLKK_R_analysis/ASCT_max_Live_Cell_Fraction/"$maxLCF_file

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

#outputName=$experimentName"."$experimentID"."$experimentDate"/"$experimentName"."$experimentID"."$experimentDate"."$well

# Building paths
path1=$(pwd)


#path=$path1"/"$experimentName"."$experimentID"."$experimentDate"/"$experimentName"."$experimentID"."$experimentDate"."$well
path=$path1"/"$experimentName"."$experimentID"."$experimentDate"_"$well



# List list of R script that will be run
#path_ofR6=$path1"/6_R_KillingFeatures_Pop_Trk.R"
#path_ofR6=$path1"/6_R_KillingFeautres_PopOnly_v11_scicore_20220107_FINAL.R"
#path_ofR6=$path1"/6_R_KillingFeatures_Pop_Trk_20230504_FINAL.R"
#path_ofR6=$path1"/6_R_KillingFeatures_Pop_Trk_Thrs_20230609_FINAL_NormCorr_v2_20230711_FINAL.R"
path_ofR6=$path1"/"$Rscript6


# Copy necessary files into specific well directory
cp $path_ofR6 $path/
cd $path


#----- SECTION 2.1 : Copying CC and maxLCF file into well directory

# Copying CC file to well dir

cp $cc_file_dir"/"$experimentName"."$experimentID"."$experimentDate"_CC.csv" $path/.
cp $maxLCF_dir $path/.



# Copy pooled Tracking analysis results
echo "Starting part ASCT part IV: killing features analysis"
echo $experimentName"."$experimentID"."$experimentDate"_"$well


#----- SECTION 3: Applying R analysis script in well directories

Rscript "$Rscript6"

echo "Script 6 killing features analysis is complete"


# remove CC file
rm $experimentName"."$experimentID"."$experimentDate"_CC.csv"

echo "CC file deleted from well directory"


echo "CONGRATULATIONS part IV analysis (killing features) is complete"


#----- SECTION 4: Pooling results in results sub directory


echo "Pooling killing features results"

cp $path"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_Results"/*_KF.Trkv2_longformat.csv $path_main"/PerWell-KF.Results"/.


cp $path"/"$experimentName"."$experimentID"."$experimentDate"_"$well"_Results"/*_KF.Trkv2_sanity.check_longformat.csv $path_main"/PerWell-KFsc.Results"/.


echo "Tracking time-kill curves cell numbers results is pooled"



