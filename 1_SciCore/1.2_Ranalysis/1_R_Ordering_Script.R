

#-----
# Authour: Alexander Jovanović
# Date: 2022.08.17
#-----

# Aim ----- 
# Once the image analysis is complete this script will go into the directories directory in sciore and copy all the necessary csv files into  new directory. The new directory will be in a common TLCI_Ranalysis directory. The recently copied files which represent field data (MOC, POC , artefacts and metadata) will then order in  in their respective wells. 




#2) Loading packages

library(ggpubr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(RColorBrewer)
library(platetools)
library(directlabels) 
library(MESS)
library(tictoc)



# Section 1: Deleting ilastik raw table results  -----
# Deleting files that are not necessary

# Delete files which are not used in analysis

# Delete all wiles with the _table.csv ( these are the files with all the ilastik ouput data)
delete.csv.files <- list.files(path= ".", 
                               pattern="_table.csv",
                               all.files=FALSE,
                               full.names=FALSE)


delete.csv.files

unlink( c(delete.csv.files))

#Houskeeping
rm(delete.csv.files)


# IMPORTANT: Please edit lines the two variables below accoring to experiment

# Please update the path to where the field data is found and the respective experiment ID.
# fields.data.dir <- c("/Users/Alex/Documents/PhD/PhD_project/PhD_project-Tolerance/SCICORE_Ranalysis/SCORE_orderingBOC.MOC.POC.BFin_in_parallel/Data_Dir/TLKK.83.20230114_directories_FLcorr")

# Section 2: Finding the experiment and well coordinate ID -----

experiment.ID <- list.files(path=".", 
                            pattern="OCsimple_POC_",
                            all.files=FALSE,
                            full.names=FALSE)

experiment.ID <- gsub("OCsimple_POC_","",experiment.ID)

experiment.ID <- gsub(".csv","",experiment.ID)

print("Starting ordering script")




#Create directories and move the data in their respective directores
all.csv.files <- list.files(path=".", 
                            pattern="OCsimple_POC_",
                            all.files=FALSE,
                            full.names=FALSE)



# need to find a all the wells present

finding.list.of.wells <-  list.files(path=".",
                                     pattern="OCsimple_MOC",
                                     all.files=FALSE,
                                     full.names=FALSE)

finding.list.of.wells <- gsub(
  "OCsimple_MOC_",
  "",
  finding.list.of.wells )

finding.list.of.wells <- gsub(
  ".csv",
  "",
  finding.list.of.wells )

finding.list.of.wells.df <- as.data.frame(finding.list.of.wells)

#houskeeping
rm(finding.list.of.wells)

finding.list.of.wells.df <- separate(finding.list.of.wells.df , col = finding.list.of.wells,
                                     into = c("ExpID",
                                              "Well",
                                              "Field"),
                                     sep="_")%>%
  select(-Field)%>%
  distinct()%>%
  mutate(Loop_ExpID_well = paste(ExpID,
                                 Well,
                                 sep="_"))


HOME.dir <- getwd()


#Looping here based on ExpID_well pattern


i <- finding.list.of.wells.df$Loop_ExpID_well[1]
i

# Section 3: Creating experiment and well specific sub directories and file moving -----

for (i in finding.list.of.wells.df$Loop_ExpID_well) {  
  
  directory.df <- finding.list.of.wells.df %>%
    filter(Loop_ExpID_well == i)
  
  
  print("1")
  directory.df <- directory.df%>%
    mutate(ExpID_Well = paste(ExpID,
                              "_",
                              Well,
                              sep=""))%>%
    mutate(Artefacts = paste(ExpID_Well,
                             "_Artefacts",
                             sep=""))%>%
    dplyr::mutate(MOC = paste(ExpID_Well,
                              "_OCsimple.moc",
                              sep=""))%>%
    dplyr::mutate(POC = paste(ExpID_Well,
                              "_OCsimple.poc",
                              sep=""))%>%
    dplyr::mutate(Metadata = paste(ExpID_Well,
                                   "_Metadata",
                                   sep=""))%>%
    dplyr::mutate(Boc = paste(ExpID_Well,
                              "_OCsimple.boc",
                              sep=""))%>%
    dplyr::mutate(Bfint = paste(ExpID_Well,
                                "_OCsimple.bfint",
                                sep=""))%>%
    dplyr::mutate(PerWell.res = paste(
      "PerWell-Results",
      sep=""))%>%
    dplyr::mutate(Results = paste(ExpID_Well,
                                  "_Results",
                                  sep=""))
  
  print("2")
  
  #create directory
  
  current.wd <- HOME.dir 
  
  wellDir <- directory.df$ExpID_Well
  

  real.well.Dir <- paste(current.wd,
                
                         sep="")
  
  
  print("3")
  # Create a sub directory within sub directory
  
  
  well.Artefacts.Dir  <- directory.df$Artefacts
  well.Metadata.Dir <- directory.df$Metadata
  well.MOC.Dir <- directory.df$MOC
  well.POC.Dir <- directory.df$POC
  well.RES.Dir <- directory.df$Results
  
  #-- BOC and Bfint
  well.BOC.Dir <- directory.df$Boc
  well.BFint.Dir <- directory.df$Bfint
  
  #-- PerWell results
  well.PerWellres.Dir <- directory.df$PerWell.res
  
  # Creating the respective directories
  dir.create(file.path(real.well.Dir,well.Artefacts.Dir))
  dir.create(file.path(real.well.Dir,well.Metadata.Dir ))
  dir.create(file.path(real.well.Dir,well.MOC.Dir ))
  dir.create(file.path(real.well.Dir,well.POC.Dir ))
  
  dir.create(file.path(real.well.Dir,well.BOC.Dir )) 
  dir.create(file.path(real.well.Dir,well.BFint.Dir )) 
  
  dir.create(file.path(real.well.Dir,well.RES.Dir ))
  dir.create(file.path(real.well.Dir,well.PerWellres.Dir ))
  
  
  print("4")
  #Move all files to their well directory
  
  
  well.file.pattern <- directory.df$Well
  
  well.file.pattern <- paste("_",
                             well.file.pattern,
                             "_",
                             sep="")
  
  # find all files in current directory
  well.csv.files <- list.files(path=".", 
                               pattern= well.file.pattern ,
                               all.files=FALSE,
                               full.names=FALSE)
  
  moving.function <- function(well.csv.files  = well.csv.files ) {
    file.rename( from = file.path(current.wd, well.csv.files) ,
                 to = file.path(real.well.Dir, well.csv.files) )
    
  }
  
  print("5")
  # apply the function to all files
  lapply(well.csv.files, moving.function)
  
  print("6")
  
  # In Well directory move all files to the subdirectory
  
  in.well.Dir <- paste(current.wd,
                      
                       sep="")
  setwd(in.well.Dir)
  
  
  #--- Metadata files
  # find all files in current directory
  metadata.csv.files <- list.files(path=".", 
                                   pattern= "_Metadata.csv",
                                   all.files=FALSE,
                                   full.names=FALSE)
  
  
  well.Metadata.Dir <- paste(in.well.Dir,
                             "/",
                             well.Metadata.Dir,
                             sep="")
  
  moving.function.metadata.files <- function(metadata.csv.files   = metadata.csv.files  ) {
    file.rename( from = file.path(in.well.Dir, metadata.csv.files) ,
                 to = file.path(well.Metadata.Dir , metadata.csv.files) )
    
  }
  
  
  # apply the function to all files
  lapply(metadata.csv.files , moving.function.metadata.files  )
  
  
  #--- Artefacts files
  # find all files in current directory
  qc.csv.files <- list.files(path=".", 
                             pattern= "_Artefacts.csv",
                             all.files=FALSE,
                             full.names=FALSE)
  
  
  well.Aretefacts.Dir <- paste(in.well.Dir,
                               "/",
                               well.Artefacts.Dir,
                               sep="")
  
  moving.function.qc.files <- function(qc.csv.files    = qc.csv.files   ) {
    file.rename( from = file.path(in.well.Dir,qc.csv.files ) ,
                 to = file.path(well.Aretefacts.Dir , qc.csv.files ) )
    
  }
  
  
  # apply the function to all files
  lapply(qc.csv.files  , moving.function.qc.files   )
  
  
  
  #--- POC files
  # find all files in current directory
  poc.csv.files <- list.files(path=".", 
                              pattern= "OCsimple_POC_",
                              all.files=FALSE,
                              full.names=FALSE)
  
  
  well.POC.Dir <- paste(in.well.Dir,
                        "/",
                        well.POC.Dir,
                        sep="")
  
  moving.function.poc.files <- function(poc.csv.files     = poc.csv.files   ) {
    file.rename(from = file.path(in.well.Dir,poc.csv.files  ) ,
                to = file.path(well.POC.Dir, poc.csv.files ) )
    
  }
  
  # apply the function to all files
  lapply(poc.csv.files,moving.function.poc.files )
  
  
  #--- MOC files
  # find all files in current directory
  moc.csv.files <- list.files(path=".", 
                              pattern= "OCsimple_MOC",
                              all.files=FALSE,
                              full.names=FALSE)
  
  
  well.MOC.Dir <- paste(in.well.Dir,
                        "/",
                        well.MOC.Dir,
                        sep="")
  
  moving.function.moc.files <- function(moc.csv.files     = moc.csv.files   ) {
    file.rename(from = file.path(in.well.Dir,moc.csv.files  ) ,
                to = file.path(well.MOC.Dir, moc.csv.files ) )
    
  }
  
  # apply the function to all files
  lapply(moc.csv.files,moving.function.moc.files)
  
  #--- Bf int files
  # find all files in current directory
  bfint.csv.files <- list.files(path=".", 
                                pattern= "OCsimple_BFbckint_",
                                all.files=FALSE,
                                full.names=FALSE)
  
  
  well.BFint.Dir <- paste(in.well.Dir,
                          "/",
                          well.BFint.Dir,
                          sep="")
  
  moving.function.bfint.files <- function(bfint.csv.files = bfint.csv.files   ) {
    file.rename(from = file.path(in.well.Dir,bfint.csv.files  ) ,
                to = file.path(well.BFint.Dir, bfint.csv.files ) )
    
  }
  
  
  # apply the function to all files
  lapply(bfint.csv.files,moving.function.bfint.files)
  #----------
  
  
  #--- BOC files
  # find all files in current directory
  boc.csv.files <- list.files(path=".", 
                              pattern= "OCsimple_BOC",
                              all.files=FALSE,
                              full.names=FALSE)
  
  
  well.BOC.Dir <- paste(in.well.Dir,
                        "/",
                        well.BOC.Dir,
                        sep="")
  
  moving.function.boc.files <- function(boc.csv.files = boc.csv.files   ) {
    file.rename(from = file.path(in.well.Dir,boc.csv.files  ) ,
                to = file.path(well.BOC.Dir, boc.csv.files ) )
    
  }
  
  
  
  # apply the function to all files
  lapply(boc.csv.files,moving.function.boc.files)
  
 
  setwd(HOME.dir)
  
}


# Section 4 : Check for BaSic default POC files -----

# Step 1: List all CSV files in the current directory
csv_files <- list.files(pattern = "\\.csv$")

# Step 2: Identify files with "OCsimple_POCbSd" in their filenames
matching_files <- grepl("OCsimple_POCbSd", csv_files)



if (any(matching_files)) {
  # Step 3: Create a new directory if it doesn't exist
  new_dir <- paste(wellDir,
                   "_",
                   "OCsimple.pocbSd",
                   sep="")
  
  if (!file.exists(new_dir)) {
    dir.create(new_dir)
  }
  
  # Step 4: Move the matching CSV files to the new directory
  matching_files_names <- csv_files[matching_files]
  file.rename(matching_files_names, file.path(new_dir, matching_files_names))
  
  cat("CSV files with 'OCsimple_POCbSd' in their filenames have been moved to the 'OCsimple.pocbSd' directory.\n")
} else {
  cat("No CSV files with 'OCsimple_POCbSd' in their filenames were found.\n")
}


print("CONGRATULATIONS file organisation complete ")









