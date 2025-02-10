


# Section 0  :  Defining and functions  --------------------------

#----- Functions
#Aim: Create functions for code that is repeating itself 

# -Generate function for setting variables; making directories; 
# -Loading Ocsimple; Experiment ID variables
# -Loading Metadata
# -Loading Qctrl data
# -Function saving csv files 


print("<<<<<<<<<<<<>>>>>>>>>> Running main poupulation analysis script >>>>>>><<<<<<<<>>>>>><<<")

# Loading library
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(RColorBrewer)
library(platetools)
library(directlabels) 
library(MESS)



#1) OC.df Function: Loading MOC.df 


loading.MOC.df <- function (moc.simpleDir){
  
  #Creacting vector with the list of file names
  OC.filenames <- list.files(path = moc.simpleDir,
                             pattern = "*.csv",
                             full.names = FALSE)
  
  #setting path to where all the csv files due to be processed are
  setwd(moc.simpleDir)
  tic("Vector file list complete")
  toc()
  
  tic("Loading a list of all the csv files")
  toc()
  OC.file.list <- lapply(OC.filenames,
                         read.csv)
  
  
  tic("csv file list loaded")
  toc()
  
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub(".csv","",OC.filenames)
  names(OC.file.list)<- Exp.Names
  
  #Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
  
  OC.df <- bind_rows(OC.file.list, .id = "Ilastik_Classification_ExpFile_Well_Field")
  
  OC.df <- tidyr::separate(OC.df, col = Ilastik_Classification_ExpFile_Well_Field,
                    into = c("Ilastik",
                             "Classification",
                             "ExpFile",
                             "Well_coordinate",
                             "Field"),
                    sep="_")
  
  OC.df <- OC.df %>%
    dplyr::mutate(MOC.Predicted.Class = Predicted.Class)%>%
    dplyr::select(-Ilastik,
                  -Classification,
                  -Predicted.Class)
  
  
  # Previously When running the above code on larger number of csv files (7000+) f0 ,f1, f2 columns were being produced which had NAs as values. We therefore want to remove these colunms from the OC.df
  OC.df <- OC.df %>% select_if(~!all(is.na(.))) 
  
  # Filtering out all False Cell (FC) labeled objects
  OC.df <- OC.df %>%
    dplyr::filter(MOC.Predicted.Class != "FC")%>%
    dplyr::select(-Mean.Intensity)
  
  
  OC.df$Field <- gsub("p","",OC.df$Field)
  
  OC.df$Field <- as.numeric(OC.df$Field)
  
  
  # #House keeping
  rm(OC.file.list,
     OC.filenames,
     Exp.Names
  )
  
  return(OC.df)
}






loading.POC.df <- function (poc.simpleDir){
  
  #Creacting vector with the list of file names
  OC.filenames <- list.files(path = poc.simpleDir,
                             pattern = "*.csv",
                             full.names = FALSE)
  
  #setting path to where all the csv files due to be processed are
  setwd(poc.simpleDir)
 
  OC.file.list <- lapply(OC.filenames,
                         read.csv)
  
  
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub(".csv","",OC.filenames)
  names(OC.file.list)<- Exp.Names
  
  #Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
  
  OC.df <- bind_rows(OC.file.list, .id = "Ilastik_Classification_ExpFile_Well_Field")
  
  OC.df <- tidyr::separate(OC.df, col = Ilastik_Classification_ExpFile_Well_Field,
                    into = c("Ilastik",
                             "Classification",
                             "ExpFile",
                             "Well_coordinate",
                             "Field"),
                    sep="_")
  
  OC.df <- OC.df %>%
    dplyr::mutate(POC.Predicted.Class = Predicted.Class)%>%
    dplyr::select(-Ilastik,
                  -Classification,
                  -Predicted.Class)
  
  
  OC.df <- OC.df %>% select_if(~!all(is.na(.))) 
  
  OC.df$Field <- gsub("p","",OC.df$Field)
  
  OC.df$Field <- as.numeric(OC.df$Field)
  
  
  #House keeping
  rm(OC.file.list,
     OC.filenames,
     Exp.Names
  )
  
  return(OC.df)
}



loading.POCbasic.df <- function (pocBasic.simpleDir){
  
  #Creacting vector with the list of file names
  OC.filenames <- list.files(path = pocBasic.simpleDir,
                             pattern = "*.csv",
                             full.names = FALSE)
  
  #setting path to where all the csv files due to be processed are
  setwd(pocBasic.simpleDir)
  
  OC.file.list <- lapply(OC.filenames,
                         read.csv)
  
  
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub(".csv","",OC.filenames)
  names(OC.file.list)<- Exp.Names
  
  #Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
  OC.df <- bind_rows(OC.file.list, .id = "Ilastik_Classification_ExpFile_Well_Field")
  
  OC.df <- tidyr::separate(OC.df, col = Ilastik_Classification_ExpFile_Well_Field,
                    into = c("Ilastik",
                             "Classification",
                             "ExpFile",
                             "Well_coordinate",
                             "Field"),
                    sep="_")
  
  OC.df <- OC.df %>%
    dplyr::mutate(POC.BaSic.Predicted.Class = Predicted.Class)%>%
    dplyr::mutate(Mean.IntensityBaSic = Mean.Intensity)%>%
    dplyr::select(-Ilastik,
                  -Classification,
                  -Predicted.Class,
                  -Mean.Intensity)
  
  # Previously When running the above code on larger number of csv files (7000+) f0 ,f1, f2 columns were being produced which had NAs as values. We therefore want to remove these colunms from the OC.df
  OC.df <- OC.df %>% select_if(~!all(is.na(.))) 
  OC.df$Field <- gsub("p","",OC.df$Field)
  OC.df$Field <- as.numeric(OC.df$Field)
  
  
  #House keeping
  rm(OC.file.list,
     OC.filenames,
     Exp.Names)
  
  return(OC.df)
}


#4) Loading QC dataframe
loading.QC.df <- function(qcDir, Exp.ID) {
  
  #Creacting vector with the list of file names
  QC.filenames <- list.files(path = qcDir,
                             pattern = "*.csv",
                             full.names = FALSE)
  #setting path to where all the csv files due to be processed are
  setwd(qcDir)
  tic("Vector file list complete")
  toc()
  
  tic("Loading a list of all the csv files")
  
  QC.file.list <- lapply(QC.filenames,
                         read.csv)
  tic("csv file list loaded")
  toc()
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub("_QC_Artefacts.csv","",QC.filenames)
  names(QC.file.list)<- Exp.Names
  
    #Creating QC dataframe
  QCntrl.df <- bind_rows(QC.file.list, .id = "ExpFile_Well_field")
  
  QCntrl.df <- QCntrl.df %>% mutate_all(~replace(., is.na(.),0))
  
  
  QCntrl.df <- tidyr::separate(QCntrl.df, col = ExpFile_Well_field,
                        into = c("ExpFile",
                                 "Well_coordinate",
                                 "Field"),
                        sep="_")
  
  
  QCntrl.df <- QCntrl.df %>%
  dplyr::select(-Slice)%>%
  dplyr::mutate(timestep = paste(1))
  
  QCntrl.df$timestep <- as.numeric(QCntrl.df$timestep)
  
  QCntrl.df <- QCntrl.df %>%
    dplyr::group_by(Field)%>%
    dplyr::mutate(timestep = cumsum(timestep)-1)%>%
    ungroup()%>%
    dplyr::mutate(Art.Area.percent = X.Area) %>%
    dplyr::mutate(n.Artefacts = Count)%>%
    dplyr::mutate(Mean.Int.Artefact = Mean)%>%
    select(-X.Area,
           -Count,
           -Mean)
  
  QCntrl.df$Mean.Int.Artefact <- round(QCntrl.df$Mean.Int.Artefact, digits = 2)
  QCntrl.df$Average.Size <- round(QCntrl.df$Average.Size, digits = 2)
  
  
  #Housekeeping
  rm(QC.file.list,
     QC.filenames,
     Exp.Names)
  
  
  return(QCntrl.df)
}

print("~ All functions are now loaded ~")
#------ All functions are now loaded ---

loading.metadata.df <- function (metadataDir,
                                 Exp.ID,
                                 expected.N.timestep
) {
  metaData.filenames <- list.files(path = metadataDir,
                                   pattern = "*.csv",
                                   full.names = FALSE)
  #setting path to where all the csv files due to be processed are
  setwd(metadataDir)
  
  
  metaData.file.list <- lapply(metaData.filenames,
                               read.delim2)
  
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub(".csv","",metaData.filenames)
  names(metaData.file.list ) <- Exp.Names
  
  #Creating the empty Metadata (MD) dataframe # HERE
  MD.df <- bind_rows(metaData.file.list, .id = "ExpFile_Well_Field")
  MD.df <- MD.df %>% select_if(~!all(is.na(.)))
  MD.df$ExpFile_Well_Field <- stringr::str_replace_all(MD.df$ExpFile_Well_Field, "_Metadata","")
  
  #Houskeeping
  rm(metaData.file.list,
     metaData.filenames)
  
  #creating appropriate colunm can colunm labels
  MD.df <- MD.df %>% separate(sep..,
                              sep =",",
                              c("TimeStamp", "Time [s]", "Zpos"))
  
  # In MD.df remove Zpos colunm 
  
  MD.df <- MD.df %>%
    select(-Zpos)
  
  #MD.df <- MD.df[!(MD.df$Zpos =="Zpos"),]
  
  MD.df <- MD.df %>%
    mutate(Time_sec = `Time [s]` )%>%
    mutate(timestep = TimeStamp) %>%
    select(-`Time [s]`,
           -TimeStamp)
  
  MD.df$timestep <- stringr::str_replace_all(MD.df$timestep, "TimeStamp ","")
  
  #chamging value type
  MD.df$Time_sec <- as.numeric(MD.df$Time_sec)
  MD.df$timestep <- as.numeric(MD.df$timestep)
  
  
  
  
  #MD.df$Zpos <- round(MD.df$Zpos,digits = 2)
  MD.df$Time_sec <- round(MD.df$Time_sec, digits = 2)
  
  
  # The code below will account for instances wheere when have two concantenated time lapses and we need MD of the timeframe of the first loopA to be added to the rest of the MD of LoopBs

  MD.df <- MD.df %>%
    drop_na()%>%
    group_by(ExpFile_Well_Field)%>%
    mutate(Time_frame_Order = 1,
           Time_frame_Order = cumsum(Time_frame_Order))
  

  
  MD.df<- MD.df %>%
    group_by(ExpFile_Well_Field)%>%
    mutate(new_time_sec2 = if_else(timestep != Time_frame_Order, 
                                   Time_sec + last(Time_sec[timestep == Time_frame_Order]), 
                                   Time_sec))
  
  
  MD.df <-  MD.df %>%
    mutate(Time_sec = new_time_sec2,
           timestep = Time_frame_Order)%>%
    select(ExpFile_Well_Field,
           Time_sec,
           timestep)
  
  
  

  
  MD.df <- MD.df %>%
    drop_na()%>% # removing rows with NAs
    mutate(timestep = timestep-1) #Replicating output generated by Ilastik where timestep stars from 0
  
  
  
  # Splitting the colum ExpFile_Well_Field to 3 colunms. Allows us to merge with MS.res.df
  MD.df <- tidyr::separate(MD.df, col = ExpFile_Well_Field,
                    into = c("ExpFile",
                             "Well_coordinate",
                             "Field"),
                    sep="_")
  
  MD.df <- MD.df %>%
    dplyr::select(Well_coordinate,
                  Field,
                  Time_sec,
                  timestep)
  
  # We could base the MD of a single field within a well e.g. field 18 which is in the middle of the well. However this may be a problem if by chance field 18 is missing. 
  MD.df <- MD.df %>%
    dplyr::ungroup()%>%
    dplyr::group_by(Well_coordinate,
                    timestep)%>%
    dplyr::mutate(Time_Hrs = mean(Time_sec)/3600) # MEan time across each field at a single timestep
  
  # Removing duplicate times. We expect the same number of real time (Time_Hrs) as we have timesteps
  MD.df <- MD.df %>%
    dplyr::ungroup()%>%
    select(Well_coordinate,
           timestep,
           Time_Hrs)%>%
    distinct()
  
  
  #Housekeeping
  rm(expected.N.timestep.96h.MD,
     MD.96h.df)
  
  return(MD.df)
}




#5) Writing csv files

# Define saving a df as a csv file. In this function you define the dataframe and the filename the csv will recieve. By default the csv file will always have TLCI.##.YYYYMMMDD_

tlci.writing.csv <- function(df, name.file.csv) {
  
  
  filename <- paste ("/",
                     RscriptUsed,
                     "_",
                     Exp.ID,
                     "_",
                     well.being.analysed,
                     sep ="")
  
  filename
  filename <- paste(filename,
                    "_",
                    name.file.csv,
                    sep="")
  
  filename <- paste (resDir,
                     filename, 
                     sep ="")
  
  write.csv(df,
            filename,
            row.names = FALSE)
  
  rm(filename)
  
}



# __ --------------------------

# Population analysis  --------------------------



# -----  ~Population analysis of drug tolerance  ~ 
# By Alexander Jovanovic
# PhD student
# Pulmonary Infection Biology
# Department of Biomedicine
# University Hospital Basel
# Hebelstrasse 20 
# 4031 Basel
# Date: Thursday 16th June 2022




# Aim  --------------------------

# Aim: Populational analysis of drug tolerance in Mycobacterium abscessus
# Evaluate time kill kinetics across clinical isolates treated with different antibiotics and across different antibiotic concentrations above MIC. This R script will run on a well basis.
# This R script is made up of a series of steps, generating a well-based result table. These general steps are listed below
#1) Quality Control: screening for image artefacts, gel-detachment, contamination and missing-image-frames (due to loss-of-focus thus no objects observations)
#2) Growth: evaluate growth  
#3) Object statistics: total number of single cells (SC), v-snapps (VS)
#4) Live cell fraction over time (Time-kill curve) across different definitions

#  Section 1 : loading libraries libraries   --------------------------

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


#Loading All functions


#  Section 2 : defining paths  --------------------------

#Identifying well to analyse


# Find all files in the current directory this Rscipt is found in
list.of.folders <- list.files(path=".", 
                             pattern="ASCT",
                              all.files=FALSE,
                              full.names=FALSE)

Exp.Well.to.Analyse <- as.data.frame(list.of.folders)
Exp.Well.to.Analyse <- Exp.Well.to.Analyse[!grepl("maxLCF", Exp.Well.to.Analyse$list.of.folders), ]

Exp.Well.to.Analyse <- as.data.frame(Exp.Well.to.Analyse)
Exp.Well.to.Analyse <- Exp.Well.to.Analyse %>%
  mutate(list.of.folders = Exp.Well.to.Analyse)%>%
  select(-Exp.Well.to.Analyse)

#Housekeeping
rm(list.of.folders)

#Creating a dataframe with directory details
Exp.Well.to.Analyse <- separate(Exp.Well.to.Analyse, col = list.of.folders,
                                into = c("ExpID",
                                         "Well",
                                         "Data"),
                                sep="_")


RscriptUsed <- list.files(path=".", 
                          pattern="2_R_Population.Analysis.R",
                          all.files=FALSE,
                          full.names=FALSE)

well.being.analysed <- unique(Exp.Well.to.Analyse$Well)



#1) IMPORANT: user input vairables

# ----- USER INPUT VARIABLE HERE 
#---- User Input variables
expected.N.fields <- 1:9
expected.N.timestep <- 0:29
analyse.Metadata <- c("yes") # should the Metadata data be analysed ?"yes" or "no" 
analyse.GelDetachment <- c("yes")
analyse.LC.oldThrs <- c("yes")
analyse.Contamination <-c("yes")


#1) Settign up path variables

#---

#Creating an Experiment variable e.g. TLCI.13.2021125
expFile <- unique(Exp.Well.to.Analyse$ExpID)



#--Setting working directory
genDir <- getwd() 

#General variables
genDir <- getwd() 

#---Artefacts data directory: qcDir
artefacts.df.dir <- Exp.Well.to.Analyse %>%
  filter(Data == "Artefacts")%>%
  mutate(artefacts = paste(ExpID,
                           Well,
                           Data,
                           sep="_"))

artefacts <- c(artefacts.df.dir$artefacts)
qcDir <- paste(genDir,"/",artefacts ,sep="")

#Houskeeping
rm(artefacts,
   artefacts.df.dir)
#--- Metadata directory
metadata.df.dir <- Exp.Well.to.Analyse %>%
  filter(Data == "Metadata")%>%
  mutate(metadata = paste(ExpID,
                          Well,
                          Data,
                          sep="_"))

metadata <- c(metadata.df.dir$metadata)
metadataDir <- paste(genDir,"/", metadata, sep="")

#Houskeeping
rm(metadata,
   metadata.df.dir)

#--- Results directory
res.df.dir <- Exp.Well.to.Analyse %>%
  filter(Data == "Results")%>%
  mutate(results = paste(ExpID,
                         Well,
                         Data,
                         sep="_"))

results<- c(res.df.dir$results)
resDir <- paste(genDir,"/",results, sep="")

#Houskeeping
rm(results,
   res.df.dir)

#----- Morphology Object classification (MOC) directory
moc.df.dir <- Exp.Well.to.Analyse %>%
  filter(Data == "OCsimple.moc")%>%
  mutate(OCsimple.moc = paste(ExpID,
                              Well,
                              Data,
                              sep="_"))

moc <- c(moc.df.dir$OCsimple.moc)

moc.simpleDir <- paste(genDir,
                       "/",
                       moc,
                       sep="")

#Houskeeping
rm(moc.df.dir,
   moc)


#----- Propidium-iodide Object classification (POC) directory
poc.df.dir <- Exp.Well.to.Analyse %>%
  filter(Data == "OCsimple.poc")%>%
  mutate(OCsimple.poc = paste(ExpID,
                              Well,
                              Data,
                              sep="_"))

poc <- c(poc.df.dir$OCsimple.poc)

poc.simpleDir <- paste(genDir,
                       "/",
                       poc,
                       sep="")
#Houskeeping
rm(poc.df.dir,
   poc)

#-----  BaSic default settings Propidium-iodide Object classification (POC) directory

pocBasic.simpleDir <- gsub("OCsimple.poc",
                           "OCsimple.pocbSd",
                           poc.simpleDir)


#Houskeeping
rm(pocBaSic.df.dir,
   pocBasic,
   Exp.Well.to.Analyse)

#Setting work directory
wdDir <- genDir
setwd(wdDir)



# Section 3  : Loading MOC and POC data  --------------------------

#2) Loading and preparing OCsimple.df 

#Using custom made ASCT_pop_perWell_mocpoc_v0_functions.Rmd we will load the different data sets 
#---1---
MOC.df <- loading.MOC.df(moc.simpleDir = moc.simpleDir)

POC.df <- loading.POC.df(poc.simpleDir = poc.simpleDir)

POC.BaSicdef <-  loading.POCbasic.df(pocBasic.simpleDir= pocBasic.simpleDir)

# POC.df  <- left_join(POC.df,
#                                      POC.BaSicdef,
#                                      by = join_by(ExpFile, Well_coordinate, Field, timestep, labelimage_oid, Center.of.the.object_0, Center.of.the.object_1))
POC.df <- dplyr::left_join(POC.df,
                                 POC.BaSicdef,
                                 by= c("ExpFile", "Well_coordinate", "Field", "timestep", "labelimage_oid", "Center.of.the.object_0", "Center.of.the.object_1"))
rm(POC.BaSicdef)

Exp.ID <- expFile

# Section 4  :  Zero bytes  --------------------------

#3) M/POC-Fields missing: Zero bytes

# When analysis times out on scicore there are instances where a M/POC-OCsimple file is written but has no lines (Zero bytes). Before creating either MOC or POC directors we must before start of the  R analysis delete all Zero bytes files otherwise the data cannot be loaded. 
# Knowing this we need t determine whether there are Zero byte files in the Well and if yes which files. Otherwise we say NA
# aim: Detemrine whether there are missing MOC files

# Find fields missing from the expected fields value
expected.Fields <- expected.N.fields

# Fields present in analysis
actual.moc.Fields <- unique(MOC.df$Field)

#Converting vector in to single string vector of 1 element
missing.moc.Fields <- expected.Fields[!(expected.Fields %in% actual.moc.Fields)]

missing.moc.Fields <- as.character(missing.moc.Fields)
missing.moc.Fields <- paste(missing.moc.Fields, collapse = "_")
missing.moc.Fields

# Creating a  Master results dataframe which will be the main results dataframe
MS.res.df <- MOC.df %>%
  select(ExpFile,
         Well_coordinate)%>%
  distinct()

#Adding a new colunm varibale to results dataframe
MS.res.df <- cbind(MS.res.df,
                   missing.moc.Fields)

MS.res.df <- MS.res.df %>%
  dplyr::mutate(MOC.zeroBytes.fields = missing.moc.Fields)%>%
  dplyr::select(-missing.moc.Fields)

#Houskeeping
rm(actual.moc.Fields,
   expected.Fields,
   missing.moc.Fields)

# If there were no missing fields convert all empty cell into NA
MS.res.df <- MS.res.df %>%
  mutate_all(na_if,"")


#----------------------------------- POC 
# Find fields missing from the expected fields value
expected.Fields <- expected.N.fields

# Fields present in analysis
actual.poc.Fields <- unique(POC.df$Field)

#Converting vector in to single string vector of 1 element
missing.poc.Fields <- expected.Fields[!(expected.Fields %in% actual.poc.Fields)]

missing.poc.Fields <- as.character(missing.poc.Fields)
missing.poc.Fields <- paste(missing.poc.Fields, collapse = "_")
missing.poc.Fields

# Creating a  Master results dataframe which will be the main results dataframe


#Adding a new colunm varibale to results dataframe
MS.res.df <- cbind(MS.res.df,
                   missing.poc.Fields)

MS.res.df <- MS.res.df %>%
  dplyr::mutate(POC.zeroBytes.fields = missing.poc.Fields)%>%
  dplyr::select(-missing.poc.Fields)

#Houskeeping
rm(actual.poc.Fields,
   expected.Fields,
   missing.poc.Fields)

# If there were no missing fields convert all empty cell into NA
MS.res.df <- MS.res.df %>%
  dplyr::mutate_all(na_if,"")



#4) MS.res.df: Expected timestep colunm

# Expected number of image frame in experiment. This is defined by the user before the start of the analysis

MS.res.df <- cbind(MS.res.df,
                   expected.N.timestep)


MS.res.df <- MS.res.df %>%
  mutate(timestep = expected.N.timestep)%>%
  select(-expected.N.timestep)




# Section 5   : Metadata  --------------------------

#4.1) Loading Metadata and appeding to Mater.result dataframe

if ( analyse.Metadata == "yes") {
  
  MD.df <- loading.metadata.df(metadataDir = metadataDir,
                               Exp.ID = Exp.ID,
                               expected.N.timestep =expected.N.timestep)
  
  MS.res.df <- left_join(MS.res.df,
                         MD.df)
  
  #Houskeeping 
  rm(MD.df)
} else {
  
  Time_Hrs <- c(NA)
  
  MS.res.df <- cbind( MS.res.df , 
                      Time_Hrs)
  
  rm(Time_Hrs)
  # Metadata: Time_hrs colunm NA to shouw it was not assessed
  (print ("Metadata integration step skipped"))}






#5) Merging MOC and POC dataframes

#Merging Ilastik MOC and POC data using common variables by = c("ExpFile", "Well_coordinate", "Field", "timestep", "labelimage_oid", "Center.of.the.object_0", "Center.of.the.object_1")
moc.poc.df <- dplyr::left_join(MOC.df,
                               POC.df)

#Huskeeping
rm(MOC.df,
   POC.df)


# Section 6   : Fluorescent artefact objects  --------------------------

#6) QC: Load quality control data (QCntrl) and writing csv

# Certain fields may have large flourescent artefacts (e.g. undissolved-gel) . We screen all fields and determine whether large fluorescent artefacts are present in the first and second image frame (timestep 0 and 1). We assume so early on it is unlikely to be growth and should reflect and artefact.

# Load all the Artefacts.csv files into a Quality Analysis dataframe
QA.df <- loading.QC.df(qcDir = qcDir,
                       Exp.ID = Exp.ID)

#Creating a Well_field coorindate variable
QA.df$Field <- sub(".*p", "", QA.df$Field) # Extracxting element after delimiter
QA.df$Field <- as.numeric(QA.df$Field)


QA.df <- QA.df %>%
  dplyr::mutate(Field_coordinate = paste(Well_coordinate,
                                         "_p",
                                         Field,
                                         sep=""
  ))

QA.df <- QA.df %>%
  select(Well_coordinate,
         Field_coordinate,
         Field,
         Total.Area,
         Average.Size,
         timestep,
         Art.Area.percent,
         n.Artefacts,
         Mean.Int.Artefact)


#6.1) QC: assessing field quality

# We are identifying fields with large fluorescence artefacts by investigating the Artifact area. 
# 6.1.1) Isolating  first two image frames of each fields with at least 1.5 % of artefact area
QA.ArtefactFields <- QA.df %>%
  ungroup()%>%
  dplyr::filter(timestep <= 1)%>%
  dplyr::filter(Art.Area.percent >= 1.5) # Feel this covers the large which are clearly visible in the first 2 frames.


# NOTIC: NEED TO CODE IF 20% of fields within well shows artefacts omit  entire WELL 

#If at least 20% of total fields in the well have an artifact quality of well is declared poor. Otherwise its declared good.
#Finding number of fields with poor quality
n.of.Artefact.fields <- QA.ArtefactFields %>%
  dplyr::select(Field)%>%
  dplyr::distinct()


n.of.Artefact.fields <- nrow(n.of.Artefact.fields)

good.or.poor.quality <- (n.of.Artefact.fields /max(expected.N.fields))*100

QC_Well.image.quality <- ifelse (good.or.poor.quality >= 20,"poor", "good" )

MS.res.df <- cbind(MS.res.df,
                   QC_Well.image.quality)

#Houskeeping
rm(n.of.Artefact.fields,
   good.or.poor.quality,
   QC_Well.image.quality)

#6.1.2) Identify fields with artefacts if there are none Label as NA 
QC_artefact.fields <- QA.ArtefactFields %>%
  dplyr::select(Field)%>%
  dplyr::distinct()

QC_artefact.fields <- QC_artefact.fields$Field

QC_artefact.fields <-  paste(  QC_artefact.fields, collapse = "_")

MS.res.df <- cbind(MS.res.df,
                   QC_artefact.fields) 

#Houskeeping
rm(QC_artefact.fields)

#ERROR HERE 
# If there were no artefact fields convert all empty cell into NA
#MS.res.df <- MS.res.df %>%
 # dplyr::mutate_all(na_if,"")
 #dplyr::mutate_all(~ifelse(. == "", "NA", .))


MS.res.df[, 8][MS.res.df[, 8] == ""] <- NA


#6.1.3) If there are fields that were reject how many out of total number of fields were rejected?
n.fields.rejected  <- QA.ArtefactFields %>%
  dplyr::select(Field)%>%
  dplyr::distinct()

n.fields.rejected <- nrow(n.fields.rejected )


n.fields.rejected <- paste(n.fields.rejected,
                           "/",
                           max(expected.N.fields),
                           sep="")

MS.res.df <- cbind(MS.res.df,
                   n.fields.rejected)

MS.res.df <- MS.res.df %>%
  dplyr::mutate(QC_n.Artefact.fields.rejected = n.fields.rejected)%>%
  dplyr::select(-n.fields.rejected)

#Houskeeping
rm(n.fields.rejected)

# Filter out fields with artefacts from the main moc.poc dataframe
#df with a single colunm showing field coordinates identified as fields with artefacts
QA.ArtefactFields <- QA.ArtefactFields %>%
  dplyr::ungroup() %>%
  dplyr::select(Field_coordinate)%>%
  dplyr::distinct()

#Vector with fields to omit
QA_Fields_to_Omitt <- QA.ArtefactFields$Field_coordinate

# Need a function that will filter  moc.poc.df if QA_Fields_to_Omitt >0 i.e has at least one field to filter by or not filter at all
# call it moc.poc.qc.df
for (length in length(QA_Fields_to_Omitt)) {
  if  (length == 0) {
    moc.poc.qc.df <- moc.poc.df}
  
  else 
    moc.poc.qc.df<- moc.poc.df %>%
      mutate(Field_coordinate = paste(Well_coordinate,
                                      "_p",
                                      Field,
                                      sep=""))%>%
      ungroup()%>%
      dplyr::filter(Field_coordinate != QA_Fields_to_Omitt)%>%
      select(-Field_coordinate)
}

#Delete the unfilted  (wrt to artefacts) moc.poc.df
#Houskeeping
rm(moc.poc.df)


#Houskeeping
rm(QA.ArtefactFields,
   QA.df)

# Section 6.1  : Object numbers  --------------------------

# 6.1) Assessing minimum number of objects at 0th frame

# Well with too few cells i.e less than 500 should be flagged for having too few cells. Thats too for statistics. 
# Aim: Evaluate whether well has less than or more than 500 cells. If less than evaluate number of cells if more than state how many objects 

# Stating whether we have at least 500 SC+ VS in frame 0 
moc.N.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  filter(timestep == 0)%>%
  filter(MOC.Predicted.Class == "SC" |
           MOC.Predicted.Class == "VS" )%>%
  mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  mutate(Sum.Objects = sum(Objects))%>%
  ungroup()%>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         Sum.Objects)%>%
  distinct()%>%
  mutate(QC_minimum.N.cells.tp0 = dplyr::if_else( Sum.Objects > 500, "sufficient.N.cells", "insufficient.N.cells"))%>%
  select(-Sum.Objects,
         -timestep)


MS.res.df <- left_join(MS.res.df,
                       moc.N.df)

#Houskeeping
rm(moc.N.df)



#--- Total number (N) of OFFC

moc.N.offc.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "OFFC")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(QC_n.OFFC = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                QC_n.OFFC)%>%
  dplyr::distinct()
  
MS.res.df <- left_join(MS.res.df,
                       moc.N.offc.df)

#Houskeeping
rm(moc.N.sc.df)

#--- Total number (N) of SC
moc.N.sc.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "SC")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(QC_n.SC = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                QC_n.SC)%>%
  dplyr::distinct()

MS.res.df <- left_join(MS.res.df,
                       moc.N.sc.df)

#Houskeeping
rm(moc.N.sc.df)

#--- Total number (N) of VS
moc.N.vs.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "VS")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(QC_n.VS = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                QC_n.VS)%>%
  dplyr::distinct()

MS.res.df <- left_join(MS.res.df,
                       moc.N.vs.df)

#Hosukeeping
rm(moc.N.vs.df)


#--- Total number (N) of C2
moc.N.c2.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "C2")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(QC_n.C2 = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                QC_n.C2 )%>%
  dplyr::distinct()

MS.res.df <- left_join(MS.res.df,
                       moc.N.c2.df)

#Hosukeeping
rm(moc.N.c2.df )


#--- Total number (N) of C5
moc.N.c5.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "C5")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(QC_n.C5 = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                QC_n.C5 )%>%
  dplyr::distinct()

MS.res.df <- left_join(MS.res.df,
                       moc.N.c5.df )

#Hosukeeping
rm(moc.N.c5.df )


#--- Total number (N) of C20
moc.N.c20.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "C20")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(QC_n.C20 = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                QC_n.C20 )%>%
  dplyr::distinct()

MS.res.df <- left_join(MS.res.df,
                       moc.N.c20.df  )

#Hosukeeping
rm(moc.N.c20.df )

#--- Total number (N) of SC + VS
moc.N.sc.vs.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class == "VS")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(QC_n.SC.VS.cells = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                QC_n.SC.VS.cells)%>%
  dplyr::distinct()

MS.res.df <- left_join(MS.res.df,
                       moc.N.sc.vs.df)

#Hosukeeping
rm(moc.N.sc.vs.df)


#--- Total number (N) of C2 C5 C20
moc.N.c2c5c20.df <- moc.poc.qc.df %>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "C2" | MOC.Predicted.Class == "C5" |  MOC.Predicted.Class == "C20")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(QC_n.c2.c5.c20 = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                QC_n.c2.c5.c20)%>%
  dplyr::distinct()

MS.res.df <- left_join(MS.res.df,
                       moc.N.c2c5c20.df )

#Hosukeeping
rm(moc.N.c2c5c20.df,
   moc.N.offc.df)

# Section 6.1.1 : Total object area --------------------------

# Generate Total.area  SC, VS ,C2, C5 and C20
# group by different classes and calulate total area >> if it works apply this calulation above

# Generate Total.area  SC   
tot.sc.Obj.Area.df <- moc.poc.qc.df %>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                Size.in.pixels,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "SC")%>%
  dplyr::group_by(MOC.Predicted.Class,
                  timestep)%>%
  mutate(QC_TotArea.sqrPx.sc = sum(Size.in.pixels))%>%
  select(-Size.in.pixels,
         -MOC.Predicted.Class)%>%
  distinct()


MS.res.df <- left_join(MS.res.df,
                       tot.sc.Obj.Area.df )

#Hosukeeping
rm(tot.sc.Obj.Area.df)

# Generate Total.area  VS   
tot.vs.Obj.Area.df <- moc.poc.qc.df %>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                Size.in.pixels,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "VS")%>%
  dplyr::group_by(MOC.Predicted.Class,
                  timestep)%>%
  dplyr::mutate(QC_TotArea.sqrPx.vs = sum(Size.in.pixels))%>%
  dplyr::ungroup()%>%
  dplyr::select(-Size.in.pixels,
                -MOC.Predicted.Class)%>%
  dplyr::distinct()


MS.res.df <- left_join(MS.res.df,
                       tot.vs.Obj.Area.df )

#Hosukeeping
rm(tot.vs.Obj.Area.df)



# Generate Total.area  OFFC   
tot.offc.Obj.Area.df <- moc.poc.qc.df %>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                Size.in.pixels,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "OFFC")%>%
  dplyr::group_by(MOC.Predicted.Class,
                  timestep)%>%
  dplyr::mutate(QC_TotArea.sqrPx.offc = sum(Size.in.pixels))%>%
  dplyr::ungroup()%>%
  dplyr::select(-Size.in.pixels,
                -MOC.Predicted.Class)%>%
  dplyr::distinct()


MS.res.df <- left_join(MS.res.df,
                       tot.offc.Obj.Area.df )

#Hosukeeping
rm(tot.offc.Obj.Area.df)

# Generate Total.area  C2   
tot.c2.Obj.Area.df <- moc.poc.qc.df %>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                Size.in.pixels,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "C2")%>%
  dplyr::group_by(MOC.Predicted.Class,
                  timestep)%>%
  dplyr::mutate(QC_TotArea.sqrPx.c2 = sum(Size.in.pixels))%>%
  dplyr::ungroup()%>%
  dplyr::select(-Size.in.pixels,
                -MOC.Predicted.Class)%>%
  dplyr::distinct()


MS.res.df <- left_join(MS.res.df,
                       tot.c2.Obj.Area.df )

#Hosukeeping
rm(tot.c2.Obj.Area.df)


# Generate Total.area  C5  
tot.c5.Obj.Area.df <- moc.poc.qc.df %>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                Size.in.pixels,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "C5")%>%
  dplyr::group_by(MOC.Predicted.Class,
                  timestep)%>%
  dplyr::mutate(QC_TotArea.sqrPx.c5 = sum(Size.in.pixels))%>%
  dplyr::ungroup()%>%
  dplyr::select(-Size.in.pixels,
                -MOC.Predicted.Class)%>%
  dplyr::distinct()


MS.res.df <- left_join(MS.res.df,
                       tot.c5.Obj.Area.df )

#Hosukeeping
rm(tot.c5.Obj.Area.df)


# Generate Total.area  C5  
tot.c20.Obj.Area.df <- moc.poc.qc.df %>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                Size.in.pixels,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "C20")%>%
  dplyr::group_by(MOC.Predicted.Class,
                  timestep)%>%
  dplyr::mutate(QC_TotArea.sqrPx.c20 = sum(Size.in.pixels))%>%
  dplyr::ungroup()%>%
  dplyr::select(-Size.in.pixels,
                -MOC.Predicted.Class)%>%
  dplyr::distinct()


MS.res.df <- left_join(MS.res.df,
                       tot.c20.Obj.Area.df )

#Hosukeeping
rm(tot.c20.Obj.Area.df)


# Section 7 : Growth  --------------------------


# 7) Growth Evaluation: Total area / Total object numbers 

# Certain isolates despite treatment are able to grow. Such isolates we assume are resistant. During growth single cells and vsnapping (dividing cells) transition in terms of their morphology to form clusters.
# Aim: Evaluate whether  growth is detected. IF growth is detected find timestep threshold is surpassed and caulatate  timepoint before this instance.

# Calculating the total object  (TO) area per timestep (frame)
moc.TO.area.df <- moc.poc.qc.df %>%
  ungroup()%>%
  dplyr::select(Well_coordinate,
                timestep,
                Size.in.pixels)%>%
  dplyr:: group_by(
    Well_coordinate,
    timestep)%>%
  dplyr::mutate(TotObj.area = sum(Size.in.pixels))%>%
  dplyr::select(-Size.in.pixels)%>%
  distinct()

# Total number of Objects per image frame
moc.TnO.df <- moc.poc.qc.df%>%
  ungroup()%>%
  dplyr::select(Well_coordinate,
                timestep,
                Size.in.pixels)%>%
  group_by(timestep,
           Well_coordinate)%>%
  summarise(TotObj.n = n()) 

# total Object Area Number. df 
moc.tot.OAN.df <- left_join(moc.TO.area.df ,
                            moc.TnO.df)

#Houkeeping
rm(moc.TnO.df,
   moc.TO.area.df)

# Calulating the corrected object area over time 


#Calulating corrected Area
moc.tot.OAN.df <- moc.tot.OAN.df %>%
  ungroup()%>%
  group_by(Well_coordinate)%>%
  arrange(timestep)%>%
  dplyr::mutate(Corr.N.tot.ObjArea = (TotObj.area/TotObj.n)/(dplyr::first(TotObj.area)/dplyr::first(TotObj.n)) )

# Append Corrected object Area to MS.res.df file
moc.tot.OAN.df <- moc.tot.OAN.df %>%
  dplyr::ungroup()%>%
  dplyr::mutate(GT_TotObj.area = TotObj.area ,
                GT_TotObj.n = TotObj.n,
                GT_CorrObjArea = Corr.N.tot.ObjArea)%>%
  select(Well_coordinate,
         timestep,
         GT_TotObj.area ,
         GT_TotObj.n ,
         GT_CorrObjArea)



# Threshold for growth (GT) to be detected is set at 6.5 
GT <- 6.5

#Evaluate growth threshold using a binary score. If growth threshold (GT) is surpassed a score of 1 is given if below 0.
# No isntance where threshold is surpassed
#Single instance where threshold is surpassed
moc.tot.OAN.df  <- moc.tot.OAN.df %>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate)%>%
  dplyr::mutate(GT_Growth.thrs.6.5 = dplyr::if_else(GT_CorrObjArea >= GT,
                                                    1,
                                                    0))%>% # Assessing whether threshold has been surpassed 
  dplyr::mutate(GT_Growth.eval.in.well = dplyr::if_else( sum(GT_Growth.thrs.6.5) >0, "growth-detected","no-growth-detected" ))%>% # if the we have a
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep)%>%
  dplyr::mutate(GT_min.tp.GT = dplyr::if_else(GT_Growth.thrs.6.5 ==1 , as.numeric(paste(timestep-1)),as.numeric("NA")))%>%
  dplyr::ungroup()



GT_min.tp.GT <-  moc.tot.OAN.df  %>%
  dplyr::ungroup()%>%
  dplyr::select(GT_min.tp.GT)%>%
  drop_na()%>%
  mutate(GT_min.tp.GT = paste(min(GT_min.tp.GT)))%>%
  dplyr::distinct()

GT_eval.min.tp.GT <- GT_min.tp.GT$GT_min.tp.GT

# Some wells may not be flagged for growth therefore will not have a minimum timepoint for growth. In that case we check the length of  the  GT_min.tp.GT vector. If it is 0 we know there was no growth and thus concatenate (NA). If the length is > 0 we know we have a value and thus shall do nothing. This is mimporant because we need to cbind again the moc.tot.OAN.df. If no growth was found we would have length vector 0 and thus cannot cbind against a dataframe with 13 rows.
if ( length( GT_eval.min.tp.GT) == 0) {
  
  GT_eval.min.tp.GT <- c(NA) } else { GT_eval.min.tp.GT <-  GT_eval.min.tp.GT } 

# append acoliunm
moc.tot.OAN.df <- cbind(moc.tot.OAN.df , GT_eval.min.tp.GT)


#Houskeeping
rm(GT_min.tp.GT)

moc.tot.OAN.df <- moc.tot.OAN.df %>%
  select(-GT_min.tp.GT)

# Append to MS.res.df   
MS.res.df <- left_join(MS.res.df,
                       moc.tot.OAN.df)






# Section 8 : Contamination --------------------------

# 8) Contamination


# To idnetify wells which are potentially contaminated we will calculate the ratio of SC numbers over time. In the past we observed some contamination with Staph.aureus, which has a similar  rod shape to M.abscessus. As a result the machine learning algorithm defines foreign bacteria  as a single cell during classification. So far we have not observed contamination during timestep 0 but instead single cells appearing over time, thus we will measure the increase in ratio of single cells over time.

# Calulate the change in ratio of total single cell wrt to first total single cell number. If the number of SC increase by to ratio > 2 we assume contamination
moc.Contamination.df <- moc.poc.qc.df %>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                MOC.Predicted.Class)%>%
  dplyr::filter(MOC.Predicted.Class == "SC" |   MOC.Predicted.Class == "OFFC")%>%
  dplyr::mutate(Objects = 1)%>%
  dplyr::group_by( Well_coordinate,
                   timestep)%>%
  dplyr:: mutate(Contamination_Sum.sc = sum(Objects))%>%
  dplyr::ungroup()%>%
  dplyr::select(ExpFile,
                Well_coordinate,
                timestep,
                Contamination_Sum.sc)%>%
  dplyr::distinct()%>%
  dplyr::ungroup()%>%
  mutate(Contamination_ratio = Contamination_Sum.sc/dplyr::first(Contamination_Sum.sc))%>%
  mutate(Contamination_ratio.thrs = dplyr::if_else(Contamination_ratio >2, 1,0))%>%
  mutate(Contamination_ratio.thrs  = cumsum(Contamination_ratio.thrs))%>%
  mutate(Contamination_ratio.thrs  = dplyr::if_else(Contamination_ratio.thrs >=1 , 1, 0 ))%>% # We now found  at what timepoint the contamination start appearing (1) but we want to keep all the informattion before contamination crops up. So we now calculate the timestep before contaminaotion is detected  
  mutate(Contamination_eval = dplyr::if_else(sum(Contamination_ratio.thrs) >=1 , "Contamination", "No-contamination"))%>%
  mutate(Contamination_timestep = dplyr::if_else(Contamination_ratio.thrs ==1, paste(timestep-1), "NA"))

cnt.timestep  <- moc.Contamination.df$Contamination_timestep

cnt.timestep <- as.numeric(cnt.timestep )

cnt.timestep <- na.omit(cnt.timestep)

if (length(cnt.timestep > 0 ) ){
  
  cnt.timestep <- min(cnt.timestep)
} else {cnt.timestep <- NA}

moc.Contamination.df <- moc.Contamination.df %>%
  mutate(Contamination_timestep = cnt.timestep )



#Houskeeping 
rm(cnt.timestep )

MS.res.df <- left_join(MS.res.df,
                       moc.Contamination.df)


MS.res.df$Contamination_timestep <- as.numeric(MS.res.df$Contamination_timestep)


#Houskeeping
rm(moc.Contamination.df)





# Section 9 : Assessing Gel detachment --------------------------

#9) Assessing Gel detachment


if ( analyse.GelDetachment == "yes") {
  #---- STICK geld-det ANALYSIS CODE HERE 
  
  
  # In the event of gel detachment from the plate bottom,  there is a rapid decrease in total object number and  total object area from one frame to another.
  
  # Aim: 1) Evaluate whether gel-detachment is occurring in the well. 2) If gel detachment is evident find the instance of detachment.
  moc.GD.df <- moc.tot.OAN.df %>%
    dplyr::ungroup()%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  GT_TotObj.area,
                  GT_TotObj.n)%>%
    mutate(GD_Tot.Area = GT_TotObj.area,
           GD_Tot.n = GT_TotObj.n)%>%
    select(-GT_TotObj.area,
           - GT_TotObj.n)
  
  #Houskeeping
  rm(moc.tot.OAN.df)
  
  moc.GD.df <- moc.GD.df %>%
    dplyr::ungroup()%>%
    dplyr::group_by(Well_coordinate)%>%
    arrange(timestep)%>%
    mutate(GD_Ratio.A = GD_Tot.Area / first(GD_Tot.Area ))%>%
    mutate(GD_Ratio.N = GD_Tot.n  / first(GD_Tot.n ))
  
  #--- Ratio N 
  #rolling window of 3 frames from 0-48h then roll window of single frame + rbind and find lowest
  
  #  "N" 48h step I (sI): calculate rolling window for 0-48 by 3 frames
  moc.GD.tl.df.rollN <- moc.GD.df %>% #time lapse (tl) 
    ungroup()%>%
    dplyr::filter(timestep < (max(timestep)))%>% # assessing everything before the last image frame which was the 96th hour image. time between 11th and 12th frame is equivalent to 48h-96h img frame
    group_by(Well_coordinate)%>%
    arrange(timestep)%>%
    mutate(GD_Roll.Ratio.N.tl = lead(GD_Tot.n, n=3)/(GD_Tot.n))%>%  
    ungroup()%>%
    mutate(GD_Roll.Ratio.N.tl = if_else(is.na(GD_Roll.Ratio.N.tl), 0,GD_Roll.Ratio.N.tl))%>%
    ungroup()
  
  moc.GD.tl.df.rollA <- moc.GD.df %>% #time lapse (tl) 
    ungroup()%>%
    dplyr::filter(timestep < (max(timestep)))%>%
    group_by(Well_coordinate)%>%
    arrange(timestep)%>%
    mutate(GD_Roll.Ratio.A.tl = lead(GD_Tot.Area, n=3)/(GD_Tot.Area))%>%  
    ungroup()%>%
    mutate(GD_Roll.Ratio.A.tl = if_else(is.na(GD_Roll.Ratio.A.tl), 0,GD_Roll.Ratio.A.tl))%>%
    ungroup()%>%
    select(Well_coordinate,
           timestep,
           GD_Roll.Ratio.A.tl)
  
  moc.GD.tl.df.rollNA.48h <- left_join(moc.GD.tl.df.rollN,
                                       moc.GD.tl.df.rollA)
  
  
  
  #Houskeeping
  rm(moc.GD.tl.df.rollN,moc.GD.tl.df.rollA)
  
  moc.GD.tl.df.rollNA.48h <- moc.GD.tl.df.rollNA.48h  %>%
    filter(GD_Roll.Ratio.N.tl != 0)
  
  
  # ---- single image 48-96h  
  # N 48-96h step II (sI): caluclate rolling window for 48-96 by 1 frames
  moc.GD.tl.df.rollN.96h <- moc.GD.df %>%
    ungroup()%>%
    dplyr::filter(timestep >= (max(timestep))-1)%>%
    group_by(Well_coordinate)%>%
    arrange(timestep)%>%
    mutate(GD_Roll.Ratio.N.96h = lead(GD_Tot.n, n=1)/(GD_Tot.n))%>%  # calling it Roll ratio Roll.Ratio.A.1 image frame
    ungroup()%>%
    mutate(GD_Roll.Ratio.N.96h = if_else(is.na(GD_Roll.Ratio.N.96h), 0,GD_Roll.Ratio.N.96h))
  
  # A 48-96h step II (sI): caluclate rolling window for 48-96 by 1 frames
  moc.GD.tl.df.rollA.96h <-  moc.GD.df %>%
    ungroup()%>%
    dplyr::filter(timestep >= (max(timestep))-1)%>%
    group_by(Well_coordinate)%>%
    arrange(timestep)%>%
    mutate(GD_Roll.Ratio.A.96h = lead(GD_Tot.Area, n=1)/(GD_Tot.Area))%>%  # calling it Roll ratio Roll.Ratio.A.1 image frame
    ungroup()%>%
    mutate(GD_Roll.Ratio.A.96h = if_else(is.na(GD_Roll.Ratio.A.96h), 0,GD_Roll.Ratio.A.96h))
  
  moc.GD.tl.df.rollNA.96h <- left_join(moc.GD.tl.df.rollN.96h,
                                       moc.GD.tl.df.rollA.96h)
  
  
  #Houskeeping
  rm(moc.GD.tl.df.rollN.96h,
     moc.GD.tl.df.rollA.96h)
  
  #-- finding the minimum Rolling ration in the moc.GD.tl.df.rollNA.48h dataframe
  moc.GD.tl.df.rollNA.48h <- moc.GD.tl.df.rollNA.48h %>%
    dplyr::ungroup()%>%
    mutate(GD_Roll.Ratio.N.48h.max.decline = min(GD_Roll.Ratio.N.tl),
           GD_Roll.Ratio.A.48h.max.decline = min (GD_Roll.Ratio.A.tl))%>%
    select(Well_coordinate,
           GD_Roll.Ratio.N.48h.max.decline,
           GD_Roll.Ratio.A.48h.max.decline)%>%
    distinct()
  
  
  moc.GD.tl.df.rollNA.96h <- moc.GD.tl.df.rollNA.96h %>%
    dplyr::ungroup()%>%
    filter(GD_Roll.Ratio.N.96h != 0)%>%
    mutate(GD_Roll.Ratio.N.96h.max.decline = min(GD_Roll.Ratio.N.96h),
           GD_Roll.Ratio.A.96h.max.decline = min (GD_Roll.Ratio.A.96h))%>%
    select(Well_coordinate,
           GD_Roll.Ratio.N.96h.max.decline,
           GD_Roll.Ratio.A.96h.max.decline)%>%
    distinct()
  
  moc.GD.tl.df.rollNA.48.96h <- left_join(moc.GD.tl.df.rollNA.48h,
                                          moc.GD.tl.df.rollNA.96h)
  
  #Houskeeping
  rm(moc.GD.tl.df.rollNA.48h,
     moc.GD.tl.df.rollNA.96h)
  
  
  moc.GD.tl.df.rollNA.48.96h <- moc.GD.tl.df.rollNA.48.96h %>%
    ungroup()%>%
    # drop_na()%>%
    group_by(Well_coordinate)%>%
    mutate(GD_Eval.0.48h = if_else( GD_Roll.Ratio.A.48h.max.decline <= 0.8 & GD_Roll.Ratio.N.48h.max.decline <= 0.9, 1,0))%>%
    mutate(GD_Eval.48.96h = if_else (GD_Roll.Ratio.A.96h.max.decline <= 0.6 &  GD_Roll.Ratio.N.96h.max.decline <= 0.6, 1,0)) %>%
    mutate(GD_eval.gelDet= GD_Eval.0.48h + GD_Eval.48.96h )%>%
    mutate(GD_eval.gelDet = if_else(GD_eval.gelDet >0 , "Det", "No-det"))
  
  
  moc.GD.df <- left_join(moc.GD.df,
                         moc.GD.tl.df.rollNA.48.96h)
  
  #houskeeping
  rm(moc.GD.tl.df.rollNA.48.96h)


  x <-  moc.GD.df$GD_eval.gelDet[1]
  if ( x == "Det" && moc.GD.df$GD_Eval.0.48h == 1 ) {
    #------ Instance of Gel detachment if Gel detachment is detcted  
    
    print("1")
    # Determine whether instance of detachment within 0-48h. We define gel detachment to be a decrease in total object area by at least 5% from one frame to another
    # if (moc.GD.df$GD_Eval.0.48h == 1) {
    
    moc.GD.df.Eval.inst.0.48h.df <- moc.GD.df %>%
      dplyr::ungroup()%>%
      select(Well_coordinate,
             timestep,
             GD_Tot.Area)%>%
      dplyr::group_by(Well_coordinate)%>%
      dplyr::arrange(timestep)%>%
      dplyr::mutate(GD_inst_Roll.Ratio.A.0.48h = GD_Tot.Area/first(GD_Tot.Area))%>%# 1.WORKS
      dplyr::mutate(Per.change = ((lead(GD_Tot.Area, n=1)/GD_Tot.Area)*100)-100)%>% # v1 roll 1 change in ratio over time
      drop_na()%>%
      mutate(Per.change = if_else(is.na(Per.change), 0,Per.change))%>% #v1
      ungroup()
    
    TP.moc.GD.df.Eval.inst.0.48h.df <- moc.GD.df.Eval.inst.0.48h.df %>%
      ungroup()%>%
      group_by(Well_coordinate,
               timestep)%>%
      dplyr::filter(Per.change <= -5)%>% #v2
      ungroup()%>%
      group_by(Well_coordinate)%>%
      arrange(timestep)%>%
      mutate(GD_timestep.of.det = (first(timestep)))%>% # The first instance where we have a deacrease by at least 5 % 
      select(Well_coordinate,
             GD_timestep.of.det)%>%
      distinct()
    
    moc.GD.df <- left_join(moc.GD.df, 
                           TP.moc.GD.df.Eval.inst.0.48h.df)
    
    #Houskeeping
    rm(moc.GD.df.Eval.inst.0.48h.df,
       TP.moc.GD.df.Eval.inst.0.48h.df)
    
    #Or
    
  }
  if ( x == "Det" && moc.GD.df$GD_Eval.0.48h == 0 ) {
    print("2")
    
    # between 48 to 96h. In TLCI experiments we know that the 11th frame is the 48h timepoint while the 12th frame the 96th h timepoint. Detachment detected here can only be in the 11th frame. 
    moc.GD.df.Eval.inst.48.96h.df <- moc.GD.df %>%
      dplyr::ungroup()%>%
      select(Well_coordinate,
             timestep)%>%
      group_by(Well_coordinate)%>%
      mutate(GD_timestep.of.det = timestep)%>%
      filter(timestep == max(timestep)-1)%>%
      select(-timestep)
    
    
    moc.GD.df <- left_join(moc.GD.df,
                           moc.GD.df.Eval.inst.48.96h.df)
    
    #Houskeeping
    rm(moc.GD.df.Eval.inst.48.96h.df)
    
  } 
  if ( x == "No-det"  ) {
    
    print("3")
    moc.GD.df <- moc.GD.df %>%
      mutate(GD_timestep.of.det = NA)
    
    
  }
  
  #Houskeeping
  rm(x)
  
  MS.res.df <- dplyr::left_join(MS.res.df,
                                moc.GD.df)
  
  #Housekeeping
  rm(moc.GD.df)
  
} else{ print("Gel detachmented-skipped")}


#9) Calculating live cell fraction POC
#left join poc.df 



# __ --------------------------
# Section 10.1 : LCF poc SC VS  Basic ff9 df0  --------------------------

poc.moc.LC.neg.sc.vs <- moc.poc.qc.df %>%
  dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class  == "VS")%>%
  dplyr::select(Well_coordinate,
                timestep,
                MOC.Predicted.Class,
                POC.Predicted.Class)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep,
                  POC.Predicted.Class) %>%
  summarise(n = n()) %>%
  dplyr::mutate(fraction = n/sum(n)) %>%
  dplyr::filter(POC.Predicted.Class== "piNEG")%>%
  dplyr::mutate(LC_piNEG.n.sc.vs = n)%>%
  dplyr::select(-n)%>%
  dplyr::mutate(LC_piNeg.fraction.sc.vs = fraction)%>%
  dplyr::select(-fraction,
                -POC.Predicted.Class)


poc.moc.LC.pos.sc.vs <- moc.poc.qc.df %>%
  dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class  == "VS")%>%
  dplyr::select(Well_coordinate,
                timestep,
                MOC.Predicted.Class,
                POC.Predicted.Class)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep,
                  POC.Predicted.Class) %>%
  summarise(n = n()) %>%
  dplyr::mutate(fraction = n/sum(n)) %>%
  dplyr::filter(POC.Predicted.Class== "piPOS")%>%
  dplyr::mutate(LC_piPOS.n.sc.vs = n)%>%
  dplyr::select(-n)%>%
  dplyr::mutate(LC_piPOS.fraction.sc.vs = fraction)%>%
  dplyr::select(-fraction,
                -POC.Predicted.Class)


MS.res.df <- left_join(MS.res.df,
                       poc.moc.LC.neg.sc.vs )


MS.res.df <- left_join(MS.res.df,
                       poc.moc.LC.pos.sc.vs)

rm( poc.moc.LC.neg.sc.vs,
    poc.moc.LC.pos.sc.vs )


# Section 10.1 : LCF poc SC VS  Basic default  --------------------------

poc.moc.LC.neg.sc.vs.BaSic <- moc.poc.qc.df %>%
  dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class  == "VS")%>%
  dplyr::select(Well_coordinate,
                timestep,
                MOC.Predicted.Class,
                POC.BaSic.Predicted.Class)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep,
                  POC.BaSic.Predicted.Class) %>%
  summarise(n = n()) %>%
  dplyr::mutate(fraction = n/sum(n))%>%
  dplyr::filter(POC.BaSic.Predicted.Class== "piNEG")%>%
  dplyr::mutate(LC_piNEG.n.sc.vs.BaSic = n)%>%
  dplyr::select(-n)%>%
  dplyr::mutate(LC_piNeg.fraction.sc.vs.BaSic = fraction)%>%
  dplyr::select(-fraction,
                -POC.BaSic.Predicted.Class)

poc.moc.LC.pos.sc.vs.BaSic <- moc.poc.qc.df %>%
  dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class  == "VS")%>%
  dplyr::select(Well_coordinate,
                timestep,
                MOC.Predicted.Class,
                POC.BaSic.Predicted.Class)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep,
                  POC.BaSic.Predicted.Class) %>%
  summarise(n = n()) %>%
  dplyr::mutate(fraction = n/sum(n)) %>%
  dplyr::filter(POC.BaSic.Predicted.Class== "piPOS")%>%
  dplyr::mutate(LC_piPOS.n.sc.vs.BaSic = n)%>%
  dplyr::select(-n)%>%
  dplyr::mutate(LC_piPOS.fraction.sc.vs.BaSic = fraction)%>%
  dplyr::select(-fraction,
                -POC.BaSic.Predicted.Class)




MS.res.df <- left_join(MS.res.df,
                       poc.moc.LC.neg.sc.vs.BaSic )


MS.res.df <- left_join(MS.res.df,
                       poc.moc.LC.pos.sc.vs.BaSic)

rm( poc.moc.LC.neg.sc.vs.BaSic,
    poc.moc.LC.pos.sc.vs.BaSic )



# __ --------------------------
# Section 10.2 : LCF poc  SC Basic ff9 df0   --------------------------

#------------ Live cell fraction of SC 
poc.moc.LC.neg.sc <- moc.poc.qc.df %>%
  dplyr::filter(MOC.Predicted.Class == "SC")%>%
  dplyr::select(Well_coordinate,
                timestep,
                MOC.Predicted.Class,
                POC.Predicted.Class)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep,
                  POC.Predicted.Class) %>%
  summarise(n = n()) %>%
  dplyr::mutate(fraction = n/sum(n)) %>%
  dplyr::filter(POC.Predicted.Class== "piNEG")%>%
  dplyr::mutate(LC_piNEG.n.sc = n)%>%
  dplyr::select(-n)%>%
  dplyr::mutate(LC_piNeg.fraction.sc = fraction)%>%
  dplyr::select(-fraction,
                -POC.Predicted.Class)

poc.moc.LC.pos.sc <- moc.poc.qc.df %>%
  dplyr::filter(MOC.Predicted.Class == "SC")%>%
  dplyr::select(Well_coordinate,
                timestep,
                MOC.Predicted.Class,
                POC.Predicted.Class)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep,
                  POC.Predicted.Class) %>%
  summarise(n = n()) %>%
  dplyr::mutate(fraction = n/sum(n)) %>%
  dplyr::filter(POC.Predicted.Class== "piPOS")%>%
  dplyr::mutate(LC_piPOS.n.sc = n)%>%
  dplyr::select(-n)%>%
  dplyr::mutate(LC_piPOS.fraction.sc = fraction)%>%
  dplyr::select(-fraction,
                -POC.Predicted.Class)



MS.res.df <- left_join(MS.res.df,
                       poc.moc.LC.neg.sc )


MS.res.df <- left_join(MS.res.df,
                       poc.moc.LC.pos.sc)

rm( poc.moc.LC.neg.sc,
    poc.moc.LC.pos.sc)

# Section 10.2 : LCF poc  SC  Basic default  --------------------------

poc.moc.LC.neg.sc.Basic <- moc.poc.qc.df %>%
  dplyr::filter(MOC.Predicted.Class == "SC")%>%
  dplyr::select(Well_coordinate,
                timestep,
                MOC.Predicted.Class,
                POC.BaSic.Predicted.Class)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep,
                  POC.BaSic.Predicted.Class) %>%
  summarise(n = n()) %>%
  dplyr::mutate(fraction = n/sum(n)) %>%
  dplyr::filter(POC.BaSic.Predicted.Class == "piNEG")%>%
  dplyr::mutate(LC_piNEG.n.sc.BaSic  = n)%>%
  dplyr::select(-n)%>%
  dplyr::mutate(LC_piNeg.fraction.sc.BaSic  = fraction)%>%
  dplyr::select(-fraction,
                -POC.BaSic.Predicted.Class)

poc.moc.LC.pos.sc.Basic <- moc.poc.qc.df %>%
  dplyr::filter(MOC.Predicted.Class == "SC")%>%
  dplyr::select(Well_coordinate,
                timestep,
                MOC.Predicted.Class,
                POC.BaSic.Predicted.Class)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  timestep,
                  POC.BaSic.Predicted.Class) %>%
  summarise(n = n()) %>%
  dplyr::mutate(fraction = n/sum(n)) %>%
  dplyr::filter(POC.BaSic.Predicted.Class== "piPOS")%>%
  dplyr::mutate(LC_piPOS.n.sc.BaSic = n)%>%
  dplyr::select(-n)%>%
  dplyr::mutate(LC_piPOS.fraction.sc.BaSic = fraction)%>%
  dplyr::select(-fraction,
                -POC.BaSic.Predicted.Class)



MS.res.df <- left_join(MS.res.df,
                       poc.moc.LC.neg.sc.Basic )


MS.res.df <- left_join(MS.res.df,
                       poc.moc.LC.pos.sc.Basic)

rm( poc.moc.LC.neg.sc.Basic ,
    poc.moc.LC.pos.sc.Basic )



# ___ --------------------------

# Section 11.1 : Threshold LCF poc Basic ff9 df0  --------------------------

#left join poc.df 
x <-  analyse.LC.oldThrs

if ( analyse.LC.oldThrs == "yes" ) {
  
  Death_threshold <- 700
  
  poc.moc.LC.neg.sc.vs <- moc.poc.qc.df %>%
    dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class  == "VS")%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  Mean.Intensity)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(POC.Predicted.Class.thrs = dplyr::if_else(Mean.Intensity >= Death_threshold, "piPOS", "piNEG"))%>%
    dplyr::group_by(Well_coordinate,
                    timestep,
                    POC.Predicted.Class.thrs) %>%
    summarise(n = n()) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::filter(POC.Predicted.Class.thrs == "piNEG")%>%
    dplyr::mutate(LC_piNEG.n.sc.vs.thrs = n)%>%
    dplyr::select(-n)%>%
    dplyr::mutate(LC_piNeg.fraction.sc.vs.thrs = fraction)%>%
    dplyr::select(-fraction,
                  -POC.Predicted.Class.thrs)
  
  
  poc.moc.LC.pos.sc.vs <- moc.poc.qc.df %>%
    dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class  == "VS")%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  Mean.Intensity)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(POC.Predicted.Class.thrs = dplyr::if_else(Mean.Intensity >= Death_threshold, "piPOS", "piNEG"))%>%
    dplyr::group_by(Well_coordinate,
                    timestep,
                    POC.Predicted.Class.thrs) %>%
    summarise(n = n()) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::filter(POC.Predicted.Class.thrs == "piPOS")%>%
    dplyr::mutate(LC_piNEG.n.sc.vs.thrs = n)%>%
    dplyr::select(-n)%>%
    dplyr::mutate(LC_piNeg.fraction.sc.vs.thrs = fraction)%>%
    dplyr::select(-fraction,
                  -POC.Predicted.Class.thrs)
  
  
  MS.res.df <- left_join(MS.res.df,
                         poc.moc.LC.neg.sc.vs )
  
  
  MS.res.df <- left_join(MS.res.df,
                         poc.moc.LC.pos.sc.vs)
  
  
  #------------ Live cell fraction of SC 
  
  poc.moc.LC.neg.sc <- moc.poc.qc.df %>%
    dplyr::filter(MOC.Predicted.Class == "SC")%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  Mean.Intensity)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(POC.Predicted.Class.thrs = dplyr::if_else(Mean.Intensity >= Death_threshold, "piPOS", "piNEG"))%>%
    dplyr::group_by(Well_coordinate,
                    timestep,
                    POC.Predicted.Class.thrs) %>%
    summarise(n = n()) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::filter(POC.Predicted.Class.thrs == "piNEG")%>%
    dplyr::mutate(LC_piNEG.n.sc.thrs= n)%>%
    dplyr::select(-n)%>%
    dplyr::mutate(LC_piNeg.fraction.sc.thrs = fraction)%>%
    dplyr::select(-fraction,
                  -POC.Predicted.Class.thrs)
  
  poc.moc.LC.pos.sc <- moc.poc.qc.df %>%
    dplyr::filter(MOC.Predicted.Class == "SC")%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  Mean.Intensity)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(POC.Predicted.Class.thrs = dplyr::if_else(Mean.Intensity >= Death_threshold, "piPOS", "piNEG"))%>%
    dplyr::group_by(Well_coordinate,
                    timestep,
                    POC.Predicted.Class.thrs) %>%
    summarise(n = n()) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::filter( POC.Predicted.Class.thrs == "piPOS")%>%
    dplyr::mutate(LC_piPOS.n.sc.thrs = n)%>%
    dplyr::select(-n)%>%
    dplyr::mutate(LC_piPOS.fraction.sc.thrs = fraction)%>%
    dplyr::select(-fraction,
                  -POC.Predicted.Class.thrs)
  
  
  
  MS.res.df <- left_join(MS.res.df,
                         poc.moc.LC.neg.sc )
  
  
  MS.res.df <- left_join(MS.res.df,
                         poc.moc.LC.pos.sc)
  
}

# Section 11.1 : Threshold LCF poc Basic Default --------------------------
if ( analyse.LC.oldThrs == "yes" ) {
  
  Death_threshold <- 700
  
  poc.moc.LC.neg.sc.vs.BaSicThr <- moc.poc.qc.df %>%
    dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class  == "VS")%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  Mean.IntensityBaSic)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(POC.BaSic.Predicted.Class.thrs = dplyr::if_else(Mean.IntensityBaSic >= Death_threshold, "piPOS", "piNEG"))%>%
    dplyr::group_by(Well_coordinate,
                    timestep,
                    POC.BaSic.Predicted.Class.thrs) %>%
    summarise(n = n()) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::filter(POC.BaSic.Predicted.Class.thrs == "piNEG")%>%
    dplyr::mutate(LC_piNEG.n.sc.vs.thrs.BaSic = n)%>%
    dplyr::select(-n)%>%
    dplyr::mutate(LC_piNeg.fraction.sc.vs.thrs.BaSic = fraction)%>%
    dplyr::select(-fraction,
                  -POC.BaSic.Predicted.Class.thrs)
  
  
  poc.moc.LC.pos.sc.vs.BaSicThr  <- moc.poc.qc.df %>%
    dplyr::filter(MOC.Predicted.Class == "SC" | MOC.Predicted.Class  == "VS")%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  Mean.IntensityBaSic)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(POC.BaSic.Predicted.Class.thrs = dplyr::if_else(Mean.IntensityBaSic >= Death_threshold, "piPOS", "piNEG"))%>%
    dplyr::group_by(Well_coordinate,
                    timestep,
                    POC.BaSic.Predicted.Class.thrs) %>%
    summarise(n = n()) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::filter(POC.BaSic.Predicted.Class.thrs == "piPOS")%>%
    dplyr::mutate(LC_piNEG.n.sc.vs.thrs = n)%>%
    dplyr::select(-n)%>%
    dplyr::mutate(LC_piNeg.fraction.sc.vs.thrs.BaSic = fraction)%>%
    dplyr::select(-fraction,
                  -POC.BaSic.Predicted.Class.thrs)
  
  
  MS.res.df <- left_join(MS.res.df,
                         poc.moc.LC.neg.sc.vs.BaSicThr  )
  
  
  MS.res.df <- left_join(MS.res.df,
                         poc.moc.LC.pos.sc.vs.BaSicThr  )
  
  
  #------------ Live cell fraction of SC 
  
  poc.moc.LC.neg.sc.BaSicThr  <- moc.poc.qc.df %>%
    dplyr::filter(MOC.Predicted.Class == "SC")%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  Mean.IntensityBaSic)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(POC.BaSic.Predicted.Class.thrs = dplyr::if_else(Mean.IntensityBaSic >= Death_threshold, "piPOS", "piNEG"))%>%
    dplyr::group_by(Well_coordinate,
                    timestep,
                    POC.BaSic.Predicted.Class.thrs) %>%
    summarise(n = n()) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::filter(POC.BaSic.Predicted.Class.thrs == "piNEG")%>%
    dplyr::mutate(LC_piNEG.n.sc.thrs.BaSic= n)%>%
    dplyr::select(-n)%>%
    dplyr::mutate(LC_piNeg.fraction.sc.thrs.BaSic = fraction)%>%
    dplyr::select(-fraction,
                  -POC.BaSic.Predicted.Class.thrs)
  
  poc.moc.LC.pos.sc.BaSicThr  <- moc.poc.qc.df %>%
    dplyr::filter(MOC.Predicted.Class == "SC")%>%
    dplyr::select(Well_coordinate,
                  timestep,
                  Mean.IntensityBaSic)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(POC.BaSic.Predicted.Class.thrs = dplyr::if_else(Mean.IntensityBaSic >= Death_threshold, "piPOS", "piNEG"))%>%
    dplyr::group_by(Well_coordinate,
                    timestep,
                    POC.BaSic.Predicted.Class.thrs) %>%
    summarise(n = n()) %>%
    dplyr::mutate(fraction = n/sum(n)) %>%
    dplyr::filter( POC.BaSic.Predicted.Class.thrs == "piPOS")%>%
    dplyr::mutate(LC_piPOS.n.sc.thrs.BaSic = n)%>%
    dplyr::select(-n)%>%
    dplyr::mutate(LC_piPOS.fraction.sc.thrs.BaSic = fraction)%>%
    dplyr::select(-fraction,
                  -POC.BaSic.Predicted.Class.thrs)
  
  
  
  MS.res.df <- left_join(MS.res.df,
                         poc.moc.LC.neg.sc.BaSicThr  )
  
  
  MS.res.df <- left_join(MS.res.df,
                         poc.moc.LC.pos.sc.BaSicThr )
  
}

rm(poc.moc.LC.neg.sc,
   poc.moc.LC.neg.sc.BaSicThr,
   poc.moc.LC.neg.sc.vs,
   poc.moc.LC.pos.sc.vs.BaSicThr,
   poc.moc.LC.pos.sc.Basic,
   poc.moc.LC.pos.sc.BaSicThr,
   poc.moc.LC.neg.sc.vs.BaSicThr,
   poc.moc.LC.pos.sc,
   poc.moc.LC.pos.sc.vs)

# ___ --------------------------

# Section 12 : Time frame of earliest event (growth, geldetachment or contamination ) --------------------------

# Image frame of timestep to ultimately filter by

# We have evaluated different partameters growth, geldetachment and contamination. We need to find the minimum timestep from these readouts by which we should filter the data

timestep.to.filter.by <- c(MS.res.df$GT_eval.min.tp.GT, # growth threhosld
                           MS.res.df$GD_timestep.of.det, # gel detachment 
                           MS.res.df$ Contamination_timestep)# Contamination


# Remove instance with NA while it is still character format
timestep.to.filter.by <- timestep.to.filter.by[!is.na(timestep.to.filter.by)]
timestep.to.filter.by

timestep.to.filter.by <- as.numeric(timestep.to.filter.by)

timestep.to.filter.by


# HERE REMOVE NA did not work
timestep.to.filter.by

x <- length(timestep.to.filter.by)
x 
if (x == 0 ){
  
  print("This well passed all quality criteria")
  
  Filterting.timestep <- as.numeric(NA)
  
  MS.res.df <- MS.res.df %>%
    mutate(Filterting.timestep =Filterting.timestep )
  
  MS.res.df$Filterting.timestep <- as.numeric(MS.res.df$Filterting.timestep )
}

timestep.to.filter.by

if(x > 0) {
  
  print("This well has to be filtered")
  
  Filterting.timestep <- unique(min(timestep.to.filter.by))
  
  MS.res.df <- MS.res.df %>%
    mutate(Filterting.timestep =Filterting.timestep )
  
  MS.res.df$Filterting.timestep <- as.numeric(MS.res.df$Filterting.timestep )
}



#Housekeeping
rm(x)

setwd(wdDir)
#We also need to know hy what variable/criteria are we filtering by.Was it growth? gel detachment or contaminaiton? what event happended the soonest? (if it did)

Variable.for.filtering <- MS.res.df %>%
  filter(timestep ==Filterting.timestep)


x <-nrow(Variable.for.filtering)
x
if ( x  == 0) {
  print("This well passed all quality criteria")

  timestep.Filtering.Event <- NA

  MS.res.df <- MS.res.df %>%
    dplyr::mutate(timestep.Filtering.Event = timestep.Filtering.Event )

  MS.res.df$timestep.Filtering.Event <- as.character(MS.res.df$timestep.Filtering.Event)

}

if( x >= 1 && analyse.GelDetachment == "yes" ) {

  timestep.to.filter.by <- unique(timestep.to.filter.by)

  timestep.to.filter.by <- min(timestep.to.filter.by)
  
  timestep.to.filter.by

  variable.for.filtering.df <- MS.res.df %>%
    dplyr::ungroup()%>%
    filter(timestep == unique(timestep.to.filter.by))
  
  variable.for.filtering.df$timestep <- as.character( variable.for.filtering.df$timestep)
  variable.for.filtering.df$Contamination_timestep <- as.character(variable.for.filtering.df$Contamination_timestep)
  variable.for.filtering.df$GT_eval.min.tp.GT <- as.character(variable.for.filtering.df$GT_eval.min.tp.GT)
  variable.for.filtering.df$GD_timestep.of.det  <- as.character(variable.for.filtering.df$GD_timestep.of.det )
  
  variable.for.filtering.df <-  variable.for.filtering.df %>%
  select(timestep,
         GT_Growth.eval.in.well,
         GT_eval.min.tp.GT,
         Contamination_timestep,
         Contamination_eval,
         GD_timestep.of.det,
         GD_eval.gelDet)
  

  
  variable.for.filtering.df <-variable.for.filtering.df %>%
    mutate_if(is.character, ~replace_na(.,"NA"))
  
  variable.for.filtering.df <-variable.for.filtering.df %>%
  mutate(timestep.Filtering.Event = dplyr::if_else(timestep == GT_eval.min.tp.GT,
                                                   paste(GT_Growth.eval.in.well),
                                                   dplyr::if_else(timestep == GD_timestep.of.det,
                                                                  paste(GD_eval.gelDet),
                                                                  dplyr::if_else(timestep == Contamination_timestep,
                                                                                 paste(Contamination_eval), "error"))))
  
  
  
  variable.for.filtering <- variable.for.filtering.df$timestep.Filtering.Event

  MS.res.df <- MS.res.df %>%
    mutate(timestep.Filtering.Event = variable.for.filtering)

  #Houskeeping
  rm( variable.for.filtering,
      variable.for.filtering.df )

}



if( x >= 1 && analyse.GelDetachment == "no" ) {
  
  timestep.to.filter.by <- unique(timestep.to.filter.by)
  
  timestep.to.filter.by <- min(timestep.to.filter.by)
  
  timestep.to.filter.by
  
  variable.for.filtering.df <- MS.res.df %>%
    dplyr::ungroup()%>%
    filter(timestep == unique(timestep.to.filter.by))
  
  variable.for.filtering.df$timestep <- as.character( variable.for.filtering.df$timestep)
  variable.for.filtering.df$Contamination_timestep <- as.character(variable.for.filtering.df$Contamination_timestep)
  variable.for.filtering.df$GT_eval.min.tp.GT <- as.character(variable.for.filtering.df$GT_eval.min.tp.GT)
 # variable.for.filtering.df$GD_timestep.of.det  <- as.character(variable.for.filtering.df$GD_timestep.of.det )
  
  variable.for.filtering.df <-  variable.for.filtering.df %>%
    select(timestep,
           GT_Growth.eval.in.well,
           GT_eval.min.tp.GT,
           Contamination_timestep,
           Contamination_eval#,
        #   GD_timestep.of.det,
         #  GD_eval.gelDet
        )
  
  
  
  variable.for.filtering.df <-variable.for.filtering.df %>%
    mutate_if(is.character, ~replace_na(.,"NA"))
  
  variable.for.filtering.df <-variable.for.filtering.df %>%
    mutate(timestep.Filtering.Event = dplyr::if_else(timestep == GT_eval.min.tp.GT,
                                                     paste(GT_Growth.eval.in.well),
                                                                    dplyr::if_else(timestep == Contamination_timestep,
                                                                                   paste(Contamination_eval), "error")))
  
  
  
  variable.for.filtering <- variable.for.filtering.df$timestep.Filtering.Event
  
  MS.res.df <- MS.res.df %>%
    mutate(timestep.Filtering.Event = variable.for.filtering)
  
  #Houskeeping
  rm( variable.for.filtering,
      variable.for.filtering.df)
  
}

#Housekeeping
rm(Variable.for.filtering)

# Section 13 : Live cell fractions can never increase  --------------------------

# Here we will make a copy of the raw  live cell fraction values for each prameter and we will correct for instances the live cell fraction goes up.
# We will consider the live cell fraction of the previous frame if the LCF increases.


variables_to_select <- names(MS.res.df)
variables_to_select <- variables_to_select[grep("LC", variables_to_select)]

#By using " {col}.RAW", you are adding ".RAW" at the end of each selected column name.
MS.res.df <-MS.res.df %>%
  mutate(across(all_of(variables_to_select), ~ ., .names = "{col}.RAW"))


rm(variables_to_select)
variables_to_select <- names(MS.res.df)
variables_to_select <- variables_to_select[grep("LC", variables_to_select)]
variables_to_select <- variables_to_select[!grepl("RAW", variables_to_select)]
variables_to_select <- variables_to_select[grepl("LC_piNeg", variables_to_select)]

LCF.res.df <-MS.res.df %>%
  select(ExpFile,
         Well_coordinate,
         timestep,
     #    Time_Hrs,
         all_of(variables_to_select))

rm(variables_to_select)

library(reshape2)


# Correct for increasing LCF fraction. If there is an incease from one frame to another then take result of previous frame.
LCF.res.df.melted  <- melt(LCF.res.df, 
                           id.vars = c("ExpFile",
                                       "Well_coordinate",
                                       "timestep"),
                           variable.name = "LCF.def",
                           value = "LCF")

LCF.res.df.melted <- LCF.res.df.melted %>%
  ungroup() %>%
  group_by(ExpFile, Well_coordinate, LCF.def) %>%
  mutate(
    LCF.def.Increase.corrected = ifelse(
      value > lag(value, default = first(value)),
      lag(value, default = first(value)),
      value
    )
  )

# Assuming your dataframe is named LCF.res.df.melted
max_iterations <- max(MS.res.df$timestep)-1  # Set a maximum number of iterations to avoid infinite loops. Increases can occure multiple time. Therefore we apply this correction for as many loops time frames

for (iteration in 1:max_iterations) {
  LCF.res.df.melted <- LCF.res.df.melted %>%
    arrange(ExpFile, Well_coordinate, LCF.def,     timestep) %>%
    group_by(ExpFile, Well_coordinate, LCF.def) %>%
    mutate(
      LCF.def.Increase.corrected = ifelse(
        LCF.def.Increase.corrected > lag(    LCF.def.Increase.corrected, default = first(    LCF.def.Increase.corrected)),
        lag(    LCF.def.Increase.corrected, default = first(    LCF.def.Increase.corrected)),
        LCF.def.Increase.corrected
      )
    ) %>%
    ungroup()
  
  # Check if there are any increases left, break the loop if none
  if (!any(LCF.res.df.melted$    LCF.def.Increase.corrected> lag(LCF.res.df.melted$LCF.def.Increase.corrected, default = LCF.res.df.melted$    LCF.def.Increase.corrected[1]))) {
    break
  }
}

# RUN THIS CODE TO CHECK IF THE ABPVE CODE IS CORRECT
# LCF.res.df.melted <- LCF.res.df.melted %>%
#   group_by(ExpFile, Well_coordinate, LCF.def) %>%
#   arrange(Time_Hrs)%>%
#    mutate(Finindg.Inc = value - lead(value))%>% ## Result should be either positive and zero. If it is negative, it reflects an increase
#   mutate(Finindg.Inc.Correded =    LCF.def.Incorrected - lead(   LCF.def.Incorrected))


LCF.res.df.melted <- LCF.res.df.melted %>%
  ungroup()%>%
  select(ExpFile,
         Well_coordinate,
       #  Time_Hrs,
       timestep,
         LCF.def,
         LCF.def.Increase.corrected)

# ERROR I think the line below does not work in scicore
 LCF.res.df.spread <- tidyr::spread(LCF.res.df.melted, key = "LCF.def",
                             value ="LCF.def.Increase.corrected")

#
#LCF.res.df.spread <- tidyr::pivot_wider(LCF.res.df.melted, names_from = LCF.def, values_from = LCF.def.Increase.corrected)


variables_to_select <- names(LCF.res.df.spread)
variables_to_select <- variables_to_select[grep("LC", variables_to_select)]

# Dropping the uncorrected variables from main results table
MS.res.df <- MS.res.df %>%
  select(-all_of(variables_to_select))


# Plugging back corrected results
MS.res.df <- left_join(MS.res.df,
                       LCF.res.df.spread,
                       by = c("ExpFile",
                              "Well_coordinate",
                              "timestep"))

# # Checking results in a plot
# check.melted <-MS.res.df %>%
#   select(Time_Hrs,
#          LC_piNeg.fraction.sc.thrs,
#          LC_piNeg.fraction.sc.thrs.RAW,
#          LC_piNeg.fraction.sc.thrs.BaSic,
#          LC_piNeg.fraction.sc.thrs.BaSic.RAW,
#          LC_piNeg.fraction.sc.vs,
#          LC_piNeg.fraction.sc.vs.RAW ,
#          LC_piNeg.fraction.sc.vs.BaSic,
#          LC_piNeg.fraction.sc.vs.BaSic.RAW )
# 
# 
# check.melted <- melt(check.melted , id.vars = c("Time_Hrs"))
# 
# check.melted %>%
#   ggplot(aes(x = Time_Hrs,
#              y = value,
#              group = variable,
#              colour = variable))+
#   #geom_point()+
#   ylim(0,1)+
#   geom_line()

# Section 14 : Export table of results --------------------------

# Saving results table with R script as a numaing variable

MS.res.df <- MS.res.df %>%
  select(-MOC.Predicted.Class)

filename <- gsub(".R","",RscriptUsed)

filename <- paste(filename, 
                  "_",
                  Exp.ID,
                  "_",
                  unique(MS.res.df$Well_coordinate),
"_results.csv",
sep="")


filename.path <- paste(resDir,
                       "/",
                       filename,
                       sep="")

write.csv(MS.res.df, filename.path,
          row.names = FALSE)


print ("--- COGRATULATIONS! The analysis is complete  ---")










