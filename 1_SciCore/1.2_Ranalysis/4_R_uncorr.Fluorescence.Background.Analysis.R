


#1) User input variables

#----  <<<<<<<<<<>>>>>> ---- 
#User Input variables
expected.N.fields <- 1:9
expected.N.timestep <- 0:29

#----  <<<<<<<<<<>>>>>> ---- 


# Loading libraries 
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(RColorBrewer)
library(platetools)
library(directlabels) 
library(MESS)

#Identifying well to analyze

# Find all files in the current directory this Rscipt is found in
list.of.folders <- list.files(path=".", 
                              pattern="OCsimple.boc",
                              all.files=FALSE,
                              full.names=FALSE)

print(list.of.folders)
Exp.Well.to.Analyse <- as.data.frame(list.of.folders)

#Housekeeping
rm(list.of.folders)
#Creating a dataframe with directory details
Exp.Well.to.Analyse <- tidyr::separate(Exp.Well.to.Analyse, col = list.of.folders,
                                into = c("ExpID",
                                         "Well",
                                         "Data"),
                                sep="_")

RscriptUsed <- list.files(path=".", 
                          pattern="4_R_uncorr.Fluorescence.Background.Analysis.R",
                          all.files=FALSE,
                          full.names=FALSE)


well.being.analysed <- unique(Exp.Well.to.Analyse$Well)


#Creating an Experiment variable e.g. TLCI.13.2021125
expID<- unique(Exp.Well.to.Analyse$ExpID)



#--- Setting directories

# working directory
genDir <- getwd() 

# result directory
resDir <- paste(genDir,
                "/",
                expID,
                "_",
                well.being.analysed,
                "_Results",
                sep="")

# input data directory
fluncorr.inputDir <-paste(genDir,
                             "/",
                             expID,
                             "_",
                             well.being.analysed,
                             "_OCsimple.boc",
                             sep="")

#----- <<<<<>>>>>>>> HM of Fluorescence background intensity <<<<<>>>>>>>>----------



# Loading data 

  OC.filenames <- list.files(path = fluncorr.inputDir,
                             pattern = "*.csv",
                             full.names = FALSE)

#setting path to where all the csv files due to be processed are
setwd(fluncorr.inputDir)

OC.file.list <- lapply(OC.filenames,
                       read.csv)

setwd(genDir)
#Naming the elements of my list with the vector list of my file names
Exp.Names <- gsub(".csv","",OC.filenames)
Exp.Names <- gsub("OCsimple_BOC_","", Exp.Names)
names(OC.file.list)<- Exp.Names

#Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 

FL.bck.Int.OC.df <- bind_rows(OC.file.list, .id = "ExpFile_Well_Field")

#Housekeeping
rm(OC.file.list)

FL.bck.Int.OC.df <- separate(FL.bck.Int.OC.df, col =ExpFile_Well_Field,
                  into = c("ExpFile",
                           "Well_coordinate",
                           "Field"),
                  sep="_")


FL.bck.Int.OC.df  <- FL.bck.Int.OC.df  %>%
  dplyr::filter(Predicted.Class == "Background")%>%
  select(-labelimage_oid)
  

# Selecting background only and calculating mean backgroubnd bright field intensity per field and per well
FL.bck.Int.OC.df  <- FL.bck.Int.OC.df  %>%
  group_by(ExpFile,
           Well_coordinate,
           Field,
           timestep)%>%
  mutate(Mean.FLuncorr.Background.Int = mean(Mean.Intensity))%>%
  select(ExpFile,
         Well_coordinate,
         Field,
         timestep,
         Mean.FLuncorr.Background.Int )%>%
  distinct()


# Exporting results
# Saving results table with R script as a naming variable
experiment.well.filename  <- paste("ASCT_FLuncorr_results",
                          unique(FL.bck.Int.OC.df$ExpFile),
                          unique(FL.bck.Int.OC.df$Well_coordinate),
                          sep="_")

experiment.well.filename <- paste(experiment.well.filename,
                                  ".csv",
                                  sep="")


experiment.well.filename 

file.name.path <- paste(resDir,
                        "/",
                        experiment.well.filename,
                        sep="")

file.name.path 

write.csv(FL.bck.Int.OC.df,
                 file.name.path,
          row.names = FALSE)

print ("--- COGRATULATIONS! Fluorescence background anlysis is complete ---")









