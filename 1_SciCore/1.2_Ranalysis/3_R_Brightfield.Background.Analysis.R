
#1) IMPORANT: user input vairables

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
                              pattern="OCsimple.bfint",
                              all.files=FALSE,
                              full.names=FALSE)


Exp.Well.to.Analyse <- as.data.frame(list.of.folders)

#Housekeeping
rm(list.of.folders)

#Creating a dataframe with directory details
Exp.Well.to.Analyse <- separate(Exp.Well.to.Analyse, col = list.of.folders,
                                into = c("ExpID",
                                         "Well",
                                         "Data"),
                                sep="_")


RscriptUsed <- list.files(path=".", 
                          pattern="3_R_Brightfield.Background.Analysis.R",
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
bf.inputDir <-paste(genDir,
                             "/",
                             expID,
                             "_",
                             well.being.analysed,
                             "_OCsimple.bfint",
                             sep="")

#----- <<<<<>>>>>>>> Analysis <<<<<>>>>>>>> ----------

# Loading data 

  OC.filenames <- list.files(path = bf.inputDir,
                             pattern = "*.csv",
                             full.names = FALSE)

#setting path to where all the csv files due to be processed are
setwd(bf.inputDir)

OC.file.list <- lapply(OC.filenames,
                       read.csv)

setwd(genDir)
#Naming the elements of my list with the vector list of my file names
Exp.Names <- gsub(".csv","",OC.filenames)
Exp.Names <- gsub("OCsimple_BFbckint_","", Exp.Names)
names(OC.file.list)<- Exp.Names

#Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 

BF.bck.Int.OC.df <- bind_rows(OC.file.list, .id = "ExpFile_Well_Field")

#Housekeeping
rm(OC.file.list)

BF.bck.Int.OC.df <- separate(BF.bck.Int.OC.df, col =ExpFile_Well_Field,
                  into = c("ExpFile",
                           "Well_coordinate",
                           "Field"),
                  sep="_")


BF.bck.Int.OC.df  <- BF.bck.Int.OC.df  %>%
  filter(Predicted.Class == "Background")%>%
  select(-labelimage_oid)
  
BF.bck.Int.OC.df  <- BF.bck.Int.OC.df  %>%
  dplyr::ungroup()%>%
  dplyr::group_by(Well_coordinate,
                  Field,
                  timestep)%>%
  filter(Size.in.pixels == max(Size.in.pixels))

BF.bck.Int.OC.df  <- BF.bck.Int.OC.df  %>%
  dplyr::ungroup()
  
# Selecting background only and calculating mean backgroubnd bright field intensity per field and per well
BF.bck.Int.OC.df  <- BF.bck.Int.OC.df  %>%
  group_by(ExpFile,
           Well_coordinate,
           Field,
           timestep)%>%
  mutate(Mean.Brightfield.Background.Int = mean(Mean.Intensity))%>%
  select(ExpFile,
         Well_coordinate,
         Field,
         timestep,
         Mean.Brightfield.Background.Int )%>%
  distinct()


# Exporting results
# Saving results table with R script as a naming variable
experiment.well.filename  <- paste("ASCT_BFint_results",
                          unique(BF.bck.Int.OC.df$ExpFile),
                          unique(BF.bck.Int.OC.df$Well_coordinate),
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

write.csv(BF.bck.Int.OC.df,
                 file.name.path,
          row.names = FALSE)



print ("--- COGRATULATIONS! Brightfield object anlysis is complete ---")









