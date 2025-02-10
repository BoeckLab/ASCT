
# Aim ------------------------


# This script will go in to the tracking labeloid directory which is Ahmads tracking results and merge the data against the Basic FF9DF0 (old background correction) and new BaSic default background correction.
# It will then find the LCF across these definitions.

# Section 1: defining paths and directories  ------------------------

#General variables
genDir <- getwd() 

# Single Master dataframe to pool results across each BaSic definition
MASTER.df.Tracking_pi <- data.frame()
MASTER.df.Tracking.N <- data.frame()



#Creacting vector with the list of file names
filenames.list <- list.files(path = genDir,
                             pattern = "Tracking_LC",
                             full.names = FALSE)


filenames.list <- gsub("_Tracking_LC","",filenames.list)



res <- c("PerWell-Results")
exp.res <- c("Experimental-Results")

tracking.res <- c("Tracking_LC") 

tracking.res.export <- c(tracking.res)

#Directories
wdDir <- genDir

perWell.resDir <- paste(genDir,"/",res,sep="")
exp.resDir <- paste(genDir,"/", exp.res,sep="")

tracking.resDir <- paste(genDir, "/",tracking.res , sep="" )
#Housekeeping
rm(mocSimple,
   pocSimple,
   qc,
   metaData,
   res)

#Setting work directory
setwd(wdDir)


# Section 2: Loading libraries  ------------------------

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
library(gghighlight)
library(reshape2)



#--- Start Tracking data ---

trk.res <- paste(filenames.list,
                 "_Tracking_labeloid/",
                 sep="")


ila.res.poc <- paste(filenames.list,
                     "_OCsimple.poc",
                     sep="")

ila.res.poc.BaSic <- paste(filenames.list,
                           "_OCsimple.pocbSd",
                           sep="")

tracking.res.export <- paste(filenames.list,
                             "_Tracking_LC",
                             sep="")

#Directories
wdDir <- genDir

trk.resDir <- paste(genDir,"/", trk.res,sep="") # Tracking data
ila.res.poc.resDir <- paste(genDir,"/", ila.res.poc,sep="") # POC old baSic
ila.res.poc.BaSic.resDir <- paste(genDir,"/",ila.res.poc.BaSic,sep="") # POC new BaSic
tracking.res.export <- paste(genDir,"/",tracking.res.export,sep="")# Export directory


# Section 2.1 : Detecting whether BaSic default data should be considered ----
# Define the directory path
directory_path <- ila.res.poc.BaSic.resDir


# Check if the directory exists
if (dir.exists(directory_path)) {
  # List files in the directory
  files <- list.files(directory_path)
  
  # Check if there are files present
  if (length(files) > 0) {
    cat("Files are present in the directory:", directory_path, "\n")
    cat("List of files:\n")
    #  cat(files, sep = "\n")
    
    list.of.POC.result.dir <- c(ila.res.poc.resDir, ila.res.poc.BaSic.resDir)
  } else {
    cat("Directory is present but no files found in:", directory_path, "\n")
    
    list.of.POC.result.dir <- c(ila.res.poc.resDir)
    
    
  }
} else {
  cat("Directory does not exist:", directory_path, "\n")
  
  list.of.POC.result.dir <- c(ila.res.poc.resDir)
}

list.of.POC.result.dir

rm(files,
   directory_path)



# Section 2.2: Checking which tracking versions are available  ----

# Define the parent directory where you want to start the search
parent_directory <- genDir

# List all directories in the parent directory
subdirectories <- list.dirs(parent_directory, full.names = TRUE)

# Filter directories containing "Tracking_labeloid" in their name
matching_directories <- subdirectories[grep("Tracking_labeloid", subdirectories)]

# Loop through matching directories and list files
for (directory_path in matching_directories) {
  cat("Files in directory:", directory_path, "\n")
  Tracking.files <- list.files(directory_path)
  
  if (length(  Tracking.files) > 0) {
    cat("List of files:\n")
    cat(  Tracking.files,
          sep = "\n")
    
  } else {
    cat("No files found in the directory.\n")
  }
  cat("\n")
}

rm(parent_directory)



# here we select only the data from old tracking algorithm
Tracking.files <- Tracking.files[grep("Trkv2", Tracking.files )]

Tracking.version <- gsub(".csv","",gsub(".*-","",Tracking.files))


#poc.director.of.interest  <- list.of.POC.result.dir[1]

# Section  Loop ------

progress <- 0
# list.of.POC.result.dir <- list.of.POC.result.dir[2]
#i <- list.of.POC.result.dir[1]


for ( i in list.of.POC.result.dir) {
  
  progress <- progress + 1
  print(i)
  poc.director.of.interest  <- i
  
  list.of.POC.result.dir.subset <- i
  
  BaSic.parameter <- sub(".*_OCsimple\\.", "", list.of.POC.result.dir.subset)
  
  print(BaSic.parameter )
  
  #   #Creacting vector with the list of file names
  perWell.filenames <- list.files(path = trk.resDir ,
                                  pattern =Tracking.files  ,
                                  full.names = FALSE)
  
  #setting path to where all the csv files due to be processed are
  setwd(trk.resDir)
  
  OC.file.list <- lapply(perWell.filenames,
                         read.csv)
  
  
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub(".csv",
                    "",
                    perWell.filenames)
  
  names(OC.file.list)<- Exp.Names
  
  #Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
  Trk.perWell.df <- rapply(
    OC.file.list,
    as.character,
    how = "replace"
  )
  
  
  Trk.perWell.df  <- bind_rows(Trk.perWell.df  , .id = "Filename")
  
  
  # Section 3.2 : wrangling tracking dataframe  ----
  
  Trk.perWell.df <- Trk.perWell.df %>%
    mutate(TrackID = TRK_ID,
           Labelimage_oid = ilastik_ID,
           timestep = Time)%>%
    select(-TRK_ID)
  
  
  #Trk.perWell.df$Filename <- gsub("Dynamic_Table_corrected_","",Trk.perWell.df$Filename)
  
  Trk.perWell.df$Filename <- gsub(".*corrected_","",Trk.perWell.df$Filename)
  
  
  
  Trk.perWell.df <- separate(Trk.perWell.df, col = well_coord,
                             into = c("Exp",
                                      "Well_coordinate",
                                      "Field"),
                             sep = "_")
  
  
  
  
  
  perWell.trk.pi.df <-  Trk.perWell.df
  
  #HOuskeeping
  rm(Trk.perWell.df)
  
 
  
  #Creacting vector with the list of file names
  OC.filenames <- list.files(path = list.of.POC.result.dir.subset,
                             pattern = "*.csv",
                             full.names = FALSE)
  
  
  
  #setting path to where all the csv files due to be processed are
  setwd(list.of.POC.result.dir.subset)
  
  OC.file.list <- lapply(OC.filenames,
                         read.csv)
  
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub(".csv","",OC.filenames)
  names(OC.file.list)<- Exp.Names
  
  #Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
  
  OC.df <- bind_rows(OC.file.list, .id = "Ilastik_Classification_ExpFile_Well_Field")
  
  OC.df <- OC.df %>%
    mutate(Ila.Mean.Intensity = Mean.Intensity )%>%
    select(-Mean.Intensity)
  
  
  OC.df$Ilastik_Classification_ExpFile_Well_Field <- gsub(".*POC_","",OC.df$Ilastik_Classification_ExpFile_Well_Field )
  OC.df$Ilastik_Classification_ExpFile_Well_Field <- gsub(".*POCbSd_","",OC.df$Ilastik_Classification_ExpFile_Well_Field )
  
  
  
  OC.df <- separate(OC.df, col = Ilastik_Classification_ExpFile_Well_Field,
                    into = c(#"Ilastik",
                      #"Classification",
                      "ExpFile",
                      "Well_coordinate",
                      "Field"),
                    sep="_")
  
  OC.df <- OC.df %>%
    dplyr::mutate(POC.Predicted.Class = Predicted.Class)%>%
    dplyr::select(#-Ilastik,
      #-Classification,
      -Predicted.Class)
  
  
  # Previously When running the above code on larger number of csv files (7000+) f0 ,f1, f2 columns were being produced which had NAs as values. We therefore want to remove these colunms from the OC.df
  OC.df <- OC.df %>% select_if(~!all(is.na(.))) 
  
  
  #House keeping
  rm(OC.file.list,
     OC.filenames,
     Exp.Names
  )
  
  
  
  
  
  
  
  perWell.trk.pi.df$Time <- as.numeric(perWell.trk.pi.df$Time)
  
  perWell.trk.pi.df <- perWell.trk.pi.df %>%
    mutate(timestep = Time-1)%>%
    select(-Time)%>%
    mutate(labelimage_oid = Labelimage_oid)%>%
    select(-Labelimage_oid)%>%
    mutate(ExpFile = Exp)%>%
    select(-Exp)
  
  #Changing class of colunms of df 
  perWell.trk.pi.df$labelimage_oid <- as.numeric(perWell.trk.pi.df$labelimage_oid)
  OC.df$labelimage_oid <- as.numeric(OC.df$labelimage_oid)
  
  perWell.trk.pi.df  <- left_join(perWell.trk.pi.df,
                                  OC.df)
  
  perWell.trk.pi.df  <-   perWell.trk.pi.df %>%
    select(-Center.of.the.object_0,
           -Center.of.the.object_1,
           -labelimage_oid,
           -ilastik_ID)
  
  #Housekeeping
  rm(OC.df)
  # Section 3.2.1 : Bridging the gap of missing frames  ----
  
  # Ahmad, has defined in his tracking algorithm that if a trackID is missing in more than 4 frames, it will not consider the ID in analysis.
  # When I now account for this by adopting the last value of the previous POC predicted results or Mean intenity i need to run this code at least 4 times in the event that a trackID has missing information in 4 sequential image frames.
  
  
  perWell.trk.pi.df <- perWell.trk.pi.df%>%
    ungroup() %>%
    group_by(Well_coordinate, Field, TrackID) %>%
    mutate(NAs_in_Ila_Mean_Intensity = sum(is.na(Ila.Mean.Intensity)))%>%
  mutate(Consecutive_NAs = ifelse(is.na(Ila.Mean.Intensity) & 
                                         (is.na(lag(Ila.Mean.Intensity)) | 
                                            is.na(lead(Ila.Mean.Intensity))), 
                                       1,0))%>%
    mutate(Sum.of.Consecutive_NAs = sum(Consecutive_NAs))%>%
    ungroup()%>%
    filter(Sum.of.Consecutive_NAs <= 2)%>%
    filter(NAs_in_Ila_Mean_Intensity <= 6)%>%
    
    ungroup()%>%
    group_by(Well_coordinate,
             Field,
             TrackID)%>% # MF abbreviation for Missing Frames
    mutate(Ila.Mean.Intensity.MF1 = if_else(is.na(Ila.Mean.Intensity), lag(Ila.Mean.Intensity), Ila.Mean.Intensity),
           Ila.Mean.Intensity.MF2 = if_else(is.na(Ila.Mean.Intensity.MF1 ), lag(Ila.Mean.Intensity.MF1 ), Ila.Mean.Intensity.MF1),
           Ila.Mean.Intensity.MF3 = if_else(is.na(Ila.Mean.Intensity.MF2 ), lag(Ila.Mean.Intensity.MF2 ), Ila.Mean.Intensity.MF2),
           Ila.Mean.Intensity.MF4 = if_else(is.na(Ila.Mean.Intensity.MF3 ), lag(Ila.Mean.Intensity.MF3 ), Ila.Mean.Intensity.MF3),
           Ila.Mean.Intensity.MF5 = if_else(is.na(Ila.Mean.Intensity.MF4 ), lag(Ila.Mean.Intensity.MF4 ), Ila.Mean.Intensity.MF4),
           Ila.Mean.Intensity.MF6 = if_else(is.na(Ila.Mean.Intensity.MF5 ), lag(Ila.Mean.Intensity.MF5 ), Ila.Mean.Intensity.MF5))%>%
    mutate(POC.Predicted.Class.MF1 = if_else(is.na(POC.Predicted.Class), lag(POC.Predicted.Class), POC.Predicted.Class),
           POC.Predicted.Class.MF2 = if_else(is.na(POC.Predicted.Class.MF1 ), lag(POC.Predicted.Class.MF1 ), POC.Predicted.Class.MF1),
           POC.Predicted.Class.MF3 = if_else(is.na(POC.Predicted.Class.MF2 ), lag(POC.Predicted.Class.MF2 ), POC.Predicted.Class.MF2),
           POC.Predicted.Class.MF4 = if_else(is.na(POC.Predicted.Class.MF3 ), lag(POC.Predicted.Class.MF3 ), POC.Predicted.Class.MF3),
           POC.Predicted.Class.MF5 = if_else(is.na(POC.Predicted.Class.MF4 ), lag(POC.Predicted.Class.MF4 ), POC.Predicted.Class.MF4),
           POC.Predicted.Class.MF6 = if_else(is.na(POC.Predicted.Class.MF5 ), lag(POC.Predicted.Class.MF5 ), POC.Predicted.Class.MF5))%>%
    select(Filename,
           Well_coordinate,
           Field,
           TrackID,
           timestep,
           ExpFile,
           Ila.Mean.Intensity,
           Ila.Mean.Intensity.MF6,
           POC.Predicted.Class,
           POC.Predicted.Class.MF6)%>%
    mutate(Ila.Mean.Intensity = Ila.Mean.Intensity.MF6,
           POC.Predicted.Class = POC.Predicted.Class.MF6)%>%
    select(-Ila.Mean.Intensity.MF6,
           -POC.Predicted.Class.MF6)%>%
    mutate(Filename  = gsub(paste("-",Tracking.version,sep=""),"",Filename))
  
  
 
  #Housekeeping
  rm(OC.df)
  
  # Section 3.3 : Defining live cell definition across all tracking definitions  ----
  
  
  # Calculate the live cell fraction  across all tracking definitions 
  # Thresholding: 1) >= 700 mean intesity & once positive always positive  2) >=700 mean intensity & twice psisive always positive
  
  # Section 3.3.1 : Thresholding 1 pos   ----
  
  # Thrs 1 : Once positive always positive
  perWell.trk.pi.df <- perWell.trk.pi.df %>%
    dplyr::ungroup()%>%
    group_by(Filename,
             Well_coordinate,
             Field,
             TrackID)%>%
    mutate(Trk_Thrs.1pos_LCfraction = dplyr::if_else(Ila.Mean.Intensity >= 700, 1, 0),
           Trk_Thrs.1pos_LCfraction = cumsum(Trk_Thrs.1pos_LCfraction),
           Trk_Thrs.1pos_LCfraction = if_else(Trk_Thrs.1pos_LCfraction >= 1, "piPOS", "piNEG"))
  
  
  # Section 3.3.2 : Thresholding 2 pos   ----
  
  # --- Thrs 2 : Twice positive always positive always positive
  
  perWell.trk.pi.df <- perWell.trk.pi.df %>%
    dplyr::ungroup()%>%
    group_by(Filename,
             Well_coordinate,
             Field,
             TrackID)%>%
    mutate(Trk_Thr.2pos_LCfraction = dplyr::if_else(Ila.Mean.Intensity >= 700, 1, 0))%>%
    mutate(Trk_Thr.2.2pos_LCfraction = if_else(lead(Trk_Thr.2pos_LCfraction) == 1 & Trk_Thr.2pos_LCfraction == 1,1,0))%>%
    mutate(Trk_Thr.2.2.2pos_LCfraction = cumsum(Trk_Thr.2.2pos_LCfraction))%>%
    mutate(Trk_Thr.test.2pos_LCfraction  = if_else(Trk_Thr.2.2.2pos_LCfraction >= 1, "piPOS", "piNEG" ))%>%
    ungroup()%>%
    group_by(Well_coordinate,
             Field,
             TrackID,
             timestep)%>%
    mutate(Trk_Thr.test.2.1pos_LCfraction = if_else( timestep == max(timestep) && is.na(Trk_Thr.test.2pos_LCfraction), "piPOS", paste(Trk_Thr.test.2pos_LCfraction) ))
  
  
  perWell.trk.pi.df <- perWell.trk.pi.df %>%
    mutate(Trk_Thrs.2pos_LCfraction = Trk_Thr.test.2.1pos_LCfraction )%>%
    ungroup()%>%
    select(
      -Trk_Thr.2.2pos_LCfraction ,
      -Trk_Thr.2.2.2pos_LCfraction,
      -Trk_Thr.test.2pos_LCfraction,
      -Trk_Thr.test.2.1pos_LCfraction)
  
  
  # Section 3.3.3 : Ilastik 1 pos   ----
  
  #----- Ilastik 1: Once positive always postive
  perWell.trk.pi.df <- perWell.trk.pi.df%>%
    dplyr::ungroup()%>%
    group_by(Filename,
             Well_coordinate,
             Field,
             TrackID)%>%
    mutate(Def.1.POC.Predicted.Class = if_else(POC.Predicted.Class == "piPOS", 1, 0))%>%
    mutate(Def.1.POC.Predicted.Class = cumsum(Def.1.POC.Predicted.Class))%>%
    mutate(Def.1.POC.Predicted.Class = if_else(Def.1.POC.Predicted.Class >= 1, "piPOS", "piNEG"))
  
  # Section 3.3.4 : Ilastik 2 pos   ----
  
  #----- Ilastik 2: Twice positive always postive
  perWell.trk.pi.df <- perWell.trk.pi.df %>%
    ungroup()%>%
    group_by(ExpFile,
             Well_coordinate,
             Field,
             TrackID)%>%
    mutate(Def.2.POC.Predicted.Class = if_else(POC.Predicted.Class == "piPOS", 1, 0))%>%
    mutate(Def.2.2.POC.Predicted.Class = if_else(lead(Def.2.POC.Predicted.Class) == 1 & Def.2.POC.Predicted.Class == 1,1,0))%>%
    mutate(Def.2.2.2.POC.Predicted.Class = cumsum(Def.2.2.POC.Predicted.Class))%>%
    mutate(Def.2.test.POC.Predicted.Class = if_else( Def.2.2.2.POC.Predicted.Class >= 1, "piPOS", "piNEG" ))%>%
    ungroup()%>%
    group_by(Well_coordinate,
             Field,
             TrackID,
             timestep)%>%
    mutate(Def.2.1.test.POC.Predicted.Class = if_else( timestep == max(timestep) && is.na(Def.2.test.POC.Predicted.Class), "piPOS", paste(Def.2.test.POC.Predicted.Class) ))
  
  
  
  perWell.trk.pi.df <- perWell.trk.pi.df %>%
    mutate(Def.2.POC.Predicted.Class = Def.2.1.test.POC.Predicted.Class )%>%
    ungroup()%>%
    select(
      -Def.2.2.POC.Predicted.Class,
      -Def.2.2.2.POC.Predicted.Class,
      -Def.2.test.POC.Predicted.Class,
      -Def.2.1.test.POC.Predicted.Class)
  
  
  # Section 3.4 : Calculating LCF and numbers across all tracking definitions  ----
  
  # ----- Calculating live cell fractions 
  
  # Section 3.4.1 : Thresholding 1 pos  ----
  
  #---- Thresholding 1 positive Calculate the live cell fraction across image frames
  Tracking_Thr.1pos <- perWell.trk.pi.df %>%
    mutate(Exp = Filename)%>%
    dplyr::select(
      Exp,
      TrackID,
      Well_coordinate,
      timestep,
      Trk_Thrs.1pos_LCfraction)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(
      Exp,
      Well_coordinate,
      timestep,
      Trk_Thrs.1pos_LCfraction) %>%
    summarise(n = n()) %>%
    dplyr::mutate(Trk_Thrs.1pos_LCF = n/sum(n)) %>%
    dplyr::filter( Trk_Thrs.1pos_LCfraction == "piNEG")%>%
    dplyr::mutate(Trk_Thrs.1pos_piNEG.n = n)%>%
    dplyr::select(-n)%>%
    ungroup()
  

  #----- PiNeg and Pos Numbers
  Tracking_Thr.1pos.N <-   perWell.trk.pi.df%>%
    mutate(Exp = Filename)%>%
    dplyr::select(
      Exp,
      TrackID,
      Well_coordinate,
      timestep,
      Trk_Thrs.1pos_LCfraction)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(
      Exp,
      Well_coordinate,
      timestep,
      Trk_Thrs.1pos_LCfraction) %>%
    summarise(n = n()) %>%
    dplyr::mutate(Trk_Thrs.1pos_LCF = n/sum(n)) %>%
    dplyr::mutate(Trk_Thrs.1pos.N = n)%>%
    dplyr::select(-n)%>%
    ungroup()
  
  #Incase we ever hit 0 live cell fractions 
  piNEG <- perWell.trk.pi.df %>%
    ungroup()%>%
    mutate(Exp = Filename)%>%
    select(timestep,
           Exp,
           Well_coordinate)%>%
    distinct()
  
  Tracking_Thr.1pos <- left_join(piNEG,
                                 Tracking_Thr.1pos)
  
  Tracking_Thr.1pos <-  Tracking_Thr.1pos %>%
    mutate(Trk_Thrs.1pos= "piNEG")%>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
  Tracking_Thr.1pos <-   Tracking_Thr.1pos %>%
    mutate(Trk_Thrs.1pos_LCfraction = Trk_Thrs.1pos_LCF)%>%
    select(-Trk_Thrs.1pos_LCF,
           -Trk_Thrs.1pos_piNEG.n)
  
  
  # Section 3.4.2 : Thresholding 2 pos  ----
  
  #---- Thresholding 2 positive Calculate the live cell fraction across image frames
  Tracking_Thr.2pos <- perWell.trk.pi.df %>%
    mutate(Exp = Filename)%>%
    dplyr::select(
      Exp,
      TrackID,
      Well_coordinate,
      timestep,
      Trk_Thrs.2pos_LCfraction)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(
      Exp,
      Well_coordinate,
      timestep,
      Trk_Thrs.2pos_LCfraction) %>%
    summarise(n = n()) %>%
    dplyr::mutate(Trk_Thrs.2pos_LCF = n/sum(n)) %>%
    dplyr::filter( Trk_Thrs.2pos_LCfraction == "piNEG")%>%
    dplyr::mutate(Trk_Thrs.2pos_piNEG.n = n)%>%
    dplyr::select(-n)%>%
    ungroup()
  
  
  #----- PiNeg and Pos Numbers
  Tracking_Thr.2pos.N <- perWell.trk.pi.df %>%
    mutate(Exp = Filename)%>%
    dplyr::select(
      Exp,
      TrackID,
      Well_coordinate,
      timestep,
      Trk_Thrs.2pos_LCfraction)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(
      Exp,
      Well_coordinate,
      timestep,
      Trk_Thrs.2pos_LCfraction) %>%
    summarise(n = n()) %>%
    dplyr::mutate(Trk_Thrs.2pos_LCF = n/sum(n)) %>%
    dplyr::mutate(Trk_Thrs.2pos.N = n)%>%
    dplyr::select(-n)%>%
    ungroup()
  
  #Incase we ever hit 0 live cell fractions 
  piNEG <- perWell.trk.pi.df %>%
    ungroup()%>%
    mutate(Exp = Filename)%>%
    select(timestep,
           Exp,
           Well_coordinate)%>%
    distinct()
  
  Tracking_Thr.2pos <- left_join(piNEG,
                                 Tracking_Thr.2pos)
  
  Tracking_Thr.2pos <-  Tracking_Thr.2pos %>%
    mutate(Trk_Thrs.2pos= "piNEG")%>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
  Tracking_Thr.2pos <-   Tracking_Thr.2pos %>%
    mutate(Trk_Thrs.2pos_LCfraction = Trk_Thrs.2pos_LCF)%>%
    select(-Trk_Thrs.2pos_LCF,
           -Trk_Thrs.2pos_piNEG.n)
  
  
  # Section 3.4.3 : Ilastik 1 pos  ----
  
  # ---- Calculate Ilastik live cell fraction dentition 1
  Tracking_pi.def1 <- perWell.trk.pi.df %>%
    mutate(Exp = Filename)%>%
    dplyr::select(
      Exp,
      TrackID,
      Well_coordinate,
      timestep,
      Def.1.POC.Predicted.Class)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(
      Exp,
      Well_coordinate,
      timestep,
      Def.1.POC.Predicted.Class) %>%
    summarise(n = n()) %>%
    dplyr::mutate(Trk_Ila.Def.1_LCfraction = n/sum(n)) %>%
    dplyr::filter(Def.1.POC.Predicted.Class == "piNEG")%>%
    dplyr::mutate(Trk_Ila.Def.1_LC_piNEG.n = n)%>%
    dplyr::select(-n)%>%
    ungroup()
  
  
  #----- PiNeg and Pos Numbers
  Tracking_pi.def1.N <- perWell.trk.pi.df %>%
    mutate(Exp = Filename)%>%
    dplyr::select(
      Exp,
      TrackID,
      Well_coordinate,
      timestep,
      Def.1.POC.Predicted.Class)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(
      Exp,
      Well_coordinate,
      timestep,
      Def.1.POC.Predicted.Class) %>%
    summarise(n = n()) %>%
    dplyr::mutate(Trk_Ila.Def.1_LCfraction = n/sum(n)) %>%
    dplyr::mutate(Trk_Ila.Def.1_LC.N = n)%>%
    dplyr::select(-n)%>%
    ungroup()
  
  
  Tracking_pi.def1 <- left_join(piNEG,
                                Tracking_pi.def1)
  
  Tracking_pi.def1 <- Tracking_pi.def1 %>%
    mutate(Def.1.POC.Predicted.Class = "piNEG")%>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
  
  
  
  Tracking_pi.def1 <- Tracking_pi.def1 %>%
    select(-Trk_Ila.Def.1_LC_piNEG.n,
           -Def.1.POC.Predicted.Class)
  
  # Section 3.4.5 : Ilastik 2 pos  ----
  
  # ---- Calculate Ilastik live cell fraction defition 2
  Tracking_pi.def2 <- perWell.trk.pi.df %>%
    mutate(Exp = Filename)%>%
    dplyr::select(
      Exp,
      TrackID,
      Well_coordinate,
      timestep,
      Def.2.POC.Predicted.Class)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(
      Exp,
      Well_coordinate,
      timestep,
      Def.2.POC.Predicted.Class) %>%
    summarise(n = n()) %>%
    dplyr::mutate(Trk_Ila.Def.2_LCfraction = n/sum(n)) %>%
    dplyr::filter(Def.2.POC.Predicted.Class == "piNEG")%>%
    dplyr::mutate(Trk_Ila.Def.2_LC_piNEG.n = n)%>%
    dplyr::select(-n)%>%
    ungroup()
  
  
  #--- Pi POS Neg numbers
  Tracking_pi.def2.N <- perWell.trk.pi.df %>%
    mutate(Exp = Filename)%>%
    dplyr::select(
      Exp,
      TrackID,
      Well_coordinate,
      timestep,
      Def.2.POC.Predicted.Class)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(
      Exp,
      Well_coordinate,
      timestep,
      Def.2.POC.Predicted.Class) %>%
    summarise(n = n()) %>%
    dplyr::mutate(Trk_Ila.Def.2_LCfraction = n/sum(n)) %>%
    dplyr::mutate(Trk_Ila.Def.2_LC_piNEG.N = n)%>%
    dplyr::select(-n)%>%
    ungroup()
  
  
  Tracking_pi.def2 <- left_join(piNEG,
                                Tracking_pi.def2)
  
  Tracking_pi.def2 <- Tracking_pi.def2 %>%
    mutate(Def.2.POC.Predicted.Class = "piNEG")%>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
  
  Tracking_pi.def2 <- Tracking_pi.def2 %>%
    select(-Trk_Ila.Def.2_LC_piNEG.n,
           -Def.2.POC.Predicted.Class)
  
  
  # Section 4 : Pooling tracking results  ----
  
  #-- Single Live cell fraction tracking dataframe
  
  
  Tracking_pi <- left_join( Tracking_pi.def1,
                           Tracking_pi.def2)
  
  Tracking_pi <- left_join(Tracking_pi,
                           Tracking_Thr.1pos)
  
  Tracking_pi <- left_join(Tracking_pi,
                           Tracking_Thr.2pos)
  
  rm(perWell.trk.pi.df,
     Tracking_pi.def1,
     Tracking_pi.def2,
     Tracking_Thr.1pos,
     Tracking_Thr.2pos)
  
  
  #Dropping irrelevant colunms
  Tracking_pi <- Tracking_pi %>%
    select(Exp,
           Well_coordinate,
           timestep,
           Trk_Ila.Def.1_LCfraction,
           Trk_Ila.Def.2_LCfraction,
           Trk_Thrs.1pos_LCfraction,
           Trk_Thrs.2pos_LCfraction)
  


  Tracking_pi.melted <- melt(Tracking_pi, 
                             id.vars =  c("Exp",
                                          "Well_coordinate",
                                          "timestep"),
                             variable.name = "Killing.Def")
  
  Tracking_pi.melted <-   Tracking_pi.melted %>%
    mutate(Exp = gsub(paste("-",Tracking.version,sep=""),"",Exp))%>%
    mutate(Killing.Def = gsub("Trk_",paste(Tracking.version,"_",sep=""),Killing.Def),
           BaSic.version = if_else(BaSic.parameter == "poc", "",
                                   if_else(BaSic.parameter == "pocbSd", "BaSic","UNKWN")))
  
  Tracking_pi.melted <-   Tracking_pi.melted %>%
    mutate(Killing.Def = if_else( BaSic.version == "BaSic", paste(Killing.Def,
                                                                  BaSic.version,
                                                                  sep="."), Killing.Def))%>%
    select(-BaSic.version)
  
  Tracking_pi <- spread(Tracking_pi.melted, key = "Killing.Def",
                        value = "value")
  
  
  
  rm(Tracking_pi.melted)
  
  if ( progress == 1 ) {
    MASTER.df.Tracking_pi <- Tracking_pi }
  
  if ( progress  > 1 ) {
    
    MASTER.df.Tracking_pi <- left_join(    MASTER.df.Tracking_pi,
                                           Tracking_pi)
    
  }
  
  
  # Section 5 : Export LCF results  ----
  
  # --- Export loop
  filename <- unique(Tracking_pi$Exp)
  
  filename <- paste(filename,
                    "_",
                    Tracking.version,
                    "_LC.csv",
                    sep="")
  
  filename <- paste(tracking.res.export,
                    "/",
                    filename,sep="")
  
  
  
  
  # Section 6 : Cell numbers across viability states  ----
  
  #------- saving tracking live.dead Numbers
  
  # corecting instance wheren cell numbers are low
  timestep.length <- max(Tracking_pi$timestep)
  timestep.length <- 0:timestep.length
  timestep.length <- rep(timestep.length, each =2)
  
  Tracking_Trk.N.corr <- data.frame( timestep =  as.numeric( timestep.length))
  
  
  
  #--- Ila.1pos
  
  Tracking_pi.def1.N.corr <- data.frame( timestep =  as.numeric( timestep.length))
  
  Tracking_pi.def1.N.corr <-   Tracking_pi.def1.N.corr  %>%
    mutate(Exp = unique(Tracking_pi$Exp),
           Well_coordinate =  unique(Tracking_pi$Well_coordinate),
           timestep =  as.numeric( timestep.length),
           timestep2 = 1,
           timestep2 = cumsum(timestep2))%>%
    mutate(Def.1.POC.Predicted.Class = if_else(timestep2 %% 2 == 1, "piNEG", "piPOS"))%>%
    select(-timestep2)
  
  
  Tracking_pi.def1.N.corr <- left_join(  Tracking_pi.def1.N.corr,
                                         Tracking_pi.def1.N  )
  
  Tracking_pi.def1.N <-   Tracking_pi.def1.N.corr%>%
    mutate(Trk_Ila.Def.1_LC.N = as.numeric(Trk_Ila.Def.1_LC.N))%>%
    group_by(timestep) %>%
    mutate(
      Trk_Ila.Def.1_LCfraction = if_else(is.na(Trk_Ila.Def.1_LCfraction) & Def.1.POC.Predicted.Class == "piNEG", 1 - Trk_Ila.Def.1_LCfraction[Def.1.POC.Predicted.Class == "piPOS"], Trk_Ila.Def.1_LCfraction),
      Trk_Ila.Def.1_LCfraction = if_else(is.na(Trk_Ila.Def.1_LCfraction) & Def.1.POC.Predicted.Class == "piPOS", 0, Trk_Ila.Def.1_LCfraction),
      Trk_Ila.Def.1_LC.N = ifelse( is.na( Trk_Ila.Def.1_LC.N), as.numeric(0), Trk_Ila.Def.1_LC.N)) %>%
    ungroup()
  
  
  
  #----- ila 2 pos
  Tracking_pi.def2.N.corr <- data.frame( timestep =  as.numeric( timestep.length))
  
  Tracking_pi.def2.N.corr <-   Tracking_pi.def2.N.corr  %>%
    mutate(Exp = unique(Tracking_pi$Exp),
           Well_coordinate =  unique(Tracking_pi$Well_coordinate),
           timestep =  as.numeric( timestep.length),
           timestep2 = 1,
           timestep2 = cumsum(timestep2))%>%
    mutate(Def.2.POC.Predicted.Class = if_else(timestep2 %% 2 == 1, "piNEG", "piPOS"))%>%
    select(-timestep2)
  
  
  Tracking_pi.def2.N.corr <- left_join(  Tracking_pi.def2.N.corr,
                                         Tracking_pi.def2.N  )
  
  Tracking_pi.def2.N <-   Tracking_pi.def2.N.corr%>%
    mutate(Trk_Ila.Def.2_LC_piNEG.N = as.numeric(Trk_Ila.Def.2_LC_piNEG.N))%>%
    group_by(timestep) %>%
    mutate(
      Trk_Ila.Def.2_LCfraction = if_else(is.na(Trk_Ila.Def.2_LCfraction) & Def.2.POC.Predicted.Class == "piNEG", 1 - Trk_Ila.Def.2_LCfraction[Def.2.POC.Predicted.Class == "piPOS"], Trk_Ila.Def.2_LCfraction),
      Trk_Ila.Def.2_LCfraction = if_else(is.na(Trk_Ila.Def.2_LCfraction) & Def.2.POC.Predicted.Class == "piPOS", 0, Trk_Ila.Def.2_LCfraction),
      Trk_Ila.Def.2_LC.N = ifelse( is.na( Trk_Ila.Def.2_LC_piNEG.N), as.numeric(0), Trk_Ila.Def.2_LC_piNEG.N)) %>%
    ungroup()%>%
    select(-Trk_Ila.Def.2_LC_piNEG.N)
  
  
  #---- Thr 1pos N 
  Tracking_Thr.1pos.N.corr <- data.frame( timestep =  as.numeric( timestep.length))
  
  Tracking_Thr.1pos.N.corr <-    Tracking_Thr.1pos.N.corr %>%
    mutate(Exp = unique(Tracking_pi$Exp),
           Well_coordinate =  unique(Tracking_pi$Well_coordinate),
           timestep =  as.numeric( timestep.length),
           timestep2 = 1,
           timestep2 = cumsum(timestep2))%>%
    mutate(Trk_Thrs.1pos_LCfraction = if_else(timestep2 %% 2 == 1, "piNEG", "piPOS"))%>%
    select(-timestep2)
  
  
  Tracking_Thr.1pos.N.corr<- left_join(Tracking_Thr.1pos.N.corr,
                                       Tracking_Thr.1pos.N )
  
  Tracking_Thr.1pos.N<-  Tracking_Thr.1pos.N.corr%>%
    mutate(Trk_Thrs.1pos.N = as.numeric(Trk_Thrs.1pos.N))%>%
    group_by(timestep) %>%
    mutate(
      Trk_Thrs.1pos_LCF = if_else(is.na(Trk_Thrs.1pos_LCF) & Trk_Thrs.1pos_LCfraction == "piNEG", 1 - Trk_Thrs.1pos_LCF[Trk_Thrs.1pos_LCfraction  == "piPOS"], Trk_Thrs.1pos_LCF),
      Trk_Thrs.1pos_LCF = if_else(is.na(Trk_Thrs.1pos_LCF) & Trk_Thrs.1pos_LCfraction  == "piPOS", 0, Trk_Thrs.1pos_LCF),
      Trk_Thrs.1pos.N = ifelse( is.na( Trk_Thrs.1pos.N), as.numeric(0), Trk_Thrs.1pos.N)) %>%
    ungroup()
  
  
  
  #---- Thr 2pos N 
  Tracking_Thr.2pos.N.corr <- data.frame( timestep =  as.numeric( timestep.length))
  
  Tracking_Thr.2pos.N.corr <-    Tracking_Thr.2pos.N.corr %>%
    mutate(Exp = unique(Tracking_pi$Exp),
           Well_coordinate =  unique(Tracking_pi$Well_coordinate),
           timestep =  as.numeric( timestep.length),
           timestep2 = 1,
           timestep2 = cumsum(timestep2))%>%
    mutate(Trk_Thrs.2pos_LCfraction = if_else(timestep2 %% 2 == 1, "piNEG", "piPOS"))%>%
    select(-timestep2)
  
  
  Tracking_Thr.2pos.N.corr<- left_join(Tracking_Thr.2pos.N.corr,
                                       Tracking_Thr.2pos.N )
  
  Tracking_Thr.2pos.N <-  Tracking_Thr.2pos.N.corr%>%
    mutate(Trk_Thrs.2pos.N = as.numeric(Trk_Thrs.2pos.N))%>%
    
    group_by(timestep) %>%
    mutate(
      Trk_Thrs.2pos_LCF = if_else(is.na(Trk_Thrs.2pos_LCF) & Trk_Thrs.2pos_LCfraction == "piNEG", 1 - Trk_Thrs.2pos_LCF[Trk_Thrs.2pos_LCfraction  == "piPOS"], Trk_Thrs.2pos_LCF),
      Trk_Thrs.2pos_LCF = if_else(is.na(Trk_Thrs.2pos_LCF) & Trk_Thrs.2pos_LCfraction  == "piPOS", 0, Trk_Thrs.2pos_LCF),
      Trk_Thrs.2pos.N = if_else( is.na( Trk_Thrs.2pos.N), as.numeric(0), Trk_Thrs.2pos.N)) %>%
    ungroup()
  
  # Houskeeping
  
  
  rm(Tracking_Trk.N.corr,
     Tracking_pi.def1.N.corr,
     Tracking_pi.def2.N.corr,
     Tracking_Thr.1pos.N.corr,
     Tracking_Thr.2pos.N.corr)
 
  #-- Ila 1 and 2
  Tracking_pi.def1.N <-  Tracking_pi.def1.N %>%
    mutate(live.dead.def = Def.1.POC.Predicted.Class) %>%
    select(-Def.1.POC.Predicted.Class)
  
  
  Tracking_pi.def2.N <-  Tracking_pi.def2.N %>%
    mutate(live.dead.def = Def.2.POC.Predicted.Class) %>%
    select(-Def.2.POC.Predicted.Class)
  
  #-- Thrs 1 and 2 
  Tracking_Thr.1pos.N<- Tracking_Thr.1pos.N %>%
    mutate(live.dead.def = Trk_Thrs.1pos_LCfraction) %>%
    select(-Trk_Thrs.1pos_LCfraction)
  
  Tracking_Thr.2pos.N<- Tracking_Thr.2pos.N %>%
    mutate(live.dead.def = Trk_Thrs.2pos_LCfraction) %>%
    select(-Trk_Thrs.2pos_LCfraction)
  
  
 
  
  Tracking.N <- left_join(   Tracking_pi.def1.N ,#Ila1
                            Tracking_pi.def2.N)# Ila2
  
  
  Tracking.N <- left_join(Tracking.N ,
                          Tracking_Thr.1pos.N)# Thr.1
  
  Tracking.N <- left_join(Tracking.N ,
                          Tracking_Thr.2pos.N)# Thr.2
  
  
  #Housekeeping
  rm(#Tracking_Trk.N,
     Tracking_pi.def1.N,
     Tracking_pi.def2.N,
     Tracking_Thr.1pos.N,
     Tracking_Thr.2pos.N)
  
  
  Tracking.N <- Tracking.N%>%
    distinct()
  
  
  
  Tracking.N.melted <- melt(Tracking.N, 
                            id.vars =  c("timestep",
                                         "Exp",
                                         "Well_coordinate",
                                      #   "Trk_LCfraction",
                                         #"piDyn.N",
                                         "live.dead.def"),
                            variable.name = "Killing.Def")
  
  Tracking.N.melted <-   Tracking.N.melted %>%
    mutate(Killing.Def = gsub("Trk_",paste(Tracking.version,"_",sep=""),Killing.Def),
           BaSic.version = if_else(BaSic.parameter == "poc", "",
                                   if_else(BaSic.parameter == "pocbSd", "BaSic","UNKWN")))
  
  Tracking.N.melted <-   Tracking.N.melted %>%
    mutate(Killing.Def = if_else( BaSic.version == "BaSic", paste(Killing.Def,
                                                                  BaSic.version,
                                                                  sep="."), Killing.Def)
           
    )%>%
    select(-BaSic.version)
  
  Tracking.N <- spread(Tracking.N.melted, key = "Killing.Def",
                       value = "value")
  
  
  rm(Tracking.N.melted)
  
  
  if ( progress == 1 ) {
    MASTER.df.Tracking.N <- Tracking.N }
  
  if (  progress  > 1 ) {
    
    MASTER.df.Tracking.N <- left_join(        MASTER.df.Tracking.N,
                                              Tracking.N)
    
  }
  
  
  
  
  filename2 <- unique(Tracking.N$Exp)
  
  
  
  filename2 <- paste(filename2,
                     "_",
                     Tracking.version,
                     "_LC.N.csv",
                     sep="")
  
  filename2 <- paste(tracking.res.export,
                     "/",
                     filename2,sep="")
  
  #Houskeeping
  rm(piNEG,
     Tracking_pi,
     Tracking.N)
  
}



write.csv(MASTER.df.Tracking_pi,
          filename,
          row.names = FALSE)



write.csv(MASTER.df.Tracking.N,
          filename2,
          row.names = FALSE)



# Statment of file completed
filename <- unique(MASTER.df.Tracking_pi$Exp)
filename <- paste(filename,
                  " tracking  live-cell fraction is now COMPLETE!",
                  sep="")


#}