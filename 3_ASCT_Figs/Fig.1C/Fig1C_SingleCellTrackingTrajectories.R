
# Aim ------------------------


# This script will go in to the tracking labeloid directory which is Ahmads tracking results and merge the data against the Basic FF9DF0 (old background correction) and new BaSic default background correction.
# It will then find the LCF across these definitions.
# The code is however rigid. Essentially I need to duplicate alot of lines to accomodate the results of BaSic correction. Unless I duplicate this R script and

# Section 1: defining paths and directories  ------------------------

#General variables
genDir <- getwd() 


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
library(gghighlight)
library(reshape2)
library(ggprism)


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

list.of.POC.result.dir <- list.of.POC.result.dir[grep("OCsimple.pocbSd", list.of.POC.result.dir)]
i <- list.of.POC.result.dir[1]

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
                    into = c(
                      "ExpFile",
                      "Well_coordinate",
                      "Field"),
                    sep="_")
  
  OC.df <- OC.df %>%
    dplyr::mutate(POC.Predicted.Class = Predicted.Class)%>%
    dplyr::select(
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
  
  # Exporting Cascading Colourmaps plot
  
setwd(genDir)
  #Housekeeping
  rm(OC.df)
  
  # Section 3.3 : Defining live cell definition across all tracking definitions  ----
  
  
  # Calculate the live cell fraction  across all tracking definitions 
  # Thresholding: 1) >= 700 mean intesity & once positive always positive  2) >=700 mean intensity & twice psisive always positive
  

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
  
 if( Tracking.version == "Trkv2") {
   
   
   # ---- Calculate Ilastik live cell fraction defition 2
   Tracking_pi.def2 <- perWell.trk.pi.df %>%
     dplyr::select(
       ExpFile,
       TrackID,
       Well_coordinate,
       timestep,
       Def.2.POC.Predicted.Class)%>%
     dplyr::ungroup()%>%
     dplyr::group_by(
       ExpFile,
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
     dplyr::select(
       ExpFile,
       TrackID,
       Well_coordinate,
       timestep,
       Def.2.POC.Predicted.Class)%>%
     dplyr::ungroup()%>%
     dplyr::group_by(
       ExpFile,
       Well_coordinate,
       timestep,
       Def.2.POC.Predicted.Class) %>%
     summarise(n = n()) %>%
     dplyr::mutate(Trk_Ila.Def.2_LCfraction = n/sum(n)) %>%
     dplyr::mutate(Trk_Ila.Def.2_LC_piNEG.N = n)%>%
     dplyr::select(-n)%>%
     ungroup()
   
 
   
   perWell.trk.pi.df <- perWell.trk.pi.df %>%
     mutate(Trkv2_Ila2BaSic_POCpredictedClass = Def.2.POC.Predicted.Class)%>%
     select(ExpFile,
            Well_coordinate,
            Field,
            TrackID,
            timestep,
            Trkv2_Ila2BaSic_POCpredictedClass,
            Ila.Mean.Intensity)
   
 
 }
 
  
  
}




Raw.TKC.df <- Tracking_pi.def2 %>%
  mutate(timestep = timestep + 1)%>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         LCF = Trk_Ila.Def.2_LCfraction) # Live cell fraction

Raw.TKC.df.0 <- data.frame(ExpFile = unique(Raw.TKC.df$ExpFile),
                           Well_coordinate = unique(Raw.TKC.df$Well_coordinate),
                           timestep = 0,
                           LCF = 1)
Raw.TKC.df <- rbind(Raw.TKC.df,
                    Raw.TKC.df.0 )
rm(Tracking_pi.def2.N,
   Raw.TKC.df.0,
   Tracking_pi.def2)


# loading time metadata information from populational analysis output
md.path <- paste(genDir,
                 "/",
                 unique(perWell.trk.pi.df$ExpFile),
                 "_",
                 unique(perWell.trk.pi.df$Well_coordinate),
                 "_Results/",
                 sep="")

Time.df <- list.files(path = md.path,
                             pattern = "*_results.csv",
                             full.names = TRUE)

Time.df <- read.csv(Time.df)

Time.df <- Time.df %>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         Time_Hrs)%>%
  distinct()


perWell.trk.pi.df <- perWell.trk.pi.df %>%
  left_join(Time.df)

Time.Raw.TKC.df <- Time.df %>%
  mutate(timestep = timestep + 1)

Time.Raw.TKC.df.0 <- data.frame(ExpFile = unique(Raw.TKC.df$ExpFile),
             Well_coordinate = unique(Raw.TKC.df$Well_coordinate),
             timestep = 0,
             Time_Hrs = 0)

Time.Raw.TKC.df <- rbind(Time.Raw.TKC.df,
                         Time.Raw.TKC.df.0)

Raw.TKC.df <- left_join(Raw.TKC.df,
                        Time.Raw.TKC.df,
                        by = c("ExpFile",
                               "Well_coordinate",
                               "timestep"))

rm(Time.df,
   Time.Raw.TKC.df,
   Time.Raw.TKC.df.0)



# ------- Rest of script -----


ID.df <-  perWell.trk.pi.df %>%
  ungroup()%>%
  mutate(Exp_Well_Field_uniTrackID = paste(ExpFile,
                                           Well_coordinate,
                                           Field,
                                           TrackID, sep="_"))%>%
  select(Exp_Well_Field_uniTrackID)%>%
  distinct()%>%
  mutate(UniqueTrackID = 1,
         UniqueTrackID = cumsum(UniqueTrackID))
           
           
       
perWell.trk.pi.df <- perWell.trk.pi.df %>%
  ungroup()%>%
  mutate(Exp_Well_Field_uniTrackID = paste(ExpFile,
                                           Well_coordinate,
                                           Field,
                                           TrackID, sep="_"))%>%
  left_join(ID.df )
  




# Rank the TrackIDs based on when they turned positive

# piNEG to piPOS
# piNEG only
# piPOS only

library(dplyr)

perWell.trk.pi.df <- perWell.trk.pi.df

# Convert Trkv2_Ila2BaSic_POCpredictedClass to a factor for easier comparison
perWell.trk.pi.df$Trkv2_Ila2BaSic_POCpredictedClass <- as.factor(perWell.trk.pi.df $Trkv2_Ila2BaSic_POCpredictedClass)



# Finding the number of piPos flags for each TrackID;
perWell.trk.pi.df <- perWell.trk.pi.df  %>%
  group_by(UniqueTrackID) %>%
  arrange(timestep)%>%
  mutate(Total.piPos.labels = sum(Trkv2_Ila2BaSic_POCpredictedClass == "piPOS"))
  

# Ranking trackID, by ordering by two variables, the total number of piPOS labels and the UniqueTrackID. Next we define a ranking definitons
perWell.trk.pi.df.RANK <- perWell.trk.pi.df %>%
  select(UniqueTrackID,
         Total.piPos.labels) %>%
  distinct()%>%
  ungroup()%>%
  arrange((Total.piPos.labels),
          UniqueTrackID) %>%
  ungroup()%>%
  mutate(Rank = 1)%>%
  mutate(Rank = cumsum(Rank))

perWell.trk.pi.df <- perWell.trk.pi.df %>%
  left_join(perWell.trk.pi.df.RANK )

rm(perWell.trk.pi.df.RANK)
plot.path <- paste(genDir,
                   "/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_Results/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_",
                   "TrackIntPlot_Ranked",
                   ".pdf", sep = "")


pdf(plot.path)

perWell.trk.pi.df  %>%
  ggplot(
    aes(x = Time_Hrs,
        y = Rank,
        colour = Ila.Mean.Intensity,
        group = UniqueTrackID
    )
  )+
  geom_line()+
  scale_color_gradient(low = "black",
                       
                       high = "red") +
  labs(x = "Time (Hrs)", 
       y = "Ranked trackIDs",
       color = "Pi Intensity")+
  scale_x_continuous(breaks = c(0, 24, 48, 72))+
  scale_y_continuous(breaks = seq(0, max(perWell.trk.pi.df$Rank), by = 1000))+
  theme_prism()
 

dev.off()


# ---- Log transformed Ranked plot
# NB: many PI values are 0 (piNeg); when you log10 transform 0 you get -INF. therefore for all 0 values I change the result to 1.
plot.path <- paste(genDir,
                   "/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_Results/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_",
                   "TrackIntPlot_Ranked_LOG",
                   ".pdf", sep = "")

pdf(plot.path)

perWell.trk.pi.df  %>%
  mutate(Ila.Mean.Intensity = if_else(Ila.Mean.Intensity <= 1, 1, Ila.Mean.Intensity))%>%
  
  ggplot(
    aes(x = Time_Hrs,
        y = Rank,
        colour = log10(Ila.Mean.Intensity),
        group = Rank
    )
  )+
  geom_line()+
  scale_color_gradient(low = "black",
                       high = "red") +
  
  labs(x = "Time (Hrs)", 
       y = "Ranked trackIDs",
       color = "Pi Intensity")+
  scale_x_continuous(breaks = c(0,
                                24,
                                48,
                                72))+
  scale_y_continuous(breaks = seq(0,
                                  max(perWell.trk.pi.df$Rank),
                                  by = 1000))+
  theme_prism()


dev.off()

# ---- Load Tracking Time kill curves
file.name <- paste(genDir,
                   "/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_Results/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   "Fitnorm",
                   "_Trkv2_Flags.csv",
                   sep="")

Time.kill.curve <- read.csv(file.name )

Time.kill.curve <- Time.kill.curve%>%
  filter(Well_coordinate == unique(perWell.trk.pi.df$Well_coordinate))%>%
  filter(Time.Kill.Definitions == "Trkv2_Ila2")



plot.path <- paste(genDir,
                   "/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_Results/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_",
                   "TrackIntPlot_Ranked_LOG-dualAxis",
                   ".pdf", sep = "")
Norm.factor <- max(perWell.trk.pi.df$Rank)

pdf(plot.path)

perWell.trk.pi.df  %>%
  mutate(Ila.Mean.Intensity = if_else(Ila.Mean.Intensity <= 1, 1, Ila.Mean.Intensity))%>%
  mutate(Ila.Mean.Intensity = log10(Ila.Mean.Intensity) )%>%
  ggplot(
    aes(x = Time_Hrs,
        y = Rank,
        colour = Ila.Mean.Intensity,
        group = Rank))+
  geom_line()+
  geom_line(data = Time.kill.curve,
            aes( x = Time_Hrs,
                 y = LC.fraction.corr*Norm.factor,
                 group =Well_coordinate), 
            color = "blue")+ # Second axis
  scale_color_gradient( name = "Log mean Pi intensity",

                       low = "black",
                       
                       high = "red") +
  labs(color = "Log mean Pi intensity")+

  labs(x = "Time (Hrs)", 
       y = "TrackIDs",
       caption =  paste(unique(Time.kill.curve$Abx.con),
                         " ",
                         unique(Time.kill.curve$Isolate),
                         " well ",
                         unique(perWell.trk.pi.df$Well_coordinate),
                         ": ",
                         
                         Norm.factor,
                        " cells",
                        
                         sep= ""))+
  scale_x_continuous(breaks = c(0,
                                24,
                                48,
                                72))+
  scale_y_continuous(breaks = seq(0,
                                  Norm.factor , by = 1000),
                     sec.axis = sec_axis(~./Norm.factor*100, name = "Live cell fraction [au]", breaks = seq(0,  1, by = 0.1)))+
  theme_prism()+
  theme( legend.position = "right",
        axis.line.y.right = element_line(color = "blue"), 
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"), 
        axis.title.y.right = element_text(color = "blue"))



dev.off()




# ---- Adjusted colour scheme: Diverging colour scheme


# 
plot.path <- paste(genDir,
                   "/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_Results/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_",
                   "TrackIntPlot_Ranked_LOG-dualAxis-AdjColour",
                   ".pdf", sep = "")
Norm.factor <- max(perWell.trk.pi.df$Rank)

#---
CHECK2 <- perWell.trk.pi.df  %>%
  #filter(UniqueTrackID <= 500)%>%
  ungroup()
 

x <-ggplot(CHECK2,
       aes(x = Ila.Mean.Intensity))+
  geom_histogram(binwidth = 1)

library(plotly)

ggplotly(x)
#---

CHECK <- perWell.trk.pi.df  %>%
  ungroup()%>%
  mutate(Ila.Mean.Intensity = if_else(Ila.Mean.Intensity <= 1, 1, Ila.Mean.Intensity))%>%
  
  mutate(Ila.Mean.Intensity = log10(Ila.Mean.Intensity) )

x2 <-  ggplot(CHECK,
         aes(x = Ila.Mean.Intensity))+
  geom_histogram(binwidth = 0.01)
  
ggplotly(x2)
min(perWell.trk.pi.df$Ila.Mean.Intensity)
log10(1)

log10(0.008849558)
log10(100)

Check.max <- max(CHECK$Ila.Mean.Intensity)
mid <- 2

# ---- GPT 4

plot.path <- paste(genDir,
                   "/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_Results/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_",
                   "TrackIntPlot_Ranked_LOG-dualAxis-AdjColour-final",
                   ".pdf", sep = "")

pdf(plot.path)

perWell.trk.pi.df  %>%
  #filter(UniqueTrackID <= 500)%>%
  mutate(Ila.Mean.Intensity = if_else(Ila.Mean.Intensity <= 1, 1, Ila.Mean.Intensity)) %>%
  mutate(Ila.Mean.Intensity = log10(Ila.Mean.Intensity)) %>%
  ggplot(
    aes(
      x = Time_Hrs,
      y = Rank,
      colour = Ila.Mean.Intensity,
      group = Rank
    )) +
  geom_line() +
  geom_line(
    data = Time.kill.curve,
    aes(
      x = Time_Hrs,
      y = LC.fraction.corr*Norm.factor,
      group = Well_coordinate
    ),
    color = "steelblue4"
  ) + # Second axis
  
  
  scale_colour_gradientn(
    name = "Log mean Pi intensity",
    colors = c("black",
               "black",
               "red"),
    values = scales::rescale(c(0, 2.0001, 3.5)),
    limits = c(0, 3.5),
    space = "Lab"
  ) +
  
  labs(color = "Log mean Pi intensity") +
  labs(
    x = "Time (Hrs)",
    y = "TrackIDs",
    caption = paste(
      unique(Time.kill.curve$Abx.con),
      " ",
      unique(Time.kill.curve$Isolate),
      " well ",
      unique(perWell.trk.pi.df$Well_coordinate),
      ": ",
      Norm.factor,
      " cells",
      sep = ""
    )
  ) +
  scale_x_continuous(breaks = c(0,
                                12,
                                24,
                                36,
                                48,
                                60,
                                72)) +
  scale_y_continuous(
    breaks = seq(0,
                Norm.factor,
                 by = 1000),
    sec.axis = sec_axis(
      ~ .*100 /Norm.factor,
      name = "Live cells [%]",
      breaks = seq(0, 100, by = 10)
    )
  ) +
 
theme_pubr()+
  theme(
    legend.position = "right",
    axis.line.y.right = element_line(color = "steelblue4"),
    axis.ticks.y.right = element_line(color = "steelblue4"),
    axis.text.y.right = element_text(color = "steelblue4"),
    axis.title.y.right = element_text(color = "steelblue4")
  )

dev.off()



# 
plot.path <- paste(genDir,
                   "/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_Results/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_",
                   "TrackIntPlot_Ranked_LOG-dualAxis-AdjColour",
                   ".pdf", sep = "")
Norm.factor <- max(perWell.trk.pi.df$Rank)

#---
CHECK2 <- perWell.trk.pi.df  %>%
  ungroup()
 

x <-ggplot(CHECK2,
       aes(x = Ila.Mean.Intensity))+
  geom_histogram(binwidth = 1)

library(plotly)

ggplotly(x)
#---

CHECK <- perWell.trk.pi.df  %>%
  ungroup()%>%
  mutate(Ila.Mean.Intensity = if_else(Ila.Mean.Intensity <= 1, 1, Ila.Mean.Intensity))%>%
  
  mutate(Ila.Mean.Intensity = log10(Ila.Mean.Intensity) )

x2 <-  ggplot(CHECK,
         aes(x = Ila.Mean.Intensity))+
  geom_histogram(binwidth = 0.01)
  
ggplotly(x2)
min(perWell.trk.pi.df$Ila.Mean.Intensity)
log10(1)

log10(0.008849558)
log10(100)

Check.max <- max(CHECK$Ila.Mean.Intensity)
#mid <- mean(CHECK$Ila.Mean.Intensity)-1*sd(CHECK$Ila.Mean.Intensity)
mid <- 2


plot.path <- paste(genDir,
                   "/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_Results/",
                   unique(perWell.trk.pi.df$ExpFile),
                   "_",
                   unique(perWell.trk.pi.df$Well_coordinate),
                   "_",
                   "TrackIntPlot_Ranked_LOG-dualAxis-AdjColour-final-RAW",
                   ".pdf", sep = "")

pdf(plot.path)

perWell.trk.pi.df  %>%
  mutate(Ila.Mean.Intensity = if_else(Ila.Mean.Intensity <= 1, 1, Ila.Mean.Intensity)) %>%
  mutate(Ila.Mean.Intensity = log10(Ila.Mean.Intensity)) %>%
  ggplot(
    aes(
      x = Time_Hrs,
      y = Rank,
      colour = Ila.Mean.Intensity,
      group = Rank
    )) +
  geom_line() +
  
  
  geom_line(
    data = Raw.TKC.df,
    aes(
      x = Time_Hrs,
      y = LCF*Norm.factor,
      group = Well_coordinate
    ),
    color = "blue"
  ) + # Second axis "
  
  
  scale_colour_gradientn(
    name = "Log mean Pi intensity",
    colors = c("black",
               "black",
               "red"),
    values = scales::rescale(c(0, 2.0001, 3.5)),
    limits = c(0, 3.5),
    space = "Lab"
  ) +
  
  labs(color = "Log mean Pi intensity") +
  labs(
    x = "Time (Hrs)",
    y = "TrackIDs",
    caption = paste(
      unique(Time.kill.curve$Abx.con),
      " ",
      unique(Time.kill.curve$Isolate),
      " well ",
      unique(perWell.trk.pi.df$Well_coordinate),
      ": ",
      Norm.factor,
      " cells",
      sep = ""
    )
  ) +
  scale_x_continuous(breaks = c(0,
                                12,
                                24,
                                36,
                                48,
                                60,
                                72)) +
  scale_y_continuous(
    breaks = seq(0,
                Norm.factor,
                 by = 1000),
    sec.axis = sec_axis(
      ~ .*100 /Norm.factor,
      name = "Live cells [%]",
      breaks = seq(0, 100, by = 10)
    )
  ) +
  
theme_pubr()+
  theme(
    legend.position = "right",
    axis.line.y.right = element_line(color = "steelblue4"),
    axis.ticks.y.right = element_line(color = "steelblue4"),
    axis.text.y.right = element_text(color = "steelblue4"),
    axis.title.y.right = element_text(color = "steelblue4")
  )

dev.off()




