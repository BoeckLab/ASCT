

#---- Aim: 
# Considering tracking and population data assess different killing features e.g Area Under the Curve, MDK and live cell fraction at distinct timpoints

# KEY input variables


# frequency of imaging (e.g loop every 2h)
freq.loop <- (2)

# How long was this timelapse experiment? Was a 24, 48 ,72 or 96h timelapse
timelapse <- c(72)


# Based on the duration of the time-lapse, what real-time live cell fractions values are you interested in? 
# E.g. In the old TLCI experiments we had considered 3, 6, 9 12,14,48 and 96h interesting real time values. Relevant for KF 3 and KF5 (see chunks)
rt.values.interest <- c(3,6,9,12,24,36,48,60,72 )


# Ranges of AUC of interest e.g. from 0-3h or from 12-24h.
auc.of.interest <- c("0-3","0-6","0-9","0-24", "0-36", "0-48", "0-60", "0-72", "48-72")



# What type input data are we dealing with Population analysis only?  Or Both Population + tracking data: variable options are: Population or Population_&_Tracking
analysis.data <- c("Population_&_Tracking")
tracking.versions <- c("Trkv2")

#1) Settign up path variables


#General variables
genDir <- getwd() 

#res <- c("PerWell-Results") 
#tracking.res <- c("Results-LCdef")

exp.res <- list.files(path = ".", pattern = "_Results", all.files = FALSE)
trk.exp.res <- list.files(path = ".", pattern = "_Tracking_LC", all.files = FALSE)

#Directories
wdDir <- genDir

perWell.resDir <- paste(genDir,"/",exp.res,sep="")
exp.resDir <- paste(genDir,"/", exp.res,sep="")
trk.resDir <- paste(genDir,"/",trk.exp.res ,sep="")


#Housekeeping
rm(res)
#Setting work directory
setwd(wdDir)



#2) Loading libraries

library(stringr)
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




perWell.resDir

trk.resDir



# Loading population analysis data
loading.perWell.df <- function (perWell.resDir){
  
  
  
  #Creacting vector with the list of file names
  perWell.filenames <- list.files(path = perWell.resDir,
                                  pattern = "2_pop",
                                  full.names = FALSE)
  
  #setting path to where all the csv files due to be processed are
  setwd(perWell.resDir)
  
  
  OC.file.list <- lapply(perWell.filenames,
                         read.csv)
  
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub(".csv",
                    "",
                    perWell.filenames)
  
  names(OC.file.list)<- Exp.Names
  
  #Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
  perWell.df <- rapply(
    OC.file.list,
    as.character,
    how = "replace"
  )
  
  
  perWell.df <- bind_rows(perWell.df , .id = "Filename")
  
  
  #Houskeeping
  rm(OC.file.list)
  
  return(perWell.df)
}


# 3) Loading per Well results

if(analysis.data == "Population_&_Tracking") {
  perWell.df <- loading.perWell.df(perWell.resDir)
  

  #Creacting vector with the list of file names
  perWell.trk.filenames <- list.files(path = trk.resDir,
                                      pattern = "*_Trkv2_LC.csv",
                                      full.names = FALSE)
  
  #setting path to where all the csv files due to be processed are
  setwd(trk.resDir)
  
  
  OC.file.list <- lapply(perWell.trk.filenames,
                         read.csv)
  
  #Naming the elements of my list with the vector list of my file names
  Exp.Names <- gsub(".csv",
                    "",
                    perWell.trk.filenames)
  
  names(OC.file.list)<- Exp.Names
  
  #Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
  perWell.trk.df <- rapply(
    OC.file.list,
    as.character,
    how = "replace"
  )
  
  
  perWell.trk.df <- bind_rows(perWell.trk.df , .id = "Filename")
  
  
  #Houskeeping
  rm(OC.file.list)
  
  
}else {perWell.df <- loading.perWell.df(perWell.resDir)}


setwd(genDir)

#4) Load Condition code file

cc.filename <- unique(perWell.df$ExpFile)

cc.filename <- paste(cc.filename,
                     "_CC.csv",
                     sep="")


CC  <- read.csv(file = cc.filename) 



perWell.df <- perWell.df %>%
  mutate(Exp.ID =ExpFile)

perWell.df <- left_join(perWell.df,
                        CC)


perWell.df <- perWell.df %>%
  select(-Exp.ID)

#Houskeeping
rm(CC)


#-------------- Target well for KF analysis --------------
#finding the well in question
specific.well.df <- list.files(path = ".", pattern = "_Results", all.files = FALSE)

# 
specific.well.df <- gsub(paste(unique(perWell.df$ExpFile),
                               "_",
                               sep=""),"",specific.well.df)

specific.well.df <- gsub("_Results","",specific.well.df )



# Idenitfy all the LCF piNeg results
variables_to_select <- names(perWell.df)
variables_to_select <- variables_to_select[grepl("LC_piNeg.fraction", variables_to_select)]
variables_to_select <- variables_to_select[!grepl("RAW", variables_to_select)]


# Converting Main dataframes to long format


TKC.corr <- perWell.df %>%
  select(ExpFile,
         Condition,
         Well_coordinate,
         timestep,
         Time_Hrs,
         Filterting.timestep,
         all_of(variables_to_select))

rm(variables_to_select)


if (analysis.data == "Population_&_Tracking") {
  
  
  TKC.trk.corr <- perWell.trk.df%>%
    mutate(ExpFile = gsub("_.*$", "",Exp))%>%
    select(-Filename,
           -Exp)
  
  
  TKC.corr <- left_join(TKC.corr,
                        TKC.trk.corr)
}



#Houskeeping
rm(TKC.trk.corr,
   perWell.df,
   perWell.trk.df)

TKC.corr <-  melt( TKC.corr,
                   id.vars = c("ExpFile",
                               "Condition",
                               "Well_coordinate",
                               "timestep",
                               "Time_Hrs",
                               "Filterting.timestep"
                   ),
                   variable.name =  "Killing.Features",
                   value.name = "LC.fraction")


TKC.corr <- TKC.corr %>%
  mutate(timestep = as.numeric(timestep),
         Time_Hrs = as.numeric(Time_Hrs),
         Filterting.timestep = as.numeric(Filterting.timestep),
         LC.fraction = as.numeric(LC.fraction))



#Creating a dataframe with directory details
TKC.corr <- separate(TKC.corr, col = Condition,
                     into = c("Abx",
                              "Concentration",
                              "Isolate"),
                     sep="_")


TKC.corr$Isolate <- gsub("Iso.", "",TKC.corr$Isolate)


TKC.corr <- TKC.corr %>%
  select(-Filterting.timestep)%>%
  drop_na()

TKC.corr <- TKC.corr %>%
  filter(Killing.Features != "LC_piNeg.fraction.vs")%>%
  filter( Killing.Features != "LC_piNeg.fraction.vs.thrs")


#-------------- New LCF normalization --------------
# load table with max live cell fraction of each isolate and add it to the maind TKC.corr dataframe

# Load maxLCF table
files <- list.files(full.names = TRUE)

dates <- substr(files, nchar(files) - 11, nchar(files) - 4)
dates <- as.Date(dates, format = "%Y%m%d")

dates <- na.omit(dates)
latest_date <- max(dates)


latest_filename <- paste0("_Population_&_Tracking_maxLCF_", format(latest_date, "%Y%m%d"), ".csv")
latest_filename <-  list.files(pattern = latest_filename)
maxLCF.df <- read.csv(latest_filename)




rm(files,
   dates,
   latest_date,
   latest_filename)

# Select Time kill definitions based on type of analysis

defintions.considered.in.analysis <- unique((TKC.corr$Killing.Features))
defintions.considered.in.analysis <- as.character(defintions.considered.in.analysis)

defintions.considered.in.analysis <- c("Isolate",defintions.considered.in.analysis)

# Select isolate 
maxLCF.df <- maxLCF.df %>%
  select(all_of(defintions.considered.in.analysis))


rm(defintions.considered.in.analysis)
# left join such that I have a Max.LCF.of.Iso.at.tp0

maxLCF.df <- maxLCF.df %>%
  mutate(Isolate = gsub("Iso.","", Isolate))%>%
  filter(Isolate == unique(TKC.corr$Isolate))

maxLCF.df <- melt(  maxLCF.df,
                    id.vars =  c("Isolate"),
                    variable.name = "Killing.Features",
                    value.name = "maxLCF",
)

#6) Live cell fraction corrected wrt to highest live cell fraction of isolate

# We need to correct the live cell fraction across isolates. We will first find the maximum fraction of cells alive across the 6 LC calulations. That will be our fraction of 1. We then correct the rest 1/max(Isolate)@timepoint0 * LCf.fraction. We will then add an artificial 0th hour timepoint where we assume all cells are alive.

#Finding the maximum live cell fraction at timepoint 0 
TKC.corr.tp0 <- TKC.corr %>%
  dplyr::ungroup()%>%
  filter(timestep == 0)%>%
  left_join(maxLCF.df)%>%
  dplyr::group_by(ExpFile,
                  Isolate,
                  Killing.Features)%>%
  mutate(Max.LCF.of.Iso.at.tp0 = as.numeric(maxLCF))%>%
  select(-maxLCF)

# Houskeeping 
rm(maxLCF.df)

TKC.corr <- left_join(TKC.corr,
                      TKC.corr.tp0 )
#Housekeeping
rm( TKC.corr.tp0 )

TKC.corr <- TKC.corr %>%
  dplyr::ungroup()%>%
  dplyr::group_by(ExpFile,
                  Abx,
                  Concentration,
                  Well_coordinate,
                  Isolate,
                  Killing.Features)%>%
  mutate(LC.fraction = (1/first(Max.LCF.of.Iso.at.tp0))*LC.fraction)%>%
  ungroup()%>%
  select(ExpFile,
         Abx,
         Concentration,
         Well_coordinate,
         Isolate,
         timestep,
         Time_Hrs,
         Killing.Features,
         LC.fraction)%>%
  mutate(timestep = timestep + 1)


# Add timestep 0 and Time_hrs 0 and Live cell fraction 1 for all isolates
TKC.corr.0th.tp <- TKC.corr %>%
  select(ExpFile,
         Abx,
         Concentration,
         Well_coordinate,
         Isolate,
         Killing.Features)%>%
  distinct()%>%
  mutate(timestep = 0,
         Time_Hrs = 0,
         LC.fraction= 1)

TKC.corr <- rbind(TKC.corr,
                  TKC.corr.0th.tp)

#Housekeeping
rm(TKC.corr.0th.tp)




#--------------  End of new Normalisation 
#Plot the corrected and filtered as per event live cell fraction over time ** ADJUST

TKC.corr <- separate( TKC.corr, col = ExpFile,
                      into = c("Experiment.Type",
                               "Exp.N",
                               "Date"))

TKC.corr <-  TKC.corr %>%
  ungroup()%>%
  mutate(Exp.Iso = paste(Experiment.Type,
                         Exp.N,
                         "_",
                         "Iso.",
                         Isolate,
                         sep=""))%>%
  mutate(Abx.con = paste(Abx,
                         Concentration,
                         sep="_"))%>%
  mutate(
    Killing.Features = gsub("LC_piNeg.fraction.sc.thrs","Pop.Thrs.LCfraction.sc", Killing.Features),
    Killing.Features = gsub("LC_piNeg.fraction.sc.vs.thrs","Pop.Thrs.LCfraction.sc.vs", Killing.Features),
    Killing.Features = gsub("LC_piNeg.fraction.sc","Pop.PiClass.LCfraction.sc", Killing.Features),
    Killing.Features = gsub("LC_piNeg.fraction.sc.vs","Pop.PiClass.LCfraction.sc.vs", Killing.Features) )


# The arguments to spread():
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
TKC.corr <- spread(   TKC.corr  ,  Killing.Features, LC.fraction)



#------------------------ Killing features  ------------------------
print ( "Killing features")



# Nice AMK killing  Thresholding dont reach MDK50 or more
TKC.corr <- TKC.corr %>%
  ungroup()%>%
  dplyr::filter(Well_coordinate == specific.well.df)%>% 
  distinct()%>%
  drop_na()

#--------------  Killing feature 1: Live cell fraction at image frames  --------------
#--- <<<>>> Killing features ---<<<>>>
#Soft-coded: Killing feature 1: Live cell fraction at image frames

# Rational: I cannot anticipate the rate of imaging so hardcoding distinct time points is dangerous. I dont know whether the 12th frame will be 48, 72 hour time point- Infact I dont even know whether there will be 12 image frames at all. Therefore I need to account for this eventuality

# Find first 3 time step
img.frames <- TKC.corr %>%
  dplyr::filter(timestep >0)%>%
  dplyr::select(timestep)%>%
  dplyr::distinct()

img.frames.First.three  <- img.frames %>%
  top_n(-3)

img.frames.First.three <- img.frames.First.three$timestep

# find middle timestep
img.frames.Middle.frame <- median(img.frames$timestep)
img.frames.Middle.frame  <- round(img.frames.Middle.frame  ,digits = 0)


# find the last 2 timesteps
img.frames.last.two.frames  <- img.frames %>%
  top_n(2)

img.frames.last.two.frames <- img.frames.last.two.frames$timestep

# List of iimage frames which reflect the beginning, middle and end of the timelapse
list.of.image.frames <- c(img.frames.First.three,img.frames.Middle.frame,img.frames.last.two.frames)
list.of.image.frames <- sort(list.of.image.frames)


#Houskeeping
rm(img.frames.First.three,
   img.frames.Middle.frame,
   img.frames.last.two.frames ,
   img.frames)




KF.df.1 <-  TKC.corr %>%
  filter(timestep > 0)


KF.df.1.Long <- melt(KF.df.1,
                     id.vars = c("Experiment.Type",
                                 "Exp.N",
                                 "Date",
                                 "Isolate",
                                 "Well_coordinate",
                                 "Abx",
                                 "Concentration",
                                 "Abx.con",
                                 "Exp.Iso",
                                 "timestep",
                                 "Time_Hrs"
                     ),
                     variable.name =  "Killing.Features")



KF.df.1.Long <- KF.df.1.Long %>%
  mutate(LC.fraction = value)%>%
  select(-value)

#list of image frames to loop through
list.of.image.frames

KF.df.1  <- data.frame()


for ( i in list.of.image.frames) {
  
  approx.real.time <- i*freq.loop
  
  approx.real.time <- paste(approx.real.time,
                            "h",
                            sep="")
  
  KF.df.1.Long.sub <- KF.df.1.Long %>%
    filter(timestep == i)%>%
    mutate(KF.1.var = paste(Killing.Features,
                            "imageframe",
                            timestep,
                            "approx",
                            approx.real.time,
                            sep="."
    ))
  
  
  
  KF.df.1   <- rbind(KF.df.1  ,  KF.df.1.Long.sub)
  
  rm(  approx.real.time,
       KF.df.1.Long.sub )
  
}

KF.df.1   <-KF.df.1    %>%
  select(-Abx.con,
         -Exp.Iso,
         -timestep,
         -Time_Hrs,
         -Killing.Features)

# The arguments to spread():
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
KF.df.1  <- spread(KF.df.1  , KF.1.var, LC.fraction)

#Housekeeping
rm(KF.df.1.Long,
   i,
   list.of.image.frames)


#-------------- Killing feature 2  Minimum Duration to Kill (MDK) 25, 50, 75 and 90 % across definitions--------------

# ---
# Soft-coded: Killing feature 2  Minimum Duriation to Kill (MDK) 25, 50, 75 and 90 % across definitons


# CHECK WHAT HAPPENS IF WE DONT KILL THAT MUCH: SHould put NA 
# First detemine whether we are able to kill 25% of population:
# IF YES find between which image frames this occures in and what real time . Find the equation of the straight line and determine X value when Y is e.g. 0.75.
# IF NO declare NA-beyond scope of analysis
#Same calulation across MDK50 75 and 90


# mdk.df <- TKC.corr %>%
#   select(-Filterting.timestep)

mdk.df <- TKC.corr


mdk.df.Long <- melt(mdk.df,
                    id.vars = c("Experiment.Type",
                                "Exp.N",
                                "Date",
                                "Isolate",
                                "Well_coordinate",
                                "Abx",
                                "Concentration",
                                "Abx.con",
                                "Exp.Iso",
                                "timestep",
                                "Time_Hrs"
                    ),
                    variable.name =  "Killing.Features")
mdk.df.Long <- mdk.df.Long%>%
  mutate(LC.fraction =value)%>%
  select(-value)

# List of killing defitions to loop by
list.of.killing.def <- (unique(as.character(mdk.df.Long$Killing.Features)))
list.of.killing.def

# List of MDK values of interest
list.of.mdk.val <- c(0.75,0.5,0.25,0.10)

KF.df.2  <- data.frame()
for( i in list.of.killing.def) {
  
  mdk.df.Longsub <- mdk.df.Long %>%
    filter(Killing.Features == i)
  
  print(paste("MDK Assessing of killing Def:", i))
  
  # Lowest live cell fraction of this particular definition
  x <- min(mdk.df.Longsub$LC.fraction)
  
  
  for ( j in list.of.mdk.val ) {
    
    print(paste("For killing Def:", i, "calculating the MDK" , j ,"fraction of the population"))
    
    if( x < j) { 
      
      x1 <- mdk.df.Longsub  %>%
        mutate(drop_timestep0 = if_else(LC.fraction[1] == 1, TRUE, FALSE)) %>% # Accounts for instances where the live cell fraction for both timestep 0 and 1 are equal to 1
        filter(!drop_timestep0 | timestep > 0) %>%
        select(-drop_timestep0)%>%
        mutate(Time.point.Before = LC.fraction -j)%>%
        arrange(timestep)%>%
        filter(Time.point.Before > 0 & lead(Time.point.Before) < 0 ) %>% # Finding the first instance when a number changes from positive to negative
        filter(timestep == min(timestep))
      
      
      if ( nrow( x1 ) == 0) {
        # if nrow is 0 it means two things; for  Time.point.Before
        #  if all numbers are negative; this means the the fraction before it crossed the corresponding MDK value was 1.0 at time 0h
        # If first couple of of numebrs are negative then positive it means again the instance before the mdk was corresed was the time 0h. Usually we expect to have positive numebrs followed by sequence of negative numbers
        
        x1 <-  mdk.df.Longsub  %>%
          filter(Time_Hrs == 0)
        
        
      }
      
      
      y1 <- x1 %>%
        select(LC.fraction)
      
      x1 <- x1 %>%
        select(#timestep,
          Time_Hrs)
      
      x2 <-  mdk.df.Longsub %>%
        mutate(drop_timestep0 = if_else(LC.fraction[1] == 1, TRUE, FALSE)) %>% # Accounts for instances where the live cell fraction for both timestep 0 and 1 are equal to 1
        filter(!drop_timestep0 | timestep > 0) %>%
        select(-drop_timestep0)%>%
        filter(LC.fraction  <= j)%>%
        arrange(timestep)%>%
        mutate(Time.point.Before = LC.fraction  - j)%>% # Finding the first instance when a number changes to negative
        filter(Time.point.Before < 0) %>%
        filter(timestep == min(timestep))
      
      
      y2 <- x2 %>%
        select(LC.fraction )
      
      
      
      x2 <- x2 %>%
        select(
          Time_Hrs)
      
      #-- Creating vectors 
      
      #Equation of straight line
      x1 <- x1$Time_Hrs
      x2 <- x2$Time_Hrs
      
      y1 <- y1$LC.fraction 
      y2 <- y2$LC.fraction 
      
      #gradient of straight lin
      m <- (y2-y1)/(x2-x1)
      
      # Intercept
      c <- y1-(m*x1)
      
      
      minimum.duration.to.kill <- (j-c)/m
      
      
      mdk.df.Longsub.res  <-   mdk.df.Longsub %>%
        mutate(KF.2.mdk = paste(minimum.duration.to.kill))%>%
        ungroup()%>%
        select(-timestep,
               -Time_Hrs,
               -LC.fraction)%>%
        distinct()%>%
        mutate(MDK.var =  dplyr::if_else(j == 0.75,"MDK25",
                                         dplyr::if_else(j == 0.50, "MDK50",
                                                        dplyr::if_else(j == 0.25, "MDK75",
                                                                       dplyr::if_else(j == 0.1, "MDK90","MDK-ERROR")))))
      
      KF.df.2 <- rbind(KF.df.2,mdk.df.Longsub.res)
      
    } else {     
      
      # If we never reach the MDK value then the MDK value lies beyond the scope of this timelapse
      mdk.df.Longsub.res  <-   mdk.df.Longsub %>%
        select(-timestep,
               -Time_Hrs,
               -LC.fraction)%>%
        distinct()%>%
        mutate(KF.2.mdk = paste(">",
                                timelapse,
                                sep=""),
               MDK.var = dplyr::if_else(j == 0.75,"MDK25",
                                        dplyr::if_else(j == 0.50, "MDK50",
                                                       dplyr::if_else(j == 0.25, "MDK75",
                                                                      dplyr::if_else(j == 0.1, "MDK90","MDK-ERROR")))))
      
      KF.df.2 <- rbind(   KF.df.2 ,mdk.df.Longsub.res)
      
      
    }
    
  }
  
  
}

KF.df.2 <- KF.df.2 %>%
  mutate(MDK.res = paste(Killing.Features,
                         ".",
                         MDK.var,
                         "hrs",
                         sep=""))%>%
  select(-Killing.Features,
         -MDK.var)%>%
  distinct()

# For downstream Sanity checks (SC)
SC.KF.df.2 <- KF.df.2


# The arguments to spread():
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
# Reshape to wide format
KF.df.2 <- spread(   KF.df.2 ,  MDK.res, KF.2.mdk)

#Housekeeping
rm(mdk.df,
   mdk.df.Long,
   mdk.df.Longsub.res,
   mdk.df.Longsub,
   c,
   i,
   j,
   list.of.killing.def,
   list.of.mdk.val,
   m,
   minimum.duration.to.kill,
   x,
   x1,
   x2,
   y1,
   y2)


#-------------- Killing feature 3: real time live cell fraction --------------

# ---
#Soft-coded: Killing feature 3: real time live cell fraction

#Using similar strategy as MDK, when calculating the live cell fraction of a particular time. Find the first time-point in our data just before and just after. From the equation of the curve I solver or m and c and find the respectinve y-value at desired time
#Think about how to deal with cases that the real metadarte doesnt go to desired timepoint e.g. 48hrs then finde the next closes time 

#Killing Feaures df: Live cell fraction at image frames 2, 3,4,7, 12 and 13
#indexing goes from 1 onwards
#Live cell fraction after 3h of treatment

#Same calulation across MDK50 75 and 90
# rt.df <- TKC.corr %>%
#   select(-Filterting.timestep)

rt.df <- TKC.corr


rt.df.Long <- melt(rt.df,
                   id.vars = c("Experiment.Type",
                               "Exp.N",
                               "Date",
                               "Isolate",
                               "Well_coordinate",
                               "Abx",
                               "Concentration",
                               "Abx.con",
                               "Exp.Iso",
                               "timestep",
                               "Time_Hrs"
                   ),
                   variable.name =  "Killing.Features")
rt.df.Long <- rt.df.Long%>%
  mutate(LC.fraction =value)%>%
  select(-value)

# Creating empty results dataframe
KF.df.3 <- data.frame()

# Loop through killing defition
list.of.killing.def <- unique(as.character(rt.df.Long$Killing.Features))

# list of real time values we are interested in
list.of.real.time.values <- rt.values.interest 

#j <- list.of.real.time.values[1]



for ( i in list.of.killing.def ) {
  
  print(paste("Assessing real-time live cell fraction of killing Def:", i))
  
  # selecting specific killing feature
  rt.df.Long.sub <- rt.df.Long %>%
    dplyr::ungroup()%>%
    filter(Killing.Features == i)
  
  for (j in list.of.real.time.values) {
    print(paste("For killing Def:", i, "calculating the live cell fraction at" , j ,"hrs."))
    
    
    
    if (between(j, min(rt.df.Long.sub$Time_Hrs), max(rt.df.Long.sub$Time_Hrs))){
      
      
      # WITHIN 
      print( paste("Time of interest",j, "lies within timelapse" ,sep =" "))
      rt.df.Long.sub.rt.BEFORE <- rt.df.Long.sub %>%
        dplyr::ungroup()%>%
        dplyr::group_by(Well_coordinate)%>%
        dplyr::filter(timestep >0)%>%
        dplyr::arrange(timestep)%>%
        dplyr::mutate(Time_hrs_before = signif(j-Time_Hrs, digits = 3))%>%
        dplyr::arrange(timestep)%>%
        dplyr::filter(Time_hrs_before > 0 & lead(Time_hrs_before) < 0 ) %>% # Finding the first instance when a number changes from positive to negative
        dplyr::filter(timestep == min(timestep))
      
      
      
      rt.df.Long.sub.rt.AFTER<- rt.df.Long.sub %>%
        dplyr::ungroup()%>%
        dplyr::group_by(Well_coordinate)%>%
        filter(timestep >0)%>%
        arrange(timestep)%>%
        mutate(Time_hrs_before = signif(j-Time_Hrs, digits = 3))%>%
        filter(Time_hrs_before < 0) %>%
        filter(timestep == min(timestep))
      
      rt.df.Long.sub.rt <- rbind(rt.df.Long.sub.rt.BEFORE ,
                                 rt.df.Long.sub.rt.AFTER  )
      
      
      #Houskeeping
      rm(rt.df.Long.sub.rt.BEFORE ,
         rt.df.Long.sub.rt.AFTER)
      
      
      # Interpolate live cell fraction found within the equation of the line
      rt.df.Long.sub.rt <- rt.df.Long.sub.rt %>%
        mutate(LC.fraction.m =   (nth(LC.fraction,2)-nth(LC.fraction,1)) /(nth(Time_Hrs,2) - nth(Time_Hrs,1)),
               LC.fraction.c = (nth(LC.fraction,1)-(LC.fraction.m*nth(Time_Hrs,1))),
               LC.fraction.rt = (j*LC.fraction.m)+LC.fraction.c )%>% 
        ungroup()%>%
        mutate(KF.3.var = paste(Killing.Features,
                                ".",
                                j,
                                "h",
                                sep=""))%>%
        select(-LC.fraction.m ,
               -LC.fraction.c)%>%
        select(Experiment.Type,
               Exp.N,
               Date,
               Isolate,
               Well_coordinate,
               Abx,
               Concentration,
               LC.fraction.rt,
               KF.3.var)%>%
        distinct()
      
      
    } else { 
      
      
      print(paste("Time of interest",j, "lies BEYOND timelapse"))
      
      rt.df.Long.sub.rt <- rt.df.Long.sub %>%
        ungroup()%>%
        dplyr::group_by(Well_coordinate)%>%
        dplyr::arrange(timestep)%>%
        top_n(2, wt = timestep)
      
      # Interpolate live cell fraction found within the equation of the line
      rt.df.Long.sub.rt  <- rt.df.Long.sub.rt %>%
        ungroup()%>%
        group_by(Well_coordinate)%>%
        arrange(timestep)%>%
        mutate(LC.fraction.m =   (nth(LC.fraction,2)-nth(LC.fraction,1)) /(nth(Time_Hrs,2) - nth(Time_Hrs,1)),
               LC.fraction.c = (nth(LC.fraction,1)-(LC.fraction.m*nth(Time_Hrs,1))),
               LC.fraction.rt = (j*LC.fraction.m)+LC.fraction.c )%>% 
        ungroup()%>%
        mutate(KF.3.var = paste(Killing.Features,
                                ".",
                                j,
                                "h",
                                sep=""))%>%
        select(-LC.fraction.m ,
               -LC.fraction.c)%>%
        select(Experiment.Type,
               Exp.N,
               Date,
               Isolate,
               Well_coordinate,
               Abx,
               Concentration,
               LC.fraction.rt,
               KF.3.var)%>%
        distinct()
      
      
      
    }
    
    
    KF.df.3 <- rbind(KF.df.3,rt.df.Long.sub.rt )
    
  }
  
}

KF.df.4.to.be.used.for.KF.5 <- KF.df.3 
KF.df.4.to.be.used.for.KF.5.1 <- KF.df.3
# The arguments to spread(Data,key,value):
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
# Reshape to wide format
KF.df.3 <- spread(   KF.df.3 ,  KF.3.var, LC.fraction.rt)


#Houskeeping
rm(rt.df,
   rt.df.Long,
   rt.df.Long.sub,
   rt.df.Long.sub.rt,
   i,
   j,
   list.of.real.time.values,
   list.of.killing.def)




#--------------  Killing features 4: AUC-image frame --------------

# ---
# Soft-coded: Killing features 4: AUC-image frame

# AUC have to be calculated using the corrected Live cell fractions
# Find the calculations
#11.1) AUC using the trapezium rule: SC threshold

# A = 0.5(b+a)*h . Where h = 1 ( interval in time)
#Calculate the area under the curve between each time interval. When we sum up the intermediary trapeziums we get an estimate of the area under the curve

# auc.imgfrm.df <- TKC.corr %>%
#   select(-Filterting.timestep)

auc.imgfrm.df <- TKC.corr

auc.imgfrm.df.Long <- melt(auc.imgfrm.df,
                           id.vars = c("Experiment.Type",
                                       "Exp.N",
                                       "Date",
                                       "Isolate",
                                       "Well_coordinate",
                                       "Abx",
                                       "Concentration",
                                       "Abx.con",
                                       "Exp.Iso",
                                       "timestep",
                                       "Time_Hrs"
                           ),
                           variable.name =  "Killing.Features")


auc.imgfrm.df.Long <- auc.imgfrm.df.Long%>%
  mutate(LC.fraction =value)%>%
  select(-value)


KF.df.4.multiple.trapz <- data.frame()


# Loop through defintions
list.of.killing.def <- unique(as.character(auc.imgfrm.df.Long$Killing.Features))
list.of.killing.def
#i <- list.of.killing.def[1]


# Loop through image frames
list.of.img.frames <- unique(auc.imgfrm.df.Long$timestep)
list.of.img.frames <- sort(list.of.img.frames)
list.of.img.frames
#j <- list.of.img.frames[1]

max(list.of.img.frames)

for ( i in list.of.killing.def ) {
  
  print(paste("AUC across frame of killing Def:", i))
  
  # selecting specific killing feature
  auc.imgfrm.df.Long.sub <- auc.imgfrm.df.Long %>%
    dplyr::ungroup()%>%
    filter(Killing.Features == i)
  
  
  
  for (j in list.of.img.frames ) {
    
    h <- ifelse( j == max(list.of.img.frames), 0 , 1)
    h
    
    next.frame <- j +1
    
    
    auc.imgfrm.df.Long.sub.sub <- auc.imgfrm.df.Long.sub %>%
      ungroup()%>%
      filter(timestep >= j)%>%
      filter(timestep <=  (j+1))%>%
      group_by(Well_coordinate) %>%
      arrange(timestep)%>%
      mutate(# sum LC_piNeg is equviavalent to the sum of the (a+b) in the equation of a trapezium
        AUC.trapz = (sum(LC.fraction) /2)*h)%>%
      mutate(AUC.frame.int = paste("AUC.frameInterval.",
                                   j,
                                   "-",
                                   next.frame,
                                   sep =""))%>%
      select(-timestep,
             -Time_Hrs,
             -LC.fraction)%>%
      distinct()
    
    
    KF.df.4.multiple.trapz <- rbind(KF.df.4.multiple.trapz,  auc.imgfrm.df.Long.sub.sub)
    
  }
  
}


KF.df.4.multiple.trapz <- KF.df.4.multiple.trapz %>%
  ungroup()%>%
  mutate(Frame.Int = as.numeric(sub(".*-", "",AUC.frame.int)))%>%
  filter(Frame.Int  < max(Frame.Int))

#---- Calculating the AUC across timestep intervals
# Derive image frame inter
KF.df.4.multiple.trapz.interval <- data.frame()

set.of.intervals <- KF.df.4.multiple.trapz$Frame.Int

set.of.intervals <- as.numeric(set.of.intervals)
set.of.intervals <- sort(unique(set.of.intervals))
set.of.intervals

#From 0 to the third image frame ~ initial level of tolerance
Int.0.3 <- KF.df.4.multiple.trapz  %>%
  select(Frame.Int)%>%
  distinct()%>%
  top_n(-3)

Int.0.3 <- as.numeric(max(Int.0.3$Frame.Int))
Int.0.3

#From 0  to middle of the timelapse ~ level of tolerance half way through experiment
Int.0.middle.of.timelapse <- round(median(set.of.intervals), digits = 0)
Int.0.middle.of.timelapse

#From 0 till the end ~ level of tolerance of of the isolate
Int.0.end  <- round(max(set.of.intervals), digits = 0)

# from (last.frame-1) till end ~ level of persisters
Int.beforeEnd  <-  KF.df.4.multiple.trapz  %>%
  select(Frame.Int)%>%
  distinct()%>%
  top_n(2)%>%
  filter(Frame.Int == min(Frame.Int))


Int.beforeEnd <- Int.beforeEnd$Frame.Int
Int.beforeEnd

# List of killing definitions
# k <- list.of.killing.def[1]
# k
# Define the relevant image frames we aim to consider
relevant.auc.img.frame.intervals <- c(Int.0.3,Int.0.middle.of.timelapse ,Int.beforeEnd,  Int.0.end)
# l <- relevant.auc.img.frame.intervals[1]
# l
# Now we will sum the many tiny trapeziums across the relevant time intervals
KF.df.4.multiple.trapz.interval <- data.frame()

for ( k in list.of.killing.def ) {
  
  KF.df.4.multiple.trapz.sub <- KF.df.4.multiple.trapz %>%
    filter(Killing.Features == k)
  
  print(paste("AUC killing Def:", k))
  
  for ( l in relevant.auc.img.frame.intervals) {
    
    print(paste("AUC time interval:", l))
    
    if (l == Int.beforeEnd ) {
      
      KF.df.4.multiple.trapz.sub.sub <-  KF.df.4.multiple.trapz.sub %>%
        filter(  Frame.Int >= l)%>%
        mutate(KF.4.auc.sum = sum(AUC.trapz),
               KF.4.Intervals = paste(min(Frame.Int),
                                      "_",
                                      max(Frame.Int),
                                      sep=""))
      
      KF.df.4.multiple.trapz.interval <- rbind( KF.df.4.multiple.trapz.interval,KF.df.4.multiple.trapz.sub.sub )
      
      rm( KF.df.4.multiple.trapz.sub.sub)
      
    } else {
      
      KF.df.4.multiple.trapz.sub.sub <-  KF.df.4.multiple.trapz.sub %>%
        filter( Frame.Int <= l)%>%
        mutate(KF.4.auc.sum = sum(AUC.trapz),
               KF.4.Intervals = paste(min(Frame.Int),
                                      "_",
                                      max(Frame.Int),
                                      sep=""))
      
      KF.df.4.multiple.trapz.interval <- rbind(KF.df.4.multiple.trapz.interval, KF.df.4.multiple.trapz.sub.sub)
      
      rm( KF.df.4.multiple.trapz.sub.sub)
      
    }
    
    
  }
  
}

KF.df.4.multiple.trapz.interval <- KF.df.4.multiple.trapz.interval  %>%
  select(-AUC.trapz,
         -AUC.frame.int,
         -Frame.Int)%>%
  distinct()


KF.df.4.multiple.trapz.interval <- KF.df.4.multiple.trapz.interval  %>%
  mutate(KF.var = paste(Killing.Features,
                        "AUCimgframe",
                        KF.4.Intervals,
                        sep="."))%>%
  select(-KF.4.Intervals,
         -Killing.Features)

# The arguments to spread(Data,key,value):
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
# Reshape to wide format
KF.df.4 <- spread( KF.df.4.multiple.trapz.interval,
                   KF.var,
                   KF.4.auc.sum)


#Houskeeping
rm(KF.df.4.multiple.trapz,
   KF.df.4.multiple.trapz.interval,
   KF.df.4.multiple.trapz.sub,
   h,
   i,
   Int.0.3,
   Int.0.end,
   Int.0.middle.of.timelapse,
   Int.beforeEnd,
   j,
   k,
   l,
   list.of.killing.def,
   list.of.img.frames,
   next.frame,
   relevant.auc.img.frame.intervals,
   set.of.intervals,
   auc.imgfrm.df,
   auc.imgfrm.df.Long,
   auc.imgfrm.df.Long.sub,
   auc.imgfrm.df.Long.sub.sub)



#--------------  Killing features 5: AUC linear-scale  across real time invervals --------------

#---
#Soft-coded: Killing features 5: AUC across real time invervals

# auc.rt.df <- TKC.corr %>%
#   select(-Filterting.timestep)

auc.rt.df <- TKC.corr

auc.rt.df.Long <- melt(auc.rt.df,
                       id.vars = c("Experiment.Type",
                                   "Exp.N",
                                   "Date",
                                   "Isolate",
                                   "Well_coordinate",
                                   "Abx",
                                   "Concentration",
                                   "Abx.con",
                                   "Exp.Iso",
                                   "timestep",
                                   "Time_Hrs"
                       ),
                       variable.name =  "Killing.Features")


auc.rt.df.Long <- auc.rt.df.Long%>%
  mutate(LC.fraction =value)%>%
  select(-value)%>%
  select(-timestep)


#KF.df.4.to.be.used.for.KF.5.OG <- KF.df.4.to.be.used.for.KF.5

#KF.df.4.to.be.used.for.KF.5 <- KF.df.4.to.be.used.for.KF.5.OG

KF.df.4.to.be.used.for.KF.5$Time_Hrs <- sub(".*\\.", "", KF.df.4.to.be.used.for.KF.5$KF.3.var)


# Need to add the real time live cell fraction values to the auc.rt.df.Long dataframe. Preparing dataframe for rbind
KF.df.4.to.be.used.for.KF.5 <- KF.df.4.to.be.used.for.KF.5 %>%
  ungroup()%>%
  mutate(  Time_Hrs = gsub("h","",Time_Hrs),
           Time_Hrs = as.numeric(Time_Hrs))%>%
  mutate(Killing.Features = KF.3.var,
         String.to.remove = paste(".",
                                  Time_Hrs,
                                  "h",
                                  sep=""))%>%
  group_by(KF.3.var)%>%
  mutate(Killing.Features =gsub(String.to.remove,
                                "",
                                Killing.Features))%>%
  ungroup()

KF.df.4.to.be.used.for.KF.5 <- KF.df.4.to.be.used.for.KF.5 %>%
  mutate(LC.fraction = LC.fraction.rt)%>%
  select(-KF.3.var,
         -String.to.remove,
         -LC.fraction.rt)

KF.df.4.to.be.used.for.KF.5 <- KF.df.4.to.be.used.for.KF.5 %>%
  mutate(Killing.Features = as.factor(Killing.Features))


auc.rt.df.Long <- auc.rt.df.Long%>%
  select(-Exp.Iso,
         -Abx.con)

auc.rt.df.Long <- rbind(auc.rt.df.Long ,
                        KF.df.4.to.be.used.for.KF.5)

#need the pooled real time and measured time results for killing feature 6
To.be.used.in.KF.6 <- auc.rt.df.Long
#Houskeeping
rm( KF.df.4.to.be.used.for.KF.5)


#---- Now we will calculate the auc between 2 timepoints across the timelapse. Essentially splitting the time-kill curve into many little trapeziums. Similar strategy used in KF4 will be used here.


KF.df.5.multiple.trapz <- data.frame()


# Loop through defintions
list.of.killing.def <- unique(as.character(auc.rt.df.Long$Killing.Features))
list.of.killing.def <- sort(list.of.killing.def)
#  i <- list.of.killing.def[1]


# Loop through image frames
list.of.hrs.frames <- unique(auc.rt.df.Long$Time_Hrs)
list.of.hrs.frames <- sort(list.of.hrs.frames)
list.of.hrs.frames
#  j <- list.of.hrs.frames[1]

max(list.of.hrs.frames)

for ( i in list.of.killing.def ) {
  
  print(paste("AUC across frame of killing Def:", i)) 
  
  # selecting specific killing feature
  auc.rt.df.Long.sub <- auc.rt.df.Long %>%
    dplyr::ungroup()%>%
    filter(Killing.Features == i)
  
  # i <- list.of.killing.def[1]
  
  for (j in  list.of.hrs.frames ) {
    #j <- list.of.hrs.frames[1]
    
    auc.rt.df.Long.sub.sub <-auc.rt.df.Long.sub %>%
      ungroup()%>%
      filter(Time_Hrs >= j)%>%
      top_n(-2,Time_Hrs)%>%
      group_by(Well_coordinate) %>%
      arrange(Time_Hrs)%>%
      mutate(# difference in lead and current time_hrs gives us h
        AUC.trapz = (sum(LC.fraction) /2)*(lead(Time_Hrs)-Time_Hrs)) # We will have NAs with this calculation because It is applied across each tim_hrs 
    
    #Important to define the next time interval which indecates the range of AUC
    next.time <-  round( max(auc.rt.df.Long.sub.sub$Time_Hrs), digits = 3)
    current.time <- round(j, digits =3)
    
    auc.rt.df.Long.sub.sub <-   auc.rt.df.Long.sub.sub %>%
      mutate(AUC.realtime.int = paste("AUC.TimeInterval_",
                                      current.time, 
                                      "-",
                                      next.time,
                                      sep =""))%>%
      select(
        -Time_Hrs,
        -LC.fraction)%>%
      distinct()%>%
      drop_na() 
    
    
    KF.df.5.multiple.trapz <- rbind(KF.df.5.multiple.trapz,   auc.rt.df.Long.sub.sub )
    
    #Houskeeping
    rm(next.time,
       current.time)
    
  }
  
}

#Houskeeping
rm(auc.rt.df.Long.sub,
   auc.rt.df.Long.sub.sub)


KF.df.5.multiple.trapz <- KF.df.5.multiple.trapz  %>%
  ungroup()%>%
  mutate(Time.Int = as.numeric(sub(".*-", "",AUC.realtime.int)))

#---- Calculating the AUC across timestep intervals
# Derive image frame inter
KF.df.5.multiple.trapz.interval <- data.frame()

set.of.intervals <- KF.df.5.multiple.trapz$Time.Int 


# List of killing definitions
list.of.killing.def <- unique(as.character(KF.df.5.multiple.trapz$Killing.Features))
list.of.killing.def 

# Define the relevant image frames we aim to consider
relevant.auc.img.frame.intervals <- auc.of.interest

# Now we will sum the many tiny trapeziums across the relevant time intervals

for ( k in list.of.killing.def ) {
  
  KF.df.5.multiple.trapz.sub <-KF.df.5.multiple.trapz %>%
    filter(Killing.Features == k)
  
  print(paste("AUC killing Def:", k))
  
  for ( l in relevant.auc.img.frame.intervals) {
    
    print(paste("AUC time interval:", l))
    
    # intervals are given by e.g. 0-3. We therfore find the numeber before and after the delimiter - to define the range of AUCs of interest
    from <- as.numeric(gsub( "-.*", "", l))
    #  from
    to <- as.numeric(sub( ".*-", "", l))
    #  to
    
    KF.df.5.multiple.trapz.sub.sub <-  KF.df.5.multiple.trapz.sub %>%
      filter(  Time.Int >= from)%>%
      filter(  Time.Int <= to)%>%
      mutate(KF.5.auc.sum = sum(AUC.trapz),
             KF.5.Intervals = paste(from,
                                    "_",
                                    to ,
                                    sep=""))
    
    KF.df.5.multiple.trapz.interval <- rbind( KF.df.5.multiple.trapz.interval,KF.df.5.multiple.trapz.sub.sub )
    
    rm( KF.df.5.multiple.trapz.sub.sub)
    
  }
  
}

KF.df.5.multiple.trapz.interval <- KF.df.5.multiple.trapz.interval  %>%
  ungroup()%>%
  select(-Time.Int,
         -AUC.trapz,
         -AUC.realtime.int)%>%
  distinct()




KF.df.5.multiple.trapz.interval  <-     KF.df.5.multiple.trapz.interval  %>%
  mutate(KF.var = paste(Killing.Features,
                        "AUCrealtimeHrs",
                        KF.5.Intervals,
                        sep="."))%>%
  select(-KF.5.Intervals,
         -Killing.Features)


# Downstream santiy checks
SC.KF.df.5 <- KF.df.5.multiple.trapz.interval 

# The arguments to spread(Data,key,value):
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
# Reshape to wide format
KF.df.5 <- spread( KF.df.5.multiple.trapz.interval,
                   KF.var,
                   KF.5.auc.sum)

#Housekeeping
rm(auc.rt.df,
   auc.rt.df.Long,
   KF.df.5.multiple.trapz,
   KF.df.5.multiple.trapz.interval,
   KF.df.5.multiple.trapz.sub,
   j,
   k,l,i,
   list.of.hrs.frames,
   list.of.killing.def,
   from,
   to)



#-------------- Killing features 5.1: AUC log-scale  across real time invervals --------------

#---
#Soft-coded: Killing features 5: AUC across real time invervals


auc.rt.df <- TKC.corr

auc.rt.df.Long <- melt(auc.rt.df,
                       id.vars = c("Experiment.Type",
                                   "Exp.N",
                                   "Date",
                                   "Isolate",
                                   "Well_coordinate",
                                   "Abx",
                                   "Concentration",
                                   "Abx.con",
                                   "Exp.Iso",
                                   "timestep",
                                   "Time_Hrs"
                       ),
                       variable.name =  "Killing.Features")


auc.rt.df.Long <- auc.rt.df.Long%>%
  mutate(LC.fraction =log10(value))%>%
  select(-value)%>%
  select(-timestep)


#KF.df.4.to.be.used.for.KF.5.OG <- KF.df.4.to.be.used.for.KF.5

#KF.df.4.to.be.used.for.KF.5 <- KF.df.4.to.be.used.for.KF.5.OG

KF.df.4.to.be.used.for.KF.5.1$Time_Hrs <- sub(".*\\.", "", KF.df.4.to.be.used.for.KF.5.1$KF.3.var)


# Need to add the real time live cell fraction values to the auc.rt.df.Long dataframe. Preparing dataframe for rbind
KF.df.4.to.be.used.for.KF.5.1 <- KF.df.4.to.be.used.for.KF.5.1 %>%
  ungroup()%>%
  mutate(  Time_Hrs = gsub("h","",Time_Hrs),
           Time_Hrs = as.numeric(Time_Hrs))%>%
  mutate(Killing.Features = KF.3.var,
         String.to.remove = paste(".",
                                  Time_Hrs,
                                  "h",
                                  sep=""))%>%
  group_by(KF.3.var)%>%
  mutate(Killing.Features =gsub(String.to.remove,
                                "",
                                Killing.Features))%>%
  ungroup()

KF.df.4.to.be.used.for.KF.5.1 <- KF.df.4.to.be.used.for.KF.5.1 %>%
  mutate(LC.fraction = LC.fraction.rt)%>%
  select(-KF.3.var,
         -String.to.remove,
         -LC.fraction.rt)

KF.df.4.to.be.used.for.KF.5.1 <- KF.df.4.to.be.used.for.KF.5.1 %>%
  mutate(Killing.Features = as.factor(Killing.Features))


auc.rt.df.Long <- auc.rt.df.Long%>%
  select(-Exp.Iso,
         -Abx.con)

auc.rt.df.Long <- rbind(auc.rt.df.Long ,
                        KF.df.4.to.be.used.for.KF.5.1)

#need the pooled real time and measured time results for killing feature 6
To.be.used.in.KF.6.1 <- auc.rt.df.Long
#Houskeeping
rm( KF.df.4.to.be.used.for.KF.5.1)


#---- Now we will calculate the auc between 2 timepoints across the timelapse. Essentially splitting the time-kill curve into many little trapeziums. Similar strategy used in KF4 will be used here.


KF.df.5.1.multiple.trapz <- data.frame()


# Loop through defintions
list.of.killing.def <- unique(as.character(auc.rt.df.Long$Killing.Features))
list.of.killing.def <- sort(list.of.killing.def)
#  i <- list.of.killing.def[1]


# Loop through image frames
list.of.hrs.frames <- unique(auc.rt.df.Long$Time_Hrs)
list.of.hrs.frames <- sort(list.of.hrs.frames)
list.of.hrs.frames
#  j <- list.of.hrs.frames[1]

max(list.of.hrs.frames)

for ( i in list.of.killing.def ) {
  
  print(paste("AUC across frame of killing Def:", i)) 
  
  # selecting specific killing feature
  auc.rt.df.Long.sub <- auc.rt.df.Long %>%
    dplyr::ungroup()%>%
    filter(Killing.Features == i)
  
  # i <- list.of.killing.def[1]
  
  for (j in  list.of.hrs.frames ) {
    #j <- list.of.hrs.frames[1]
    
    auc.rt.df.Long.sub.sub <-auc.rt.df.Long.sub %>%
      ungroup()%>%
      filter(Time_Hrs >= j)%>%
      top_n(-2,Time_Hrs)%>%
      group_by(Well_coordinate) %>%
      arrange(Time_Hrs)%>%
      mutate(# difference in lead and current time_hrs gives us h
        AUC.trapz = (sum(LC.fraction) /2)*(lead(Time_Hrs)-Time_Hrs)) # We will have NAs with this calculation because It is applied across each tim_hrs 
    
    #Important to define the next time interval which indecates the range of AUC
    next.time <-  round( max(auc.rt.df.Long.sub.sub$Time_Hrs), digits = 3)
    current.time <- round(j, digits =3)
    
    auc.rt.df.Long.sub.sub <-   auc.rt.df.Long.sub.sub %>%
      mutate(AUC.realtime.int = paste("AUC.TimeInterval.LOG_",
                                      current.time, 
                                      "-",
                                      next.time,
                                      sep =""))%>%
      select(
        -Time_Hrs,
        -LC.fraction)%>%
      distinct()%>%
      drop_na() 
    
    
    KF.df.5.1.multiple.trapz <- rbind(KF.df.5.1.multiple.trapz,   auc.rt.df.Long.sub.sub )
    
    #Houskeeping
    rm(next.time,
       current.time)
    
  }
  
}

#Houskeeping
rm(auc.rt.df.Long.sub,
   auc.rt.df.Long.sub.sub)


KF.df.5.1.multiple.trapz <- KF.df.5.1.multiple.trapz  %>%
  ungroup()%>%
  mutate(Time.Int = as.numeric(sub(".*-", "",AUC.realtime.int)))

#---- Calculating the AUC across timestep intervals
# Derive image frame inter
KF.df.5.1.multiple.trapz.interval <- data.frame()

set.of.intervals <- KF.df.5.1.multiple.trapz$Time.Int 



list.of.killing.def <- unique(as.character(KF.df.5.1.multiple.trapz$Killing.Features))
list.of.killing.def 

# Define the relevant image frames we aim to consider
relevant.auc.img.frame.intervals <- auc.of.interest
#l <- relevant.auc.img.frame.intervals[1]
# l
# Now we will sum the many tiny trapeziums across the relevant time intervals

for ( k in list.of.killing.def ) {
  
  KF.df.5.1.multiple.trapz.sub <-KF.df.5.1.multiple.trapz %>%
    filter(Killing.Features == k)
  
  print(paste("AUC killing Def:", k))
  
  for ( l in relevant.auc.img.frame.intervals) {
    
    print(paste("AUC time interval:", l))
    
    # intervals are given by e.g. 0-3. We therfore find the numeber before and after the delimiter - to define the range of AUCs of interest
    from <- as.numeric(gsub( "-.*", "", l))
    #  from
    to <- as.numeric(sub( ".*-", "", l))
    #  to
    
    KF.df.5.1.multiple.trapz.sub.sub <-  KF.df.5.1.multiple.trapz.sub %>%
      filter(  Time.Int >= from)%>%
      filter(  Time.Int <= to)%>%
      mutate(KF.5.auc.sum = sum(abs(AUC.trapz)),
             KF.5.Intervals = paste(from,
                                    "_",
                                    to ,
                                    sep=""))
    
    KF.df.5.1.multiple.trapz.interval <- rbind( KF.df.5.1.multiple.trapz.interval,KF.df.5.1.multiple.trapz.sub.sub )
    
    rm( KF.df.5.1.multiple.trapz.sub.sub)
    
  }
  
}

KF.df.5.1.multiple.trapz.interval <- KF.df.5.1.multiple.trapz.interval  %>%
  ungroup()%>%
  select(-Time.Int,
         -AUC.trapz,
         -AUC.realtime.int)%>%
  distinct()


KF.df.5.1.multiple.trapz.interval  <-     KF.df.5.1.multiple.trapz.interval  %>%
  mutate(KF.var = paste(Killing.Features,
                        "AUCrealtimeHrsLOG",
                        KF.5.Intervals,
                        sep="."))%>%
  select(-KF.5.Intervals,
         -Killing.Features)


# Downstream santiy checks
SC.KF.df.5.1 <- KF.df.5.1.multiple.trapz.interval 

# The arguments to spread(Data,key,value):
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
# Reshape to wide format
KF.df.5.1 <- spread( KF.df.5.1.multiple.trapz.interval,
                     KF.var,
                     KF.5.auc.sum)

#Housekeeping
rm(auc.rt.df,
   auc.rt.df.Long,
   KF.df.5.1.multiple.trapz,
   KF.df.5.1.multiple.trapz.interval,
   KF.df.5.1.multiple.trapz.sub,
   j,
   k,l,i,
   list.of.hrs.frames,
   list.of.killing.def,
   from,
   to)


#--------------  Killing features 6: Crossing time: Persister live cell fraction and y-intercept --------------

#---
#Soft-coded: Killing features 6: Crossing time: Persister live cell fraction and y-intercept; STRATEGY 2 new Killing features using Linear regression

# Consider the first 4 points when generating the equation of the first line and the last 2 points to consider the last equation of the line. The point where the tow curves intersect will approximately show us the inflection/crossing time, which could indicate the time we start killing the persister cell population.

# The five summaries required for calculating the slope and intercept of the best fitting line are:

# The mean of x
# The mean of y
# The sd of x
# The sd of y
# The r between x and y
# We can work out the five summarises using built-in functions in R. 1) formula for the slope m, response variable first m = r(sd(x) / sd(y)) . 2) # formula for intercept as above
# b = mean(y) - (m * mean(x))

KF.df.6.EqN.linear <- To.be.used.in.KF.6


KF.df.6.Eq1.linearReg <- KF.df.6.EqN.linear%>%
  filter(Time_Hrs > 0)%>%
  ungroup()%>%
  group_by(Killing.Features)%>%
  top_n(-4,Time_Hrs)%>%
  mutate(Eq1.m = cor(LC.fraction,Time_Hrs)*sd(LC.fraction),
         Eq1.c = mean(LC.fraction)- (Eq1.m*mean(Time_Hrs)))%>%
  ungroup()%>%
  select(-Time_Hrs,
         -LC.fraction)%>%
  distinct()


KF.df.6.Eq2.linearReg <- KF.df.6.EqN.linear %>%
  ungroup()%>%
  group_by(Killing.Features)%>%
  top_n(2,Time_Hrs)%>%
  # mutate(Eq2.m = cor(LC.fraction,Time_Hrs)*sd(LC.fraction),
  #        Eq2.c = mean(LC.fraction)- (Eq2.m*mean(Time_Hrs)) )%>%
  mutate(Eq2.m = (lead(LC.fraction)-LC.fraction) /(lead(Time_Hrs)- Time_Hrs))%>%
  mutate(Eq2.c = LC.fraction - Eq2.m*Time_Hrs)%>%
  ungroup()%>%
  drop_na()%>%
  select(-Time_Hrs,
         -LC.fraction)%>%
  distinct()

KF.df.6.Eq1.Eq2.linearReg <- left_join(KF.df.6.Eq1.linearReg,
                                       KF.df.6.Eq2.linearReg)

#Houskeeping
rm(KF.df.6.Eq1.linearReg,
   KF.df.6.Eq2.linearReg,
   KF.df.6.EqN.linear)

KF.df.6.Eq1.Eq2.linearReg <- KF.df.6.Eq1.Eq2.linearReg %>%
  group_by(Killing.Features)%>%
  mutate(X.crossing_time.LinReg = (Eq2.c-Eq1.c)/ (Eq1.m-Eq2.m))%>% # Crossing time i.e time it takes till we treat most of the tolerant population
  mutate(Y.intercept.LinReg = Eq2.c )%>%# Initial fraction of persiter population 
  mutate(Persister.fraction.LinReg= Eq2.m*X.crossing_time.LinReg + Eq2.c)%>% # Fraction of persister population at crossing time
  select(-Eq1.m,
         -Eq1.c,
         -Eq2.m,
         -Eq2.c)

# Linear Crossing time (ctR)results
KF.df.6.ctr <- KF.df.6.Eq1.Eq2.linearReg %>%
  select(-Y.intercept.LinReg,
         -Persister.fraction.LinReg)%>%
  mutate(Killing.Features = as.character(paste(Killing.Features, 
                                               "Crossing.time.Hrs",sep=".")))%>%
  mutate(LC.fraction = X.crossing_time.LinReg )%>%
  select(-X.crossing_time.LinReg)


# Linear Persiter fraction (pf)
KF.df.6.pf <- KF.df.6.Eq1.Eq2.linearReg %>%
  select(-X.crossing_time.LinReg,
         -Y.intercept.LinReg)%>%
  mutate(Killing.Features = as.character(paste(Killing.Features, 
                                               "persister.fraction.at.crossingtime",sep=".")))%>%
  mutate(LC.fraction = Persister.fraction.LinReg )%>%
  select(-Persister.fraction.LinReg)


# Linear Y intercept
KF.df.6.yInt <- KF.df.6.Eq1.Eq2.linearReg %>%
  select(-X.crossing_time.LinReg,
         -Persister.fraction.LinReg)%>%
  mutate(Killing.Features = as.character(paste(Killing.Features, 
                                               "initial.persister.fraction",sep=".")))%>%
  mutate(LC.fraction = Y.intercept.LinReg )%>%
  select(-Y.intercept.LinReg)


KF.df.6  <- rbind(KF.df.6.yInt,
                  KF.df.6.pf)

KF.df.6 <- rbind(KF.df.6,
                 KF.df.6.ctr)


SC.KF.df.6 <- KF.df.6
# The arguments to spread(Data,key,value):
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
# Reshape to wide format
KF.df.6 <- spread( KF.df.6,
                   Killing.Features,
                   LC.fraction)

#Houskeeping
rm(KF.df.6.ctr,
   KF.df.6.Eq1.Eq2.linearReg,
   KF.df.6.pf,
   KF.df.6.yInt)


#-------------- Exporting Killing Features --------------

#---
#Master Killing features df: long format
#  
# MS.KF.df.Wide <- left_join(KF.df.1,
#                            KF.df.2)

MS.KF.df.Wide <- left_join(KF.df.2,
                           KF.df.3)


MS.KF.df.Wide <- left_join(MS.KF.df.Wide,
                           KF.df.4)

MS.KF.df.Wide <- left_join(MS.KF.df.Wide,
                           KF.df.5)

MS.KF.df.Wide <- left_join(MS.KF.df.Wide,
                           KF.df.5.1)

# Dropping persister and crossing time
MS.KF.df.Wide <- left_join(MS.KF.df.Wide,
                           KF.df.6)


# Adding the crossing time definition, to the main dataframe 2 with linear Regression. 
# MS.KF.df.Wide <- left_join(MS.KF.df.Wide,
#                            KF.df.6)



MS.KF.df.Wide <- MS.KF.df.Wide %>%
  mutate(Abx.con = paste(Abx,
                         Concentration,
                         sep="_"))

MS.KF.df.Long <- melt(MS.KF.df.Wide,
                      id.vars = c("Experiment.Type",
                                  "Exp.N",
                                  "Date",
                                  "Exp.Iso",
                                  "Isolate",
                                  "Well_coordinate",
                                  "Abx",
                                  "Concentration",
                                  "Abx.con"),
                      variable.name =  "Killing.Features")

#---
#Exporting files

filename.csv.KF <- paste(exp.resDir,
                         "/",
                         
                         MS.KF.df.Wide$Experiment.Type,
                         MS.KF.df.Wide$Exp.N,
                         ".",
                         MS.KF.df.Wide$Date,
                         
                         "_",
                         MS.KF.df.Wide$Well_coordinate,
                         "_",
                         MS.KF.df.Wide$Abx,
                         "_",
                         MS.KF.df.Wide$Concentration,
                         "_",
                         MS.KF.df.Wide$Isolate,
                         "_",
paste("KF.",tracking.versions,"_longformat.csv",sep=""),

                         sep="")


write.csv(MS.KF.df.Long,filename.csv.KF,row.names = FALSE)
filename.csv.KF



#--------------  Sanity check --------------
#--------------  SC 1: LCF always increasing  --------------

# grouping by the defintions and arrangeing time varibale in ascending order compare the current and leadning value
SC.KF3.Long.df <- To.be.used.in.KF.6 %>%
  ungroup()%>%
  group_by(Killing.Features)%>%
  arrange(Time_Hrs)%>%
  glimpse()%>%
  mutate(SC.KF3.diff = lead(LC.fraction)-LC.fraction )%>%
  drop_na()%>%
  mutate(SC.LCF.always.decreasing = if_else(SC.KF3.diff < 0, # We expect fractions to always be decreasing as time goes by
                                            0,#  0 means Yes its below zero
                                            1))%>%  # 1 means no its its above zero
  mutate(SC.LCF.always.decreasing = cumsum(SC.LCF.always.decreasing),
         SC.LCF.always.decreasing = if_else(sum(SC.LCF.always.decreasing) > 0 ,"no","yes"))%>%
  select(-SC.KF3.diff,
         # -Killing.Features,
         -Time_Hrs,
         #  - Real.time.fractions,
         -LC.fraction)%>%
  distinct()


SC.KF.3.df.checked <- SC.KF3.Long.df 

#Houskeeping
rm(SC.KF3.Long.df)

#--------------  SC 2: MDK  --------------
#Sanity check MDK: MDK values can only increase as a greater proportion of the population is killed

# Assess MDK values and test whether MDK values increase as a greater proportion of the population is being killed.
# MDK25 had to be shorter/less than MDK50  same true for MDK25 ,MDK50 and MDK75 and 90
# If MDK value is >96 for MDK50 then it has to also be >96 for MDK75




# Change MDK fraction variable to characher
SC.KF2.Long.df <- SC.KF.df.2 %>%
  mutate(MDK.res = as.character(MDK.res))


# Create variable with the kill curve definitnons
SC.KF2.Long.df <- SC.KF2.Long.df %>%
  mutate(Def = MDK.res)%>%
  mutate(Def = sub(".MDK.*", "",Def ))


# create a new variable with  MDK percentages only 
SC.KF2.Long.df <- SC.KF2.Long.df %>%
  mutate(MDK.perc = MDK.res)%>%
  mutate(MDK.perc = gsub(pattern = "[^0-9.-]","", MDK.perc))%>%
  mutate(MDK.perc = gsub("[[:punct:]]", "",MDK.perc))


# Change MDK values >96 to only 96 
time <- paste(">",timelapse,sep="")

time <- as.character(time)
timelapse2 <- as.character(timelapse)



SC.KF2.Long.df <- SC.KF2.Long.df %>%
  mutate(KF.2.mdk = dplyr::if_else(KF.2.mdk ==  time,timelapse2, KF.2.mdk))%>%
  mutate(KF.2.mdk = as.numeric(KF.2.mdk))

#Housekeeping
rm(time,
   timelapse2)

# Group by kill curve defnintions and test whther MDK KF.2.mdks are always increasing
SC.KF2.Long.df <- SC.KF2.Long.df %>%
  ungroup()%>%
  group_by(Def)%>%
  arrange(MDK.perc)%>%
  mutate(MDK.time.increasing = if_else(lead(KF.2.mdk) >= KF.2.mdk, 0, # 0 means it is increasing
                                       1)) # 1 means KF.2.mdks did not increase
# 
# # Check what happens when we have mutliple>96 :: make is such that when we have Low levels of killing and MDKe.g. 50 never is reached (therefore >96h) the subsequent MDK75 will also be 96h.  IF this conditon is met then we still classifiy it as a n increases

SC.KF2.Long.df <- SC.KF2.Long.df %>%
  ungroup()%>%
  drop_na()

SC.KF2.Long.df <- SC.KF2.Long.df %>%
  select(-MDK.perc,
         -KF.2.mdk)%>%
  distinct()%>%
  mutate(MDK.time.increasing = if_else(MDK.time.increasing  == 0 , "yes", "no" ))

SC.KF.2.df.checked <- SC.KF2.Long.df %>%
  select(-MDK.res)%>%
  distinct()

#Houskeeping
rm(SC.KF2.Long.df)

#--------------  SC 2: AUC linear --------------

#---
# Check AUC: area under the curve can only increases or stay the same but not decrease and sub sequentially increase

# AUC can never decrease over time. A decrease would suggest an increase and subsequent decrease in live cell fraction. Here we aim to detect whether there is a decrease in AUC over time.


SC.KF5.Long.df <- SC.KF.df.5 %>%
  ungroup()%>%
  mutate(KF.var = as.character(KF.var))%>%
  mutate(Def = KF.var)


SC.KF5.Long.df <- SC.KF5.Long.df %>%
  mutate(Def = sub("\\.AUCrealtimeHrs.*", "", Def))

SC.KF5.Long.df <- SC.KF5.Long.df %>%
  mutate(Time_Hrs = KF.var)%>%
  mutate(Time_Hrs = sub(".*AUCrealtimeHrs.0.", "", Time_Hrs))%>%
  filter(Time_Hrs == "3" | Time_Hrs == "6" | Time_Hrs =="9" | Time_Hrs == "12" | Time_Hrs == "18" |  Time_Hrs == "24" | Time_Hrs == "36" |  Time_Hrs == "48"  | Time_Hrs == "60"  |Time_Hrs == "72" ) %>%
  select(-KF.var)%>%
  mutate(Time_Hrs = as.numeric(Time_Hrs))


# Group by kill curve definitions and test whether AUC values are always increasing
SC.KF5.Long.df <- SC.KF5.Long.df %>%
  ungroup()%>%
  group_by(Def)%>%
  arrange(Time_Hrs)%>%
  mutate(AUC.time.increasing = if_else(lead(KF.5.auc.sum) >= KF.5.auc.sum, 0, # 0 means it is increaseing
                                       1)) # 1 means values did not increase

SC.KF5.Long.df <- SC.KF5.Long.df %>%
  ungroup()%>%
  drop_na()%>%
  select(-KF.5.auc.sum,
         -Time_Hrs)%>%
  distinct()%>%
  mutate(AUCrealtime.always.increasing = if_else(AUC.time.increasing  == 0 , "yes", "no" ))

SC.KF.5.df.checked <- SC.KF5.Long.df


#Houskeeping
rm(SC.KF5.Long.df)


#--------------  SC 2.1: AUC logarithmic --------------

# Check AUC: area under the curve can only increases or stay the same but not decrease and sub sequentially increase

# AUC  can never decrease over time. A decrease would suggest an increase and subsequent decrease in live cell fraction. Here we aim to detect whether there is a decrease in AUC over time.


SC.KF5.1.Long.df <- SC.KF.df.5.1 %>%
  ungroup()%>%
  mutate(KF.var = as.character(KF.var))%>%
  mutate(Def = KF.var)


SC.KF5.1.Long.df <- SC.KF5.1.Long.df %>%
  mutate(Def = sub("\\.AUCrealtimeHrsLOG.*", "", Def))

SC.KF5.1.Long.df <- SC.KF5.1.Long.df %>%
  mutate(Time_Hrs = KF.var)%>%
  mutate(Time_Hrs = sub(".*AUCrealtimeHrsLOG.0.", "", Time_Hrs))%>%
  filter(Time_Hrs == "3" | Time_Hrs == "6" | Time_Hrs =="9" | Time_Hrs == "12" | Time_Hrs == "18" |  Time_Hrs == "24" | Time_Hrs == "36" |  Time_Hrs == "48"  | Time_Hrs == "60"  |Time_Hrs == "72" ) %>%
  select(-KF.var)%>%
  mutate(Time_Hrs = as.numeric(Time_Hrs))


# Group by kill curve definitions and test whether AUC values are always increasing
SC.KF5.1.Long.df <- SC.KF5.1.Long.df %>%
  ungroup()%>%
  group_by(Def)%>%
  arrange(Time_Hrs)%>%
  mutate(AUC.time.increasing = if_else(lead(KF.5.auc.sum) >= KF.5.auc.sum, 0, # 0 means it is increaseing
                                       1)) # 1 means values did not increase

SC.KF5.1.Long.df <- SC.KF5.1.Long.df %>%
  ungroup()%>%
  drop_na()%>%
  select(-KF.5.auc.sum,
         -Time_Hrs)%>%
  distinct()%>%
  mutate(AUCrealtimeLOG.always.increasing = if_else(AUC.time.increasing  == 0 , "yes", "no" ))

SC.KF.5.1.df.checked <- SC.KF5.1.Long.df


#Houskeeping
rm(SC.KF5.1.Long.df)


#--------------  SC 3: Crossing time --------------

#---
# Check Crossing time: Live-cell intercept : Crossing time : Persister live cell fraction

SC.KF6.Long.df <- SC.KF.df.6 %>%
  ungroup()%>%
  mutate(Abx.con = paste(Abx,
                         Concentration,
                         sep="_"))%>%
  select(-Abx,
         -Concentration)%>%
  mutate(Crossing.Time.results = Killing.Features)%>%
  select(-Killing.Features)

# Initial Persister fraction
SC.KF6.Long.initial.pers.fraction.df <- SC.KF6.Long.df %>%
  mutate(Initial.Persister = as.character(Crossing.Time.results),
         Initial.Persister.eval = grepl("initial.persister.fraction",Initial.Persister ))%>% #grepl(needle, haystack, fixed=TRUE)
  filter(Initial.Persister.eval == TRUE )%>%
  mutate(Initial.Persister = as.numeric(LC.fraction),
         Initial.Persister.above.1 = if_else(Initial.Persister > 1, "yes","no"))

SC.KF6.Long.initial.pers.fraction.df <- SC.KF6.Long.initial.pers.fraction.df %>%
  select(-LC.fraction,
         -Initial.Persister.eval)%>%
  mutate(Initial.persister.fraction = Crossing.Time.results)

# Persister fraction at crossing time: Can never be >96 
SC.KF6.Long.fraction.at.crossingtime.df <- SC.KF6.Long.df %>%
  mutate(fraction.at.crossingtime = as.character(Crossing.Time.results),
         fraction.at.crossingtime.eval = grepl("persister.fraction.at.crossingtime",fraction.at.crossingtime))%>% #grepl(needle, haystack, fixed=TRUE)
  filter(fraction.at.crossingtime.eval == TRUE )%>%
  mutate(fraction.at.crossingtime.above.1 = if_else(LC.fraction >1 ,"yes","no"))


SC.KF6.Long.fraction.at.crossingtime.df <- SC.KF6.Long.fraction.at.crossingtime.df%>%
  select(-Crossing.Time.results)


# Crossing time 
SC.KF6.Long.crossingtime.df <- SC.KF6.Long.df %>%
  mutate(Crossingtime = as.character(Crossing.Time.results),
         Crossingtime.eval = grepl("Crossing.time.Hrs",Crossingtime))%>% #grepl(needle, haystack, fixed=TRUE)
  filter(Crossingtime.eval == TRUE )%>%
  mutate(Crossingtime.above.end = if_else(LC.fraction > max( TKC.corr$Time_Hrs) ,"yes","no"))

#Houskeeping
rm(SC.KF6.Long.df)




SC.KF.df.6.checked <- SC.KF6.Long.crossingtime.df %>%
  select(-Crossingtime,
         -Crossing.Time.results,
         -Crossingtime.eval,
         -LC.fraction)

rm(SC.KF6.Long.crossingtime.df)

# --------------------- Exporting SC ------------


SC.KF.df.6.checked <- left_join(SC.KF.df.6.checked,SC.KF6.Long.fraction.at.crossingtime.df)

rm(SC.KF6.Long.fraction.at.crossingtime.df)



SC.KF.6.df.checked <- left_join(SC.KF.df.6.checked,SC.KF6.Long.initial.pers.fraction.df )


rm(SC.KF6.Long.initial.pers.fraction.df )



print("All sanity checks complete!")

#---
# Combine all sanity checks

SC.POOLED <- left_join(SC.KF.2.df.checked,
                       SC.KF.3.df.checked)

SC.POOLED  <- left_join(SC.POOLED ,SC.KF.5.df.checked)

SC.POOLED  <- left_join(SC.POOLED ,SC.KF.5.1.df.checked)


SC.POOLED  <- left_join(SC.POOLED ,SC.KF.6.df.checked)

#Housekeeping
rm(SC.KF.2.df.checked,
   SC.KF.3.df.checked,
   SC.KF.5.df.checked,
   SC.KF.6.df.checked)

SC.POOLED <- SC.POOLED %>%
  distinct()%>%
  select(-LC.fraction,
         -Initial.Persister)%>%
  distinct()%>%
  
  select(Experiment.Type,
         Exp.N,
         Date,
         Isolate,
         Well_coordinate,
         Killing.Features,
         MDK.time.increasing,
         AUCrealtime.always.increasing,
         AUCrealtimeLOG.always.increasing,
         SC.LCF.always.decreasing)%>% # Dropped the Persister time and crossing time
  distinct()


SC.POOLED <- SC.POOLED %>%
  distinct()%>%
  ungroup()



#---
# Exporting sanity checks

filename.csv.KF <- paste(exp.resDir,
                         "/",
                         
                         MS.KF.df.Wide$Experiment.Type,
                         MS.KF.df.Wide$Exp.N,
                         ".",
                         MS.KF.df.Wide$Date,
                         
                         "_",
                         MS.KF.df.Wide$Well_coordinate,
                         "_",
                         MS.KF.df.Wide$Abx,
                         "_",
                         MS.KF.df.Wide$Concentration,
                         "_",
                         MS.KF.df.Wide$Isolate,
                         "_",paste("KF.",tracking.versions,"_sanity.check_longformat.csv",sep=""),
                    #    "KF_sanity.check_longformat.csv",
                         sep="")


filename.csv.KF

write.csv(SC.POOLED,filename.csv.KF,row.names = FALSE)
filename.csv.KF





