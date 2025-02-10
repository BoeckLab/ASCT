## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------------
# SECTION 1: defining variables ----

old.or.new.tracking <-c("Trkv2")
tracking.availablitiy <- ("Yes") # options Yes or No. If yes in SECTION 8, tracking killing feautres will be loaded

Two.BaSic.def <- ("Yes") # options Yes or No. If yes then the ordering of killing features export table will also include the two basic outputs.

analysis.normalisation <- "Top2norm"
#General variables
genDir <- getwd() 

# EDIT the time scale axis that should be used 
plot.xaxis.time.scale <- c(0,12,24,36,48,60,72)

plot.xaxis.time.scale.timestep <- c(0:29)

# Are we using a 1526 well plate or 384 well plate? If 1536 well plate variable has to be set as
type.of.wellplate <- c(1,48,1)

#Chose between these two variables; type.of.data c("Pop.only")  ; type.of.data c("Pop_&_Trk")
type.of.data <- c("Pop.only")
path.to.minimum.LCF.data <-c("/Users/Alex/Documents/PhD/PhD_project/PhD_project-Tolerance/R-ANALYSIS/Ranalysis_perwell/Ranalysis_perWell_TLKK/ASCT_Analysis_AJ_orginal/ASCT_EXPERIMENTS_v2-critical/ASCT_Initial_LCF/Experimental-Results/ASCT.Avg_Isolates_minLCF_20231001.csv")

path.to.isol.too.few.numbers.LCF.data <-c("/Users/Alex/Documents/PhD/PhD_project/PhD_project-Tolerance/R-ANALYSIS/Ranalysis_perwell/Ranalysis_perWell_TLKK/ASCT_Analysis_AJ_orginal/ASCT_EXPERIMENTS_v2-critical/ASCT_Initial_LCF/Experimental-Results/ASCT.Avg_Isolates_too.few.numbers_20231001.csv")


maxLCF.df.datapath <-c("/Users/Alex/Documents/PhD/PhD_project/PhD_project-Tolerance/R-ANALYSIS/Ranalysis_perwell/Ranalysis_perWell_TLKK/ASCT_Analysis_AJ_orginal/ASCT_EXPERIMENTS_v2-critical/ASCT_Initial_LCF/Experimental-Results/ASCT.Avg_Population_&_Tracking_maxLCF_20231001.csv")

# ---- CC file information
# Navigate back one directory (assuming you want to go up one level)
cc_dir <- file.path(genDir, "..")

# Enter the ASCT.Ranalysis_Exp_Info directory
cc_dir <- file.path(cc_dir, "ASCT.Ranalysis_Exp_Info")

# Enter the ASCT_CC directory
cc_dir <- file.path(cc_dir, "ASCT_CC")


# ---- label directory file information to be used in SECTION 
# Navigate back one directory (assuming you want to go up one level)
label_dir <- file.path(genDir, "..")

# Enter the ASCT.Ranalysis_Exp_Info directory
label_dir <- file.path(label_dir , "ASCT.Ranalysis_Exp_Info")

# Enter the ASCT_CC directory
label_dir <- file.path(label_dir , "ASCT_Label_directory")

  
  
#Pop data variables
res <- c("ASCT_Data/PerWell-Pop.Results")

#Tracking data
res.tracking.data <-c("ASCT_Data/PerWell-Trk.Results")
res.tracking.data.numbers <-c("ASCT_Data/PerWell-Trk.n.Results")

exp.res <- c("Experimental-Results")


#Directories
wdDir <- genDir

perWell.resDir <- paste(genDir,"/",res,sep="")
res.tracking.data.Dir <-  paste(genDir,"/",res.tracking.data,sep="")
res.tracking.data.N.Dir <-  paste(genDir,"/",res.tracking.data.numbers ,sep="")


exp.resDir <- paste(genDir,"/", exp.res,sep="")


#Housekeeping
rm(mocSimple,
   pocSimple,
   qc,
   metaData,
   res)
#Setting work directory
setwd(wdDir)


#SECTION 2:  Loading packages ----

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
library(scales)


# Capture the start time
start_time <- Sys.time()



# SECTION 3: Loading population data -----------
loading.perWell.df <- function (perWell.resDir){
  
#Creacting vector with the list of file names
perWell.filenames <- list.files(path = perWell.resDir,
                            pattern = "*.csv",
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


perWell.df <- loading.perWell.df(perWell.resDir)




# SECTION 4:  Export for Ahmad zeroBytes files ----------------------------------------------------------------------------------------------------------------------------------------------------------

ZB.df <- perWell.df %>%
  select(ExpFile,
         Well_coordinate,
         MOC.zeroBytes.fields)

zb.exp <- unique(ZB.df$ExpFile)

experiment.id <- zb.exp


ZB.df <- ZB.df %>%
  drop_na()%>%
  distinct()




 # SECTION 5: Load CC condiiton code  ----------------------------------------------------------------------------------------------------------------------------------------------------------

CC.filename <- unique(perWell.df$ExpFile)
CC.filename <- paste(CC.filename,
                     "_CC.csv",
                     sep="")

CC.filename <- file.path(cc_dir , CC.filename)
CC  <- read.csv(file = CC.filename) 

perWell.df <- left_join(perWell.df,
                        CC)

#Houskeeping
rm(CC,
   CC.filename)

  Conditions <- perWell.df %>%
  select(ExpFile,
         Well_coordinate,
         Condition)%>%
  distinct()


  Conditions <- separate(Conditions, col = Condition,
                      into = c("Abx",
                               "Concentration",
                               "Isolate"),
                      sep="_")

  Conditions <- Conditions %>%
  mutate(ExpFile = paste(ExpFile,
                         Well_coordinate,
                         sep ="_"))


  perWell.df <- perWell.df %>%
  group_by(Well_coordinate) %>%
  filter(!is.na(Time_Hrs)) %>%
  ungroup()


  
# SECTION 6:  Events ----
  ## 6.1 Psudo Growth curves: ---- 
  
# An overview of growth curves
Growth.df <- perWell.df %>%
  select(
         Well_coordinate,
         timestep,
         Time_Hrs,
         GT_CorrObjArea,
         GT_Growth.eval.in.well)%>%
  mutate(timestep = as.numeric(timestep),
         Time_Hrs = as.numeric(Time_Hrs),
         GT_CorrObjArea = as.numeric(GT_CorrObjArea),
        GT_Growth.eval.in.well = GT_Growth.eval.in.well)%>%
  left_join(Conditions , by = c("Well_coordinate"))%>%
  mutate(ExpFile = sub("_.*","", ExpFile))%>%
  mutate(Abx.con = paste(Abx,
                         Concentration,sep="_"))

# Plot psudo growth curve
# PLOTTING TIME KILL CURVES
plot.title.growth <- paste(exp.resDir,"/",unique(Growth.df$ExpFile),"_Psudo_GrowthCurves.pdf", sep = "")
#LC_plottitle <- paste ( "Time kill curves")
gg <-Growth.df%>%
  ggplot(aes(x = Time_Hrs,
             y = GT_CorrObjArea,
             group = Well_coordinate))+
  geom_line(size = 0.2,
            alpha = 0.5)+
  #    scale_x_continuous( breaks = plot.xaxis.time.scale)+
  geom_dl(aes(label = Well_coordinate), method = list(dl.combine( "last.points"),cex =0.1) ) +
  theme(legend.position="none")+
  theme_bw()+
     theme(aspect.ratio = 1)+
 labs(title = "Psudo growth curves",
       x = "Time [Hrs]",
       y = "Ratio tot.obj.Area/total obj.N [a.u]",
     #  color = "Isolate",
      subtitle = "Rel line is the growth threshold Calculation; Corr.N.tot.ObjArea = (TotObj.area/TotObj.n)/(dplyr::first(TotObj.area)/dplyr::first(TotObj.n)) ")+
   theme(legend.position="none",
        axis.text.x = element_text(size = 3),
        axis.text.y = element_text(size = 3),
        plot.subtitle = element_text(size = 2)
          )+
  geom_hline(yintercept = 6.5, color = "red", linetype = "dashed", size = 0.2, alpha = 0.5) +  # Add 
  
    scale_y_continuous(breaks = seq(2, 12, by = 2))+

    coord_cartesian(ylim = c(0, 12))+  # Set Y-axis limits to 0-1.5
 facet_wrap_paginate(~ Abx.con + Isolate,
                      ncol =7,
                      nrow =7,
                      page =1)
n <- n_pages(gg)

pdf(plot.title.growth,paper = "a4", width = 20 , height = 15 )
for(i in 1:n){
    print(gg + facet_wrap_paginate(~ Abx.con + Isolate,
                      ncol =7,
                      nrow =7, page = i)) ## WORKS with 240 wells ## I think Error: Cannot create zero-length unit vector ("unit" subsetting) means it does not know how to plot the remianing 2 wells on the last page given ncol 4 and nrow 6
}
dev.off()

#Houskeeping
rm(gg,
   LiveCells)


rm(plot.title.growth,
   Growth.df)





### 6.1.1 ROC analysis based growth evaluation ----
# Susceptibiltiy, resistance and heteroresistance evaluation


Growth.ResHR.df <- perWell.df


Growth.ResHR.df  <-Growth.ResHR.df  %>%
      mutate(Exp.Well = paste(ExpFile,
                          Well_coordinate,
                          sep="_"))%>%
  mutate(      timestep = as.numeric(timestep),
               Time_Hrs = as.numeric(      Time_Hrs))%>%
       select(
         ExpFile,
         Exp.Well,
         Well_coordinate,
         timestep,
         Time_Hrs,
         GT_TotObj.area,
    QC_n.SC.VS.cells,
    QC_n.SC)




Master.PerWell.df.12_72h <- Growth.ResHR.df

# Find closest value  before and after 24hrs in each
list.of.Exp.Well <- Master.PerWell.df.12_72h %>%
  select(Exp.Well)%>%
  distinct()

total_iterations <- nrow(list.of.Exp.Well)

list.of.Exp.Well <- unique(list.of.Exp.Well$Exp.Well)

#i <- c("ASCT.22.20230814_D13")

#timepoints.of.interest <- c(0,72)
timepoints.of.interest <- c(12,72)


loop.progress <- 0
Results.df <- data.frame()


### 6.2.1 Considering 12-72h ----
### 6.2.3 Long Analaysis: Readouts where max time.hrs is within timelapse ----

Results.df.within <- data.frame()
Results.df.beyond <-  data.frame()
    
Master.PerWell.df.subset <- melt(   Master.PerWell.df.12_72h , id.vars = c("ExpFile",
                                                                      "Well_coordinate",
                                                                      "Exp.Well" ,
                                                                     # "Well_coordinate",
                                                                      "timestep", "Time_Hrs" ),
                                variable.name = "Feature")


# Identify the well coordinates with the last time hrs before 72h

  Master.PerWell.df.subset <-   Master.PerWell.df.subset  %>%
    dplyr:: mutate(value = as.numeric(value))%>%
    dplyr::mutate(value = if_else(is.na(value), 1,value)) %>%
    ungroup()%>%
    group_by(Exp.Well)%>%
    mutate(Within.or.beyond = if_else(max(Time_Hrs) <= 72 , "Beyond", "Within"))%>%
    ungroup()%>%
    filter(Within.or.beyond == "Within")


    Master.PerWell.df.subset$Feature <- as.character(  Master.PerWell.df.subset$Feature)


     Master.PerWell.df.subset.subset <- Master.PerWell.df.subset %>%
      ungroup()%>%
      distinct()%>%
      dplyr:: mutate(value = as.numeric(value))%>%
      dplyr::mutate(value = if_else(is.na(value), 1,value)) #%>%
   #  filter( Feature == element)

     for ( tp in timepoints.of.interest ) {

       print(tp)
  # tp <- timepoints.of.interest[1]


# If it lies within timelapse
    Master.PerWell.df.subset.subset.subset <-    Master.PerWell.df.subset.subset


    BEFORE <- Master.PerWell.df.subset.subset.subset%>%
            dplyr::ungroup()%>%
        dplyr::group_by(ExpFile,
                        Well_coordinate,
                        Feature)%>%
        dplyr::filter(timestep >0)%>%
        dplyr::arrange(timestep)%>%
        dplyr::mutate(Time_hrs_before = signif(tp -Time_Hrs, digits = 3))%>%
        dplyr::arrange(timestep)%>%
        dplyr::filter(Time_hrs_before > 0 & lead(Time_hrs_before) < 0 ) %>% # Finding the first instance when a number changes from positive to negative
        dplyr::filter(timestep == min(timestep))



    AFTER <- Master.PerWell.df.subset.subset.subset %>%
        dplyr::ungroup()%>%
           dplyr::group_by(ExpFile,
                        Well_coordinate,
                        Feature)%>%
        filter(timestep >0)%>%
        arrange(timestep)%>%
        mutate(Time_hrs_before = signif(tp -Time_Hrs, digits = 3))%>%
        filter(Time_hrs_before < 0) %>%
        filter(timestep == min(timestep))

  BEFORE.AFTER <- rbind(  BEFORE  ,
                                 AFTER  )


      #Houskeeping
      rm(  BEFORE ,
           AFTER )

          # Interpolate live cell fraction found within the equation of the line
  BEFORE.AFTER <-  BEFORE.AFTER %>%
    ungroup()%>%
    group_by(Exp.Well,
             Feature)%>%
        mutate(value.m =   (nth(value,2)-nth(value,1)) /(nth(Time_Hrs,2) - nth(Time_Hrs,1)),
               value.c = (nth(value,1)-(value.m*nth(Time_Hrs,1))),
               value.rt = (tp*value.m)+value.c )%>%

        ungroup()%>%

        distinct()

    BEFORE.AFTER.results <-   BEFORE.AFTER %>%
      select(-timestep)%>%
      mutate(Feature.val =   value.rt )%>%
      mutate(Time_Hrs =tp )%>%
      select(ExpFile,
             Well_coordinate,
             Exp.Well,
             Time_Hrs,
           Feature,
           Feature.val)%>%
      distinct()

    # Saving results in Results table
      Results.df.within  <- rbind(Results.df.within ,
                               BEFORE.AFTER.results)

     }





#### 6.2.4 Long Analaysis: Beyond timelapse  ---------------
Master.PerWell.df.subset <- melt(   Master.PerWell.df.12_72h , id.vars = c("ExpFile",
                                                                      "Well_coordinate",
                                                                      "Exp.Well" ,
                                                                     # "Well_coordinate",
                                                                      "timestep", "Time_Hrs" ),
                                variable.name = "Feature")


# Identify the well coordinates with the last time hrs before 72h

  Master.PerWell.df.subset <-   Master.PerWell.df.subset  %>%
    dplyr:: mutate(value = as.numeric(value))%>%
    dplyr::mutate(value = if_else(is.na(value), 1,value)) %>%
    ungroup()%>%
    group_by(Exp.Well)%>%
    mutate(Within.or.beyond = if_else(max(Time_Hrs) <= 72 , "Beyond", "Within"))%>%
    ungroup()%>%
    filter(Within.or.beyond == "Beyond")


    Master.PerWell.df.subset$Feature <- as.character(  Master.PerWell.df.subset$Feature)


    Master.PerWell.df.subset <- Master.PerWell.df.subset %>%
      ungroup()%>%
      distinct()%>%
      dplyr:: mutate(value = as.numeric(value))%>%
      dplyr::mutate(value = if_else(is.na(value), 1,value)) #%>%
   #  filter( Feature == element)



       list.of.Exp.Well <- Master.PerWell.df.subset%>%
            select(Exp.Well) %>%
            distinct()

        list.of.Exp.Well  <-        list.of.Exp.Well$Exp.Well


#i <- list.of.Exp.Well[1]


  for ( i in list.of.Exp.Well ) {
#  loop.progress
    loop.progress <- loop.progress +1
      progress_percentage <- round((loop.progress / total_iterations) * 100,digits = 1)


    if (    progress_percentage %% 2 == 0) {
  print(paste(progress_percentage, "% complete", sep = ""))
      }


      Master.PerWell.df.subset.subset <- Master.PerWell.df.subset %>%
        ungroup()%>%
        filter(Exp.Well == i)

     for ( tp in timepoints.of.interest ) {



       Master.PerWell.df.subset.subset


# If it lies within timelapse
     if (between(tp, min(  Master.PerWell.df.subset.subset$Time_Hrs), max(  Master.PerWell.df.subset.subset$Time_Hrs))){

       Master.PerWell.df.subset.subset.subset <-    Master.PerWell.df.subset.subset

    BEFORE <- Master.PerWell.df.subset.subset.subset%>%
            dplyr::ungroup()%>%
        dplyr::group_by(ExpFile,
                        Well_coordinate,
                        Feature)%>%
        dplyr::filter(timestep >0)%>%
        dplyr::arrange(timestep)%>%
        dplyr::mutate(Time_hrs_before = signif(tp -Time_Hrs, digits = 3))%>%
        dplyr::arrange(timestep)%>%
        dplyr::filter(Time_hrs_before > 0 & lead(Time_hrs_before) < 0 ) %>% # Finding the first instance when a number changes from positive to negative
        dplyr::filter(timestep == min(timestep))



    AFTER <- Master.PerWell.df.subset.subset.subset %>%
        dplyr::ungroup()%>%
           dplyr::group_by(ExpFile,
                        Well_coordinate,
                        Feature)%>%
        filter(timestep >0)%>%
        arrange(timestep)%>%
        mutate(Time_hrs_before = signif(tp -Time_Hrs, digits = 3))%>%
        filter(Time_hrs_before < 0) %>%
        filter(timestep == min(timestep))

  BEFORE.AFTER <- rbind(  BEFORE,
                                 AFTER  )


      #Houskeeping
      rm(  BEFORE ,
           AFTER )

          # Interpolate live cell fraction found within the equation of the line
  BEFORE.AFTER <-  BEFORE.AFTER %>%
    ungroup()%>%
    group_by(Exp.Well,
             Feature)%>%
        mutate(value.m =   (nth(value,2)-nth(value,1)) /(nth(Time_Hrs,2) - nth(Time_Hrs,1)),
               value.c = (nth(value,1)-(value.m*nth(Time_Hrs,1))),
               value.rt = (tp*value.m)+value.c )%>%

        ungroup()%>%

        distinct()

    BEFORE.AFTER.results <-   BEFORE.AFTER %>%
      select(-timestep)%>%
      mutate(Feature.val =   value.rt )%>%
      mutate(Time_Hrs =tp )%>%
      select(ExpFile,
             Well_coordinate,
             Exp.Well,
             Time_Hrs,
           Feature,
           Feature.val)%>%
      distinct()

    # Saving results in Results table
      Results.df.beyond <- rbind(    Results.df.beyond ,
                               BEFORE.AFTER.results)
     }else {


   #   print(paste("Time of interest",j, "lies BEYOND timelapse"))

     BEFORE.AFTER <- Master.PerWell.df.subset.subset.subset %>%
        ungroup()%>%
        dplyr::group_by(ExpFile,
                        Well_coordinate,
                        Feature)%>%
        dplyr::arrange(timestep)%>%
        top_n(2, wt = timestep)

      # Interpolate live cell fraction found within the equation of the line
          BEFORE.AFTER   <-      BEFORE.AFTER  %>%
        ungroup()%>%
       dplyr::group_by(ExpFile,
                        Well_coordinate,
                        Feature)%>%
        arrange(timestep)%>%
       group_by(Feature)%>%
        mutate(value.m =   (nth(value,2)-nth(value,1)) /(nth(Time_Hrs,2) - nth(Time_Hrs,1)),
               value.c = (nth(value,1)-(value.m*nth(Time_Hrs,1))),
               value.rt = (tp*value.m)+value.c )%>%

        distinct()

    BEFORE.AFTER.results <-   BEFORE.AFTER %>%
      select(-timestep)%>%
      mutate(Feature.val =   value.rt )%>%
      mutate(Time_Hrs =tp )%>%
      select(ExpFile,
             Well_coordinate,
             Exp.Well,
             Time_Hrs,
           Feature,
           Feature.val)%>%
      distinct()

    # Saving results in Results table
     Results.df.beyond <- rbind(    Results.df.beyond  ,
                               BEFORE.AFTER.results)




  }

  }
}



  

   rm(BEFORE.AFTER,
   BEFORE.AFTER.results,
   Results.df,
   Master.PerWell.df.subset.subset.subset,
   Master.PerWell.df.subset.subset,
   Master.PerWell.df.subset)
   

rm(Master.PerWell.df.subset)

rm(BEFORE.AFTER,
   BEFORE.AFTER.results)

Results.df <- rbind(Results.df.within,
                    Results.df.beyond)

rm(Results.df.beyond,
   Results.df.within)



Master.PerWell.df.12_72h <- melt(  Master.PerWell.df.12_72h , id.vars = c("ExpFile",
                                                                      "Well_coordinate",
                                                                      "Exp.Well" ,
                                                                      "timestep", 
                                                                      "Time_Hrs" ),
                                variable.name = "Feature")

Master.PerWell.df.12_72h <- Master.PerWell.df.12_72h %>%
  ungroup()%>%
  distinct()%>%
  mutate(Feature = as.character(Feature))%>%
  mutate(Feature.val = as.numeric(value))%>%
  select(-value)


Master.PerWell.df <- Master.PerWell.df.12_72h %>%
  select(-timestep,
         -Exp.Well)


Results.df<- Results.df %>%
  select(-Exp.Well)

Master.PerWell.df <- rbind(Master.PerWell.df,
                           Results.df)

# --- 

rm(Master.PerWell.df.12_72h)


# Considering new Ratio from 6 to 72 hours
Master.PerWell.df <- Master.PerWell.df %>%
  ungroup()%>%
  filter(Time_Hrs >= 12)%>%
    filter(Time_Hrs <= 72)%>%
    dplyr::mutate(Feature.val = if_else(is.na(Feature.val), 1,Feature.val))

#Houskeeping
rm(i,
   progress_percentage,
   total_iterations,
   tp,
   loop.progress)

Growth.ResHR.df <- Master.PerWell.df

rm(Master.PerWell.df)


### 6.2.5 Growth detection Feature 1: MinRatio;n.SC;c:none -----

GrowthFeature.1 <-Growth.ResHR.df %>%
  ungroup()%>%
  filter( Feature == "QC_n.SC")%>%
  group_by(Well_coordinate,
           Feature)%>%
  arrange(Time_Hrs)%>%
  mutate(Ratio = Feature.val/ first(Feature.val),
         MaxRatio= max(Ratio),
         MinRatio = min(Ratio))%>%
  mutate(Min.Ratio_n.SC_corr.None = MinRatio)%>%
  ungroup()%>%
  select(ExpFile,
         Well_coordinate,
       Min.Ratio_n.SC_corr.None)%>%
  distinct()


### 6.2.6 Growth detection Feature 2: MaxRatio;TotalArea;c:n.SC_VS ----

GrowthFeature.2 <- Growth.ResHR.df

GrowthFeature.2  <- spread(GrowthFeature.2  , key = Feature,
                     value =  Feature.val)


GrowthFeature.2  <- GrowthFeature.2  %>%
  group_by(ExpFile,
           Well_coordinate)%>%
  arrange(Time_Hrs)%>%
  mutate( TotArea_corr.nSCVS =  (GT_TotObj.area /QC_n.SC.VS.cells) / (dplyr::first(GT_TotObj.area)/dplyr::first(QC_n.SC.VS.cells) ))%>%
  ungroup()%>%
  select(ExpFile,
         Well_coordinate,
         Time_Hrs,
         TotArea_corr.nSCVS)
  
GrowthFeature.2 <- melt( GrowthFeature.2  , id.vars =  c ("ExpFile",
                                              "Well_coordinate",
                                              "Time_Hrs"),
                      variable.name = "Feature",
                   value.name =  "Feature.val")

GrowthFeature.2 <- GrowthFeature.2 %>%
        ungroup()%>%
  group_by(Well_coordinate,
           Feature)%>%
  arrange(Time_Hrs)%>%
  mutate(
         MaxRatio= max(Feature.val),
         MinRatio = min(Feature.val))%>%
  ungroup()%>%
  mutate(Max.R_TotArea_corr.nSCVS = MaxRatio)%>%
  select(ExpFile,
         Well_coordinate,
        Max.R_TotArea_corr.nSCVS)%>%
  distinct()





### 6.2.7 Merging feature 1 & feature 2 ----------------------------------------------------------------------------------------------------------------------------------------------------------
GrowthFeature = left_join(GrowthFeature.1,
                          GrowthFeature.2)


rm(GrowthFeature.1,
   GrowthFeature.2,
   Growth.ResHR.df)



### 6.2.8 Applying thresholds ----

# growth.threshold 1 
# Min.Ratio_n.SC_corr.None
AI.scDensity.threshold <- log10(0.286630) # Derived from Groundtruth pooled analysis MIC growth evaluation

# Max.R_TotArea_corr.nSCVS
AI.avgArea.threshold <- log10(1.2208106)

Manual.scDensity.threshold <- log10(0.316227766) # -0.5
Manual.avgArea.threshold <- log10(3.1622776602) # 0.5


GrowthFeature <- GrowthFeature %>%
  ungroup()%>%
  group_by(ExpFile,
           Well_coordinate)%>%
  mutate(Min.Ratio_n.SC_corr.None = log10(Min.Ratio_n.SC_corr.None),
         Max.R_TotArea_corr.nSCVS = log10(Max.R_TotArea_corr.nSCVS))%>%
  mutate(AI.Growth.Eval = if_else( Max.R_TotArea_corr.nSCVS <= AI.avgArea.threshold & Min.Ratio_n.SC_corr.None >= AI.scDensity.threshold, "S" ,
                                if_else(Max.R_TotArea_corr.nSCVS >= AI.avgArea.threshold  & Min.Ratio_n.SC_corr.None <= AI.scDensity.threshold  , "R",
                                        if_else(Max.R_TotArea_corr.nSCVS >= AI.avgArea.threshold  & Min.Ratio_n.SC_corr.None >= AI.scDensity.threshold  ,"HR", "S"))))%>%
  
   mutate(Manual.Growth.Eval = if_else( Max.R_TotArea_corr.nSCVS <= Manual.avgArea.threshold & Min.Ratio_n.SC_corr.None >= Manual.scDensity.threshold, "S" ,
                                if_else(Max.R_TotArea_corr.nSCVS >= Manual.avgArea.threshold & Min.Ratio_n.SC_corr.None <= Manual.scDensity.threshold , "R",
                                        if_else(Max.R_TotArea_corr.nSCVS >= Manual.avgArea.threshold & Min.Ratio_n.SC_corr.None >= Manual.scDensity.threshold ,"HR", "S"))))%>%
  ungroup()


#### 6.2.9 Export new results of new growth threshold --------
# Export 

write.csv(GrowthFeature,file = paste(exp.resDir,
                          "/",
                          unique(GrowthFeature$ExpFile),
                          "_GrowthEval_ROCbased.csv",
                          sep=""),
          row.names = FALSE)





#### 6.3.0 Merge Growth eval to main perWell.df ----------------------------------------------------------------------------------------------------------------------------------------------------------
# KEY STEP everything > 0.5 Max.R_TotArea_corr.nSCVS is defined as growth
GrowthFeature.merge <- GrowthFeature %>%
  mutate(GT_Growth.eval.in.well = if_else(Manual.Growth.Eval == "S" | Manual.Growth.Eval == "HR" , "no-growth-detected", "growth-detected"))%>%
  select(ExpFile,
         Well_coordinate,
         GT_Growth.eval.in.well)

perWell.df <- perWell.df %>%
  select(-GT_Growth.eval.in.well)

perWell.df <- left_join(perWell.df,
                             GrowthFeature.merge)

rm(GrowthFeature.merge)




## Geldetachment re-assessment  ----------------------------------------------------------------------------------------------------------------------------------------------------------
  Geldetachment.df <- perWell.df %>%
   select(ExpFile,
        Well_coordinate,
        timestep,
      #  Time_Hrs,
        QC_n.OFFC,
        QC_n.SC.VS.cells,
        QC_n.SC,
        QC_n.VS,
        QC_n.c2.c5.c20,
        QC_TotArea.sqrPx.sc,
        QC_TotArea.sqrPx.vs,
        QC_TotArea.sqrPx.c2,
        QC_TotArea.sqrPx.offc,
        QC_TotArea.sqrPx.c5,
        QC_TotArea.sqrPx.c20)%>%
    mutate(QC_n.OFFC = as.numeric(QC_n.OFFC))%>%
    mutate(QC_n.SC.VS.cells = as.numeric(QC_n.SC.VS.cells))%>%
    mutate(QC_n.SC = as.numeric(QC_n.SC))%>%
    mutate(QC_n.VS = as.numeric(QC_n.VS))%>%
    mutate(QC_n.c2.c5.c20 = as.numeric(QC_n.c2.c5.c20))%>%
    mutate(QC_n.OFFC = as.numeric(QC_n.OFFC))%>%
    mutate(QC_n.SC.VS.cells = as.numeric(QC_n.SC.VS.cells))%>%
    mutate(QC_n.SC = as.numeric(QC_n.SC))%>%
    mutate(QC_n.VS = as.numeric(QC_n.VS))%>%
    mutate(QC_n.c2.c5.c20 = as.numeric(QC_n.c2.c5.c20))%>%
    mutate(QC_TotArea.sqrPx.sc = as.numeric(QC_TotArea.sqrPx.sc))%>%
    mutate(QC_TotArea.sqrPx.vs = as.numeric(QC_TotArea.sqrPx.vs))%>%
    mutate(QC_TotArea.sqrPx.c2 = as.numeric(QC_TotArea.sqrPx.c2))%>%
    mutate(QC_TotArea.sqrPx.offc = as.numeric(QC_TotArea.sqrPx.offc))%>%
      mutate( QC_TotArea.sqrPx.c5 = as.numeric( QC_TotArea.sqrPx.c5))%>%
      mutate(  QC_TotArea.sqrPx.c20 = as.numeric(  QC_TotArea.sqrPx.c20))%>%

  mutate(QC_TotArea.sqrPx.c5 = as.numeric(QC_TotArea.sqrPx.c5))%>%
   mutate(QC_TotArea.sqrPx.c20 = as.numeric(QC_TotArea.sqrPx.c20))%>%
   replace(is.na(.), 0)%>%
#  filter(timestep == "0")%>%
  mutate(timestep = as.numeric(timestep),
         Number.of.C2.C5.C20 = as.numeric( QC_n.c2.c5.c20),
        Number.of.Out.Of.Focus.cells = as.numeric(QC_n.OFFC),
        Number.of.SingleCells = as.numeric(QC_n.SC),
         Number.of.Vsnapp = as.numeric(QC_n.VS),
       Number.of.SingleCells.Vsnaps = as.numeric(QC_n.SC.VS.cells),
       QC_TotArea.sqrPx.sc = as.numeric(QC_TotArea.sqrPx.sc),
       QC_TotArea.sqrPx.vs = as.numeric(QC_TotArea.sqrPx.vs),
       QC_TotArea.sqrPx.c2 = as.numeric(QC_TotArea.sqrPx.c2),
       QC_TotArea.sqrPx.c5 = as.numeric(QC_TotArea.sqrPx.c5),
       QC_TotArea.sqrPx.c20 = as.numeric(QC_TotArea.sqrPx.c20))%>%
  select(-QC_n.OFFC,
        -QC_n.SC.VS.cells,
        -QC_n.SC,
        -QC_n.VS)%>%
  group_by(ExpFile,
    Well_coordinate,
           timestep)%>%
  mutate(Tot.Number.of.Objects = sum(Number.of.C2.C5.C20,
                                    Number.of.Out.Of.Focus.cells,
                                     Number.of.SingleCells.Vsnaps),
         Tot.Area.of.Objects = sum(QC_TotArea.sqrPx.sc,
                                  QC_TotArea.sqrPx.vs,
                                   QC_TotArea.sqrPx.c2,
                                   QC_TotArea.sqrPx.c5,
                                    QC_TotArea.sqrPx.c20)) %>%
  select(ExpFile,
        Well_coordinate,
        timestep,
        Tot.Number.of.Objects,
        Tot.Area.of.Objects )%>%
  mutate(Norm.A =  Tot.Area.of.Objects /Tot.Number.of.Objects)%>%
  ungroup()%>%
  group_by(ExpFile,
           Well_coordinate)%>%
  arrange(timestep)%>%
  mutate(Ratio.Area = Norm.A/ first(Norm.A) )
# 
#   Geldetachment.df %>%
#   ggplot(aes(x = timestep,
#              y = Ratio.Area,
#              group = Well_coordinate))+
#   geom_line(alpha = 0.2)+
#   ylim(0, 3)+
#   theme_bw()+
#   scale_x_continuous(limits = c(min(plot.xaxis.time.scale.timestep), max(plot.xaxis.time.scale.timestep)), breaks = seq(0, max(plot.xaxis.time.scale.timestep), by = 1)) 

  Geldetachment.df.Min <- Geldetachment.df %>%
  ungroup()%>%
  group_by(ExpFile,
           Well_coordinate)%>%
  mutate(Minimum.Ratio.in.timelapse = min(Ratio.Area))%>%
  select(ExpFile,
         Well_coordinate,
         Minimum.Ratio.in.timelapse)%>%
  distinct()%>%
  mutate(GelDet.eval.2 = if_else(  Minimum.Ratio.in.timelapse < (1/3), "Det", "No-det"))

  
perWell.df <- perWell.df %>%
  left_join(Geldetachment.df.Min)


perWell.df <- perWell.df %>%
  mutate(GD_eval.gelDet = GelDet.eval.2)%>%
  select(-Minimum.Ratio.in.timelapse,
         -GelDet.eval.2)


#HOuskeeping
rm(Geldetachment.df,
   Geldetachment.df.Min)


## Contamination  -----------------------------

Contamination.df <- perWell.df %>%
  select(ExpFile,
         Well_coordinate,
         Condition,
         timestep,
         Contamination_ratio,
         Contamination_eval)%>%
  drop_na()%>%
  mutate(Contamination_ratio = as.numeric(Contamination_ratio),
         timestep = as.numeric(timestep))


## Quality control ----


### Plot heatmap: Population and tracking cell numbers  ----
QC.n <- perWell.df %>%
  select(ExpFile,
        Well_coordinate,
        timestep,
        Time_Hrs,
        QC_n.OFFC,
        QC_n.SC.VS.cells,
        QC_n.SC,
        QC_n.VS)%>%
  
   mutate(QC_n.OFFC = as.numeric(QC_n.OFFC))%>%
  mutate(QC_n.SC.VS.cells = as.numeric(QC_n.SC.VS.cells))%>%
  mutate(QC_n.SC = as.numeric(QC_n.SC))%>%
    mutate(QC_n.VS = as.numeric(QC_n.VS))%>%
   replace(is.na(.), 0)%>%
#  filter(timestep == "0")%>%
  mutate(
        Number.of.Out.Of.Focus.cells = as.numeric(QC_n.OFFC),
        Number.of.SingleCells = as.numeric(QC_n.SC),
         Number.of.Vsnapp = as.numeric(QC_n.VS),
       Number.of.SingleCells.Vsnaps = as.numeric(QC_n.SC.VS.cells))


write.csv(QC.n,file = paste(exp.resDir,
                          "/",
                          unique(QC.n$ExpFile),
                          "_cell.N.allframes.csv",
                          sep=""),
          row.names = FALSE)


QC.n <- QC.n %>%
  filter(timestep == "0")%>%
  select(-Time_Hrs)%>%
  select(-QC_n.SC.VS.cells,
         -QC_n.OFFC,
         -QC_n.SC,
         -QC_n.VS)



QC.N <- melt(QC.n,
           id.vars = c("ExpFile",
                                  "Well_coordinate",
                       "timestep"),
                      variable.name =  "QC.n")


QC.N <- QC.N %>%
  mutate(value = as.numeric(value))

write.csv(QC.n,file = paste(exp.resDir,
                          "/",
                          unique(QC.n$ExpFile),
                          "_cell.N.csv",
                          sep=""),
          row.names = FALSE)

 max(QC.N$value)
 min(QC.N$value)
 
 #Converting Well Names from A01 to A1
QC.N$Well_coordinate <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", QC.N$Well_coordinate, perl = TRUE)

QC.N$Well_coordinate <- sub('(?<![0-9])0*(?=[0-9])', '',QC.N$Well_coordinate , perl=TRUE)

qc.Condition <- Conditions %>%
  mutate(ExpFile = sub("_.*","", ExpFile ))

QC.N <- left_join(QC.N,
                  qc.Condition)

#Housekeeping
rm(qc.Condition)


if ( type.of.data == "Pop.only") {
  
  pdf(paste(exp.resDir,"/",unique(QC.N$ExpFile),"_QC.MOCn.pdf",
          sep=""))
  print(
  QC.N %>%
  mutate(Isolate = sub("Iso.","",Isolate))%>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  mutate(Exp.Colunm = paste(ExpFile,
                            Colunm,
                            sep="_"))%>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  mutate(Exp.Well = paste(ExpFile,
                          Well_coordinate,
                          sep="_"))%>%
  ggplot(aes(x = Colunm,
             y= Row,
             label = Well_coordinate))+
  geom_tile(aes(fill = value),width=0.7, height=0.7,color="black")+
  geom_text(aes(label = Well_coordinate), color = "black", size = 0.2, nudge_y = 0.2)+
  geom_text(aes(label = value), color = "black", size = 0.2 )+
  geom_text(aes(label = paste(Isolate,sep=".")), color = "black", size = 0.2, nudge_y = -0.2 )+
     labs(title = paste(unique(QC.N$Abx),
                        ".",unique(QC.N$Concentration),
                        " cell numbers",sep=""))+
  coord_equal()+
 
  theme_bw()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+
  scale_y_discrete(limits=rev)+
  scale_x_continuous( breaks = seq(type.of.wellplate ),
                      limits=c(0,NA))+
  scale_fill_distiller(type = "seq", palette = "RdYlBu",
                       limits = c( min( QC.N$value), max( QC.N$value)))+
  theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+  
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))+
  facet_wrap(~QC.n)
  )
  dev.off()


}


### Plot growth, artefact, gel-detachment and contamination  ----

if ( type.of.data == "Pop.only") {
  
  QC <- perWell.df %>%
  select(ExpFile,
         Well_coordinate,
         QC_Well.image.quality,
         GT_Growth.eval.in.well,
         Contamination_eval,
         GD_eval.gelDet)%>%
  distinct()
  
  QC <- QC %>%
    mutate(QC_Well.image.quality = if_else(QC_Well.image.quality  == "good", "Keep", "Omit"),
           GT_Growth.eval.in.well = if_else(GT_Growth.eval.in.well == "no-growth-detected", "Keep", "Omit"),
          Contamination_eval = if_else(Contamination_eval == "No-contamination", "Keep", "Omit"),
          GD_eval.gelDet = if_else(GD_eval.gelDet == "No-det", "Keep", "Omit"))
  
    QC <- QC %>%
      mutate(Fluorescence.Artefact = QC_Well.image.quality,
             Growth =GT_Growth.eval.in.well,
             Gel.detachment = GD_eval.gelDet,
             Contamination = Contamination_eval)%>%
      select(-QC_Well.image.quality ,
             -GT_Growth.eval.in.well,
             -Contamination_eval,
             -GD_eval.gelDet)

    QC <- QC %>%
      ungroup()%>%
      group_by(ExpFile,
               Well_coordinate)%>%
     # mutate(Overall.assessment = if_else(Growth == "Keep" & Fluorescence.Artefact == "Keep" & Gel.detachment == "Keep" & Contamination == "Keep", 0,1))%>%
      mutate(Overall.assessment = if_else(Growth == "Keep" & Gel.detachment == "Keep" & Contamination == "Keep", 0,1))%>%

      mutate(Overall.assessment = if_else(Overall.assessment == 0 , "Keep", "Omit" ))

    QC <- melt(QC,
           id.vars = c("ExpFile",
                                  "Well_coordinate"),
           
                      variable.name =  "QC.readouts")
    
    
    qc.Condition <- Conditions %>%
  mutate(ExpFile = sub("_.*","", ExpFile ))
    
    
        QC <- QC %>%
          left_join(qc.Condition)
        
        
        # Adding info of population analysis i.e cell numbers whether we have at least 1k cells or not
        
        QC.n.thrs <- QC.n %>%
          select(ExpFile,
                 Well_coordinate,
                 Number.of.SingleCells.Vsnaps)%>%
          mutate(value = if_else(Number.of.SingleCells.Vsnaps > 1000, "Keep", "Omit"))%>%
          left_join(qc.Condition)%>%
          mutate(QC.readouts = "1000 Single Cell and Vsnapp population analysis")%>%
          select(-Number.of.SingleCells.Vsnaps)
        
      QC <- rbind(QC,QC.n.thrs )
      
      QC.subset <- pivot_wider(QC, names_from = QC.readouts, values_from = value)
        
        QC.subset <-       QC.subset%>%
        ungroup()%>%
        group_by(Well_coordinate)%>%
        mutate(Experiment.quality = if_else(Overall.assessment == "Keep" & `1000 Single Cell and Vsnapp population analysis` == "Keep", "Keep","Omit" ))%>%
          select(-Overall.assessment)

            QC.subset <-       QC.subset%>%
              drop_na()

              QC.subset<- QC.subset %>%
                mutate(Fluorescence.Artefact  = as.character(Fluorescence.Artefact),
                       Growth = as.character(Growth),
                       Gel.detachment = as.character(Gel.detachment),
                       Contamination = as.character(Contamination),
                       `1000 Single Cell and Vsnapp population analysis`= as.character(`1000 Single Cell and Vsnapp population analysis`),
                       Experiment.quality = as.character(Experiment.quality))


       QC<- melt(QC.subset,
           id.vars = c("ExpFile",
                                  "Well_coordinate",
                      "Abx",
                      "Concentration",
                      "Isolate"),
           
                      variable.name =  "QC.readouts")
  
  
    QC <- QC %>%
      mutate(value = factor(value, levels = c("Keep","Omit")))
        
        QC <- QC %>%
          mutate(Isolate = sub("Iso.","",Isolate))
        
        
        gg.experiment.title <- paste(unique(QC$ExpFile),
                                     " ",
                                     unique(QC$Abx),
                                     ".",
                                     unique(QC$Concentration),
                                     sep="") 
          
          
  Flagged.Analysis <- QC
          

}




### Plot heatmap ratio of sc and vs   ---------------------------------------------------------------------------------------------------------------------------------------------------------
QC <- perWell.df %>%
  select(ExpFile,
         Well_coordinate,
       timestep,
       QC_n.SC,
       QC_n.VS)%>%
  filter(timestep == "0")%>%
  mutate(QC_n.SC = as.numeric(QC_n.SC),
           QC_n.VS = as.numeric(QC_n.VS),
        Ratio.of.SingleCells.Vsnaps = QC_n.SC / QC_n.VS)%>%
  mutate(Ratio.of.SingleCells.Vsnaps.thrs.5x = if_else( Ratio.of.SingleCells.Vsnaps >= 5, 5,  Ratio.of.SingleCells.Vsnaps))%>%
  mutate(Ratio.of.SingleCells.Vsnaps.thrs.10x = if_else( Ratio.of.SingleCells.Vsnaps >= 10, 10,  Ratio.of.SingleCells.Vsnaps))%>%
  select(-QC_n.SC,
         -QC_n.VS,
         -timestep)



QC <- melt(QC,
           id.vars = c("ExpFile",
                                  "Well_coordinate"),
                      variable.name =  "QC.n")


#QC$Well_coordinate <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", QC$Well_coordinate, perl = TRUE)

Conditions2 <- Conditions %>%
  ungroup()%>%
  rowwise()%>%
  mutate(ExpFile = sub(paste("_",Well_coordinate,sep=""), "", ExpFile))

QC <- QC %>%
  left_join(Conditions2,by = c("ExpFile","Well_coordinate"))


pdf(paste(exp.resDir,"/",unique(QC$ExpFile),"_QC.Ratio.pdf",
          sep=""))
  QC %>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  mutate(Exp.Colunm = paste(ExpFile,
                            Colunm,
                            sep="_"))%>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  mutate(Exp.Well = paste(ExpFile,
                          Well_coordinate,
                          sep="_"))%>%
  ggplot(aes(x = Colunm,
             y= Row))+
    ggtitle(paste(
          "Ration sc/vs at t0",          sep=""))+
 geom_tile(aes(fill =  value),width=0.7, height=0.7,color="black")+
  coord_equal()+
  theme_bw()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+
  geom_text(aes(label = round(value,digits = 3)), color = "black", size = 0.5) +
     geom_text(aes(label = Well_coordinate), color = "black", size = 0.2, nudge_y = 0.2)+
   geom_text(aes(label = paste(Isolate,sep=".")), color = "black", size = 0.2, nudge_y = -0.2 )+
   scale_y_discrete(limits=rev)+
   scale_x_continuous( breaks = seq(type.of.wellplate ),
                      limits=c(0,NA))+
    scale_fill_distiller(
                       type = "seq", palette = "Spectral",
                      # limits = c(0,5),
                       guide ="colourbar",
                      limits = c(0,10),
                       breaks = c(0,
                                1,2,3,4,5,6,7,8,9,10))+
                             #  1,2,3,4,5)
                     
                   #  )+
   theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
 theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+  
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))+
    facet_wrap(~QC.n,nrow = 3,ncol = 1)
dev.off()

#Houskeeping
#rm(QC)


#HOuskeeping
Conditions <- Conditions2

rm(Conditions2)


## Brightfield background intensity  ----------------------------------------------------------------------------------------------------------------------------------------------------------
perWell.resDir.BFbck <- gsub("PerWell-Pop.Results","PerWell-BFInt.Results",perWell.resDir)


perWell.resDir.BFbck
#Creacting vector with the list of file names
perWell.filenames <- list.files(path = perWell.resDir.BFbck ,
                            pattern = "*.csv",
                            full.names = FALSE)

#setting path to where all the csv files due to be processed are
setwd(perWell.resDir.BFbck )


OC.file.list <- lapply(perWell.filenames,
                      read.csv)

setwd(genDir )
#Naming the elements of my list with the vector list of my file names
Exp.Names <- gsub(".csv",
                  "",
                  perWell.filenames)

Exp.Names <- gsub("ASCT_BFint_results_",
                  "",
                 Exp.Names)

names(OC.file.list)<- Exp.Names

#Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list
perWell.BFbck.df <- rapply(
  OC.file.list,
  as.character,
  how = "replace"
)


perWell.BFbck.df<- bind_rows(perWell.BFbck.df , .id = "Filename")


#Houskeeping
rm(OC.file.list)

  perWell.BFbck.df <-  perWell.BFbck.df %>%
    mutate(Mean.Brightfield.Background.Int = as.numeric(Mean.Brightfield.Background.Int))%>%
  ungroup()%>%
  mutate(Mean.Brightfield.Background.Int = as.numeric(Mean.Brightfield.Background.Int))%>%
  mutate(Norm.BFbck = (Mean.Brightfield.Background.Int/ 4095)*100)%>%
  group_by(Well_coordinate)%>%
  mutate(Max.NormBf = max(Norm.BFbck),
         Min.NormBf = min (Norm.BFbck),
         Diff.Max.Min.NormBF =Max.NormBf-Min.NormBf )%>%
    mutate(QC_BFbck = if_else(Max.NormBf > 75, "Omit", 
                              if_else(Min.NormBf < 25, "Omit", "Keep")))
  
  Rate.of.Change.BFbck <- perWell.BFbck.df 

pdf(paste(exp.resDir,"/",unique(perWell.df$ExpFile),"_QC_BFbck.pdf",
          sep=""))

   perWell.BFbck.df %>%
  distinct()%>%
  mutate(timestep = as.numeric(timestep))%>%
  select(-timestep)%>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  ggplot(aes(x = Colunm,
             y= Row))+
 geom_tile(aes(fill =QC_BFbck),width=0.7, height=0.7,color="black")+
  coord_equal()+
  theme_pubclean()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+
   scale_y_discrete(limits=rev)+
   scale_x_continuous( breaks = seq(1,48,1),
                      limits=c(0,NA))+
   theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
  ggtitle("Qualtiy control Brightfield background assessment")+
 theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))
dev.off()

##### per Field coordinate X-Y space  ----

perWell.BFbck.df <- perWell.BFbck.df %>%
      dplyr::mutate(Field = gsub("[^0-9.-]", "", Field))%>%
    dplyr::mutate(Field = as.numeric(Field))%>%
    dplyr::mutate(Field_x = if_else(Field == 1 | Field == 6 | Field == 7, 1,
                                 if_else(Field == 2 | Field == 5 | Field == 8  , 2, 
                                         if_else(Field == 3 | Field == 4 | Field == 9 , 3,0))))%>%
    dplyr::mutate(Field_y = if_else(Field == 7 | Field == 8 | Field == 9, 1 ,
                                 if_else(Field == 6 | Field == 5 | Field == 4 , 2,
                                         if_else( Field == 1 | Field == 2 | Field== 3 ,3 ,0))))%>%
  ungroup()%>%
   group_by(Well_coordinate,
            Field)%>%
  mutate(Mean.NormBf = mean(Norm.BFbck))%>%
  select(Filename,
         ExpFile,
         Well_coordinate,
         Field,
         Field_x,
         Field_y,
         Mean.NormBf)%>%
  distinct()


pdf(paste(exp.resDir,"/",unique(perWell.df$ExpFile),"_QC_HM_BF.field.bckInt.pdf",
          sep=""))  

    perWell.BFbck.df %>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.character(Row))%>%
  mutate(Col = as.numeric(gsub("[^0-9.-]", "",Well_coordinate)))%>%
  ggplot(aes(x = Field_x,
             y= Field_y,
             group = Well_coordinate))+
   geom_tile(aes(fill = Mean.NormBf))+
  theme_classic()+
  geom_text(aes(label = Field), color = "black", size = 0.5) +
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"))+
  theme(legend.position = "right")+
  scale_fill_distiller(type = "seq", palette = "RdYlBu",
                      limits = c(0,100))+

  scale_x_continuous( breaks = seq(1,4,1),
                      limits=c(0,NA))+
   scale_x_continuous( breaks = seq(1,4,1),
                      limits=c(0,NA))+
   theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
   theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(legend.position="top")+
    theme(strip.text.x = element_text(size = 5))+
     theme(strip.text.y = element_text(size = 5))+
  labs(title = "Heatmap mean Brightfield background Int. per field over time")+
     coord_equal()+
    theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank())+
      theme(legend.text=element_text(size=3))+
    theme(strip.background = element_blank())+
      theme(plot.background = element_rect(fill = "lightgrey"))+
      theme(strip.text.y = element_text(size = 3))+
  facet_grid(cols= vars(Col),rows = vars(Row))

dev.off()

#HOuskeeping
rm(perWell.BFbck.df)


## Fluorescence (raw) background intensity  -----

perWell.resDir.FLuncorr <- gsub("PerWell-Pop.Results","PerWell-FLInt.Results",perWell.resDir)

perWell.resDir.FLuncorr


#Creacting vector with the list of file names
perWell.filenames <- list.files(path = perWell.resDir.FLuncorr ,
                            pattern = "*.csv",
                            full.names = FALSE)

#setting path to where all the csv files due to be processed are
setwd(perWell.resDir.FLuncorr )


OC.file.list <- lapply(perWell.filenames,
                      read.csv)

setwd(genDir )
#Naming the elements of my list with the vector list of my file names
Exp.Names <- gsub(".csv",
                  "",
                  perWell.filenames)

Exp.Names <- gsub("ASCT_FLuncorr_results_",
                  "",
                 Exp.Names)

names(OC.file.list)<- Exp.Names

#Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list
perWell.FLuncorr.df <- rapply(
  OC.file.list,
  as.character,
  how = "replace"
)


perWell.FLuncorr.df <- bind_rows(perWell.FLuncorr.df , .id = "Filename")


  perWell.FLuncorr.df <- perWell.FLuncorr.df %>%
    mutate(Mean.FLuncorr.Background.Int = as.numeric(Mean.FLuncorr.Background.Int))%>%
  ungroup()%>%
  mutate(Mean.FLuncorr.Background.Int = as.numeric(Mean.FLuncorr.Background.Int))%>%
  mutate(Norm.FLbck = (Mean.FLuncorr.Background.Int/ 4095)*100)%>%
  group_by(Well_coordinate)%>%
  mutate(Max.NormFL = max(Norm.FLbck),
         Min.NormFL = min (Norm.FLbck),
         Diff.Max.Min.NormFL =Max.NormFL-Min.NormFL )%>%
    mutate(QC_FLbck = if_else(Max.NormFL >= 30, "Omit",
                              if_else(Min.NormFL >= 25, "Omit", "Keep")))

  pdf(paste(exp.resDir,"/",unique(perWell.df$ExpFile),"_QC_FLbck.pdf",
          sep=""))

  perWell.FLuncorr.df %>%
  distinct()%>%
  mutate(timestep = as.numeric(timestep))%>%
  select(-timestep)%>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  ggplot(aes(x = Colunm,
             y= Row))+
 geom_tile(aes(fill =QC_FLbck),width=0.7, height=0.7,color="black")+
  coord_equal()+
  theme_pubclean()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+
   scale_y_discrete(limits=rev)+
   scale_x_continuous( breaks = seq(1,48,1),
                      limits=c(0,NA))+

   theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
  ggtitle("Quality control Fluorescence background assessment")+
 theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))
dev.off()



##### per Field -----
  perWell.FLuncorr.df  <-   perWell.FLuncorr.df  %>%
      dplyr::mutate(Field = gsub("[^0-9.-]", "", Field))%>%
    dplyr::mutate(Field = as.numeric(Field))%>%
    dplyr::mutate(Field_x = if_else(Field == 1 | Field == 6 | Field == 7, 1,
                                 if_else(Field == 2 | Field == 5 | Field == 8  , 2,
                                         if_else(Field == 3 | Field == 4 | Field == 9 , 3,0))))%>%
    dplyr::mutate(Field_y = if_else(Field == 7 | Field == 8 | Field == 9, 1 ,
                                 if_else(Field == 6 | Field == 5 | Field == 4 , 2,
                                         if_else( Field == 1 | Field == 2 | Field== 3 ,3 ,0))))%>%
  ungroup()%>%
   group_by(Well_coordinate,
            Field)%>%
  mutate(Mean.NormFL = mean(Norm.FLbck))%>%
  select(Filename,
         ExpFile,
         Well_coordinate,
         Field,
         Field_x,
         Field_y,
         Mean.NormFL)%>%
  distinct()



pdf(paste(exp.resDir,"/",unique(perWell.df$ExpFile),"_QC_HM_FL.field.bckInt.pdf",
          sep=""))

  perWell.FLuncorr.df  %>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.character(Row))%>%
  mutate(Col = as.numeric(gsub("[^0-9.-]", "",Well_coordinate)))%>%
  ggplot(aes(x = Field_x,
             y= Field_y,
             group = Well_coordinate))+
   geom_tile(aes(fill = Mean.NormFL))+
   theme_classic()+
    geom_text(aes(label = Field), color = "black", size = 0.5) +
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.5, "cm"),
  legend.key.width = unit(0.5,"cm"))+
  theme(legend.position = "right")+
  scale_fill_distiller(type = "seq", palette = "RdYlBu",
                      # limits = c(minInt,maxInt))+
                      limits = c(0,100))+

  scale_x_continuous( breaks = seq(1,4,1),
                      limits=c(0,NA))+
   scale_x_continuous( breaks = seq(1,4,1),
                      limits=c(0,NA))+
   theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
   theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(legend.position="top")+
    theme(strip.text.x = element_text(size = 5))+
     theme(strip.text.y = element_text(size = 5))+
  labs(title = "Heatmap mean Fluorescence background Int. per field over time")+
     coord_equal()+
    theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank())+
      theme(legend.text=element_text(size=3))+
    theme(strip.background = element_blank())+
      theme(plot.background = element_rect(fill = "lightgrey"))+
      theme(strip.text.y = element_text(size = 3))+
  facet_grid(cols= vars(Col),rows = vars(Row))

dev.off()



#### Plot out-of-focus  morphology class  ----

  perWell.MOC.pred.fraction <- perWell.df %>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         QC_n.OFFC,
         QC_n.SC,
         QC_n.VS,
         QC_n.c2.c5.c20)%>%
  mutate(QC_n.OFFC = as.numeric(QC_n.OFFC))%>%
  mutate(QC_n.SC = as.numeric(QC_n.SC))%>%
    mutate(QC_n.VS = as.numeric(QC_n.VS))%>%
    mutate(QC_n.c2.c5.c20 = as.numeric(QC_n.c2.c5.c20))%>%
 replace(is.na(.), 0)%>%
  mutate(QC_n.OFFC = as.numeric(QC_n.OFFC),
         QC_n.SC = as.numeric(QC_n.SC),
         QC_n.VS = as.numeric(QC_n.VS),
          QC_n.c2.c5.c20 = as.numeric(QC_n.c2.c5.c20))%>%
  mutate(timestep = as.numeric(timestep))

perWell.MOC.pred.fraction  <-perWell.MOC.pred.fraction %>%
  group_by(Well_coordinate,
           timestep)%>%
  mutate(Sum.MOC.pred = QC_n.OFFC + QC_n.SC + QC_n.VS + QC_n.c2.c5.c20,
         MOC.class.OFFC = round((QC_n.OFFC/Sum.MOC.pred)*100, digits = 3),
         MOC.class.SC =  round((QC_n.SC/Sum.MOC.pred)*100,digits = 3),
        MOC.class.VS =  round((QC_n.VS/Sum.MOC.pred)*100,digits = 3),
        MOC.class.c2.c5.c20 =  round((QC_n.c2.c5.c20/Sum.MOC.pred)*100, digits = 3))%>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         MOC.class.SC,
         MOC.class.VS,
         MOC.class.OFFC,
         MOC.class.c2.c5.c20)


LiveCells <- paste(exp.resDir,"/",unique(perWell.MOC.pred.fraction$ExpFile),"_QC_HM_MOC.percent.pdf", sep = "")

gg <-perWell.MOC.pred.fraction  %>%
  mutate(timestep = as.numeric(timestep))%>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  ggplot(aes(x = Colunm,
             y= Row))+
   geom_tile(aes(fill =  MOC.class.OFFC),width=0.7, height=0.7,color="black")+
  geom_text(aes(label = Well_coordinate), color = "black", size = 0.5)+

  coord_equal()+
  theme_bw()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+

   scale_y_discrete(limits=rev)+
   scale_x_continuous( breaks = type.of.wellplate,
                      limits=c(0,NA))+
    scale_fill_distiller(
                       type = "seq", palette = "Spectral",
                       limits = c(0,100),
                       guide ="colourbar",
                       breaks = c(0,
                               0,20,40,60,80,100),oob=scales::squish)+
   theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
  ggtitle("Heatmap of the percentage of out-of-focus cells per well")+
 theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+  
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))+
    facet_wrap_paginate(~timestep,
                      ncol =3,
                      nrow =4,
                      page =1)
n <- n_pages(gg)

pdf(LiveCells ,paper = "a4", width = 20 , height = 15 )
for(i in 1:n){
    print(gg + facet_wrap_paginate( ~ timestep, ncol =3, nrow = 4, page = i)) ## WORKS with 240 wells ## I think Error: Cannot create zero-length unit vector ("unit" subsetting) means it does not know how to plot the remianing 2 wells on the last page given ncol 4 and nrow 6
}
dev.off()


#Houskeeping
#rm( perWell.MOC.pred.fraction)

## Rate of change in out-of-focus morphology class in experiment----

  Rate.Change.OFFC.df <- perWell.MOC.pred.fraction %>%
  ungroup()%>%
  group_by(ExpFile)%>%
  filter(timestep == 0 | timestep == max(perWell.MOC.pred.fraction$timestep))%>%
  group_by(ExpFile,
           timestep)%>%
  mutate(Mean.OFFC.perentage = mean(MOC.class.OFFC))%>%
  ungroup()%>%
  select(ExpFile,
         timestep,
        Mean.OFFC.perentage )%>%
  distinct()%>%
  ungroup()%>%
  group_by(ExpFile)%>%
  mutate(Change.OFFC = (lead(Mean.OFFC.perentage) - Mean.OFFC.perentage) / last(timestep))%>%
  ungroup()%>%
  select(ExpFile,
         Change.OFFC)%>%
  drop_na()


## Plot Cell numbers over time  ----------------------------------------------------------------------------------------------------------------------------------------------------------
QC.n.Overtime <- perWell.df %>%
  select(ExpFile,
        Well_coordinate,
        timestep,
        Time_Hrs,
        QC_n.OFFC,
        QC_n.SC.VS.cells,
        QC_n.SC,
        QC_n.VS)%>%
   mutate(QC_n.OFFC = as.numeric(QC_n.OFFC))%>%
  mutate(QC_n.SC.VS.cells = as.numeric(QC_n.SC.VS.cells))%>%
  mutate(QC_n.SC = as.numeric(QC_n.SC))%>%
    mutate(QC_n.VS = as.numeric(QC_n.VS))%>%
 replace(is.na(.), 0)%>%
#  filter(timestep == "0")%>%
  mutate(
        Number.of.Out.Of.Focus.cells = as.numeric(QC_n.OFFC),
        Number.of.SingleCells = as.numeric(QC_n.SC),
         Number.of.Vsnapp = as.numeric(QC_n.VS),
       Number.of.SingleCells.Vsnaps = as.numeric(QC_n.SC.VS.cells))

QC.n.Overtime <- QC.n.Overtime %>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         Number.of.SingleCells.Vsnaps,
           Number.of.Out.Of.Focus.cells)%>%
  left_join(Conditions, by = c("Well_coordinate"))
  

LiveCells <- paste(exp.resDir,"/",unique(perWell.df$ExpFile),"_Cell.n.pdf", sep = "")
LC_plottitle <- paste ( "Total Number of SC and VS per frame")


gg <- QC.n.Overtime  %>%
  mutate(Abx.con = paste(Abx,
                         Concentration))%>%
  mutate(timestep = as.numeric(timestep))%>%
  ggplot(aes(x = timestep,
             y = Number.of.SingleCells.Vsnaps,
             group = Well_coordinate,
             label = Well_coordinate,
             colour = Abx.con))+
  geom_line( colour = "black",
             alpha = 0.5,
             size = 0.1)+

 geom_hline(aes(yintercept = 1000, colour ="1000 cells")) +
  geom_hline(aes(yintercept = 100, colour ="100 cells")) +
  ggtitle(LC_plottitle)+
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 1000))+ 
  ggtitle(LC_plottitle)+

  theme(legend.key.size = unit(0.1, "cm"),
  legend.key.width = unit(0.1,"cm"))+
  
scale_x_continuous(limits = c(min(plot.xaxis.time.scale.timestep), max(plot.xaxis.time.scale.timestep)), breaks = seq(0, max(plot.xaxis.time.scale.timestep), by = 1))+ 
   labs(x = "Image Frame [a.u]",
         y = "Number of SC and VS")+

 
  theme_bw()+
    theme(aspect.ratio = 1)+
theme(axis.text = element_text(size = 3)) +
    labs(color=NULL)+
  facet_wrap_paginate(~Abx.con,
                      ncol =1,
                      nrow =1,
                      page =1)
n <- n_pages(gg)

pdf(LiveCells ,paper = "a4", width = 20 , height = 15 )
for(i in 1:n){
    print(gg + facet_wrap_paginate(~Abx.con,
                      ncol =1,
                      nrow =1, page = i)) ## WORKS with 240 wells ## I think Error: Cannot create zero-length unit vector ("unit" subsetting) means it does not know how to plot the remianing 2 wells on the last page given ncol 4 and nrow 6
}
dev.off()

#Houskeeping
rm(gg,
   LiveCells,
   LC_plottitle)



## Plot of sc , vs and offc, overview of out-of-focus cells----------------------------------------------------------------------------------------------------------------------------------------------------------
  
  LiveCells <- paste(exp.resDir,"/",unique(perWell.df$ExpFile),"_Cell.n_SC.VS.OFFC.pdf", sep = "")
  LC_plottitle <- paste ( "Total Number of SC + VS+ OFFC  per frame")


  gg <- QC.n.Overtime  %>%
  mutate(Abx.con = paste(Abx,
                         Concentration))%>%
  mutate(timestep = as.numeric(timestep))%>%
  group_by(Well_coordinate,
           timestep)%>%
  mutate(Tot.SC.VS.OFFC = (Number.of.SingleCells.Vsnaps+
                           Number.of.Out.Of.Focus.cells  ))%>%
  ungroup()%>%
  ggplot(aes(x = timestep,
             y = Tot.SC.VS.OFFC,
             group = Well_coordinate,
             label = Well_coordinate,
             colour = Abx.con))+
  geom_line( colour = "black",
             alpha = 0.5,
             size = 0.1)+

 geom_hline(aes(yintercept = 1000, colour ="1000 cells")) +
  geom_hline(aes(yintercept = 100, colour ="100 cells")) +
  ggtitle(LC_plottitle)+
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 1000))+ 
  ggtitle(LC_plottitle)+

  theme(legend.key.size = unit(0.1, "cm"),
  legend.key.width = unit(0.1,"cm"))+
  
scale_x_continuous(limits = c(min(plot.xaxis.time.scale.timestep), max(plot.xaxis.time.scale.timestep)), breaks = seq(0, max(plot.xaxis.time.scale.timestep), by = 1))+ 
   labs(x = "Image Frame [a.u]",
         y = "Number of SC and VS")+

 
  theme_bw()+
    theme(aspect.ratio = 1)+
theme(axis.text = element_text(size = 3)) +
    labs(color=NULL)+
  facet_wrap_paginate(~Abx.con,
                      ncol =1,
                      nrow =1,
                      page =1)
n <- n_pages(gg)

pdf(LiveCells ,paper = "a4", width = 20 , height = 15 )
for(i in 1:n){
    print(gg + facet_wrap_paginate(~Abx.con,
                      ncol =1,
                      nrow =1, page = i)) ## WORKS with 240 wells ## I think Error: Cannot create zero-length unit vector ("unit" subsetting) means it does not know how to plot the remianing 2 wells on the last page given ncol 4 and nrow 6
}
dev.off()

#Houskeeping
rm(gg,
   LiveCells,
   LC_plottitle)


#Houskeeping
rm(QC.n.thrs,
   QC.n.Overtime,
   QC.N,
   QC.n,
   qc.Condition,
   QC,
   perWell.FLuncorr.df,
   OC.file.list,
   n,
   i)


# SECTION 7: Live cell fractions ----
## 7.1 Live cell fraction normalisation ----
Pop.TKC <- perWell.df %>%
  select(ExpFile,
         Condition,
         Well_coordinate,
         timestep,
         Time_Hrs,
     matches("LC_piNeg.fraction"))%>%
    select(-matches("RAW"))%>%
  drop_na()


Pop.TKC <- melt(Pop.TKC,
           id.vars = c("ExpFile",
                       "Condition",
                                  "Well_coordinate",
                       "timestep",
                       "Time_Hrs"),
                      variable.name =  "Time.Kill.Definitions")


Pop.TKC <- Pop.TKC %>%
  mutate(timestep = as.numeric(timestep),
         Time_Hrs = as.numeric(Time_Hrs),
         value = as.numeric(value))

# Loading external input data for maximum live cell fraction
maxLCF.df <- read.csv(maxLCF.df.datapath )

maxLCF.df <-maxLCF.df %>%
  mutate(Isolate = if_else(Isolate =="ATc.19979", "Iso.ATc.19979", Isolate))

# Finding killing defitnions
defintions.considered.in.analysis <- unique((Pop.TKC$Time.Kill.Definitions))
defintions.considered.in.analysis <- as.character(defintions.considered.in.analysis)

defintions.considered.in.analysis <- c("Isolate",defintions.considered.in.analysis)


# Select isolate 
maxLCF.df <- maxLCF.df %>%
  select(all_of(defintions.considered.in.analysis))

  
  rm(defintions.considered.in.analysis)
# left join such that I have a Max.LCF.of.Iso.at.tp0
  
  maxLCF.df <- maxLCF.df %>%
   mutate(Isolate = gsub("Iso.","", Isolate))
  

#Creating a dataframe with directory details
Pop.TKC <- separate(Pop.TKC, col = Condition,
                      into = c("Abx",
                               "Concentration",
                               "Isolate"),
                      sep="_")


Pop.TKC$Isolate <- gsub("Iso.", "",Pop.TKC$Isolate)

maxLCF.df$Isolate  <- gsub("Iso.", "",maxLCF.df$Isolate)

  maxLCF.df <- melt(  maxLCF.df,
                      id.vars =  c("Isolate"),
                      variable.name = "Time.Kill.Definitions",
                      value.name = "maxLCF",
                      )
#Finding the maximum live cell fraction at timepoint 0 
Pop.TKC.tp0 <- Pop.TKC %>%
  dplyr::ungroup()%>%
   left_join(maxLCF.df)%>%
  dplyr::group_by(ExpFile,
                  Isolate,
                  Time.Kill.Definitions)%>%
  mutate(Max.LCF.of.Iso.at.tp0 = maxLCF)%>%
    select(-maxLCF)
  
Pop.TKC.corr <- left_join(Pop.TKC,
                                Pop.TKC.tp0  )

#Houskeeping


Pop.TKC.corr <- Pop.TKC.corr %>%
    dplyr::group_by(ExpFile,
                   Abx,
                   Concentration,
                   Well_coordinate,
                  Isolate,
                  Time.Kill.Definitions)%>%
  mutate(LC.fraction.corr = (1/first(Max.LCF.of.Iso.at.tp0))*value)%>%
  mutate(timestep = timestep + 1)%>%
  select(-value)



Pop.TKC.corr.0th.tp <- Pop.TKC.corr %>%
  ungroup()%>%
   select(ExpFile,
                   Abx,
                   Concentration,
                   Well_coordinate,
                  Isolate,
                  Time.Kill.Definitions)%>%
    distinct()%>%
    mutate(timestep = 0,
           Time_Hrs = 0,
          LC.fraction.corr = 1)


Pop.TKC.corr <- rbind(Pop.TKC.corr,
                  Pop.TKC.corr.0th.tp)



Pop.TKC.corr <- Pop.TKC.corr %>%
  select(-Max.LCF.of.Iso.at.tp0)


#Housekeeping
rm(Pop.TKC.corr.0th.tp,
   Pop.TKC,
   Pop.TKC.tp0)



write.csv(Pop.TKC.corr, paste(exp.resDir,
                                 "/",
                                 unique(Pop.TKC.corr$ExpFile),
                                 "_Pop.PiClass_MultiBaSic.csv",
                                 sep=""),
          row.names = FALSE)

rm(QC.thrs.n,
   QC.n)

rm(
  maxLCF.df,
  Rate.Change.OFFC.df,
  Rate.of.Change.BFbck,
   QC.subset,
   perWell.MOC.pred.fraction,
  i,
  n)


# ..... POPULATION ANALYSIS ----
# Here we aim to assess the overall experimental quality regarding events (e.g. growth) and reproducibility and evaluate which isolate passes all quality control to be considered in downstream analysis. Three different killing features, i.e. the live cell fraction at comparable time points, the area under the curve (AUC) and minimum duration of killing (MDK), are aggregated and exported in a final table of results across the different time kill definitions.

# SECTION 8: Kill features ----
##  Loading killing features data (population & tracking ) 
# Here we will load the population analysis KF and also the tracking killing features if available
### Loading population Killing Features

kf.dir <- c("ASCT_Data/PerWell-KF.Results_Top2_MeanNorm")
kf.dir <- paste(genDir,
                "/",
                kf.dir,sep="")
  
#Creacting vector with the list of file names
perKF.filenames <- list.files(path = kf.dir,
                            pattern = "*_KF_longformat.csv",
                            full.names = FALSE)

#setting path to where all the csv files due to be processed are
setwd(kf.dir)


OC.file.list <- lapply(perKF.filenames,
                      read.csv)


# Convert all columns of each data frame into character
OC.file.list <- lapply(OC.file.list, function(df) {
  as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
})


KF.df.Pop <- bind_rows(OC.file.list, .id = "Exp.N")



#Houskeeping
rm(OC.file.list)

KF.df.Pop <- KF.df.Pop%>%
  mutate(Killing.Def = Killing.Features )


##  If-else statment to determine whther to also load tracking killing features ----------------------------------------------------------------------------------------------------------------------------------------------------------
 if ( tracking.availablitiy == "Yes") {
  # _KF.Trkv1_longformat.csv
   kf.pattern <- paste("*_KF.",old.or.new.tracking,"_longformat.csv",sep="")
   #Creacting vector with the list of file names
perKF.filenames <- list.files(path = kf.dir,
                            pattern =    kf.pattern,
                            full.names = FALSE)

   
 
#setting path to where all the csv files due to be processed are
setwd(kf.dir)


OC.file.list <- lapply(perKF.filenames,
                      read.csv)




# Convert all columns of each data frame into character
OC.file.list <- lapply(OC.file.list, function(df) {
  as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
})

   

#Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
KF.df.Trk  <- rapply(
  OC.file.list,
  as.character,
  how = "replace"
)

KF.df.Trk <- bind_rows(OC.file.list, .id = "Exp.N")



#Houskeeping
rm(OC.file.list)
  rm(   kf.pattern)

KF.df.Trk <- KF.df.Trk %>%
  mutate(Killing.Def = Killing.Features )


   
 } else ( print("SECTION 8: tracking killing feautres unavailable "))


# Merging data

 if ( tracking.availablitiy == "Yes") {
KF.df <- rbind(KF.df.Pop,
               KF.df.Trk)

rm(KF.df.Pop,
   KF.df.Trk)
 } else ( KF.df <-KF.df.Pop )



KF.df <- KF.df %>%
  distinct()

KF.df <- KF.df %>%
 mutate(ExpFile = zb.exp)%>%
  select( ExpFile ,
         Isolate,
         Well_coordinate,
         Abx.con,
         Killing.Features,
         Killing.Def,
         value) %>%
  distinct()



## Relabeling killing features  ----------------------------------------------------------------------------------------------------------------------------------------------------------

# The variable names I made while generating the killing features were horrible. I need to make it cleaner before downstream analysis. This will overall help condense the code because I have allot of lines that infact correct the killing features headings
label.kf.old.new <- paste(label_dir,
                          "/",
                          "ASCT_KF_Old_to_New_MultiBaSic_MultiTrk_v3.csv",sep="")
neo.labels <- read.csv(label.kf.old.new)

neo.labels <-neo.labels %>%
   mutate(Pop.Trk.new = if_else(Pop.Trk == "Pop.PiClass.",
                                   "P.Ila",
                                   if_else(Pop.Trk == "Pop.Thrs.", 
                                           "P.Thr", 
                                           if_else(Pop.Trk == paste(old.or.new.tracking,"_",sep=""),
                                                  old.or.new.tracking,
                                                  Pop.Trk))))%>%
  mutate(Killing.Features.new = paste(Pop.Trk.new,
                                      "_",
                                  Organised_Killing_Features,
                                  sep=""))%>%
  mutate(Killing.Features = paste(Pop.Trk,
                                  Current_Killing_Features,
                                  sep=""))%>%
  select(Killing.Features,
        Killing.Features.new )



KF.df <- KF.df%>%
  left_join(neo.labels, by = c("Killing.Features"))

KF.df <- KF.df %>%
  ungroup()%>%
  mutate(ExpFile = zb.exp)%>%
    select(ExpFile ,
         Isolate,
         Well_coordinate,
         Abx.con,
         Killing.Features,
         Killing.Features.new,
         value)%>%
  mutate(Killing.Features.new = gsub("_Drug_","_" ,Killing.Features.new))

# Dropping killing features per image frame,intinal persister fraction and crossing time. The will have NA as Killing feature because they were not included in renaming script

KF.df <- KF.df %>%
 group_by(Well_coordinate) %>%
  filter(!is.na(Killing.Features.new)) %>%
  ungroup()

KF.df <- KF.df %>%
  mutate(Killing.Features = Killing.Features.new)

rm(neo.labels)


KF.df <-  separate(KF.df, col = Killing.Features.new,
                      into = c("Analysis.Str",
                               "Main.KF",
                               "Sub.KF",
                               "cell.type"),
                   sep= "_")

KF.df <- KF.df %>%
  mutate(Killing.Def = paste(Analysis.Str,
                             Main.KF,
                             cell.type,
                             sep="_") )
  Main.KF.df <- KF.df
  Main.KF.df.TRACKING <- KF.df


  
# SECTION 9: Detecting outliers ----
###  Assigning replicate numbers to isolates ----

Iso.Rep.N <- Conditions %>%
  ungroup()%>%
  mutate(Isolate = if_else(Isolate == "ATc.19979", "Iso.ATc.19979", Isolate))%>%
  group_by(ExpFile,
           Isolate)%>%
  mutate(rep = 1,
         rep = cumsum(rep))


  # Killing features of OutLiers Coefficient of Variation
  KF.OL.CV.df  <- KF.df %>%
    mutate(Isolate = paste("Iso.",Isolate,sep=""))%>%
  filter(Main.KF == "LCF" & Analysis.Str == "P.Ila")%>%
  mutate(value = as.numeric(value),
         ExpID = ExpFile)%>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  ungroup()%>%
    select(ExpFile,
           ExpID,
           Abx.con,
           Isolate,
           Well_coordinate,
           Colunm,
           Killing.Features,
           Killing.Def,
           Main.KF,
           value)
  
 
 

    KF.OL.CV.df <-   KF.OL.CV.df  %>%
    left_join(Iso.Rep.N)%>%
    drop_na()
  

  # Houskeeping
  rm(Iso.Rep.N)



### Coefficient of variation of triplicates ----
  
KF.OL.CV.df <- KF.OL.CV.df 


# KF.OL.3cv.Mean.perRep sHould be changed to mean.perISo
KF.OL.3cv.Mean.perRep <- KF.OL.CV.df %>%
  ungroup()%>%
  mutate(rep = as.numeric(rep))%>%
  mutate(value = if_else(value > 1 , 1, value))%>%
  mutate(value.log10 = signif(log10(value), digits = 3))%>%
  group_by(Abx.con,
           Isolate,
           Killing.Features)%>%
  mutate(CV.triplicate = sd(value)/mean(value))%>%
  mutate(CV.triplicate.lg10 = sd(value.log10)/mean(value.log10))%>%
  ungroup()%>%
   group_by(Abx.con,
           Isolate,
           Killing.Def)%>%
  mutate(CV.triplicate.Mean.perRep = mean(CV.triplicate),
         CV.triplicate.Mean.perRep.lg10 = mean(CV.triplicate.lg10))



KF.OL.3cv.Mean.perRep <- KF.OL.3cv.Mean.perRep %>%
  ungroup()%>%
  select(-Killing.Features,
         -value,
         -Main.KF,
         -Colunm,
         -value.log10,
         -CV.triplicate,
         -CV.triplicate.lg10)%>%
  distinct()%>%
  drop_na()
  

### Coefficient of variation of duplicates LONG LOOP ----

CV.combo.res <- data.frame()

  KF.OL.CV.df <- KF.OL.CV.df %>%
  mutate(loop.var = paste(ExpFile,
                          Abx.con,
                          Isolate,sep="_"))%>%
    
  mutate(value.lg10 = signif(log10(value), digits = 3))


loop.list <-   unique(  KF.OL.CV.df$loop.var)


#x <- loop.list[1]
#x
for ( x in loop.list) {
  
  
  KF.OL.CV.df.subset <- KF.OL.CV.df %>%
  filter(loop.var == x)
  
  a <- unique((KF.OL.CV.df.subset$Well_coordinate))
  a
  b <- unique((KF.OL.CV.df.subset$Well_coordinate))
  b
  combo <- crossing(a,b)
  
  
  #Housleeping
  rm(a,
     b)
  
  
  combo <- combo %>%
    mutate(Identical.combo = if_else (a == b , "Omit", "keep"))%>%
    filter( Identical.combo == "keep")
  
  
    combo <- combo %>%
    mutate(Well_rep_combinations = paste(a,b,sep="_"))
    
    
    list.of.duplicates <- combo$Well_rep_combinations


    for ( i in     list.of.duplicates  ) {
 
      # before and after the _ 
      i.Well_rep1 <- sub("_.*", "", i)  
      i.Well_rep2 <- sub(".*_", "", i)  
  
      
      Subset.combo <-KF.OL.CV.df.subset %>%
        filter(Well_coordinate == i.Well_rep1 | Well_coordinate ==  i.Well_rep2)%>%
        mutate(Combo = paste(i))%>%
        ungroup()%>%
        group_by(Killing.Features)%>%
        mutate(CV.combo2 = sd(value)/mean(value))%>%
        mutate(CV.combo2.lg10 = sd(value.lg10)/mean(value.lg10))%>%
        select(-value,
               -value.lg10,
               -Colunm,
               )%>%
        distinct()
      
        CV.combo.res <- rbind(  CV.combo.res ,Subset.combo)
        
          rm(i,
           Subset.combo)
    
      
      
    }
    # break
    }
  
  #Houskeeping
  rm(combo,
     KF.OL.CV.df.subset,
        Subset.combo,
     x,
    
     i.Well_rep1,
     i.Well_rep2)
  

 



# Capture the end time
end_time <- Sys.time()

# Calculate the time difference in minutes
execution_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

# Print the execution time
cat("Script execution time:", execution_time, "minutes\n")


###  Creating dataframe for best CV  ----
# Finding the average of the best duplicates
  KF.OL.2cv.perRep.Mean.perRep <- CV.combo.res %>%
  ungroup()%>%
  select(ExpID,
         Abx.con,
         Killing.Def,
        Isolate,
         Combo,
         CV.combo2,
         CV.combo2.lg10)%>%
  group_by(ExpID,
         Abx.con,
         Killing.Def,
         Isolate,
         Combo)%>%
  mutate(CV.duplicate.Mean.perRep = mean(CV.combo2),
         CV.duplicate.Mean.perRep.lg10 = mean(CV.combo2.lg10))%>%
  ungroup()%>%
  select(-CV.combo2.lg10,
         -CV.combo2)%>%
  distinct()%>%
   group_by(ExpID,
         Abx.con,
         Killing.Def,
         Isolate)%>%
  mutate(CV.duplicate.Mean.perRep.BestDup = min(CV.duplicate.Mean.perRep),# Finding the best duplicate combinaiton, by looking at the combinations witht he smalles CV
         CV.duplicate.Mean.perRep.BestDup.lg10 = min( CV.duplicate.Mean.perRep.lg10))%>%
   mutate(Best.combo.Dup = if_else(CV.duplicate.Mean.perRep.BestDup == CV.duplicate.Mean.perRep , "Best-combo","Poor-combo"))%>%
   mutate(Best.combo.Dup.lg10 = if_else(CV.duplicate.Mean.perRep.BestDup.lg10 == CV.duplicate.Mean.perRep.lg10, "Best-combo","Poor-combo"))%>%
  ungroup()%>%
  distinct()


### Merging triplicate and duplicate CV ----------------------------------------------------------------------------------------------------------------------------------------------------------
KF.OL.3cv.2cv.Mean.perRep <- left_join(KF.OL.3cv.Mean.perRep,
                           KF.OL.2cv.perRep.Mean.perRep)


Check.CV <- KF.OL.CV.df
#Houskeeping
rm(KF.OL.3cv.Mean.perRep,
    KF.OL.2cv.perRep.Mean.perRep)


## Distribution of coeffient of variation: 3cv and 2cv  ----------------------------------------------------------------------------------------------------------------------------------------------------------
Outlier.df <- KF.OL.3cv.2cv.Mean.perRep %>%
  ungroup()%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Def,
         CV.triplicate.Mean.perRep,
         CV.duplicate.Mean.perRep.BestDup)%>%
  ungroup()%>%
  distinct()%>%
  drop_na() # SOLUTION HERE 20230618


Outlier.df <- melt(Outlier.df,
                                       id.vars = c("ExpFile",
                                                   "Abx.con",
                                                   "Isolate",
                                                   "Killing.Def"))
  Outlier.df <- Outlier.df %>%
  group_by(Killing.Def,
           variable)%>%
  mutate(Mean3cv = mean(value))%>%
  mutate(Mean.3stdev =3*sd(value) + Mean3cv)
  
  
  CV.duplicate.thrs <- Outlier.df %>%
  ungroup()%>%
  filter(variable == "CV.duplicate.Mean.perRep.BestDup" & Killing.Def ==  "P.Ila_LCF_scvsBaSic")%>%
  select(Mean.3stdev)%>%
  distinct()

CV.duplicate.thrs <- CV.duplicate.thrs$Mean.3stdev



CV.triplicate.thrs <-Outlier.df %>%
  ungroup()%>%
  filter(variable == "CV.triplicate.Mean.perRep" & Killing.Def ==  "P.Ila_LCF_scvsBaSic")%>%
  select(Mean.3stdev)%>%
  distinct()

CV.triplicate.thrs <- CV.triplicate.thrs$Mean.3stdev


  
library(ggprism)
  plot.path <- paste(exp.resDir,
                     "/",
                     unique(perWell.df$ExpFile),
  "_Reproducibility_CoV.pdf",sep="")


# HOWKRING HERE


pdf(  plot.path)
  Outlier.df %>%
    filter(Killing.Def == "P.Ila_LCF_scvsBaSic")%>%
  ggplot(aes(x = value,
           colour = variable,
           fill = variable))+
    
     geom_histogram(aes(y = ..density..),
             #   colour = 1,
                 alpha = 0.25) +
  geom_density(lwd = 1,  alpha = 0.25)+
      xlim(0, 0.75) +


  xlab("Coefficient of variation [a.u]")+
  ylab("Frequency [a.u]")+
  ggtitle(paste("Histogram CoV",
                zb.exp,
                unique(Outlier.df$Abx.con),sep=" "))+
    labs( caption = paste("Mean +3sd : Triplicate ", signif(CV.triplicate.thrs, digits = 4) , " | Duplicate:", signif(CV.duplicate.thrs, digits = 4) ,sep="") )+
  geom_vline(aes(xintercept = Mean3cv),
               colour = "Black")+
  geom_vline(aes(xintercept = Mean.3stdev),
               colour = "Grey")+
    theme(
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal" # Display legend items horizontally
    
    ) +

  theme(aspect.ratio =1)+

  facet_wrap(Killing.Def~variable)
dev.off()
#--


### Exporting outliers data ----------------------------------------------------------------------------------------------------------------------------------------------------------

Outlier.Data.export <- Outlier.df %>%
  filter(Killing.Def == "P.Ila_LCF_scvsBaSic")%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Def,
         variable,
         value,
         Mean3cv,
         Mean.3stdev)%>%
  distinct()

write.csv(Outlier.Data.export,file = paste(exp.resDir,
                          "/",
                          unique(Outlier.Data.export$ExpFile),
                          "_Outlier_CoV.Values.csv",
                          sep=""),
          row.names = FALSE)



Outlier.Data.export <- Outlier.df %>%
  filter(Killing.Def == "P.Ila_LCF_scvsBaSic")%>%
  select(ExpFile,
         Abx.con,
         Killing.Def,
         variable,
         Mean3cv,
         Mean.3stdev)%>%
  distinct()

write.csv(Outlier.Data.export,file = paste(exp.resDir,
                          "/",
                          unique(Outlier.Data.export$ExpFile),
                          "_Outlier_CoV.Values_Mean.Stdev.csv",
                          sep=""),
          row.names = FALSE)


rm(Outlier.Data.export)



### Flagging outliers----------------------------------------------------------------------------------------------------------------------------------------------------------
Outlier.df.intermed <- KF.OL.3cv.2cv.Mean.perRep %>%
  distinct() %>%
  filter(Killing.Def == "P.Ila_LCF_scvsBaSic")%>%
  group_by(ExpFile,
           Abx.con,
           Isolate)%>%
  mutate(Dup.Cal = (100/CV.triplicate.Mean.perRep)*CV.duplicate.Mean.perRep.BestDup)%>% # How much lower is the duplicate cv compared to the triplicate cv
  mutate(Dup.cal.test  =  if_else( Dup.Cal <= 50, "Keep-Best-Duplicate" ,"consider-triplcate"))%>%
  mutate(Outlier.Evaluation = if_else(signif(CV.triplicate.Mean.perRep, digits = 2)  > CV.triplicate.thrs  & signif(CV.duplicate.Mean.perRep.BestDup, digits = 2)  <= CV.duplicate.thrs , "Keep-Best-Duplicate" ,
                                      if_else( signif(CV.triplicate.Mean.perRep, digits= 2) <= CV.triplicate.thrs , "Keep-Triplicate", "Omit-Isolate")))%>%
  ungroup()%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Outlier.Evaluation)%>%
  distinct()

Outlier.df <- KF.OL.3cv.2cv.Mean.perRep %>%
  left_join(Outlier.df.intermed)

rm(  Outlier.df.intermed)


# Build a dataframe of all the Well coordinates to keep and the ones flagged as outliers.
list.of.outlier.evaluation <-c ("Keep-Triplicate", "Keep-Best-Duplicate" )

Outlier.wells <- data.frame()

for ( i in list.of.outlier.evaluation ) {
  print(i)
  if ( i == "Keep-Triplicate" ) {
    
     Outlier.df.subset <-   Outlier.df %>%
    filter(Outlier.Evaluation  == "Keep-Triplicate")%>%
    select(ExpID,
           Well_coordinate,
           Abx.con,
           Isolate)%>%
  distinct()%>%
    mutate(Outlier.Evaluation = "Reproducible")
     
     Outlier.wells <- rbind(    Outlier.wells  ,
                                 Outlier.df.subset)
  }
  
  
  if ( i == "Keep-Best-Duplicate" ) {
    
     Outlier.df.subset <-   Outlier.df %>%
    filter(Outlier.Evaluation  == "Keep-Best-Duplicate")%>%
       ungroup()%>%
    select(ExpID,
           Well_coordinate,
           Abx.con,
           Isolate,
           Outlier.Evaluation,
           Combo,
           Best.combo.Dup)%>%
    distinct()%>%
    filter(Best.combo.Dup == "Best-combo")%>%
    mutate(Outlier.Evaluation = "Reproducible")%>%
       select(-Outlier.Evaluation)
     
    Well_coordinate <- unlist(paste(Outlier.df.subset$Combo,collapse="_"))
    Well_coordinate <- strsplit( Well_coordinate , split = "_")
    Well_coordinate <- Well_coordinate[[1]]
    Well_coordinate <-  unique(Well_coordinate)
      
  
     Outlier.df.subset.subset <- data.frame( Well_coordinate =  Well_coordinate,
                                      Outlier.Evaluation = "Reproducible"
                                      )
       
          Outlier.df.subset <- left_join(   Outlier.df.subset,   Outlier.df.subset.subset , by = ("Well_coordinate"))
          
           Outlier.df.subset <-  Outlier.df.subset %>%
             drop_na()%>%
               select(ExpID,
           Well_coordinate,
           Abx.con,
           Isolate,
            Outlier.Evaluation)%>%
  distinct()
             
     Outlier.wells <- rbind(    Outlier.wells  ,
                                 Outlier.df.subset)
     
     #Houskepeing
     rm(Outlier.df.subset)
  }
 
}
 



### Merging time kill curve data with Outlier data  ----------------------------------------------------------------------------------------------------------------------------------------------------------
  exp.name <- unique(perWell.df$ExpFile)


  Outlier.wells <- Outlier.wells %>%
  mutate(ExpFile = exp.name)

  Pop.TKC.corr.OL1.df  <- Pop.TKC.corr %>%
  mutate(Time.Kill.Definitions = gsub("LC_piNeg.fraction.", "P.Ila_LCF_", Time.Kill.Definitions))%>% 
     mutate(Time.Kill.Definitions = gsub("P.Ila_LCF_sc.thrs", "P.Thr_LCF_sc", Time.Kill.Definitions))%>%
      mutate(Time.Kill.Definitions = gsub("P.Ila_LCF_sc.vs.thrs", "P.Thr_LCF_sc.vs", Time.Kill.Definitions))%>%
  

  mutate(Isolate = paste("Iso.",Isolate,sep=""))%>%
  mutate(Abx.con = paste(Abx,
                         Concentration,
                         sep="_"))


  Pop.TKC.corr.OL2.df <- left_join(Pop.TKC.corr.OL1.df, 
                                   Outlier.wells, by = c("ExpFile","Well_coordinate","Isolate","Abx.con"))


  Pop.TKC.corr.OL2.df <- Pop.TKC.corr.OL2.df %>% 
  replace_na(list(Outlier.Evaluation = 'Outlier'))%>%
    mutate(ExpID = ExpFile)
  
  Pop.TKC.corr.OL2.df$Outlier.Evaluation%>% replace_na('Outlier')
  
  



rm(CV.combo.res,
   Check.CV,
   KF.OL.3cv.2cv.Mean.perRep,
   KF.OL.CV.df,
   Outlier.df.subset.subset,
   Outlier.df,
   Pop.TKC.corr.OL1.df,
   i,
   i.Well_rep1,
   i.Well_rep2,
   list.of.duplicates,
   list.of.drugs,
   list.of.outlier.evaluation,
   loop.list,
   n,
   x)


## SECTION 10:  Final check ; qc outliers and sanity checks ----------------------------------------------------------------------------------------------------------------------------------------------------------
# We will find the difference between 48 and 72h mark of experiment
library(stringr)


drug <- Pop.TKC.corr.OL2.df %>%
  ungroup()%>%
  select(Abx.con)%>%
  drop_na()

drug  <- unique(drug$Abx.con)



LCF.Inc.All.timepoints.df.II <- Pop.TKC.corr.OL2.df %>%
  ungroup()%>%
 select(ExpFile,
        Well_coordinate,
         Isolate,
         Abx.con,
        Time.Kill.Definitions,
        timestep,
        Time_Hrs,
        LC.fraction.corr)%>%
  mutate(Abx.con = drug )%>%
  filter( Time.Kill.Definitions == "P.Ila_LCF_sc.vs.BaSic")%>%
  ungroup()%>%
  group_by(Well_coordinate)%>%
  filter(timestep == max(timestep))



LCF.Inc.All.timepoints.df.I <- Pop.TKC.corr.OL2.df %>%
  ungroup()%>%
 select(ExpFile,
        Well_coordinate,
         Isolate,
         Abx.con,
        Time.Kill.Definitions,
        timestep,
        Time_Hrs,
        LC.fraction.corr)%>%
    mutate(Abx.con = drug )%>%
  filter( Time.Kill.Definitions == "P.Ila_LCF_sc.vs.BaSic")%>%
  ungroup()%>%
  group_by(Well_coordinate)%>%
  filter(timestep < max(timestep))%>%
  filter(LC.fraction.corr == min(LC.fraction.corr))


LCF.Inc.df <- rbind(LCF.Inc.All.timepoints.df.I,
                    LCF.Inc.All.timepoints.df.II)

  rm(LCF.Inc.All.timepoints.df.I,
                    LCF.Inc.All.timepoints.df.II)

  LCF.Inc.df <- LCF.Inc.df %>%
  select(-timestep)%>%
  ungroup()%>%
  group_by(
           Abx.con,
           Isolate,
           Well_coordinate)%>%
  arrange(Time_Hrs)%>%
  mutate(LCF.diff = (lead(LC.fraction.corr)- LC.fraction.corr ) *100)%>%
  ungroup()%>%
   filter(!is.na(LCF.diff)) %>%
  drop_na()%>%
  mutate(LCF.increase = if_else(signif(LCF.diff, digits =3) >= 5 , "Increasing" , "Decreasing"))%>%
  select(-LC.fraction.corr)%>%
  select(-Time_Hrs)


# Reproducible time kill curbes
LCF.Inc.df.REPRODUCIBLE.allIso <- Pop.TKC.corr.OL2.df %>%
  ungroup()%>%
  select(ExpFile,
         Abx.con,
         Isolate,
        Well_coordinate,
        timestep,
        Time_Hrs,
        Time.Kill.Definitions,
         Outlier.Evaluation,
        LC.fraction.corr)%>%
    mutate(Abx.con = drug)


# Isolating all the Reproducible time-kill curves

LCF.Inc.df.REPRODUCIBLE.allIso <- dplyr::left_join(LCF.Inc.df.REPRODUCIBLE.allIso,
                                            LCF.Inc.df,
                                            by = c("ExpFile","Abx.con", "Isolate","Well_coordinate"))

LCF.Inc.df.REPRODUCIBLE <- LCF.Inc.df.REPRODUCIBLE.allIso%>%
  filter(Outlier.Evaluation == "Reproducible")


LCF.Inc.df.REPRODUCIBLE <- LCF.Inc.df.REPRODUCIBLE %>%
  mutate(Outlier_LCFin = paste(Outlier.Evaluation,
                               LCF.increase,
                               sep="_"))


# Houskeeping
rm(drug,
   LCF.Inc.df,
   LCF.Inc.df.REPRODUCIBLE)

LCF.Inc.df.REPRODUCIBLE.allIso <- LCF.Inc.df.REPRODUCIBLE.allIso %>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Well_coordinate,
         LCF.increase)%>%
  distinct()

### Event: Growth, artefact, gel-detachment, contamination, late-occuring contamination (LOC), reproduicble ----------------------------------------------------------------------------------------------------------------------------------------------------------

if ( type.of.data == "Pop.only") {
  
  QC <- perWell.df %>%
  select(ExpFile,
         Well_coordinate,
         QC_Well.image.quality,
         GT_Growth.eval.in.well,
         Contamination_eval,
         GD_eval.gelDet)%>%
  distinct()
  
  # Adding information for LOC 
  LOC.df <- LCF.Inc.df.REPRODUCIBLE.allIso %>%
      mutate(LOC = LCF.increase,
             LOC = if_else(LOC == "Increasing", "Contamination", "No-contamination"))%>%
    select(-LCF.increase)
  
  QC <- left_join(QC,
                  LOC.df)
  


#  here we consider instance 
#   When an isolate grows despite treatment, we assume it is resistant. We aim to identify isolates for growth because it does not make sense to assess killing in conditions where there is growth.
# Sometimes, two out of three replicates are flagged for growth in a few isolates. The third replicate is not flagged because it is shy of hitting the growth threshold . Here we consider when 2 out of 3 replicates are flagged for growth, and then the 3rd curve/ replicate should also be flagged for growth. E.g. In ASCT.03 Iso.159

  QC <- QC %>%
  group_by(ExpFile,
           Isolate) %>%
  mutate(Growth_Count = sum(GT_Growth.eval.in.well == "growth-detected")) %>%
  mutate(GT_Growth.eval.in.well = ifelse(Growth_Count >= 2 & GT_Growth.eval.in.well == "no-growth-detected",
                                          "growth-detected", GT_Growth.eval.in.well)) %>%
  select(-Growth_Count)
  
  QC <- QC %>%
    drop_na()%>%
    mutate(Contamination_eval2 = if_else(Contamination_eval == "No-contamination" & LOC == "Contamination", "Contamination",
                                         if_else(Contamination_eval == "Contaminaiton", "Contamination",
                                                 if_else(Contamination_eval == "No-contamination" & LOC == "No-contamination", "No-contamination", 
                                                            if_else(Contamination_eval == "Contamination" & LOC == "Contamination", "Contamination", 
                                                         "Undefined")))))%>%
    mutate(Contamination_eval = Contamination_eval2)%>%
    select(-Contamination_eval2,
           -LOC)
  
  QC <- QC %>%
    mutate(QC_Well.image.quality = if_else(QC_Well.image.quality  == "good", "Keep", "Omit"),
           GT_Growth.eval.in.well = if_else(GT_Growth.eval.in.well == "no-growth-detected", "Keep", "Omit"),
          Contamination_eval = if_else(Contamination_eval == "No-contamination", "Keep", "Omit"),
          GD_eval.gelDet = if_else(GD_eval.gelDet == "No-det", "Keep", "Omit"))
  
    QC <- QC %>%
      mutate(Fluorescence.Artefact = QC_Well.image.quality,
             Growth =GT_Growth.eval.in.well,
             Gel.detachment = GD_eval.gelDet,
             Contamination = Contamination_eval)%>%
      select(-QC_Well.image.quality ,
             -GT_Growth.eval.in.well,
             -Contamination_eval,
             -GD_eval.gelDet)

    QC <- QC %>%
      ungroup()%>%
      group_by(Well_coordinate)%>%
   #   mutate(Overall.assessment = if_else(Growth == "Keep" & Fluorescence.Artefact == "Keep" & Gel.detachment == "Keep" & Contamination == "Keep", 0,1))%>%
         mutate(Overall.assessment = if_else(Growth == "Keep"  & Gel.detachment == "Keep" & Contamination == "Keep", 0,1))%>%
      mutate(Overall.assessment = if_else(Overall.assessment == 0 , "Keep", "Omit" ))

    QC <- melt(QC,
           id.vars = c("ExpFile",
                                  "Well_coordinate",
                       "Abx.con",
                       "Isolate"),
           
                      variable.name =  "QC.readouts")
    
    
    qc.Condition <- Conditions %>%
  mutate(ExpFile = sub("_.*","", ExpFile ))
    
    
        QC <- QC %>%
          left_join(qc.Condition)
        
        
        QC.n <- perWell.df %>%
  select(ExpFile,
        Well_coordinate,
        timestep,
        Time_Hrs,
        QC_n.OFFC,
        QC_n.SC.VS.cells,
        QC_n.SC,
        QC_n.VS)%>%
  mutate(QC_n.OFFC = as.numeric(QC_n.OFFC))%>%
  mutate(QC_n.SC.VS.cells = as.numeric(QC_n.SC.VS.cells))%>%
  mutate(QC_n.SC = as.numeric(QC_n.SC))%>%
    mutate(QC_n.VS = as.numeric(QC_n.VS))%>%
   replace(is.na(.), 0)%>%
  mutate(
        Number.of.Out.Of.Focus.cells = as.numeric(QC_n.OFFC),
        Number.of.SingleCells = as.numeric(QC_n.SC),
         Number.of.Vsnapp = as.numeric(QC_n.VS),
       Number.of.SingleCells.Vsnaps = as.numeric(QC_n.SC.VS.cells))
        
        
        QC.n <- QC.n %>%
  filter(timestep == "0")%>%
  select(-Time_Hrs)%>%
  select(-QC_n.SC.VS.cells,
         -QC_n.OFFC,
         -QC_n.SC,
         -QC_n.VS)

      
QC <- QC %>%
  select(-Abx,
         -Concentration)
        
        # Adding info of population analysis i.e cell numbers whether we have at least 1k cells or not
        
        QC.n.thrs <- QC.n %>%
          select(ExpFile,
                 Well_coordinate,
                 Number.of.SingleCells.Vsnaps)%>%
          mutate(value = if_else(Number.of.SingleCells.Vsnaps > 1000, "Keep", "Omit"))%>%
          left_join(qc.Condition)%>%
          mutate(Isolate = if_else(Isolate == "ATc.19979" ,"Iso.ATc.19979",Isolate))%>%
          mutate(QC.readouts = "1000 Single Cell and Vsnapp population analysis")%>%
          select(-Number.of.SingleCells.Vsnaps)%>%
          mutate(Abx.con = paste(Abx,
                                 Concentration,
                                 sep="_"))%>%
          ungroup()%>%
          select(-Abx,
                 -Concentration)

      QC <- rbind(QC,QC.n.thrs )
      
      QC.subset <- pivot_wider(QC, names_from = QC.readouts, values_from = value)
        
        QC.subset <-       QC.subset%>%
        ungroup()%>%
        group_by(Well_coordinate)%>%
        mutate(Experiment.quality = if_else(Overall.assessment == "Keep" & `1000 Single Cell and Vsnapp population analysis` == "Keep", "Keep","Omit" ))%>%
          select(-Overall.assessment)

            QC.subset <-       QC.subset%>%
              drop_na()

              QC.subset<- QC.subset %>%
                mutate(Fluorescence.Artefact  = as.character(Fluorescence.Artefact),
                       Growth = as.character(Growth),
                       Gel.detachment = as.character(Gel.detachment),
                       Contamination = as.character(Contamination),
                       `1000 Single Cell and Vsnapp population analysis`= as.character(`1000 Single Cell and Vsnapp population analysis`),
                       Experiment.quality = as.character(Experiment.quality))


       QC<- melt(QC.subset,
           id.vars = c("ExpFile",
                                  "Well_coordinate",
                       "Abx.con",
                      "Isolate"),
           
                      variable.name =  "QC.readouts")
  
  
    QC <- QC %>%
      mutate(value = factor(value, levels = c("Keep","Omit")))
        
        QC <- QC %>%
          mutate(Isolate = sub("Iso.","",Isolate))
        
        
        gg.experiment.title <- paste(unique(QC$ExpFile),
                                     " ",
                                     unique(QC$Abx),
                                     ".",
                                     unique(QC$Concentration),
                                     sep="") 
          
          
  Flagged.Analysis <- QC
          
pdf(paste(exp.resDir,"/",unique(perWell.df$ExpFile),"_QC.readouts_LOC.pdf",
          sep=""))
  print(
  QC %>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  mutate(Exp.Colunm = paste(ExpFile,
                            Colunm,
                            sep="_"))%>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  mutate(Exp.Well = paste(ExpFile,
                          Well_coordinate,
                          sep="_"))%>%
  ggplot(aes(x = Colunm,
             y= Row,
             label = Well_coordinate))+
  geom_tile(aes(fill = value),width=0.7, height=0.7)+
  scale_fill_manual(values = c("SteelBlue","coral2"))+
  geom_text(aes(label = Well_coordinate), color = "black", size = 0.3, nudge_y = 0.1)+
  geom_text(aes(label = paste(Isolate,sep=".")), color = "white", size = 0.2, nudge_y = -0.2 )+
  coord_equal()+
  theme_bw()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+
   scale_y_discrete(limits=rev)+
   scale_x_continuous( breaks = seq(type.of.wellplate),
                      limits=c(0,NA))+
   labs(title = gg.experiment.title)+
      theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
 theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+  
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))+
   
   facet_wrap(~QC.readouts,
               ncol = 2,
              nrow = 4)
)
dev.off()

rm(gg.experiment.title)
  
}


if (type.of.data == "Pop_&_Trk") {
   QC <- perWell.df %>%
  select(ExpFile,
         Well_coordinate,
         QC_Well.image.quality,
         GT_Growth.eval.in.well,
         Contamination_eval,
         GD_eval.gelDet)%>%
  distinct()
  
  QC <- QC %>%
    mutate(QC_Well.image.quality = if_else(QC_Well.image.quality  == "good", "Keep", "Omit"),
           GT_Growth.eval.in.well = if_else(GT_Growth.eval.in.well == "no-growth-detected", "Keep", "Omit"),
          Contamination_eval = if_else(Contamination_eval == "No-contamination", "Keep", "Omit"),
          GD_eval.gelDet = if_else(GD_eval.gelDet == "No-det", "Keep", "Omit"))
  
    QC <- QC %>%
      mutate(Fluorescence.Artefact = QC_Well.image.quality,
             Growth =GT_Growth.eval.in.well,
             Gel.detachment = GD_eval.gelDet,
             Contamination = Contamination_eval)%>%
      select(-QC_Well.image.quality ,
             -GT_Growth.eval.in.well,
             -Contamination_eval,
             -GD_eval.gelDet)

    QC <- QC %>%
      ungroup()%>%
      group_by(Well_coordinate)%>%
      mutate(Overall.assessment = if_else(Growth == "Keep" & Gel.detachment == "Keep" & Contamination == "Keep", 0,1))%>%
      mutate(Overall.assessment = if_else(Overall.assessment == 0 , "Keep", "Omit" ))

    QC <- melt(QC,
           id.vars = c("ExpFile",
                                  "Well_coordinate"),
           
                      variable.name =  "QC.readouts")
    
    
    qc.Condition <- Conditions %>%
  mutate(ExpFile = sub("_.*","", ExpFile ))
    
    
        QC <- QC %>%
          left_join(qc.Condition)
        
        
        # Adding info of population analysis i.e cell numbers whether we have at least 1k cells or not
        QC.n.thrs <- QC.n %>%
          select(ExpFile,
                 Well_coordinate,
                 Number.of.SingleCells.Vsnaps)%>%
          mutate(value = if_else(Number.of.SingleCells.Vsnaps > 1000, "Keep", "Omit"))%>%
          left_join(qc.Condition)%>%
          mutate(QC.readouts = "1000 Single Cell and Vsnapp population analysis")%>%
          select(-Number.of.SingleCells.Vsnaps)
        
      QC <- rbind(QC,QC.n.thrs )
      
      QC.subset <- pivot_wider(QC, names_from = QC.readouts, values_from = value)
        
        QC.subset <-       QC.subset%>%
        ungroup()%>%
        group_by(Well_coordinate)%>%
        mutate(Experiment.quality = if_else(Overall.assessment == "Keep" & `1000 Single Cell and Vsnapp population analysis` == "Keep", "Keep","Omit" ))%>%
          select(-Overall.assessment)



       QC<- melt(QC.subset,
           id.vars = c("ExpFile",
                                  "Well_coordinate",
                      "Abx",
                      "Concentration",
                      "Isolate"),
           
                      variable.name =  "QC.readouts")
  
  
    QC <- QC %>%
      mutate(value = factor(value, levels = c("Keep","Omit")))
        
        QC <- QC %>%
          mutate(Isolate = sub("Iso.","",Isolate))
        
        
        gg.experiment.title <- paste(unique(QC$ExpFile),
                                     " ",
                                     unique(QC$Abx),
                                     ".",
                                     unique(QC$Concentration),
                                     sep="") 
          
          
  Flagged.Analysis <- QC
          
pdf(paste(exp.resDir,"/",unique(perWell.df$ExpFile),"_QC.readouts.pdf",
          sep=""))
print(
  QC %>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  mutate(Exp.Colunm = paste(ExpFile,
                            Colunm,
                            sep="_"))%>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  mutate(Exp.Well = paste(ExpFile,
                          Well_coordinate,
                          sep="_"))%>%
  ggplot(aes(x = Colunm,
             y= Row,
             label = Well_coordinate))+
  geom_tile(aes(fill = value),width=0.7, height=0.7)+
  scale_fill_manual(values = c("SteelBlue","coral2"))+
  geom_text(aes(label = Well_coordinate), color = "black", size = 0.3, nudge_y = 0.1)+
  geom_text(aes(label = paste(Isolate,sep=".")), color = "white", size = 0.2, nudge_y = -0.2 )+
  coord_equal()+
  theme_bw()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+
   scale_y_discrete(limits=rev)+
   scale_x_continuous( breaks = seq(type.of.wellplate),
                      limits=c(0,NA))+
   labs(title = gg.experiment.title)+
      theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
 theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+  
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))+
   
   facet_wrap(~QC.readouts,
               ncol = 2,
              nrow = 4)
)
dev.off()

rm(gg.experiment.title)

}




### Experiment quality control:  outliers and  sanity checks ---------------------------------------------------------------------------------------------------------------------------------------------------------

# Outlier detection results
experiment.drug <- Pop.TKC.corr.OL2.df %>%
  ungroup()%>%
  select(Abx.con)%>%
  drop_na()

experiment.drug <- unique(experiment.drug$Abx.con)

  
FC.df <- Pop.TKC.corr.OL2.df %>%
  ungroup()%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Well_coordinate,
         Outlier.Evaluation)%>%
  distinct()%>%
  mutate(Abx.con = experiment.drug)


# -- Experiemnt quality contorl i.e growth, contamination ,artefacts
Flagged.Analysis.overall <- Flagged.Analysis %>%
  mutate(Isolate = paste("Iso.",Isolate,sep=""))%>%
  filter(QC.readouts == "Experiment.quality" )%>%
  mutate(Experiment.Quality.control = value)%>%
   select(
         -value,
         -QC.readouts)


FC.df <- FC.df %>%
  left_join(Flagged.Analysis.overall)





### Loading sanity check data  ----------------------------------------------------------------------------------------------------------------------------------------------------------


#------------------------------------- Loading Killing features
kfsc.dir <- c("ASCT_Data/PerWell-KFsc.Results_Top2_MeanNorm")
kfsc.dir <- paste(genDir,
                "/",
                kfsc.dir,sep="")
  
#Creacting vector with the list of file names
perKFsc.filenames <- list.files(path = kfsc.dir,
                            pattern = "*_KF_sanity.check_longformat.csv",
                            full.names = FALSE)

#setting path to where all the csv files due to be processed are
setwd(kfsc.dir)


OC.file.list <- lapply(perKFsc.filenames,
                      read.csv)


# Convert all columns of each data frame into character
OC.file.list <- lapply(OC.file.list, function(df) {
  as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
})


KF.sc.df <- bind_rows(OC.file.list, .id = "Exp.N")

KF.sc.df <- KF.sc.df %>%
  select(
      #   -X,
         -Experiment.Type,
         -Exp.N)





### Load as well for tracking  Sanity Check  ----------------------------------------------------------------------------------------------------------------------------------------------------------
if ( tracking.availablitiy == "Yes") {
  
   kf.pattern <- paste("*_KF.",old.or.new.tracking,"_sanity.check_longformat.csv",sep="")
   
# ----- Load  one of each set of tracking data
perWell.filenames <- list.files(path = kfsc.dir,
                            pattern = kf.pattern,
                            full.names = FALSE)




#setting path to where all the csv files due to be processed are
setwd(kfsc.dir)


OC.file.list <- lapply(perWell.filenames,
                      read.csv)


#Naming the elements of my list with the vector list of my file names
Exp.Names <- gsub(".csv",
                  "",
                  perWell.filenames)

Exp.Names <- gsub(paste("_",
                        old.or.new.tracking,
                        "_",
                        "KF_sanity.check_longformat",
                        sep=""),
                  "",
                  Exp.Names)


names(OC.file.list)<- Exp.Names

#Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
perWell.Trk.SC.df <- rapply(
  OC.file.list,
  as.character,
  how = "replace"
)


perWell.Trk.SC.df <- bind_rows(perWell.Trk.SC.df , .id = "Filename")

rm(OC.file.list,
   Exp.Names,
   perWell.filenames)

perWell.Trk.SC.df<- perWell.Trk.SC.df%>%
  select(-Filename,
   #      -X,
         -Experiment.Type,
         -Exp.N)


 }


## Merging pooulation and tracking  ----------------------------------------------------------------------------------------------------------------------------------------------------------

if ( tracking.availablitiy == "Yes") {
  
KF.sc.df <- rbind(KF.sc.df,
                       perWell.Trk.SC.df)

KF.sc.df<- KF.sc.df %>%
  distinct()

rm( perWell.Trk.SC.df,
    OC.file.list)

} 



KF.sc.df <- KF.sc.df %>%
  group_by(Isolate) %>%
  filter(Killing.Features %in% duplicated(Killing.Features) | MDK.time.increasing == "yes")%>%
    filter(Killing.Features %in% duplicated(Killing.Features) | AUCrealtime.always.increasing == "yes")%>%
        filter(Killing.Features %in% duplicated(Killing.Features) |   AUCrealtimeLOG.always.increasing == "yes")%>%
ungroup()%>%
  distinct()




### Adjusting sanity check data labels ----------------------------------------------------------------------------------------------------------------------------------------------------------
# The variable names I made while generating the killing features were horrible. I need to make it cleaner before downstream analysis. This will overall help condense the code because I have alot of lines that in fact rename the killing features headings

label.kf.old.new <- paste(label_dir,
                          "/",
                          "ASCT_KFsc_relabel_MultiBaSic_MultiTrk.csv",sep="")
neo.labels <- read.csv(label.kf.old.new)



KF.sc.df <- KF.sc.df %>%
  mutate(ExpFile = experiment.id,
         Abx.con = experiment.drug,
         Isolate = paste("Iso.",
                         Isolate,
                         sep=""))%>%
    mutate(Killing.Def = Killing.Features)%>%
  left_join(neo.labels)%>%
  
mutate(Killing.Def = Killing.Def.new)%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Well_coordinate,
         Killing.Def,
         MDK.time.increasing,
         AUCrealtime.always.increasing,
             AUCrealtimeLOG.always.increasing,
         SC.LCF.always.decreasing)
  

KF.sc.df.assesmment <- KF.sc.df %>%
  ungroup()%>%
  group_by(Isolate,
           Well_coordinate,
           Killing.Def)%>%
  mutate(Sanity.Check = if_else(MDK.time.increasing == "yes" & AUCrealtime.always.increasing == "yes", "passed", "failed" ))%>%

 group_by(ExpFile,
          Isolate,
           Well_coordinate) %>%
  select(ExpFile,
         Well_coordinate,
         Isolate,
         Killing.Def,
         Sanity.Check)%>%
  ungroup()%>%
  distinct()
  
FC.complete <- left_join(KF.sc.df.assesmment,
                         FC.df)

#Houskeeping
rm(FC.df,
   KF.sc.df.assesmment)


# Rearrange order of variables,
FC.complete <- FC.complete %>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Well_coordinate,
         Killing.Def,
         Experiment.Quality.control,
         Outlier.Evaluation,
         Sanity.Check)
#---
# Add element for Increasing live cell fraction 

FC.complete <- left_join(FC.complete,
                         LCF.Inc.df.REPRODUCIBLE.allIso)


#-- 
 FC.complete.AS <- FC.complete %>%
   mutate(Overall.Assessment = if_else(Experiment.Quality.control == "Keep" & Outlier.Evaluation == "Reproducible" & Sanity.Check == "passed"  & LCF.increase == "Decreasing" , "Keep", "Omit"))%>%
   select(-Killing.Def)%>%
   distinct()
 

###  Overall Plot label ----------------------------------------------------------------------------------------------------------------------------------------------------------

# Here we aim to comprehensively label each well coordinate with regard to the overall analysis label
Plot.label.df <- QC


Plot.label.df <- pivot_wider( data = Plot.label.df,
                        names_from = QC.readouts,
                        values_from = value)


Plot.label.df <- Plot.label.df %>%
  select(-Experiment.quality)%>%
  select(-Fluorescence.Artefact)%>%
  mutate(Growth = if_else(Growth == "Keep", "", "Growth"),
         Gel.detachment = if_else(Gel.detachment == "Keep", "", "GD"),
         `1000 Single Cell and Vsnapp population analysis` = if_else(`1000 Single Cell and Vsnapp population analysis` == "Keep", "", "<1K"),
         Contamination = if_else(Contamination == "Keep", "", "Contam"))%>%
  mutate(Isolate= paste("Iso.",
                        Isolate,
                        sep=""))

Plot.label.df.subset <- FC.complete %>%
  select(ExpFile, 
         Abx.con,
         Well_coordinate,
         Isolate,
         Outlier.Evaluation) %>%
  mutate(Outlier.Evaluation =  if_else(Outlier.Evaluation == "Reproducible", "","OT"))


Plot.label.df <- Plot.label.df %>%
  left_join(Plot.label.df.subset)%>%
  mutate(Flag.label = paste(Growth,
                            Gel.detachment,
                             `1000 Single Cell and Vsnapp population analysis`,
                             Contamination ,
                           Outlier.Evaluation,
                            sep=" "))%>%
  select(ExpFile,
         Abx.con,
         Well_coordinate,
         Isolate,
         Flag.label)%>%
  distinct()

Plot.label.df <- Plot.label.df %>%
 mutate(Flag.label = str_replace_all(Flag.label, "\\s+", " "))


#
rm(Plot.label.df.subset )



## SECTION 11: KEY FIG  Visualising All time kill curves flagged and not-flagged **  ----------------------------------------------------------------------------------------------------------------------------------------------------------
All.iso.df <- Pop.TKC.corr.OL2.df %>%
  select(ExpFile,
         Abx,
         Concentration,
         Isolate,
         Well_coordinate,
         timestep,
         Time_Hrs,
         Time.Kill.Definitions,
         LC.fraction.corr)

All.iso.df <- left_join(All.iso.df ,
                        Plot.label.df )

All.iso.df <- left_join(All.iso.df ,
                        FC.complete.AS,
                        by = join_by("ExpFile", "Isolate", "Well_coordinate", "Abx.con"))


  write.csv(All.iso.df,
          paste(exp.resDir,
                "/",
                experiment.id,
                "_LCF_flags.csv",
                sep=""))
    
  AS.All.iso.df <- All.iso.df %>%
        ungroup()%>%
    select(ExpFile,
           Abx.con,
           Isolate,
           Well_coordinate,
           timestep,
           Time_Hrs,
           Flag.label,
           Outlier.Evaluation,
           Overall.Assessment)%>%
    distinct()%>%
    left_join(ZB.df)
  
  
  write.csv(AS.All.iso.df,
          paste(exp.resDir,
                "/",
                experiment.id,
                "_AS_MD_flags.csv",
                sep=""),
          row.names = FALSE)
          
          


## ## 11.1.1 Plot  Time kill curves with flags log scale ----------------------------------------------------------------------------------------------------------------------------------------------------------
# PLOTTING TIME KILL CURVES
LiveCells <- paste(exp.resDir,"/",exp.name,"_Pop_TimeKillCurves_Experiment_Flags.pdf", sep = "")
#LC_plottitle <- paste ( "Time kill curves")
library(ggrepel)
gg <- All.iso.df  %>%
   mutate(Abx.con = paste(Abx,
                          Concentration,
                          sep="_"))%>%
 # mutate(Outlier_LCFin = as.factor(Outlier_LCFin))%>%
  mutate(Overall.Assessment = as.factor(Overall.Assessment))%>%
  filter(Time.Kill.Definitions =="P.Ila_LCF_sc.vs.BaSic")%>%
  mutate(Time.Kill.Definitions = gsub("LC_piNeg.fraction.","Pop.",Time.Kill.Definitions))%>%
  ggplot(aes(x = Time_Hrs,
             y = log10(LC.fraction.corr),
             group = Well_coordinate,
             label = Well_coordinate,
             colour = Overall.Assessment))+
  geom_line( 
             alpha = 0.5,
             size = 0.2)+
geom_label_repel(data = . %>% group_by(Well_coordinate) %>% filter(timestep == 3) %>% 
                      distinct(),
                  aes(label = Flag.label),
                 fill = "transparent",
                   label.size = NA, label.padding = unit(0, "lines"), label.color = NA,
                 size = 0.2,
                  nudge_y = -0.3,
                   segment.color = "transparent", point.color = NA, arrow = arrow(length = unit(0.03, "npc"))) +
  theme_bw()+
   scale_y_continuous(

                    breaks = c(0, -0.999,  -1.999),
                    limits = c(-3,0),
                     labels = c(100, 10, 1))+
    scale_colour_manual(values=c(Keep="steelblue",
                                 Omit="red"))+
    geom_text(data = .  %>% group_by(Well_coordinate) %>% top_n(n = 1, wt = Time_Hrs),
    aes(label = Well_coordinate), size = 0.2)+
    theme(axis.text = element_text(size = 3))     +
    theme(strip.text.x = element_text(size = 3))+
    theme(legend.key.size = unit(0.1, "cm"),
    legend.key.width = unit(0.1,"cm"))+
    scale_x_continuous( breaks = plot.xaxis.time.scale)+
    labs(title = paste(experiment.id,
                       unique(All.iso.df$Abx.con),
                           sep=" "),

       x = "Time [Hrs]",
       y = "Percentage of cells alive [%]",
       color = "Condition evaluation",
       subtitle = "Contam = Contamination, Artefacts = AF, Gel-detachment = GD, Less 1K sc & vs = <1K, Outlier = OT")+
    geom_text(x=48, y=1,
            size = 1)+ # Plot label
  theme(legend.position="top")+
   theme(aspect.ratio = 1)+
   theme(plot.subtitle = element_text(size = 5))+

  facet_wrap_paginate(Abx.con~Isolate,
                      ncol =12,
                      nrow =12,
                      page =1)
n <- n_pages(gg)

pdf(LiveCells ,paper = "a4", width = 20 , height = 15 )
for(i in 1:n){
    print(gg + facet_wrap_paginate(Abx.con~Isolate,
                      ncol =12,
                      nrow =12, page = i))
}
dev.off()

#Houskeeping
rm(gg,
   LiveCells)


### 11.1.1 Plot  Time kill curves with flags linear ---------------------------------------------------------------------------------------------------------------------------------------------------------

# PLOTTING TIME KILL CURVES
LiveCells <- paste(exp.resDir,"/",exp.name,"_Pop_TimeKillCurves_Experiment_Flags_LINEAR.pdf", sep = "")
library(ggrepel)
gg <- All.iso.df  %>%
   mutate(Abx.con = paste(Abx,
                          Concentration,
                          sep="_"))%>%
  mutate(Overall.Assessment = as.factor(Overall.Assessment))%>%
  filter(Time.Kill.Definitions =="P.Ila_LCF_sc.vs.BaSic")%>%
  mutate(Time.Kill.Definitions = gsub("LC_piNeg.fraction.","Pop.",Time.Kill.Definitions))%>%
  ggplot(aes(x = Time_Hrs,
             y = LC.fraction.corr,
             group = Well_coordinate,
             label = Well_coordinate,
             colour = Overall.Assessment))+
  geom_line( 
             alpha = 0.5,
             size = 0.2)+
geom_label_repel(data = . %>% group_by(Well_coordinate) %>% filter(timestep == 3) %>% 
                      distinct(),
                  aes(label = Flag.label),
                 fill = "transparent",
                   label.size = NA, label.padding = unit(0, "lines"), label.color = NA,
                 size = 0.2,
                  nudge_y = -0.3,
                   segment.color = "transparent", point.color = NA, arrow = arrow(length = unit(0.03, "npc"))) +
theme_bw()+
 scale_y_continuous(

                    breaks = c(0, 0.2,0.4,0.6,0.8,1),
                    limits = c(0,1),
                     labels = c(0, 0.2,0.4,0.6,0.8, 1))+
    scale_colour_manual(values=c(Keep="steelblue",
                                 Omit="red"))+
    geom_text(data = .  %>% group_by(Well_coordinate) %>% top_n(n = 1, wt = Time_Hrs),
    aes(label = Well_coordinate), size = 0.2)+ 
    theme(axis.text = element_text(size = 3))     +
    theme(strip.text.x = element_text(size = 3))+
    theme(legend.key.size = unit(0.1, "cm"),
    legend.key.width = unit(0.1,"cm"))+
    scale_x_continuous( breaks = plot.xaxis.time.scale)+
    labs(title = paste(experiment.id,
                       unique(All.iso.df$Abx.con),
                           sep=" "),

       x = "Time [Hrs]",
       y = "Percentage of cells alive [%]",
       color = "Condition evaluation",
       subtitle = "Contam = Contamination, Artefacts = AF, Gel-detachment = GD, Less 1K sc & vs = <1K, Outlier = OT")+
    geom_text(x=48, y=1,
            size = 1)+ # Plot label
  theme(legend.position="top")+
   theme(aspect.ratio = 1)+
   theme(plot.subtitle = element_text(size = 5))+

  facet_wrap_paginate(Abx.con~Isolate,
                      ncol =12,
                      nrow =12,
                      page =1)
n <- n_pages(gg)

pdf(LiveCells ,paper = "a4", width = 20 , height = 15 )
for(i in 1:n){
    print(gg + facet_wrap_paginate(Abx.con~Isolate,
                      ncol =12,
                      nrow =12, page = i)) 
}
dev.off()

#Houskeeping
rm(gg,
   LiveCells)


## 11.2. KEY STEP exporting intermediary KF table and flags----
#### Excluding instances of two different events  ----
# Preparing Main. Killing dataframe for left_join
experiment.ID <- unique(Pop.TKC.corr$ExpFile)


# REMOVE LCF in the middle here 
  Main.KF.df.To.merge <- Main.KF.df %>%
  ungroup()%>%
  mutate(ExpFile = experiment.ID )%>%
  mutate(Isolate = paste("Iso.",
                         Isolate,
                         sep=""))%>%
  mutate(Merge.def = paste(Analysis.Str,
                           cell.type,
                           sep="_"))%>%
  select(ExpFile,
         Abx.con,
         Well_coordinate,
         Isolate,
         Merge.def,
         Main.KF,
         Killing.Def,
         Killing.Features,
         value)
  

###### 11.2.1 Instance of different events in across two replicates ---------------------------------------------------------------------------------------------------------------------------------------------------------
  FC.complete.intermediary <- FC.complete %>%
    mutate(Merge.def = gsub("_LCF_","_", Killing.Def))%>%
    mutate(Merge.def = gsub("P.Thrs", "P.Thr", Merge.def))%>%
    ungroup()%>%
    select(ExpFile,
           Abx.con,
           Isolate,
           Well_coordinate,
           Experiment.Quality.control,
           Outlier.Evaluation,
           Merge.def)%>%
    distinct()%>%
    mutate(Multiple.Events.intermediary = if_else(Experiment.Quality.control == "Keep" & Outlier.Evaluation == "Reproducible", 0 , 1))%>%
    group_by(ExpFile,
             Isolate,
             Merge.def)%>%
    mutate(Multiple.Events = sum(Multiple.Events.intermediary),
           Multiple.Events.Check = if_else(Multiple.Events > 1, "Omit", "Keep"))
  
  
  FC.complete <- FC.complete %>%
    mutate(Merge.def = gsub("_LCF_","_", Killing.Def))%>%
    mutate(Merge.def = gsub("P.Thrs", "P.Thr", Merge.def))%>%
    select(-Killing.Def)
  
  
  FC.complete <- left_join(  FC.complete ,
                             FC.complete.intermediary  )
  
    
       FC.complete <-FC.complete %>%
         mutate(Multiple.Events = Multiple.Events.Check)%>%
         select(-Multiple.Events.intermediary,
                - Multiple.Events.Check)
  
filename <- paste(exp.resDir,
                    "/",
                  unique(FC.complete$ExpFile),
                    "_",
                    experiment.drug,
                    "_Event(multi)_SanityCheck_Outlier.csv",
                    sep="")


  # Saving ASCT Clinical isoalte data
    write.csv(FC.complete.intermediary,
               filename,
              row.names=FALSE)
    
    
    rm(filename)
  FC.complete <- left_join(  FC.complete,
                        Main.KF.df.To.merge )

# Exporting Raw KF table  ----
# Exporting Killing feature population results ---- 
  filename.rawKF <- paste(exp.resDir,
                    "/",
                    unique(KF.df$ExpFile),
                    "_",
                    experiment.drug,
                    "_PoprawKF.csv",
                    sep="")
  
  kf.export <- KF.df %>%
    filter( Main.KF != "MDK")
  
  # Saving ASCT Clinical isoalte data
  write.csv(  kf.export,
            filename.rawKF,
            row.names=FALSE)
  
  rm(  kf.export,
       filename.rawKF )

##### Applying Exclusion ----------------------------------------------------------------------------------------------------------------------------------------------------------
#Houskeeping
  rm(experiment.ID,
   Main.KF.df.To.merge)

All.Iso.Analysis <- FC.complete %>%
  ungroup()%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Def,
         Main.KF,
         Killing.Features)%>%
  distinct()
  
All.Iso.Analysis.intermediary.CHECK <-All.Iso.Analysis

Average.Result.perIso.PASSED.all.Criteria <- FC.complete %>%
  filter(Multiple.Events == "Keep")%>%
  filter(Experiment.Quality.control == "Keep" ,
         Outlier.Evaluation == "Reproducible" ,
         Sanity.Check == "passed") %>%
  mutate( value = if_else(value == ">72", "72",value))%>%
  mutate(value = as.numeric(value))%>%
  ungroup()%>%
    group_by(ExpFile,
             Isolate,
             Abx.con,
             Main.KF,
             Killing.Def,
             Killing.Features) %>%
  summarise(Avg_Killing_Features = mean(value))%>%
  ungroup()%>%
  distinct()


Average.Result.perIso.PASSED.all.Criteria.INTERMEDIARY <- FC.complete %>%
  filter(Multiple.Events == "Keep")%>%
  filter(Experiment.Quality.control == "Keep" ,
         Outlier.Evaluation == "Reproducible" ,
         Sanity.Check == "passed") %>%
  mutate( value = if_else(value == ">72", "72",value))%>%
  mutate(value = as.numeric(value))%>%
  ungroup()

filename <- paste(exp.resDir,
                    "/",
                  unique(FC.complete$ExpFile),
                    "_",
                    experiment.drug,
                    "_KFintermediary.csv",
                    sep="")


  # Saving ASCT Clinical isoalte data
    write.csv(Average.Result.perIso.PASSED.all.Criteria.INTERMEDIARY,
               filename,
              row.names=FALSE)
    
    rm(Average.Result.perIso.PASSED.all.Criteria.INTERMEDIARY)
    



All.Iso.Analysis<- left_join(All.Iso.Analysis,
                              Average.Result.perIso.PASSED.all.Criteria)


#HOuskeeping
rm(FC.complete.AS)
 
# ---


All.iso <- All.Iso.Analysis  %>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Def,
         Main.KF,
         Killing.Features)

#---
All.Iso.Analysis  <- All.Iso.Analysis  %>%
  ungroup()%>%
  drop_na()%>%
  ungroup()%>%
  mutate(Killing.Features = as.factor(Killing.Features))%>%
  mutate(Avg_Killing_Features = Avg_Killing_Features)%>%
  mutate(Avg_Killing_Featureslg10 = log10(Avg_Killing_Features))%>%
  group_by(ExpFile,
           Abx.con,
           Killing.Features)%>%
  mutate(Min.Per.Feature.lg10 = min(Avg_Killing_Featureslg10))%>%
  mutate(Max.Per.Feature.lg10 = max(Avg_Killing_Featureslg10))%>%
  mutate(Norm_Avg_Killing_Features.lg10 = (Avg_Killing_Featureslg10 - Min.Per.Feature.lg10) / (Max.Per.Feature.lg10 -Min.Per.Feature.lg10) )%>%
    mutate(Min.Per.Feature = min(Avg_Killing_Features))%>%
  mutate(Max.Per.Feature = max(Avg_Killing_Features))%>%
  mutate(Norm_Avg_Killing_Features = (Avg_Killing_Features - Min.Per.Feature) / (Max.Per.Feature -Min.Per.Feature) )%>%
## Here make two sets of calculations a log10 transformed average and a raw average
  ungroup()%>%
    select(-Min.Per.Feature,
         -Max.Per.Feature,
         -Min.Per.Feature.lg10,
         -Max.Per.Feature.lg10)


# Previously I saw a NaN error when caluclating MDK90%. In cases when this is never reached across all isolates, when we normalize b (x-xmin)/(xmin-xmax) we will get a final calculation zero / zero. This will throw a NaN error in R. To avoid this, I made an if_else statment which will consider For all MDK values across each isolate, if the Norm_Avg_Killing_Features.lg10 is a NaN value, it will be changed to 1 otherwise it remains the value it was defined.
All.Iso.Analysis  <- All.Iso.Analysis %>%
  ungroup()%>%
  group_by(ExpFile,
           Abx.con,
           Isolate,
           Killing.Features)%>%
  mutate(Norm_Avg_Killing_Features.lg10 =  if_else(is.nan(Norm_Avg_Killing_Features.lg10) & Main.KF == "MDK" ,
                                                   1,
                                                   Norm_Avg_Killing_Features.lg10))%>%
  
 mutate(Norm_Avg_Killing_Features =  if_else(is.nan(Norm_Avg_Killing_Features) & Main.KF == "MDK" ,
                                                   72,
                                                   Norm_Avg_Killing_Features))
  
All.Iso.Analysis <- left_join(All.iso,
                              All.Iso.Analysis)




### Reordering axis ----------------------------------------------------------------------------------------------------------------------------------------------------------
time.kill.def.part_I <- ("P.Ila")
time.kill.def.part_II <- ("scvs")


order.x.axis.LC.fraction <- c("LCF_3h",
                              "LCF_6h",
                              "LCF_9h",
                              "LCF_12h",
                              "LCF_24h",
                              "LCF_36h",
                              "LCF_48h",
                              "LCF_60h",
                              "LCF_72h")

order.x.axis.AUC.realtime <- c("AUCrt_0.3h",
                               "AUCrt_0.6h",
                               "AUCrt_0.9h",
                               "AUCrt_0.24h",
                               "AUCrt_0.36h",
                               "AUCrt_0.48h",
                               "AUCrt_0.60h",
                               "AUCrt_0.72h",
                                "AUCrt_48.72h")

order.x.axis.MDK <- c("MDK_25pct",
                      "MDK_50pct",
                      "MDK_75pct",
                      "MDK_90pct")

order.x.axis.overall.order <- c(order.x.axis.LC.fraction,
                               order.x.axis.AUC.realtime,
                               order.x.axis.MDK)

order.x.axis.overall.order <- paste(time.kill.def.part_I  ,
                                    order.x.axis.overall.order,
                                    time.kill.def.part_II,
                                    sep="_")
order.x.axis.overall.order

#Houskeeping
rm(order.x.axis.LC.fraction,
   order.x.axis.AUC.realtime,
   order.x.axis.MDK)

Killing.features.unique <- unique(All.Iso.Analysis$Killing.Features)


# Difine a time-kill curve definiton

All.Iso.Analysis <- All.Iso.Analysis %>%
  mutate(Time.Kill.curve.defintion = Killing.Def)

All.Iso.Analysis <- separate(All.Iso.Analysis, col =Time.Kill.curve.defintion,
                                  into = c("Analysis",
                                           "Feature",
                                           "Cell.type"),
                                  sep="_")

All.Iso.Analysis<- All.Iso.Analysis %>%
  mutate(Killing.Def  = paste(Analysis,
                              Cell.type,
                              sep="_"))%>%
  select(-Analysis,
         -Feature,
         -Cell.type)



# SECTION 12:  ASCT Heatmap of killing features ----
### Ordering as per phylogenetic tree facet_grid ----------------------------------------------------------------------------------------------------------------------------------------------------------


colors.gpt <- colorRampPalette(c("Darkblue", "white", "red"))(10)


exp.name <- unique(Pop.TKC.corr$ExpFile)

exp.name <- paste("/",exp.name,
                  sep="")


LiveCells <- paste(exp.resDir,exp.name,"_Pop_Heatmap_lg10_Avg_Norm_KillFeatures_OrderdPhyloTree.pdf", sep = "")

gg <- All.Iso.Analysis  %>%
  filter(Killing.Def == paste(time.kill.def.part_I,
                              "_",
                              time.kill.def.part_II,
                              sep=""))%>%
  mutate(Isolate.label = sub("Iso.", "", Isolate))
 

gg$Killing.Features <- factor(gg$Killing.Features, levels = order.x.axis.overall.order)

library(gtools)


setwd(genDir)
  Clinical.Iso.List <- read.csv("ClinicalIsolates_LookUpMasterList_BWupdated.csv")
  Clinical.Iso.List <- Clinical.Iso.List  %>%
  mutate( Isolate = paste("Iso.",Sample_ID_Full,
                         sep=""))%>%
  select(Isolate,
         Subspecies,
         OrderingVector)%>%
  mutate(Subspecies = if_else(Subspecies == "", "unknown-subspecies" ,Subspecies))



gg <- left_join(gg,
                Clinical.Iso.List)

gg <- gg %>%
  mutate(Subspecies = ifelse(is.na(Subspecies), "unknown-subspecies", Subspecies))



gg <-gg %>%
 mutate(Isolate = paste(Subspecies,
                        " ",
                        Isolate,
                        sep =""))

gg <- gg %>%
  ggplot(aes(x = Killing.Features, 
             y = reorder(Isolate, desc(OrderingVector)), 
             fill = Norm_Avg_Killing_Features.lg10,
             label = Norm_Avg_Killing_Features.lg10)) +
   geom_tile() +
  ggtitle("lg10 Normalised time-kill features") +
  ylab("Isolates") +
  xlab("Time Kill features") +
      theme_bw()+
   geom_text(aes(label = signif(Norm_Avg_Killing_Features.lg10, digits = 3)), color = "black", size = 0.2 )+
  theme(
    panel.background = element_rect(fill = "black"),
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 7),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 2),
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    strip.text = element_text(size = 5),
    strip.background = element_blank(), #remove background for facet labels
    panel.border = element_rect(colour = "black", fill = NA), #add black border
    panel.spacing = unit(0.05, "lines"))+
  scale_fill_gradientn(colors = colors.gpt, breaks = seq(0.0, 1.0, by = 0.1), na.value = "gray") +
  theme(panel.spacing = unit(0.25, "lines"))+
   facet_grid_paginate(Subspecies ~., 
                       scales = "free", 
                       space = "free" ,
                      ncol =1,
                      nrow =5,
                      page =1)
n <- n_pages(gg)

pdf(LiveCells ,paper = "a4", width = 20 , height = 15 )
for(i in 1:n){
    print(gg + facet_grid_paginate(Subspecies ~., 
                                   scales = "free", 
                               space = "free" ,
                      ncol =1,
                      nrow =5, page = i)) 
}
dev.off()


#Houskeeping
rm(gg,
   gg.mtrx,
   hc,
   df_ordered,
   dist_mat,
   Well_coordinate,
   ZB.df,
   QC,
   qc.Condition,
   QC.n,
   QC.n.thrs,
   QC.subset,
   Contamination.df,
   AS.All.iso.df,
   Average.Result.perIso.PASSED.all.Criteria,
      Contamination.df,
     FC.complete,
     Flagged.Analysis,
     Flagged.Analysis.overall,
     KF.df,
     KF.per.drug,
     Main.KF.df,
     Outlier.wells)


# New ordering strategy


# SECTION 13: Exporting killing features table ----


## Average killing features result ----
expID <- unique( Pop.TKC.corr$ExpFile)

    library(stringr)

# Export log10 transformed and raw values

    All.Iso.Analysis.WIDER.subset <- All.Iso.Analysis %>%
      ungroup()%>%
      mutate(Avg_Killing_Features = as.character(Avg_Killing_Features),
             Avg_Killing_Features = if_else( Avg_Killing_Features == "72" & Main.KF == "MDK", ">72", Avg_Killing_Features))
    
    # -----Order of colunms
  time.kill.def.part_I <- ("P.Ila")
  time.kill.def.part_II <- ("scvs")

  

  order.x.axis.LC.fraction <- c("LCF_3h",
                              "LCF_6h",
                              "LCF_9h",
                              "LCF_12h",
                              "LCF_24h",
                              "LCF_36h",
                              "LCF_48h",
                              "LCF_60h",
                              "LCF_72h")

  order.x.axis.AUC.realtime <- c("AUCrt_0.3h",
                               "AUCrt_0.6h",
                               "AUCrt_0.9h",
                               "AUCrt_0.24h",
                               "AUCrt_0.36h",
                               "AUCrt_0.48h",
                               "AUCrt_0.60h",
                               "AUCrt_0.72h",
                                "AUCrt_48.72h")
  
  

  order.x.axis.AUC.realtimeLOG <- c("AUCrtLOG_0.3h",
                               "AUCrtLOG_0.6h",
                               "AUCrtLOG_0.9h",
                               "AUCrtLOG_0.24h",
                               "AUCrtLOG_0.36h",
                               "AUCrtLOG_0.48h",
                               "AUCrtLOG_0.60h",
                               "AUCrtLOG_0.72h",
                                "AUCrtLOG_48.72h")
  


  order.x.axis.MDK <- c("MDK_25pct",
                      "MDK_50pct",
                      "MDK_75pct",
                      "MDK_90pct")

  order.x.axis.overall.order <- c(order.x.axis.LC.fraction,
                               order.x.axis.AUC.realtime,
                                order.x.axis.AUC.realtimeLOG,
                               order.x.axis.MDK)

  # Ordering all time kill cuvre defitions
  order.x.axis.overall.order.P.Ila_scvs <- paste(time.kill.def.part_I  ,
                                    order.x.axis.overall.order,
                                    time.kill.def.part_II,
                                    sep="_")
  
  
  order.x.axis.overall.order.P.Ila_sc <- gsub(  "vs",
                                                "",
                                                order.x.axis.overall.order.P.Ila_scvs)
  
  order.x.axis.overall.order.P.Thr_scvs <- gsub(  "Ila",
                                                "Thr",
                                                order.x.axis.overall.order.P.Ila_scvs)
  
  order.x.axis.overall.order.P.Thr_sc <- gsub(  "vs",
                                                "",
                                                   order.x.axis.overall.order.P.Thr_scvs )
  
    if(Two.BaSic.def == "Yes") {
      
       order.x.axis.overall.order.P.Ila_scBaSic <- paste(    order.x.axis.overall.order.P.Ila_sc,
                                                            "BaSic",sep="")
       
         order.x.axis.overall.order.P.Ila_scvsBaSic <- paste(    order.x.axis.overall.order.P.Ila_scvs,
                                                            "BaSic",sep="")
         
         order.x.axis.overall.order.P.Thr_scBaSic <- paste(      order.x.axis.overall.order.P.Thr_sc,
                                                            "BaSic",sep="")
         
           order.x.axis.overall.order.P.Thr_scvsBaSic <- paste(      order.x.axis.overall.order.P.Thr_scvs,
                                                            "BaSic",sep="")
  
      
      
    }


    if(Two.BaSic.def == "Yes") {
      
      # If multiple Basic outputs are present then we order Basic output aswell
      general.order <- c( order.x.axis.overall.order.P.Ila_sc,
                              order.x.axis.overall.order.P.Ila_scBaSic,
                              order.x.axis.overall.order.P.Ila_scvs,
                              order.x.axis.overall.order.P.Ila_scvsBaSic,
                              order.x.axis.overall.order.P.Thr_sc,
                              order.x.axis.overall.order.P.Thr_scBaSic,
                              order.x.axis.overall.order.P.Thr_scvs,
                              order.x.axis.overall.order.P.Thr_scvsBaSic)
      
      
    } else (   
      # Otherwise, we 
      general.order <- c( order.x.axis.overall.order.P.Ila_sc,
                               order.x.axis.overall.order.P.Ila_scvs,
                               order.x.axis.overall.order.P.Thr_sc,
                               order.x.axis.overall.order.P.Thr_scvs) 
      )
  
      general.order <- paste(gsub("_",".",experiment.drug),
                             general.order,
                             experiment.id,
                             sep="_")
   
      general.order <- c(#"ExpFile",
                         #"Abx.con",
                         "Isolate",
                       #  "Killing.Def",
                         general.order)
      
    #  general.order
      
        All.Iso.Analysis.WIDER.subset <- All.Iso.Analysis.WIDER.subset %>%
        filter(grepl("P.Ila|P.Thr", Killing.Def))%>%
        select(ExpFile,
                 Abx.con,
                 Isolate,
                 Killing.Features,
                 Avg_Killing_Features,
                Avg_Killing_Featureslg10)
         
      All.Iso.Analysis.WIDER.subset <-   All.Iso.Analysis.WIDER.subset %>%
      mutate(Abx.con = gsub("_",".",Abx.con))
             
             
             All.Iso.Analysis.WIDER.subset <-   All.Iso.Analysis.WIDER.subset %>%
               mutate(Killing.Features = paste(Abx.con,
                                               Killing.Features,
                                               ExpFile,
                                               sep="_"))
               
            # grabbing the lg10 data need to add reference to heading   
           All.Iso.Analysis.WIDER.subset.lg10 <-  All.Iso.Analysis.WIDER.subset %>%
            select(Isolate,
                 Killing.Features,
                  Avg_Killing_Featureslg10)
               
            All.Iso.Analysis.WIDER.subset <- All.Iso.Analysis.WIDER.subset %>%
               select(Isolate,
                      Killing.Features,
                      Avg_Killing_Features)
    
  All.Iso.Analysis.WIDER.subset <- pivot_wider( All.Iso.Analysis.WIDER.subset,
                                                  names_from = Killing.Features,
                                                values_from = Avg_Killing_Features)
  
  All.Iso.Analysis.WIDER.subset.Organised.KF <-   All.Iso.Analysis.WIDER.subset[,general.order ]  
  
  rm(All.Iso.Analysis.WIDER.subset)
  
# ---- Generating output of log10 data
  
           All.Iso.Analysis.WIDER.subset.lg10 <- separate(   All.Iso.Analysis.WIDER.subset.lg10 , col = Killing.Features,
                                                             into = c("c1",
                                                                      "c2",
                                                                      "c3",
                                                                      "c4",
                                                                      "c5",
                                                                      "c6"),
                                                             sep = "_")
           
           
               All.Iso.Analysis.WIDER.subset.lg10 <- All.Iso.Analysis.WIDER.subset.lg10 %>%   
                 mutate(Killing.Features = paste(c1,
                                                 c2,
                                                 c3,
                                                 paste(c4,
                                                       ".lg10",
                                                       sep=""),
                                                 c5,
                                                 c6,
                                                 sep="_"))%>%
                 select(Isolate,
                        Killing.Features,
                        Avg_Killing_Featureslg10)
  
    # ---- Defining Order of colunms
               
  time.kill.def.part_I <- ("P.Ila")
  time.kill.def.part_II <- ("scvs")


  order.x.axis.LC.fraction <- c("LCF_3h.lg10",
                              "LCF_6h.lg10",
                              "LCF_9h.lg10",
                              "LCF_12h.lg10",
                              "LCF_24h.lg10",
                              "LCF_36h.lg10",
                              "LCF_48h.lg10",
                              "LCF_60h.lg10",
                              "LCF_72h.lg10")

  order.x.axis.AUC.realtime <- c("AUCrt_0.3h.lg10",
                               "AUCrt_0.6h.lg10",
                               "AUCrt_0.9h.lg10",
                               "AUCrt_0.24h.lg10",
                               "AUCrt_0.36h.lg10",
                               "AUCrt_0.48h.lg10",
                               "AUCrt_0.60h.lg10",
                               "AUCrt_0.72h.lg10",
                                "AUCrt_48.72h.lg10")
  
  
   order.x.axis.AUC.realtimeLOG <- c("AUCrtLOG_0.3h.lg10",
                               "AUCrtLOG_0.6h.lg10",
                               "AUCrtLOG_0.9h.lg10",
                               "AUCrtLOG_0.24h.lg10",
                               "AUCrtLOG_0.36h.lg10",
                               "AUCrtLOG_0.48h.lg10",
                               "AUCrtLOG_0.60h.lg10",
                               "AUCrtLOG_0.72h.lg10",
                                "AUCrtLOG_48.72h.lg10")

  order.x.axis.MDK <- c("MDK_25pct.lg10",
                      "MDK_50pct.lg10",
                      "MDK_75pct.lg10",
                      "MDK_90pct.lg10")

  order.x.axis.overall.order <- c(order.x.axis.LC.fraction,
                               order.x.axis.AUC.realtime,
                                  order.x.axis.AUC.realtimeLOG,
                               order.x.axis.MDK)

  # Ordering all time kill cuvre defitions
  order.x.axis.overall.order.P.Ila_scvs <- paste(time.kill.def.part_I  ,
                                    order.x.axis.overall.order,
                                    time.kill.def.part_II,
                                    sep="_")
  
  
  order.x.axis.overall.order.P.Ila_sc <- gsub(  "vs",
                                                "",
                                                order.x.axis.overall.order.P.Ila_scvs)
  
  order.x.axis.overall.order.P.Thr_scvs <- gsub(  "Ila",
                                                "Thr",
                                                order.x.axis.overall.order.P.Ila_scvs)
  
  order.x.axis.overall.order.P.Thr_sc <- gsub(  "vs",
                                                "",
                                                   order.x.axis.overall.order.P.Thr_scvs )
  
  
    if(Two.BaSic.def == "Yes") {
      
       order.x.axis.overall.order.P.Ila_scBaSic <- paste(    order.x.axis.overall.order.P.Ila_sc,
                                                            "BaSic",sep="")
       
         order.x.axis.overall.order.P.Ila_scvsBaSic <- paste(    order.x.axis.overall.order.P.Ila_scvs,
                                                            "BaSic",sep="")
         
         order.x.axis.overall.order.P.Thr_scBaSic <- paste(      order.x.axis.overall.order.P.Thr_sc,
                                                            "BaSic",sep="")
         
           order.x.axis.overall.order.P.Thr_scvsBaSic <- paste(      order.x.axis.overall.order.P.Thr_scvs,
                                                            "BaSic",sep="")
  
      
      
    }


    if(Two.BaSic.def == "Yes") {
      
      # If multiple Basic outputs are present then we order Basic output aswell
      general.order <- c( order.x.axis.overall.order.P.Ila_sc,
                              order.x.axis.overall.order.P.Ila_scBaSic,
                              order.x.axis.overall.order.P.Ila_scvs,
                              order.x.axis.overall.order.P.Ila_scvsBaSic,
                              order.x.axis.overall.order.P.Thr_sc,
                              order.x.axis.overall.order.P.Thr_scBaSic,
                              order.x.axis.overall.order.P.Thr_scvs,
                              order.x.axis.overall.order.P.Thr_scvsBaSic)
      
      
    } else (   
      # Otherwise, we 
      general.order <- c( order.x.axis.overall.order.P.Ila_sc,
                               order.x.axis.overall.order.P.Ila_scvs,
                               order.x.axis.overall.order.P.Thr_sc,
                               order.x.axis.overall.order.P.Thr_scvs) 
      )
  
      
     
      general.order <- paste(gsub("_",".",experiment.drug),
                             general.order,
                             experiment.id,
                             sep="_")
   
      general.order <- c(
                         "Isolate",
                      
                         general.order)
      
#  general.order
  
   All.Iso.Analysis.WIDER.subset.lg10 <- pivot_wider( All.Iso.Analysis.WIDER.subset.lg10,
                                                  names_from = Killing.Features,
                                                values_from = Avg_Killing_Featureslg10)
  
  All.Iso.Analysis.WIDER.subset.Organised.KF.lg10 <-   All.Iso.Analysis.WIDER.subset.lg10[,general.order ]  
  


### Exporting  lg10 Average killing features result  ----------------------------------------------------------------------------------------------------------------------------------------------------------
All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- left_join(All.Iso.Analysis.WIDER.subset.Organised.KF,
                                                                  All.Iso.Analysis.WIDER.subset.Organised.KF.lg10,
                                                                  by = "Isolate")

    
    All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- data.frame(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED )
    
    
      All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED  <- All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED%>%
      mutate_all(as.character)
    
      All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- melt(      All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                                                  id.vars = "Isolate")
      
            All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <-       All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
              group_by(Isolate)%>%
              mutate(Col.order = 1,
                     Col.order = cumsum(Col.order))%>%
              ungroup()
            
          
# -- Editing column heading to reflect Drug_Analysis_AnalysisDef_MainKF_Specific.variable.KF_ExperimentID and defining the order  of colunms.
      
            All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- separate(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,variable, into = c("Drug", 
                                                                                                                                               "Analysis",
                                                                                                                                              
                                                                                                                                               "Main.KF",
                                                                                                                                               "Specific.var.KF",
                                                                                                                                                "Analysis.Def",
                                                                                                                                               "Exp.ID"),
                                                                               sep = "_")
            
       All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED<-        All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
                   mutate(Specific.var.KF = gsub(".lg10", "-lg10", Specific.var.KF))%>%
                   mutate(Full.Killing.Feature = paste(Drug,
                                                       Analysis,
                                                       Analysis.Def,
                                                       Main.KF,
                                                       Specific.var.KF,
                                                       Exp.ID,sep="_"))%>%
                   select(Isolate,
                          Full.Killing.Feature,
                          value,
                          Col.order)%>%
         ungroup()

# Defining colunm order
       Colunm.order.def <-  All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
         ungroup()%>%
         select(Full.Killing.Feature,
                Col.order)%>%
         distinct()
       
             Colunm.order.def <-Colunm.order.def$Full.Killing.Feature
               Colunm.order.def <- c("Isolate",
                                     Colunm.order.def)
       

                 All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED<-All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
                   ungroup()%>%
                   select(-Col.order)

  All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED<- spread(  All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                                                                                     Full.Killing.Feature,
                                                value)

        All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <-  All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
          select(Colunm.order.def)


## SECTION 14: Removing isolates with low LCF  & Low numbers to determine LCF ----------------------------------------------------------------------------------------------------------------------------------------------------------
# We will finally only consider the isolates which across all definitions had an inital live cell fraction fo at least 80%. Otherwise it is omitted. This is considered from the control experiment or considering all the initial LCF in the ASCT experiment per isolate.

minLCF.df <- read.csv(path.to.minimum.LCF.data)

min.LCF.df.subset <-  All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
  filter(Isolate %in% minLCF.df$Isolate )%>%
   mutate(across(-Isolate, ~ as.character(NA)))

All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED  <- All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
  filter(!Isolate %in% minLCF.df$Isolate )

All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- rbind(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                                                                min.LCF.df.subset )


rm(
   min.LCF.df.subset)



#--- Too few numbers  for max LCF ----


too.few.LCF.df <- read.csv(path.to.isol.too.few.numbers.LCF.data)

too.few.LCF.df.subset <-  All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
  filter(Isolate %in% too.few.LCF.df$Isolate )%>%
   mutate(across(-Isolate, ~ as.character(NA)))

All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED  <- All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
  filter(!Isolate %in% too.few.LCF.df$Isolate )

All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- rbind(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                                                               too.few.LCF.df.subset )


rm(
   too.few.LCF.df.subset )





## Exporting iTol compatible files ----------------------------------------------------------------------------------------------------------------------------------------------------------
filename2 <- paste(exp.resDir,
                    "/",
                    expID,
                    "_",
                    experiment.drug,
                    "_",
                    "Pop.Clinial.Isolates.KillingFeatures.csv",
                    sep="")
                    
  # Saving ASCT Clinical isoalte data
    write.csv(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                file =filename2,
              row.names=FALSE)


    All.Iso.Analysis.WIDER <- All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED
    
    #Houskeeping
     rm(    All.Iso.Analysis.WIDER.subset.Organised.KF,
           All.Iso.Analysis.WIDER.subset.Organised.KF.lg10,
           All.Iso.Analysis.WIDER.subset,
           All.Iso.Analysis.WIDER.subset.lg10,
              KF.labs.df)
  

All.Iso.Analysis.LONG<- pivot_longer(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED, 
                        cols = -Isolate)


All.Iso.Analysis.LONG <- All.Iso.Analysis.LONG %>%
  mutate(Organised_Killing_Features = name)%>%
  select(Isolate,
         Organised_Killing_Features,
         value)



# Creating iTol file directory
path.iTol.dir <- paste(exp.resDir,
                       "/",
                                 analysis.normalisation,

                       "_iTol-Results",
                       sep="")
dir.create(path.iTol.dir)
#-- 

setwd(genDir)

label.Clinical.Iso <- paste(label_dir,
                            "/ClinicalIsolates_LookUpMasterList_BWupdated.csv",
                            sep="")


iTol.Clinical.Iso.List <- read.csv(label.Clinical.Iso)
iTol.Clinical.Iso.List <-iTol.Clinical.Iso.List %>%
  mutate( Isolate = paste("Iso.",Sample_ID_Full,
                         sep=""))%>%
  select(Isolate,
         Subspecies,
         Lane)%>%
  mutate(Subspecies = if_else(Subspecies == "", "unknown-subspecies" ,Subspecies))


iTol <- left_join(All.Iso.Analysis.LONG,
                  iTol.Clinical.Iso.List)


itol.list.features.loop <- unique(iTol$Organised_Killing_Features)


for ( i in itol.list.features.loop) {
  
  
  iTol.subset <- iTol %>%
    filter(Organised_Killing_Features ==i)
  
  #---- NORMALISE values here HERE 
  
  All.Iso.Analysis.no.NAs <- iTol.subset  %>%
  ungroup()%>%
  drop_na()%>%
  ungroup()%>%
  mutate(Organised_Killing_Features = as.factor(Organised_Killing_Features))%>%
  mutate(value = if_else ( value == ">72", "72", value))%>%
  mutate(value =as.numeric(value))%>%
  group_by(
           Organised_Killing_Features)%>%
  mutate(Min.value = min(value))%>%
  mutate(Max.value= max(value))%>%
  mutate(Norm_value = (value - Min.value) / (Max.value -Min.value) )%>%
    ungroup()%>%
    select(-Min.value,
            -Max.value)%>%
    mutate(value = as.character(value))

  
    iTol.subset <- left_join(  iTol.subset,
                                All.Iso.Analysis.no.NAs)
  

  #Housekeeping
    rm(All.Iso.Analysis.no.NAs)
  
  
  iTol.subset <-  iTol.subset %>%
    drop_na()%>%
       filter(!is.na(Lane))%>%
    select(Lane,
           Norm_value)%>%

    mutate(Lane_Killing_feature = paste(Lane,
                                        " ",
                                       Norm_value,
                                        sep=""))%>%
    select(Lane_Killing_feature)
  
  txt <- c("DATASET_GRADIENT",
         "SEPARATOR SPACE",
     
        paste("DATASET_LABEL ", i,sep=""),
         "COLOR #ff0000",
         "STRIP_WIDTH 200",
         "USE_MID_COLOR 1",
         "COLOR_MIN #0000ff",
         "COLOR_MAX #ff0000",
         "COLOR_MID #aaaaaa",
         "DATA",
         paste("Lane " ,i, ")",sep="") )


  txt.data <- c(  iTol.subset$Lane_Killing_feature)
  
  
  iTol.txt <- c(txt,
         txt.data )
  
  
  itol.filename <- paste(
                       "iTol",
                       "_",  
                       analysis.normalisation,
                       "_",
                       i,
                       ".txt",
                       sep="")
  path.iTol.dir.results <- paste(path.iTol.dir ,
                       "/",
                       itol.filename,
                       sep="")
  
  fileConn <- file(path.iTol.dir.results)
  writeLines(iTol.txt , fileConn)
  close(fileConn)

  #Houskeeping
  rm(path.iTol.dir.results,
   itol.filename,
   iTol.txt)
  

  
}
  
#Houskeeping
rm(txt,
   iTol,
   iTol.Clinical.Iso.List,
   path.iTol.dir,Clinical.Iso.List,
   iTol.subset,LOC.df,
   neo.labels)


message("Population analysis has now ended")

# Capture the end time
end_time.pop <- Sys.time()

# Calculate the time difference in minutes
execution_time.pop <- as.numeric(difftime(end_time.pop , start_time, units = "mins"))

# Print the execution time
cat("Script execution time:", execution_time.pop , "minutes\n")

start_time <- Sys.time()

#.... TRACKING ANALYSIS ----
# SECTION 14: Loading tracking data ----
## Live cell fraction data ----

trk.file.pattern <- paste("*_",old.or.new.tracking,"_LC.csv",sep="")

perWell.filenames <- list.files(path = res.tracking.data.Dir,
                            pattern = trk.file.pattern,
                            full.names = FALSE)




#setting path to where all the csv files due to be processed are
setwd(res.tracking.data.Dir)


OC.file.list <- lapply(perWell.filenames,
                      read.csv)



#Naming the elements of my list with the vector list of my file names
Exp.Names <- gsub(".csv",
                  "",
                  perWell.filenames)

Exp.Names <- gsub(paste("_",old.or.new.tracking,
                        "_LC",sep=""),
                  "",
                  Exp.Names)


names(OC.file.list)<- Exp.Names

#Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
perWell.Trk.df <- rapply(
  OC.file.list,
  as.character,
  how = "replace"
)


perWell.Trk.df <- bind_rows(perWell.Trk.df , .id = "Filename")


#Houskeeping
rm(OC.file.list,
   trk.file.pattern)


perWell.Trk.df <- perWell.Trk.df %>%
  mutate(ExpFile = experiment.id)%>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         matches("Trkv")) # Grabbing all the colunms that have Trkv in string name


perWell.Trk.df <- left_join(perWell.Trk.df ,
                            Plot.label.df, 
                            by = c("ExpFile",
                                   "Well_coordinate"))


# Adding time metadata to tracking data
MD.data <- perWell.df %>%
  select(ExpFile,
         #Abx.con,
         Well_coordinate,
         timestep,
         Time_Hrs)

perWell.Trk.df <- left_join(perWell.Trk.df ,
                           MD.data, 
                            by = c("ExpFile",
                                   "Well_coordinate",
                                   "timestep"))

# Need to edit the column names of the different time kill definitions
perWell.Trk.df.melted <- melt(perWell.Trk.df, id.vars = c("ExpFile",
                                                          "Well_coordinate",
                                                          "Abx.con",
                                                          "Isolate",
                                                          "Flag.label",
                                                          "timestep",
                                                          "Time_Hrs"),variable.name = "Time.Kill.Definitions")

perWell.Trk.df.melted <-perWell.Trk.df.melted %>%
  mutate( Time.Kill.Definitions =  Time.Kill.Definitions)%>%
  # Improving time kill defintions names
  mutate(
        Time.Kill.Definitions = gsub("_Ila.Def.1_LCfraction.BaSic", "_Ila1BaSic",Time.Kill.Definitions),
        Time.Kill.Definitions = gsub("_Ila.Def.2_LCfraction.BaSic", "_Ila2BaSic",Time.Kill.Definitions),
        
        Time.Kill.Definitions = gsub("_Ila.Def.1_LCfraction", "_Ila1",Time.Kill.Definitions),
        Time.Kill.Definitions = gsub("_Ila.Def.2_LCfraction", "_Ila2",Time.Kill.Definitions),
        
        Time.Kill.Definitions = gsub("_Thrs.1pos_LCfraction.BaSic", "_Thr1BaSic",Time.Kill.Definitions),
        Time.Kill.Definitions = gsub("_Thrs.2pos_LCfraction.BaSic", "_Thr2BaSic",Time.Kill.Definitions),
        
        Time.Kill.Definitions = gsub("_Thrs.1pos_LCfraction", "_Thr1",Time.Kill.Definitions),
        Time.Kill.Definitions = gsub("_Thrs.2pos_LCfraction", "_Thr2",Time.Kill.Definitions),
        
        Time.Kill.Definitions = gsub("_LCfraction_PiDyn", "_piDyn",Time.Kill.Definitions)
        )

perWell.Trk.df <- spread(perWell.Trk.df.melted ,
                              key = Time.Kill.Definitions , 
                              value = value) 



if ( old.or.new.tracking == "Trkv2") {
  perWell.Trk.df <-perWell.Trk.df %>%
  select(-matches("_piDyn"))
  
}


 perWell.Trk.df <-perWell.Trk.df %>%
  select(-matches("_piDyn.BaSic"))

rm(perWell.Trk.df.melted)
                                
                              


setwd(genDir)
### Tracking numbers data ----------------------------------------------------------------------------------------------------------------------------------------------------------


trk.file.pattern <- paste("*_",old.or.new.tracking,"_LC.N.csv",sep="")


#Creacting vector with the list of file names
perWell.filenames <- list.files(path = res.tracking.data.N.Dir,
                            pattern = trk.file.pattern,
                            full.names = FALSE)

#setting path to where all the csv files due to be processed are
setwd(res.tracking.data.N.Dir)


OC.file.list <- lapply(perWell.filenames,
                      read.csv)

#Naming the elements of my list with the vector list of my file names
Exp.Names <- gsub(".csv",
                  "",
                  perWell.filenames)

Exp.Names <- gsub(paste("_",old.or.new.tracking,"_LC.N",sep=""),
                  "",
                  Exp.Names)


names(OC.file.list)<- Exp.Names

#Irrespective of the number of rows bind all the elemtents of my list using a unique identifier "well_coord" which is based on the name of each element of the list 
perWell.Trk.n.df <- rapply(
  OC.file.list,
  as.character,
  how = "replace"
)


perWell.Trk.n.df <- bind_rows(perWell.Trk.n.df , .id = "Filename")


#Houskeeping
rm(OC.file.list)

perWell.Trk.n.df <- perWell.Trk.n.df %>%
  mutate(ExpFile = experiment.id)%>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         live.dead.def,
         matches("_Ila.Def.2_LC.N"))%>%
  select(-matches(".N.BaSic"))



perWell.Trk.n.df <- melt(perWell.Trk.n.df, id.vars = c("ExpFile",
                                                          "Well_coordinate",
                                                          "timestep",
                                                       "live.dead.def"),
                         variable.name = "Time.Kill.Definitions")

perWell.Trk.n.df <-perWell.Trk.n.df %>%
  mutate(Tracking.N = "Tracking.N")%>%
  select(-Time.Kill.Definitions)


perWell.Trk.n.df <- spread(perWell.Trk.n.df, key = c("Tracking.N"),
                           value = value)
  

perWell.Trk.n.df <- perWell.Trk.n.df %>%
  group_by(ExpFile,
           Well_coordinate,
          timestep)%>%
  mutate(Tracking.N = as.numeric(Tracking.N))%>%
  mutate(T_Tot.N = sum(Tracking.N))%>%
  mutate(ExpFile = experiment.id)%>%
  ungroup()%>%
  select(ExpFile,
         Well_coordinate,
         timestep,
         T_Tot.N)%>%
  distinct()%>%
  ungroup()%>%
  group_by(ExpFile,
           Well_coordinate)%>%
  mutate(T_Tot.N = max(T_Tot.N))%>% ### KEY STEP 
  ungroup()%>%
  distinct()%>%
  mutate(T_Eval_Trck.Numbers = if_else(T_Tot.N > 1000, ">1K", "<1K"))


setwd(genDir)
#### Merging tracking fractions and numbers ----------------------------------------------------------------------------------------------------------------------------------------------------------
perWell.Trk.df <- left_join(perWell.Trk.df,
                            perWell.Trk.n.df)


rm(  perWell.Trk.n.df)



## SECTION 15:  Identfiying QC QC flages from Pop analysis----------------------------------------------------------------------------------------------------------------------------------------------------------

perWell.Trk.df  <- perWell.Trk.df %>%
  mutate(Contamination = if_else(grepl("Contam", Flag.label), "Contamination", "No-Contamination"))%>%
  mutate(Growth = if_else(grepl("Growth", Flag.label), "Growth", "No-Growth"),
         Artefact = if_else(grepl("AF",  Flag.label),"Artefact", "No-Artefact"))



### Ignore contaminaiton flags when using tracking----------------------------------------------------------------------------------------------------------------------------------------------------------
# Given contmainaition occures late in experiment and tracking corrects for that we will ingnore the conaminaiton flages generated in population analysis
perWell.Trk.df <- perWell.Trk.df %>%
  mutate(Contamination = "No-Contamination")



### Ignore Artefacts flags when using tracking----------------------------------------------------------------------------------------------------------------------------------------------------------
perWell.Trk.df <- perWell.Trk.df %>%
  mutate(Artefact = "No-Artefact")



## SECTION 16:  Normalisation based on Top2norm  ----------------------------------------------------------------------------------------------------------------------------------------------------------

T.TKC <- perWell.Trk.df %>%
  select(-T_Eval_Trck.Numbers,
         -Flag.label,
                       -Contamination,
                       -Growth,
                       -Artefact,
         -T_Eval_Trck.Numbers) # NEW


T.TKC <- melt(T.TKC,
           id.vars = c("ExpFile",
                       "Abx.con",
                                  "Well_coordinate",
                       "Isolate",
                       "timestep",
                       "Time_Hrs"),
                      variable.name =  "Time.Kill.Definitions")


T.TKC <- T.TKC %>%
  mutate(timestep = as.numeric(timestep),
         Time_Hrs = as.numeric(Time_Hrs),
         value = as.numeric(value))

# Loading external input data for maximum live cell fraction
maxLCF.df <- read.csv(maxLCF.df.datapath )

maxLCF.df <-maxLCF.df %>%
  mutate(Isolate = if_else(Isolate =="ATc.19979", "Iso.ATc.19979", Isolate))


# Finding killing defitnions
defintions.considered.in.analysis <- unique((T.TKC$Time.Kill.Definitions))
defintions.considered.in.analysis <- as.character(defintions.considered.in.analysis)

defintions.considered.in.analysis <- c("Isolate",defintions.considered.in.analysis)
defintions.considered.in.analysis <- subset(defintions.considered.in.analysis, defintions.considered.in.analysis != "T_Tot.N")

# Select isolate 


# Change column name Trk1
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_LCfraction_PiDyn.BaSic"] <- "Trkv1_piDynBaSic"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_LCfraction_PiDyn"] <- "Trkv1_piDyn"

# BaSic colunm names
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_Ila.Def.1_LCfraction.BaSic"] <- "Trkv1_Ila1BaSic"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_Ila.Def.2_LCfraction.BaSic"] <- "Trkv1_Ila2BaSic"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_Thrs.1pos_LCfraction.BaSic"] <- "Trkv1_Thr1BaSic"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_Thrs.2pos_LCfraction.BaSic"] <- "Trkv1_Thr2BaSic"

colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_Ila.Def.1_LCfraction"] <- "Trkv1_Ila1"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_Ila.Def.2_LCfraction"] <- "Trkv1_Ila2"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_Thrs.1pos_LCfraction"] <- "Trkv1_Thr1"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv1_Thrs.2pos_LCfraction"] <- "Trkv1_Thr2"

# Change column name Trk2
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv2_Ila.Def.1_LCfraction.BaSic"] <- "Trkv2_Ila1BaSic"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv2_Ila.Def.2_LCfraction.BaSic"] <- "Trkv2_Ila2BaSic"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv2_Thrs.1pos_LCfraction.BaSic"] <- "Trkv2_Thr1BaSic"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv2_Thrs.2pos_LCfraction.BaSic"] <- "Trkv2_Thr2BaSic"

colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv2_Ila.Def.1_LCfraction"] <- "Trkv2_Ila1"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv2_Ila.Def.2_LCfraction"] <- "Trkv2_Ila2"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv2_Thrs.1pos_LCfraction"] <- "Trkv2_Thr1"
colnames(maxLCF.df)[colnames(maxLCF.df) == "Trkv2_Thrs.2pos_LCfraction"] <- "Trkv2_Thr2"

if ( old.or.new.tracking == "Trkv1") {
  

  maxLCF.df <- maxLCF.df %>%
  select(-matches("Trkv1_PiDynBaSic"))%>%
  select(all_of(defintions.considered.in.analysis))

  
}


  
  rm(defintions.considered.in.analysis)
# left join such that I have a Max.LCF.of.Iso.at.tp0
  maxLCF.df <- melt(  maxLCF.df,
                      id.vars =  c("Isolate"),
                      variable.name = "Time.Kill.Definitions",
                      value.name = "maxLCF",
                      )

#Finding the maximum live cell fraction at timepoint 0 
T.TKC.tp0 <- T.TKC %>%
 dplyr::ungroup()%>%
   left_join(maxLCF.df)%>%
  dplyr::group_by(ExpFile,
                  Isolate,
                  Time.Kill.Definitions)%>%
  mutate(Max.LCF.of.Iso.at.tp0 = maxLCF)%>%
    select(-maxLCF)

T.TKC.corr <- left_join(T.TKC,
                                T.TKC.tp0  )

#Houskeeping


T.TKC.corr <- T.TKC.corr %>%
    dplyr::group_by(ExpFile,
                   Abx.con,
                   Well_coordinate,
                  Isolate,
                  Time.Kill.Definitions)%>%
  mutate(LC.fraction.corr = (1/first(Max.LCF.of.Iso.at.tp0))*value)%>%
  mutate(timestep = timestep + 1)%>%
  select(-value)


T.TKC.corr.0th.tp <- T.TKC.corr %>%
  ungroup()%>%
   select(ExpFile,
                   Abx.con,
                  
                   Well_coordinate,
                  Isolate,
                  Time.Kill.Definitions)%>%
    distinct()%>%
    mutate(timestep = 0,
           Time_Hrs = 0,
          LC.fraction.corr = 1)


T.TKC.corr <- rbind(T.TKC.corr,
                  T.TKC.corr.0th.tp)



T.TKC.corr <- T.TKC.corr %>%
  select(-Max.LCF.of.Iso.at.tp0)


#Housekeeping
rm(T.TKC.corr.0th.tp,
   T.TKC,
   T.TKC.tp0)



write.csv(T.TKC.corr, paste(exp.resDir,
                                 "/",
                                 unique(T.TKC.corr$ExpFile),
                            "_",
                            old.or.new.tracking,
                                 "_MultiBaSic.csv",
                                 sep=""),
          row.names = FALSE)

QC.T <- perWell.Trk.df %>%
  ungroup()%>%
  select(ExpFile,
         Well_coordinate,
         Isolate,
          Contamination,
         Growth,
         Artefact,
         T_Tot.N,
         T_Eval_Trck.Numbers)%>%
  distinct()



T.TKC.corr <- left_join(T.TKC.corr,
                        QC.T)



## SECTION 17: Detecting Outliers for tracking using Ilastik once positive classification ----
##  Assigning replicate numbers to isolates ------------

  KF.df <- Main.KF.df.TRACKING %>%
  filter(Analysis.Str == old.or.new.tracking & cell.type == "Ila2" & Main.KF == "LCF")
  

  Iso.Rep.N <- Conditions %>%
  ungroup()%>%
  mutate(Isolate = if_else(Isolate == "ATc.19979", "Iso.ATc.19979", Isolate))%>%
  group_by(ExpFile,
           Isolate)%>%
  mutate(rep = 1,
         rep = cumsum(rep))%>%
    mutate(Abx.con = paste(Abx,
                           Concentration,
                           sep="_"))%>%
    select(ExpFile,
           Abx.con,
           Isolate,
           rep)%>%
    distinct()


  # Killing features of OutLiers Coefficient of Variation
  KF.OL.CV.df  <- KF.df %>%
    mutate(value = if_else(value == ">72", "72", value))%>%
  mutate(value = as.numeric(value),
         ExpID = ExpFile)%>%
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
    mutate(Isolate = paste("Iso.",Isolate,
                           sep=""))%>%
  ungroup()%>%
    
    select(ExpFile,
           ExpID,
           Abx.con,
           Isolate,
           Well_coordinate,
           Colunm,
           Killing.Features,
           Killing.Def,
           Main.KF,
           value)
  
 
    KF.OL.CV.df <-   KF.OL.CV.df  %>%
    left_join(Iso.Rep.N)%>%
    drop_na()
  

  # Houskeeping
  rm(Iso.Rep.N)




#### Coefficient of variation of triplicates  ----------------------------------------------------------------------------------------------------------------------------------------------------------

  KF.OL.CV.df <- KF.OL.CV.df 


# KF.OL.3cv.Mean.perRep sHould be changed to mean.perISo
KF.OL.3cv.Mean.perRep <- KF.OL.CV.df %>%
  ungroup()%>%
 # filter(Isolate =="Iso.857")%>%
  mutate(rep = as.numeric(rep))%>%
  mutate(value = if_else(value > 1 , 1, value))%>%
  mutate(value.log10 = signif(log10(value), digits = 3))%>%
  group_by(Abx.con,
           Isolate,
         #  rep,
           Killing.Features)%>%
  mutate(CV.triplicate = sd(value)/mean(value))%>%
  mutate(CV.triplicate.lg10 = sd(value.log10)/mean(value.log10))%>%
  ungroup()%>%
   group_by(Abx.con,
           Isolate)%>%
  mutate(CV.triplicate.Mean.perRep = mean(CV.triplicate),
         CV.triplicate.Mean.perRep.lg10 = mean(CV.triplicate.lg10))



KF.OL.3cv.Mean.perRep <- KF.OL.3cv.Mean.perRep %>%
  ungroup()%>%
  select(-Killing.Features,
         -value,
         -Main.KF,
         -Colunm,
         -value.log10,
         -CV.triplicate,
         -CV.triplicate.lg10)%>%
  distinct()
  


#### Coefficient of variation of duplicates  ----------------------------------------------------------------------------------------------------------------------------------------------------------


CV.combo.res <- data.frame()

  KF.OL.CV.df <- KF.OL.CV.df %>%
  mutate(loop.var = paste(KF.OL.CV.df$Abx.con,KF.OL.CV.df$Isolate,sep="_"))%>%
  mutate(value.lg10 = signif(log10(value), digits = 3))


loop.list <- unique(paste(KF.OL.CV.df$Abx.con,KF.OL.CV.df$Isolate,sep="_"))


#x <- loop.list[1]
#x
for ( x in loop.list) {
  
  
  KF.OL.CV.df.subset <- KF.OL.CV.df %>%
  filter(loop.var == x)
  
  # list.var <- unique(Main.KF.df.TRACKING$Killing.Features)
  a <- unique((KF.OL.CV.df.subset$Well_coordinate))
  a
  b <- unique((KF.OL.CV.df.subset$Well_coordinate))
  b
  combo <- crossing(a,b)
  
  
  #Housleeping
  rm(a,
     b)
  
  
  combo <- combo %>%
    mutate(Identical.combo = if_else (a == b , "Omit", "keep"))%>%
    filter( Identical.combo == "keep")
  
  
    combo <- combo %>%
    mutate(Well_rep_combinations = paste(a,b,sep="_"))
    
    
    list.of.duplicates <- combo$Well_rep_combinations

    
 # for ( i in list.of.combo.a  ) {
    for ( i in     list.of.duplicates  ) {
 #   print(i)
    
      # before and after the _ 
      i.Well_rep1 <- sub("_.*", "", i)  
       i.Well_rep1 
      i.Well_rep2 <- sub(".*_", "", i)  
    # for ( j in list.of.combo.b ) {
    #     print(j)
    #   
      
      Subset.combo <-KF.OL.CV.df.subset %>%
        filter(Well_coordinate == i.Well_rep1 | Well_coordinate ==  i.Well_rep2)%>%
        mutate(Combo = paste(i))%>%
        ungroup()%>%
        group_by(Killing.Features)%>%
        mutate(CV.combo2 = sd(value)/mean(value))%>%
        mutate(CV.combo2.lg10 = sd(value.lg10)/mean(value.lg10))%>%
        select(-value,
               -value.lg10,
               -Colunm,
               )%>%
        distinct()
      
        CV.combo.res <- rbind(  CV.combo.res ,Subset.combo)
        
          rm(i,
           Subset.combo)
    
      
      
    }
    # break
    }
  
  #Houskeeping
  rm(combo,
     KF.OL.CV.df.subset)
  

 



#### Creating dataframe for best CoV ----------------------------------------------------------------------------------------------------------------------------------------------------------

KF.OL.2cv.perRep.Mean.perRep <- CV.combo.res %>%
  ungroup()%>%
  select(ExpID,
     #    Exp.Iso,
         Abx.con,
         Killing.Def,
        Isolate,
         Combo,
         CV.combo2,
         CV.combo2.lg10)%>%
  group_by(ExpID,
      #   Exp.Iso,
         Abx.con,
         Killing.Def,
         Isolate,
         Combo)%>%
  mutate(CV.duplicate.Mean.perRep = mean(CV.combo2),
         CV.duplicate.Mean.perRep.lg10 = mean(CV.combo2.lg10))%>%
  ungroup()%>%
  select(-CV.combo2.lg10,
         -CV.combo2)%>%
  distinct()%>%
   group_by(ExpID,
      #   Exp.Iso,
         Abx.con,
         Killing.Def,
         Isolate)%>%
  mutate(CV.duplicate.Mean.perRep.BestDup = min(CV.duplicate.Mean.perRep),
         CV.duplicate.Mean.perRep.BestDup.lg10 = min( CV.duplicate.Mean.perRep.lg10))%>%
   mutate(Best.combo.Dup = if_else(CV.duplicate.Mean.perRep.BestDup == CV.duplicate.Mean.perRep , "Best-combo","Poor-combo"))%>%
   mutate(Best.combo.Dup.lg10 = if_else(CV.duplicate.Mean.perRep.BestDup.lg10 == CV.duplicate.Mean.perRep.lg10, "Best-combo","Poor-combo"))%>%
  ungroup()%>%
  #select(-Well_coordinate)%>%
  distinct()


#### Merging triplicate and duplicate CV ---------------------------------------------------------------------------------------------------------------------------------------------------------
KF.OL.3cv.2cv.Mean.perRep <- left_join(KF.OL.3cv.Mean.perRep,
                           KF.OL.2cv.perRep.Mean.perRep)


Check.CV <- KF.OL.CV.df

#Houskeeping
rm(KF.OL.3cv.Mean.perRep,
    KF.OL.2cv.perRep.Mean.perRep)


### Distribution of coeffient of variation: 3cv and 2cv  ----------------------------------------------------------------------------------------------------------------------------------------------------------
Outlier.df <- KF.OL.3cv.2cv.Mean.perRep %>%
  ungroup()%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Def,
         CV.triplicate.Mean.perRep,
         CV.duplicate.Mean.perRep.BestDup)%>%
  ungroup()%>%
  distinct()%>%
  drop_na() # SOLUTION HERE 20230618


Outlier.df <- melt(Outlier.df,
                                       id.vars = c("ExpFile",
                                                   "Abx.con",
                                                   "Isolate",
                                                   "Killing.Def"))
  Outlier.df <- Outlier.df %>%
  group_by(Killing.Def,
           variable)%>%
  mutate(Mean3cv = mean(value))%>%
  mutate(Mean.3stdev =3*sd(value) + Mean3cv)
  
  
  CV.duplicate.thrs <- Outlier.df %>%
  ungroup()%>%
  filter(variable == "CV.duplicate.Mean.perRep.BestDup" & Killing.Def ==  "Trkv2_LCF_Ila2")%>%
      #filter(variable == "CV.duplicate.Mean.perRep.BestDup" & Killing.Def ==  "P.Ila_LCF_scvsBaSic")%>%

  select(Mean.3stdev)%>%
  distinct()

CV.duplicate.thrs <- CV.duplicate.thrs$Mean.3stdev



CV.triplicate.thrs <-Outlier.df %>%
  ungroup()%>%
  filter(variable == "CV.triplicate.Mean.perRep" & Killing.Def ==  "Trkv2_LCF_Ila2")%>%
  select(Mean.3stdev)%>%
  distinct()

CV.triplicate.thrs <- CV.triplicate.thrs$Mean.3stdev


  
library(ggprism)
  plot.path <-  paste(exp.resDir,
                       "/",
                       unique(All.iso$ExpFile),
                       "_",
                       old.or.new.tracking,
                       "_Reproducibility_Cov.pdf",sep="")



pdf(  plot.path)
  Outlier.df %>%
    filter(Killing.Def == "Trkv2_LCF_Ila2")%>%
  ggplot(aes(x = value,
           colour = variable,
           fill = variable))+
    
     geom_histogram(aes(y = ..density..),
             #   colour = 1,
                 alpha = 0.25) +
  geom_density(lwd = 1,  alpha = 0.25)+
      xlim(0, 0.75) +

  # geom_density(bw = 0.1,
  #              size = 1)+
  #     geom_histogram(binwidth = 0.001, fill = "lightblue", 
  #                    bins = 50) +

  xlab("Coefficient of variation [a.u]")+
  ylab("Frequency [a.u]")+
  ggtitle(paste("Histogram CoV",
                zb.exp,
                unique(Outlier.df$Abx.con),sep=" "))+
    labs( caption = paste("Mean +3sd : Triplicate ", signif(CV.triplicate.thrs, digits = 4) , " | Duplicate:", signif(CV.duplicate.thrs, digits = 4) ,sep="") )+

  #theme_prism()+
   # theme_pubr()+
#  theme(text = element_text(size =4))  +
  geom_vline(aes(xintercept = Mean3cv),
               colour = "Black")+
  geom_vline(aes(xintercept = Mean.3stdev),
               colour = "Grey")+
    theme(
    legend.position = "bottom",  # Move legend to the bottom
    legend.box = "horizontal" # Display legend items horizontally
    
    ) +

  theme(aspect.ratio =1)+

  facet_wrap(Killing.Def~variable)
dev.off()

#-- WORK HERE !!!! ----



Outlier.Data.export <- Outlier.df %>%
  filter(Killing.Def == "Trkv2_LCF_Ila2")%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Def,
         variable,
         value,
         Mean3cv,
         Mean.3stdev)%>%
  distinct()

write.csv(Outlier.Data.export,file = paste(exp.resDir,
                                           "/",
                                           unique(Outlier.Data.export$ExpFile),
                                           "_",
                                           old.or.new.tracking,
                                           "_",
                                           analysis.normalisation,
                                           "_Outlier_CoV.Values.csv",
                                           sep=""),
          row.names = FALSE)



Outlier.Data.export <- Outlier.df %>%
  filter(Killing.Def == "Trkv2_LCF_Ila2")%>%
  select(ExpFile,
         Abx.con,
         #    Isolate,
         Killing.Def,
         variable,
         #    value,
         Mean3cv,
         Mean.3stdev)%>%
  distinct()

write.csv(Outlier.Data.export,file = paste(exp.resDir,
                                           "/",
                                           unique(Outlier.Data.export$ExpFile),
                                           "_",
                                           old.or.new.tracking,
                                           "_",
                                           analysis.normalisation,
                                           "_Outlier_CoV.Values_Mean.Stdev.csv",
                                           sep=""),
          row.names = FALSE)


rm(Outlier.Data.export)

### Flagging outliers based on CoV ----------------------------------------------------------------------------------------------------------------------------------------------------------
# 0.1395
Outlier.df.intermed.tracking <- KF.OL.3cv.2cv.Mean.perRep %>%
  ungroup()%>%
  distinct() %>%
  filter(Killing.Def == "Trkv2_LCF_Ila2")%>%
  group_by(ExpFile,
           Abx.con,
           Isolate)%>%
  
  mutate(Dup.Cal = (100/CV.triplicate.Mean.perRep)*CV.duplicate.Mean.perRep.BestDup)%>% # How much lower is the duplicate cv compared to the triplicate cv
  mutate(Dup.cal.test  =  if_else( Dup.Cal <= 50, "Keep-Best-Duplicate" ,"consider-triplcate"))%>%
  
  mutate(Outlier.Evaluation = if_else(signif(CV.triplicate.Mean.perRep, digits = 2)  > CV.triplicate.thrs  & signif(CV.duplicate.Mean.perRep.BestDup, digits = 2)  <= CV.duplicate.thrs , "Keep-Best-Duplicate" ,
                                      if_else( signif(CV.triplicate.Mean.perRep, digits= 2) <= CV.triplicate.thrs , "Keep-Triplicate", "Omit-Isolate")))





Outlier.df <- KF.OL.3cv.2cv.Mean.perRep %>%
  left_join(  Outlier.df.intermed.tracking)

# Build a dataframe of all the Well corrdiantes to kee and the ones flagges as outliers
list.of.outlier.evaluation <-c ("Keep-Triplicate", "Keep-Best-Duplicate" )

Outlier.wells <- data.frame()

for ( i in list.of.outlier.evaluation ) {
  print(i)
  if ( i == "Keep-Triplicate" ) {
    
    print("here")
     Outlier.df.subset <-   Outlier.df %>%
    filter(Outlier.Evaluation  == "Keep-Triplicate")%>%
    select(ExpID,
      #     Exp.Iso,
           Well_coordinate,
           Abx.con,
           Isolate)%>%
  distinct()%>%
    mutate(Outlier.Evaluation = "Reproducible")
     
     Outlier.wells <- rbind(    Outlier.wells  ,
                                 Outlier.df.subset)
  }
  
  
  if ( i == "Keep-Best-Duplicate" ) {
    
        print("here Dups")
     Outlier.df.subset <-   Outlier.df %>%
    filter(Outlier.Evaluation  == "Keep-Best-Duplicate")%>%
       ungroup()%>%
    select(ExpID,
           Well_coordinate,
           Abx.con,
           Isolate,
           Outlier.Evaluation,
           Combo,
           Best.combo.Dup)%>%
    distinct()%>%
    filter(Best.combo.Dup == "Best-combo")%>%
    mutate(Outlier.Evaluation = "Reproducible")%>%
       select(-Outlier.Evaluation)
     
    Well_coordinate <- unlist(paste(Outlier.df.subset$Combo,collapse="_"))
    Well_coordinate <- strsplit( Well_coordinate , split = "_")
    Well_coordinate <- Well_coordinate[[1]]
    Well_coordinate <-  unique(Well_coordinate)
      
  
     Outlier.df.subset.subset <- data.frame( Well_coordinate =  Well_coordinate,
                                      Outlier.Evaluation = "Reproducible"
                                      )
       
          Outlier.df.subset <- left_join(   Outlier.df.subset,   Outlier.df.subset.subset , by = ("Well_coordinate"))
          
           Outlier.df.subset <-  Outlier.df.subset %>%
             drop_na()%>%
               select(ExpID,
           Well_coordinate,
           Abx.con,
           Isolate,
            Outlier.Evaluation)%>%
  distinct()
             
     Outlier.wells <- rbind(    Outlier.wells  ,
                                 Outlier.df.subset)
     
     #Houskepeing
     rm(Outlier.df.subset)
  }
 
  
 

  
}
 





## SECTION 18: Merging time kill curve data with outlier (CoV) data ------------------------------------------------------------------------------------------------------------------------------------------------------

  exp.name <- unique(perWell.df$ExpFile)

  Outlier.wells <- Outlier.wells %>%
  mutate(ExpFile = exp.name)

  Pop.TKC.corr.OL1.df  <- T.TKC.corr 

  Pop.TKC.corr.OL2.df <- left_join(Pop.TKC.corr.OL1.df, 
                                   Outlier.wells, by = c("ExpFile","Well_coordinate","Isolate","Abx.con"))


  Pop.TKC.corr.OL2.df <- Pop.TKC.corr.OL2.df %>% 
  replace_na(list(Outlier.Evaluation = 'Outlier'))%>%
    mutate(ExpID = ExpFile)
  
  Pop.TKC.corr.OL2.df$Outlier.Evaluation%>% replace_na('Outlier')
  
rm(CV.combo.res,
   Check.CV,
   KF.OL.3cv.2cv.Mean.perRep,
   KF.OL.CV.df,
   Outlier.df.subset.subset,
   Outlier.df,
   Pop.TKC.corr.OL1.df,
   i,
   i.Well_rep1,
   i.Well_rep2,
   list.of.duplicates,
   list.of.drugs,
   list.of.outlier.evaluation,
   loop.list,
   n,
   x)


# SECTION 19: Visualising tracking data ----
## Plot  heatmap of wells with less than 1000 cells ----

QC.plot <- QC.T %>%
  mutate( Contamination  = if_else(Contamination == "No-Contamination", "Keep", "Omit"),
          Growth = if_else(Growth == "No-Growth", "Keep", "Omit"),
          Artefact = if_else(Artefact == "No-Artefact", "Keep", "Omit"),
          A.Thousand.cells.tracked = if_else(T_Eval_Trck.Numbers == ">1K", "Keep", "Omit"))%>%
  select(-T_Eval_Trck.Numbers)


QC.plot <- melt(QC.plot,
           id.vars = c("ExpFile",
                                  "Well_coordinate",
                       "Isolate"
             ),
                      variable.name =  "QC.readouts")

gg.experiment.title <- paste("Threshold of 1k cells", expID, experiment.drug, sep =" ")

pdf(paste(exp.resDir,"/",unique(perWell.Trk.df$ExpFile),
          "_",
          old.or.new.tracking,

          "_Trk.N_QC.readouts.pdf",
          sep=""))

QC.plot%>%
    ungroup()%>%
      filter(QC.readouts == "A.Thousand.cells.tracked" )%>%

  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  mutate(Exp.Colunm = paste(ExpFile,
                            Colunm,
                            sep="_"))%>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  mutate(Exp.Well = paste(ExpFile,
                          Well_coordinate,
                          sep="_"))%>%
  ggplot(aes(x = Colunm,
             y= Row,
             label = Well_coordinate))+
  geom_tile(aes(fill = value),width=0.7, height=0.7)+
  scale_fill_manual(values = c("SteelBlue","coral2"))+
  geom_text(aes(label = Well_coordinate), color = "black", size = 0.3, nudge_y = 0.1)+
  geom_text(aes(label = paste(Isolate,sep=".")), color = "white", size = 0.2, nudge_y = -0.2 )+
  coord_equal()+
  theme_bw()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+
   scale_y_discrete(limits=rev)+
   scale_x_continuous( breaks = seq(type.of.wellplate),
                      limits=c(0,NA))+
   labs(title = gg.experiment.title)+
      theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
 theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))+

   facet_wrap(~QC.readouts,
               ncol = 1,
              nrow = 1)


dev.off()

rm(gg.experiment.title)




### Plot heatmap of number of tracked cells  ----------------------------------------------------------------------------------------------------------------------------------------------------------
  QC.T.n <- QC.T %>%
  ungroup()%>%
  select(ExpFile,
            Well_coordinate,
           Isolate,
           T_Tot.N)%>%
  mutate(value = T_Tot.N)%>%
  distinct()%>%
  select(ExpFile,
         Well_coordinate,
         Isolate,
        value)


  gg.experiment.title <- paste("Heatmap of cell numbers from tracking", expID, experiment.drug, sep =" ")


 pdf(paste(exp.resDir,"/",unique(  QC.T.n$ExpFile),
            "_",
          old.or.new.tracking,
          "_T_QC.MOCn.pdf",sep=""))
   QC.T.n %>%
   
  mutate(Colunm = as.numeric(gsub("[^0-9.-]", "", Well_coordinate)))%>%
  mutate(Exp.Colunm = paste(ExpFile,
                            Colunm,
                            sep="_"))%>%
  mutate(Row = gsub('[[:digit:]]+', '', Well_coordinate))%>%
  mutate(Row= as.factor(Row))%>%
  mutate(Exp.Well = paste(ExpFile,
                          Well_coordinate,
                          sep="_"))%>%
  ggplot(aes(x = Colunm,
             y= Row,
             label = Well_coordinate))+
  geom_tile(aes(fill = value),width=0.7, height=0.7,color="black")+
  geom_text(aes(label = Well_coordinate), color = "black", size = 0.2, nudge_y = 0.2)+
  geom_text(aes(label = value), color = "black", size = 0.2 )+
    geom_text(aes(label = paste(Isolate,sep=".")), color = "black", size = 0.2, nudge_y = -0.2 )+
  coord_equal()+
  theme_bw()+
  theme(aspect.ratio = 0.5)+
  theme(legend.key.size = unit(0.3, "cm"),
  legend.key.width = unit(0.3,"cm"))+
  theme(legend.position = "top")+
  scale_y_discrete(limits=rev)+
    ggtitle(  gg.experiment.title)+
  scale_x_continuous( breaks = seq(type.of.wellplate ),
                      limits=c(0,NA))+
  scale_fill_distiller(type = "seq", palette = "RdYlBu",
                       limits = c( min(   QC.T.n$value), max(   QC.T.n$value)))+
  theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
 theme(axis.text.y = element_text(size = 3),
         legend.text=element_text(size=2))+
     labs(fill = "Total number of tracked cells")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 3))

 dev.off()


 rm(  QC.T.n,
      QC.plot)


### Plot time kill curve and flag labels ---------------------------------------------------------------------------------------------------------------------------------------------------------
All.iso.T.df <-Pop.TKC.corr.OL2.df



Plot.label.df <- All.iso.T.df%>%
  filter(Time.Kill.Definitions!="T_Tot.N")%>%
  mutate( Contamination  = if_else(Contamination == "No-Contamination", "Keep", "Omit"),
          Growth = if_else(Growth == "No-Growth", "Keep", "Omit"),
           Fluorescence.Artefact = if_else(Artefact == "No-Artefact", "Keep", "Omit"),
          `1000 Single Cell and Vsnapp population analysis` = if_else(T_Eval_Trck.Numbers == ">1K", "Keep", "Omit"))%>%
  mutate(Growth = if_else(Growth == "Keep", "", "Growth"),
         Fluorescence.Artefact = if_else(Fluorescence.Artefact == "Keep", "", "AF"),
         `1000 Single Cell and Vsnapp population analysis` = if_else( `1000 Single Cell and Vsnapp population analysis` == "Keep", "", "<1K"),
         Contamination = if_else(Contamination == "Keep", "", "Contam"))%>%
    mutate(Outlier.Evaluation =  if_else(Outlier.Evaluation == "Reproducible", "","OT"))




Plot.label.df.subset <- Plot.label.df %>%
  select(ExpFile,
         Abx.con,
         Well_coordinate,
         Isolate,
         Outlier.Evaluation) %>%
  ungroup()%>%
  mutate(Outlier.Evaluation =  if_else(Outlier.Evaluation == "Reproducible", "","OT"))


Plot.label.df <- Plot.label.df %>%
    # mutate(Overall.Assessment = if_else(Contamination == "" & Growth =="" & Fluorescence.Artefact=="" & `1000 Single Cell and Vsnapp population analysis` == "" & Outlier.Evaluation == "", "Keep", "Omit"))%>%

     mutate(Overall.Assessment = if_else(Contamination == "" & Growth =="" &  `1000 Single Cell and Vsnapp population analysis` == "" & Outlier.Evaluation == "", "Keep", "Omit"))%>%

  mutate(Flag.label = paste(Growth,
                             `1000 Single Cell and Vsnapp population analysis`,
                             Contamination ,
                           Outlier.Evaluation,
                            sep=" "))%>%

  select(ExpFile,
         Abx.con,
         Well_coordinate,
         Isolate,
         timestep,
         Time_Hrs,
         LC.fraction.corr,
         Flag.label,
         Time.Kill.Definitions,
         Overall.Assessment)%>%
  distinct()

Plot.label.df <- Plot.label.df %>%
  mutate(Flag.label = str_replace_all(Flag.label, "\\s+", " "))%>%
  mutate(Flag.label = gsub("Reproducible", "", Flag.label))

export.filename <-  paste(exp.resDir,"/",exp.name,
                    "_",
          old.or.new.tracking,
          "_Flags.csv", sep = "")
write.csv(Plot.label.df, export.filename,
          row.names = FALSE)
rm(Plot.label.df.subset ,
   export.filename)
  


### Ilastik twice positive ------------------------------------------------------------------------------------------------------------------------------------------------------

# PLOTTING TIME KILL CURVES
LiveCells <- paste(exp.resDir,"/",exp.name,
                    "_",
          old.or.new.tracking,
          "_TimeKillCurves_Experiment_Flags.pdf", sep = "")


library(ggrepel)
gg <- Plot.label.df  %>%
  filter(Time.Kill.Definitions != "T_Tot.N" )%>%
 filter(Time.Kill.Definitions == paste(old.or.new.tracking,"_Ila2",sep=""))%>%

  distinct()%>%
  mutate(Overall.Assessment = as.factor(Overall.Assessment))%>%

  ggplot(aes(x = Time_Hrs,
             y = log10(LC.fraction.corr),
             group = Well_coordinate,
             label = Well_coordinate,
             colour = Overall.Assessment))+
  geom_line( #colour = "black",
             alpha = 0.5,
             size = 0.2)+
geom_label_repel(data = . %>% group_by(Well_coordinate) %>% filter(timestep == 3) %>% #slice_tail(n = 1) %>%
                      distinct(),
                  aes(label = Flag.label),
                 fill = "transparent",
                 #color = "transparent",
                   label.size = NA, label.padding = unit(0, "lines"), label.color = NA,
                 size = 0.2,
                  nudge_y = -0.3,
                   segment.color = "transparent", point.color = NA, arrow = arrow(length = unit(0.03, "npc"))) +


  theme_bw()+

   scale_y_continuous(

                    breaks = c(0, -0.999,  -1.999),
                    limits = c(-3,0),
                     labels = c(100, 10, 1))+
    scale_colour_manual(values=c(Keep="steelblue",
                                 Omit="red"))+
    geom_text(data = .  %>% group_by(Well_coordinate) %>% top_n(n = 1, wt = Time_Hrs),
      aes(label = Well_coordinate), size = 0.2)+

    theme(axis.text = element_text(size = 3))     +
    theme(strip.text.x = element_text(size = 3))+
    theme(legend.key.size = unit(0.1, "cm"),
    legend.key.width = unit(0.1,"cm"))+
    scale_x_continuous( breaks = plot.xaxis.time.scale)+
    labs(title = paste(experiment.id,
                       unique(All.iso.df$Abx.con),
                           sep=" "),

       x = "Time [Hrs]",
       y = "Percentage of cells alive [%]",
       color = "Condition evaluation",
       subtitle = "Contam = Contamination, Artefacts = AF, Gel-detachment = GD, Less 1K sc & vs = <1K, Outlier = OT")+
    geom_text(x=48, y=1,
            size = 1)+ # Plot label
  theme(legend.position="top")+
   theme(aspect.ratio = 1)+
   theme(plot.subtitle = element_text(size = 5))+

  facet_wrap_paginate(Abx.con+Time.Kill.Definitions~Isolate,
                      ncol =12,
                      nrow =12,
                      page =1)
n <- n_pages(gg)

pdf(LiveCells ,paper = "a4", width = 20 , height = 15 )
for(i in 1:n){
    print(gg + facet_wrap_paginate(Abx.con+Time.Kill.Definitions~Isolate,
                      ncol =12,
                      nrow =12, page = i))
}
dev.off()

#Houskeeping
rm(gg,
   LiveCells,
   g)



All.Iso.Analysis.KEEP <- All.Iso.Analysis


## SECTION 20: merging tracking killing features with experiment QC --------------------------------------------------------------------------------------------------------------------------------------------------------
# REMOVE LCF in the middle here  
  All.Iso.Analysis.Tracking <- Main.KF.df.TRACKING %>%
  ungroup()%>%
  mutate(Merge.def = paste(Analysis.Str,
                           cell.type,
                           sep="_"))%>%
  select(ExpFile,
         Abx.con,
         Well_coordinate,
         Isolate,
         Merge.def,
         Main.KF,
         Killing.Def,
         Analysis.Str,
         Killing.Features,
         value)%>%
  filter(Analysis.Str == old.or.new.tracking )
  
  Quality.Eval <- Plot.label.df %>%
  ungroup()%>%
  select(ExpFile,
         Abx.con,
         Well_coordinate,
         Flag.label,
         Overall.Assessment)%>%
  distinct()
 
  All.Iso.Analysis.Tracking <- left_join(All.Iso.Analysis.Tracking,
                                           Quality.Eval)


    #Houskeeping
    rm(Quality.Eval)

# INTERMEDIARY ----
##### 20.2.1 Instance of different events in across two replicates and ----

    
    Average.Result.perIso.PASSED.all.Criteria.INTERMEDIARY.tracking <- All.Iso.Analysis.Tracking %>%
      select(ExpFile,
             Abx.con,
             Isolate,
             Well_coordinate,
             Overall.Assessment,
             Flag.label,
             Merge.def)%>%
      ungroup()%>%
      distinct()%>%
      mutate(Outlier.Evaluation = ifelse(grepl("OT", Flag.label), "Outlier", "Reproducible"))%>%
      mutate(Multiple.Events.intermediary = if_else(Overall.Assessment == "Keep" & Outlier.Evaluation == "Reproducible", 0 , 1))%>%
            group_by(ExpFile,
               Isolate,
               Merge.def)%>%
      mutate(Multiple.Events = sum(Multiple.Events.intermediary),
             Multiple.Events.Check = if_else(Multiple.Events > 1, "Omit", "Keep"))%>%
      mutate(Multiple.Events  = Multiple.Events.Check)%>%
      ungroup()%>%
      select(ExpFile,
             Abx.con,
             Isolate,
             Well_coordinate,
             Merge.def,
             Multiple.Events)%>%
      distinct()
    
    
    
    
    
    
    filename <- paste(exp.resDir,"/",exp.name,
                      "_",
                      old.or.new.tracking,
                      "_KFintermediary_.csv",
                      sep = "")
    
    
# Exporting raw KF tracking  ----
    
    
    filename.rawKF <-  paste(exp.resDir,"/",exp.name,
          "_",
          old.or.new.tracking,
          "_rawKF.csv",
          sep = "")
    
    
    kf.export <- Main.KF.df.TRACKING %>%
      filter(cell.type == "Ila2BaSic")%>%
      filter( Main.KF != "MDK")
    
    # Saving ASCT Clinical isoalte data
    write.csv(  kf.export,
                filename.rawKF,
                row.names=FALSE)
    
    rm(  kf.export,
         filename.rawKF )
    
    
    # Saving ASCT Clinical isoalte data
    write.csv(Average.Result.perIso.PASSED.all.Criteria.INTERMEDIARY.tracking ,
              filename,
              row.names=FALSE)
    
    
    
    
    Average.Result.perIso.PASSED.all.Criteria <-    All.Iso.Analysis.Tracking %>%
      ungroup()%>%
      left_join(Average.Result.perIso.PASSED.all.Criteria.INTERMEDIARY.tracking)%>%
      filter(Overall.Assessment == "Keep") %>%
      filter(  Multiple.Events == "Keep") %>%
      
      mutate( value = if_else(value == ">72", "72",value))%>%
      mutate(value = as.numeric(value))
    
    
  
  Average.Result.perIso.PASSED.all.Criteria <- Average.Result.perIso.PASSED.all.Criteria %>%
  ungroup()%>%
  group_by(ExpFile,
             Isolate,
             Abx.con,
             Main.KF,
             Killing.Def,
             Killing.Features) %>%
  summarise(Avg_Killing_Features = mean(value))%>%
  ungroup()%>%
  distinct()%>%
  mutate(Isolate = paste("Iso.",
                         Isolate,
                         sep=""))%>%
  select(-Killing.Def,
         -Main.KF)


  Variable.defintion <- All.Iso.Analysis %>%
  ungroup()%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Def,
         Main.KF,
         Killing.Features)

  All.Iso.Analysis<- All.Iso.Analysis %>% 
     filter(grepl(old.or.new.tracking, Killing.Def ))%>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Features)



  All.Iso.Analysis <- left_join(All.Iso.Analysis ,
                          Average.Result.perIso.PASSED.all.Criteria,
                         by = join_by(ExpFile, Abx.con, Isolate, Killing.Features))


All.Iso.Analysis <- left_join(All.Iso.Analysis,
                                   Variable.defintion )

# need to add back the Killing Def, Main.KF and Killing Features

All.Iso.Analysis <- All.Iso.Analysis %>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Main.KF,
         Killing.Def,
         Killing.Features,
         Avg_Killing_Features)

# Housekeeping
rm(Variable.defintion)




All.iso <- All.Iso.Analysis  %>%
  select(ExpFile,
         Abx.con,
         Isolate,
         Killing.Def,
         Main.KF,
         Killing.Features)



All.Iso.Analysis  <- All.Iso.Analysis  %>%
  ungroup()%>%
  drop_na()%>%
  ungroup()%>%
  mutate(Killing.Features = as.factor(Killing.Features))%>%
  mutate(Avg_Killing_Features = Avg_Killing_Features)%>%
  mutate(Avg_Killing_Featureslg10 = log10(Avg_Killing_Features))%>%
  group_by(ExpFile,
           Abx.con,
           Killing.Features # IMPORTANT FOR TRACKING
           )%>%
  mutate(Min.Per.Feature.lg10 = min(Avg_Killing_Featureslg10))%>%
  mutate(Max.Per.Feature.lg10 = max(Avg_Killing_Featureslg10))%>%
  mutate(Norm_Avg_Killing_Features.lg10 = (Avg_Killing_Featureslg10 - Min.Per.Feature.lg10) / (Max.Per.Feature.lg10 -Min.Per.Feature.lg10) )%>%
    mutate(Min.Per.Feature = min(Avg_Killing_Features))%>%
  mutate(Max.Per.Feature = max(Avg_Killing_Features))%>%
  mutate(Norm_Avg_Killing_Features = (Avg_Killing_Features - Min.Per.Feature) / (Max.Per.Feature -Min.Per.Feature) )%>%
## Here make two sets of calculations a log10 transfromed average and a raw average
  ungroup()%>%
    select(-Min.Per.Feature,
         -Max.Per.Feature,
         -Min.Per.Feature.lg10,
         -Max.Per.Feature.lg10)

# Previously I saw a NaN error when caluclating MDK90%. In cases when this is never reached across all isolates, when we normalize b (x-xmin)/(xmin-xmax) we will get a final calculation zero / zero. This will throw a NaN error in R. To avoid this, I made an if_else statement which will consider For all MDK values across each isolate, if the Norm_Avg_Killing_Features.lg10 is a NaN value, it will be changed to 1 otherwise it remains the value it was defined.
All.Iso.Analysis  <- All.Iso.Analysis %>%
  ungroup()%>%
  group_by(ExpFile,
          #Killing.Def,
           Abx.con,
           Isolate)%>%
  mutate(Norm_Avg_Killing_Features.lg10 =  if_else(is.nan(Norm_Avg_Killing_Features.lg10) & Main.KF == "MDK" ,
                                                   1,
                                                   Norm_Avg_Killing_Features.lg10))%>%
  mutate(Norm_Avg_Killing_Features =  if_else(is.nan(Norm_Avg_Killing_Features) & Main.KF == "MDK" ,
                                                   72,
                                                   Norm_Avg_Killing_Features))
  


All.Iso.Analysis <- left_join(All.iso,
                              All.Iso.Analysis)




### Reordering labels  ------------------------------------------------------------------------------------------------------------------------------------------------------

Killing.features.unique <- unique(All.Iso.Analysis$Killing.Features)


# Define a time-kill curve definition
All.Iso.Analysis <- All.Iso.Analysis %>%
  mutate(Time.Kill.curve.defintion = Killing.Def)



All.Iso.Analysis <- separate(All.Iso.Analysis, col =Time.Kill.curve.defintion,
                                  into = c("Analysis",
                                           "Feature",
                                           "Cell.type"),
                                  sep="_")

All.Iso.Analysis<- All.Iso.Analysis%>%
  mutate(Killing.Def  = paste(Analysis,
                              Feature,
                              sep="_"))%>%
  select(-Analysis,
         -Feature,
         -Cell.type)




# SECTION 21:  Exporting tracking results ----
## Defining order of killing feautre variables ----


rm(general.order)
expID <- unique( Pop.TKC.corr$ExpFile)

    library(stringr)

# Export log10 transformed and raw values

    All.Iso.Analysis.WIDER.subset <- All.Iso.Analysis %>%
      ungroup()%>%
      mutate(Avg_Killing_Features = as.character(Avg_Killing_Features),
             Avg_Killing_Features = if_else( Avg_Killing_Features == "72" & Main.KF == "MDK", ">72", Avg_Killing_Features))
    
    # -----Order of colunms
    
    
  time.kill.def.part_I <- (old.or.new.tracking)
  time.kill.def.part_II <- ("Ila1")


  order.x.axis.LC.fraction <- c("LCF_3h",
                              "LCF_6h",
                              "LCF_9h",
                              "LCF_12h",
                              "LCF_24h",
                              "LCF_36h",
                              "LCF_48h",
                              "LCF_60h",
                              "LCF_72h")

  order.x.axis.AUC.realtime <- c("AUCrt_0.3h",
                               "AUCrt_0.6h",
                               "AUCrt_0.9h",
                               "AUCrt_0.24h",
                               "AUCrt_0.36h",
                               "AUCrt_0.48h",
                               "AUCrt_0.60h",
                               "AUCrt_0.72h",
                                "AUCrt_48.72h")
  
  order.x.axis.AUC.realtimeLOG <- c("AUCrtLOG_0.3h",
                               "AUCrtLOG_0.6h",
                               "AUCrtLOG_0.9h",
                               "AUCrtLOG_0.24h",
                               "AUCrtLOG_0.36h",
                               "AUCrtLOG_0.48h",
                               "AUCrtLOG_0.60h",
                               "AUCrtLOG_0.72h",
                                "AUCrtLOG_48.72h")
  

   

  order.x.axis.MDK <- c("MDK_25pct",
                      "MDK_50pct",
                      "MDK_75pct",
                      "MDK_90pct")

  order.x.axis.overall.order <- c(order.x.axis.LC.fraction,
                               order.x.axis.AUC.realtime,
                               order.x.axis.AUC.realtimeLOG,
                               order.x.axis.MDK)
  


# Changing colunm nammes
  order.x.axis.overall.order.T_Ila1 <- paste(time.kill.def.part_I  ,
                                    order.x.axis.overall.order,
                                    time.kill.def.part_II,
                                    sep="_")
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Ila1BaSic <- gsub(  "Ila1",
                                                "Ila1BaSic",
                                                 order.x.axis.overall.order.T_Ila1)
  
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Ila2BaSic <- gsub(
                                                "Ila1BaSic",
                                                 "Ila2BaSic",
                                                  order.x.axis.overall.order.T_Ila1BaSic)
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Ila2 <- gsub(
                                                "Ila1",
                                                 "Ila2",
                                                    order.x.axis.overall.order.T_Ila1 )
  
  
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Thr1BaSic <- gsub(   "Ila1BaSic",
                                                "Thr1BaSic",
                                                  order.x.axis.overall.order.T_Ila1BaSic)
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Thr2BaSic <- gsub(  
                                                "Thr1BaSic",
                                                  "Thr2BaSic",
                                                  order.x.axis.overall.order.T_Thr1BaSic )
  
  
    
# Changing colunm nammes
  order.x.axis.overall.order.T_Thr1 <- gsub(            "Thr1BaSic",
                                                "Thr1",
                                                    order.x.axis.overall.order.T_Thr1BaSic)
  
# Changing colunm nammes
  # Changing colunm nammes
  order.x.axis.overall.order.T_Thr2 <- gsub(            "Thr1BaSic",
                                                "Thr2",
                                                    order.x.axis.overall.order.T_Thr1BaSic)
  
  
  
  if ( old.or.new.tracking == "Trkv1") {
    
      
   # Ordering all time kill cuvre defitions
  order.x.axis.overall.order.T_piDyn <- gsub(  "Ila1",
                                                "piDyn",
                                                 order.x.axis.overall.order.T_Ila1)
  
      general.order <- c(
  order.x.axis.overall.order.T_piDyn,
  order.x.axis.overall.order.T_Ila1,
  order.x.axis.overall.order.T_Ila2,
  order.x.axis.overall.order.T_Ila1BaSic,
  order.x.axis.overall.order.T_Ila2BaSic,
  order.x.axis.overall.order.T_Thr1BaSic,
  order.x.axis.overall.order.T_Thr2BaSic,
  order.x.axis.overall.order.T_Thr1,
  order.x.axis.overall.order.T_Thr2
  )
  
    
  } else (      general.order <- c(
  order.x.axis.overall.order.T_Ila1,
  order.x.axis.overall.order.T_Ila2,
  order.x.axis.overall.order.T_Ila1BaSic,
  order.x.axis.overall.order.T_Ila2BaSic,
  order.x.axis.overall.order.T_Thr1BaSic,
  order.x.axis.overall.order.T_Thr2BaSic,
  order.x.axis.overall.order.T_Thr1,
  order.x.axis.overall.order.T_Thr2
  )
  )
  


      general.order <- paste(gsub("_",".",experiment.drug),
                             general.order,
                             experiment.id,
                             sep="_")
   
      general.order <- c(#"ExpFile",
                         #"Abx.con",
                         "Isolate",
                       #  "Killing.Def",
                         general.order)
     
      
         All.Iso.Analysis.WIDER.subset <- All.Iso.Analysis.WIDER.subset %>%
      filter(grepl(old.or.new.tracking, Killing.Def)) %>%
          select(ExpFile,
                 Abx.con,
                 Isolate,
                 Killing.Features,
                 Avg_Killing_Features,
                  Avg_Killing_Featureslg10)
         
             All.Iso.Analysis.WIDER.subset <-   All.Iso.Analysis.WIDER.subset %>%
      mutate(Abx.con = gsub("_",".",Abx.con))
             
             
             All.Iso.Analysis.WIDER.subset <-   All.Iso.Analysis.WIDER.subset %>%
               mutate(Killing.Features = paste(Abx.con,
                                               Killing.Features,
                                               ExpFile,
                                               sep="_"))
               
            # grabbing the lg10 data need to add reference to heading   
           All.Iso.Analysis.WIDER.subset.lg10 <-  All.Iso.Analysis.WIDER.subset %>%
                   select(
                 Isolate,
                 Killing.Features,
                  Avg_Killing_Featureslg10)
               
            All.Iso.Analysis.WIDER.subset <- All.Iso.Analysis.WIDER.subset %>%
               select(Isolate,
                      Killing.Features,
                      Avg_Killing_Features)
    
  All.Iso.Analysis.WIDER.subset <- pivot_wider( All.Iso.Analysis.WIDER.subset,
                                                  names_from = Killing.Features,
                                                values_from = Avg_Killing_Features)
  
 
  
  All.Iso.Analysis.WIDER.subset.Organised.KF <-   All.Iso.Analysis.WIDER.subset[,general.order ]
  
  
  
  


#### A check if variables are missing ----------------------------------------------------------------------------------------------------------------------------------------------------------
# List of column names from general.order
columns_to_check <- general.order

# Get the column names that do not exist in the dataset
missing_columns <- setdiff(columns_to_check, colnames(All.Iso.Analysis.WIDER.subset))

# Print the missing columns
print(missing_columns)



### Defining order of log10 variables  ----------------------------------------------------------------------------------------------------------------------------------------------------
rm(All.Iso.Analysis.WIDER.subset)
  
# ---- Generating output of log10 data
  
           All.Iso.Analysis.WIDER.subset.lg10 <- separate(   All.Iso.Analysis.WIDER.subset.lg10 , col = Killing.Features,
                                                             into = c("c1",
                                                                      "c2",
                                                                      "c3",
                                                                      "c4",
                                                                      "c5",
                                                                      "c6"),
                                                             sep = "_")
           
           
               All.Iso.Analysis.WIDER.subset.lg10 <- All.Iso.Analysis.WIDER.subset.lg10 %>%   
                 mutate(Killing.Features = paste(c1,
                                                 c2,
                                                 c3,
                                                 paste(c4,
                                                       ".lg10",
                                                       sep=""),
                                                 c5,
                                                 c6,
                                                 sep="_"))%>%
                 select(Isolate,
                        Killing.Features,
                        Avg_Killing_Featureslg10)
  
   
  time.kill.def.part_I <- (old.or.new.tracking)
  time.kill.def.part_II <- ("Ila1")


  order.x.axis.LC.fraction <- c("LCF_3h.lg10",
                              "LCF_6h.lg10",
                              "LCF_9h.lg10",
                              "LCF_12h.lg10",
                              "LCF_24h.lg10",
                              "LCF_36h.lg10",
                              "LCF_48h.lg10",
                              "LCF_60h.lg10",
                              "LCF_72h.lg10")

  order.x.axis.AUC.realtime <- c("AUCrt_0.3h.lg10",
                               "AUCrt_0.6h.lg10",
                               "AUCrt_0.9h.lg10",
                               "AUCrt_0.24h.lg10",
                               "AUCrt_0.36h.lg10",
                               "AUCrt_0.48h.lg10",
                               "AUCrt_0.60h.lg10",
                               "AUCrt_0.72h.lg10",
                                "AUCrt_48.72h.lg10")
  
  order.x.axis.AUC.realtimeLOG <- c("AUCrtLOG_0.3h.lg10",
                               "AUCrtLOG_0.6h.lg10",
                               "AUCrtLOG_0.9h.lg10",
                               "AUCrtLOG_0.24h.lg10",
                               "AUCrtLOG_0.36h.lg10",
                               "AUCrtLOG_0.48h.lg10",
                               "AUCrtLOG_0.60h.lg10",
                               "AUCrtLOG_0.72h.lg10",
                                "AUCrtLOG_48.72h.lg10")
  

   

  order.x.axis.MDK <- c("MDK_25pct.lg10",
                      "MDK_50pct.lg10",
                      "MDK_75pct.lg10",
                      "MDK_90pct.lg10")

  order.x.axis.overall.order <- c(order.x.axis.LC.fraction,
                               order.x.axis.AUC.realtime,
                               order.x.axis.AUC.realtimeLOG,
                               order.x.axis.MDK)
  


# Changing colunm nammes
  order.x.axis.overall.order.T_Ila1 <- paste(time.kill.def.part_I  ,
                                    order.x.axis.overall.order,
                                    time.kill.def.part_II,
                                    sep="_")
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Ila1BaSic <- gsub(  "Ila1",
                                                "Ila1BaSic",
                                                 order.x.axis.overall.order.T_Ila1)
  
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Ila2BaSic <- gsub(
                                                "Ila1BaSic",
                                                 "Ila2BaSic",
                                                  order.x.axis.overall.order.T_Ila1BaSic)
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Ila2 <- gsub(
                                                "Ila1",
                                                 "Ila2",
                                                    order.x.axis.overall.order.T_Ila1 )
  
  
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Thr1BaSic <- gsub(   "Ila1BaSic",
                                                "Thr1BaSic",
                                                  order.x.axis.overall.order.T_Ila1BaSic)
  
# Changing colunm nammes
  order.x.axis.overall.order.T_Thr2BaSic <- gsub(  
                                                "Thr1BaSic",
                                                  "Thr2BaSic",
                                                  order.x.axis.overall.order.T_Thr1BaSic )
  
  
    
# Changing colunm nammes
  order.x.axis.overall.order.T_Thr1 <- gsub(            "Thr1BaSic",
                                                "Thr1",
                                                    order.x.axis.overall.order.T_Thr1BaSic)
  
# Changing colunm nammes
  # Changing colunm nammes
  order.x.axis.overall.order.T_Thr2 <- gsub(            "Thr1BaSic",
                                                "Thr2",
                                                    order.x.axis.overall.order.T_Thr1BaSic)
  
  
  
  if ( old.or.new.tracking == "Trkv1") {
    
      
   # Ordering all time kill cuvre defitions
  order.x.axis.overall.order.T_piDyn <- gsub(  "Ila1",
                                                "piDyn",
                                                 order.x.axis.overall.order.T_Ila1)
  
      general.order <- c(
          order.x.axis.overall.order.T_piDyn,
        order.x.axis.overall.order.T_Ila1,
  order.x.axis.overall.order.T_Ila2,
  order.x.axis.overall.order.T_Ila1BaSic,
  order.x.axis.overall.order.T_Ila2BaSic,
  order.x.axis.overall.order.T_Thr1BaSic,
  order.x.axis.overall.order.T_Thr2BaSic,
  order.x.axis.overall.order.T_Thr1,
  order.x.axis.overall.order.T_Thr2
  )
  
    
  } else (      general.order <- c(
    # order.x.axis.overall.order.T_piDyn,
        order.x.axis.overall.order.T_Ila1,
  order.x.axis.overall.order.T_Ila2,
  order.x.axis.overall.order.T_Ila1BaSic,
  order.x.axis.overall.order.T_Ila2BaSic,
  order.x.axis.overall.order.T_Thr1BaSic,
  order.x.axis.overall.order.T_Thr2BaSic,
  order.x.axis.overall.order.T_Thr1,
  order.x.axis.overall.order.T_Thr2
  )
  )
  
     
      general.order <- paste(gsub("_",".",experiment.drug),
                             general.order,
                             experiment.id,
                             sep="_")
   
      general.order <- c(#"ExpFile",
                         #"Abx.con",
                         "Isolate",
                       #  "Killing.Def",
                         general.order)
      
      
 # general.order
  
   All.Iso.Analysis.WIDER.subset.lg10 <- pivot_wider( All.Iso.Analysis.WIDER.subset.lg10,
                                                  names_from = Killing.Features,
                                                values_from = Avg_Killing_Featureslg10)
  
  All.Iso.Analysis.WIDER.subset.Organised.KF.lg10 <-   All.Iso.Analysis.WIDER.subset.lg10[,general.order ]  
  



All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- left_join(All.Iso.Analysis.WIDER.subset.Organised.KF,
                                                                  All.Iso.Analysis.WIDER.subset.Organised.KF.lg10,
                                                                  by = "Isolate")
 #   rm(  All.Iso.Analysis.WIDER.subset.OrganisedKF.POOLED )

    
    All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- data.frame(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED )
    
    
      All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED  <- All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED%>%
      mutate_all(as.character)
    
      
      # -------------- Renaming colunm headings and respecting colunm order
        All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- melt(      All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                                                  id.vars = "Isolate")
      
                 All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <-              All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
              group_by(Isolate)%>%
              mutate(Col.order = 1,
                     Col.order = cumsum(Col.order))%>%
              ungroup()
            
          
# -- Editing column heading to reflect Drug_Analysis_AnalysisDef_MainKF_Specific.variable.KF_ExperimentID and defining the order  of colunms.
      
                   All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- separate(        All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,variable, into = c("Drug", 
                                                                                                                                               "Analysis",
                                                                                                                                              
                                                                                                                                               "Main.KF",
                                                                                                                                               "Specific.var.KF",
                                                                                                                                                "Analysis.Def",
                                                                                                                                               "Exp.ID"),
                                                                               sep = "_")
            
          All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED<-              All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED%>%
                   mutate(Specific.var.KF = gsub(".lg10", "-lg10", Specific.var.KF))%>%
                   mutate(Full.Killing.Feature = paste(Drug,
                                                       Analysis,
                                                       Analysis.Def,
                                                       Main.KF,
                                                       Specific.var.KF,
                                                       Exp.ID,sep="_"))%>%
                   select(Isolate,
                          Full.Killing.Feature,
                          value,
                          Col.order)%>%
         ungroup()

      # Defining colunm order
       Colunm.order.def <-          All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
         ungroup()%>%
         select(Full.Killing.Feature,
                Col.order)%>%
         distinct()
       
             Colunm.order.def <-Colunm.order.def$Full.Killing.Feature
               Colunm.order.def <- c("Isolate",
                                     Colunm.order.def)
       

                      All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <-        All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
                   ungroup()%>%
                   select(-Col.order)

        All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- spread(          All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                                                                                     Full.Killing.Feature,
                                                value)

          All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <-   All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
          select(Colunm.order.def)
                  
          


### Removing isolates with  low live cell fraction ----------------------------------------------------------------------------------------------------------------------------------------------
# We will finally only consider the isolates which across all defitnions had an inital live cell fraction fo at least 80%. Otherisw it is omitted. This is considered from the control experiment or considering all the intiial LCF in the ASCT experiemnt per isolate.

minLCF.df <- read.csv(path.to.minimum.LCF.data)

min.LCF.df.subset <-  All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
  filter(Isolate %in% minLCF.df$Isolate )%>%
   mutate(across(-Isolate, ~ as.character(NA)))

All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED  <- All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
  filter(!Isolate %in% minLCF.df$Isolate )

All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- rbind(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                                                                min.LCF.df.subset )


rm(
   min.LCF.df.subset)




#--- Too few numbers  for max LCF


too.few.LCF.df <- read.csv(path.to.isol.too.few.numbers.LCF.data)

too.few.LCF.df.subset <-  All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
  filter(Isolate %in% too.few.LCF.df$Isolate )%>%
   mutate(across(-Isolate, ~ as.character(NA)))

All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED  <- All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED %>%
  filter(!Isolate %in% too.few.LCF.df$Isolate )

All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED <- rbind(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                                                               too.few.LCF.df.subset )


rm(
   too.few.LCF.df.subset )




### Exporting tacking killing feature results -----------------------------------------------------------------------------------------------------------------------------------------------------
# -- End of edinting colunm headings  
   filename2 <- paste(exp.resDir,
                    "/",
                    expID,
                    "_",
                    experiment.drug,
                    "_",
                          
          old.or.new.tracking,
                    "Clinial.Isolates.Tracking_KillingFeatures.csv",
                    sep="")
                    
  # Saving ASCT Clinical isoalte data
    write.csv(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED,
                file =filename2,
              row.names=FALSE)

  
    All.Iso.Analysis.WIDER <- All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED
    
    #Houskeeping
     rm(    All.Iso.Analysis.WIDER.subset.Organised.KF,
           All.Iso.Analysis.WIDER.subset.Organised.KF.lg10,
           All.Iso.Analysis.WIDER.subset,
           All.Iso.Analysis.WIDER.subset.lg10,
              KF.labs.df)
  


### Exporting iTol compatible files ----------------------------------------------------------------------------------------------------------------------------------------------------
All.Iso.Analysis.LONG<- pivot_longer(All.Iso.Analysis.WIDER.subset.Organised.KF.POOLED, 
                        cols = -Isolate)


All.Iso.Analysis.LONG <- All.Iso.Analysis.LONG %>%
  mutate(Organised_Killing_Features = name)%>%
  select(Isolate,
         Organised_Killing_Features,
         value)


# Creating iTol file directory
path.iTol.dir <- paste(exp.resDir,
                       "/",
                                
          old.or.new.tracking,
                       "_",
          analysis.normalisation,
                       "_iTol-Tracking-Results",
                       sep="")
dir.create(path.iTol.dir)
#-- 

setwd(genDir)
iTol.Clinical.Iso.List <- read.csv("ClinicalIsolates_LookUpMasterList_BWupdated.csv")
iTol.Clinical.Iso.List <-iTol.Clinical.Iso.List %>%
  mutate( Isolate = paste("Iso.",Sample_ID_Full,
                         sep=""))%>%
  select(Isolate,
         Subspecies,
         Lane)%>%
  mutate(Subspecies = if_else(Subspecies == "", "unknown-subspecies" ,Subspecies))


iTol <- left_join(All.Iso.Analysis.LONG,
                  iTol.Clinical.Iso.List)


itol.list.features.loop <- unique(iTol$Organised_Killing_Features)


for ( i in itol.list.features.loop) {
  
  
  iTol.subset <- iTol %>%
    filter(Organised_Killing_Features == i)
  
  #---- NORMALISE values here HERE 
  
  All.Iso.Analysis.no.NAs <- iTol.subset  %>%
  ungroup()%>%
  drop_na()%>%
  ungroup()%>%
  mutate(Organised_Killing_Features = as.factor(Organised_Killing_Features))%>%
  mutate(value = if_else ( value == ">72", "72", value))%>%
  mutate(value =as.numeric(value))%>%
  group_by(
           Organised_Killing_Features)%>%
  mutate(Min.value = min(value))%>%
  mutate(Max.value= max(value))%>%
  mutate(Norm_value = (value - Min.value) / (Max.value -Min.value) )%>%
    ungroup()%>%
    select(-Min.value,
            -Max.value)%>%
    mutate(value = as.character(value))

  
    iTol.subset <- left_join(  iTol.subset,
                                All.Iso.Analysis.no.NAs)
  
  #Housekeeping
    rm(All.Iso.Analysis.no.NAs)
  
  
  iTol.subset <-  iTol.subset %>%
      filter(!is.na(Lane))%>%
    select(Lane,
           Norm_value)%>%
    mutate(Lane_Killing_feature = paste(Lane,
                                        " ",
                                       Norm_value,
                                        sep=""))%>%
    select(Lane_Killing_feature)
  
  txt <- c("DATASET_GRADIENT",
         "SEPARATOR SPACE",
     
        paste("DATASET_LABEL ", i,sep=""),
         "COLOR #ff0000",
         "STRIP_WIDTH 200",
         "USE_MID_COLOR 1",
         "COLOR_MIN #0000ff",
         "COLOR_MAX #ff0000",
         "COLOR_MID #aaaaaa",
         "DATA",
         paste("Lane " ,i, ")",sep="") )


  txt.data <- c(  iTol.subset$Lane_Killing_feature)
  
  
  iTol.txt <- c(txt,
         txt.data )
  
  
  itol.filename <- paste(
                       "iTol_Tracking",
                       "_",
                       analysis.normalisation,
                       "_",
                       i,
                       ".txt",
                       sep="")
  path.iTol.dir.results <- paste(path.iTol.dir ,
                       "/",
                       itol.filename,
                       sep="")
  
  fileConn <- file(path.iTol.dir.results)
  writeLines(iTol.txt , fileConn)
  close(fileConn)

  #Houskeeping
  rm(path.iTol.dir.results,
   itol.filename,
   iTol.txt)
  

  
}
  
#Houskeeping
rm(txt,
   iTol,
   iTol.Clinical.Iso.List,
   path.iTol.dir,
   Clinical.Iso.List,
   iTol.subset,LOC.df,
   
   neo.labels)


## Rename files based on analysis strategy at the end  --------------------------------------------------------------------------------------------------------------------------------------------
# Define the directory you want to search in
directory_path <- exp.resDir

dir.create(paste(exp.resDir,
                 "/",
                 analysis.normalisation,
                 sep=""))

# List all files in the specified directory (not in subdirectories)
files_in_directory <- c(list.files(path = directory_path, full.names = FALSE, recursive = FALSE , pattern =  "csv"),
                        list.files(path = directory_path, full.names = FALSE, recursive = FALSE , pattern =  "pdf")
                        )


old.files_in_directory  <- paste(exp.resDir,
                                 "/",
                                 files_in_directory,
                                 sep="")
                                 
new.files_in_directory  <- gsub(expID,
                                paste(expID,
                                      "_",
                                      analysis.normalisation,sep=""),
                                files_in_directory)

new.files_in_directory <-  paste(exp.resDir,
                 "/",
                 analysis.normalisation,
                 "/",
                                new.files_in_directory,
                                 sep="")

                                 

file.rename(old.files_in_directory,
new.files_in_directory )


# End of script --------------

end_time <- Sys.time()
# Calculate the time difference in minutes
execution_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

# Print the execution time
cat("Script execution time:", execution_time, "minutes\n")


message(paste("CONGRATULATIONS ASCT  Tracking analsis of", experiment.drug, experiment.id, "is complete", sep = " "))




