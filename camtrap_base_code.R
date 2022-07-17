#Preparation####
#install exiftool for the first time
library("exiftoolr")

#load more library
library("camtrapR")  #camtrap related function
library("tidyverse") #modifying data frame 
library("lubridate") #modifying date and time format

#run only for the first time the code below, to download latest version in your computer
#install_exiftool("D:/Ryan/camtrapR/Exif") #change path to your directory
#after download finish, rename exiftool(-K) to exiftool.exe
#directing R to run exiftool in your specified folder
exiftoolPath("D:/Ryan/camtrapR/Exif/win_exe")
Sys.which("exiftool") #if works, it will shows your directory

#set to your folder of camera trap
wd_raw <- "D:/Ryan/camtrapR/Trial/SMTB/raw_vid"

##Load effort dataset####
effort <- read.csv("data/Station.csv")

##Read folder to make species detection####
#work for most jpg format (bushnell and reconyx)
rec_table1 <- recordTable(inDir  = wd_raw,
                          IDfrom = "directory",
                          timeZone = "Asia/Jakarta",
                          video  = list(file_formats = c("jpg", "mp4"),
                                        dateTimeTag  = "QuickTime:CreateDate")
)

#work for avi (bushnell and cuddeback)
rec_table2 <- recordTable(inDir  = wd_raw,
                          IDfrom = "directory",
                          timeZone = "Asia/Jakarta",
                          video  = list(file_formats = "avi",
                                        dateTimeTag  = "FileModifyDate")
)

#combine data
detection <- rbind(rec_table1, rec_table2)

#write detection dataset to csv
detection %>%
  select(Station, Species, Date, Time, FileName) %>%
  write.csv("data/detection.csv")


###optional features####
#add copyright
#copyrightTagToAdd <- "FFI`s IP & BKSDA Kalbar"
#addCopyrightTag(inDir = wd_raw,
#                copyrightTag = copyrightTagToAdd,
 #               askFirst = FALSE)
#read all metadata in your folder
#metadatatable <- exifTagNames(inDir = wd_raw)

#add more metadata into file
#rec_table3 <- recordTable(inDir  = wd_raw,
#                          IDfrom = "directory",
#                         additionalMetadataTags = c("EXIF:Model", "EXIF:Make"),
#                        timeZone = "Asia/Jakarta",
#                       video  = list(file_formats = c("jpg", "mp4"),
#                                        dateTimeTag  = "QuickTime:CreateDate")
)


#Generate report####
#generate ops matrix
camop <- cameraOperation(CTtable = effort,
                         stationCol = "Station",
                         setupCol = "Setup_date",
                         retrievalCol = "End_date",
                         writecsv = FALSE,
                         hasProblems = FALSE,
                         dateFormat = "%d-%b-%y"
)

##Effective trap days (exclude broken cams) ####
camdays_effective <- rowSums(camop, na.rm = T) %>%
  as.data.frame() %>%
  colSums(.)

##Calculate independent event for 30m####
detection_30m <- detection %>%
  camtrapR:::assessTemporalIndependence(deltaTimeComparedTo = "lastIndepentRecord",
                                        columnOfInterest = "Species",
                                        stationCol = "Station",
                                        minDeltaTime = 30,
                                        camerasIndependent = FALSE) %>%
  drop_na(Species) %>%
  select(Station, Species, DateTimeOriginal, FileName, n_images) %>%
  mutate(IE = 1) %>%
  filter(!Species=='StartEnd') #remove people who set cameras


##Generate summary table for each species####
sum_table1 <- detection_30m %>% group_by(Species) %>% 
  summarise_at(.vars = c("n_images", "IE"),
               .funs = c("sum")) %>%
  mutate(RAI = IE/(camdays_effective/100))

sum_table2 <- detection_30m %>%
  group_by(Species) %>%
  summarise(n_locs = length(unique(Station))) %>%
  mutate(`Psi Naiive` = n_locs/5)

sum_table <- left_join(sum_table1, sum_table2, by="Species")

#Generate Species presence maps####
##Generate Species occurence tables####
Sp_occs <- detectionMaps(CTtable = effort,
                         recordTable = detection_30m,
                         Xcol = "Longitude",
                         Ycol = "Latitude",
                         stationCol = "Station",
                         speciesCol = "Species",
                         printLabels = TRUE,
                         richnessPlot = TRUE, # by setting this argument TRUE
                         speciesPlots = FALSE,
                         addLegend = TRUE
)

##Export to shapefiles####
#load library
library(sf)

# define shapefile name
shapefileName <- "SpeciesOccurences"

# set projection: WGS 84  = EPSG:4326
# see: https://spatialreference.org
shapefileProjection <- 4326

#set output
sf_out <- "D:\\Ryan\\camtrapR\\camtrap_base\\sf_output"

# run detectionMaps with shapefile creation
Sp_occs<- detectionMaps(CTtable = effort,
                        recordTable = detection_30m,
                        Xcol = "Longitude",
                        Ycol = "Latitude",
                        stationCol = "Station",
                        speciesCol = "Species",
                        richnessPlot = FALSE, # no richness plot
                        speciesPlots = FALSE, # no species plots
                        writeShapefile = TRUE, # but shaepfile creation
                        shapefileName = shapefileName,
                        shapefileDirectory = sf_out, # change this in your scripts!
                        shapefileProjection = shapefileProjection
)

detections_sf <- st_read(dsn = sf_out,
                         layer = shapefileName)
##Inspect map ####
library(mapview)
mapview(detections_sf)
mapview(detections_sf, zcol = "n_specs")

#Generate Detection for Occupancy models####
sp_list <- unique(detection_30m$Species)

list_dh <- lapply(sp_list, FUN = detectionHistory, 
                  recordTable     = detection,
                  camOp           = camop,
                  occasionLength  = 14,
                  day1            = "survey",
                  timeZone        = "Asia/Jakarta",
                  includeEffort   = FALSE, 
                  scaleEffort     = FALSE
)
names(list_dh) <- sp_list

#Generate species activity####
#all species at once

po <- "D:\\Ryan\\camtrapR\\camtrap_base\\plot_output"

activityDensity(recordTable = detection_30m,
                allSpecies  = TRUE,
                writePNG    = TRUE,
                plotDirectory = po,
                plotR       = TRUE,
                add.rug     = TRUE)

#overlap between species
activityOverlap (recordTable = detection_30m,
                 speciesA = "Neofelis diardi",
                 speciesB = "Tragulus kanchil",
                 writePNG    = TRUE,
                 plotDirectory = po,
                 plotR = TRUE,
                 pngMaxPix = 1000,
                 linecol = c("black", "blue"),
                 linewidth = c(5,3),
                 linetype = c(1, 2),
                 olapcol = "darkgrey",
                 add.rug = TRUE,
                 extend = "lightgrey",
                 ylim = c(0, 0.25),
                 main = paste("Activity overlap: ", "Neofelis diardi", "-", "Tragulus kanchil")
)
