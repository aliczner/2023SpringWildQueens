#Code for wild gyne 2021 analyses. 
#Wild gynes were caught on flowers where observed and from a nest located at Rare

#########################
### Data preparation  ###

#Creating a datafile for just wild gynes
tags<-read.csv("tagsTriangulate.csv") #tag datasheet after tower processing and triangulation
#includes pesticide bees, triangulation in GynePesticideCode script

wild<-subset(tags, Treatment == "wild")

library(dplyr)
wildnodups<-wild %>% #checking for duplications, no duplicates were found
  distinct() 

write.csv(wild, "WildGyneTriangulated.csv")

### Making landcover raster for RARE. First classification done by hand in QGIS

## Function to calculate mode
library(sf)
landcover<- st_read("Fall2021Landclass.shp")
rasterDimX = round(extent(landcover)[2]-extent(landcover)[1], 0 ) ## Rows per raster
rasterDimY = round(extent(landcover)[4]-extent(landcover)[3], 0)  ## columns per raster
empty<-raster(nrows= rasterDimX,  ## create an empty raster with the number of rows based on extent
              ncols = rasterDimY, ## create an empty raster with the number of cols based on extent
              ext = extent(landcover), 
              crs = crs(landcover))
landcover$ClassFactor <- as.factor(landcover$landclass) ## switch classes to numbers alphabetically
rastland <- rasterize(landcover, empty, "ClassFactor") ## convert polygon to raster based on class as a number
rastland2 <- focal(rastland, matrix(1,21,21), fun=modal, NAonly=T,
                   na.rm=T)  ## Take most common land cover within 10 m of NA (for NA values only)

plot(rastland2)
writeRaster(rastland2, "landcoverRaster.tif", overwrite=T)


#need to remove any points that are outside of the towers detection range
wild<-read.csv("WildGyneTriangulated.csv")
towers<-read.csv("TowerActualLocations.csv")

landcover<-raster("landcoverRaster.tif") #Rare landcover

coordinates(towers)<-~X+Y
proj4string(towers) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
towersproj<-spTransform(towers, crs(landcover)) #reprojecting CRS to utm matching landcover
towersdf<-as.data.frame(towersproj) 

#making buffer around towers
perim<-buffer(towersproj, width=220, dissolve=T)

#removing points not in the polygon, using datafrom below
pointsIN<- tagsproj[!is.na(over(tagsproj, perim)),]
#2018 points before polygon, 1998 within polygon

pointsIN.df<-data.frame(pointsIN)
detsum<-aggregate(pointsIN.df$AnimalID, by=list(pointsIN.df$AnimalID), FUN=length)
detsum

write.csv(pointsIN.df, "TriangulatTagsPerim.csv")
tagsproj<-pointsIN

nums<- pointsIN.df %>% 
  mutate(Var_X = as.numeric(Var_X)) %>% 
  mutate(Var_Y=as.numeric(Var_Y)) %>% 
  filter(!is.na(Var_X)) %>% 
  filter(!is.infinite(Var_X)) %>% 
  filter(!is.na(Var_Y)) %>% 
  filter(!is.infinite(Var_Y))

nums %>% 
  group_by(AnimalID) %>% 
  summarise(meanX = mean(Var_X), meanY=mean(Var_Y), meanAng=mean(AngleDiff))

nums %>% 
  summarise(meanX = mean(Var_X), meanY=mean(Var_Y), meanAng=mean(AngleDiff))

###################################################################################
######  Pre-specifying Breakpoints and Dealing with Zero-inflated Variables ####
##################################################################################

#knowledge of animal behaviour may allow pre-specification of breakpoints
#can deal with a large number of zeros by lumping all zeros of a given variable into a single bin

##  Prepare the data (using lost of methods from above)

#need to start with tracks dataset (copied from above)

library(bayesmove)
library(raster)
library(tidyr)
library(lubridate)

wild<-read.csv("WildGyneTriangulated.csv")
landcover<-raster("landcoverRaster.tif") #Rare landcover

coordinates(wild)<-~x+y
proj4string(wild) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
tagsproj<-spTransform(wild, crs(landcover)) #reprojecting CRS to utm matching landcover
tagsdf<-as.data.frame(tagsproj) 

## how much landcover in the study area?

studyarea <- rgdal::readOGR ("studyArea.shp")
landclip <- mask(landcover, studyarea)
landvalue<-(ncell(landclip)) - (freq(landclip, value=NA))

landvalues <- raster::extract(landcover, studyarea)
landtable <- table(landvalues)
land.df <- data.frame(landtable)

library(dplyr)

land.df2 <- land.df %>% 
  mutate (prop.land = Freq/landvalue) %>% 
  mutate (perc.land = prop.land * 100)



#date is in the wrong format
hseptags<- separate(tagsdf %>% filter(TowerID=="handheldGPS"), GPSdate, into = c("month", "day", "year"), "_")
hseptags <- hseptags %>% mutate(day = ifelse(as.numeric(day)<10, paste0("0",day), day),
                                month = ifelse(as.numeric(month)<10, paste0("0",month), month))
tseptags<- separate(tagsdf %>% filter(TowerID!="handheldGPS"),GPSdate, into = c("day", "month", "year"), "_")
septags<-rbind(hseptags, tseptags)
septags$goodDate<-paste0("20", septags$year, "-", septags$month, "-", septags$day)

septags$dtime <- ymd_hms(paste0(septags$goodDate, septags$GPStime))
septagsNoNA <- septags %>% filter(!is.na(dtime))

tagsgdate<-septagsNoNA %>% 
  dplyr::select(x, y, dtime, AnimalID, Triangulate)
names(tagsgdate)[3] <- "date"
names(tagsgdate)[4] <- "id"
tagsgdate<-tagsgdate %>% arrange(id, date)

#calculating step length, turning angle and time interval
tracks<-prep_data(dat=tagsgdate, coord.names=c("x","y"), id="id") 

head(tracks)
unique (tracks$id) #21 unique track ids 

library(ggplot2)

ggplot() +
  geom_path(data = tracks, aes(x, y, color = id, group = id), size=0.75) +
  geom_point(data = tracks %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_d("id", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))


#boxplot(log10(tracks$dt))
#tracks <- tracks %>% filter(dt < 10000 & dt > 1)

#round time steps to specified interval
hist(tracks$dt, breaks=200, xaxt="n")
axis(side=1, at=seq(50, 10000, 50)) 

tracks<-round_track_time(dat=tracks, id="id", int=250, tol=3600, time.zone="UTC", units="secs")

#will break up the data into  rest and non-rest

tracks <- tracks %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )
#now rest column where 1 = rest and 2 = non-rest

# Create list from data frame where each element is a different track
tracks.list<- df_to_list(dat = tracks, ind = "id")

# Filter observations to time interval
tracks_filt.list<- filter_time(dat.list = tracks.list, int = 250)

#### Discretize data streams (put into bins)

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)
angle.bin.lims

dist.bin.lims=quantile(tracks[tracks$dt==250,]$step,
                       c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), na.rm=T) 
dist.bin.lims

library(tidyverse)
# Assign bins to observations
tracks_disc.list<- map(tracks_filt.list,
                       discrete_move_var,
                       lims = list(dist.bin.lims, angle.bin.lims),
                       varIn = c("step", "angle"),
                       varOut = c("SL", "TA"))
# Since 0s get lumped into bin 1 for SL, need to add a 8th bin to store only 0s
tracks_disc.list2 <- tracks_disc.list %>%
  map(., ~mutate(.x, SL = SL + 1)) %>%  #to shift bins over
  map(., ~mutate(.x, SL = case_when(step == 0 ~ 1,  #assign 0s to bin 1
                                    step != 0 ~ SL)  #otherwise keep the modified SL bin
  ))

##pre-specify breakpoints

# Only retain id, discretized step length (SL), turning angle (TA), and rest columns
tracks.list2<- map(tracks_disc.list2,
                   subset,
                   select = c(id, SL, TA, rest))

## Drop bees without non-rest state (i.e., no "2"s)
beesToDrop <- sapply(tracks.list2, function(x){
  ifelse(sum(x$rest) == nrow(x) *2 || sum(x$rest) == nrow(x), FALSE, TRUE)
})
tracks.list3 <- tracks.list2[beesToDrop]

tracks_disc.list3<-tracks_disc.list2 [beesToDrop]

# Pre-specify breakpoints based on 'rest'
breaks<- purrr::map(tracks.list3, ~find_breaks(dat = ., ind = "rest"))

### run the segmentation model

set.seed(1)

# Define hyperparameter for prior distribution
alpha<- 1

# Set number of iterations for the Gibbs sampler
ngibbs<- 10000

# Set the number of bins used to discretize each data stream
nbins<- c(8,8,2)  #SL, TA, rest (in the order from left to right in tracks.list2)

progressr::handlers(progressr::handler_progress(clear = FALSE))
future::plan(future::multisession, workers = 3)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.res<- segment_behavior(data = tracks.list3, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha, breakpt = breaks)

future::plan(future::sequential)  #return to single core

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res, type = "nbrks")

# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res, type = "LML")

## Determine MAP for selecting breakpoints (maximum a posteriori)
MAP.est<- get_MAP(dat = dat.res$LML, nburn = 5000)
MAP.est

brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

# How many breakpoints estimated per ID?
apply(brkpts[,-1], 1, function(x) length(purrr::discard(x, is.na)))

plot_breakpoints(data = tracks_disc.list2, as_date = FALSE, var_names = c("step","angle","rest"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)", "Resting"), brkpts = brkpts)

tracks.seg<- assign_tseg(dat = tracks_disc.list3, brkpts = brkpts)
head(tracks.seg)

# Select only id, tseg, SL, TA, and rest columns
tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA","rest")]

# Summarize observations by track segment
nbins<- c(8,8,2)
obs<- summarize_tsegs(dat = tracks.seg2, nbins = nbins)
head(obs)

#### Running the clustering model

set.seed(1)

# Prepare for Gibbs sampler
ngibbs <- 20000  #number of MCMC iterations for Gibbs sampler
nburn <- ngibbs/2  #number of iterations for burn-in
nmaxclust <- max(nbins) - 1  #one fewer than max number of bins used for data streams
ndata.types <- length(nbins)  #number of data types

# Priors
gamma1 <- 0.1
alpha <- 0.1

# Run LDA model
res <- cluster_segments(dat=obs, gamma1=gamma1, alpha=alpha,
                        ngibbs=ngibbs, nmaxclust=nmaxclust,
                        nburn=nburn, ndata.types=ndata.types)

plot(res$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")

### Determine the number of likely behaviour states

# Extract proportions of behaviors per track segment
theta.estim <- extract_prop(res = res, ngibbs = ngibbs, nburn = nburn, nmaxclust = nmaxclust)

# Calculate mean proportions per behaviour
(theta.means<- round(colMeans(theta.estim), digits = 3))

# Calculate cumulative sum
cumsum(theta.means)
#first 2 states comprise 98% of all observations

# Convert to data frame for ggplot2
theta.estim_df<- theta.estim %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:nmaxclust, names_to = "behaviour", values_to = "prop") %>%
  modify_at("behaviour", factor)
levels(theta.estim_df$behaviour)<- 1:nmaxclust

# Plot results
ggplot(theta.estim_df, aes(behaviour, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehaviour", y="Proportion of Total Behaviour\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#most are in the first two behaviour states, a little bit in 3

#### Classify the states as behaviours

# Extract bin estimates from phi matrix
behav.res<- get_behav_hist(dat = res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Step Length","Turning Angle","Resting"))

# Plot histograms of proportion data
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c("#21908CFF","#440154FF","#FDE725FF",
                               "grey35","grey35","grey35","grey35"), guide = "none") +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")

#first behaviour: mainly not resting, some resting. Some zero step length and many at max step length
# all turning angles except directly straight, foraging or area-restricted search
#Second behaviour: Only resting. No step length, facing straight

theta.estim.long<- expand_behavior(dat = tracks.seg, theta.estim = theta.estim, obs = obs,
                                   nbehav = 2, behav.names = c("ARS", "Rest"),
                                   behav.order = c(1,2))

# Plot results
ggplot(theta.estim.long) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behaviour\n") +
  scale_fill_viridis_d("Behaviour") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, nrow = 2)


##### Assign behavioural states to tracks

# Merge results with original data
tracks.out<- assign_behavior(dat.orig = tracks,
                             dat.seg.list = df_to_list(tracks.seg, "id"),
                             theta.estim.long = theta.estim.long,
                             behav.names = c("ARS", "Rest"))

# Map dominant behavior for all IDs
ggplot() +
  geom_path(data = tracks.out, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = tracks.out, aes(x, y, fill=behav), size=1.5, pch=21, alpha=0.7) +
  geom_point(data = tracks.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks.out %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id, scales = "free", ncol = 2)

# Map of all IDs
ggplot() +
  geom_path(data = tracks.out, aes(x, y, color = behav, group = id), size=1.25) +
  geom_point(data = tracks.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks.out %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_d("Behaviour", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Proportion ARS
ggplot() +
  geom_path(data = tracks.out, aes(x, y, color = ARS, group = id), size=1.25, alpha=0.7) +
  geom_point(data = tracks.out %>%
               slice(which(row_number() == 1)),
             aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks.out %>%
               slice(which(row_number() == n())),
             aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_distiller("Proportion\nARS", palette = "Spectral", na.value = "grey50") +
  labs(x = "Easting", y = "Northing", title = "ARS") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Proportion resting
ggplot() +
  geom_path(data = tracks.out, aes(x, y, color = Rest, group = id), size=1.25, alpha=0.7) +
  geom_point(data = tracks.out %>%
               slice(which(row_number() == 1)),
             aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks.out %>%
               slice(which(row_number() == n())),
             aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_distiller("Proportion\nRest", palette = "Spectral", na.value = "grey50") +
  labs(x = "Easting", y = "Northing", title = "Resting") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

## conditions edits
## Bees to drop is really be bees to keep. So we need the inverse
beesToRevise <- names(beesToDrop)[! beesToDrop]

## Revise Rest == 1 for all observations
tracks.outTemp <- tracks.out
tracks.outTemp[tracks.outTemp$id %in% beesToRevise &   tracks.outTemp$rest == 1 & tracks.outTemp$dt == 250 &
                 tracks.outTemp$angle > -1 & tracks.outTemp$angle < 1, "behav" ] <- "Rest"
tracks.outTemp <- tracks.outTemp %>% filter(!is.na(behav))

## join processed dataset
tracks.outWithMissing <- tracks.out %>% 
  dplyr::select(names(tracks.outTemp)) %>% 
  rbind(tracks.outTemp)  
tracks.out2<-tracks.outWithMissing

#####################
### Flight Paths  ###
####################

###adding a column to the track.out2 dataset that indicates if the wild gyne was found at the nest
library (dplyr)

nestgynes <-c("G082901", "G082902", "G082903", "G082001", "G083002", "G083003", "G083101", "G083102",
              "G090202", "G090701") #gynes that were tagged at the nest
gynetracks<- tracks.out2 %>% 
  mutate(FromNest = ifelse(id %in% nestgynes, "nest", "notnest"))

library(ggplot2)
library(leaflet)
library(shiny)
library(dplyr)
library(htmltools)

## set leaflet CRS
UTMtoLatLon <- function(dfToConvert){
  dfIn <- dfToConvert
  coordinates(dfIn) <- ~x+y
  proj4string(dfIn) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m"
  dfIntLatLon <- spTransform(dfIn, CRS("+proj=longlat +datum=WGS84"))
  dfLatLon <- data.frame(dfIntLatLon)
  return(dfLatLon)
}

## CSS code to style title in leaflet map
tag.map.title <- htmltools::tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(10%,20%);
    position: fixed !important;
    left: 70%;
    bottom: 10%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 16px;
  }
"))



i="G083001"
j=2
allIds <- unique(gynetracks$id)
for( i in allIds) {
  beeSingular <- gynetracks %>% filter(id == i) %>% UTMtoLatLon(.) 
  ## need to find out number of steps before second loop
  ## maximum 65 frames for very long steps
  nSteps <- ifelse(nrow(beeSingular) > 65, 65, nrow(beeSingular))
  for( j in 2:nSteps) { ## need to start at 2 to have a line
    beeStep<-beeSingular[1:j,]
    
    colorSelect <- colorRampPalette(c("blue","orange"))
    beeStep$colour <- colorSelect(j)
    
    ## HTML code to add title - uses CSS styles from above
    title <- htmltools::tags$div(
      tag.map.title, HTML(paste0(beeStep[j, "date"]))
    )  
    
    leafletMap <- leaflet() %>% addTiles() %>% 
      addProviderTiles('Esri.WorldImagery') %>% 
      setView(-80.352, 43.377, zoom = 16) %>% 
      addControl(title, position = "topleft", className="map-title")
    for(k in 1:nrow(beeStep)){
      lineStart <- k
      lineEnd <- k+1
      leafletMap <- addPolylines(leafletMap, 
                                 data = beeStep[lineStart:lineEnd,], 
                                 lng = ~x, lat = ~y, color = ~colour)
    }
    leafletMap
    
    ## Make directory
    beeDir <- paste0("./scratch/",i)
    dir.create(beeDir, showWarnings = FALSE)
    ## save file
    mapview::mapshot(leafletMap, file = paste0(beeDir,"/", i, "_", j, ".png"))
    print(paste0(i, "-", j)) ## prints iteration
  }
}

## Need to rename all the pngs to have equal padding of numbers so
## that smaller numbers dont come before biggers ones
## e.g., 3 > 21
allFiles <- list.files("scratch", recursive = T, full.names = T)
listLengths <- strsplit(allFiles, "_")
fileID <- do.call(rbind, listLengths)[,2]
fileID_addZeros <- ifelse(nchar(fileID) == 5, paste0("00",fileID), 
                          ifelse(nchar(fileID) == 6, paste0("0",fileID), fileID))
renameFiles <- paste0(do.call(rbind, listLengths)[,1], "_",fileID_addZeros)
file.rename(allFiles, renameFiles)

## function to save as GIF
makeGIFwithPNG <- function(directory, savePath, FPS) {
  require(magick)
  ## list file names and read in
  imgs <- list.files(directory, full.names = TRUE)
  maxImages <- ifelse(length(imgs) > 50, 50, length(imgs))
  imgs <- imgs[1:maxImages] ## cap the max number of images
  img_list <- lapply(imgs, image_read)
  
  ## join the images together
  img_joined <- image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated <- image_animate(img_joined, fps = FPS)
  
  ## save to disk
  image_write(image = img_animated,
              path = savePath)
}

processedBeeMaps <- list.dirs("./scratch")[-1]
for(k in 10:length(processedBeeMaps)){
  makeGIFwithPNG(processedBeeMaps[k], paste0(processedBeeMaps[k], ".gif"), 2)
  print(k)
  print(Sys.time())
}

# Making flight paths for each bee
library(leaflet)

allIds <- unique(gynetracks$id)
for( i in allIds) {
  beeSingular <- gynetracks %>% filter(id == i)
  beePlot <- ggplot(data = beeSingular, aes(x, y, color = date, label = date)) +
    geom_path(size=0.7) +
    geom_text() +
    scale_color_distiller(palette = "Spectral")
  labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  beePlot
  ggsave( paste0("flightpaths/",i,".pdf"), width = 11, height = 9)
  print(i)
}


##make a figure of flight paths
nottracks<-subset(gynetracks, id %in% c("G081302", "G081601", "G082701", "G090702", "G091001"))
nesttracks<-subset(gynetracks, id %in% c("G082903", "G083001", "G083002", "G083101","G090202"))

library(raster)

nottracksGPS <- nottracks
coordinates(nottracksGPS) <- ~x+y
proj4string(nottracksGPS) <- crs(landcover)
nottracksGPS <- spTransform(nottracksGPS, "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +towgs84=0,0,0")
nottracksGPS <- data.frame(nottracksGPS)

library(leaflet)

colourDF1 <- data.frame(id = unique(nottracks$id), colours = c("#FFFBFF","#FB5012","#FFD23F",
                                                               "#01FDF6", "#3772FF"))
nottracksGPS  <- nottracksGPS %>% select(-colours) %>% left_join(colourDF1)

leaflet() %>% addTiles() %>% 
  addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addPolylines(data = nottracksGPS %>% filter(id == "G081302"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1.4) %>% 
  addPolylines(data = nottracksGPS %>% filter(id == "G081601"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1.25) %>% 
  addPolylines(data = nottracksGPS %>% filter(id == "G082701"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1) %>% 
  addPolylines(data = nottracksGPS %>% filter(id == "G090702"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1) %>% 
  addPolylines(data = nottracksGPS %>% filter(id == "G091001"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1) %>% 
  addLegend("topright", color=unique(nottracksGPS$colours), labels=unique(nottracksGPS$id), 
            opacity = 1)

nesttracksGPS <- nesttracks
coordinates(nesttracksGPS) <- ~x+y
proj4string(nesttracksGPS) <- crs(landcover)
nesttracksGPS <- spTransform(nesttracksGPS, "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +towgs84=0,0,0")
nesttracksGPS <- data.frame(nesttracksGPS)

library(leaflet)

colourDF2 <- data.frame(id = unique(nesttracks$id), colours = c("#FFFBFF","#FB5012","#FFD23F",
                                                                "#01FDF6", "#3772FF"))
nesttracksGPS <- left_join(nesttracksGPS, colourDF2)

leaflet() %>% addTiles() %>% 
  addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addPolylines(data = nesttracksGPS %>% filter(id == "G082903"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1) %>% 
  addPolylines(data = nesttracksGPS %>% filter(id == "G083001"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1) %>% 
  addPolylines(data = nesttracksGPS %>% filter(id == "G083002"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1) %>% 
  addPolylines(data = nesttracksGPS %>% filter(id == "G083101"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1) %>% 
  addPolylines(data = nesttracksGPS %>% filter(id == "G090202"),
               ~x, ~y, label = ~id, group = ~id, color = ~colours, weight = 1) %>% 
  addLegend("topright", color=unique(nesttracksGPS$colours), labels=unique(nesttracksGPS$id), 
            opacity = 1)

##to do stats with landcover need to add landcover data for true and random points
library(amt)
gynetracks.sp<-gynetracks
coordinates(gynetracks.sp)<-~x+y
proj4string(gynetracks.sp) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs" 

tracksRandom <- spsample(gynetracks.sp, 3161, "random")

rand_sl = random_numbers(make_exp_distr(), n = 1e+05)
rand_ta = random_numbers(make_unif_distr(), n = 1e+05)

getRandomCoords <- function(df, N) { ## needs dataframe and number of random steps per observation
  rand_sl =  random_numbers(fit_distr(df$step, "gamma"), N) ## set random steps based on dataframe
  rand_ta = random_numbers(fit_distr(df$angle, "vonmises"), N) # set random angles based on dataframe
  randomizedCoords <- data.frame() ## empty dataframe to fill with randomized points
  for(i in 1:nrow(df) ){
    x2 <- df[i,"x"] + (rand_sl * cos(rand_ta)) ## calculate change in x-coordinate based on origin, step, and angle
    y2 <- df[i,"y"] + (rand_sl * sin(rand_ta)) ## calculate change in y-coordinate based on origin, step, and angle
    tempDf <- data.frame(id = df[i, "id"], case = "random", ## create dataframe with original bee ID, change in points, and movement info
                         x1= df[i,"x"], y1= df[i,"y"], x2, y2, 
                         sl = rand_sl, ta = rand_ta)
    randomizedCoords <- rbind(randomizedCoords, tempDf) ## output all randomizations per observation
  }
  return(randomizedCoords)
  
}
randomCoords <- getRandomCoords(gynetracks, 10) ## calculate all random points
plot(randomCoords$x2, randomCoords$y2)

getTrueCoords <- function(df){ ## revise the true steps to match the same as random
  revisedCoordsDF <- df %>% 
    group_by(id) %>%  ## select by bee identifier
    mutate(x2 = lead(x), y2 = lead(y)) %>%  ## create a new column for the end travel point
    filter(!is.na(x2)) %>%  ## drop all initial x-y that don't have an end point
    mutate(case = "true") %>%  ## add a column to identify true
    select(case, x1 = x, y1 = y, x2, y2, sl = step, ta = angle) ## match same data structure
  return(revisedCoordsDF)
}
trueCoords <- getTrueCoords(gynetracks) ## revise data to get all true points
gynetracks.tf <- rbind(randomCoords,trueCoords ) 

#landcover raster needs to be a rasterstack
landStack <- stack()
for(i in 1:7) { 
  tempRaster <- landcover
  tempRaster[tempRaster != i] <- 0
  tempRaster[tempRaster == i] <- 1
  names(tempRaster) <- paste0("landcover",i)
  landStack <- stack(landStack, tempRaster)
}

#extract landcover data
gynetracks.tfs <- gynetracks.tf
coordinates(gynetracks.tfs)<-~x2+y2
proj4string(gynetracks.tfs) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs" 
library(raster)
landex<- raster::extract(landStack, gynetracks.tfs)

#add covariates from landStack to the dataframe for SSF
gynetracks.tfs$agriculture <- landex[,1]
gynetracks.tfs$developed <- landex[,2]
gynetracks.tfs$forest <- landex[,3]
gynetracks.tfs$highFloral <- landex[,4]
gynetracks.tfs$lowFloral <- landex[,5]
gynetracks.tfs$modFloral <- landex[,6]
gynetracks.tfs$wetland<- landex[,7]
head(gynetracks.tfs)

gynetracks.tfs$presence <- ifelse(gynetracks.tfs$case == "random", FALSE, TRUE)
  
### Step-selection functions

library(survival)
SSF1<-clogit(presence ~ highFloral + lowFloral + modFloral,
             method= "approximate", data=gynetracks.tfs) 
SSF1
summary(SSF1)
library(sjPlot)
plot_model(SSF1) #selects for floral areas

SSF2<-clogit(presence ~ 1, method="approximate", data=gynetracks.tfs) #not the best model

gynetracks.tfs.NONA <- data.frame(gynetracks.tfs)
gynetracks.tfs.NONA[is.na(gynetracks.tfs.NONA)] <- 0
SSF6<-clogit(presence~crop+developed+forest+highFloral+lowFloral+modFloral+water+wetland,
             method="approximate", na.action="na.fail", data=gynetracks.tfs.NONA) #significant forest, avoids
plot_model(SSF6)
library(MuMIn)

D1<-dredge(SSF6)

# explain what happened above: could not apriori determine which variables to include
#realized I could use dredge on the full model for variable selection
#there are multiple best models:
#high low and mod floral
#forest, mod floral and wetland
#forest, mod floral
#forest high and mod floral
#high mod and low floral and wetland
#forest, high mod and low floral
#forest, high mod and wetland
#forest high

SSF3 <- clogit(presence ~ forest +modFloral +wetland, method="approximate", data=gynetracks.tfs)
plot_model(SSF3)
SSF4 <- clogit(presence ~ forest + modFloral, 
               method ="approximate", data=gynetracks.tfs)
plot_model(SSF4) #selects for both
SSF5 <- clogit(presence ~ forest+highFloral+modFloral,
               method="approximate", data=gynetracks.tfs)
plot_model(SSF5)
SSF7 <- clogit(presence ~ forest + highFloral + modFloral +lowFloral + wetland, 
                method="approximate", data=gynetracks.tfs)
plot_model(SSF7)
SSF8 <- clogit(presence ~ forest + highFloral + modFloral +lowFloral,
               method="approximate", data=gynetracks.tfs) #marginally selecting for everything
plot_model(SSF8)

nest<-data.frame(y=540940.972081,	x=1286341.239233)
coordinates(nest)<-~x+y
proj4string(nest) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"

nest.ras<-distanceFromPoints(landcover, nest)
dnest.ras<-raster::extract(nest.ras, gynetracks.tfs)
gynetracks.tfs$dnest <- dnest.ras

##going to test if distance from nest predicts anything for the nest gynes

#lost the nest/notnest column need to add it back
nestgynes <-c("G082901", "G082902", "G082903", "G083001", "G083002", "G083003", "G083101", "G083102",
              "G090202", "G090701") #gynes that were tagged at the nest
gynetracks.tfs2<- as.data.frame(gynetracks.tfs) %>% 
  mutate(FromNest = ifelse(id %in% nestgynes, "nest", "notnest"))

gynetracks.nest2 <- data.frame(gynetracks.nest)
gynetracks.nest2[is.na(gynetracks.nest2)] <- 0
SSF12<-clogit(presence~crop+developed+forest+highFloral+lowFloral+modFloral+wetland, 
              method="approximate", data=gynetracks.nest2, na.action="na.fail")
d2<-dredge(SSF12) #the null model has the lowest AIC

mean(gynetracks.nest2$sl)
gynetracks.nest2 %>% group_by(case) %>% summarise(mean(sl))
gynetracks.nest2 %>% group_by(case) %>% summarise(mean(dnest))

###################################
### Step Lengths and turn angles ###
###################################

library(pastecs)
library (summarytools)#for descriptive statistics
#summary stats for step lengths
stat.desc(gynetracks)
gynetracks.abs <- abs(gynetracks$angle)
stat.desc(gynetracks.abs)

## Step length x time of day
library(dplyr)
library(tidyr)
library(emmeans)

#Separating dataset to make a daytime column
gynetracks.period <- separate(gynetracks, date, into = c("day","time"), sep=" ")
gynetracks.period <-separate(gynetracks.period, time, sep=":", into=c("hour","minute","sec"))

gynetracks.period[,"daytime"] <- ifelse( as.numeric(gynetracks.period$hour) < 5
                                         | as.numeric(gynetracks.period$hour) > 20, "night", "day")

t1<-lm(step ~ daytime, data=gynetracks.period)
summary(t1)
anova(t1) #no significant difference in step lengths between time of day

t2 <- lm(angle~daytime, data=gynetracks.period)
summary(t2)
anova(t2) #no significant different in turning angle between time of day

t3<- lm(abs(angle)~daytime, data=gynetracks.period)
anova(t3) #not significant

stby(data= gynetracks.period, INDICES=gynetracks.period$daytime, FUN=descr, stats="common")

## adding landcover as categorical variable
gynetracks.ps <- gynetracks.period
coordinates(gynetracks.ps)<-~x+y
proj4string(gynetracks.ps) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"

landex<- raster::extract(landcover, gynetracks.ps) #makes a vector of landcover numbers
gynetracks.ps$landtype <- landex
gynetracks.psf<-as.data.frame(gynetracks.ps)
landClasses <- data.frame(landtype = 1:7, landTypeNew = c("agriculture","developed","forest", "highFloral", "modFloral", "lowFloral", "wetland"))
gynetracks.psf <- gynetracks.psf %>% left_join(landClasses) %>% select(-landtype) %>% rename(landtype = landTypeNew)

gynetracks.psf<- gynetracks.psf %>% 
  filter(!is.na(behav)) %>% 
  mutate(ARS = ifelse(behav == "ARS", 1, 0),
         Rest = ifelse(behav == "Rest", 1, 0))
gynetracks.psf %>% select(-Forage)#accidentally added forage from previous code

s1<-lm(step~landtype, data=gynetracks.psf) #significant
car::Anova(s1, type=2)
s1.a<-emmeans(s1, pairwise~landtype, data=gynetracks.psf) 
s2 <-lm(step~landtype*daytime, data=gynetracks.psf) #landtype x daytime significant
library(emmeans)
s2.a<-emmeans(s2, pairwise~landtype*daytime, data=gynetracks.psf)
s2.apostDFSignificant <- data.frame(s2.a$contrasts) %>% filter(p.value < 0.05)
write.csv(s2.apostDFSignificant, "s2aPost.csv", row.names= F)

s3<-lm(step~behav*daytime, data=gynetracks.psf) #not significant interaction
s3.a<-emmeans(s3, pairwise~behav*daytime, data=gynetracks.psf)

s4<-lm(abs(angle) ~ landtype, data=gynetracks.psf)
anova(s4)
s4.a <-emmeans(s4, pairwise ~landtype, data=gynetracks.psf)

s5 <- lm(angle ~landtype * daytime, data=gynetracks.psf) #significant interaction
s5.a <-emmeans(s5, pairwise ~ landtype *daytime, data=gynetracks.psf)

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))
stepplotdat <-gynetracks.psf %>% group_by (landtype) %>% 
  summarize (avgStep = mean(step, na.rm=T), errorstep = se(step))

stepplotdat2 <-stepplotdat %>% filter(!is.na(landtype))

library(ggplot2)
ggplot(data=stepplotdat2, aes(x=landtype, y=avgStep)) +
  theme_classic() +
  geom_bar(stat="identity", fill="#858786") +
    geom_errorbar(aes(x = landtype, ymin = avgStep - errorstep+1, ymax = avgStep + errorstep+1), 
                  width = 0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 14)) +
  theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=12)) +
  scale_y_continuous(name= "average step length") +
  scale_x_discrete(limits=c("agriculture", "forest", "lowFloral", "modFloral", "highFloral"),
                   labels=c("agriculture", "forest", "low floral", "moderate floral", "high floral"))

turnplotdat <-gynetracks.psf %>% group_by (landtype) %>% 
  summarize (avgTurn = mean(abs(angle), na.rm=T), errorturn = se(angle))

turnplotdat2 <-turnplotdat %>% filter(!is.na(landtype))

ggplot(data=turnplotdat2, aes(x=landtype, y=avgTurn)) +
  theme_classic() +
  geom_bar(stat="identity", fill="#858786") +
  geom_errorbar(aes(x = landtype, ymin = avgTurn - errorturn, ymax = avgTurn + errorturn), 
                width = 0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 14)) +
  theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=12)) +
  scale_y_continuous(name= "average turn angle (radians)") +
  scale_x_discrete(limits=c("agriculture", "forest", "lowFloral", "modFloral", "highFloral"),
                   labels=c("agriculture", "forest", "low floral", "moderate floral", "high floral"))


### behaviour as a response variable

b1<-glm(ARS ~ landtype, family="binomial", data=gynetracks.psf) 
anova(b1, test="Chisq") #significant
b1.a<-emmeans(b1, pairwise~landtype, data=gynetracks.psf)

b2<-glm(ARS~landtype*daytime, family="binomial", data=gynetracks.psf)
anova(b2, test="Chisq")
b2.a<-emmeans(b2, pairwise~landtype*daytime, data=gynetracks.psf)

b3<- glm(ARS~daytime, family="binomial", data=gynetracks.psf)
anova(b3, test="Chisq") 
b3.a<-emmeans(b3, pairwise~daytime, data=gynetracks.psf)

b4<-glm(Rest~landtype, family="binomial", data=gynetracks.psf)
anova(b4, test="Chisq") #significant
b4.a <-emmeans(b4, pairwise~landtype, data=gynetracks.psf)

b5<-glm(Rest~daytime, family="binomial", data=gynetracks.psf)
anova(b5, test="Chisq")
b5.a<-emmeans(b5, pairwise~daytime, data=gynetracks.psf)

b6<-glm(Rest~daytime*landtype, family="binomial", data=gynetracks.psf)
anova(b6, test="Chisq")
b6.a<-emmeans(b6, pairwise~landtype*daytime, data=gynetracks.psf)

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

aggtab <- gynetracks.psf[,c("landtype", "daytime" ,"ARS", "Rest")]
behavPlot<-aggtab %>% group_by (landtype, daytime) %>% 
  summarise(ARS = mean(ARS),  rest = mean(Rest)) %>%
  tidyr::gather(behav, value, ARS:rest )

ggplot(behavPlot, aes(x= landtype, y = value, fill=behav)) + theme_classic() + 
  facet_grid(daytime ~ behav)+ geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#858786", "#858786")) +
  scale_x_discrete(name="", limits=c("agriculture", "forest", "highFloral", "modFloral", "lowFloral"),
                   labels=c("agriculture", "forest", "high floral", "moderate floral", "low floral")) +
  scale_y_continuous(name="proportion of behaviour") +
  theme(strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y  = element_text(size=14)) +
  guides(fill=FALSE)

####################
### Total Flight ###
####################

gynetracks.psf2<-gynetracks.psf %>% 
  group_by(id, daytime, landtype, FromNest, behav) %>% 
  summarise(totalFlight=sum(step, na.rm=T)) %>% 
  data.frame()

sumflight<-gynetracks.psf %>% 
  group_by(id) %>% 
  summarise(totalFlight=sum(step, na.rm=T)) %>% 
  data.frame()

library(pastecs)
stat.desc(sumflight)

#is there a difference in total flight by time of day
f1<-lm(totalFlight~daytime, data=gynetracks.psf2)#not sig
 
#is there a difference in total flight by landtype
f2<- lm(totalFlight~landtype, data=gynetracks.psf2) #not sig

#is there a difference in total flight by time of day and landtype
f3<- lm(totalFlight~landtype*daytime, data=gynetracks.psf2) #not sig

##################################
### non-track based stats     ###
##################################

#tagsgdate dataframe (original "occurrence points") 2018 obs
#need to add daytime, nest status and landcover to the dataframe

#adding nest status
nestgynes <-c("G082901", "G082902", "G082903", "G083001", "G083002", "G083003", "G083101", "G083102",
              "G090202", "G090701") #gynes that were tagged at the nest

gyne.pts<- as.data.frame(tagsgdate) %>% 
  mutate(FromNest = ifelse(id %in% nestgynes, "nest", "notnest"))

#adding time of day
gyne.pts <- separate(gyne.pts, date, into = c("day","time"), sep=" ")
gyne.pts <-separate(gyne.pts, time, sep=":", into=c("hour","minute","sec"))
gyne.pts[,"daytime"] <- ifelse( as.numeric(gyne.pts$hour) < 5
                                         | as.numeric(gyne.pts$hour) > 20, "night", "day")

#adding landcover
gyne.pts.sp<-gyne.pts
  coordinates(gyne.pts.sp)<-~x+y
proj4string(gyne.pts.sp) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"

landex.pts<- raster::extract(landcover, gyne.pts.sp) #makes a vector of landcover numbers
gyne.pts$landtype <- landex.pts
landClasses <- data.frame(landtype = 1:7, landTypeNew = c("agriculture","developed","forest", "highFloral", "modFloral", "lowFloral", "wetland"))
gyne.pts <- gyne.pts %>% left_join(landClasses) %>% select(-landtype) %>% rename(landtype = landTypeNew)

###calculating the number of days flying
daysFlown <- gyne.pts %>%
  group_by(id) %>% 
  summarise(firstDay = as.Date(min(day)), lastDay = as.Date(max(day))) %>% 
  mutate(daysFlown = lastDay - firstDay)

daysFlown.df<-data.frame(daysFlown)

#extracting the last detection point

lastDay<-gyne.pts %>% 
  group_by(id) %>% 
  slice(which.max(as.Date(day)))

firstDay<-gyne.pts %>% 
  group_by(id) %>% 
  slice(which.min(as.Date(day)))

library(ggplot2)

#plot of landcover types associated with last occurrence point
ggplot(data=lastDay, aes(x=landtype)) +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 14)) +
  theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=12)) +
  geom_bar(fill="#858786", stat="count") +
  scale_y_continuous(breaks=seq(0, 15, 3), name= "count of last occurrences") +
  scale_x_discrete(limits=c("forest", "highFloral", "modFloral", "lowFloral"),
                   labels=c("forest", "high floral", "low floral", "moderate floral"))

#making a map of where the last points are
lastDay.sp<-lastDay
coordinates(lastDay.sp)<-~x+y
proj4string(lastDay.sp) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"
lastDay.sp <- spTransform(lastDay.sp, "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +towgs84=0,0,0")

set.seed(11) ## makes the plot repeatable
randomAdjustment <- function(nObs, adjRange){
  return(runif(nObs, -0.5, 0.5)/adjRange)
}
## Then run this
lastDay.sp$xAdj <- lastDay.sp$x + randomAdjustment(21, 1000)
lastDay.sp$yAdj <- lastDay.sp$y + randomAdjustment(21, 1000)

library(leaflet)

pal<-colorFactor(c("#660000", "#C40808", "#FF0000", "#DA6A08", "#FF7F00", "#FDAD5E", "#DCB91E",
                   "#FFFF00", "#FFFF99", "#CC005F", "#FF3355", "#F17B8B", "#025033", "#2EB14F",
                   "#90EE90", "#090984", "#1500FF", "#2A98CB", "#520A85", "#7E1FC2", "#C588EB"),
                 lastDay.sp$id)

leaflet(lastDay.sp) %>% addTiles() %>%
  addProviderTiles('Esri.WorldImagery') %>%
  setView(-80.352, 43.377, zoom = 16) %>%
  addCircleMarkers(~xAdj, ~yAdj, stroke=FALSE, color= ~pal(id), fillOpacity=0.9) %>%
  addLegend("topright", pal=pal, values= ~id, opacity=1)

#distance from nest on last day
lastDayNest <-subset(lastDay, FromNest=="nest")
coordinates(lastDayNest)<-~x+y
proj4string(lastDayNest) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"
last.pts<- raster::extract(nest.ras, lastDayNest) #makes a vector of landcover numbers
lastDayNest$dNest <- last.pts
as.data.frame(lastDayNest)

library(pastecs)
stat.desc(lastDayNest)

### Distance from nest for gynes from the nest

gynetracks.nest.true<-subset(gynetracks.tfs2, case=="true")
gynetracks.nest.true<-subset(gynetracks.nest, FromNest=="nest")

dnest.gyne <- gynetracks.nest.true %>%
  group_by(id) %>% 
  summarise(maxDis = (max(dnest)), meanDis=(mean(dnest)), SDDis=(sd(dnest)))
dnest.gyne

##################################
### Minimum Convex Polygon  ###
##################################

library(adehabitatHR)
library(dplyr)
library(tidyr)

gyne.pts5 <- gyne.pts %>% group_by(id) %>% mutate(nHits = length(id)) %>% 
  filter(nHits > 5)
gyne.pts5<-gyne.pts5[, c("id", "x", "y")] #mcp needs only 3 columns
coordinates(gyne.pts5) <- ~x+y #specifying coordinates
proj4string(gyne.pts5) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs"

gyne.mcp <- mcp(gyne.pts5, percent = 100, unout="km2")
gyne.mcp
names(gyne.mcp)[1]<-"GyneID"
gyne.mcp.df<-as.data.frame(gyne.mcp)

tagsmcpDF<-septagsNoNA
tagsmcpDFoutlierRemoved <-tags.mcp[tags.mcp$area < max(tags.mcp$area),]
gyne.mcp.t <- spTransform(gyne.mcp, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

pal2<-colorFactor(c("#660000", "#C40808", "#FF0000", "#DA6A08", "#FF7F00", "#FDAD5E", "#DCB91E",
                    "#FFFF00", "#FFFF99", "#CC005F", "#FF3355", "#F17B8B", "#025033", "#2EB14F",
                    "#90EE90", "#090984", "#1500FF", "#2A98CB", "#520A85", "#7E1FC2", "#C588EB"),
                  gyne.mcp.t$GyneID)

library(leaflet)
## CSS code to style title in leaflet map
tag.map.title <- htmltools::tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(10%,20%);
    position: fixed !important;
    left: 70%;
    bottom: 10%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 16px;
  }
"))
allGynes <- gyne.mcp.t$GyneID

for(i in allGynes){
## HTML code to add title - uses CSS styles from above
title <- htmltools::tags$div(
  tag.map.title, HTML(paste0(i))
)  
gyneSingular <- subset(gyne.mcp.t, GyneID ==i)

leafletMap <- leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addPolygons(data = gyneSingular, weight = 5,
              col="orange",
             fillOpacity=0.3)    %>% 
      addControl(title, position = "topleft", className="map-title")

mapview::mapshot(leafletMap, file = paste0("MCP_",i,".png"))
}

gyne.mcp.tReproj <- spTransform(gyne.mcp.t, crs(landcover))
outDF <- data.frame()
for(i in allGynes){
  gyneSingular <- subset(gyne.mcp.tReproj, GyneID ==i)
  gyneland <- mask(landcover, gyneSingular)
  gynetable <- data.frame(freq(gyneland))
  gynetableNoNA <- gynetable[!is.na(gynetable$value),]
  gyneout <- gynetableNoNA
  gyneout[,"Proportion"] <- gyneout$count/sum(gyneout$count)
  gyneout[,"GyneID"] <- i
  outDF <- rbind(outDF, gyneout)
}

library(tidyr)
MCPdata <- read.csv ("MCPLancoverEach.csv")
MCPdata2 <- spread(MCPdata, value, Proportion)
MCPdata2[is.na(MCPdata2)] <-0

p1 <- glm (Area ~ agriculture +forest + highFloral + modFloral + lowFloral , family = "binomial", 
           data = MCPdata2, na.action="na.fail",)
p2 <- dredge(p1)

library(ggplot2)

#plot of landcover types associated with last occurrence point

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

mcptab <- MCPdata[,c("value", "Proportion")]
MCPPlot<-mcptab %>% group_by (value) %>% 
  summarise(Proportion = mean(Proportion), errorprop = se(Proportion)) 

ggplot(data=MCPPlot, aes(x=value, y = Proportion)) +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size = 14)) +
  theme(axis.title.y=element_text(size=14), axis.text.y=element_text(size=12)) +
  geom_bar(fill="#858786", stat="identity") +
  geom_errorbar(aes(x = value, ymin = Proportion - errorprop, ymax = Proportion + errorprop), 
                width = 0) +
  scale_y_continuous(name= "proportion within MCP") +
  scale_x_discrete(limits=c("agriculture", "forest", "highFloral", "modFloral", "lowFloral"),
                   labels=c("agriculture", "forest", "high floral", "moderate floral", "low floral"))

turnplotdat <-gynetracks.psf %>% group_by (landtype) %>% 
  summarize (avgTurn = mean(abs(angle), na.rm=T), errorturn = se(angle))

turnplotdat2 <-turnplotdat %>% filter(!is.na(landtype))

ggplot(data=turnplotdat2, aes(x=landtype, y=avgTurn)) +
  theme_classic() +
  geom_bar(stat="identity", fill="#858786") +
  geom_errorbar(aes(x = landtype, ymin = avgTurn - errorturn, ymax = avgTurn + errorturn), 
                width = 0) +

###########################
### Kernel Density  #######
#############################

library(tidyr)
library(MASS)
library(raster)
library (ggplot2)
library(dplyr)

landcover2 <- aggregate(landcover, fact=5, fun = modal) #reduce the resolution to increase speed

kd<-kde2d(x=gyne.pts$x, y=gyne.pts$y)
kdr<-raster(kd)

#need to change crs

crs(kdr) <-crs(landcover)
kdr.p<-projectRaster(kdr, crs=crs(gyne.mcp.t))


library(leaflet)
library(colorspace)

#set up colour palette
cPal<-sequential_hcl(6, palette = "OrYel", rev=TRUE)
cval = c(0,seq(0, maxValue(kdr.p), by = maxValue(kdr.p)/5))
cval = cval *100000 #numbers were very small, added this to clean up the legend
cval <- round(cval, 2)
cpal = c("#FFFFFF00", cPal)

leaflet() %>% addProviderTiles('Esri.WorldImagery') %>% 
  setView(-80.352, 43.377, zoom = 16) %>% 
  addRasterImage(kdr.p, colors = cpal, opacity = 0.85) %>% 
  addLegend(colors=cpal[-c(1,2)], labels = cval[-c(1,2)])

#extracting landcover variables
landcoverPoly <- rasterToPoints(landcover2)
landpoint<-as.data.frame(landcoverPoly)
coordinates(landpoint)<- ~x + y
proj4string(landpoint)<-crs(landcover2)

gynevalues <- raster::extract(kdr.p, landpoint)
gyneex<-cbind(gynevalues, data.frame(landpoint))

gyneex$landcoverRaster <-factor(gyneex$landcoverRaster)

newgyne <- gyneex %>% filter(!(landcoverRaster %in% c(7,8)))

se <- function(x) sd(x, na.rm =T) / sqrt(length(x[!is.na(x)]))

summarizedLandcover <- newgyne %>% 
  group_by(landcoverRaster) %>% 
  summarize(meanGyne = mean(gynevalues *10000, na.rm = T), errorval = se(gynevalues *10000))

ggplot(summarizedLandcover, aes(x = landcoverRaster, y = meanGyne)) + 
  geom_bar(stat = "identity", fill ="#858786") +
  geom_errorbar(aes(x = landcoverRaster, ymin = meanGyne - errorval, ymax = meanGyne + errorval), 
                 width = 0) +
  theme_classic() +
scale_x_discrete(limits=c("1", "2", "3", "4", "6", "5", "9"),
                 labels=c("agriculture", "developed", "forest", "high floral",  "moderate floral", "low floral",
                          "wetland"), name="") +
theme(axis.text.x = element_text (size = 12),
      axis.title.y = element_text(size = 14)) +
  scale_y_continuous(name="mean kernel density") 

kd1<- lm(gynevalues ~ as.factor(landcoverRaster), data=newgyne)
car::Anova(kd1, type=2)
emmeans(kd1, pairwise~as.factor(landcoverRaster), data=newgyne)

## getting random points to compare to check for bias in landcover abundance

newgyne.sp <- newgyne
coordinates(newgyne.sp)<-~x+y
proj4string(newgyne.sp) <- "+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" 

randKern <- spsample(newgyne.sp, 77758, "random")

rand.kern<-kde2d(x=randKern$x, y=randKern$y)
rand.kern.ras<-raster(rand.kern)

randvalues <- raster::extract(rand.kern.ras, landpoint)
randbind<-cbind(randvalues, data.frame(landpoint))

randbind$landcoverRaster <-factor(randbind$landcoverRaster)

kd2<- lm(randvalues ~ as.factor(landcoverRaster), data=randbind)
car::Anova(kd2, type=2)
emmeans(kd2, pairwise~as.factor(landcoverRaster), data=randbind)

anova(kd1, kd2) #kernel density model with bee data is significantly different from random


