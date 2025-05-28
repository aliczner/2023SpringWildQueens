####    code for triangulating the radio tower and motus 

# first need to remove any duplicate values from the radio tower data
library(tidyverse)
library(lubridate)
library(sf)


## Load all the triangulation functions
source("helperFunctions.r")

tag1 <- read.csv("Spring2023TagsProcessed.csv") #3354 observations
tag2 <- tag1 %>%  distinct() #no duplicates

###   triangulation code needs spatial data

## Fix the issue with the Rare tower without ID
tower1Rare <- data.frame(x = 1286438, y = 541122.7)
tower1Rare <- st_as_sf(tower1Rare, coords=c("x","y"), crs = "ESRI:102002")
tower1Rare <- tower1Rare %>% st_transform("epsg:4326")
tag1[tag1$TowerID=="Tower1",c("Easting")] <- st_coordinates(tower1Rare)[1]
tag1[tag1$TowerID=="Tower1",c("Northing")] <- st_coordinates(tower1Rare)[2]
tag1[grep("Tower1", tag1$TowerID),"TowerID"] <- "TowerR01"

#fix a typo in an AnimalID name
tag1$AnimalID[tag1$AnimalID == "EF051913"] <- "RF051913"

## Separate rare and wallace
tag1$SiteName <- ifelse(grepl("W",substr(tag1$AnimalID, 1,1)), "Wallace","Rare")
tag1$GPSdatetime <- dmy_hms(paste(tag1$GPSdate, tag1$GPStime))
tag1 <- tag1 %>%
  filter(!is.na(Easting)) %>%
  mutate(GroupTime = binTime(GPStime, 10))

#create the uniqueID column
tag1$uniqueID <- paste0(tag1$AnimalID, tag1$GPSdatetime, tag1$TowerID)

rare.tri <- tag1 %>%
  filter(SiteName == "Rare")

wallace.tri <- tag1 %>%
  filter(SiteName == "Wallace")

# For 'Rare' sites, keep only rows where the 6th character of TowerID is 'R'
rare.tri <- rare.tri %>%
  filter(substr(TowerID, 6, 6) == "R") #284 observations

# For 'Wallace' sites, keep only rows where the 6th character of TowerID is 'W'
wallace.tri <- wallace.tri %>%
  filter(substr(TowerID, 6, 6) == "W")

# Convert to sf object
rare.triSF <- st_as_sf(rare.tri, coords = c("Easting", "Northing"), crs = 4326)
wallace.triSF <- st_as_sf(wallace.tri, coords = c("Easting", "Northing"), crs = 4326)

# triangulation code needs UTM
rareUTM <- st_transform(rare.triSF, crs = "ESRI:102002")
wallaceUTM <- st_transform(wallace.triSF, crs = "ESRI:102002")

#Extract X and Y from geometry
coordsRare <- st_coordinates(rareUTM)
coordsWallace <- st_coordinates(wallaceUTM)

# Combine with other columns
rareUTM <- cbind(st_drop_geometry(rareUTM), longitude = coordsRare[,1], 
                  latitude = coordsRare[,2])
wallaceUTM <- cbind(st_drop_geometry(wallaceUTM), longitude = coordsWallace[,1], 
                 latitude = coordsWallace[,2])

## Rare
rare.receiver <- as.receiver(rareUTM)
rareLocate <- locate(rare.receiver)

## Wallace
wallace.receiver <- as.receiver(wallaceUTM)
wallaceLocate <- locate(wallace.receiver)



for(i in 1:nrow(rareLocate)) {
  rareUTM[rareLocate$uniqueID[i] == rareUTM$uniqueID,"longitude"] <- rareLocate[i,"X"]
  rareUTM[rareLocate$uniqueID[i] == rareUTM$uniqueID,"latitude"] <- rareLocate[i,"Y"]
}

write.csv(rareUTM, "2023SpringRareTriangulated.csv")

rare.wallace <- rbind (rareUTM, wallaceUTM)

write.csv(rare.wallace, "2023SpringBothSites.csv")

