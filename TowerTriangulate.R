####    code for triangulating the radio tower and motus 

# first need to remove any duplicate values from the radio tower data
library(tidyverse)

tag1 <- read.csv("Spring2023TagsProcessed.csv") #3354 observations
tag2 <- tag1 %>%  distinct() #no duplicates

library(dplyr)

binTime <- function(timeColumn, interval = 20) {
  # Convert string to POSIXct (today's date + time)
  time <- as.POSIXct(timeColumn, format = "%H:%M:%S", tz = "UTC")
  
  # Get seconds since midnight
  secs <- as.numeric(difftime(time, floor_date(time, "day"), units = "secs"))
  
  # Round seconds to nearest interval
  intervalSecs <- interval * 60
  roundedSecs <- round(secs / intervalSecs) * intervalSecs
  
  # Add rounded seconds back to midnight
  rounded_time <- as.POSIXct(floor_date(time, "day") + roundedSecs, 
                             origin = "1970-01-01", tz = "UTC")
  
  # Return time in HH:MM format
  format(rounded_time, "%H:%M:%S")
}

tag1$SiteName <- ifelse(grepl("W",substr(tag1$TowerID, 6,6)), "Wallace","Rare")
tag1$GroupTime <- binTime(tag1$GPStime, 10) %>% filter(!is.na(Easting))
tag1$uniqueID <- paste0(tag1$AnimalID,tag1$GPSdate,tag1$GroupTime)  
rare.tri <- tag1 %>% filter(SiteName == "Rare") #447
wallace.tri <-tag1 %>% filter(SiteName == "Wallace") #2907

###   triangulation code needs spatial data
library(sf)


tower1Rare <- data.frame(x = 1286438, y = 541122.7)
tower1Rare <- st_as_sf(tower1Rare, coords=c("x","y"), crs = "ESRI:102002")
tower1Rare <- tower1Rare %>% st_transform("epsg:4326")
rare.tri[rare.tri$TowerID=="Tower1",c("Easting")] <- st_coordinates(tower1Rare)[1]
rare.tri[rare.tri$TowerID=="Tower1",c("Northing")] <- st_coordinates(tower1Rare)[2]
rare.tri[grep("Tower1", rare.tri$TowerID),"TowerID"] <- "TowerR01"

# Convert to sf object
rare.triSF <- st_as_sf(rare.tri, coords = c("Easting", "Northing"), crs = 4326)

# triangulation code needs UTM
rareUTM <- st_transform(rare.triSF, crs = "ESRI:102002")

#Extract X and Y from geometry
coords <- st_coordinates(rareUTM)

# Combine with other columns
rareUTM <- cbind(st_drop_geometry(rareUTM), longitude = coords[,1], 
                  latitude = coords[,2])
x<-rareUTM

####    Triangulation Function    ####

as.receiver <-function(x) {
  class(x) <-"receiver"
  return(x)
}

findintersects <- function(x) {
  
  ## Convert compass Angle to standard radian measure (Lenth 1981)
  x$theta = (pi/180*(90-x$Angle)) 
  
  ## Calculate additional variables
  x$sin = sin(x$theta)
  x$cos = cos(x$theta)
  x$tan = tan(x$theta)
  
  ## Create data frame to store coordinates of all intersections
  results <-data.frame(uniqueID=NA,X=NA,Y=NA) 
  
  ## Keep track of total number of intersections in the data
  counter=1; 
  for (group in unique(x$uniqueID)) {
    
    ## Create subdata to store data of the single UniqueID
    subdata=data.frame(X=x$longitude[x$uniqueID==group],Y=x$latitude[x$uniqueID==group],
                       TAN=x$tan[x$uniqueID==group])
    
    ## Calculate number of Angles in the single UniqueID
    subdata_length=nrow(subdata)
    
    ## Cycle through all combinations of Angles in the single UniqueID 
    for (i in 1:subdata_length) {
      for (j in 1:subdata_length) {
        
        ## Do not calculate intersections of Angles with themselves
        if (i !=j ) {
          
          ## Import UniqueID variable
          results[counter,1]=group
          
          ## Calculate Angle intersect
          results[counter,3]=((subdata[j,'X']-subdata[i,'X'])*subdata[i,'TAN']*subdata[j,'TAN']-subdata[j,'Y']*subdata[i,'TAN']+subdata[i,'Y']*subdata[j,'TAN'])/(subdata[j,'TAN']-subdata[i,'TAN'])
          results[counter,2]=(results[counter,3]-subdata[j,'Y'])/subdata[j,'TAN']+subdata[j,'X']
          
          ## Increment total number of calculated intersects
          counter=counter+1
        }
      }
    }
  }
  
  class(results)<-"intersect"
  return(results)
}


locateDF <- data.frame()

locate <-function(x) {

## Define function to solve system of MLE equations
solver <-function(par) {
  
  ## Isolate portion of data corresponding to the unique grouping
  subdata=data.frame(X=x$longitude[x$uniqueID==group],Y=x$latitude[x$uniqueID==group],
                     SIN=x$sin[x$uniqueID==group],COS=x$cos[x$uniqueID==group])
  
  ## Create variables to hold transmitter location estimates
  loc_est = numeric(length(par))
  transmitter_x <- par[1]
  transmitter_y <- par[2]
  
  ## Define system of MLE equations (Lenth 1981)
  loc_est[1] = -sum((transmitter_y-subdata$Y)*(subdata$SIN*(transmitter_x-subdata$X)-
                    subdata$COS*(transmitter_y-subdata$Y))/(((transmitter_x-subdata$X)^2+
                    (transmitter_y-subdata$Y)^2)^0.5)^3)
  loc_est[2] = sum((transmitter_x-subdata$X)*(subdata$SIN*(transmitter_x-subdata$X)-
                  subdata$COS*(transmitter_y-subdata$Y))/(((transmitter_x-subdata$X)^2+
                  (transmitter_y-subdata$Y)^2)^0.5)^3)
  
  ## Return transmitter location estimate
  loc_est  
}

## Convert compass Angle to standard radian measure (Lenth 1981)
x$theta = (pi/180*(90-x$Angle)) 

## Calculate necessary variables
x$sin = sin(x$theta)
x$cos = cos(x$theta)
x$tan = tan(x$theta)

## Calculates intersection of all Angles in each grouping
intersections <-findintersects(x) 

## Create data frame to store calculated transmitter locations
transmitter<-data.frame(X=NA,Y=NA,BadPoint=NA,Var_X=NA,Var_Y=NA,Cov_XY=NA,AngleDiff=NA,
                        Date=NA,Time=NA)

## Declare buffer (in meters) for ensuring location validity  
buffer<-0.00005

## Calculate transmitter location for each grouping
for (group in unique(intersections$uniqueID)) {
  tryCatch({
    
    ## Calculate starting point for solver function based on the mean of Angle intersections
    start_point_x<-mean(intersections$X[intersections$uniqueID==group], na.rm = T) 
    start_point_y<-mean(intersections$Y[intersections$uniqueID==group], na.rm = T) 
  
    
    ## Temporarily store results of solver function
    
    triangulation_results<- nleqslv::nleqslv(c(start_point_x,start_point_y),solver) 
    
    ## Withdraws transmitter location from list of solver function results
    location<-triangulation_results$x
    
    ## Ensure location validity prior to recording
    valid=TRUE
    # if (location[1]<min(intersections$X[intersections$uniqueID==group]-buffer)) valid=FALSE
    # if (location[1]>max(intersections$X[intersections$uniqueID==group]+buffer)) valid=FALSE
    # if (location[2]<min(intersections$Y[intersections$uniqueID==group]-buffer)) valid=FALSE
    # if (location[2]>max(intersections$Y[intersections$uniqueID==group]+buffer)) valid=FALSE
    
    ## Set default color scheme for transmitter ploting
    row.names(transmitter) <- group
    transmitter[group,3]<-0
    
    ## Transfer calculated locations into results variable
    if (valid) {
      transmitter[group,1]<-location[1]
      transmitter[group,2]<-location[2]
    } else {
      warning("Bad point detected in Grouping ",group)
      transmitter[group,1]<-start_point_x
      transmitter[group,2]<-start_point_y
      ## Set color scheme to red for transmitter ploting of bad points
      transmitter[group,3]<-1
    }
    
    ## Create error data frame
    errors=data.frame(X=x$longitude[x$uniqueID==group],Y=x$latitude[x$uniqueID==group],
                      theta=(pi/180*(90-x$Angle[x$uniqueID==group])),si=x$sin[x$uniqueID==group],
                      ci=x$cos[x$uniqueID==group])
    errors$X_hat<-transmitter[group,1]
    errors$Y_hat<-transmitter[group,2]
    errors$di<-((errors$X_hat-errors$X)^2+(errors$Y_hat-errors$Y)^2)^(1/2)
    errors$si_star<-(errors$Y_hat-errors$Y)/(errors$di)^3
    errors$ci_star<-(errors$X_hat-errors$X)/(errors$di)^3
    #errors$mew_i<-atan((errors$Y_hat-errors$Y)/(errors$X_hat-errors$X))
    
    ## Calculate Angle Angles to transmitter location and append to error data frame
    errors$mew_i<-NA
    for (e in 1:nrow(errors)) {
      if ((errors$X[e]==errors$X_hat[e])&&(errors$Y[e]<errors$Y_hat[e])) {errors$mew_i[e]<-0}
      if ((errors$X[e]==errors$X_hat[e])&&(errors$Y[e]>errors$Y_hat[e])) {errors$mew_i[e]<-180}
      if ((errors$Y[e]==errors$Y_hat[e])&&(errors$X[e]<errors$X_hat[e])) {errors$mew_i[e]<-90}
      if ((errors$Y[e]==errors$Y_hat[e])&&(errors$X[e]>errors$X_hat[e])) {errors$mew_i[e]<-270}
      if ((errors$X[e]<errors$X_hat[e])&&(errors$Y[e]<errors$Y_hat[e])) {
        errors$mew_i[e]<-180/pi*atan(abs(errors$X_hat[e]-errors$X[e])/abs(errors$Y_hat[e]-errors$Y[e]))
      }
      if ((errors$X[e]<errors$X_hat[e])&&(errors$Y[e]>errors$Y_hat[e])) {
        errors$mew_i[e]<-180/pi*atan(abs(errors$Y_hat[e]-errors$Y[e])/abs(errors$X_hat[e]-errors$X[e]))+90
      }
      if ((errors$X[e]>errors$X_hat[e])&&(errors$Y[e]>errors$Y_hat[e])) {
        errors$mew_i[e]<-180/pi*atan(abs(errors$X_hat[e]-errors$X[e])/abs(errors$Y_hat[e]-errors$Y[e]))+180
      }
      if ((errors$X[e]>errors$X_hat[e])&&(errors$Y[e]<errors$Y_hat[e])) {
        errors$mew_i[e]<-180/pi*atan(abs(errors$Y_hat[e]-errors$Y[e])/abs(errors$X_hat[e]-errors$X[e]))+270
      }
    }
    
    ## Convert calculated Angles to standard radian measure (Lenth 1981)
    errors$mew_i = (pi/180*(90-errors$mew_i))  
    
    ## Calculate covariance entries as per Lenth (1981)
    C_bar=sum(cos(errors$theta-errors$mew_i)/nrow(errors))
    k_inv=2*(1-C_bar)+(1-C_bar)^2*(.48794-.82905*C_bar-1.3915*C_bar^2)/C_bar
    Q_mat=rbind(c(sum(errors$si_star*errors$si),(-0.5)*sum(errors$si_star*errors$ci+errors$ci_star*errors$si)),c((-0.5)*sum(errors$si_star*errors$ci+errors$ci_star*errors$si),sum(errors$ci_star*errors$ci)))
    Q_hat=k_inv*solve(Q_mat)
    if (transmitter$BadPoint==1) transmitter[group,4]<-0 else transmitter[group,4]<-Q_hat[1,1] #Var_X
    if (transmitter$BadPoint==1) transmitter[group,5]<-0 else transmitter[group,5]<-Q_hat[2,2] #Var_Y
    if (transmitter$BadPoint==1) transmitter[group,6]<-0 else transmitter[group,6]<-Q_hat[1,2] #Cov_XY
    
    ## Calculate average angular difference
    if (transmitter$BadPoint==1) transmitter[group,7]<-0 else transmitter[group,7]=mean(abs(errors$mew_i-errors$theta))*180/pi
    
    ## Attach date and time stamps to results variable
    transmitter[group,8]<-(x$GPSdate[x$uniqueID==group])[1]
    transmitter[group,9]<-(x$GPStime[x$uniqueID==group])[1]
    locateDF <- rbind(locateDF, transmitter)
    
}, error = function(e) e, finally = print(group))}

 
## Return transmitter locations
# class(transmitter)<-"transmitter"
return(locateDF)
}

x2 <- as.receiver(x)
test <- locate(rareUTM)
test$uniqueID <- row.names(test)
test$uniqueID <- gsub("001","01", test$uniqueID)

for(i in 1:nrow(test)) {
  rareUTM[test$uniqueID[i] == rareUTM$uniqueID,"longitude"] <- test[i,"X"]
  rareUTM[test$uniqueID[i] == rareUTM$uniqueID,"latitude"] <- test[i,"Y"]
}

rareUTM ## the done one
write.csv(rareUTM, "2023SpringRareTriangulated.csv")

test <- locate(x)

test[,"uniqueID"] <- row.names(test)

newdata<-left_join(tag1.df2, test, by="uniqueID")

newdata2 <- newdata %>% 
  dplyr::mutate(X = dplyr::coalesce(X,longitude)) %>% 
  dplyr::mutate(Y = dplyr::coalesce(Y,latitude)) 

write.csv(newdata2, "Triangulate.csv")

## add handheld data to the triangulate data

Trian <- read.csv("Triangulate.csv")
hand <- read.csv("HandheldCombinedSpring2022.csv")

toget <- bind_rows(Trian, hand)

write.csv(toget, "TriangulateHandheld.csv")
