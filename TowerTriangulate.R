####    code for triangulating the radio tower and motus data

##  NOTE *** I have been triangulating periodically to help with field work only using tower data

# first need to remove any duplicate values from the radio tower data
library(tidyverse)
tag1 <- read.csv("Spring2023TagsProcessed.csv") #787 observations
tag2 <- tag1 %>%  distinct(Antenna, GPSdate, Northing, Easting, TagID, AnimalID, GPStime, 
                           .keep_all=TRUE) #5991 observations

# next, need to join motus and radio tower data - motus will have lots of false detections
library(lubridate)
motus <- read.csv("motus2022.csv") #8761 obs
motus$month <-month(motus$GPSdate) 
motus2<- subset(motus, month < 7) #7013 obs

motag <- plyr::rbind.fill (tag2, motus2)

#code to triangulate bumble bee location. Taken from package sigloc
#cannot use package directly as it has been removed from CRAN

separateTime <- function(time){
  timeDF <- stringr::str_split(time,":", simplify=T)
  timeDF <- apply(timeDF, 2, as.numeric)
  timeDF <- data.frame(timeDF)
  names(timeDF) <- c("hours","minutes","seconds")
  return(timeDF)
}
tagsTime <- separateTime(tag1$GPStime)
tag1[,"uniqueID"] <- paste0(tag1proj$AnimalID, "-",tag1proj$GPSdate,"-",
                            tagsTime$hours)

## need utm for triangulation code
regcrs <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

land2021 <- rast("2021landcoverRaster.tif") #just using this file for the CRS
tag1sp <- vect(tag1, geom=c("longitude", "latitude"), crs=regcrs, geomkeep=FALSE)
tag1proj <- project(tag1sp, land2021)


as.receiver<-function(x) {
  class(x)<-"receiver"
  return(x)
}
findintersects<-function(x) {
  
  ## Convert compass Angle to standard radian measure (Lenth 1981)
  x$theta = (pi/180*(90-x$Angle)) 
  
  ## Calculate additional variables
  x$sin = sin(x$theta)
  x$cos = cos(x$theta)
  x$tan = tan(x$theta)
  
  ## Create data frame to store coordinates of all intersections
  results<-data.frame(uniqueID=NA,X=NA,Y=NA) 
  
  ## Keep track of total number of intersections in the data
  counter=1; 
  for (group in unique(x$uniqueID)) {
    
    ## Create subdata to store data of the single UniqueID
    subdata=data.frame(X=x$Easting[x$uniqueID==group],Y=x$Northing[x$uniqueID==group],
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
x <- tag1proj

locate<-function(x) {

## Define function to solve system of MLE equations
solver<-function(par) {
  
  ## Isolate portion of data corresponding to the unique grouping
  subdata=data.frame(X=x$Easting[x$uniqueID==group],Y=x$Northing[x$uniqueID==group],
                     SIN=x$sin[x$uniqueID==group],COS=x$cos[x$uniqueID==group])
  
  ## Create variables to hold transmitter location estimates
  loc_est = numeric(length(par))
  transmitter_x<-par[1]
  transmitter_y<-par[2]
  
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
intersections<-findintersects(x) 

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
    if (location[1]<min(intersections$X[intersections$uniqueID==group]-buffer)) valid=FALSE
    if (location[1]>max(intersections$X[intersections$uniqueID==group]+buffer)) valid=FALSE
    if (location[2]<min(intersections$Y[intersections$uniqueID==group]-buffer)) valid=FALSE
    if (location[2]>max(intersections$Y[intersections$uniqueID==group]+buffer)) valid=FALSE
    
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
    errors=data.frame(X=x$Easting[x$uniqueID==group],Y=x$Northing[x$uniqueID==group],
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

test <- locate(x)


test[,"uniqueID"] <- row.names(test)

newdata<-left_join(tag1proj, test, by="uniqueID")
write.csv(newdata, "Triangulate.csv")

## add handheld data to the triangulate data

Trian <- read.csv("Triangulate.csv")
hand <- read.csv("HandheldCombinedSpring2022.csv")

toget <- bind_rows(Trian, hand)

write.csv(toget, "TriangulateHandheld.csv")
