#bayesmove for spring 2023 project

data <- read.csv("2023SpringBothSites.csv")

# Remove the first 3 columns
data <- data[, -c(1,2, 3)]

# View the updated dataset
head(data)

############################################################
###               Bayes move                            ###
#############################################################

library(bayesmove)
library(lubridate)
library(dplyr)
library(tidyr)
library(sf)

#rename columns to work with bayesmove
tagsgdate <-data %>% 
  dplyr::select(longitude, latitude, GPSdatetime, AnimalID, SiteName)
names(tagsgdate)[1] <- "x"
names(tagsgdate)[2] <- "y"
names(tagsgdate)[3] <- "date"
names(tagsgdate)[4] <- "id"
names(tagsgdate)[5] <- "SiteName"
tagsgdate <-tagsgdate %>% arrange(id, date)

#remove any points without dates
tagsNoNA <- tagsgdate %>% filter(!is.na(date))

tagsNoNA$date <- as.POSIXct(tagsNoNA$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# need to separate out the two sites

raredata <- subset(tagsNoNA, SiteName == "Rare")
wallacedata <- subset(tagsNoNA, SiteName == "Wallace")

# Rare

#calculating step length, turning angle and time interval
# there is at least 1 ID with only one occurrence that needs to be removed
# 1. Count occurrences of each ID
id_counts <- raredata %>%
  count(id)

# 2. Filter to find IDs with only one occurrence
single_occurrence_ids <- id_counts %>%
  filter(n == 1) # there is only 1 RFO51907

# remove these IDs before running prep_data
 raredata_filtered <- raredata %>%
 filter(id %in% single_occurrence_ids$id == FALSE)
 
# Then run  prep_data on the filtered data:
tracks.rare <- prep_data(dat = raredata_filtered, coord.names = c("x", "y"), id = "id")
tracks.rare <- replace (tracks.rare, is.na(tracks.rare), 0)
#can't run prep_data on any ids with just one occurrence.

head(tracks.rare)
unique (tracks.rare$id) #21 unique tracks

library(ggplot2)

boxplot(log10(tracks.rare$dt))
tracks.big <- tracks.rare %>% filter(dt < 10000 & dt > 1)
tracks.big2 <- tracks.big %>% filter(dt < 2000 & dt > 1)
#use 15 min, 1 hr, and one day as the time intervals

#round time steps to the three intervals
hist(tracks.rare$dt, breaks=100)

tracks.15 <-round_track_time(dat=tracks.rare, id="id", int=450, tol=450, 
                              time.zone="UTC", units="secs")

tracks.15 <- tracks.15 %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )


# Create list from data frame where each element is a different track

tracks.list.15 <- df_to_list(dat = tracks.15, ind = "id")

# Filter observations to time interval

tracks_filt.list.15 <- filter_time(dat.list = tracks.list.15, int = 450)

#### Discretize data streams (put into bins)
angle.bin.lims=c(-3.143, -2.356, -1.571, -0.785,  0, 
                 0.785,  1.571,  2.356,  3.143)
angle.bin.lims

dist.bin.lims = quantile(tracks.15[tracks.15$dt==450,]$step,
                       c(0, 0.1, 0.25, 0.30, 0.50, 0.75, 0.90, 1), na.rm=T) #7 bins
dist.bin.lims

library(tidyverse)
library(purrr)
# Assign bins to observations
tracks_disc.list.15 <- map(tracks_filt.list.15,
                            discrete_move_var,
                            lims = list(dist.bin.lims, angle.bin.lims),
                            varIn = c("step", "angle"),
                            varOut = c("SL", "TA"))

# Only retain id, discretized step length (SL), turning angle (TA), and rest columns
tracks.15.list2 <- map(tracks_disc.list.15,
                        subset,
                        select = c(id, SL, TA))
### run the segmentation model

set.seed(1)

# Define hyperparameter for prior distribution
alpha <- 1

# Set number of iterations for the Gibbs sampler
ngibbs <- 10000

# Set the number of bins used to discretize each data stream
nbins<- c(7,8)  #SL, TA, rest (in the order from left to right in tracks.list2)

progressr::handlers(progressr::handler_progress(clear = FALSE))
future::plan(future::multisession, workers = 3)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.15.res <- segment_behavior(data = tracks.15.list2, 
                               ngibbs = ngibbs, 
                               nbins = nbins,
                                alpha = alpha)

future::plan(future::sequential)  #return to single core

## Determine MAP for selecting breakpoints (maximum a posteriori)
get_MAP <- function(dat, nburn) {
  
  # Subset only numeric columns after burn-in and convert to matrix
  tmp <- as.matrix(dat[, (nburn + 2):ncol(dat)])  # Subset columns after burn-in and ensure matrix format
  
  # Find max LML per row (ID) after burn-in
  MAP.est <- as.integer(apply(tmp, 1, function(x) which.max(x)) + nburn)
  
  return(MAP.est)
}

MAP.est.15 <- get_MAP(dat = dat.15.res$LML[,-1], nburn = 5000)
MAP.est.15
idColumn <- dat.15.res$LML[,1]

brkpts.15 <- get_breakpts(dat = dat.15.res$brkpts, MAP.est = MAP.est.15)

# How many breakpoints estimated per ID?
apply(brkpts.15[,-1], 1, function(x) length(purrr::discard(x, is.na)))

tracks.seg.15 <- assign_tseg(dat = tracks_disc.list.15, 
                              brkpts = brkpts.15)
head(tracks.seg.15)

# Select only id, tseg, SL, TA,
tracks.seg.15.2 <- tracks.seg.15[,c("id","tseg","SL","TA")]

# Summarize observations by track segment
nbins <- c(7,8)

obs <- summarize_tsegs(dat = tracks.seg.15.2, nbins = nbins)

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
res.15 <- cluster_segments(dat=obs, gamma1=gamma1, alpha=alpha,
                            ngibbs=ngibbs, nmaxclust=nmaxclust,
                            nburn=nburn, ndata.types=ndata.types)

### Determine the number of likely behaviour states

# Extract proportions of behaviors per track segment
theta.estim.15 <- extract_prop(res = res.15, ngibbs = ngibbs, nburn = nburn,
                                nmaxclust = nmaxclust)

# Calculate mean proportions per behaviour
(theta.means.15 <- round(colMeans(theta.estim.15), digits = 3))

# Calculate cumulative sum
cumsum(theta.means.15) #first 2 states comprise 99.5% of all observations

# Convert to data frame for ggplot2
theta.estim_df.15 <- theta.estim.15 %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:nmaxclust, names_to = "behaviour", 
               values_to = "prop") %>%
  modify_at("behaviour", factor)

levels(theta.estim_df.15$behaviour) <- 1:nmaxclust

# Plot results
ggplot(theta.estim_df.15, aes(behaviour, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehaviour", y="Proportion of Total Behaviour\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#most are in the first two behaviour states, but mostly just one

#### Classify the states as behaviours

# Extract bin estimates from phi matrix
behav.res.15 <- get_behav_hist(dat = res.15, nburn = nburn, 
                                ngibbs = ngibbs, nmaxclust = nmaxclust,
                                var.names = c("Step Length","Turning Angle"))

# Plot histograms of proportion data
ggplot(behav.res.15, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c("#21908CFF","#440154FF","#FDE725FF",
                               "grey35","grey35","grey35","grey35",
                               "grey35","grey35","grey35","grey35")
                    , guide = "none") +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:7) +
  facet_grid(behav ~ var, scales = "free_x")

# behaviour 1 = Transit only large step lengths, mainly slight TA but a 
#minority are more tortuous or straight 
# behaviour 2 = explore, only large step lengths, mainly tortuous TA,
# some straight and few at less slight turn angles 

theta.estim.long.15 <- expand_behavior(dat = tracks.seg.15, 
                                       theta.estim = theta.estim.15, 
                                        obs = obs,
                                        nbehav = 2, 
                                        behav.names = c("Transit", 
                                                        "Explore"),
                                        behav.order = c(1,2))

##### Assign behavioural states to tracks

# Merge results with original data
tracks.out.15 <- assign_behavior(dat.orig = tracks.15,
                                  dat.seg.list = df_to_list(tracks.seg.15, "id"),
                                  theta.estim.long = theta.estim.long.15,
                                  behav.names = c("Transit", "Explore"))

write.csv(tracks.out.15, "tracksout_ten.csv")
