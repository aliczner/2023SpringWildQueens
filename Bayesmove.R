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

tracks.ten <-round_track_time(dat=tracks.rare, id="id", int=302400, tol=302400, 
                              time.zone="UTC", units="secs")

tracks.ten <- tracks.ten %>%
  mutate(rest = case_when(step > 0 ~ 2,
                          step == 0 ~ 1)
  )


# Create list from data frame where each element is a different track

tracks.list.ten <- df_to_list(dat = tracks.ten, ind = "id")

# Filter observations to time interval

tracks_filt.list.ten <- filter_time(dat.list = tracks.list.ten, int = 302400)

#will go forward with 10 minute, min, day. 

#### Discretize data streams (put into bins)
angle.bin.lims=c(-3.143, -2.356, -1.571, -0.785,  0, 0.785,  1.571,  2.356,  3.143)
angle.bin.lims

tracks.ten_norest <- subset(tracks.ten, rest == 2)
dist.bin.lims=quantile(tracks.ten_norest[tracks.ten_norest$dt==300,]$step,
                       c(0, 0.1, 0.25, 0.30, 0.50, 0.75, 0.90, 1), na.rm=T) #7 bins

dist.bin.lims

library(tidyverse)
library(purrr)
# Assign bins to observations
tracks_disc.list.ten <- map(tracks_filt.list.ten,
                            discrete_move_var,
                            lims = list(dist.bin.lims, angle.bin.lims),
                            varIn = c("step", "angle"),
                            varOut = c("SL", "TA"))

# Since 0s get lumped into bin 1 for SL, need to add an extra bin to store 0s
tracks_disc.list.ten2 <- tracks_disc.list.ten %>%
  map(., ~mutate(.x, SL = SL + 1)) %>%  #to shift bins over
  map(., ~mutate(.x, SL = case_when(step == 0 ~ 1,  #assign 0s to bin 1
                                    step != 0 ~ SL)  #otherwise keep the modified SL bin
  ))

# Only retain id, discretized step length (SL), turning angle (TA), and rest columns
tracks.ten.list2 <- map(tracks_disc.list.ten2,
                        subset,
                        select = c(id, SL, TA))



### run the segmentation model

set.seed(1)

# Define hyperparameter for prior distribution
alpha<- 1

# Set number of iterations for the Gibbs sampler
ngibbs<- 10000

# Set the number of bins used to discretize each data stream
nbins<- c(8,8)  #SL, TA, rest (in the order from left to right in tracks.list2)

progressr::handlers(progressr::handler_progress(clear = FALSE))
future::plan(future::multisession, workers = 3)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.ten.res <- segment_behavior(data = tracks.ten.list2, ngibbs = ngibbs, nbins = nbins,
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

MAP.est.ten <- get_MAP(dat = dat.ten.res$LML[,-1], nburn = 5000)
MAP.est.ten
idColumn <- dat.ten.res$LML[,1]

brkpts.ten <- get_breakpts(dat = dat.ten.res$brkpts, MAP.est = MAP.est.ten)

# How many breakpoints estimated per ID?
apply(brkpts.ten[,-1], 1, function(x) length(purrr::discard(x, is.na)))

tracks.seg.ten <- assign_tseg(dat = tracks_disc.list.ten2, 
                              brkpts = brkpts.ten)
head(tracks.seg.ten)

# Select only id, tseg, SL, TA,
tracks.seg.ten2 <- tracks.seg.ten[,c("id","tseg","SL","TA")]


# Summarize observations by track segment
nbins <- c(8,8)

obs <- summarize_tsegs(dat = tracks.seg.ten2, nbins = nbins)

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
res.ten <- cluster_segments(dat=obs, gamma1=gamma1, alpha=alpha,
                            ngibbs=ngibbs, nmaxclust=nmaxclust,
                            nburn=nburn, ndata.types=ndata.types)

### Determine the number of likely behaviour states

# Extract proportions of behaviors per track segment
theta.estim.ten <- extract_prop(res = res.ten, ngibbs = ngibbs, nburn = nburn,
                                nmaxclust = nmaxclust)

# Calculate mean proportions per behaviour
(theta.means.ten <- round(colMeans(theta.estim.ten), digits = 3))

# Calculate cumulative sum
cumsum(theta.means.ten) #first 2 states comprise 97.4% of all observations

# Convert to data frame for ggplot2
theta.estim_df.ten<- theta.estim.ten %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:nmaxclust, names_to = "behaviour", 
               values_to = "prop") %>%
  modify_at("behaviour", factor)

levels(theta.estim_df.ten$behaviour)<- 1:nmaxclust

# Plot results
ggplot(theta.estim_df.ten, aes(behaviour, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehaviour", y="Proportion of Total Behaviour\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
#most are in the first two behaviour states

#### Classify the states as behaviours

# Extract bin estimates from phi matrix
behav.res.ten <- get_behav_hist(dat = res.ten, nburn = nburn, 
                                ngibbs = ngibbs, nmaxclust = nmaxclust,
                                var.names = c("Step Length","Turning Angle"))

# Plot histograms of proportion data
ggplot(behav.res.ten, aes(x = bin, y = prop, fill = as.factor(behav))) +
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

# behaviour 1 = MPC, Minimal Positional Change majority zero step lengths, 
#with mainly straight but some turns
# behaviour 2 = exploratory, lots of differnt SL and lots of turns, larger SL
# behaviour 3 = rest but is very rare, not sure if it can be tested statistically.

theta.estim.long.ten <- expand_behavior(dat = tracks.seg.ten, theta.estim = theta.estim.ten, 
                                        obs = obs,
                                        nbehav = 3, 
                                        behav.names = c("MPC", "Explore", "ARS"),
                                        behav.order = c(1,2,3))

##### Assign behavioural states to tracks

# Merge results with original data
tracks.out.ten <- assign_behavior(dat.orig = tracks.ten,
                                  dat.seg.list = df_to_list(tracks.seg.ten, "id"),
                                  theta.estim.long = theta.estim.long.ten,
                                  behav.names = c("MPC", "Explore", "ARS"))

write.csv(tracks.out.ten, "tracksout_pesticideten.csv")
