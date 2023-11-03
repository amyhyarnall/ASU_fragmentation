# Community analyses for Epibenthic fish and invertebrates (minnow traps), 
# larval fish and crab megalopae (SMURFs), and fish length z-scores (DIDSON) 

# MS Title: Habitat area more consistently influences seagrass faunal communities than fragmentation per se
# Authors: Yarnall AH, Yeager LA, Lopazanski C, Poray AK, Morley J, Fodrie FJ
# Journal: Ecological monographs

# All dataset related to this manuscript are publicly available at:
# https://www.bco-dmo.org/project/714026

# Load libraries
library(lubridate)
library(dplyr)
library(DHARMa)
library(ggplot2)
library(vegan)
library(ecodist)
library(ggvegan)
library(ggrepel)
library(gridExtra)

# Minnow trap, SMURF, and DIDSON z-score data sets were organized to Site_ID by species matrices in excel 
# Read in data sets
setwd("C:/Users/Amy/Documents/UNC_current use files/Ch. 3 NSF ASU fragmentation study/GitHub_Rscripts")
minnow_matrix <- read.csv("asufrag_trapfaunalcpue_communityanalysisformat.csv")
smurf_matrix <- read.csv("asufrag_smurf_settlement_rate_communityanalysisformat.csv")
didson_matrix <- read.csv("asufrag_didson_fish_zscore_communityanalysisformat.csv")

#### Epibenthic fauna (Minnow traps): Data formatting ####
# R-friendly data formatting 
minnow_matrix$Date_Out <- as.Date(minnow_matrix$Date_Out, format = "%m/%d/%Y") # convert to date format
minnow_matrix$month <- month(minnow_matrix$Date_Out, label = T, abbr = F) # create month column

# Remove October samples
minnow_matrix <- minnow_matrix[!minnow_matrix$month == "October",]

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
minnow_matrix$Per_cov100 <- minnow_matrix$Per_cov
minnow_matrix$Per_cov <- minnow_matrix$Per_cov100/100
minnow_matrix <- minnow_matrix[,!colnames(minnow_matrix) == "Per_cov100"]

# Remove A and B designation from Sites 69%-059 A and B
# Later CPUE will be average for these sites
minnow_matrix$Site_ID <- ifelse(grepl("60-0.59", minnow_matrix$Site_ID), 
                             "60-0.59", minnow_matrix$Site_ID)

# There was inconsistency in how Penaeid shrimp were recorded
shrimp.names <- c('Penaeid.Shrimp','Brown.shrimp','White.shrimp')
shrimp.sum <- c(sum(minnow_matrix$Penaeid.Shrimp),sum(minnow_matrix$Brown.shrimp),sum(minnow_matrix$White.shrimp))
names(shrimp.sum) <- shrimp.names; shrimp.sum

# Combine these shrimp columns into 1
minnow_matrix$Penaeid.shrimp <- apply(minnow_matrix[,c(which(colnames(minnow_matrix) %in% shrimp.names))], 1, sum)
# Remove the original shrimp cols from the df
minnow_matrix <- minnow_matrix[,!colnames(minnow_matrix) %in% shrimp.names]

# Remove one column of unknown shrimp ID
minnow_matrix <- minnow_matrix[,!colnames(minnow_matrix) == "Unknown.Shrimp"]

# Create metadata dataframe
minnow.meta <- minnow_matrix[colnames(minnow_matrix) %in% c("Site_ID", "Per_cov", "Frag", "Date_Out", "month", "Trap_class")]

# remove unneeded meta columns and reorder
minnow_matrix <- minnow_matrix[,c(1:4,42,5:40,43)]
minnow_matrix_reduced <- minnow_matrix[,!colnames(minnow_matrix) %in% c("Cell_coord", "Cell_class")]

# Remove minnow traps (rows) that had no catch
minnow_matrix_reduced <- minnow_matrix_reduced[apply(minnow_matrix_reduced[,7:40], 1, function(x) !all(x==0)),]

#### Epibenthic fauna (Minnow traps): Create community matrix (Hellinger transform and Bray Curtis dist) ####

# Vector of all the spp names
minnow.spp.names <- names(minnow_matrix_reduced[7:40])

# Average catch across traps for each landscape across all sampling periods
minnow.df <- minnow_matrix_reduced[,-c(4:6)] %>% group_by(Site_ID, Per_cov, Frag) %>%
  summarise_all(funs(mean))
minnow.comm <- minnow.df[,colnames(minnow.df) %in% minnow.spp.names] #Pull out community data only

# Remove species (columns) that were not caught
minnow.comm <- minnow.comm[,colSums(minnow.comm) > 0]
rownames(minnow.comm) <- minnow.df$Site_ID

# Exploration of raw count means across landscape parameter treatment levels
p.treat <- minnow.df[,-c(1,3)] %>% group_by(Per_cov) %>% summarise_all(funs(mean))
f.treat <- minnow.df[,-c(1,2)] %>% group_by(Frag) %>% summarise_all(funs(mean))

# This matrix can be used to plot vectors for percov and frag over ordination plot 
minnow.env <- minnow.df[,1:3]

# 'Hellinger' transformation - helps standardize for a few abundant spp and many rare spp

# Hellinger y'ij = sqrt(yij/yi+); where I'm using yi+ to indicate the sample total count over all j=1,.,m species, for the ith sample.
minnow.comm_hel <- decostand(minnow.comm, method = 'hellinger')
# Bray curtis distance with ecodist fxs
minnow.comm_hel.bcd <- distance(as.matrix(minnow.comm_hel), method='bray')
# Extended bray curtis distance
minnow.comm_hel.xbcd <- stepacross(minnow.comm_hel.bcd, path = 'shortest', toolong = 1, trace = TRUE)

#### Epibenthic fauna (Minnow traps): NMDS plots  and PERMANOVA ####

# Start community analysis
# First plot the data to see what is looks like
NMDS.minnow <-metaMDS(comm = minnow.comm_hel, distance = 'bray', k = 2, try = 20, 
                      engine = 'monoMDS', 
                      autotransform = FALSE, 
                      noshare = TRUE, 
                      stress = 1, 
                      wascores = TRUE, expand = TRUE, 
                      trace = FALSE, plot = FALSE)
#NMDS.minnow
NMDS.minnow$stress # k = 2; stress < 0.2
stressplot(NMDS.minnow)

#ordiplot(NMDS.minnow, type = "p")
ordiplot(NMDS.minnow, type = "t")

minnow.ord.fit <- envfit(NMDS.minnow ~ Per_cov + Frag, data=minnow.env, perm=999); minnow.ord.fit

plot(NMDS.minnow, type = "t", display = "sites")
plot(minnow.ord.fit)

# Start plotting with ggvegan and ggplot2
#autoplot(NMDS.minnow) #basic plot

minnow.fort <- fortify(NMDS.minnow)
minnow.env.scores <- as.data.frame(scores(minnow.ord.fit, "vectors")) #extracts relevant scores from envfit
minnow.env.scores_labels <- cbind(env.variables = c('Percent cover','Fragmentation')) #and then gives them their names
options(ggrepel.max.overlaps = Inf)

NMDSsite.minnow.plot <- ggplot() + theme_classic() + ylim(-1.5,1) + xlim(-1.5,1) +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', col = 'gray', size = 0.8) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', col = 'gray', size = 0.8) +
  geom_point(data = subset(minnow.fort, Score == 'sites'),
             mapping = aes(x = NMDS1, y = NMDS2), col = 'black', alpha = 0.5) +
  geom_segment(data = minnow.env.scores, size = 1.5,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, col = row.names(minnow.env.scores)), 
               arrow = arrow(length = unit(0.02, "npc"))) + #arrows for envfit
  scale_color_manual(values = c('darkblue','darkgreen'))+
  geom_text_repel(data = subset(minnow.fort, Score == 'sites'), size = 3,
                  mapping = aes(label = Label, x = NMDS1, y = NMDS2)) +
  geom_text(data = NMDS.minnow, x = 0.75, y = 1, col = 'red', size = 4,
            aes(label = paste0("Stress = ", round(NMDS.minnow$stress, digits = 3)))) + 
  geom_text(aes(label = "Percent cover"), x = -0.7, y = 0.1, col = 'darkgreen') +
  geom_text(aes(label = "Fragmentation"), x = -0.1, y = 0.5, col = 'darkblue') +
  theme(legend.position = 'none')

NMDSspecies.minnow.plot <- ggplot() + theme_classic() + ylim(-1.5,1) + xlim(-1.5,1) +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', col = 'gray', size = 0.8) + 
  labs(x = 'NMDS1', y = 'NMDS2')+
  geom_vline(aes(xintercept = 0), linetype = 'dashed', col = 'gray', size = 0.8) +
  geom_segment(data = subset(minnow.fort, Score == 'species'), 
               mapping = aes(x = 0, y= 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.015, 'npc'), type = "closed"), col = 'darkgray', size = 0.8) +
  geom_text_repel(data = subset(minnow.fort, Score == 'species'), size = 3,
                  mapping = aes(label = Label, x = NMDS1, y = NMDS2))+
  geom_segment(data = minnow.env.scores, size = 1.5,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, col = row.names(minnow.env.scores)), 
               arrow = arrow(length = unit(0.02, "npc"))) + #arrows for envfit
  scale_color_manual(values = c('darkblue','darkgreen'))+
  geom_text(aes(label = "Percent cover"), x = -0.7, y = 0.1, col = 'darkgreen') +
  geom_text(aes(label = "Fragmentation"), x = -0.1, y = 0.5, col = 'darkblue') +
  theme(legend.position = 'none')

grid.arrange(NMDSsite.minnow.plot, NMDSspecies.minnow.plot, nrow = 1)

# Run Permanova (adonis2 for Type II SS)
# Treat frag and %cover as numeric 
# Test effect of Per_cov and Frag (and reverse order)
# Note: adonis wants the sample distance matrix and the metadata in different dataframes

adonis.minnowa <- adonis2(minnow.comm_hel.xbcd ~ Per_cov * Frag, data = minnow.env);adonis.minnowa # Nothing is sig
adonis.minnowb <- adonis2(minnow.comm_hel.xbcd ~ Frag * Per_cov, data = minnow.env);adonis.minnowb # Nothing is sig

#### Larval fishes and crab megalopae (SMURFs): Data formatting ####
# R-friendly data formatting 
smurf_matrix$Date_Out <- as.Date(smurf_matrix$Date_Out, format = "%m/%d/%Y") # convert to date format
smurf_matrix$month <- month(smurf_matrix$Date_Out, label = T, abbr = F) # create month column

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
smurf_matrix$Per_cov100 <- smurf_matrix$Per_cov
smurf_matrix$Per_cov <- smurf_matrix$Per_cov100/100
smurf_matrix <- smurf_matrix[,!colnames(smurf_matrix) == "Per_cov100"]

# Remove A and B designation from Sites 69%-059 A and B
# Later counts will be average for these sites
smurf_matrix$Site_ID <- ifelse(grepl("60-0.59", smurf_matrix$Site_ID), 
                                "60-0.59", smurf_matrix$Site_ID)

# Remove columns of unknown fish and crab IDs
smurf_matrix <- smurf_matrix[,!colnames(smurf_matrix) %in% c("Unidentified.crab", "Unided.fish")]

# Brachyuran crabs were not reliably identifiable - pool them
Brachyuran.names <- c("Brachyuran.megalopa","Brachyuran.zoea","Callinectes.sp.")
Brachyuran.sum <- c(sum(smurf_matrix$Brachyuran.megalopa),sum(smurf_matrix$Brachyuran.zoea),sum(smurf_matrix$Callinectes.sp.))
names(Brachyuran.sum) <- Brachyuran.names; Brachyuran.sum

# Combine these Brachyuran columns into 1
smurf_matrix$Brachyura <- apply(smurf_matrix[,c(which(colnames(smurf_matrix) %in% Brachyuran.names))], 1, sum)
# Remove the original Brachyuran cols from the df
smurf_matrix <- smurf_matrix[,!colnames(smurf_matrix) %in% Brachyuran.names]

# reorder meta columns  
smurf_matrix <- smurf_matrix[,c(1:2,21,3:20,22)]

# Remove smurfs (rows) that had no catch - a lot of AUGUST SMURFS
smurf_matrix_reduced <- smurf_matrix[apply(smurf_matrix[,6:22], 1, function(x) !all(x==0)),]

# Create metadata dataframe
smurf.meta <- smurf_matrix_reduced[colnames(smurf_matrix_reduced) %in% c("Site_ID", "Per_cov", "Frag", "Date_Out", "month")]

# Remove species (columns) that were not caught
smurf_matrix_reduced <- smurf_matrix_reduced[,colSums(smurf_matrix_reduced[6:22]) > 0]

#### Larval fishes and crab megalopae (SMURFs): Create community matrix (Hellinger transform and Bray Curtis dist) ####

# Vector of all the spp names
smurf.spp.names <- names(smurf_matrix_reduced[6:22])

# Average catch across traps for each landscape across all sampling periods
smurf.df <- smurf_matrix_reduced[,-c(2:3)] %>% group_by(Site_ID, Per_cov, Frag) %>%
  summarise_all(funs(mean))
smurf.comm <- smurf.df[,colnames(smurf.df) %in% smurf.spp.names] #Pull out community data only

# Remove species (columns) that were not caught
smurf.comm <- smurf.comm[,colSums(smurf.comm) > 0]
rownames(smurf.comm) <- smurf.df$Site_ID

# Exploration of raw count means across landscape parameter treatment levels
p.treat <- smurf.df[,-c(1,3)] %>% group_by(Per_cov) %>% summarise_all(funs(mean))
f.treat <- smurf.df[,-c(1,2)] %>% group_by(Frag) %>% summarise_all(funs(mean))

# This matrix can be used to plot vectors for percov and frag over ordination plot 
smurf.env <- smurf.df[,1:3]

# 'Hellinger' transformation - helps standardize for a few abundant spp and many rare spp

# Hellinger y'ij = sqrt(yij/yi+); where I'm using yi+ to indicate the sample total count over all j=1,.,m species, for the ith sample.
smurf.comm_hel <- decostand(smurf.comm, method = 'hellinger')
# Bray curtis distance with ecodist fxs
smurf.comm_hel.bcd <- distance(as.matrix(smurf.comm_hel), method='bray')
# Extended bray curtis distance
smurf.comm_hel.xbcd <- stepacross(smurf.comm_hel.bcd, path = 'shortest', toolong = 1, trace = TRUE)

#### Larval fishes and crab megalopae (SMURFs): NMDS plots and PERMANOVA ####

# Start community analysis
# First plot the data to see what is looks like
NMDS.smurf <-metaMDS(comm = smurf.comm_hel, distance = 'bray', k = 2, try = 20, 
                      engine = 'monoMDS', 
                      autotransform = FALSE, 
                      noshare = TRUE, 
                      stress = 1, 
                      wascores = TRUE, expand = TRUE, 
                      trace = FALSE, plot = FALSE)
#NMDS.smurf
NMDS.smurf$stress # k = 2; stress < 0.2
stressplot(NMDS.smurf)

#ordiplot(NMDS.smurf, type = "p")
ordiplot(NMDS.smurf, type = "t")

smurf.ord.fit <- envfit(NMDS.smurf ~ Per_cov + Frag, data=smurf.env, perm=999); smurf.ord.fit

plot(NMDS.smurf, type = "t", display = "sites")
plot(smurf.ord.fit)

# Start plotting with ggvegan and ggplot2
#autoplot(NMDS.smurf) #basic plot

smurf.fort <- fortify(NMDS.smurf)
smurf.env.scores <- as.data.frame(scores(smurf.ord.fit, "vectors")) #extracts relevant scores from envfit
smurf.env.scores_labels <- cbind(env.variables = c('Percent cover','Fragmentation')) #and then gives them their names
options(ggrepel.max.overlaps = Inf)

NMDSsite.smurf.plot <- ggplot() + theme_classic() + ylim(-1,1.5) + xlim(-1.5,1) +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', col = 'gray', size = 0.8) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', col = 'gray', size = 0.8) +
  geom_point(data = subset(smurf.fort, Score == 'sites'),
             mapping = aes(x = NMDS1, y = NMDS2), col = 'black', alpha = 0.5) +
  geom_segment(data = smurf.env.scores, size = 1.5,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, col = row.names(smurf.env.scores)), 
               arrow = arrow(length = unit(0.02, "npc"))) + #arrows for envfit
  scale_color_manual(values = c('darkblue','darkgreen'))+
  geom_text_repel(data = subset(smurf.fort, Score == 'sites'), size = 3,
                  mapping = aes(label = Label, x = NMDS1, y = NMDS2)) +
  geom_text(data = NMDS.smurf, x = 0.75, y = 1.5, col = 'red', size = 4,
            aes(label = paste0("Stress = ", round(NMDS.smurf$stress, digits = 3)))) + 
  geom_text(aes(label = "Percent cover"), x = -0.1, y = -0.7, col = 'darkgreen') +
  geom_text(aes(label = "Fragmentation"), x = -0.8, y = 0.2, col = 'darkblue') +
  theme(legend.position = 'none')

NMDSspecies.smurf.plot <- ggplot() + theme_classic() + ylim(-1,1.5) + xlim(-1.5,1) +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', col = 'gray', size = 0.8) + 
  labs(x = 'NMDS1', y = 'NMDS2')+
  geom_vline(aes(xintercept = 0), linetype = 'dashed', col = 'gray', size = 0.8) +
  geom_segment(data = subset(smurf.fort, Score == 'species'), 
               mapping = aes(x = 0, y= 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.015, 'npc'), type = "closed"), col = 'darkgray', size = 0.8) +
  geom_text_repel(data = subset(smurf.fort, Score == 'species'), size = 3,
                  mapping = aes(label = Label, x = NMDS1, y = NMDS2)) +
  geom_segment(data = smurf.env.scores, size = 1.5,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, col = row.names(smurf.env.scores)), 
               arrow = arrow(length = unit(0.02, "npc"))) + #arrows for envfit
  scale_color_manual(values = c('darkblue','darkgreen'))+
  geom_text(aes(label = "Percent cover"), x = -0.1, y = -0.7, col = 'darkgreen') +
  geom_text(aes(label = "Fragmentation"), x = -0.8, y = 0.2, col = 'darkblue') +
  theme(legend.position = 'none')

grid.arrange(NMDSsite.smurf.plot, NMDSspecies.smurf.plot, nrow = 1)

# Run Permanova (adonis2 for Type II SS)
# Treat frag and %cover as numeric 
# Test effect of Per_cov and Frag (and reverse order)
# Note: adonis wants the sample distance matrix and the metadata in different dataframes

adonis.smurfa <- adonis2(smurf.comm_hel.xbcd ~ Per_cov * Frag, data = smurf.env);adonis.smurfa # Nothing is sig
adonis.smurfb <- adonis2(smurf.comm_hel.xbcd ~ Frag * Per_cov, data = smurf.env);adonis.smurfb # Nothing is sig

#### Create ggplot themes and plotting functions for efficient figure creation ####

# top plot theme
toptheme <- theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.y = element_text(size=15),
        axis.ticks = element_line(color = "black", size  =  1),
        axis.line = element_line(color = "black", size = 1))

# side plot theme
sidetheme <- theme_classic() + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size=15, angle = -90, vjust = -0.02),
        axis.title.x = element_text(size=15),axis.ticks = element_line(color = "black", size  =  1),
        axis.line = element_line(color = "black", size = 1)) 

# interaction plot aes and theme
siteaes <- list(geom_point(size = 35,pch = 22, color = "black"),
                scale_fill_continuous(low = "deepskyblue1", high = "navyblue"),
                scale_x_continuous(breaks = c(10,22.5,35,47.5,60), limits = c(5,65)),
                scale_y_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06)),
                labs(y ="Fragmentation per se treatment", x ="Percent cover treatment"),
                coord_fixed(ratio = 60/0.6))

sitetheme <- theme_classic() + 
  theme(axis.text = element_text(size = 15, color = "black", lineheight = 0.9),
        axis.ticks = element_line(color = "black", size  =  1),
        axis.title = element_text(size = 15, color = "black", margin = margin(0, 10, 0, 0)),
        axis.line = element_line(color = "black", size = 1))

# Legend for interaction plot 
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}

# ggpredict function
ggpred <- function(ggpreddata){
  list(geom_line(data = ggpreddata, aes(x = x, y = predicted), col = 'black'),
       geom_ribbon(data = ggpreddata, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey'))
}
######## COMPARE SPECIES RICHNESS & PIELOU'S EVENNESS BETWEEN MINNOW TRAPS AND SMURFS ####

# Minnow trap data
# Examine sampling across sites and months
table(minnow_matrix_reduced$Site_ID, droplevels(minnow_matrix_reduced$month)) # uneven sampling across locations so take the mean per month then sum?

# subset trap data to the same months that SMURFs were used
minnow_matrix_reduced.srj <- minnow_matrix_reduced[minnow_matrix_reduced$month %in% c("June", "July", "August"),]

lpnpip.momean <- minnow_matrix_reduced.srj %>%
  group_by(month, Site_ID, Per_cov, Frag) %>%
  summarise_if(.predicate = function(x) is.numeric(x),
               .funs = funs(mean="mean"))
lpnpip.momean_sum<- lpnpip.momean %>%
  group_by(Site_ID, Per_cov, Frag) %>%
  summarise_if(.predicate = function(x) is.numeric(x),
               .funs = funs(sum="sum"))

trap.H <- diversity(lpnpip.momean_sum[,4:37])
trap.J <- trap.H/log(specnumber(lpnpip.momean_sum[,4:37])) # PIELOU'S EVENNESS
trap.SR <- specnumber(lpnpip.momean_sum[,4:37]) # SPECIES RICHNESS

# Create dataframe of Trap monthly SR and J
trap_SR.J <- lpnpip.momean_sum[,1:3]
trap_SR.J$J <- trap.J
trap_SR.J$SR <- trap.SR


# SMURF data
smurf.momean <- smurf_matrix_reduced %>%
  group_by(month, Site_ID, Per_cov, Frag) %>%
  summarise_if(.predicate = function(x) is.numeric(x),
               .funs = funs(mean='mean'))
smurf.momean_sum <- smurf.momean %>%
  group_by(Site_ID, Per_cov, Frag) %>%
  summarise_if(.predicate = function(x) is.numeric(x),
               .funs = funs(sum='sum'))

smurf.H <- diversity(smurf.momean_sum[,4:20])
smurf.J <- smurf.H/log(specnumber(smurf.momean_sum[,4:20])) # PIELOU'S EVENNESS
smurf.SR <- specnumber(smurf.momean_sum[,4:20]) # SPECIES RICHNESS

# Create dataframe of smurf monthly SR and J
smurf_SR.J <- smurf.momean_sum[,1:3]
smurf_SR.J$J <- smurf.J
smurf_SR.J$SR <- smurf.SR

# Linear regressions for minnow trap and SMURF species richness and evenness

# Minnow traps
# LARGEST PATCH, NEAR PATCH, INTERPATCH DATA
hist(trap_SR.J$SR)
hist(trap_SR.J$J)

# SPECIES RICHNESS
summary(trap.sr.mod <-lm(SR ~ Per_cov*Frag, data = trap_SR.J))# NS
plot(simulateResiduals(fittedModel = trap.sr.mod, n = 250), asFactor = F)#good fit
# PIELOU'S EVENNESS
summary(trap.j.mod <-lm(J ~ Per_cov*Frag, data = trap_SR.J))# NS
plot(simulateResiduals(fittedModel = trap.j.mod, n = 250), asFactor = F)#good fit

# SMURFs
hist(smurf_SR.J$SR)
hist(smurf_SR.J$J)
smurf_SR.J$J[smurf_SR.J$J == "NaN"] <- NA

# SPECIES RICHNESS
summary(smurf.sr.mod <- lm(SR ~ Per_cov*Frag, data = smurf_SR.J))#NS
plot(simulateResiduals(fittedModel = smurf.sr.mod, n = 250), asFactor = F)#good fit
# PIELOU'S EVENNESS
summary(smurf.j.mod <- lm(J ~ Per_cov*Frag, data = smurf_SR.J[!is.na(smurf_SR.J$J),]))#Near sig frag
plot(simulateResiduals(fittedModel = smurf.j.mod, n = 250), asFactor = F)#good fit


# Plot results for Epibenthic faunal SR and J vs Per_cov*Frag 
trap.SR_J.pc <- trap_SR.J %>% group_by(Per_cov) %>% summarise(mean.SR = mean(SR), se.SR = sd(SR)/sqrt(length(SR[!is.na(SR)])), 
                                                                   mean.J = mean(J), se.J = sd(J)/sqrt(length(J[!is.na(J)])))
trap.SR_J.frag <- trap_SR.J %>% group_by(Frag) %>% summarise(mean.SR = mean(SR), se.SR = sd(SR)/sqrt(length(SR[!is.na(SR)])), 
                                                                  mean.J = mean(J), se.J = sd(J)/sqrt(length(J[!is.na(J)])))

new.trap.SR.pc.scatplot <- ggplot() +  ylim(4, 10) +
  geom_errorbar(data = trap.SR_J.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= mean.SR-se.SR, ymax = mean.SR+se.SR), col = "black", width = 0, size = 1)+
  geom_point(data = trap.SR_J.pc, aes(x= as.numeric(as.character(Per_cov)), y= mean.SR), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Epibenthic\nfaunal SR", x = NULL) + toptheme
new.trap.SR.f.scatplot <- ggplot() + ylim(4, 10) +
  geom_errorbar(data = trap.SR_J.frag, aes(x= as.numeric(as.character(Frag)), ymin= mean.SR-se.SR, ymax = mean.SR+se.SR),  col = "black", width = 0, size = 1)+
  geom_point(data = trap.SR_J.frag, aes(x = as.numeric(as.character(Frag)), y= mean.SR), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  labs(y = "Epibenthic\nfaunal SR", x = NULL) +  sidetheme + coord_flip()
new.trap.SR.plot <- ggplot(data = trap_SR.J, aes(x=Per_cov*100, y=Frag, fill = SR) ) +
  labs(fill = "Epibenthic\nfaunal SR")+ siteaes + sitetheme 
new.trap.SR.l <- g_legend(new.trap.SR.plot)
grid.arrange(new.trap.SR.pc.scatplot, new.trap.SR.l, new.trap.SR.plot+theme(legend.position='none'), new.trap.SR.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

new.trap.J.pc.scatplot <- ggplot() +  ylim(0.7,0.9) +
  geom_errorbar(data = trap.SR_J.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= mean.J-se.J, ymax = mean.J+se.J), col = "black", width = 0, size = 1)+
  geom_point(data = trap.SR_J.pc, aes(x= as.numeric(as.character(Per_cov)), y= mean.J), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Epibenthic\nfaunal J", x = NULL) + toptheme
new.trap.J.f.scatplot <- ggplot() + ylim(0.7,0.9) +
  geom_errorbar(data = trap.SR_J.frag, aes(x= as.numeric(as.character(Frag)), ymin= mean.J-se.J, ymax = mean.J+se.J),  col = "black", width = 0, size = 1)+
  geom_point(data = trap.SR_J.frag, aes(x = as.numeric(as.character(Frag)), y= mean.J), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  labs(y = "Epibenthic\nfaunal J", x = NULL) +  sidetheme + coord_flip()
new.trap.J.plot <- ggplot(data = trap_SR.J, aes(x=Per_cov*100, y=Frag, fill = J) ) +
  labs(fill = "Epibenthic\nfaunal J")+ siteaes + sitetheme 
new.trap.J.l <- g_legend(new.trap.J.plot)
grid.arrange(new.trap.J.pc.scatplot, new.trap.J.l, new.trap.J.plot+theme(legend.position='none'), new.trap.J.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))



# Plot results for Settler SR and J vs Per_cov*Frag 
smurf.SR_J.pc <- smurf_SR.J %>% group_by(Per_cov) %>% summarise(mean.SR = mean(SR), se.SR = sd(SR)/sqrt(length(SR[!is.na(SR)])), 
                                                              mean.J = mean(J[!is.na(J)]), se.J = sd(J[!is.na(J)])/sqrt(length(J[!is.na(J)])))
smurf.SR_J.frag <- smurf_SR.J %>% group_by(Frag) %>% summarise(mean.SR = mean(SR), se.SR = sd(SR)/sqrt(length(SR[!is.na(SR)])), 
                                                               mean.J = mean(J[!is.na(J)]), se.J = sd(J[!is.na(J)])/sqrt(length(J[!is.na(J)])))

new.smurf.SR.pc.scatplot <- ggplot() + ylim(2, 7) +
  geom_errorbar(data = smurf.SR_J.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= mean.SR-se.SR, ymax = mean.SR+se.SR), col = "black", width = 0, size = 1)+
  geom_point(data = smurf.SR_J.pc, aes(x= as.numeric(as.character(Per_cov)), y= mean.SR), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Settler SR", x = NULL) + toptheme
new.smurf.SR.f.scatplot <- ggplot() + ylim(2, 7) +
  geom_errorbar(data = smurf.SR_J.frag, aes(x= as.numeric(as.character(Frag)), ymin= mean.SR-se.SR, ymax = mean.SR+se.SR),  col = "black", width = 0, size = 1)+
  geom_point(data = smurf.SR_J.frag, aes(x = as.numeric(as.character(Frag)), y= mean.SR), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  labs(y = "Settler SR", x = NULL) +  sidetheme + coord_flip()
new.smurf.SR.plot <- ggplot(data = smurf_SR.J, aes(x=Per_cov*100, y=Frag, fill = SR) ) +
  labs(fill = "Settler SR")+ siteaes + sitetheme 
new.smurf.SR.l <- g_legend(new.smurf.SR.plot)
grid.arrange(new.smurf.SR.pc.scatplot, new.smurf.SR.l, new.smurf.SR.plot+theme(legend.position='none'), new.smurf.SR.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

new.smurf.J.pc.scatplot <- ggplot() +  ylim(0,0.9) +
  geom_errorbar(data = smurf.SR_J.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= mean.J-se.J, ymax = mean.J+se.J), col = "black", width = 0, size = 1)+
  geom_point(data = smurf.SR_J.pc, aes(x= as.numeric(as.character(Per_cov)), y= mean.J), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Settler J", x = NULL) + toptheme
new.smurf.J.f.scatplot <- ggplot() + ylim(0,0.9) +
  geom_errorbar(data = smurf.SR_J.frag, aes(x= as.numeric(as.character(Frag)), ymin= mean.J-se.J, ymax = mean.J+se.J),  col = "black", width = 0, size = 1)+
  geom_point(data = smurf.SR_J.frag, aes(x = as.numeric(as.character(Frag)), y= mean.J), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  labs(y = "Settler J", x = NULL) +  sidetheme + coord_flip()
new.smurf.J.plot <- ggplot(data = smurf_SR.J, aes(x=Per_cov*100, y=Frag, fill = J) ) +
  labs(fill = "Settler J")+ siteaes + sitetheme 
new.smurf.J.l <- g_legend(new.smurf.J.plot)
grid.arrange(new.smurf.J.pc.scatplot, new.smurf.J.l, new.smurf.J.plot+theme(legend.position='none'), new.smurf.J.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

#### Benthopelagic fish (DIDSON) z-score: Data formatting ####
# rename the bin columns names in r
names(didson_matrix)[6:18] <- c('[-2,-1.5)','[-1.5,-1)','[-1,-0.5)','[-0.5,0)','[0,0.5)','[0.5,1)',
                                   '[1,1.5)','[1.5,2)','[2,2.5)','[2.5,3)','[3,3.5)','[3.5,4)','[4,4.5)')

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
didson_matrix$Per_cov100 <- didson_matrix$Per_cov
didson_matrix$Per_cov <- didson_matrix$Per_cov100/100
didson_matrix <- didson_matrix[,!colnames(didson_matrix) == "Per_cov100"]

# Remove A and B designation from Sites 69%-059 A and B
# Later CPUE will be average for these sites
didson_matrix$Site_ID <- ifelse(grepl("60-0.59", didson_matrix$Site_ID), 
                                "60-0.59", didson_matrix$Site_ID)
# Remove samples (rows) that had no catch
didson_matrix_reduced <- didson_matrix[apply(didson_matrix[,6:18], 1, function(x) !all(x==0)),]

#### Benthopelagic fish (DIDSON) z-score: Create community matrix (Hellinger transform and Bray Curtis dist) ####

# Average catch across videos for each landscape across all sampling periods
didson.df <- didson_matrix_reduced[,-c(4,5)] %>% group_by(Site_ID, Per_cov, Frag) %>%
  summarise_all(funs(mean))
didson.comm <- didson.df[,4:16] #Pull out community data only

# Remove species (columns) that were not caught
didson.comm <- didson.comm[,colSums(didson.comm) > 0]
rownames(didson.comm) <- didson.df$Site_ID

# This matrix can be used to plot vectors for percov and frag over ordination plot 
did.env <- didson.df[,1:3]

# 'Hellinger' transformation - helps standardize for a few abundant spp and many rare spp

# Hellinger y'ij = sqrt(yij/yi+); where I'm using yi+ to indicate the sample total count over all j=1,.,m species, for the ith sample.
didson.comm_hel <- decostand(didson.comm, method = 'hellinger')
# Bray curtis distance with ecodist fxs
didson.comm_hel.bcd <- distance(as.matrix(didson.comm_hel), method='bray')
# Extended bray curtis distance
didson.comm_hel.xbcd <- stepacross(didson.comm_hel.bcd, path = 'shortest', toolong = 1, trace = TRUE)

#### Benthopelagic fish (DIDSON) z-score: NMDS plots & PERMANOVA ####

# Start community analysis
# First plot the data to see what is looks like
NMDS.didson <-metaMDS(comm = didson.comm_hel, distance = 'bray', k = 2, try = 20, 
                        engine = 'monoMDS', 
                        autotransform = FALSE, 
                        noshare = TRUE, 
                        stress = 1, 
                        wascores = TRUE, expand = TRUE, 
                        trace = FALSE, plot = FALSE)
#NMDS.didson
NMDS.didson$stress #lower stress with k = 2; stress < 0.2
stressplot(NMDS.didson)

#ordiplot(NMDS.didson, type = "p")
ordiplot(NMDS.didson, type = "t")

didson.ord.fit <- envfit(NMDS.didson ~ Per_cov + Frag, data=did.env, perm=999); didson.ord.fit

plot(NMDS.didson, type = "t", display = "sites")
plot(didson.ord.fit)

# Start plotting with ggvegan and ggplot2
#autoplot(NMDS.didson) #basic plot
didson.fort <- fortify(NMDS.didson)
didson.env.scores <- as.data.frame(scores(didson.ord.fit, "vectors")) #extracts relevant scores from envfit
didson.env.scores_labels <- cbind(env.variables = c('Percent cover','Fragmentation')) #and then gives them their names
options(ggrepel.max.overlaps = Inf)

NMDSsite.didson.plot <- ggplot() + theme_classic() + ylim(-0.6,0.8) + xlim(-0.8,1.5) +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', col = 'gray', size = 0.8) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', col = 'gray', size = 0.8) +
  geom_point(data = subset(didson.fort, Score == 'sites'),
             mapping = aes(x = NMDS1, y = NMDS2), col = 'black', alpha = 0.5) +
  geom_segment(data = didson.env.scores, size = 1.5,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, col = row.names(didson.env.scores)), 
               arrow = arrow(length = unit(0.02, "npc"))) + #arrows for envfit
  scale_color_manual(values = c('darkblue','darkgreen'))+
  geom_text_repel(data = subset(didson.fort, Score == 'sites'), size = 3,
                  mapping = aes(label = Label, x = NMDS1, y = NMDS2)) +
  geom_text(data = NMDS.didson, x = 1.2, y = 0.8, col = 'red', size = 4,
            aes(label = paste0("Stress = ", round(NMDS.didson$stress, digits = 3)))) + 
  geom_text(aes(label = "Percent cover"), x = -0.6, y = 0.1, col = 'darkgreen') +
  geom_text(aes(label = "Fragmentation"),  x = 0, y = -0.5, col = 'darkblue') +
  theme(legend.position = 'none')

NMDSspecies.didson.plot <- ggplot() + theme_classic() + ylim(-0.6,0.8) + xlim(-0.8,1.5) +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed', col = 'gray', size = 0.8) + 
  labs(x = 'NMDS1', y = 'NMDS2')+
  geom_vline(aes(xintercept = 0), linetype = 'dashed', col = 'gray', size = 0.8) +
  geom_segment(data = subset(didson.fort, Score == 'species'), 
               mapping = aes(x = 0, y= 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.015, 'npc'), type = "closed"), col = 'darkgray', size = 0.8) +
  geom_text_repel(data = subset(didson.fort, Score == 'species'), size = 3,
                  mapping = aes(label = Label, x = NMDS1, y = NMDS2)) +
  geom_segment(data = didson.env.scores, size = 1.5,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, col = row.names(didson.env.scores)), 
               arrow = arrow(length = unit(0.02, "npc"))) + #arrows for envfit
  scale_color_manual(values = c('darkblue','darkgreen'))+
  geom_text(aes(label = "Percent cover"), x = -0.6, y = 0.1, col = 'darkgreen') +
  geom_text(aes(label = "Fragmentation"),  x = 0, y = -0.5, col = 'darkblue') +
  theme(legend.position = 'none')

grid.arrange(NMDSsite.didson.plot, NMDSspecies.didson.plot, nrow = 1)

# Run Permanova (adonis2 for Type II SS)
# Treat frag and %cover as numeric (this is what we intended in original design, otherwise we have really low replication)
# Test effect of Per_cov and Frag (and reverse order)
# Note: adonis wants the sample distance matrix and the metadata in different dataframes

# PERMANOVA
adonis.didsona <- adonis2(didson.comm_hel.xbcd ~ Per_cov * Frag, data = did.env);adonis.didsona # Nothing is sig
adonis.didsonb <- adonis2(didson.comm_hel.xbcd ~ Frag * Per_cov, data = did.env);adonis.didsonb # Nothing is sig

#### Benthopelagic fish (DIDSON) z-score: Distribution across landscape parameters ####
# Minnow trap data
# Examine sampling across sites and months
table(as.factor(didson_matrix_reduced$Site_ID), droplevels(as.factor(didson_matrix_reduced$month))) # uneven sampling across locations so take the mean per month then sum?

didson_lpnpip.momean <- didson_matrix_reduced %>%
  group_by(month, Site_ID, Per_cov, Frag) %>%
  summarise_if(.predicate = function(x) is.numeric(x),
               .funs = funs(mean="mean"))
didson_lpnpip.momean_sum<- didson_lpnpip.momean %>%
  group_by(Site_ID, Per_cov, Frag) %>%
  summarise_if(.predicate = function(x) is.numeric(x),
               .funs = funs(sum="sum"))

didson.H <- diversity(didson_lpnpip.momean_sum[,4:16])
didson.J <- didson.H/log(specnumber(didson_lpnpip.momean_sum[,4:16])) # PIELOU'S EVENNESS
didson.SR <- specnumber(didson_lpnpip.momean_sum[,4:16]) # SPECIES RICHNESS

# Create dataframe of didson monthly SR and J
didson_SR.J <- didson_lpnpip.momean_sum[,1:3]
didson_SR.J$J <- didson.J
didson_SR.J$SR <- didson.SR

# Linear regressions for minnow trap and SMURF species richness and evenness

# Minnow traps
# LARGEST PATCH, NEAR PATCH, INTERPATCH DATA
hist(didson_SR.J$SR)
hist(didson_SR.J$J)

# SPECIES RICHNESS
summary(didson.sr.mod <-lm(SR ~ Per_cov*Frag, data = didson_SR.J))# NS
plot(simulateResiduals(fittedModel = didson.sr.mod, n = 250), asFactor = F)#good fit
# PIELOU'S EVENNESS
summary(didson.j.mod <-lm(J ~ Per_cov*Frag, data = didson_SR.J))# NS
plot(simulateResiduals(fittedModel = didson.j.mod, n = 250), asFactor = F)#good fit

######## Results for MS ####

# Linear regressions for minnow trap and SMURF species richness and evenness

# Minnow traps
# SPECIES RICHNESS
summary(trap.sr.mod)# NS
grid.arrange(new.trap.SR.pc.scatplot, new.trap.SR.l, new.trap.SR.plot+theme(legend.position='none'), new.trap.SR.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# PIELOU'S EVENNESS
summary(trap.j.mod)# NS
grid.arrange(new.trap.J.pc.scatplot, new.trap.J.l, new.trap.J.plot+theme(legend.position='none'), new.trap.J.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))


# SMURFs
# SPECIES RICHNESS
summary(smurf.sr.mod)#NS
grid.arrange(new.smurf.SR.pc.scatplot, new.smurf.SR.l, new.smurf.SR.plot+theme(legend.position='none'), new.smurf.SR.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# PIELOU'S EVENNESS
summary(smurf.j.mod)#Near sig frag
grid.arrange(new.smurf.J.pc.scatplot, new.smurf.J.l, new.smurf.J.plot+theme(legend.position='none'), new.smurf.J.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))


# NMDS AND PERMANOVA 
# Minnow traps
grid.arrange(NMDSsite.minnow.plot, NMDSspecies.minnow.plot, nrow = 1)
adonis.minnowa
adonis.minnowb

# SMURFS
grid.arrange(NMDSsite.smurf.plot, NMDSspecies.smurf.plot, nrow = 1)
adonis.smurfa
adonis.smurfb

# DIDSON
grid.arrange(NMDSsite.didson.plot, NMDSspecies.didson.plot, nrow = 1)
adonis.didsona
adonis.didsonb

