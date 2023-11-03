# Fine-scale complexity influences on settler and post-settler fauna
# Compare landscape-scale vs. fine-scale habitat structure effects

# MS Title: Habitat area more consistently influences seagrass faunal communities than fragmentation per se
# Authors: Yarnall AH, Yeager LA, Lopazanski C, Poray AK, Morley J, Fodrie FJ
# Journal: Ecological monographs

# All dataset related to this manuscript are publicly available at:
# https://www.bco-dmo.org/project/714026

# Load libraries
library(lubridate)
library(dplyr)
library(MASS)
library(DHARMa)
library(MuMIn)
library(rsq)
library(ggeffects)
library(ggplot2)
library(gridExtra)

# ASU fine-scale complexity - canopy height and epiphyte biomass
sg <- read.csv('https://datadocs.bco-dmo.org/file/XYM7M14TONR8mY/asu_frag_landscape_complexity.csv')
# Larval fish and megalopae settlement rates (caught in SMURFs)
smurf <- read.csv('https://datadocs.bco-dmo.org/file/EKPKZqQCgogj25/asufrag_smurf_settlement_rate.csv')
# Epibenthic fish and invertebrate CPUE - caught in 'minnow' traps
minnow <- read.csv('https://datadocs.bco-dmo.org/file/WW8WwADCDoXgmx/asufrag_trapfaunalcpue.csv')
# Benthopelagic fish CPUE (recorded in DIDSON samples)
didson <- read.csv('https://datadocs.bco-dmo.org/file/N7870JESMjKmpp/asufrag_didsonfishcpue.csv')

#### Data formatting for ASU characteristic, SMURF, Minnow trap, and DIDSON data sets ####

### Format ASU characteristic DATA
## R-friendly data formatting 
sg$Date_collected <- as.Date(sg$Date_collected)
sg$month <- month(sg$Date_collected, label = T, abbr = F)

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
sg$Per_cov100 <- sg$Per_cov
sg$Per_cov <- sg$Per_cov100/100

# Remove samples for which the Site_ID info was lost
sg <- sg[!sg$Site_ID == "",]

# Calculate SG surface area
sg$blade.mm2 <- sg$Height*5*2 # blade surface area by multiplying height (mm) by 5 mm width * 2 sides

# Dry epiphyte biomass
sg$dry_epibio <- as.numeric(sg$Dry_wt)-as.numeric(sg$Tin_wt) # dry epiphyte biomass
sg$dry_epibio[sg$dry_epibio < 0] <- NA # remove neg numbers due to measurement errors
sg$dry_epibio <- as.numeric(sg$dry_epibio) # make sure class numeric
sg$dry_epibio_g.mm2 <- sg$dry_epibio/sg$blade.mm2 # dry epiphyte biomass per mm2

# Ash free epiphyte biomass
sg$ashfree_epibio <- as.numeric(sg$Dry_wt)-as.numeric(sg$Ash_wt)# ash-free epiphyte biomass (tin wt included in both dry and ash wt so cancels out)
sg$ashfree_epibio[sg$ashfree_epibio < 0] <- NA # remove neg numbers due to measurement errors
sg$ashfree_epibio <- as.numeric(sg$ashfree_epibio)#make sure class numeric
sg$ashfree_epibio_g.mm2 <- sg$ashfree_epibio/sg$blade.mm2 # ash free epiphyte biomass per mm2 blade SA

# Remove any samples where the ash wt ended up being more than the dry wt as this is not possible - measurement error
sg$ashfree_epibio_g.mm2[sg$ashfree_epibio_g.mm2 > sg$dry_epibio_g.mm2] <- NA

# Convert to units that are more readable
sg$epi_mg.cm2 <- sg$ashfree_epibio_g.mm2*100*1000

# Estimate biomass per blade
sg$epibio.mg.blade <- sg$ashfree_epibio_g.mm2*sg$blade.mm2*1000 
# ^ mean mg of ash free epiphyte biomass per blade

# Are canopy (blade) height and epiphyte biomass correlated? Yes 
cor.test(sg$Height, sg$epi_mg.cm2, paired = T)

par(mfrow = c(3,1))
hist(sg$Height)#pretty normal
hist(log(sg$epibio.mg.blade[!is.na(sg$epibio.mg.blade)]))#log-normal? gamma? Has two humps....
hist(log(sg$epi_mg.cm2[!is.na(sg$epi_mg.cm2)]))#log-normal? gamma? Has two humps....
par(mfrow = c(1,1))

# Some discrepancy in Sitea-IDs for 60%-0.59 sites
sg_60.59 <- sg[sg$Per_cov == 0.60 & sg$Frag == 0.59,]
table(sg_60.59$month, sg_60.59$Site_ID) # A and B were not designated for July samples
table(sg_60.59$Site_ID, sg_60.59$Cell_coord) # Let's try to use Cell_coord for separation
# ^ will need to remove  C3R7 from sg data to match up the 60%-0.59 site data

# removed C3R7 from sg dataset
sg_60.59 <- sg[sg$Per_cov == 0.60 & sg$Frag == 0.59,]
table(sg_60.59$month, sg_60.59$Site_ID) # A and B were not designated for July samples
table(sg_60.59$Site_ID, sg_60.59$Cell_coord) # Let's try to use Cell_coord for separation
# ^ will need to remove  C3R7 from sg data to match up the 60%-0.59 site data

# removed C3R7 from sg dataset
sg <- sg[!(sg$Site_ID == "60-0.59" & sg$Cell_coord == "C3R7"),] 
sg_60.59 <- sg_60.59[!(sg_60.59$Site_ID == "60-0.59" & sg_60.59$Cell_coord == "C3R7"),] 

# How does fine-scale complexity differ between 60%-0.59 A and B sites?
sg_60.59$Site_ID.new <- ifelse(sg_60.59$Site_ID == "60-0.59" & sg_60.59$Date_scraping == "2019-02-05","60-0.59B",
                              ifelse(sg_60.59$Site_ID == "60-0.59" & !sg_60.59$Date_scraping == "2019-02-05", "60-0.59A", sg_60.59$Site_ID)) 
# used context clues from metadata to infer which site was which
sg_60.59$ht_a <- ifelse(sg_60.59$Site_ID.new == "60-0.59A", sg_60.59$Height, NA)
sg_60.59$ht_b <- ifelse(sg_60.59$Site_ID.new == "60-0.59B", sg_60.59$Height, NA)
sg_60.59$eb_a <- ifelse(sg_60.59$Site_ID.new == "60-0.59A", sg_60.59$epi_mg.cm2, NA)
sg_60.59$eb_b <- ifelse(sg_60.59$Site_ID.new == "60-0.59B", sg_60.59$epi_mg.cm2, NA)
sg_60.59$hteb_a <- ifelse(sg_60.59$Site_ID.new == "60-0.59A", sg_60.59$epibio.mg.blade, NA)
sg_60.59$hteb_b <- ifelse(sg_60.59$Site_ID.new == "60-0.59B", sg_60.59$epibio.mg.blade, NA)

sg_ab <- sg_60.59 %>% group_by(month, Cell_coord) %>% summarise(hta = mean(ht_a, na.rm = T),htb = mean(ht_b, na.rm = T),
                                                                eba = mean(eb_a, na.rm = T),ebb = mean(eb_b, na.rm = T),
                                                                hteba = mean(hteb_a, na.rm = T),htebb = mean(hteb_b, na.rm = T))

par(mfrow = c(3,2))
hist(sg_ab$hta)
hist(sg_ab$htb)
hist(sg_ab$eba)
hist(sg_ab$ebb)
hist(sg_ab$hteba)
hist(sg_ab$htebb)
par(mfrow = c(1,1))

t.test(sg_ab$hta, sg_ab$htb, paired = T) #SIG
t.test(sg_ab$eba, sg_ab$ebb, paired = T) #SIG
t.test(sg_ab$hteba, sg_ab$htebb, paired = T) #SIG
# differ between sites - however faunal responses did not differ so we will still take mean for 60%-0.59 a and b

# As further background info: How many ASUs were sampled per site? 
sg_asu.count<- sg %>% group_by(Site_ID) %>% summarise(asu_num = length(unique(Cell_coord)))
sg_asu.count <- sg_asu.count[-c(26,27),] #remove the duplicate "60-0.59" sites - will through off the mean

# number of asus sampled per landscape
mean(sg_asu.count$asu_num) #6.8
sd(sg_asu.count$asu_num) # 3.3

### Format SMURF DATA
## R-friendly data formatting 
smurf$Date_Out <- as.Date(smurf$Date_Out) # convert to date format
smurf$Sci_name <- ifelse(smurf$Sci_name == "", NA, smurf$Sci_name) # make empty strings into NAs

# Create site list
smurf.sites <- unique(smurf[,c("Site_ID","Latitude","Longitude","Per_cov","Frag","Date_Out")])

# Select only just settled individuals and larvae
smurf.set <- smurf[smurf$Settler == "Y",]

# Merge with site list to add back in samples with zero settlers
smurf.settlers <- merge(smurf.sites, smurf.set, all = TRUE)
smurf.settlers$N[is.na(smurf.settlers$N)] <- 0
smurf.settlers$month <- month(smurf.settlers$Date_Out, label = T, abbr = F) # create month column

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
smurf.settlers$Per_cov100 <- smurf.settlers$Per_cov
smurf.settlers$Per_cov <- smurf.settlers$Per_cov100/100

# Remove A and B designation from Sites 69%-059 A and B
smurf.settlers$Site_ID <- ifelse(grepl("60-0.59", smurf.settlers$Site_ID), 
                                 "60-0.59", smurf.settlers$Site_ID)

# Sum settler abundance across all sampling dates for each landscape
smurf.set.sum <- aggregate(N ~ Site_ID + Per_cov + Frag + Date_Out + month, data = smurf.settlers, FUN = "sum")

# Divide the N from the 60%-0.59 site by 2 to give the mean catch of both sites on each date 
# only divide by 2 on dates where both were sampled (i.e., not on "6/8/2018")
smurf.set.sum$N <- ifelse(smurf.set.sum$Site_ID == "60-0.59" & !smurf.set.sum$Date_Out == "2018-06-08", 
                          smurf.set.sum$N/2, smurf.set.sum$N)

# Make the one non-whole number whole to maintain count distribution (poisson or neg. binomial)
smurf.set.sum$N <- ceiling(smurf.set.sum$N) # this fxn rounds up to the next whole integer


### Format Minnow trap DATA
## R-friendly data formatting 
minnow$Date_Out <- as.Date(minnow$Date_Out) # convert to date format
minnow$month <- month(minnow$Date_Out, label = T, abbr = F) # create month column
minnow$Sp_name <- ifelse(minnow$Sp_name == "", NA, minnow$Sp_name) # make empty strings into NAs
minnow$Sci_name <- ifelse(minnow$Sci_name == "", NA, minnow$Sci_name) # make empty strings into NAs

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
minnow$Per_cov100 <- minnow$Per_cov
minnow$Per_cov <- minnow$Per_cov100/100

# Subset data to only include largest patch samples prior to Hurricane Florence
minnow.lp_data <- minnow[minnow$Trap_class == "largest patch" & !minnow$month == "October",]

# Remove A and B designation from Sites 69%-059 A and B
# Later CPUE will be average for these sites
minnow.lp_data$Site_ID.new <- ifelse(grepl("60-0.59", minnow.lp_data$Site_ID), 
                                     "60-0.59", minnow.lp_data$Site_ID)

# Separate out the minnow traps that had zero catch
minnow.lp_datasub <- minnow.lp_data[!is.na(minnow.lp_data$Sci_name),]
minnow.nas.lp <- minnow.lp_data[is.na(minnow.lp_data$Sci_name),]

# Pull out traps with no catch (because length(unique(x)) would still = 1)
minnow.nas.lp <- minnow.nas.lp[,colnames(minnow.nas.lp) %in% c("Site_ID.new","Date_Out", "Per_cov", "Frag", "Sci_name")]
colnames(minnow.nas.lp)[5] <-  "Site_ID" # change Site_ID.new back to Site_ID
minnow.nas.lp$N <- rep(0, length(minnow.nas.lp$Site_ID))# create zero entries for N

# Find the (non-zero) CPUE of each site on each date
minnow.N.lp <- aggregate(Sci_name ~  Site_ID.new + Date_Out + Per_cov + Frag, data = minnow.lp_datasub,
                         FUN = function(x) length((x)))
colnames(minnow.N.lp) <- c("Site_ID", "Date_Out","Per_cov", "Frag", "N")

# Divide the N from the 60%-0.59 site by 2 to give the mean N of both sites on each date 
# Only divide by 2 on dates where both were sampled (i.e., not on "2018-06-07")
minnow.N.lp$N <- ifelse(minnow.N.lp$Site_ID == "60-0.59" & !(as.character(minnow.N.lp$Date_Out) == "2018-06-07"), 
                        minnow.N.lp$N/2, minnow.N.lp$N)
minnow.N.lp$N <- ceiling(minnow.N.lp$N) # Make the 2 non-whole numbers whole to maintain count distribution

# Add back in minnow traps that had no catch - gives full dataset
minnow.lp <- rbind(minnow.N.lp, minnow.nas.lp[colnames(minnow.nas.lp) %in% c("Site_ID", "Date_Out",  "Per_cov", "Frag","N")])

# reorder df by date the site
minnow.lp <- minnow.lp[with(minnow.lp, order(minnow.lp$Date_Out, minnow.lp$Site_ID)),]

# create month column
minnow.lp$month <- month(minnow.lp$Date_Out, label = T, abbr = F) # add back in month column

### Format DIDSON DATA
## R-friendly data formatting 
didson$Date <- as.Date(didson$Date) # convert to date format

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
didson$Per_cov100 <- didson$Per_cov
didson$Per_cov <- didson$Per_cov100/100

# Add column for fish count (1 or 0) because "no fish" were given 0 length
didson$N <- ifelse(didson$L_cm == 0, 0, 1)

# Remove MaxN frames - not needed
didson <- didson[didson$MaxN == "n",]

# Delete sample where it was noted that the wrong frames were read by one of the readers 
didson <- didson[!(didson$Month == 'July' & didson$Site_ID == '60-0.475' & didson$Class == "inter-patch" & didson$Reader == 'AY'),]

# Create a dataframe giving the mean N and length of fish per frame and in each sample (i.e., recording)
didson_FrameN <- didson %>% group_by(Class, Month, Site_ID, Per_cov, Frag, Reader, Frames_sampled, FrameSample) %>%
  summarise(FrameN = sum(N))
didson_SampleN <- didson_FrameN %>% group_by(Class, Month, Site_ID, Per_cov, Frag, Reader, Frames_sampled) %>%
  summarise(SampleN = sum(FrameN), muFrameN = mean(FrameN)) 

# Rename df for further manipulation
didson_sub <- didson_SampleN

# Remove A and B designation from Sites 69%-059 A and B
didson_sub$Site_ID <- ifelse(grepl("60-0.59", didson_Sample$Site_ID), 
                             "60-0.59", didson_Sample$Site_ID)
didson_sub <- didson_sub[!didson_sub$Month == "October",]

# Sum "Frames_sampled", "SampleN" across all sampling dates for each landscape
didson_subsum <- didson_sub %>% group_by(Class,Month,Site_ID,Per_cov,Frag)%>%
  summarise(sum_Frames_sampled=sum(Frames_sampled),sum_N=sum(SampleN))

# Divide sum_N by (sum_Frames_sampled/10) - get average number of fish seen per 10 frames
didson_subsum$N_10frames <- didson_subsum$sum_N/(didson_subsum$sum_Frames_sampled/10)

# Make the non-whole numbers whole to maintain count distribution (poisson or neg. binomial)
didson_subsum$N <- round(didson_subsum$N_10frames)

# How many samples did this change? What % of the dataset? (9 samples - 6% of the dataset)
length(didson_subsum$N[didson_subsum$N != didson_subsum$N_10frames])/length(didson_subsum$N)

# Create the dataset for all the following analyses
didson_data <- subset(didson_subsum, select = -c(sum_Frames_sampled, sum_N, N_10frames))

# Subset the dataframe for the largest patch data
didsonlp.sub <- didson_data[didson_data$Class == "largest patch",]
didsonlp.sub$month <- didsonlp.sub$Month

#### ASU characteristics vs. faunal density ####

# Subset the SMURF and Trap datasets to match the ASU sg samples
smurf.set.sub <- smurf.set.sum[smurf.set.sum$month == "July",]
minnowlp.sub <- minnow.lp[minnow.lp$Date_Out == '2018-07-06' & minnow.lp$month == "July",]

# summarize each dataset to the catch/SR per site in each month
smurf.sg.sub <- smurf.set.sub %>% group_by(month,Per_cov,Frag) %>% summarise(sm.N = sum(N)) 
minnow.sg.sub <- minnowlp.sub %>% group_by(month,Per_cov,Frag) %>% summarise(mt.N = sum(N))
didson.sg.sub <- didsonlp.sub %>% group_by(month,Per_cov,Frag) %>% summarise(did.N = sum(N))

# This line takes care of using the mean on 60%-0.59 A and B sg data
sg.sum <- sg %>% group_by(month,Per_cov,Frag) %>% summarise(mu.ht = mean(Height, na.rm = T), mu.epi.bio = mean(epi_mg.cm2, na.rm = T), 
                                                            mu.ht.epi = mean(epibio.mg.blade, na.rm = T),)

# merge all the dataframes into one
sg.didson <-  merge(sg.sum, didson.sg.sub, all = F)
sg.didson$site <- paste0(sg.didson$Per_cov, "_", sg.didson$Frag)
sg.didson.j <- sg.didson[sg.didson$month == "July",]

sg.totcatch <- merge(merge(merge(
  sg.sum[sg.sum$month == "July",],
  smurf.sg.sub, all = TRUE),
  minnow.sg.sub, all = TRUE),
  sg.didson.j, all = TRUE)# 47.5%-0.225 was not sampled by DIDSON
sg.totcatch$site <- paste0(sg.totcatch$Per_cov, "_", sg.totcatch$Frag)


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
                scale_x_continuous(breaks = c(0.10,0.225,0.35,0.475,0.60), limits = c(0.05,0.65)),
                scale_y_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06)),
                labs(y ="Fragmentation per se treatment", x ="Percent cover treatment"),
                coord_fixed(ratio = 0.6/0.6))

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

#### ASU characteristics across landscapes (Figures) ####
sg.site <- sg[sg$month == "July",] %>% group_by(Per_cov, Frag) %>% 
  summarise(epibio = mean(epi_mg.cm2, na.rm = T), se.epibio = sd(epi_mg.cm2, na.rm = T)/sqrt(length(epi_mg.cm2)), 
            ht = mean(Height), se.ht = sd(Height)/sqrt(length(Height)), 
            htepi = mean(epibio.mg.blade, na.rm = T), se.htepi = sd(epibio.mg.blade, na.rm = T)/sqrt(length(epibio.mg.blade)))
sg.pc <- sg.site %>% group_by(Per_cov) %>% 
  summarise(epibio.mu = mean(epibio, na.rm = T), se.epibio = sd(epibio, na.rm = T)/sqrt(length(epibio)),
            ht.mu = mean(ht), se.ht = sd(ht)/sqrt(length(ht)), 
            htepi.mu = mean(htepi), se.htepi = sd(htepi)/sqrt(length(htepi)))
sg.frag <- sg.site %>% group_by(Frag) %>% 
  summarise(epibio.mu = mean(epibio, na.rm = T), se.epibio = sd(epibio, na.rm = T)/sqrt(length(epibio)),
            ht.mu = mean(ht), se.ht = sd(ht)/sqrt(length(ht)), 
            htepi.mu = mean(htepi), se.htepi = sd(htepi)/sqrt(length(htepi)))

# MODELING
# #*************************************
hist(sg.site$epibio)
hist(sg.site$ht)
hist(sg.site$htepi)

epi_lm <- lm(epibio ~ Per_cov * Frag, data = sg.site)
plot(simulateResiduals(fittedModel = epi_lm, n = 250), asFactor = F)# good fit
summary(epi_lm)

ht_lm <- lm(ht ~ Per_cov * Frag, data = sg.site)
plot(simulateResiduals(fittedModel = ht_lm, n = 250), asFactor = F)# good fit
summary(ht_lm)

htepi_lm <- lm(htepi ~ Per_cov * Frag, data = sg.site)
plot(simulateResiduals(fittedModel = htepi_lm, n = 250), asFactor = F)# good fit
summary(htepi_lm)
# #*************************************

epibio.pc.plot <- ggplot(sg.pc, aes(x= as.numeric(as.character(Per_cov)), y= epibio.mu)) + 
  geom_errorbar(aes(ymin = epibio.mu-se.epibio, ymax = epibio.mu+se.epibio), position = position_dodge(0.9), width = 0, size = 1) +
  geom_point(data = sg.pc, aes(x= as.numeric(as.character(Per_cov)), y= epibio.mu), shape = 21, fill = "lightslategrey", size = 5) +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  scale_y_continuous(breaks = seq(from = 1, to = 5, by = 1)) +
  labs(y = "Epiphyte biomass", x = NULL) + theme_classic() + toptheme
epibio.f.plot <- ggplot(sg.frag, aes(x= as.numeric(as.character(Frag)), y= epibio.mu)) + 
  geom_errorbar(aes(ymin = epibio.mu-se.epibio, ymax = epibio.mu+se.epibio), position = position_dodge(0.9), width = 0, size = 1) +
  geom_point(data = sg.frag, aes(x= as.numeric(as.character(Frag)), y= epibio.mu), shape = 21, fill = "lightslategrey", size = 5) +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  scale_y_continuous(breaks = seq(from = 1, to = 5, by = 1)) +
  labs(y = "Epiphyte biomass", x = NULL) + theme_classic() + sidetheme + coord_flip()
epibio.plot <- ggplot(data = sg.site, aes(y=Per_cov, x=Frag, fill = epibio) ) +
  labs(fill = "Epiphyte\nbiomass") + scale_fill_continuous(low = "#CCFF99", high = "#006600") + 
  theme_classic() + siteaes + sitetheme 
epibio.l <- g_legend(epibio.plot)
grid.arrange(epibio.pc.plot, epibio.l, epibio.plot+theme(legend.position='none'), epibio.f.plot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

ht.pc.plot <- ggplot(sg.pc, aes(x= as.numeric(as.character(Per_cov)), y= ht.mu)) + 
  geom_errorbar(aes(ymin = ht.mu-se.ht, ymax = ht.mu+se.ht), position = position_dodge(0.9), width = 0, size = 1) +
  geom_point(data = sg.pc, aes(x= as.numeric(as.character(Per_cov)), y= ht.mu), shape = 21, fill = "lightslategrey", size = 5) +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  scale_y_continuous(breaks = seq(from = 80, to = 120, by = 10), limits = c(80,120))+
  labs(y = "Canopy height", x = NULL) + theme_classic() + toptheme
ht.f.plot <- ggplot(sg.frag, aes(x= as.numeric(as.character(Frag)), y= ht.mu)) + 
  geom_errorbar(aes(ymin = ht.mu-se.ht, ymax = ht.mu+se.ht), position = position_dodge(0.9), width = 0, size = 1) +
  geom_point(data = sg.frag, aes(x= as.numeric(as.character(Frag)), y= ht.mu), shape = 21, fill = "lightslategrey", size = 5) +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  scale_y_continuous(breaks = seq(from = 80, to = 120, by = 10), limits = c(80,120))+
  labs(y = "Canopy height", x = NULL) + theme_classic() + sidetheme + coord_flip()
ht.plot <- ggplot(data = sg.site, aes(y=Per_cov, x=Frag, fill = ht)) +
  labs(fill = "Canopy\nheight") + scale_fill_continuous(low = "#CCFF99", high = "#006600") + 
  theme_classic() + siteaes + sitetheme 
ht.l <- g_legend(ht.plot)
grid.arrange(ht.pc.plot, ht.l, ht.plot+theme(legend.position='none'), ht.f.plot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

htepi.pc.plot <- ggplot(sg.pc, aes(x= as.numeric(as.character(Per_cov)), y= htepi.mu)) + 
  geom_errorbar(aes(ymin = htepi.mu-se.htepi, ymax = htepi.mu+se.htepi), position = position_dodge(0.9), width = 0, size = 1) +
  geom_point(data = sg.pc, aes(x= as.numeric(as.character(Per_cov)), y= htepi.mu), shape = 21, fill = "lightslategrey", size = 5) +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  #scale_y_continuous(breaks = seq(from = 80, to = 120, by = 10), limits = c(80,120))+
  labs(y = "Epiphyte biomass per blade", x = NULL) + theme_classic() + toptheme
htepi.f.plot <- ggplot(sg.frag, aes(x= as.numeric(as.character(Frag)), y= htepi.mu)) + 
  geom_errorbar(aes(ymin = htepi.mu-se.htepi, ymax = htepi.mu+se.htepi), position = position_dodge(0.9), width = 0, size = 1) +
  geom_point(data = sg.frag, aes(x= as.numeric(as.character(Frag)), y= htepi.mu), shape = 21, fill = "lightslategrey", size = 5) +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  #scale_y_continuous(breaks = seq(from = 80, to = 120, by = 10), limits = c(80,120))+
  labs(y = "Epiphyte biomass per blade", x = NULL) + theme_classic() + sidetheme + coord_flip()
htepi.plot <- ggplot(data = sg.site, aes(y=Per_cov, x=Frag, fill = htepi)) +
  labs(fill = "Epiphyte biomass\nper blade") + scale_fill_continuous(low = "#CCFF99", high = "#006600") + 
  theme_classic() + siteaes + sitetheme 
htepi.l <- g_legend(htepi.plot)
grid.arrange(htepi.pc.plot, htepi.l, htepi.plot+theme(legend.position='none'), htepi.f.plot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))
# Examine canopy height and epiphyte biomass separately - abiotic vs biotic, also very different scales and ranges

#### Landscape- and fine-scale relative importance MODELS ####

# no. fauna caught in smurfs in july vs landscape-scale and fine-scale parameters 
smurf.sg <- sg.totcatch[-22,] #(exlcude outlier smurf)
hist(smurf.sg$sm.N)
smurflpsg.glmnb <- glm.nb(sm.N ~ Per_cov*Frag  + mu.ht + mu.epi.bio, data =smurf.sg)
plot(simulateResiduals(fittedModel = smurflpsg.glmnb, n = 250), asFactor = F)#good fit
summary(smurflpsg.glmnb)

# rsq.partial(smurflpsg.glmnb)# Frag(.), mu.ht(*), Per_cov:Frag(.)
# options(na.action = "na.fail")
# smurf.aicc <- dredge(smurflpsg.glmnb)# top 3 models (deltaAICc<2) include mu.ht (1 includes frag, 1 includes per_cov) 
# setwd("C:/Users/Amy/Documents/UNC_current use files/Ch. 3 NSF ASU fragmentation study/Ecological Monographs/ECM23-0015/Figures")
# write.table(smurf.aicc, file = "SMURF_AICc_landscapeVSfine.csv", row.names = T, sep = ';')

# Independent effects on epibio and ht on smurf catch
smurfht.glmnb <- glm.nb(sm.N ~ mu.ht, data =smurf.sg)
plot(simulateResiduals(fittedModel = smurfht.glmnb, n = 250), asFactor = F)#good fit
summary(smurfht.glmnb) # sig neg
smurf.sg$pred.sm.ht <- predict(smurfht.glmnb, type = "response")
# rsq(smurfht.glmnb)

smurfeb.glmnb <- glm.nb(sm.N ~ mu.epi.bio, data =smurf.sg)
plot(simulateResiduals(fittedModel = smurfeb.glmnb, n = 250), asFactor = F)#good fit
summary(smurfeb.glmnb) # sig neg
smurf.sg$pred.sm.eb <- predict(smurfeb.glmnb, type = "response")
# rsq(smurfeb.glmnb)

# mean no. fauna caught in 2 LP traps in july vs landscape-scale and fine-scale parameters
hist(sg.totcatch$mt.N)
mtlpsg.glmnb <- glm.nb(mt.N ~ Per_cov*Frag + mu.ht + mu.epi.bio, data = sg.totcatch)
plot(simulateResiduals(fittedModel = mtlpsg.glmnb, n = 250), asFactor = F)#good fit
summary(mtlpsg.glmnb)# per_cov(**), Per_cov:Frag(.)

# rsq.partial(mtlpsg.glmnb)
# options(na.action = "na.fail")
# trap.aicc <- dredge(mtlpsg.glmnb)
# setwd("C:/Users/Amy/Documents/UNC_current use files/Ch. 3 NSF ASU fragmentation study/Ecological Monographs/ECM23-0015/Figures")
# write.table(trap.aicc, file = "TRAP_AICc_landscapeVSfine.csv", row.names = T, sep = ';')

# Independent effects on epibio and ht on trap catch
mtht.glmnb <- glm.nb(mt.N ~ mu.ht, data =sg.totcatch)
plot(simulateResiduals(fittedModel = mtht.glmnb, n = 250), asFactor = F)#good fit
summary(mtht.glmnb) # sig pos
sg.totcatch$pred.mt.ht <- predict(mtht.glmnb, type = "response")
# rsq(mtht.glmnb)

mteb.glmnb <- glm.nb(mt.N ~ mu.epi.bio, data =sg.totcatch)
plot(simulateResiduals(fittedModel = mteb.glmnb, n = 250), asFactor = F)#good fit
summary(mteb.glmnb) # NS
sg.totcatch$pred.mt.eb <- predict(mteb.glmnb, type = "response")
# rsq(mteb.glmnb)

# tot fish counted per frame in july LP vids vs landscape-scale and fine-scale parameters (exlcude outlier video)
sg.didson.jsub <- sg.totcatch[-17,]
hist(sg.didson.jsub$did.N)
didsonlpsg.glmnb <- glm.nb(did.N ~ Per_cov*Frag + mu.ht + mu.epi.bio, data = sg.didson.jsub, link = log)
plot(simulateResiduals(fittedModel = didsonlpsg.glmnb, n = 250), asFactor = F)
summary(didsonlpsg.glmnb)#NS

# rsq.partial(didsonlpsg.glmnb)#per_cov is the only non-neg rsq
# options(na.action = "na.fail")
# didson.aicc <- dredge(didsonlpsg.glmnb)
# setwd("C:/Users/Amy/Documents/UNC_current use files/Ch. 3 NSF ASU fragmentation study/Ecological Monographs/ECM23-0015/Figures")
# write.table(didson.aicc, file = "DIDSON_AICc_landscapeVSfine.csv", row.names = T, sep = ';')

# Independent effects on epibio and ht on didson catch
didsonht.glmnb <- glm.nb(did.N ~ mu.ht, data =sg.didson.jsub)
plot(simulateResiduals(fittedModel = didsonht.glmnb, n = 250), asFactor = F)#good fit
summary(didsonht.glmnb) # NS
sg.didson.jsub$pred.did.ht <- predict(didsonht.glmnb, type = "response")
# rsq(didsonht.glmnb)

didsoneb.glmnb <- glm.nb(did.N ~ mu.epi.bio, data =sg.didson.jsub)
plot(simulateResiduals(fittedModel = didsoneb.glmnb, n = 250), asFactor = F)#good fit
summary(didsoneb.glmnb) # NS
sg.didson.jsub$pred.did.eb <- predict(didsoneb.glmnb, type = "response")
# rsq(didsoneb.glmnb)

#### Marginal effects model predictions figure  (fine-scale influences) ####
# FULL model predictions graphs 
# ggpred.a2 <- ggpredict(mtlpsg.glmnb, "mu.ht")
ggpred.b2 <- ggpredict(mtlpsg.glmnb, "mu.epi.bio")
# ggpred.c2 <- ggpredict(didsonlpsg.glmnb, "mu.ht")
# ggpred.d2 <- ggpredict(didsonlpsg.glmnb, "mu.epi.bio")
ggpred.e2 <- ggpredict(smurflpsg.glmnb, "mu.ht")
# ggpred.f2 <- ggpredict(smurflpsg.glmnb, "mu.epi.bio")

a2.plot <- ggplot(sg.totcatch) + geom_point(aes(x= mu.ht, y= mt.N)) + ylim(-1,10) +
  # geom_line(data = ggpred.a2, aes(x = x, y = predicted), col = 'black')+ 
  # geom_ribbon(data = ggpred.a2, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x= NULL, y= "Faunal CPUE") + theme_classic() + mytheme + theme(axis.text.x = element_blank())
b2.plot <- ggplot(sg.totcatch) + geom_point(aes(x= mu.epi.bio, y= mt.N)) + ylim(-1,10) +
  geom_line(data = ggpred.b2, aes(x = x, y = predicted), col = 'black')+ 
  geom_ribbon(data = ggpred.b2, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x= NULL, y= NULL) + theme_classic() + mytheme + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
c2.plot <- ggplot(sg.didson.jsub) + geom_point(aes(x= mu.ht, y= did.N)) + ylim(0,20) +
  # geom_line(data = ggpred.c2, aes(x = x, y = predicted), col = 'black')+ 
  # geom_ribbon(data = ggpred.c2, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x = NULL, y= "Fish CPUE") + theme_classic() + mytheme + theme(axis.text.x = element_blank())
d2.plot <- ggplot(sg.didson.jsub) + geom_point(aes(x= mu.epi.bio, y= did.N)) + ylim(0,20) +
  # geom_line(data = ggpred.d2, aes(x = x, y = predicted), col = 'black')+ 
  # geom_ribbon(data = ggpred.d2, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x = NULL, y= NULL) + theme_classic() + mytheme + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
e2.plot <- ggplot(smurf.sg) + geom_point(aes(x= mu.ht, y= sm.N)) + ylim(0,50) +
  geom_line(data = ggpred.e2, aes(x = x, y = predicted), col = 'black')+ 
  geom_ribbon(data = ggpred.e2, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x= "Mean canopy height", y= "Settlement rate") + theme_classic() + mytheme 
f2.plot <- ggplot(smurf.sg) + geom_point(aes(x= mu.epi.bio, y= sm.N)) + ylim(0,50) + 
  # geom_line(data = ggpred.f2, aes(x = x, y = predicted), col = 'black')+ 
  # geom_ribbon(data = ggpred.f2, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x= "Mean epiphyte biomass", y= NULL) + theme_classic() + mytheme + theme(axis.text.y = element_blank())

grid.arrange(a2.plot, b2.plot, c2.plot, d2.plot, e2.plot, f2.plot, ncol = 2, nrow = 3, widths=c(1.2,1), heights=c(1,1,1))

#### 1-way model predictions figure  (fine-scale influences) ####
require(ggeffects)
# Basic plot theme
mytheme <- theme(legend.position= "bottom", 
                 legend.text=element_text(size=20), 
                 legend.title = element_text(size=20), 
                 axis.text.x = element_text(colour = "black", size=20), 
                 axis.text.y = element_text(colour = "black", size=20), 
                 axis.title.x = element_text(size=20), 
                 axis.title.y = element_text(size=20), 
                 strip.text.x = element_text(colour = "black", size = 20))

# 1-way model prediction graphs 
ggpred.a <- ggpredict(mtht.glmnb, "mu.ht")
# ggpred.b <- ggpredict(mteb.glmnb, "mu.epi.bio")
# ggpred.c <- ggpredict(didsonht.glmnb, "mu.ht")
# ggpred.d <- ggpredict(didsoneb.glmnb, "mu.epi.bio")
ggpred.e <- ggpredict(smurfht.glmnb, "mu.ht")
ggpred.f <- ggpredict(smurfeb.glmnb, "mu.epi.bio")

b.plot <- ggplot(sg.totcatch) + geom_point(aes(x= mu.epi.bio, y= mt.N)) + ylim(-1,10) + ggtitle("a") +
  # geom_line(data = ggpred.b, aes(x = x, y = predicted), col = 'black')+
  # geom_ribbon(data = ggpred.b, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x= NULL, y= "Epibenthic fauna CPUE") + theme_classic() + mytheme + theme(axis.text.x = element_blank())
a.plot <- ggplot(sg.totcatch) + geom_point(aes(x= mu.ht, y= mt.N)) + ylim(-1,10) + ggtitle("b") +
  geom_line(data = ggpred.a, aes(x = x, y = predicted), col = 'black')+ 
  geom_ribbon(data = ggpred.a, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x= NULL, y= NULL) + theme_classic() + mytheme + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
d.plot <- ggplot(sg.didson.jsub) + geom_point(aes(x= mu.epi.bio, y= did.N)) + ylim(0,20) + ggtitle("c") +
  # geom_line(data = ggpred.d, aes(x = x, y = predicted), col = 'black')+ 
  # geom_ribbon(data = ggpred.d, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x = NULL, y= "Benthopelagic fish CPUE") + theme_classic() + mytheme + theme(axis.text.x = element_blank())
c.plot <- ggplot(sg.didson.jsub) + geom_point(aes(x= mu.ht, y= did.N)) + ylim(0,20) + ggtitle("d") +
  # geom_line(data = ggpred.c, aes(x = x, y = predicted), col = 'black')+ 
  # geom_ribbon(data = ggpred.c, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x = NULL, y= NULL) + theme_classic() + mytheme + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
f.plot <- ggplot(smurf.sg) + geom_point(aes(x= mu.epi.bio, y= sm.N)) + ylim(0,50) + ggtitle("e") +
  geom_line(data = ggpred.f, aes(x = x, y = predicted), col = 'black')+ 
  geom_ribbon(data = ggpred.f, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x= "Mean epiphyte biomass", y= "Settler CPUE") + theme_classic() + mytheme 
e.plot <- ggplot(smurf.sg) + geom_point(aes(x= mu.ht, y= sm.N)) + ylim(0,50) + ggtitle("f") +
  geom_line(data = ggpred.e, aes(x = x, y = predicted), col = 'black')+ 
  geom_ribbon(data = ggpred.e, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = 'lightslategrey') +
  labs(x= "Mean canopy height", y= NULL) + theme_classic() + mytheme + theme(axis.text.y = element_blank())

grid.arrange(b.plot, a.plot, d.plot, c.plot, f.plot, e.plot,  ncol = 2, nrow = 3, widths=c(1.2,1), heights=c(1,1,1))

