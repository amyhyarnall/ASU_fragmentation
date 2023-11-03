# Epibenthic fish and macroinvertebrate CPUE (catch per unit effort) in Minnow traps

# MS Title: Habitat area more consistently influences seagrass faunal communities than fragmentation per se
# Authors: Yarnall AH, Yeager LA, Lopazanski C, Poray AK, Morley J, Fodrie FJ
# Journal: Ecological monographs

# All dataset related to this manuscript are publicly available at:
# https://www.bco-dmo.org/project/714026

# Load libraries
library(lubridate)
library(dplyr)
library(ggplot2)
library(MASS)
library(lmtest) 
library(DHARMa)
library(gridExtra)
library(sjPlot)
library(ggeffects)

# Epibenthic fish and invertebrate CPUE - caught in 'minnow' traps
minnow <- read.csv('https://datadocs.bco-dmo.org/file/WW8WwADCDoXgmx/asufrag_trapfaunalcpue.csv')

#### Preliminary Epibenthic fish and invertebrate analyses ####

## R-friendly data formatting 
minnow$Date_Out <- as.Date(minnow$Date_Out) # convert to date format
minnow$month <- month(minnow$Date_Out, label = T, abbr = F) # create month column
minnow$Sp_name <- ifelse(minnow$Sp_name == "", NA, minnow$Sp_name) # make empty strings into NAs
minnow$Sci_name <- ifelse(minnow$Sci_name == "", NA, minnow$Sci_name) # make empty strings into NAs

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
minnow$Per_cov100 <- minnow$Per_cov
minnow$Per_cov <- minnow$Per_cov100/100

# Use Wilcox.test() to:
# 1. Compare minnow trap CPUE in (near-patch) Matrix vs ASU sampling locations (paired by site and month)
# 2. Compare minnow trap CPUE in Sites 69%-059 A and B on the two shoal environments (paired by site and date)

# Wilcoxon Signed-Rank Test
# This test compares two related samples, e.g., paired/matched, or repeated measures on the same samples, 
# to make inferences as to whether the mean ranks of the two related populations differ. 
# It is the non-parametric equivalent of the paired two-sample t test.

# 1. Compare CPUE in near-patch Matrix vs ASU sampling locations (Wilcox rank sum tests)
minnow$num <- ifelse(is.na(minnow$Sp_name),0,1) # denote traps that had no catch

# Create a dataframe that summarizes each individual minnow trap (exclude October -> post-Hurricane Florence)
date.trapcpue <- minnow[!minnow$month == "October",] %>% group_by(Date_Out,month,Site_ID,Per_cov,Frag,Trap_class,Cell_class) %>% 
  summarise(N = sum(num))

# Subset the dataframe for each group
asu_trapsamples <- date.trapcpue[date.trapcpue$Cell_class == "ASU",]
np_trapsamples <- date.trapcpue[date.trapcpue$Trap_class == "near-patch",]
ip_trapsamples <- date.trapcpue[date.trapcpue$Trap_class == "inter-patch",]

# Keep only the "paired" asu/mtrx sample dates - offset by ~ 1 week
asu_trapsamples <- asu_trapsamples[!as.Date(asu_trapsamples$Date_Out) %in% as.Date(c("2018-06-06","2018-06-07","2018-07-27","2018-08-22")),]
np_trapsamples <- np_trapsamples[!np_trapsamples$Date_Out == "2018-09-06",]
ip_trapsamples <- ip_trapsamples[!ip_trapsamples$Date_Out == "2018-09-06",]

# Merge paired trap samples
asu.np_trapsamples <- merge(asu_trapsamples, np_trapsamples, by = c("month", "Site_ID", "Per_cov", "Frag"), all = F)
asu.ip_trapsamples <- merge(asu_trapsamples, ip_trapsamples, by = c("month", "Site_ID", "Per_cov", "Frag"), all = F)

# Histograms
par(mfrow=c(2,1))
# CPUE asu vs near-patch
hist(asu.np_trapsamples$N.x) #asu
hist(asu.np_trapsamples$N.y) #mtrx
# CPUE asu vs inter-patch
hist(asu.ip_trapsamples$N.x) #asu
hist(asu.ip_trapsamples$N.y) #mtrx

#   What sampling location habitat type had greater minnow trap CPUE: ASU or near-patch Matrix?
wilcox.test(asu.np_trapsamples$N.x, asu.np_trapsamples$N.y, paired = T)
median(asu.np_trapsamples$N.x)/median(asu.np_trapsamples$N.y)
# 3-fold higher mean CPUE in ASU over near-patch Mtrx

#   What sampling location habitat type had greater minnow trap CPUE: ASU or inter-patch Matrix?
wilcox.test(asu.ip_trapsamples$N.x, asu.ip_trapsamples$N.y, paired = T)
median(asu.ip_trapsamples$N.x)/median(asu.ip_trapsamples$N.y)
# 2.5-fold higher mean CPUE in ASU over inter-patch Mtrx

# 2. Compare trap CPUE in Sites 69%-059 A and B on the two shoal environments
# Subset data for the comparison
minnow_60.59 <- date.trapcpue[!date.trapcpue$month == "October"& date.trapcpue$Site_ID %in% c("60-0.59A", "60-0.59B"),]
minnow_60.59 <- minnow_60.59[!minnow_60.59$Date_Out == "2018-06-07",]# remove 6/7/2018 as b was not sampled then

# Delineate and merge paired data
a_trapsamples <- minnow_60.59[minnow_60.59$Site_ID == "60-0.59A",]
b_trapsamples <- minnow_60.59[minnow_60.59$Site_ID == "60-0.59B",]
ab_trapsamples <- merge(a_trapsamples, b_trapsamples, by = c("Date_Out", "month", "Per_cov", "Frag", "Trap_class", "Cell_class"), all = F)

# Histograms
hist(ab_trapsamples$N.x) # "60-0.59A"
hist(ab_trapsamples$N.y) # "60-0.59B"
par(mfrow=c(1,1))

#    Are trap CPUEs in Sites 69%-059 A and B on the two shoal environments different
wilcox.test(ab_trapsamples$N.x, ab_trapsamples$N.y, paired = T)
# Not significantly different


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

#### LARGEST PATCH Epibenthic fish and invertebrate CPUE (Minnow Traps) ####
# Subset data to only include largest patch samples prior to Hurricane Florence
minnow.lp_data <- minnow[minnow$Trap_class == "largest patch" & !minnow$month == "October",]

# Remove A and B designation from Sites 69%-059 A and B
# Later CPUE will be average for these sites
minnow.lp_data$Site_ID.new <- ifelse(grepl("60-0.59", minnow.lp_data$Site_ID), 
                                   "60-0.59", minnow.lp_data$Site_ID)

# Are both 69%-059 A and B sampled on all the same dates?
unique(minnow.lp_data$Date_Out[minnow.lp_data$Site_ID == "60-0.59A"]) #dates for site A
unique(minnow.lp_data$Date_Out[minnow.lp_data$Site_ID == "60-0.59B"]) #dates for Site B
# Site B was not sampled on date "2018-06-07"

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

## Minnow trap CPUE (epibenthic fauna = fish + invertebrates)
par(mfrow=c(1,1))
hist(minnow.lp$N)

# GLMs for largest patch faunal CPUE across landscape %cover and frag per se
faunal.lp.glm <- glm(N ~ Per_cov * Frag, data = minnow.lp, family = "poisson")
faunal.lp.glmnb <- glm.nb(N ~ Per_cov * Frag, data = minnow.lp)
lrtest(faunal.lp.glm, faunal.lp.glmnb) # Likelihood ratio test to assess fit of models -> nb better than pois

plot(simulateResiduals(fittedModel = faunal.lp.glmnb, n = 250), asFactor = F)
# ^ there's an issue here with the lower quantile, let's check it out. 
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

simulationOutput <- simulateResiduals(fittedModel = faunal.lp.glmnb, n = 250)
testDispersion(simulationOutput) # No over dispersion
# The is issue is with the quantiles. Is if from a particular predictor? Nope
par(mfrow = c(1,2))
plotResiduals(simulationOutput, minnow.lp$Per_cov, quantreg = T) # look at within group variance for Per_cov
plotResiduals(simulationOutput, minnow.lp$Frag, quantreg = T) # look at within group variance for Frag
# It's unclear what causing this issue but the model seems a good fit to the data otherwise

# Plot the model results
new.minnow.lp <- with(minnow.lp, data.frame(Per_cov = Per_cov, Frag = Frag, N = N))
new.minnow.lp$pred.p <- predict(update(faunal.lp.glmnb, .~. -Frag-Per_cov:Frag), type = "response")
new.minnow.lp$pred.f <- predict(update(faunal.lp.glmnb, .~. -Per_cov-Per_cov:Frag), type = "response")

plot_model(faunal.lp.glmnb, type = "pred", show.ci = TRUE, terms = "Per_cov")
plot_model(faunal.lp.glmnb, type = "pred", show.ci = TRUE, terms = "Frag")

par(mfrow=c(1,1))
lp.site <- new.minnow.lp %>% group_by(Per_cov, Frag) %>% summarise(mean.N = mean(N), se = sd(N)/sqrt(length(N[!is.na(N)])))
lp.pc <- lp.site %>% group_by(Per_cov) %>% summarise(mean.siteN = mean(mean.N), se = sd(mean.N)/sqrt(length(mean.N[!is.na(mean.N)])))
lp.frag <- lp.site %>% group_by(Frag) %>% summarise(mean.siteN = mean(mean.N), se = sd(mean.N)/sqrt(length(mean.N[!is.na(mean.N)])))

# In these plots - lines are predicted by the model
summary(faunal.lp.glmnb)
trap.lp.summmary <- cbind("mod" = rep("trap lp", length(faunal.lp.glmnb)), coef(summary(faunal.lp.glmnb)), "exp.est" = exp(coef(faunal.lp.glmnb)), exp(confint(faunal.lp.glmnb)))

ggpred.lp.p <- ggpredict(faunal.lp.glmnb, "Per_cov")
ggpred.lp.f <- ggpredict(faunal.lp.glmnb, "Frag")

# Plot results for Epibenthic fish and Invert CPUE vs Per_cov*Frag 
new.minnowlp.pc.scatplot <- ggplot() + ylim(0,5) + ggpred(ggpred.lp.p) +
  geom_errorbar(data = lp.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= mean.siteN-se, ymax = mean.siteN+se), col = "black", width = 0, size = 1)+
  geom_point(data = lp.pc, aes(x= as.numeric(as.character(Per_cov)), y= mean.siteN), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Faunal CPUE", x = NULL) + toptheme
new.minnowlp.f.scatplot <- ggplot() + ylim(0,5) + #ggpred(ggpred.lp.f) +
  geom_errorbar(data = lp.frag, aes(x= as.numeric(as.character(Frag)), ymin= mean.siteN-se, ymax = mean.siteN+se),  col = "black", width = 0, size = 1)+
  geom_point(data = lp.frag, aes(x = as.numeric(as.character(Frag)), y= mean.siteN), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  labs(y = "Faunal CPUE", x = NULL) +  sidetheme + coord_flip()
lp.plot <- ggplot(data = lp.site, aes(x=Per_cov*100, y=Frag, fill = mean.N) ) +
  labs(fill = "Faunal\nCPUE")+ siteaes + sitetheme 
lp.l <- g_legend(lp.plot)
grid.arrange(new.minnowlp.pc.scatplot, lp.l, lp.plot+theme(legend.position='none'), new.minnowlp.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))


#### NEAR-PATCH Epibenthic fish and invertebrate CPUE (Minnow Traps) ####
# Subset data to only include near patch samples prior to Hurricane Florence
minnow.np_data <- minnow[minnow$Trap_class == "near-patch" & !minnow$month == "October",]

# Remove A and B designation from Sites 69%-059 A and B
# Later CPUE will be average for these sites
minnow.np_data$Site_ID.new <- ifelse(grepl("60-0.59", minnow.np_data$Site_ID), 
                                     "60-0.59", minnow.np_data$Site_ID)

# Are both 69%-059 A and B sampled on all the same dates? - yes
unique(minnow.np_data$Date_Out[minnow.np_data$Site_ID == "60-0.59A"]) #dates for site A
unique(minnow.np_data$Date_Out[minnow.np_data$Site_ID == "60-0.59B"]) #dates for Site B

# Separate out the minnow traps that had zero catch
minnow.np_datasub <- minnow.np_data[!is.na(minnow.np_data$Sci_name),]
minnow.nas.np <- minnow.np_data[is.na(minnow.np_data$Sci_name),]

# Pull out traps with no catch (because length(unique(x)) would still = 1)
minnow.nas.np <- minnow.nas.np[,colnames(minnow.nas.np) %in% c("Site_ID.new","Date_Out", "Per_cov", "Frag", "Sci_name")]
colnames(minnow.nas.np)[5] <-  "Site_ID" # change Site_ID.new back to Site_ID
minnow.nas.np$N <- rep(0, length(minnow.nas.np$Site_ID))# create zero entries for N

# Find the (non-zero) CPUE of each site on each date
minnow.N.np <- aggregate(Sci_name ~  Site_ID.new + Date_Out + Per_cov + Frag, data = minnow.np_datasub,
                         FUN = function(x) length((x)))
colnames(minnow.N.np) <- c("Site_ID", "Date_Out","Per_cov", "Frag", "N")

# Add back in minnow traps that had no catch - gives full dataset
minnow.np <- rbind(minnow.N.np, minnow.nas.np[colnames(minnow.nas.np) %in% c("Site_ID", "Date_Out",  "Per_cov", "Frag","N")])

# remove one row from the dataset - double zero sample from 60-0.59 A and B

# Divide the N from the 60%-0.59 site by 2 to give the mean N of both sites on each date 
# Only divide by 2 on dates where both were sampled (i.e., not on "2018-06-07")
minnow.np$N <- ifelse(minnow.np$Site_ID == "60-0.59", 
                      minnow.np$N/2, minnow.np$N)
minnow.np$N <- ceiling(minnow.np$N) # Make the 1 non-whole numbers whole to maintain count distribution

# reorder df by date the site
minnow.np <- minnow.np[with(minnow.np, order(minnow.np$Date_Out, minnow.np$Site_ID)),]

# remove one duplicate row from the dataset - double zero sample from 60-0.59 A and B
minnow.np <- minnow.np[!duplicated(minnow.np), ]

# create month column
minnow.np$month <- month(minnow.np$Date_Out, label = T, abbr = F) # add back in month column

## Minnow trap CPUE (epibenthic fauna = fish + invertebrates)
par(mfrow=c(1,1))
hist(minnow.np$N)

# GLMs for near patch faunal CPUE across landscape %cover and frag per se
faunal.np.glm <- glm(N ~ Per_cov * Frag, data = minnow.np, family = "poisson")
faunal.np.glmnb <- glm.nb(N ~ Per_cov * Frag, data = minnow.np)
lrtest(faunal.np.glm, faunal.np.glmnb)# nb = pois, check the fit of both

plot(simulateResiduals(fittedModel = faunal.np.glm, n = 250), asFactor = F)#good fit
plot(simulateResiduals(fittedModel = faunal.np.glmnb, n = 250), asFactor = F)#good fit

# go with NB to keep with LP model
summary(faunal.np.glmnb) # NS

# Plot the model results
par(mfrow=c(1,1))
new.minnow.np <- with(minnow.np, data.frame(Per_cov = Per_cov, Frag = Frag, N = N))
new.minnow.np$pred.p <- predict(update(faunal.np.glmnb, .~. -Frag-Per_cov:Frag), type = "response")
new.minnow.np$pred.f <- predict(update(faunal.np.glmnb, .~. -Per_cov-Per_cov:Frag), type = "response")

plot_model(faunal.np.glmnb, type = "pred", show.ci = TRUE, terms = "Per_cov")
plot_model(faunal.np.glmnb, type = "pred", show.ci = TRUE, terms = "Frag")

np.site <- new.minnow.np %>% group_by(Per_cov, Frag) %>% summarise(mean.N = mean(N), se = sd(N)/sqrt(length(N[!is.na(N)])))
np.pc <- np.site %>% group_by(Per_cov) %>% summarise(mean.siteN = mean(mean.N), se = sd(mean.N)/sqrt(length(mean.N[!is.na(mean.N)])))
np.frag <- np.site %>% group_by(Frag) %>% summarise(mean.siteN = mean(mean.N), se = sd(mean.N)/sqrt(length(mean.N[!is.na(mean.N)])))

# In these plots - lines are predicted by the model
summary(faunal.np.glmnb)
trap.np.summmary <- cbind("mod" = rep("trap np", length(faunal.np.glmnb)), coef(summary(faunal.np.glmnb)), "exp.est" = exp(coef(faunal.np.glmnb)), exp(confint(faunal.np.glmnb)))

ggpred.np.p <- ggpredict(faunal.np.glmnb, "Per_cov")
ggpred.np.f <- ggpredict(faunal.np.glmnb, "Frag")

# Plot results for Epibenthic fish and Invert CPUE vs Per_cov*Frag 
new.minnownp.pc.scatplot <- ggplot() + ylim(0,1.5) + #ggpred(ggpred.np.p) +
  geom_errorbar(data = np.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= mean.siteN-se, ymax = mean.siteN+se), col = "black", width = 0, size = 1)+
  geom_point(data = np.pc, aes(x= as.numeric(as.character(Per_cov)), y= mean.siteN), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Faunal CPUE", x = NULL) + toptheme
new.minnownp.f.scatplot <- ggplot() + ylim(0,1.5) + #ggpred(ggpred.np.f) +
  geom_errorbar(data = np.frag, aes(x= as.numeric(as.character(Frag)), ymin= mean.siteN-se, ymax = mean.siteN+se),  col = "black", width = 0, size = 1)+
  geom_point(data = np.frag, aes(x = as.numeric(as.character(Frag)), y= mean.siteN), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  labs(y = "Faunal CPUE", x = NULL) +  sidetheme + coord_flip()
np.plot <- ggplot(data = np.site, aes(x=Per_cov*100, y=Frag, fill = mean.N) ) +
  labs(fill = "Faunal\nCPUE")+ siteaes + sitetheme 
np.l <- g_legend(np.plot)
grid.arrange(new.minnownp.pc.scatplot, np.l, np.plot+theme(legend.position='none'), new.minnownp.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

#### INTERPATCH Epibenthic fish and invertebrate CPUE (Minnow Traps) ####
# Subset data to only include near patch samples prior to Hurricane Florence
minnow.ip_data <- minnow[minnow$Trap_class == "inter-patch" & !minnow$month == "October",]

# NO COMBINATION OF 60%-0.59 SITES SAMPLED IN THIS DATASET 
# 0.59 FRAG SITES DID NOT RECEIVE 'INTER-PATCH' TRAPS

# Separate out the minnow traps that had zero catch
minnow.ip_datasub <- minnow.ip_data[!is.na(minnow.ip_data$Sci_name),]
minnow.nas.ip <- minnow.ip_data[is.na(minnow.ip_data$Sci_name),]

# Pull out traps with no catch (because length(unique(x)) would still = 1)
minnow.nas.ip <- minnow.nas.ip[,colnames(minnow.nas.ip) %in% c("Site_ID","Date_Out", "Per_cov", "Frag", "Sci_name")]
minnow.nas.ip$N <- rep(0, length(minnow.nas.ip$Site_ID))# create zero entries for N

# Find the (non-zero) CPUE of each site on each date
minnow.N.ip <- aggregate(Sci_name ~  Site_ID + Date_Out + Per_cov + Frag, data = minnow.ip_datasub,
                         FUN = function(x) length((x)))
colnames(minnow.N.ip) <- c("Site_ID", "Date_Out","Per_cov", "Frag", "N")

# Add back in minnow traps that had no catch - gives full dataset
minnow.ip <- rbind(minnow.N.ip, minnow.nas.ip[colnames(minnow.nas.ip) %in% c("Site_ID", "Date_Out",  "Per_cov", "Frag","N")])

# reorder df by date the site
minnow.ip <- minnow.ip[with(minnow.ip, order(minnow.ip$Date_Out, minnow.ip$Site_ID)),]

# create month column
minnow.ip$month <- month(minnow.ip$Date_Out, label = T, abbr = F) # add back in month column

## Minnow trap CPUE (epibenthic fauna = fish + invertebrates)
par(mfrow=c(1,1))
hist(minnow.ip$N)

# GLMs for interpatch faunal CPUE across landscape %cover and frag per se
faunal.ip.glm <- glm(N ~ Per_cov * Frag, data = minnow.ip, family = "poisson")
faunal.ip.glmnb <- glm.nb(N ~ Per_cov * Frag, data = minnow.ip)
lrtest(faunal.ip.glm, faunal.ip.glmnb) # fits of pois = nb, but let's check resid of both

plot(simulateResiduals(fittedModel = faunal.ip.glm, n = 250), asFactor = F)#good fit
plot(simulateResiduals(fittedModel = faunal.ip.glmnb, n = 250), asFactor = F)#good fit 
# ^ upper quantile deviation detected - but combined adjusted quantile test NS - okay

# Go with NB to be cohesive with previous CPUE models
summary(faunal.ip.glmnb) # nearly sig interaction, sig effect of Per_cov and Frag

plot_model(faunal.ip.glmnb, type = "pred", show.ci = TRUE, terms = "Per_cov")
plot_model(faunal.ip.glmnb, type = "pred", show.ci = TRUE, terms = "Frag")

ip.site <- minnow.ip %>% group_by(Per_cov, Frag) %>% summarise(mean.N = mean(N), se = sd(N)/sqrt(length(N[!is.na(N)])))
ip.pc <- ip.site %>% group_by(Per_cov) %>% summarise(mean.siteN = mean(mean.N), se = sd(mean.N)/sqrt(length(mean.N[!is.na(mean.N)])))
ip.frag <- ip.site %>% group_by(Frag) %>% summarise(mean.siteN = mean(mean.N), se = sd(mean.N)/sqrt(length(mean.N[!is.na(mean.N)])))

# In these plots - lines are predicted by the model
summary(faunal.ip.glmnb)
trap.ip.summmary <- cbind("mod" = rep("trap ip", length(faunal.ip.glmnb)), coef(summary(faunal.ip.glmnb)), "exp.est" = exp(coef(faunal.ip.glmnb)), exp(confint(faunal.ip.glmnb)))

ggpred.ip.p <- ggpredict(faunal.ip.glmnb, "Per_cov")
ggpred.ip.f <- ggpredict(faunal.ip.glmnb, "Frag")

# Plot results for Epibenthic fish and Invert CPUE vs Per_cov*Frag 
new.minnowip.pc.scatplot <- ggplot() + ylim(0,2) + ggpred(ggpred.ip.p) +
  geom_errorbar(data = ip.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= mean.siteN-se, ymax = mean.siteN+se), col = "black",  width = 0, size = 1)+
  geom_point(data = ip.pc, aes(x= as.numeric(as.character(Per_cov)), y= mean.siteN), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Faunal CPUE", x = NULL) + toptheme
new.minnowip.f.scatplot <- ggplot() + ylim(0,2) + ggpred(ggpred.ip.f) +
  geom_errorbar(data = ip.frag, aes(x= as.numeric(as.character(Frag)), ymin= mean.siteN-se, ymax = mean.siteN+se), col = "black", width = 0, size = 1)+
  geom_point(data = ip.frag, aes(x = as.numeric(as.character(Frag)), y= mean.siteN), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  labs(y = "Faunal CPUE", x = NULL) +  sidetheme + coord_flip()
ip.plot <- ggplot(data = ip.site, aes(x=Per_cov*100, y=Frag, fill = mean.N) ) +
  labs(fill = "Faunal\nCPUE")+ siteaes + sitetheme 
ip.l <- g_legend(ip.plot)
grid.arrange(new.minnowip.pc.scatplot, ip.l, ip.plot+theme(legend.position='none'), new.minnowip.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

ggpred.ip.int <- ggpredict(faunal.ip.glmnb, c("Per_cov", "Frag"))
interpatch.pred <- ggpredict(faunal.ip.glmnb, terms = c("Per_cov", "Frag"))
ggplot(interpatch.pred, aes(x, predicted, colour = group)) + geom_line()

# Does minnow trap CPUE directly related to the interpatch distance (m)?
inter_dist <- data.frame(Site_ID = unique(minnow.ip$Site_ID), 
                           interDist_m = c(3.6,8.6,9.96,10.8,3.44,4.3,2.1,6.44,1.72,2.55,1.48,6.29,3.44,1.2,3.6,1.72,1.2,1.2,2.58,2.58))
minnow.ip_dist <- merge(minnow.ip, inter_dist, by = "Site_ID", all = T)

faunal.ip_dist.glmnb <- glm.nb(N ~ interDist_m, data = minnow.ip_dist)
plot(simulateResiduals(fittedModel = faunal.ip_dist.glmnb, n = 250), asFactor = F)#good fit 
summary(faunal.ip_dist.glmnb) #NS

#### Results for MS ####

# Preliminary analyses
# Paired Wilcoxon rank sum tests

# 1. Compare CPUE in near-patch Matrix vs ASU sampling locations (Wilcox rank sum tests)
#   What sampling location habitat type had greater minnow trap CPUE: ASU or near-patch Matrix?
wilcox.test(asu.np_trapsamples$N.x, asu.np_trapsamples$N.y, paired = T)
median(asu.np_trapsamples$N.x)/median(asu.np_trapsamples$N.y)
# 3-fold higher mean CPUE in ASU over near-patch Mtrx

#   What sampling location habitat type had greater minnow trap CPUE: ASU or inter-patch Matrix?
wilcox.test(asu.ip_trapsamples$N.x, asu.ip_trapsamples$N.y, paired = T)
median(asu.ip_trapsamples$N.x)/median(asu.ip_trapsamples$N.y)
# 2.5-fold higher mean CPUE in ASU over inter-patch Mtrx

# 2. Compare trap CPUE in Sites 69%-059 A and B on the two shoal environments

#    Are trap CPUEs in Sites 69%-059 A and B on the two shoal environments different
wilcox.test(ab_trapsamples$N.x, ab_trapsamples$N.y, paired = T)
# Not significantly different


# Epibenthic fish and invertebrate CPUE across landscape parameters (Per_cov*Frag) 

# Largest ASU patch
trap.lp.summmary
grid.arrange(new.minnowlp.pc.scatplot, lp.l, lp.plot+theme(legend.position='none'), new.minnowlp.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# Near-patch matrix
trap.np.summmary
grid.arrange(new.minnownp.pc.scatplot, np.l, np.plot+theme(legend.position='none'), new.minnownp.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# Interpatch matrix
trap.ip.summmary
grid.arrange(new.minnowip.pc.scatplot, ip.l, ip.plot+theme(legend.position='none'), new.minnowip.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

