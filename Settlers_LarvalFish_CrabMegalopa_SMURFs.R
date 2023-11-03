# Larval fish and crab megalopae settlement rates caught in SMURFs

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

# Larval fish and megalopae settlement rates (caught in SMURFs)
smurf <- read.csv('https://datadocs.bco-dmo.org/file/EKPKZqQCgogj25/asufrag_smurf_settlement_rate.csv')

####  Preliminary Settlement rate analyses ####

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

# Use Wilcox.test() to compare smurf catch in Sites 69%-059 A and B on the two shoal environments (paired by site and date)

# Subset the dataframe for each group
smurf_60.59 <- smurf.settlers[smurf.settlers$Site_ID %in% c("60-0.59A", "60-0.59B"),] %>% group_by(Site_ID, Date_Out, Per_cov, Frag) %>% summarise(N = sum(N))

# Delineate and merge paired data
a_smurfsamples <- smurf_60.59[smurf_60.59$Site_ID == "60-0.59A",]
b_smurfsamples <- smurf_60.59[smurf_60.59$Site_ID == "60-0.59B",]
ab_smurfsamples <- merge(a_smurfsamples, b_smurfsamples, by = c("Date_Out", "Per_cov", "Frag"), all = F)

# Histograms
par(mfrow=c(2,1))
hist(ab_smurfsamples$N.x) # "60-0.59A"
hist(ab_smurfsamples$N.y) # "60-0.59B"
par(mfrow=c(1,1))

#    Are SMURF settlement rates  in Sites 69%-059 A and B on the two shoal environments different
wilcox.test(ab_smurfsamples$N.x, ab_smurfsamples$N.y, paired = T)
# Not significantly different

# Because the settlement rates are not different between Sites 69%-059 A and B, 
# Pool these sites (take mean) to have 25 unique site designs for the analyses

# Are both 69%-059 A and B sampled on all the same dates?
unique(smurf.settlers$Date_Out[smurf.settlers$Site_ID == "60-0.59A"]) #dates for site A
unique(smurf.settlers$Date_Out[smurf.settlers$Site_ID == "60-0.59B"]) #dates for Site B
# Site B was not sampled on one date = "6/8/2018"

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

#### Settlement rates of larval fish and crab megalopae (SMURFs) ####
par(mfrow=c(1,1))
hist(smurf.set.sum$N)

# GLMs for settlement rate across landscape %cover and frag per se
set.glm <- glm(N ~ Per_cov * Frag, data = smurf.set.sum, family = "poisson") #pois model
set.glmnb <- glm.nb(N ~ Per_cov * Frag, data = smurf.set.sum) #neg bin model
lrtest(set.glm, set.glmnb) # Likelihood ratio test to assess fit of models -> nb is better model

plot(simulateResiduals(fittedModel = set.glmnb, n = 250), asFactor = F)

# Though not significant - there may be an outlier - the red points - let's find them
outliers(simulateResiduals(fittedModel = set.glmnb, n = 250), return = "index") # indexes of outliers
plot(set.glmnb, which = 4)
smurf.set.sum[107,] # this one is clearly the most extreme - cook's distance >0.5
smurf.set.sum[smurf.set.sum$N == 49,]# obs 107 needs to be excluded - extreme high catch
smurf.set[smurf.set$Date_Out == "2018-07-11" & smurf.set$Site_ID == "60-0.225",] # 48 megalopa + 1 fish

# This seems like an extreme outlier 98% composes of Brachyuran megalopa - remove it and compare fit and results
new.smurf.set.sum <- smurf.set.sum[!smurf.set.sum$N == 49,] # remove smurf which had 48 megalopa + 1 fish
hist(smurf.set.sum$N)
hist(new.smurf.set.sum$N)

# GLMs (EXCLUDING OUTLIER) for settlement rate across landscape %cover and frag per se
new.set.glm <- glm(N ~ Per_cov * Frag, data = new.smurf.set.sum, family = "poisson") #pois model
new.set.glmnb <- glm.nb(N ~ Per_cov * Frag, data = new.smurf.set.sum) #neg bin model
lrtest(new.set.glm, new.set.glmnb) # Likelihood ratio test to assess fit of models -> nb is still better model

plot(simulateResiduals(fittedModel = new.set.glmnb, n = 250), asFactor = F)# good fit
# still pointing out one outlier but considering the histograms and the counts - I don't think they should be removed

# FINAL MODEL
summary(new.set.glmnb) # NS

# Model predictions
new.smurf <- with(new.smurf.set.sum, data.frame(Per_cov = Per_cov, Frag = Frag, N = N))
new.smurf$pred.p <- predict(update(new.set.glmnb, .~. -Frag-Per_cov:Frag), type = "response")
new.smurf$pred.f <- predict(update(new.set.glmnb, .~. -Per_cov-Per_cov:Frag), type = "response")

plot_model(new.set.glmnb, type = "pred", show.ci = TRUE, terms = "Per_cov")
plot_model(new.set.glmnb, type = "pred", show.ci = TRUE, terms = "Frag")

# GRAPH DATA without outlier
new.set.site.mu <- new.smurf.set.sum %>% group_by(Per_cov, Frag) %>% summarise(mean.N = mean(N),se = sd(N)/sqrt(length(N[!is.na(N)])))
new.set.pc <- new.set.site.mu %>% group_by(Per_cov) %>% summarise(mean.siteN = mean(mean.N), se = sd(mean.N)/sqrt(length(mean.N[!is.na(mean.N)])))
new.set.f <- new.set.site.mu %>% group_by(Frag) %>% summarise(mean.siteN = mean(mean.N), se = sd(mean.N)/sqrt(length(mean.N[!is.na(mean.N)])))

# In these plots - lines are predicted by the model
summary(new.set.glmnb)
smurf.summmary <- cbind("mod" = rep("smurf", length(new.set.glmnb)), coef(summary(new.set.glmnb)), "exp.est" = exp(coef(new.set.glmnb)), exp(confint(new.set.glmnb)))

ggpred.p <- ggpredict(new.set.glmnb, "Per_cov")
ggpred.f <- ggpredict(new.set.glmnb, "Frag")

# Plot results for settlement rates of larval fish and megalopa  vs Per_cov*Frag 
new.set.pc.scatplot <- ggplot() + ylim(0,6) + #ggpred(ggpred.p) + 
  geom_errorbar(data = new.set.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= mean.siteN-se, ymax = mean.siteN+se),  col = "black", width = 0, size = 1)+
  geom_point(data = new.set.pc, aes(x= as.numeric(as.character(Per_cov)), y= mean.siteN), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Settlement rate", x = NULL) +  toptheme
new.set.f.scatplot <- ggplot() + ylim(0,6) + #ggpred(ggpred.f) + 
  geom_errorbar(data = new.set.f, aes(x= as.numeric(as.character(Frag)), ymin= mean.siteN-se, ymax = mean.siteN+se),  col = "black", width = 0, size = 1)+
  geom_point(data = new.set.f, aes(x = as.numeric(as.character(Frag)), y= mean.siteN), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  labs(y = "Settlement rate", x = NULL) +  sidetheme + coord_flip()
smurf.set.sum.NEW.plot <- ggplot(data = new.set.site.mu, aes(y=Frag, x=Per_cov*100, fill = mean.N)) +
  labs(fill = "Settlement\nrate")+ siteaes + sitetheme 
sm.set.l <- g_legend(smurf.set.sum.NEW.plot)
grid.arrange(new.set.pc.scatplot, sm.set.l, smurf.set.sum.NEW.plot+theme(legend.position='none'), new.set.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

#### Results for MS ####

#    Are SMURF settlement rates  in Sites 69%-059 A and B on the two shoal environments different
wilcox.test(ab_smurfsamples$N.x, ab_smurfsamples$N.y, paired = T)
# Not significantly different

# Settlement rates of larval fishes and crab megalopae (in SMURFs) across landscape parameters (Per_cov*Frag) 
smurf.summmary
grid.arrange(new.set.pc.scatplot, sm.set.l, smurf.set.sum.NEW.plot+theme(legend.position='none'), new.set.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))
