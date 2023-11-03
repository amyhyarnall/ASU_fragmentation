# Squidpop generalist consumption assays vs Post-Florence Landscape % Cover and Frag

# MS Title: Habitat area more consistently influences seagrass faunal communities than fragmentation per se
# Authors: Yarnall AH, Yeager LA, Lopazanski C, Poray AK, Morley J, Fodrie FJ
# Journal: Ecological monographs

# All dataset related to this manuscript are publicly available at:
# https://www.bco-dmo.org/project/714026

# Load libraries
library(dplyr)
library(lubridate)
library(DHARMa)
library(ggplot2)

# Generalist consumption assays - Squidpops
squidpop <- read.csv('https://datadocs.bco-dmo.org/file/6Y0Y04NCo2ZyzP/asufrag_squidpopconsumptionprob.csv')
# Pre- and Post-Hurricane Florence landscape parameters
landscape <- read.csv('https://datadocs.bco-dmo.org/file/JEPYRqPuB2GxYo/asu_frag_landscape_parameters.csv')

#### Post-Hurricane Florence Landscape parameter exploration ####

# Separate out Pre- and Post- Flo data
# format pre-flo data
landscape.pre <- landscape[landscape$Landscape_timeperiod == "original design",] # pull out pre-flo data
landscape.pre <- landscape.pre[names(landscape.pre)%in% c("Per_cov", "Frag", "Actual_cover", "Site_ID", "N_patches")]

# format post-flo data
# only consider parameters calculated when buried ASUs are excluded
landscape.post <- landscape[landscape$Buried_ASUs == "Excluded",]# pull out post-flo data
landscape.post <- landscape.post[names(landscape.post)%in% c("Per_cov", "Frag", "Actual_cover", "Site_ID", "N_patches")]
names(landscape.post)[names(landscape.post)=="N_patches"] <- "new.N_patches" #rename col
names(landscape.post)[names(landscape.post)=="Actual_cover"] <- "new.Actual_cover" #rename col

landscape.prepost <- merge(landscape.pre, landscape.post, all = T) #remerge dataframes

# If the new %cover or n.patches is 0 or "NA" make this an NA because the landscape was completely removed
#     therefore the landscape was non-existent and not sampled post-flo
landscape.prepost$new.Actual_cover <- ifelse(landscape.prepost$new.Actual_cover == 0, NA, landscape.prepost$new.Actual_cover)
landscape.prepost$new.N_patches <- ifelse(landscape.prepost$new.N_patches == 0, NA, landscape.prepost$new.N_patches)

# The treatment level of %cover may not have *exactly* matched the true percent cover (rounding involved)
# make sure this difference does not get carried over as a change post-flo 
# i.e., if the true %covers match pre- and post- so do the treatment levels
landscape.prepost$Per_cov_upd <- ifelse(landscape.prepost$new.Actual_cover == landscape.prepost$Actual_cover, 
                                        landscape.prepost$Per_cov, landscape.prepost$new.Actual_cover)

# Show the relationship between Percolation probability (Frag treatment) and number of patches
frag.npatches_mod <- lm(Frag ~ as.numeric(N_patches), data = landscape.prepost)
summary(frag.npatches_mod)

frag.patch <- ggplot(landscape.prepost, aes(y = Frag, x = as.numeric(N_patches))) + 
  geom_jitter(size = 5, shape = 21,width = 0.1, height = 0.01) +
  geom_smooth(method = "lm", se = T, col = "black") + 
  #scale_y_continuous(breaks = c(0,0.1,0.225,0.35,0.475,0.59)) + 
  scale_x_continuous(breaks = seq(2,10,2), limits = c(1,10)) + 
  labs(y = "Percolation Probability", x = "Number of Patches") + theme_classic() +
  theme(legend.position= "bottom", 
                   legend.text=element_text(size=20), 
                   legend.title = element_text(size=20), 
                   axis.text.x = element_text(colour = "black", size=20), 
                   axis.text.y = element_text(colour = "black", size=20), 
                   axis.title.x = element_text(size=20), 
                   axis.title.y = element_text(size=20), 
                   strip.text.x = element_text(colour = "black", size = 20))


# create new dataframe with pre and post flo Site info and use it to predict Frag_upd
prepostflo.sites <- data.frame(with(landscape.prepost, cbind(Site_ID,old.N_patches=N_patches,Per_cov,Frag,N_patches=new.N_patches,Per_cov_upd)))
prepostflo.sites$Frag_upd <- predict(frag.npatches_mod, type = "response", newdata = prepostflo.sites)
names(prepostflo.sites)[names(prepostflo.sites)=="N_patches"] <- "new.N_patches" #rename col

# Same as with Per_cov: The treatment level of Frag may not have *exactly* matched the true PP (rounding involved)
# make sure this difference does not get carried over as a change post-flo 
# i.e., if the true PP match pre- and post- so do the treatment levels
prepostflo.sites$Frag_upd <- ifelse((prepostflo.sites$Per_cov == prepostflo.sites$Per_cov_upd &
                                       prepostflo.sites$old.N_patches == prepostflo.sites$new.N_patches),
                                    prepostflo.sites$Frag, prepostflo.sites$Frag_upd)
prepostflo.sites <- prepostflo.sites[,c("Site_ID","Per_cov","Frag","Per_cov_upd","Frag_upd","old.N_patches","new.N_patches")] #reorder cols in df


# Are Post-flo Per_cov and Frag correlated??
# Significant negative correlation
cor.test(as.numeric(prepostflo.sites$Per_cov_upd)/100, as.numeric(prepostflo.sites$Frag_upd), paired = T)

par(mfrow=c(2,2))
hist(as.numeric(prepostflo.sites$Per_cov), 
     main=NULL, 
     xlab= "Percent cover", 
     ylab = "Number of landscapes",
     breaks = seq(0,60,5),
     ylim = c(0,5))
title("Pre-Florence",adj=0)
hist(as.numeric(prepostflo.sites$Frag), 
     main=NULL, 
     xlab= "Fragmentation per se", 
     ylab = "Number of landscapes",
     breaks = seq(0.0,0.6,0.05),
     ylim = c(0,5))
hist(as.numeric(prepostflo.sites$Per_cov_upd), 
     main=NULL, 
     xlab="Percent cover", 
     ylab = "Number of landscapes",
     breaks = seq(0,60,5),
     ylim = c(0,10))
title("Post-Florence",adj=0)
hist(as.numeric(prepostflo.sites$Frag_upd), 
     main=NULL, 
     xlab= "Fragmentation per se", 
     ylab = NULL,
     breaks = seq(0.0,0.6,0.05),
     ylim = c(0,10))

#### Squidpop data formatting ####
squidpop$Date <- as.Date(squidpop$Date)
squidpop$month <- month(squidpop$Date, label = T, abbr = F)

# Remove erroneous columns 
squidpop <- squidpop[!names(squidpop) %in% c("Postflo_Per_cov","Postflo_Frag")]

# Find the mean +/- SD of number of squidpops deployed per landscape
mean(squidpop$N) #9.18
sd(squidpop$N) #1.75

# UPDATE THE SQUIDPOP DATA TO REFLECT THE TRUE LANDSCAPE METRICS (use the buried ASUs excluded parameters)
depredation <- squidpop[!squidpop$Site_ID == "60-0.59A",] # remove "60-0.59A" - post-flo parameters were not calculated
depredation$Site_ID <- gsub("60-0.59B", "60-0.59", depredation$Site_ID) # remove the B designation
depredation$Site_ID <- gsub("-0.10", "-0.1", depredation$Site_ID) # format site IDs to match
depredation$Site_ID <- gsub("22-", "22.5-", depredation$Site_ID) # format site IDs to match
depredation <- merge(depredation, prepostflo.sites[,names(prepostflo.sites) %in% c("Per_cov_upd", "Frag_upd", "Site_ID")], all = T)

# Find the time (rounded to nearest hour) elapsed between setting and checking squidpops
depredation$h_elapsed <- signif((as.difftime(depredation$Time_Check, format = "%H:%M") - as.difftime(depredation$Time_In, format = "%H:%M")), digits = 1)
unique(depredation$h_elapsed)
depredation$h_elapsed <- ifelse(depredation$h_elapsed < 1, 1, depredation$h_elapsed) # 0.8 and 0.9 are close to 1 h, round up
depredation$h_elapsed <- ifelse(depredation$h_elapsed == 4, 3, depredation$h_elapsed) # one 4 h data point - must be rounding error, only checked up to 3 h

# Remove sites that were not sampled  because they did not exist post-flo
depredation <- depredation[!is.na(depredation$Date),] 

# Is there are relationship between post-flo percent cover and the number of squidpops used in the site? 
# This could create bias in the analysis
ggplot(depredation) + geom_point(aes(x = as.numeric(Per_cov_upd), y = N)) + 
  labs(x = " Post-Florence Percent Cover", y = "Number of squidpops used") # most received 10 squidpops 
require(car)
summary(aov(N ~ as.numeric(Per_cov_upd), data = depredation))

# Choose which check hour that has the appropriate % consumption (~30%)
# Separate out the values of num_eaten into 3 check h cols
dep_rates <- depredation
dep_rates$N_eaten1 <- ifelse(dep_rates$h_elapsed == 1, dep_rates$N_eaten, 0)
dep_rates$N_eaten2 <- ifelse(dep_rates$h_elapsed == 2, dep_rates$N_eaten, 0)
dep_rates$N_eaten3 <- ifelse(dep_rates$h_elapsed == 3, dep_rates$N_eaten, 0)

# Reduce the date to be grouped by site and month
dep_rates_site <- dep_rates %>% group_by(month, Site_ID, Per_cov_upd, Frag_upd, Time_In, N) %>%
  summarise(N_eaten1 = sum(N_eaten1), N_eaten2 = sum(N_eaten2), N_eaten3 = sum(N_eaten3))
colSums(dep_rates_site[6:9])
sum(dep_rates_site$N_eaten1)/sum(dep_rates_site$N) #33% <- use this
(sum(dep_rates_site$N_eaten1)+sum(dep_rates_site$N_eaten2))/sum(dep_rates_site$N) #46%
(sum(dep_rates_site$N_eaten1)+sum(dep_rates_site$N_eaten2)+sum(dep_rates_site$N_eaten3))/sum(dep_rates_site$N) #55%

# COVERT to 1s and 0s for each squidpop
dep_rates_site1h <- dep_rates_site[,1:7] # take meta-data cols + N_eaten1
dep_rates_site1h$N_left <- dep_rates_site1h$N - dep_rates_site1h$N_eaten1 # N_eaten are successes, N_left are failures

# make N eaten and left into bernoulli trails of successes (k=1), and failures(k=0)
new.rep1 <- dep_rates_site1h[rep(seq_len(nrow(dep_rates_site1h[,1:6])),  dep_rates_site1h$N_eaten1), ]
new.rep1$k <- rep(1, nrow(new.rep1))
new.rep0 <- dep_rates_site1h[rep(seq_len(nrow(dep_rates_site1h[,1:6])),  dep_rates_site1h$N_left), ]
new.rep0$k <- rep(0, nrow(new.rep0))

depk1 <- merge(new.rep1[,c(1:5, 9)], new.rep0[,c(1:5, 9)], all = T) #0s = not eaten and 1s =  eaten!!
depk1$Per_cov_upd <-as.numeric(depk1$Per_cov_upd)
depk1$Frag_upd <-as.numeric(depk1$Frag_upd)

# For Per_cov and Frag comparability in model results and figures divide Per_cov by 100 to put in on a similar scale to Frag
depk1$Per_cov_upd100 <- depk1$Per_cov_upd
depk1$Per_cov_upd <- depk1$Per_cov_upd/100

#### Results for MS ####
#**************************************************************************************************
#1 h checks PERCENT COVER ONLY
depk1h_p.glm <- glm(k ~ Per_cov_upd, data = depk1, family = binomial(logit))
# DIAGNOSITICS
testUniformity(simulateResiduals(depk1h_p.glm)) #GOOD
testDispersion(simulateResiduals(depk1h_p.glm)) #GOOD

summary(depk1h_p.glm) # sig effect
exp(cbind(OR = coef(depk1h_p.glm), confint(depk1h_p.glm))) # odds ratios
# increase of percent cover by 1 unit (i.e. 1%) increases odds of consumption 

require(sjPlot)
plot_model(depk1h_p.glm, type = "pred", show.ci = TRUE, terms = "Per_cov_upd [all]")

require(ggeffects)
depk1_p.plot <- ggpredict(depk1h_p.glm, "Per_cov_upd [all]") %>% plot(ci = T, line.size = 1) +  
  geom_jitter(data = depk1, aes(y = k, x = Per_cov_upd),
              shape = 21, size = 3, fill = "black", alpha = 0.2, width = 0.01, height = 0.02) + 
  labs(y = "Probability of consumption at 1 h", x= "Percent cover") + theme_classic() + 
  scale_x_continuous(breaks = c(0,0.1,0.225,0.35,0.475,0.60), limits = c(0,0.60))+ 
  theme(legend.position= "bottom", plot.title = element_blank(),
                   legend.text=element_text(size=20), 
                   legend.title = element_text(size=20), 
                   axis.text.x = element_text(colour = "black", size=20), 
                   axis.text.y = element_text(colour = "black", size=20), 
                   axis.title.x = element_text(size=20), 
                   axis.title.y = element_text(size=20), 
                   strip.text.x = element_text(colour = "black", size = 20))
depk1_p.plot
summary(depk1h_p.glm)
exp(cbind(OR = coef(depk1h_p.glm), confint(depk1h_p.glm))) # odds ratios

#http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/
#py = 1/[1 + exp(-(Bo + B1*x))]

#**************************************************************************************************
# AS A CHECK OF CONFIDENCE IN THE RESULTS - 
# LET'S MAKE SURE THAT REMOVING LANDSCAPES WITH A LOWER SAMPLE SIZE OF SQUIDPOPS
# DOES NOT CHANGE THE PATTERN SEEN FOR DEPREDATION VS PERCENT COVER

# sites that got fewer than 10 squidpops
shorted_sites <- unique(depredation$Site_ID[depredation$N <10])

# 1 h checks PERCENT COVER ONLY - remove landscapes with fewer than 10 squidpop used
dep10k1h_p.glm <- glm(k ~ Per_cov_upd, data = depk1[!depk1$Site_ID %in% shorted_sites,], family = binomial(logit))
# DIAGNOSITICS
testUniformity(simulateResiduals(dep10k1h_p.glm)) #GOOD
testDispersion(simulateResiduals(dep10k1h_p.glm)) #GOOD

summary(dep10k1h_p.glm) # sig effect
exp(cbind(OR = coef(dep10k1h_p.glm), confint(dep10k1h_p.glm))) # odds ratios
# increase of percent cover by 1 unit (i.e. 1%) increases odds of consumption 

require(sjPlot)
plot_model(dep10k1h_p.glm, type = "pred", show.ci = TRUE, terms = "Per_cov_upd [all]")

require(ggeffects)
dep10k1_p.plot <- ggpredict(dep10k1h_p.glm, "Per_cov_upd [all]") %>% plot(ci = T, line.size = 1) +  
  geom_jitter(data = depk1, aes(y = k, x = Per_cov_upd),
              shape = 21, size = 3, fill = "black", alpha = 0.2, width = 0.01, height = 0.02) + 
  labs(y = "Probability of consumption at 1 h", x= "Percent cover") + theme_classic() + 
  scale_x_continuous(breaks = c(0,0.1,0.225,0.35,0.475,0.60), limits = c(0,0.60))+ 
  theme(legend.position= "bottom", plot.title = element_blank(),
        legend.text=element_text(size=20), 
        legend.title = element_text(size=20), 
        axis.text.x = element_text(colour = "black", size=20), 
        axis.text.y = element_text(colour = "black", size=20), 
        axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), 
        strip.text.x = element_text(colour = "black", size = 20))
dep10k1_p.plot
summary(dep10k1h_p.glm)
exp(cbind(OR = coef(dep10k1h_p.glm), confint(dep10k1h_p.glm))) # odds ratios

#http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/
#py = 1/[1 + exp(-(Bo + B1*x))]

# RESULTS ARE VERY SIMILAR REGARDLESS OF LOWER SAMPLE SIZE SITES INCLUDED OR EXCLUDED
