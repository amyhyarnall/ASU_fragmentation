# Benthopelagic fish CPUE (catch per unit effort) in DIDSON samples

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
library(lme4)
library(car)

# Benthopelagic fish CPUE (recorded in DIDSON samples)
didson <- read.csv('https://datadocs.bco-dmo.org/file/N7870JESMjKmpp/asufrag_didsonfishcpue.csv')

#### Preliminary Benthopelagic fish analyses ####

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
didson_SampleL <- didson %>% group_by(Class, Month, Site_ID, Per_cov, Frag, Reader, Frames_sampled) %>%
  summarise(muL_cm = mean(L_cm))

didson_Sample <- merge(didson_SampleN, didson_SampleL, all = T)

#create sample Id column
didson_Sample$ID <- paste(didson_Sample$Class, didson_Sample$Month, didson_Sample$Site_ID, sep = "_")
length(unique(didson_Sample$ID)) # 166 samples

## Is there any reader bias? 
# Compare samples with more than one reader
reader.bias <- didson_Sample[duplicated(didson_Sample[,c(1:5,7)]) | duplicated(didson_Sample[,c(1:5,7)], fromLast = T),]

# Count the number of reader for each of these samples
N_readers <- reader.bias %>% group_by(Class,Month,Site_ID,Per_cov,Frag,Frames_sampled) %>% summarize(N_readers = length(unique(Reader)))
reader.bias <- merge(reader.bias, N_readers, all = T)

# How many samples were read by 2 and 3 readers
length(N_readers$N_readers[N_readers$N_readers == 2]) #13 read by 2 readers
length(N_readers$N_readers[N_readers$N_readers == 3]) #15 read by 3 readers

# Create individual dataframes for each reader's samples
# Reader AY
AY_samples <- didson_Sample[didson_Sample$Reader == "AY",]
AY_samples <- AY_samples[,-6]
colnames(AY_samples)[7:9] <- paste(colnames(AY_samples)[7:9],"AY", sep = "_")
# Reader MM
MM_samples <- didson_Sample[didson_Sample$Reader == "MM",]
MM_samples <- MM_samples[,-6]
colnames(MM_samples)[7:9] <- paste(colnames(MM_samples)[7:9],"MM", sep = "_")
# Reader MC
MC_samples <- didson_Sample[didson_Sample$Reader == "MC",]
MC_samples <- MC_samples[,-6]
colnames(MC_samples)[7:9] <- paste(colnames(MC_samples)[7:9],"MC", sep = "_")

# Merge the samples back together with the new tags
reader.bias_test <- merge(merge(merge(reader.bias[,1:6], AY_samples, all = T), MM_samples, all = T), MC_samples, all = T)
# remove samples with only one reader
reader.bias_test <- reader.bias_test[reader.bias_test$ID %in% reader.bias_test$ID,]

par(mfrow=c(3,1))
hist(reader.bias_test$SampleN_AY)
hist(reader.bias_test$SampleN_MM)
hist(reader.bias_test$SampleN_MC)

# Wilcoxon Signed-Rank Tests to compare each reader fish counts to the other pairwise 
wilcox.test(reader.bias_test$SampleN_AY, reader.bias_test$SampleN_MM, paired = T) # NS between AY and MM fish counts
wilcox.test(reader.bias_test$SampleN_AY, reader.bias_test$SampleN_MC, paired = T) # NS between AY and MC fish counts
wilcox.test(reader.bias_test$SampleN_MM, reader.bias_test$SampleN_MC, paired = T) # NS between MM and MC fish counts

par(mfrow=c(3,1))
hist(reader.bias_test$muL_cm_AY)
hist(reader.bias_test$muL_cm_MM)
hist(reader.bias_test$muL_cm_MC)
par(mfrow=c(1,1))

# Wilcoxon Signed-Rank Tests to compare each reader fish counts to the other pairwise 
wilcox.test(reader.bias_test$muL_cm_AY, reader.bias_test$muL_cm_MM, paired = T) # SIG DIFF between AY and MM fish lengths
wilcox.test(reader.bias_test$muL_cm_AY, reader.bias_test$muL_cm_MC, paired = T) # SIG DIFF between AY and MC fish lengths
wilcox.test(reader.bias_test$muL_cm_MM, reader.bias_test$muL_cm_MC, paired = T) # SIG DIFF between MM and MC fish lengths

# Could this be due to the cases where the are disagreements on fish counts? 

# Is there a relationship between fish count and variation of the fish count - difficult to accurately count more fish
reader.bias_fishcountvar <- reader.bias %>% group_by(ID, Class, Month, Site_ID, Per_cov, Frag, Frames_sampled) %>%
  summarise(mu_N = mean(SampleN), sd_N = sd(SampleN), mu_L = mean(muL_cm), sd_L = sd(muL_cm), readers = length(unique(Reader)))

hist(reader.bias_fishcountvar$mu_N)
summary(lm(mu_N ~ sd_N, data = )) # NS
ggplot(reader.bias_fishcountvar) + geom_point(aes(x=mu_N, y=sd_N)) # pretty scattered

# No relationship: Subset to only include frames where fish counts agreed
# Only compare frames where we counted the same number of fish
reader.bias_LFrame <- didson %>% group_by(Month, Site_ID, Class, Reader, FrameSample) %>% summarise(FrameN = sum(N)) #total fish per frame per reader
reader.bias_LFrame <- reader.bias_LFrame[reader.bias_LFrame$FrameN > 0,] # remove 0 fish frames
reader.bias_LFrame <- reader.bias_LFrame[duplicated(reader.bias_LFrame[,-4]) | duplicated(reader.bias_LFrame[,-4], fromLast = T),] # find frames on which more than 1 person got the same count
reader.bias_LFrame <- reader.bias_LFrame %>% group_by(Month, Site_ID, Class, FrameSample) %>% filter(length(unique(Reader))==3) #find the frames where all saw the same number of fish

framelist.allreaders <- unique(reader.bias_LFrame[,c("Month", "Site_ID", "Class","FrameSample")])
all.saw <- merge(framelist.allreaders, didson, all = F)
par(mfrow=c(3,1))
hist(all.saw$L_cm[all.saw$Reader == "AY"],breaks=seq(0,50,10), ylim=c(0,30))
hist(all.saw$L_cm[all.saw$Reader == "MM"],breaks=seq(0,50,10), ylim=c(0,30))
hist(all.saw$L_cm[all.saw$Reader == "MC"],breaks=seq(0,50,10), ylim=c(0,30)) # smaller than the other two

par(mfrow=c(1,1))
hist(all.saw$L_cm, breaks=seq(0,50,10))
length(unique(all.saw$FrameSample))
L_bias.mod <- Anova(lmer(L_cm ~ Reader + (1|FrameSample), data = all.saw)); L_bias.mod
# Fish Length Z-scores were calculated to account for significant differences in reader measurements!

# Compare DIDSON CPUE in Sites 69%-059 A and B on the two shoal environments
# Subset data for the comparison
didson_60.59 <- didson_Sample[!didson_Sample$Month == "October"& didson_Sample$Site_ID %in% c("60-0.59A", "60-0.59B"),]

# If samples were read by more than one reader - take average
didson_60.59sum <- didson_60.59 %>% group_by(Class,Month,Site_ID,Per_cov,Frag)%>%
  summarise(sum_Frames_sampled=sum(Frames_sampled),sum_N=sum(SampleN))

# Divide sum_N by (sum_Frames_sampled/10) - get average number of fish seen per 10 frames - round to whole number
didson_60.59sum$N <- round(didson_60.59sum$sum_N/(didson_60.59sum$sum_Frames_sampled/10))

# Delineate and merge paired data
a_didsonsamples <- didson_60.59sum[didson_60.59sum$Site_ID == "60-0.59A",]
b_didsonsamples <- didson_60.59sum[didson_60.59sum$Site_ID == "60-0.59B",]
ab_didsonsamples <- merge(a_didsonsamples, b_didsonsamples, by = c("Month", "Per_cov", "Frag", "Class"), all = F)

# Histograms
par(mfrow=c(2,1))
hist(ab_didsonsamples$N.x) # "60-0.59A"
hist(ab_didsonsamples$N.y) # "60-0.59B"

# Are DIDSON CPUEs in Sites 69%-059 A and B on the two shoal environments different
wilcox.test(ab_didsonsamples$N.x, ab_didsonsamples$N.y, paired = T)
# Not significantly different

# Format didson fish count data for analyses
didson_sub <- didson_Sample 

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

# Subset the dataframe for each group
lp_didsonsamples <- didson_data[didson_data$Class == "largest patch",]
np_didsonsamples <- didson_data[didson_data$Class == "near-patch",]
ip_didsonsamples <- didson_data[didson_data$Class == "inter-patch",]

# Merge paired didson samples
lp.np_didsonsamples <- merge(lp_didsonsamples, np_didsonsamples, by = c("Month", "Site_ID", "Per_cov", "Frag"), all = F)
lp.ip_didsonsamples <- merge(lp_didsonsamples, ip_didsonsamples, by = c("Month", "Site_ID", "Per_cov", "Frag"), all = F)

# Histograms
par(mfrow=c(2,1))
# CPUE asu vs near-patch
hist(lp.np_didsonsamples$N.x) #asu
hist(lp.np_didsonsamples$N.y) #mtrx
# CPUE asu vs inter-patch
hist(lp.ip_didsonsamples$N.x) #asu
hist(lp.ip_didsonsamples$N.y) #mtrx
par(mfrow=c(1,1))

#   What sampling location habitat type had greater didson CPUE: ASU or near-patch Matrix?
wilcox.test(lp.np_didsonsamples$N.x, lp.np_didsonsamples$N.y, paired = T) # NS

#   What sampling location habitat type had greater didson CPUE: ASU or inter-patch Matrix?
wilcox.test(lp.ip_didsonsamples$N.x, lp.ip_didsonsamples$N.y, paired = T)# NS

# Fish Length Z-scores were calculated to account for significant differences in reader measurements!
# Read in and format Fish Length Z-score data
didson_Lz <- read.csv('https://datadocs.bco-dmo.org/file/A8m8MoXhVk7Gxz/asufrag_didson_fish_zscore.csv')
didson_Lz$Per_cov100 <- didson_Lz$Per_cov
didson_Lz$Per_cov <- didson_Lz$Per_cov100/100

didson_Lz_data <- didson_Lz %>% group_by(Class, Month, Site_ID, Per_cov, Frag) %>%
  summarise(mu_z = mean(L_z, na.rm = T))

lp_didson_Lz_data <- didson_Lz_data[didson_Lz_data$Class == "largest patch",]
np_didson_Lz_data <- didson_Lz_data[didson_Lz_data$Class == "near-patch",]
ip_didson_Lz_data <- didson_Lz_data[didson_Lz_data$Class == "inter-patch",]

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



#### LARGEST PATCH Benthopelagic fish CPUE and Length z-score (DIDSON) ####

# Largest patch FISH COUNTS 
didsonlp.glm <- glm(N ~ Per_cov * Frag, data = lp_didsonsamples, family = "poisson")
didsonlp.glmnb <- glm.nb(N ~ Per_cov * Frag, data = lp_didsonsamples)
lrtest(didsonlp.glm, didsonlp.glmnb) # nb better

plot(simulateResiduals(fittedModel = didsonlp.glmnb, n = 250), asFactor = F)# ok fit - but nearly sig dispersion
testDispersion(simulateResiduals(fittedModel = didsonlp.glmnb, n = 250)) # nearly sig dispersion
hist(lp_didsonsamples$N)
plot(didsonlp.glmnb, which = 4) # extreme outlier

lp_didsonsamples[17,] # this site had huge aggregation of (forage) fish - remove point
new.lp_didsonsamples <- lp_didsonsamples[-17,]

new.didsonlp.glm <- glm(N ~ Per_cov * Frag, data = new.lp_didsonsamples, family = "poisson")
new.didsonlp.glmnb <- glm.nb(N ~ Per_cov * Frag, data = new.lp_didsonsamples)
lrtest(new.didsonlp.glm, new.didsonlp.glmnb) # nb better

plot(simulateResiduals(fittedModel = new.didsonlp.glmnb, n = 250), asFactor = F)# good fit
summary(new.didsonlp.glmnb) #NS

# Largest patch FISH CPUE
new.fish.lp <- with(new.lp_didsonsamples, data.frame(Per_cov = Per_cov, Frag = Frag, N = N))
new.fish.lp$pred.p <- predict(update(new.didsonlp.glmnb, .~. -Frag-Per_cov:Frag), type = "response")
new.fish.lp$pred.f <- predict(update(new.didsonlp.glmnb, .~. -Per_cov-Per_cov:Frag), type = "response")

plot_model(new.didsonlp.glmnb, type = "pred", show.ci = TRUE, terms = "Per_cov")
plot_model(new.didsonlp.glmnb, type = "pred", show.ci = TRUE, terms = "Frag")

didson.lp.site <- new.lp_didsonsamples %>% group_by(Per_cov,Frag) %>% summarise(mu_N = mean(N), se = sd(N)/sqrt(length(N[!is.na(N)])))
didson.lp.pc <- didson.lp.site %>% group_by(Per_cov) %>% summarise(musite_N = mean(mu_N), se = sd(mu_N)/sqrt(length(mu_N[!is.na(mu_N)])))
didson.lp.f <- didson.lp.site %>% group_by(Frag) %>% summarise(musite_N = mean(mu_N), se = sd(mu_N)/sqrt(length(mu_N[!is.na(mu_N)])))

# In these plots - lines are predicted by the model
summary(new.didsonlp.glmnb)
didson.lp.summmary <- cbind("mod" = rep("didson lp", length(new.didsonlp.glmnb)), coef(summary(new.didsonlp.glmnb)), "exp.est" = exp(coef(new.didsonlp.glmnb)), exp(confint(new.didsonlp.glmnb)))

ggpred.didlp.p <- ggpredict(new.didsonlp.glmnb, "Per_cov")
ggpred.didlp.f <- ggpredict(new.didsonlp.glmnb, "Frag")

new.didsonlp.pc.scatplot <- ggplot() + ylim(0,15) + #ggpred(ggpred.didlp.p) +
  geom_errorbar(data = didson.lp.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= musite_N-se, ymax = musite_N+se), col = "black", width = 0, size = 1)+
  geom_point(data = didson.lp.pc, aes(x= as.numeric(as.character(Per_cov)), y= musite_N), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Fish CPUE", x = NULL) + theme_classic() + toptheme 
new.didsonlp.f.scatplot <- ggplot() + ylim(0,15) + #ggpred(ggpred.didlp.f) +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  geom_errorbar(data = didson.lp.f, aes(x= as.numeric(as.character(Frag)), ymin= musite_N-se, ymax = musite_N+se),  col = "black", width = 0, size = 1)+
  geom_point(data = didson.lp.f, aes(x = as.numeric(as.character(Frag)), y= musite_N), shape = 21, fill = "lightslategrey", size = 5)  +
  labs(y = "Fish CPUE", x = NULL) + theme_classic() + sidetheme + coord_flip()
didson.lp.site.plot <- ggplot(data = didson.lp.site, aes(x=Per_cov*100, y=Frag, fill = mu_N) ) +
  labs(fill = "Fish\nCPUE") + siteaes + sitetheme  
didson.lp.l <- g_legend(didson.lp.site.plot)
grid.arrange(new.didsonlp.pc.scatplot, didson.lp.l, didson.lp.site.plot+theme(legend.position='none'), new.didsonlp.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# Largest patch LENGTH Z-SCORES
par(mfrow=c(1,1))
hist(lp_didson_Lz_data$mu_z)

didsonlpz.lm <- lm(mu_z ~ Per_cov * Frag, data = lp_didson_Lz_data[!is.na(lp_didson_Lz_data$mu_z),])

plot(simulateResiduals(fittedModel = didsonlpz.lm, n = 250), asFactor = F)# good fit
testDispersion(simulateResiduals(fittedModel = didsonlpz.lm, n = 250))#good
summary(didsonlpz.lm)# NS

#### NEAR-PATCH Benthopelagic fish CPUE and Length z-score (DIDSON) ####

# Near-patch FISH COUNTS 
didsonnp.glm <- glm(N ~ Per_cov * Frag, data = np_didsonsamples, family = "poisson")
didsonnp.glmnb <- glm.nb(N ~ Per_cov * Frag, data = np_didsonsamples)
lrtest(didsonnp.glm, didsonnp.glmnb) # nb better

par(mfrow=c(1,1))
plot(simulateResiduals(fittedModel = didsonnp.glmnb, n = 250), asFactor = F)# bad fit
testDispersion(simulateResiduals(fittedModel = didsonnp.glmnb, n = 250)) 
hist(np_didsonsamples$N)
plot(didsonnp.glmnb, which = 4) # extreme outlier

np_didsonsamples[17,] # this site had huge aggregation of (forage) fish - remove point
new.np_didsonsamples <- np_didsonsamples[-17,]

new.didsonnp.glm <- glm(N ~ Per_cov * Frag, data = new.np_didsonsamples, family = "poisson")
new.didsonnp.glmnb <- glm.nb(N ~ Per_cov * Frag, data = new.np_didsonsamples)
lrtest(new.didsonnp.glm, new.didsonnp.glmnb) # nb better

plot(simulateResiduals(fittedModel = new.didsonnp.glmnb, n = 250), asFactor = F)# good fit
summary(new.didsonnp.glmnb) #NS

# Near-patch FISH CPUE
new.fish.np <- with(new.np_didsonsamples, data.frame(Per_cov = Per_cov, Frag = Frag, N = N))
new.fish.np$pred.p <- predict(update(new.didsonnp.glmnb, .~. -Frag-Per_cov:Frag), type = "response")
new.fish.np$pred.f <- predict(update(new.didsonnp.glmnb, .~. -Per_cov-Per_cov:Frag), type = "response")

plot_model(new.didsonnp.glmnb, type = "pred", show.ci = TRUE, terms = "Per_cov")
plot_model(new.didsonnp.glmnb, type = "pred", show.ci = TRUE, terms = "Frag")

didson.np.site <- new.np_didsonsamples %>% group_by(Per_cov,Frag) %>% summarise(mu_N = mean(N), se = sd(N)/sqrt(length(N[!is.na(N)])))
didson.np.pc <- didson.np.site %>% group_by(Per_cov) %>% summarise(musite_N = mean(mu_N), se = sd(mu_N)/sqrt(length(mu_N[!is.na(mu_N)])))
didson.np.f <- didson.np.site %>% group_by(Frag) %>% summarise(musite_N = mean(mu_N), se = sd(mu_N)/sqrt(length(mu_N[!is.na(mu_N)])))

# In these plots - lines are predicted by the model
summary(new.didsonnp.glmnb)
didson.np.summmary <- cbind("mod" = rep("didson np", length(new.didsonnp.glmnb)), coef(summary(new.didsonnp.glmnb)), "exp.est" = exp(coef(new.didsonnp.glmnb)), exp(confint(new.didsonnp.glmnb)))

ggpred.didnp.p <- ggpredict(new.didsonnp.glmnb, "Per_cov")
ggpred.didnp.f <- ggpredict(new.didsonnp.glmnb, "Frag")

new.didsonnp.pc.scatplot <- ggplot() + ylim(0,35) + #ggpred(ggpred.didnp.p) +
  geom_errorbar(data = didson.np.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= musite_N-se, ymax = musite_N+se), col = "black", width = 0, size = 1)+
  geom_point(data = didson.np.pc, aes(x= as.numeric(as.character(Per_cov)), y= musite_N), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Fish CPUE", x = NULL) + theme_classic() + toptheme 
new.didsonnp.f.scatplot <- ggplot() + ylim(0,35) + #ggpred(ggpred.didnp.f) +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  geom_errorbar(data = didson.np.f, aes(x= as.numeric(as.character(Frag)), ymin= musite_N-se, ymax = musite_N+se),  col = "black", width = 0, size = 1)+
  geom_point(data = didson.np.f, aes(x = as.numeric(as.character(Frag)), y= musite_N), shape = 21, fill = "lightslategrey", size = 5)  +
  labs(y = "Fish CPUE", x = NULL) + theme_classic() + sidetheme + coord_flip()
didson.np.site.plot <- ggplot(data = didson.np.site, aes(x=Per_cov*100, y=Frag, fill = mu_N) ) +
  labs(fill = "Fish\nCPUE") + siteaes + sitetheme  
didson.np.l <- g_legend(didson.np.site.plot)
grid.arrange(new.didsonnp.pc.scatplot, didson.np.l, didson.np.site.plot+theme(legend.position='none'), new.didsonnp.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# Near-patch LENGTH Z-SCORES
par(mfrow=c(1,1))
hist(np_didson_Lz_data$mu_z)

didsonnpz.lm <- lm(mu_z ~ Per_cov * Frag, data = np_didson_Lz_data[!is.na(np_didson_Lz_data$mu_z),])

plot(simulateResiduals(fittedModel = didsonnpz.lm, n = 250), asFactor = F)# good fit
testDispersion(simulateResiduals(fittedModel = didsonnpz.lm, n = 250))#good
summary(didsonnpz.lm)# NS

#### INTERPATCH Benthopelagic fish CPUE and Length z-score (DIDSON) ####

# Interpatch FISH COUNTS 
didsonip.glm <- glm(N ~ Per_cov * Frag, data = ip_didsonsamples, family = "poisson")
didsonip.glmnb <- glm.nb(N ~ Per_cov * Frag, data = ip_didsonsamples)
lrtest(didsonip.glm, didsonip.glmnb) # nb better

par(mfrow=c(1,1))
plot(simulateResiduals(fittedModel = didsonip.glmnb, n = 250), asFactor = F)# meh fit
testDispersion(simulateResiduals(fittedModel = didsonip.glmnb, n = 250)) 
hist(ip_didsonsamples$N)
plot(didsonip.glmnb, which = 4) # extreme outlier

ip_didsonsamples[14,] # this site had huge aggregation of (forage) fish - remove point
new.ip_didsonsamples <- ip_didsonsamples[-14,]

new.didsonip.glm <- glm(N ~ Per_cov * Frag, data = new.ip_didsonsamples, family = "poisson")
new.didsonip.glmnb <- glm.nb(N ~ Per_cov * Frag, data = new.ip_didsonsamples)
lrtest(new.didsonip.glm, new.didsonip.glmnb) # nb better

plot(simulateResiduals(fittedModel = new.didsonip.glmnb, n = 250), asFactor = F)# good fit
summary(new.didsonip.glmnb) #NS

# Interpatch FISH CPUE
new.fish.ip <- with(new.ip_didsonsamples, data.frame(Per_cov = Per_cov, Frag = Frag, N = N))
new.fish.ip$pred.p <- predict(update(new.didsonip.glmnb, .~. -Frag-Per_cov:Frag), type = "response")
new.fish.ip$pred.f <- predict(update(new.didsonip.glmnb, .~. -Per_cov-Per_cov:Frag), type = "response")

plot_model(new.didsonip.glmnb, type = "pred", show.ci = TRUE, terms = "Per_cov")
plot_model(new.didsonip.glmnb, type = "pred", show.ci = TRUE, terms = "Frag")

didson.ip.site <- new.ip_didsonsamples %>% group_by(Per_cov,Frag) %>% summarise(mu_N = mean(N), se = sd(N)/sqrt(length(N[!is.na(N)])))
didson.ip.pc <- didson.ip.site %>% group_by(Per_cov) %>% summarise(musite_N = mean(mu_N), se = sd(mu_N)/sqrt(length(mu_N[!is.na(mu_N)])))
didson.ip.f <- didson.ip.site %>% group_by(Frag) %>% summarise(musite_N = mean(mu_N), se = sd(mu_N)/sqrt(length(mu_N[!is.na(mu_N)])))

# In these plots - lines are predicted by the model
summary(new.didsonip.glmnb)
didson.ip.summmary <- cbind("mod" = rep("didson ip", length(new.didsonip.glmnb)), coef(summary(new.didsonip.glmnb)), "exp.est" = exp(coef(new.didsonip.glmnb)), exp(confint(new.didsonip.glmnb)))

ggpred.didip.p <- ggpredict(new.didsonip.glmnb, "Per_cov")
ggpred.didip.f <- ggpredict(new.didsonip.glmnb, "Frag")

new.didsonip.pc.scatplot <- ggplot() + ylim(0,25) + #ggpred(ggpred.didip.p) +
  geom_errorbar(data = didson.ip.pc, aes(x= as.numeric(as.character(Per_cov)), ymin= musite_N-se, ymax = musite_N+se), col = "black", width = 0, size = 1)+
  geom_point(data = didson.ip.pc, aes(x= as.numeric(as.character(Per_cov)), y= musite_N), shape = 21, fill = "lightslategrey", size = 5)  +
  scale_x_continuous(breaks = c(0.1, 0.225, 0.35, 0.475, 0.6), limits = c(0.06, 0.63))+
  labs(y = "Fish CPUE", x = NULL) + theme_classic() + toptheme 
new.didsonip.f.scatplot <- ggplot() + ylim(0,25) + #ggpred(ggpred.didip.f) +
  scale_x_reverse(breaks = c(0.59, 0.475, 0.35, 0.225, 0.1), limits = c(0.63,0.06))+
  geom_errorbar(data = didson.ip.f, aes(x= as.numeric(as.character(Frag)), ymin= musite_N-se, ymax = musite_N+se),  col = "black", width = 0, size = 1)+
  geom_point(data = didson.ip.f, aes(x = as.numeric(as.character(Frag)), y= musite_N), shape = 21, fill = "lightslategrey", size = 5)  +
  labs(y = "Fish CPUE", x = NULL) + theme_classic() + sidetheme + coord_flip()
didson.ip.site.plot <- ggplot(data = didson.ip.site, aes(x=Per_cov*100, y=Frag, fill = mu_N) ) +
  labs(fill = "Fish\nCPUE") + siteaes + sitetheme  
didson.ip.l <- g_legend(didson.ip.site.plot)
grid.arrange(new.didsonip.pc.scatplot, didson.ip.l, didson.ip.site.plot+theme(legend.position='none'), new.didsonip.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# Interpatch LENGTH Z-SCORES
par(mfrow=c(1,1))
hist(ip_didson_Lz_data$mu_z)

didsonipz.lm <- lm(mu_z ~ Per_cov * Frag, data = ip_didson_Lz_data[!is.na(ip_didson_Lz_data$mu_z),])

plot(simulateResiduals(fittedModel = didsonipz.lm, n = 250), asFactor = F)# good fit
testDispersion(simulateResiduals(fittedModel = didsonipz.lm, n = 250))#good
summary(didsonipz.lm)# NS

# Does minnow trap CPUE directly related to the interpatch distance (m)?
inter_dist <- data.frame(Site_ID = unique(new.ip_didsonsamples$Site_ID), 
                         interDist_m = c(3.6,8.6,9.96,10.8,3.44,4.3,2.1,6.44,1.72,2.55,1.48,6.29,3.44,1.2,3.6,1.72,1.2,1.2,2.58,2.58))
didson.ip_dist <- merge(new.ip_didsonsamples, inter_dist, by = "Site_ID", all = T)

didson.ip_dist.glmnb <- glm.nb(N ~ interDist_m, data = didson.ip_dist)
plot(simulateResiduals(fittedModel = didson.ip_dist.glmnb, n = 250), asFactor = F)#good fit 
summary(didson.ip_dist.glmnb) #NS

#### Results for MS ####

# Preliminary analyses
# Paired Wilcoxon rank sum tests

# Wilcoxon Signed-Rank Tests to compare each reader fish counts to the other pairwise 
wilcox.test(reader.bias_test$SampleN_AY, reader.bias_test$SampleN_MM, paired = T) # NS between AY and MM fish counts
wilcox.test(reader.bias_test$SampleN_AY, reader.bias_test$SampleN_MC, paired = T) # NS between AY and MC fish counts
wilcox.test(reader.bias_test$SampleN_MM, reader.bias_test$SampleN_MC, paired = T) # NS between MM and MC fish counts

# fish lengths
summary(lm(mu_N ~ sd_N, data = reader.bias_fishcountvar)) # NS
L_bias.mod <- Anova(lmer(L_cm ~ Reader + (1|FrameSample), data = all.saw)); L_bias.mod
# Fish Length Z-scores were calculated to account for significant differences in reader measurements!

# Are DIDSON CPUEs in Sites 69%-059 A and B on the two shoal environments different
wilcox.test(ab_didsonsamples$N.x, ab_didsonsamples$N.y, paired = T)
# Not significantly different

# What sampling location habitat type had greater didson CPUE: ASU or near-patch Matrix?
wilcox.test(asu.np_didsonsamples$N.x, asu.np_didsonsamples$N.y, paired = T) # NS

# What sampling location habitat type had greater didson CPUE: ASU or inter-patch Matrix?
wilcox.test(asu.ip_didsonsamples$N.x, asu.ip_didsonsamples$N.y, paired = T)# NS

# Benthopelagic fish CPUE across landscape parameters (Per_cov*Frag) 

# Largest ASU patch
didson.lp.summmary
grid.arrange(new.didsonlp.pc.scatplot, didson.lp.l, didson.lp.site.plot+theme(legend.position='none'), new.didsonlp.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# Near-patch matrix
didson.np.summmary
grid.arrange(new.didsonnp.pc.scatplot, didson.np.l, didson.np.site.plot+theme(legend.position='none'), new.didsonnp.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))

# Interpatch matrix
didson.ip.summmary
grid.arrange(new.didsonip.pc.scatplot, didson.ip.l, didson.ip.site.plot+theme(legend.position='none'), new.didsonip.f.scatplot, ncol = 2, nrow = 2, widths=c(4, 1), heights=c(1, 4))


# Benthopelagic fish length z-scores across landscape parameters (Per_cov*Frag) 

# Largest ASU patch
summary(didsonlpz.lm)# NS

# Near-patch matrix
summary(didsonnpz.lm)# NS

# Interpatch matrix
summary(didsonipz.lm)# NS
