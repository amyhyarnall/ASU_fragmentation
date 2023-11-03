# Do settlement rates of larval fish and crab meglopae in SMURFs 
# correspond to the CPUE of epibenthic fish and macroinvertebrate in Minnow traps across landscapes

# MS Title: Habitat area more consistently influences seagrass faunal communities than fragmentation per se
# Authors: Yarnall AH, Yeager LA, Lopazanski C, Poray AK, Morley J, Fodrie FJ
# Journal: Ecological monographs

# All dataset related to this manuscript are publicly available at:
# https://www.bco-dmo.org/project/714026

# Load libraries
library(dplyr)

# Epibenthic fish and invertebrate CPUE - caught in 'minnow' traps
minnow <- read.csv('https://datadocs.bco-dmo.org/file/WW8WwADCDoXgmx/asufrag_trapfaunalcpue.csv')
# Larval fish and megalopae settlement rates (caught in SMURFs)
smurf <- read.csv('https://datadocs.bco-dmo.org/file/EKPKZqQCgogj25/asufrag_smurf_settlement_rate.csv')

#### Initial Data formatting ####

## R-friendly data formatting 
minnow$Date_Out <- as.Date(minnow$Date_Out) # convert to date format
minnow$Sci_name <- ifelse(minnow$Sci_name == "", NA, minnow$Sci_name) # make empty strings into NAs

## R-friendly data formatting 
smurf$Date_Out <- as.Date(smurf$Date_Out) # convert to date format
smurf$Sci_name <- ifelse(smurf$Sci_name == "", NA, smurf$Sci_name) # make empty strings into NAs

# Select only just settled individuals and larvae
smurf.set <- smurf[smurf$Settler == "Y",]

# Standardize Site_ID for smurf and minnow 
smurf.set$Site_ID <- gsub("-0.10", "-0.1", smurf.set$Site_ID)

#### SMURF/MINNOW TRAP Data merging ####
unique(minnow$Sci_name)
unique(smurf.set$Sci_name)

# Find species that were caught in both gear types
minnow$Sci_name <- ifelse(minnow$Sci_name %in% c("Callinectes similis", "Callinectes Similis", "Callinectes sapidus"), "Callinectes sp.", minnow$Sci_name)
sp.in.common <- intersect(smurf.set$Sci_name, minnow$Sci_name);sp.in.common

# Pull out smurf sampled fish in common and format for combining
smurf.incom <- smurf.set[smurf.set$Sci_name %in% sp.in.common,]
smurf.incom <- smurf.incom[,c("Date_Out","Site_ID","Per_cov", "Frag", "Sci_name", "N")]
smurf.incom$Gear <- rep("SMURF", length(smurf.incom$Sci_name))

# Pull out minnow trap sampled fish in common and format for combining
minnow.incom <- minnow[minnow$Sci_name %in% sp.in.common,]
minnow.incom <- minnow.incom %>% group_by(Date_Out,Site_ID,Per_cov,Frag,Sci_name) %>% summarise(N = length(Sci_name))
minnow.incom$Gear <- rep("Minnow Trap", length(minnow.incom$Sci_name))

smurf.minnow <- rbind(smurf.incom, data.frame(minnow.incom))
smurf.minnow$month <- as.numeric(format(as.Date(smurf.minnow$Date_Out, format = "%m/%d/%Y"), "%m"))

# Create site list
smurf.minnow.sites <- unique(smurf[,c("Per_cov", "Frag")])
smurf.minnow.sites <- smurf.minnow.sites[rep(seq_len(nrow(smurf.minnow.sites)), each = 5), ]
smurf.minnow.sites$month <- rep(6:10, 25)

#### SMURF/MINNOW TRAP CPUE comparisons ####
# what months of smurfs and minnow traps can we compare?
unique(smurf.minnow$month[smurf.minnow$Gear == "SMURF"]) #June, July, August
unique(smurf.minnow$month[smurf.minnow$Gear == "Minnow Trap"]) #June, July, August, September, October

# What was the smurf and trap effort in each site during these months
smurf.set$month <- as.numeric(format(as.Date(smurf.set$Date_Out, format = "%m/%d/%Y"), "%m"))
smurf.effort <- smurf.set %>% group_by(month, Site_ID, Per_cov, Frag) %>% summarise(smurf.effort = length(unique(Date_Out)))

minnow$month <- as.numeric(format(as.Date(minnow$Date_Out, format = "%m/%d/%Y"), "%m"))
minnow.loc.effort <- minnow[minnow$Trap_class %in% c('largest patch','near-patch','inter-patch'),] %>% group_by(month, Site_ID, Per_cov, Frag, Trap_class) %>% summarise(minnow = length(unique(Date_Out)))
minnow.effort <- minnow.loc.effort %>% group_by(month, Site_ID, Per_cov, Frag) %>% summarise(minnow.effort = sum(minnow))

gear.effort <- merge(smurf.effort, minnow.effort, all = T)
table(gear.effort$month, gear.effort$Site_ID)

# How many of each species seen in each gear?
gear.mo <- smurf.minnow %>% group_by(month, Gear, Site_ID, Per_cov, Frag, Sci_name) %>% summarise(n.tot = sum(N))
gearN.effort <- merge(gear.mo, gear.effort, all = T)
gearN.effort$n.eff <- ifelse(gearN.effort$Gear == 'Minnow Trap', gearN.effort$n.tot/gearN.effort$minnow.effort,
                             gearN.effort$n.tot/gearN.effort$smurf.effort)

N.effort <- gearN.effort[,!names(gearN.effort) %in% c('n.tot', 'smurf.effort', 'minnow.effort')]
N.effort <- N.effort[!is.na(N.effort$Gear),]

# which species are caught most frequently and in high numbers?
table(N.effort$Sci_name, N.effort$Gear) # number of times each spp was caught by the gear
N.effort.tbl <- N.effort %>% group_by(Gear, Sci_name) %>% summarise(N = mean(n.eff)) # mean CPUE
N.effort.tbl[order(N.effort.tbl$Sci_name),]

#### SMURF/MINNOW TRAP - DIRECT SPP COMPARISON w/ 2 month timelag ####

#****************************************************************
# PLANEHEAD FILEFISH (high catch in both gears)
sh.lag <- N.effort[N.effort$Sci_name == "Stephanolepis hispidus",]
sh.tbl <- sh.lag %>% group_by(Gear, month) %>% summarise(n = sum(n.eff))

# separate out the catch of each gear type
sh.lag.sm <- sh.lag[sh.lag$Gear == "SMURF",]
sh.lag.sm$smurf.n <- sh.lag.sm$n.eff

sh.lag.minnow <- sh.lag[sh.lag$Gear == "Minnow Trap",]
sh.lag.minnow$minnow.n <- sh.lag.minnow$n.eff

#     Let's try a two month lag
sh.lag2.sm <- sh.lag.sm
sh.lag2.sm$month <- sh.lag.sm$month +2 #lag of 2 month: match june smurf w/ aug Minnow Traps
sh.lag2mo <- merge(sh.lag2.sm[,-c(5:7)], sh.lag.minnow[,-c(5:7)], all = T)
sh.lag2mo <- merge(sh.lag2mo, smurf.minnow.sites, all = T)
sh.lag2mo[is.na(sh.lag2mo)] <- 0
sh.lag2mo <- sh.lag2mo[sh.lag2mo$month %in% c(8, 9, 10),]
sh2.tbl <- sh.lag2mo %>% group_by(as.factor(month)) %>% summarise(sm = sum(smurf.n), minnow = sum(minnow.n))

#   CORRELATION OF SMURFS TO TRAPS WITH 2 MONTH LAG
sh.2molag.cor <- cor.test(sh.lag2mo$smurf.n, sh.lag2mo$minnow.n, paired = T);sh.2molag.cor # POSITIVE CORRELATION!
#****************************************************************

#****************************************************************
# SNAPPERS POOLED TOGETHER (SOME catch in both gears)
lutj.lag <- N.effort[grepl("Lutjanus", N.effort$Sci_name),]
table(lutj.lag$month, lutj.lag$Gear)

# separate out the catch of each gear type
lutj.lag.sm <- lutj.lag[lutj.lag$Gear == "SMURF",]
lutj.lag.sm$smurf.n <- lutj.lag.sm$n.eff

lutj.lag.minnow <- lutj.lag[lutj.lag$Gear == "Minnow Trap",]
lutj.lag.minnow$minnow.n <- lutj.lag.minnow$n.eff

#     Let's try a two month lag
lutj.lag2.sm <- lutj.lag.sm
lutj.lag2.sm$month <- lutj.lag.sm$month +2 #lag of 2 month: match june smurf w/ august Minnow Traps
lutj.lag2mo <- merge(lutj.lag2.sm[,-c(5:7)], lutj.lag.minnow[,-c(5:7)], all = T)
lutj.lag2mo <- merge(lutj.lag2mo, smurf.minnow.sites, all = T)
lutj.lag2mo[is.na(lutj.lag2mo)] <- 0
lutj.lag2mo <- lutj.lag2mo[lutj.lag2mo$month %in% c(8, 9, 10),]
lutj2.tbl <- lutj.lag2mo %>% group_by(as.factor(month)) %>% summarise(sm = sum(smurf.n), minnow = sum(minnow.n))

#   CORRELATION OF SMURFS TO TRAPS WITH 2 MONTH LAG
lutj.2molag.cor <- cor.test(lutj.lag2mo$smurf.n, lutj.lag2mo$minnow.n, paired = T);lutj.2molag.cor  # positive correlation
#****************************************************************

#****************************************************************
# Callinectes spp. (C. silimus and C. sapidus pooled because can't tell difference among megalope) (high catch in both gears)
cs.lag <- N.effort[N.effort$Sci_name == "Callinectes sp.",]
table(cs.lag$month, cs.lag$Gear)

# separate out the catch of each gear type
cs.lag.sm <- cs.lag[cs.lag$Gear == "SMURF",]
cs.lag.sm$smurf.n <- cs.lag.sm$n.eff

cs.lag.minnow <- cs.lag[cs.lag$Gear == "Minnow Trap",]
cs.lag.minnow$minnow.n <- cs.lag.minnow$n.eff

#     Let's try a two month lag
cs.lag2.sm <- cs.lag.sm
cs.lag2.sm$month <- cs.lag.sm$month +2 #lag of 2 month: match june smurf w/ july Minnow Traps
cs.lag2mo <- merge(cs.lag2.sm[,-c(5:7)], cs.lag.minnow[,-c(5:7)], all = T)
cs.lag2mo <- merge(cs.lag2mo, smurf.minnow.sites, all = T)
cs.lag2mo[is.na(cs.lag2mo)] <- 0
cs.lag2mo <- cs.lag2mo[cs.lag2mo$month %in% c(8, 9, 10),]

#   CORRELATION OF SMURFS TO TRAPS WITH 2 MONTH LAG
cs.2molag.cor <- cor.test(cs.lag2mo$smurf.n, cs.lag2mo$minnow.n, paired = T);cs.2molag.cor #No correlation
#****************************************************************


#### Results for MS ####
# Frequency of catch
table(N.effort$Sci_name, N.effort$Gear) 

# PLANEHEAD FILEFISH (Stephanolepis hispidus)
sh.2molag.cor

# SNAPPERS (Lutjanus spp.)
lutj.2molag.cor

# BLUE CRABS (Callinectes spp.)
cs.2molag.cor

