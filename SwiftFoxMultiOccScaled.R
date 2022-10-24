## load unmarked package
library(AICcmodavg)
library(unmarked)
library(LaplacesDemon)
library(MuMIn)
library(ggplot2)
library(gtable)
library(gridExtra)

## Import data
cdata = read.csv("SwiftFoxMultiOcc.csv")

## Organize site and sampling data
S = nrow(cdata) # number of sites, = 383
J = 28 # number of secondary sampling occasions
T = 3 # number of primary sampling occasions

## detection history
swiftfox = cdata[, c("SwiftFoxOcc1_2018", "SwiftFoxOcc2_2018", "SwiftFoxOcc3_2018", "SwiftFoxOcc4_2018", "SwiftFoxOcc5_2018", 
                     "SwiftFoxOcc6_2018", "SwiftFoxOcc7_2018", "SwiftFoxOcc8_2018", "SwiftFoxOcc9_2018", "SwiftFoxOcc10_2018",
                     "SwiftFoxOcc11_2018", "SwiftFoxOcc12_2018", "SwiftFoxOcc13_2018", "SwiftFoxOcc14_2018", "SwiftFoxOcc15_2018",
                     "SwiftFoxOcc16_2018", "SwiftFoxOcc17_2018", "SwiftFoxOcc18_2018", "SwiftFoxOcc19_2018", "SwiftFoxOcc20_2018",
                     "SwiftFoxOcc21_2018", "SwiftFoxOcc22_2018", "SwiftFoxOcc23_2018", "SwiftFoxOcc24_2018", "SwiftFoxOcc25_2018",
                     "SwiftFoxOcc26_2018", "SwiftFoxOcc27_2018", "SwiftFoxOcc28_2018",
                     "SwiftFoxOcc1_2019", "SwiftFoxOcc2_2019", 
                     "SwiftFoxOcc3_2019", "SwiftFoxOcc4_2019", "SwiftFoxOcc5_2019", 
                     "SwiftFoxOcc6_2019", "SwiftFoxOcc7_2019", "SwiftFoxOcc8_2019", "SwiftFoxOcc9_2019", "SwiftFoxOcc10_2019",
                     "SwiftFoxOcc11_2019", "SwiftFoxOcc12_2019", "SwiftFoxOcc13_2019", "SwiftFoxOcc14_2019", "SwiftFoxOcc15_2019",
                     "SwiftFoxOcc16_2019", "SwiftFoxOcc17_2019", "SwiftFoxOcc18_2019", "SwiftFoxOcc19_2019", "SwiftFoxOcc20_2019",
                     "SwiftFoxOcc21_2019", "SwiftFoxOcc22_2019", "SwiftFoxOcc23_2019", "SwiftFoxOcc24_2019", "SwiftFoxOcc25_2019",
                     "SwiftFoxOcc26_2019", "SwiftFoxOcc27_2019", "SwiftFoxOcc28_2019",
                     "SwiftFoxOcc1_2020", "SwiftFoxOcc2_2020", "SwiftFoxOcc3_2020", "SwiftFoxOcc4_2020", "SwiftFoxOcc5_2020", 
                     "SwiftFoxOcc6_2020", "SwiftFoxOcc7_2020", "SwiftFoxOcc8_2020", "SwiftFoxOcc9_2020", "SwiftFoxOcc10_2020",
                     "SwiftFoxOcc11_2020", "SwiftFoxOcc12_2020", "SwiftFoxOcc13_2020", "SwiftFoxOcc14_2020", "SwiftFoxOcc15_2020",
                     "SwiftFoxOcc16_2020", "SwiftFoxOcc17_2020", "SwiftFoxOcc18_2020", "SwiftFoxOcc19_2020", "SwiftFoxOcc20_2020",
                     "SwiftFoxOcc21_2020", "SwiftFoxOcc22_2020", "SwiftFoxOcc23_2020", "SwiftFoxOcc24_2020", "SwiftFoxOcc25_2020",
                     "SwiftFoxOcc26_2020", "SwiftFoxOcc27_2020", "SwiftFoxOcc28_2020")]

## detection covariates
scent = cdata[, c("Scent_1_2018", "Scent_2_2018", "Scent_3_2018", "Scent_4_2018", "Scent_5_2018", 
                  "Scent_6_2018", "Scent_7_2018", "Scent_8_2018", "Scent_9_2018", "Scent_10_2018",
                  "Scent_11_2018", "Scent_12_2018", "Scent_13_2018", "Scent_14_2018", "Scent_15_2018",
                  "Scent_16_2018", "Scent_17_2018", "Scent_18_2018", "Scent_19_2018", "Scent_20_2018",
                  "Scent_21_2018", "Scent_22_2018", "Scent_23_2018", "Scent_24_2018", "Scent_25_2018",
                  "Scent_26_2018", "Scent_27_2018", "Scent_28_2018",
                  "Scent_1_2019", "Scent_2_2019", "Scent_3_2019", "Scent_4_2019", "Scent_5_2019", 
                  "Scent_6_2019", "Scent_7_2019", "Scent_8_2019", "Scent_9_2019", "Scent_10_2019",
                  "Scent_11_2019", "Scent_12_2019", "Scent_13_2019", "Scent_14_2019", "Scent_15_2019",
                  "Scent_16_2019", "Scent_17_2019", "Scent_18_2019", "Scent_19_2019", "Scent_20_2019",
                  "Scent_21_2019", "Scent_22_2019", "Scent_23_2019", "Scent_24_2019", "Scent_25_2019",
                  "Scent_26_2019", "Scent_27_2019", "Scent_28_2019",
                  "Scent_1_2020", "Scent_2_2020", "Scent_3_2020", "Scent_4_2020", "Scent_5_2020", 
                  "Scent_6_2020", "Scent_7_2020", "Scent_8_2020", "Scent_9_2020", "Scent_10_2020",
                  "Scent_11_2020", "Scent_12_2020", "Scent_13_2020", "Scent_14_2020", "Scent_15_2020",
                  "Scent_16_2020", "Scent_17_2020", "Scent_18_2020", "Scent_19_2020", "Scent_20_2020",
                  "Scent_21_2020", "Scent_22_2020", "Scent_23_2020", "Scent_24_2020", "Scent_25_2020",
                  "Scent_26_2020", "Scent_27_2020", "Scent_28_2020")] ## Days since Scent_ applied


year = cdata[, c("Year.1", "Year.2", "Year.3")]

doy = cdata[,c("DOY_1_2018", "DOY_2_2018", "DOY_3_2018", "DOY_4_2018", "DOY_5_2018", 
               "DOY_6_2018", "DOY_7_2018", "DOY_8_2018", "DOY_9_2018", "DOY_10_2018",
               "DOY_11_2018", "DOY_12_2018", "DOY_13_2018", "DOY_14_2018", "DOY_15_2018",
               "DOY_16_2018", "DOY_17_2018", "DOY_18_2018", "DOY_19_2018", "DOY_20_2018",
               "DOY_21_2018", "DOY_22_2018", "DOY_23_2018", "DOY_24_2018", "DOY_25_2018",
               "DOY_26_2018", "DOY_27_2018", "DOY_28_2018",
               "DOY_1_2019", "DOY_2_2019", "DOY_3_2019", "DOY_4_2019", "DOY_5_2019", 
               "DOY_6_2019", "DOY_7_2019", "DOY_8_2019", "DOY_9_2019", "DOY_10_2019",
               "DOY_11_2019", "DOY_12_2019", "DOY_13_2019", "DOY_14_2019", "DOY_15_2019",
               "DOY_16_2019", "DOY_17_2019", "DOY_18_2019", "DOY_19_2019", "DOY_20_2019",
               "DOY_21_2019", "DOY_22_2019", "DOY_23_2019", "DOY_24_2019", "DOY_25_2019",
               "DOY_26_2019", "DOY_27_2019", "DOY_28_2019",
               "DOY_1_2020", "DOY_2_2020", "DOY_3_2020", "DOY_4_2020", "DOY_5_2020", 
               "DOY_6_2020", "DOY_7_2020", "DOY_8_2020", "DOY_9_2020", "DOY_10_2020",
               "DOY_11_2020", "DOY_12_2020", "DOY_13_2020", "DOY_14_2020", "DOY_15_2020",
               "DOY_16_2020", "DOY_17_2020", "DOY_18_2020", "DOY_19_2020", "DOY_20_2020",
               "DOY_21_2020", "DOY_22_2020", "DOY_23_2020", "DOY_24_2020", "DOY_25_2020",
               "DOY_26_2020", "DOY_27_2020", "DOY_28_2020")]

dlost = cdata[,c("DaysLost_1_2018", "DaysLost_2_2018", "DaysLost_3_2018", "DaysLost_4_2018", "DaysLost_5_2018", 
                 "DaysLost_6_2018", "DaysLost_7_2018", "DaysLost_8_2018", "DaysLost_9_2018", "DaysLost_10_2018",
                 "DaysLost_11_2018", "DaysLost_12_2018", "DaysLost_13_2018", "DaysLost_14_2018", "DaysLost_15_2018",
                 "DaysLost_16_2018", "DaysLost_17_2018", "DaysLost_18_2018", "DaysLost_19_2018", "DaysLost_20_2018",
                 "DaysLost_21_2018", "DaysLost_22_2018", "DaysLost_23_2018", "DaysLost_24_2018", "DaysLost_25_2018",
                 "DaysLost_26_2018", "DaysLost_27_2018", "DaysLost_28_2018",
                 
                 "DaysLost_1_2019", "DaysLost_2_2019", "DaysLost_3_2019", "DaysLost_4_2019", "DaysLost_5_2019", 
                 "DaysLost_6_2019", "DaysLost_7_2019", "DaysLost_8_2019", "DaysLost_9_2019", "DaysLost_10_2019",
                 "DaysLost_11_2019", "DaysLost_12_2019", "DaysLost_13_2019", "DaysLost_14_2019", "DaysLost_15_2019",
                 "DaysLost_16_2019", "DaysLost_17_2019", "DaysLost_18_2019", "DaysLost_19_2019", "DaysLost_20_2019",
                 "DaysLost_21_2019", "DaysLost_22_2019", "DaysLost_23_2019", "DaysLost_24_2019", "DaysLost_25_2019",
                 "DaysLost_26_2019", "DaysLost_27_2019", "DaysLost_28_2019",
                 "DaysLost_1_2020", "DaysLost_2_2020", "DaysLost_3_2020", "DaysLost_4_2020", "DaysLost_5_2020", 
                 "DaysLost_6_2020", "DaysLost_7_2020", "DaysLost_8_2020", "DaysLost_9_2020", "DaysLost_10_2020",
                 "DaysLost_11_2020", "DaysLost_12_2020", "DaysLost_13_2020", "DaysLost_14_2020", "DaysLost_15_2020",
                 "DaysLost_16_2020", "DaysLost_17_2020", "DaysLost_18_2020", "DaysLost_19_2020", "DaysLost_20_2020",
                 "DaysLost_21_2020", "DaysLost_22_2020", "DaysLost_23_2020", "DaysLost_24_2020", "DaysLost_25_2020",
                 "DaysLost_26_2020", "DaysLost_27_2020", "DaysLost_28_2020")]

## landscape covariates
shdi = cdata$SHDI900
sgp = cdata$SGPPrp150
loamy = cdata$Loamy50
crp = cdata$CRPPrp500
ag = cdata$RowcropPrp100

#scale
shdi <- scale(shdi)
sgp <- scale(sgp)
loamy <- scale(loamy)
crp <- scale(crp)
ag <- scale(ag)

## converts site covariates into data frame that can be read by unmarkedMultiFrame
##Landscape
sitecovs = data.frame(shdi, loamy, sgp, crp, ag)

obs.stuff = list(scent = scent, doy = doy, dlost=dlost)
sescovs = list(year = year)

## assemble in unmarkedMultFrame
cd = unmarkedMultFrame(y = swiftfox, siteCovs = sitecovs, yearlySiteCovs = sescovs, numPrimary=3, obsCovs = obs.stuff)

# Quick look at data - make sure everything looks alright
summary(cd)

sfcor = cor(sitecovs, use = "na.or.complete")
sfcor

##### Detection #####
## get estimates of occupancy and detection
occ = predict(s.02, type = "psi")
occ[1,]
occ

#Null model(ie, no covariates)
p.01 = colext(~1, ~1, ~1, ~1, cd)
# Models corrected for DETECTION
p.02 = colext(~1, ~1, ~1, ~scent,cd)
p.03 = colext(~1, ~1, ~1, ~doy, cd)
p.04 = colext(~1, ~1, ~1, ~scent+doy, cd)
p.05 = colext(~1, ~1, ~1, ~dlost, cd)
p.06 = colext(~1, ~1, ~1, ~scent+dlost,cd)
p.07 = colext(~1, ~1, ~1, ~doy+dlost, cd)
p.08 = colext(~1, ~1, ~1, ~scent+doy+dlost, cd)

# set up candidate model list
Cands = list(p.01, p.02, p.03, p.04, p.05, p.06, p.07, p.08)

# assign model names
Model.names = c("p.01 Intercept", "p.02 scent", "p.03 doy", "p.04 scent+doy", "p.05 Days lost",
                "p.06 scent+dlost", "p.07 doy+dlost", "p.08 scent+doy+dlost")

ptab = aictab(cand.set = Cands, modnames = Model.names) ### Select top model based on AICc
ptab

summary(p.08)

############ Model Avg (If you have multiple competitve models)########################################################
#library(MuMIn)
#d = model.avg(Cands, cumsum(weights) <.95)
#summary(d)


#capture.output(d, file = "swiftfox_det_modavg.txt")

#######################################################################


##### Occupancy ######
#Landscape
### Null Model ###
s.01 = colext(~1, ~1, ~1, ~doy+scent+dlost,cd)
##Covariates
s.02 = colext(~crp+loamy+sgp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.03 = colext(~ag+crp+loamy+sgp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.04 = colext(~crp+loamy+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.05 = colext(~ag+crp+loamy+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.06 = colext(~crp+loamy+sgp, ~1, ~1, ~doy+scent+dlost,cd)
s.07 = colext(~ag+loamy+sgp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.08 = colext(~ag+crp+loamy+sgp, ~1, ~1, ~doy+scent+dlost,cd)
s.09 = colext(~crp+loamy, ~1, ~1, ~doy+scent+dlost,cd)
s.10 = colext(~loamy+sgp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.11 = colext(~ag+crp+loamy, ~1, ~1, ~doy+scent+dlost,cd)
s.12 = colext(~crp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.13 = colext(~ag+crp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.14 = colext(~loamy+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.15 = colext(~ag+crp+sgp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.16 = colext(~crp+sgp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.17 = colext(~ag+loamy+sgp, ~1, ~1, ~doy+scent+dlost,cd)
s.18 = colext(~ag+loamy+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.19 = colext(~ag+sgp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.20 = colext(~loamy+sgp, ~1, ~1, ~doy+scent+dlost,cd)
s.21 = colext(~sgp+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.22 = colext(~ag+shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.23 = colext(~shdi+I(shdi^2), ~1, ~1, ~doy+scent+dlost,cd)
s.24 = colext(~ag+loamy, ~1, ~1, ~doy+scent+dlost,cd)
s.25 = colext(~ag+crp, ~1, ~1, ~doy+scent+dlost,cd)
s.26 = colext(~loamy, ~1, ~1, ~doy+scent+dlost,cd)
s.27 = colext(~ag+crp+sgp, ~1, ~1, ~doy+scent+dlost,cd)
s.28 = colext(~crp, ~1, ~1, ~doy+scent+dlost,cd)
s.29 = colext(~crp+sgp, ~1, ~1, ~doy+scent+dlost,cd)
s.30 = colext(~ag+sgp, ~1, ~1, ~doy+scent+dlost,cd)
s.31 = colext(~ag, ~1, ~1, ~doy+scent+dlost,cd)
s.32 = colext(~sgp, ~1, ~1, ~doy+scent+dlost,cd)
s.33 = colext(~crp+loamy+sgp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.34 = colext(~ag+crp+loamy+sgp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.35 = colext(~crp+loamy+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.36 = colext(~ag+crp+loamy+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.37 = colext(~loamy+sgp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.38 = colext(~ag+loamy+sgp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.39 = colext(~crp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.40 = colext(~loamy+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.41 = colext(~crp+sgp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.42 = colext(~ag+crp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.43 = colext(~ag+loamy+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.44 = colext(~ag+crp+sgp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.45 = colext(~ag+sgp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.46 = colext(~sgp+shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.47 = colext(~shdi, ~1, ~1, ~doy+scent+dlost,cd)
s.48 = colext(~ag+shdi, ~1, ~1, ~doy+scent+dlost,cd)

# set up candidate model list
Cands.s=list(s.01, s.02,	s.03,	s.04,	s.05,	s.06,	s.07,	s.08,	s.09,	
             s.10, s.11,	s.12,	s.13,	s.14,	s.15,	s.16,	s.17,	s.18,	
             s.19,	s.20,	s.21,	s.22,	s.23,	s.24,	s.25,	s.26,	s.27,	
             s.28,	s.29,	s.30,	s.31,	s.32, s.33, s.34, s.35, s.36,
             s.37, s.38,s.39, s.40, s.41, s.42, s.43, s.44, s.45,
             s.46,s.47, s.48)

# Cands.s = list(s.01, s.02, s.03, s.04, s.05, s.06, s.08, s.09, s.10, s.11, s.12, s.13, s.14, s.15, s.16, s.18, s.21, s.22) 
s.Model.names = c('s.01 Null','s.02 crp+loamy+sgp+shdi+I(shdi^2)',	's.03 ag+crp+loamy+sgp+shdi+I(shdi^2)',	
                  's.04 crp+loamy+shdi+I(shdi^2)',	's.05 ag+crp+loamy+shdi+I(shdi^2)',	's.06 crp+loamy+sgp',
                  's.07 ag+loamy+sgp+shdi+I(shdi^2)',	's.08 ag+crp+loamy+sgp',	's.09 crp+loamy',	
                  's.10 loamy+sgp+shdi+I(shdi^2)',	's.11 ag+crp+loamy',	's.12 crp+shdi+I(shdi^2)',	
                  's.13 ag+crp+shdi+I(shdi^2)',	's.14 loamy+shdi+I(shdi^2)',	's.15 ag+crp+sgp+shdi+I(shdi^2)',	
                  's.16 crp+sgp+shdi+I(shdi^2)',	's.17 ag+loamy+sgp',	's.18 ag+loamy+shdi+I(shdi^2)',	
                  's.19 ag+sgp+shdi+I(shdi^2)',	's.20 loamy+sgp',	's.21 sgp+shdi+I(shdi^2)',	
                  's.22 ag+shdi+I(shdi^2)',	's.23 shdi+I(shdi^2)',	's.24 ag+loamy',	's.25 ag+crp',	
                  's.26 loamy',	's.27 ag+crp+sgp',	's.28 crp',	's.29 crp+sgp',	's.30 ag+sgp',	's.31 ag',	
                  's.32 sgp', 's.33 crp+loamy+sgp+shdi',	's.34 ag+crp+loamy+sgp+shdi',	
                  's.35 crp+loamy+shdi',	's.36 ag+crp+loamy+shdi',	's.37 loamy+sgp+shdi',	
                  's.38 ag+loamy+sgp+shdi',	's.39 crp+shdi',	's.40 loamy+shdi',	's.41 crp+sgp+shdi',
                  's.42 ag+crp+shdi',	's.43 ag+loamy+shdi',	's.44 ag+crp+sgp+shdi',	's.45 ag+sgp+shdi',	
                  's.46 sgp+shdi',	's.47 shdi',	's.48 ag+shdi')



psitab = aictab(cand.set = Cands.s, second.ord = TRUE,  modnames = s.Model.names) ### Select top model based on AICc
psitab

summary(s.02)

###################################################
#### Colonization ### 
## Null ##
g.01 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~1, ~1, ~scent+doy+dlost,cd)

#SiteCovs crp, fallow, ag, grass, aindex, contag, rain, sgp)
## Corrected colonization models ##
g.02 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+sgp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.03 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+loamy+sgp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.04 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.05 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+loamy+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.06 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+sgp,~1,~scent+doy+dlost,cd)
g.07 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+loamy+sgp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.08 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+loamy+sgp,~1,~scent+doy+dlost,cd)
g.09 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy,~1,~scent+doy+dlost,cd)
g.10 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~loamy+sgp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.11 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+loamy,~1,~scent+doy+dlost,cd)
g.12 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.13 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.14 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~loamy+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.15 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+sgp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.16 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+sgp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.17 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+loamy+sgp,~1,~scent+doy+dlost,cd)
g.18 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+loamy+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.19 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+sgp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.20 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~loamy+sgp,~1,~scent+doy+dlost,cd)
g.21 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~sgp+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.22 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.23 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~shdi+I(shdi^2),~1,~scent+doy+dlost,cd)
g.24 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+loamy,~1,~scent+doy+dlost,cd)
g.25 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp,~1,~scent+doy+dlost,cd)
g.26 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~loamy,~1,~scent+doy+dlost,cd)
g.27 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+sgp,~1,~scent+doy+dlost,cd)
g.28 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp,~1,~scent+doy+dlost,cd)
g.29 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+sgp,~1,~scent+doy+dlost,cd)
g.30 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+sgp,~1,~scent+doy+dlost,cd)
g.31 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag,~1,~scent+doy+dlost,cd)
g.32 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~sgp,~1,~scent+doy+dlost,cd)
g.33 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+sgp+shdi,~1,~scent+doy+dlost,cd)
g.34 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+loamy+sgp+shdi,~1,~scent+doy+dlost,cd)
g.35 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~1,~scent+doy+dlost,cd)
g.36 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+loamy+shdi,~1,~scent+doy+dlost,cd)
g.37 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~loamy+sgp+shdi,~1,~scent+doy+dlost,cd)
g.38 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+loamy+sgp+shdi,~1,~scent+doy+dlost,cd)
g.39 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+shdi,~1,~scent+doy+dlost,cd)
g.40 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~loamy+shdi,~1,~scent+doy+dlost,cd)
g.41 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+sgp+shdi,~1,~scent+doy+dlost,cd)
g.42 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+shdi,~1,~scent+doy+dlost,cd)
g.43 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+loamy+shdi,~1,~scent+doy+dlost,cd)
g.44 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+crp+sgp+shdi,~1,~scent+doy+dlost,cd)
g.45 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+sgp+shdi,~1,~scent+doy+dlost,cd)
g.46 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~sgp+shdi,~1,~scent+doy+dlost,cd)
g.47 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~shdi,~1,~scent+doy+dlost,cd)
g.48 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~ag+shdi,~1,~scent+doy+dlost,cd)


# set up candidate model list
Cands.g=list(g.01, g.02,	g.03,	g.04,	g.05,	g.06,	g.07,	g.08,	g.09,	
             g.10, g.11,	g.12,	g.13,	g.14,	g.15,	g.16,	g.17,	g.18,	
             g.19,	g.20,	g.21,	g.22,	g.23,	g.24,	g.25,	g.26,	g.27,	
             g.28,	g.29,	g.30,	g.31,	g.32, g.33, g.34, g.35, g.36,
             g.37, g.38,g.39, g.40, g.41, g.42, g.43, g.44, g.45,
             g.46,g.47, g.48)


g.Model.names = c('g.01 Null','g.02 crp+loamy+sgp+shdi+I(shdi^2)',	'g.03 ag+crp+loamy+sgp+shdi+I(shdi^2)',	
                  'g.04 crp+loamy+shdi+I(shdi^2)',	'g.05 ag+crp+loamy+shdi+I(shdi^2)',	'g.06 crp+loamy+sgp',
                  'g.07 ag+loamy+sgp+shdi+I(shdi^2)',	'g.08 ag+crp+loamy+sgp',	'g.09 crp+loamy',	
                  'g.10 loamy+sgp+shdi+I(shdi^2)',	'g.11 ag+crp+loamy',	'g.12 crp+shdi+I(shdi^2)',	
                  'g.13 ag+crp+shdi+I(shdi^2)',	'g.14 loamy+shdi+I(shdi^2)',	'g.15 ag+crp+sgp+shdi+I(shdi^2)',	
                  'g.16 crp+sgp+shdi+I(shdi^2)',	'g.17 ag+loamy+sgp',	'g.18 ag+loamy+shdi+I(shdi^2)',	
                  'g.19 ag+sgp+shdi+I(shdi^2)',	'g.20 loamy+sgp',	'g.21 sgp+shdi+I(shdi^2)',	
                  'g.22 ag+shdi+I(shdi^2)',	'g.23 shdi+I(shdi^2)',	'g.24 ag+loamy',	'g.25 ag+crp',	
                  'g.26 loamy',	'g.27 ag+crp+sgp',	'g.28 crp',	'g.29 crp+sgp',	'g.30 ag+sgp',	'g.31 ag',	
                  'g.32 sgp', 'g.33 crp+loamy+sgp+shdi',	'g.34 ag+crp+loamy+sgp+shdi',	
                  'g.35 crp+loamy+shdi',	'g.36 ag+crp+loamy+shdi',	'g.37 loamy+sgp+shdi',	
                  'g.38 ag+loamy+sgp+shdi',	'g.39 crp+shdi',	'g.40 loamy+shdi',	'g.41 crp+sgp+shdi',
                  'g.42 ag+crp+shdi',	'g.43 ag+loamy+shdi',	'g.44 ag+crp+sgp+shdi',	'g.45 ag+sgp+shdi',	
                  'g.46 sgp+shdi',	'g.47 shdi',	'g.48 ag+shdi')

gamtab = aictab(cand.set = Cands.g, second.ord = TRUE, modnames = g.Model.names) ### Select top model based on AICc ## set second.ord = FALSE for AIC
gamtab

summary(g.35)


############################################################################################
#top occupancy model and colonization model is added
#### Extinction #####
##  Null  ##
e.01 = colext(~~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi, ~1, ~scent+doy,cd)
### Models corrected for extinction ###
e.02 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+loamy+sgp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.03 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+loamy+sgp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.04 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+loamy+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.05 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+loamy+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.06 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+loamy+sgp,~scent+doy+dlost,cd)
e.07 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+loamy+sgp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.08 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+loamy+sgp,~scent+doy+dlost,cd)
e.09 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+loamy,~scent+doy+dlost,cd)
e.10 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~loamy+sgp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.11 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+loamy,~scent+doy+dlost,cd)
e.12 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.13 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.14 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~loamy+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.15 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+sgp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.16 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+sgp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.17 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+loamy+sgp,~scent+doy+dlost,cd)
e.18 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+loamy+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.19 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+sgp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.20 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~loamy+sgp,~scent+doy+dlost,cd)
e.21 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~sgp+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.22 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+shdi+I(shdi^2),~scent+doy+dlost,cd)
e.23 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~shdi+I(shdi^2),~scent+doy+dlost,cd)
e.24 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+loamy,~scent+doy+dlost,cd)
e.25 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp,~scent+doy+dlost,cd)
e.26 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~loamy,~scent+doy+dlost,cd)
e.27 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+sgp,~scent+doy+dlost,cd)
e.28 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp,~scent+doy+dlost,cd)
e.29 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+sgp,~scent+doy+dlost,cd)
e.30 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+sgp,~scent+doy+dlost,cd)
e.31 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag,~scent+doy+dlost,cd)
e.32 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~sgp,~scent+doy+dlost,cd)
e.33 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+loamy+sgp+shdi,~scent+doy+dlost,cd)
e.34 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+loamy+sgp+shdi,~scent+doy+dlost,cd)
e.35 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+loamy+shdi,~scent+doy+dlost,cd)
e.36 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+loamy+shdi,~scent+doy+dlost,cd)
e.37 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~loamy+sgp+shdi,~scent+doy+dlost,cd)
e.38 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+loamy+sgp+shdi,~scent+doy+dlost,cd)
e.39 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+shdi,~scent+doy+dlost,cd)
e.40 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~loamy+shdi,~scent+doy+dlost,cd)
e.41 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~crp+sgp+shdi,~scent+doy+dlost,cd)
e.42 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+shdi,~scent+doy+dlost,cd)
e.43 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+loamy+shdi,~scent+doy+dlost,cd)
e.44 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+crp+sgp+shdi,~scent+doy+dlost,cd)
e.45 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+sgp+shdi,~scent+doy+dlost,cd)
e.46 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~sgp+shdi,~scent+doy+dlost,cd)
e.47 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~shdi,~scent+doy+dlost,cd)
e.48 = colext(~crp+loamy+sgp+shdi+I(shdi^2),~crp+loamy+shdi,~ag+shdi,~scent+doy+dlost,cd)


# set up candidate model list
Cands.e=list(e.01, e.02,	e.03,	e.04,	e.05,	e.06,	e.07,	e.08,	e.09,	
             e.10, e.11,	e.12,	e.13,	e.14,	e.15,	e.16,	e.17,	e.18,	
             e.19,	e.20,	e.21,	e.22,	e.23,	e.24,	e.25,	e.26,	e.27,	
             e.28,	e.29,	e.30,	e.31,	e.32, e.33, e.34, e.35, e.36,
             e.37, e.38,e.39, e.40, e.41, e.42, e.43, e.44, e.45,
             e.46,e.47, e.48)

# Cande.s = list(e.01, e.02, e.03, e.04, e.05, e.06, e.08, e.09, e.10, e.11, e.12, e.13, e.14, e.15, e.16, e.18, e.21, e.22) 
e.Model.names = c('e.01 Null','e.02 crp+loamy+sgp+shdi+I(shdi^2)',	'e.03 ag+crp+loamy+sgp+shdi+I(shdi^2)',	
                  'e.04 crp+loamy+shdi+I(shdi^2)',	'e.05 ag+crp+loamy+shdi+I(shdi^2)',	'e.06 crp+loamy+sgp',
                  'e.07 ag+loamy+sgp+shdi+I(shdi^2)',	'e.08 ag+crp+loamy+sgp',	'e.09 crp+loamy',	
                  'e.10 loamy+sgp+shdi+I(shdi^2)',	'e.11 ag+crp+loamy',	'e.12 crp+shdi+I(shdi^2)',	
                  'e.13 ag+crp+shdi+I(shdi^2)',	'e.14 loamy+shdi+I(shdi^2)',	'e.15 ag+crp+sgp+shdi+I(shdi^2)',	
                  'e.16 crp+sgp+shdi+I(shdi^2)',	'e.17 ag+loamy+sgp',	'e.18 ag+loamy+shdi+I(shdi^2)',	
                  'e.19 ag+sgp+shdi+I(shdi^2)',	'e.20 loamy+sgp',	'e.21 sgp+shdi+I(shdi^2)',	
                  'e.22 ag+shdi+I(shdi^2)',	'e.23 shdi+I(shdi^2)',	'e.24 ag+loamy',	'e.25 ag+crp',	
                  'e.26 loamy',	'e.27 ag+crp+sgp',	'e.28 crp',	'e.29 crp+sgp',	'e.30 ag+sgp',	'e.31 ag',	
                  'e.32 sgp', 'e.33 crp+loamy+sgp+shdi',	'e.34 ag+crp+loamy+sgp+shdi',	
                  'e.35 crp+loamy+shdi',	'e.36 ag+crp+loamy+shdi',	'e.37 loamy+sgp+shdi',	
                  'e.38 ag+loamy+sgp+shdi',	'e.39 crp+shdi',	'e.40 loamy+shdi',	'e.41 crp+sgp+shdi',
                  'e.42 ag+crp+shdi',	'e.43 ag+loamy+shdi',	'e.44 ag+crp+sgp+shdi',	'e.45 ag+sgp+shdi',	
                  'e.46 sgp+shdi',	'e.47 shdi',	'e.48 ag+shdi')

eptab = aictab(cand.set = Cands.e, second.ord = TRUE, modnames = e.Model.names) ### Select top model based on AICc
#eptab = aictab(cand.set = Cands.e, second.ord = FALSE,  modnames = e.Model.names) ### Select top model based on AIC
eptab
psitab
gamtab

### our extinction models were not significant and results were not reported
summary(e.21)

####### PLOTTING CODE #####################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# Load graph package
library(ggplot2)


## look at top models (can be done with single season models as well) 
##AICc Tables
psitab
gamtab
eptab

summary(s.02)
summary(g.35)
summary(e.21)


## check range of covariates that were in top models to be backtransformed
range(loamy, na.rm = TRUE)
range(sgp, na.rm = TRUE)
range(shdi, na.rm=T)
range(crp, na.rm=T)

#backtransform scaled data 
nd.loamy <- data.frame(loamy = seq(-0.7058545,1.7388209,length=381), sgp=0, shdi=0, crp=0) ## setting scaled sequence for proportion of urban landcover, may have to play with these to ##get the correct length on graph   
nd.sgp <- data.frame(sgp = seq(-0.8260378,3.2826471, length =381), loamy=0, shdi=0, crp=0) ## setting scaled sequence for distance to water 
nd.shdi <- data.frame(shdi = seq(-3.051014,2.314599,length=381), loamy=0, sgp=0, crp=0)
nd.crp <- data.frame(crp = seq(-0.6606292,3.9564714, length =381), loamy=0, shdi=0, sgp=0) ## setting scaled sequence for distance to water 


cloamy = predict(s.02, newdata = nd.loamy, type='psi')
csgp = predict(s.02, newdata = nd.sgp, type='psi')
cshdi = predict(s.02, newdata = nd.shdi, type='psi')
ccrp = predict(s.02, newdata = nd.crp, type='psi')

cloamy = predict(g.35, newdata = nd.loamy, type='col')
cshdi = predict(g.35, newdata = nd.shdi, type='col')
ccrp = predict(g.35, newdata = nd.crp, type='col')

#FOR MOD AVERAGE
#cloamy = modavgPred(Cands.g, g.Model.names, newdata = nd.loamy, second.ord = TRUE, nobs = NULL,
#                    uncond.se = "revised", conf.level = 0.95, type = "response", c.hat = 1,
#                   parm.type = 'gamma')


#colstream = predict(g.61, newdata = nd.stream, type='col')
#extrain = predict(e.18, newdata = nd.rain, type = "ext" )
#extgraze = predict(e.18, newdata = nd.graze, type = "ext" )

## put urban into right format (unscale and put into data frame)
new.loamy <- nd.loamy*attr(loamy,"scaled:scale")+attr(loamy,"scaled:center")
new.loamy1 <- new.loamy$loamy
new.loamy2 <- data.frame(new.loamy1)

new.sgp <- nd.sgp*attr(sgp,"scaled:scale")+attr(sgp,"scaled:center")
new.sgp1 <- new.sgp$sgp
new.sgp2 <- data.frame(new.sgp1)

new.shdi <- nd.shdi*attr(shdi,"scaled:scale")+attr(shdi,"scaled:center")
new.shdi1 <- new.shdi$shdi
new.shdi2 <- data.frame(new.shdi1)

new.crp <- nd.crp*attr(crp,"scaled:scale")+attr(crp,"scaled:center")
new.crp1 <- new.crp$crp
new.crp2 <- data.frame(new.crp1)


## put model averaged predictions into data frame
predict2 <- cloamy$Predicted
predict3 <- data.frame(predict2)

predict4 <- csgp$Predicted
predict5 <- data.frame(predict4)

predict6 <- cshdi$Predicted
predict7 <- data.frame(predict6)

predict8 <- ccrp$Predicted
predict9 <- data.frame(predict8)



## uper and lower confidence limits in data frame
lower1 <- data.frame(cloamy$lower)
upper1 <- data.frame(cloamy$upper)

lower2 <- data.frame(csgp$lower)
upper2 <- data.frame(csgp$upper)

lower3 <- data.frame(cshdi$lower)
upper3 <- data.frame(cshdi$upper)

lower4 <- data.frame(ccrp$lower)
upper4 <- data.frame(ccrp$upper)

###PLOT with regular plot###
summary(g.35)
range(new.crp$crp)
par(mfrow=c(2,2))

plot(new.shdi$shdi, cshdi$Predicted, ylim=c(0,.6),xlim=c(0.8281973,2.1498729),type="n", ylab= "Predicted Occupancy", xlab="SHDIw900", font.lab=2, cex.lab=1.5)
lines(new.shdi$shdi, cshdi$Predicted, lty = 1, lwd=2.5, xlab = "", ylab = "")
lines(new.shdi$shdi, cshdi$upper, lty = 2, lwd=2.5, xlab = "", ylab = "")
lines(new.shdi$shdi, cshdi$lower, lty = 2, lwd=2.5, xlab = "", ylab = "")

plot(new.sgp$sgp, csgp$Predicted, ylim=c(0,.6),xlim=c(0,.9579828), type="n", ylab= "Predicted Occupancy", xlab="SGPw150", font.lab=2, cex.lab=1.5)
lines(new.sgp$sgp, csgp$Predicted, lty = 1, lwd=2.5, xlab = "", ylab = "")
lines(new.sgp$sgp, csgp$upper, lty = 2, lwd=2.5, xlab = "", ylab = "")
lines(new.sgp$sgp, csgp$lower, lty = 2, lwd=2.5, xlab = "", ylab = "")

plot(new.loamy$loamy, cloamy$Predicted, ylim=c(0,.6), xlim=c(0,1), type="n", ylab= "Predicted Occupancy", xlab="LTw50", font.lab=2, cex.lab=1.5)
lines(new.loamy$loamy, cloamy$Predicted, lty = 1, lwd=2.5, xlab = "", ylab = "")
lines(new.loamy$loamy, cloamy$upper, lty = 2, lwd=2.5, xlab = "", ylab = "")
lines(new.loamy$loamy, cloamy$lower, lty = 2, lwd=2.5, xlab = "", ylab = "")

plot(new.crp$crp, ccrp$Predicted, ylim=c(0,.6),xlim=c(0,.9329637), type="n", ylab= "Predicted Occupancy", xlab="CRPw500", font.lab=2, cex.lab=1.5)
lines(new.crp$crp, ccrp$Predicted, lty = 1, lwd=2.5, xlab = "", ylab = "")
lines(new.crp$crp, ccrp$upper, lty = 2, lwd=2.5, xlab = "", ylab = "")
lines(new.crp$crp, ccrp$lower, lty = 2, lwd=2.5, xlab = "", ylab = "")

#Unlist all covariate data frames
#Loamy
new.loamy2 = as.numeric( unlist(new.loamy2))
predict3 = as.numeric( unlist(predict3))
lower1 = as.numeric( unlist(lower1))
upper1 = as.numeric( unlist(upper1))

LoamyTable = data.frame(new.loamy2,predict3, lower1, upper1)
range(LoamyTable$upper1)

#SGP
new.sgp2 = as.numeric( unlist(new.sgp2))
predict5 = as.numeric( unlist(predict5))
lower2 = as.numeric( unlist(lower2))
upper2 = as.numeric( unlist(upper2))

SGPTable = data.frame(new.sgp2,predict5, lower2, upper2)
range(SGPTable$upper2)

#SHDI
new.shdi2 = as.numeric( unlist(new.shdi2))
predict7 = as.numeric( unlist(predict7))
lower3 = as.numeric( unlist(lower3))
upper3 = as.numeric( unlist(upper3))

SHDITable = data.frame(new.shdi2,predict7, lower3, upper3)
range(SHDITable$upper3)

#CRP
new.crp2 = as.numeric(unlist(new.crp2))
predict9 = as.numeric(unlist(predict9))
lower4 = as.numeric(unlist(lower4))
upper4 = as.numeric(unlist(upper4))

CRPTable = data.frame(new.crp2,predict9, lower4, upper4)
range(CRPTable$upper4)


## Start code for plotting graphs

#### OCCUPANCY ####

#define font for plot
black.bold.text <- element_text(face = "bold", color = "black", size="20")

#plot predictions
p <- ggplot(data=LoamyTable, aes(x = new.loamy1, y = predict2)) + geom_line(size=1.9, colour = "black") +
  geom_ribbon(aes(ymin = lower1, ymax = upper1), alpha = 0.2) + xlab("Proportion of Loamy Tableland Soil") + ylab("") 


#make graph pretty
p2 = p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(5,40,10,0))) + 
  theme(axis.text.y=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(10,5,0,10))) + 
  theme(axis.title = black.bold.text) + coord_cartesian(xlim = c(0, 1), ylim = c(0, .6), expand = FALSE) + theme(plot.margin = margin(30,17,15,15)) 
#plot graph
p
p2

p3 <- ggplot(data=SGPTable, aes(x = new.sgp1, y = predict4)) + geom_line(size=1.9, colour = "black") +
  geom_ribbon(aes(ymin = lower2, ymax = upper2), alpha = 0.2) + xlab("Proportion of Shortgrass Prairie") + ylab("") 


#make graph pretty
p4 = p3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(5,40,10,0))) + 
  theme(axis.text.y=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(10,5,0,10))) + 
  theme(axis.title = black.bold.text) + coord_cartesian(xlim = c(0, 1), ylim = c(0, .6), expand = FALSE) + theme(plot.margin = margin(30,17,15,15)) 
#plot graph
p3
p4

#plot predictions
p5 <- ggplot(data=SHDITable, aes(x = new.shdi1, y = predict6)) + geom_line(size=1.9, colour = "black") +
  geom_ribbon(aes(ymin = lower3, ymax = upper3), alpha = 0.2) + xlab("Habitat Diversity") + ylab("") 


#make graph pretty
p6 = p5 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(5,40,10,0))) + 
  theme(axis.text.y=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(10,5,0,10))) + 
  theme(axis.title = black.bold.text) + coord_cartesian(xlim = c(0.8281973,2.1498729), ylim = c(0, .6), expand = FALSE) + theme(plot.margin = margin(30,17,15,15)) 
#plot graph
p5
p6

p7 <- ggplot(data=CRPTable, aes(x = new.crp1, y = predict8)) + geom_line(size=1.9, colour = "black") +
  geom_ribbon(aes(ymin = lower4, ymax = upper4), alpha = 0.2) + xlab("Proportion of CRP") + ylab("") 


#make graph pretty
p8 = p7 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(5,40,10,0))) + 
  theme(axis.text.y=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(10,5,0,10))) + 
  theme(axis.title = black.bold.text) + coord_cartesian(xlim = c(0, 1), ylim = c(0, .6), expand = FALSE) + theme(plot.margin = margin(30,17,15,15)) 
#plot graph
p7
p8

gA <- ggplot_gtable(ggplot_build(p2))
gB <- ggplot_gtable(ggplot_build(p4))
gC <- ggplot_gtable(ggplot_build(p6))
gD <- ggplot_gtable(ggplot_build(p8))

grid.arrange(gC, gB, gA,  gD, ncol = 2) #here, you can fix the 'nrow' or 'ncol' to arrange the plots the way you want

#### Colonization ####

#define font for plot
black.bold.text <- element_text(face = "bold", color = "black", size="20")

#plot predictions
p <- ggplot(data=LoamyTable, aes(x = new.loamy1, y = predict2)) + geom_line(size=1.9, colour = "black") +
  geom_ribbon(aes(ymin = lower1, ymax = upper1), alpha = 0.2) + xlab("Proportion of Loamy Tableland Soil") + ylab("") 


#make graph pretty
p2 = p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(5,40,10,0))) + 
  theme(axis.text.y=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(10,5,0,10))) + 
  theme(axis.title = black.bold.text) + coord_cartesian(xlim = c(0, 1), ylim = c(0, .6), expand = FALSE) + theme(plot.margin = margin(30,17,15,15)) 
#plot graph
p
p2



#plot predictions
p5 <- ggplot(data=SHDITable, aes(x = new.shdi1, y = predict6)) + geom_line(size=1.9, colour = "black") +
  geom_ribbon(aes(ymin = lower3, ymax = upper3), alpha = 0.2) + xlab("Habitat Diversity") + ylab("") 


#make graph pretty
p6 = p5 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(5,40,10,0))) + 
  theme(axis.text.y=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(10,5,0,10))) + 
  theme(axis.title = black.bold.text) + coord_cartesian(xlim = c(0.8281973,2.1498729), ylim = c(0, .6), expand = FALSE) + theme(plot.margin = margin(30,17,15,15)) 
#plot graph
p5
p6



p7 <- ggplot(data=CRPTable, aes(x = new.crp1, y = predict8)) + geom_line(size=1.9, colour = "black") +
  geom_ribbon(aes(ymin = lower4, ymax = upper4), alpha = 0.2) + xlab("Proportion of CRP") + ylab("") 


#make graph pretty
p8 = p7 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(5,40,10,0))) + 
  theme(axis.text.y=element_text(family="sans", face="bold", colour="black", size="15", margin=margin(10,5,0,10))) + 
  theme(axis.title = black.bold.text) + coord_cartesian(xlim = c(0, 1), ylim = c(0, .6), expand = FALSE) + theme(plot.margin = margin(30,17,15,15)) 
#plot graph
p7
p8

gA <- ggplot_gtable(ggplot_build(p2))
gC <- ggplot_gtable(ggplot_build(p6))
gD <- ggplot_gtable(ggplot_build(p8))

grid.arrange(gC, gA,  gD, ncol = 2) #here, you can fix the 'nrow' or 'ncol' to arrange the plots the way you want

############### Extinction
### We didn't plot Extinction as our results were not significant, but the same process can be repeated