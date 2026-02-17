#Changes in PLant cover in the ITEX network- the role of mycorrhizal type
#Lisa Pilkinton
#16/2/26 - edit to remove juniper

#clear the brain
rm (list=ls())
#setwd("M:/") #hash out for use remotely

library(tidyverse)
library(cowplot)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)
library(broom)
library(extrafont)
library(ggeffects)
library(brms)
library(lme4)
#library(cmdstanr)
library(lmerTest)
#install.packages("cmdstanr") #this will not install on this version of R

#load("New_R_Lisa/MYC_ARCTIC_SP_AN.RData") #75470 (for comp use )
load("~/New_R_Lisa/MYC_ARCTIC_SP_AN.RData") #(for remote use)

#lets get a dataset that only has sites that include OTC's

#lets get just warming data, then we can chack how many sitesubsites have warming data- so we have  reference for later.

unique(MYC_ARCTIC_ix$TreatmentGroup)

MYC_ARCTIC_Test<-MYC_ARCTIC_ix %>% filter(TreatmentGroup !="Control") #13688
unique(MYC_ARCTIC_Test$SiteSubsite) #48
#"FAROE:SORNFELLI"               "KLUANE:PIKA"                   "ABISKO:ABISKOWET_MULTE"        "ABISKO:PEATLAND"               "WOLFCREEK:WOLFCREEK"          
# "ZACKENBERG:SALIX_SITE"         "ABISKO:lower_forest"           "ABISKO:lower_tundra"           "ABISKO:upper_forest"           "ABISKO:upper_tundra"          
# "ADVENT:ADVENT_1"               "ADVENT:ADVENT_2"               "ADVENT:ADVENT_3"               "ADVENT:ADVENT_4"               "ADVENT:ADVENT_5"              
# "ADVENT:ADVENT_6"               "ADVENT:ADVENT_7"               "KILPISJARVI:KILPISJARVI"       "LATNJA:DRY_HEATH"              "LATNJA:DRY_MEADOW"            
#"ABISKO:ABISKOWET_FLUX"         "ZACKENBERG:CASSIOPE_SITE_FLUX" "ZACKENBERG:SALIX_SITE_FLUX"    "ABISKO:ABISKODRY_PADDUS"       "SADVENT:MES_POINTFRAME"       
#"LATNJA:POOR_HEATH"             "LATNJA:RICH_MEADOW"            "SADVENT:WET_POINTFRAME"        "ZACKENBERG:CASSIOPE_SITE"      "TOOLIK:MOIST"                 
#"TOOLIK:DRY"                    "NARSARSUAQ:HIGH_ELEVATION"     "ATQASUK:AD"                    "ATQASUK:AW"                    "BARROW:BD"                    
#"BARROW:BW"                     "ALEXFIORD:MEADOW"              "ALEXFIORD:WILLOW"              "ALEXFIORD:CASSIOPE"            "ALEXFIORD:FERT"               
#"ALEXFIORD:DOMEDOLOMITE"        "ALEXFIORD:DOMEGRANITE"         "ALEXFIORD:DRYAS"               "ENDALEN:DRY-L"                 "ENDALEN:CAS-L"                
#"ENDALEN:BIS-L"                 "AUDKULUHEIDI:BETULAHEATH"      "THINGVELLIR:MOSS HEATH" 


MYC_ARCTIC_TREAT <- MYC_ARCTIC_ix %>% #27104
  group_by(SiteSubsite, YEAR) %>%
  filter(all(c("Warming", "Control") %in% TreatmentGroup)) %>%
  ungroup()

unique(MYC_ARCTIC_TREAT$SiteSubsite) #48- Great

#I will get a duration of monitoring before filtering out all observation but the last year.

# duration of monitoring # we have plots from 1 year of treatment to 24 year of treatment!
duration_per_plot <- MYC_ARCTIC_TREAT %>%
  group_by(SiteSubsitePlot) %>%
  summarize(
    start_year = min(YEAR),
    end_year   = max(YEAR),
    duration   = end_year - start_year
  ) %>%
  ungroup()



# Add to the main dataset
MYC_ARCTIC_TREAT <- MYC_ARCTIC_TREAT %>%
  left_join(duration_per_plot, by = "SiteSubsitePlot")

#Remove plots with duration < 5 years, consistent with CTL data too # 25276
MYC_ARCTIC_TREAT <- MYC_ARCTIC_TREAT %>%
  filter(duration >= 5)

#keep only the final year #3978 ( we will keep all years for the interaction model...so hash this out for 2 and 3)
MYC_ARCTIC_TREAT<- MYC_ARCTIC_TREAT %>%
 group_by(SiteSubsite) %>%
 filter(YEAR == max(YEAR, na.rm = TRUE)) %>%
ungroup()

#add adjusted year relative for interaction model
# use the middle year as the relative instead of start, so that 
# each interecept of each subsite is centered around the middle of the survey
# Same as in Bjorkman 2018
MYC_ARCTIC_TREAT <- MYC_ARCTIC_TREAT%>%
  group_by(SiteSubsitePlot) %>%
  mutate(
    YEAR_relative_2 = YEAR - mean(range(YEAR)),        # or use your existing YEAR_relative
    # scale YEAR_relative globally for numerical stability:
    YEAR_relative_sc = as.numeric(scale(YEAR_relative_2)) #Lisa editetd to YEAR_relative_2, rather than YEAR_relative
  )


#create propertion

MYC_ARCTIC_TREAT<-mutate(MYC_ARCTIC_TREAT,PropCov=Cover/100)#make a proportion
MYC_ARCTIC_TREAT$GMYCPFUNC <- as.factor(MYC_ARCTIC_TREAT$GMYCPFUNC)
MYC_ARCTIC_TREAT$MOISTURE <- as.factor(MYC_ARCTIC_TREAT$MOISTURE)


#### Checking the proportional values

sum(MYC_ARCTIC_TREAT$PropCov>= 1)
# There is only 11 100% in the whole dataset, but they need 
# there own "one model" which increase the complexity a bunch
# I suggest turning ones into the biggest non-one cover value
(new_val <- max(MYC_ARCTIC_TREAT$PropCov[MYC_ARCTIC_TREAT$PropCov < 1]))
MYC_ARCTIC_TREAT$PropCov <- ifelse(MYC_ARCTIC_TREAT$PropCov>= 1,new_val,MYC_ARCTIC_TREAT$PropCov)
sum(MYC_ARCTIC_TREAT$PropCov>= 1) # worked like a charm


#remove the group EcMAM:SEVER #25262
MYC_ARCTIC_TREAT<-MYC_ARCTIC_TREAT %>% filter(GMYCPFUNC !="EcMAM:SEVER") #13688
unique(MYC_ARCTIC_Test$SiteSubsite)

#remove for treat_6 without am:SDECI
MYC_ARCTIC_TREAT<-MYC_ARCTIC_TREAT %>% filter(GMYCPFUNC !="AM:SDECI") #13688
unique(MYC_ARCTIC_Test$SiteSubsite)

## I suggest setting more explciit prior, this will also help the model to fit
my_prior <- set_prior("normal(-2,0.5)",class = "Intercept",dpar = "mu")# Cover Intercept
my_prior <- c(my_prior,set_prior("normal(-2,1)",class = "Intercept",dpar = "zi")) # absence intercept

my_prior <- c(my_prior,prior(gamma(0.01, 0.02), class = phi))
my_prior <- c(prior(gamma(0.01, 0.02), class = phi))


#hash this out, but use again if we try the removal of EcMAMSEVER - Juniper
#TREAT_MODEL_5 <- brm( #I use bf() to set up several formula
# the cover model
 # bf(PropCov ~ 0 + TreatmentGroup * GMYCPFUNC + (1 + TreatmentGroup | SiteSubsite:GMYCPFUNC) + (1|SiteSubsitePlot:GMYCPFUNC)+ (1 + TreatmentGroup|SPECIES_NAME:GMYCPFUNC),
# the absence presence model
  #  zi~ 0 + TreatmentGroup * GMYCPFUNC + (1 + TreatmentGroup | SiteSubsite:GMYCPFUNC) + (1|SiteSubsitePlot:GMYCPFUNC)+ (1 + TreatmentGroup |SPECIES_NAME:GMYCPFUNC) ),
# data = MYC_ARCTIC_TREAT,
 # family = zero_inflated_beta(),
  #prior = my_prior,  
#  chains = 3,
 # iter = 1250, #in the exploratory phase Less iterations and 3 chains is enough
  #warmup = 250,# again, in exploratory phase less burnin
#  cores = 3, # A core is used per chain, compute three time faster !
 # threads = threading(4), # hom many core per core (confusing I know) to speed up even more, use less core that parallel::detectCores()
# control = list(adapt_delta = 0.99), # this increase the computation time like crazy, and often better model specification can fix the issues
#  backend = "cmdstanr", # you can comment this or install cmdstanr, more stable API <- wouldn't install on my R
# file = "TREAT_MODEL_5"   # <-- autosaves every few minutes
#)


#now run with all years and interaction with year

#TREAT_MODEL_3 <- brm(
#  bf(PropCov ~ 0 + TreatmentGroup * YEAR_relative_2 * GMYCPFUNC + (1 + TreatmentGroup | SiteSubsite:GMYCPFUNC) + (1|SiteSubsitePlot:GMYCPFUNC)+ (1 + TreatmentGroup|SPECIES_NAME:GMYCPFUNC),
 #    #the absence presence model
  #   zi~ 0 + TreatmentGroup * GMYCPFUNC + (1 + TreatmentGroup | SiteSubsite:GMYCPFUNC) + (1|SiteSubsitePlot:GMYCPFUNC)+ (1 + TreatmentGroup |SPECIES_NAME:GMYCPFUNC) ),
  #data = MYC_ARCTIC_TREAT,
  #family = zero_inflated_beta(),
#  prior = my_prior,  
 # chains = 3,
  #iter = 1250, #in the exploratory phase Less iterations and 3 chains is enough
#  warmup = 250,# again, in exploratory phase less burnin
 # cores = 3, # A core is used per chain, compute three time faster !
  #threads = threading(4), # hom many core per core (confusing I know) to speed up even more, use less core that parallel::detectCores()
#  control = list(adapt_delta = 0.99), # this increase the computation time like crazy, and often better model specification can fix the issues
  #  backend = "cmdstanr", # you can comment this or install cmdstanr, more stable API <- wouldn't install on my R
 # file = "TREAT_MODEL_3"   # <-- autosaves every few minutes
#)



#hash this out, but use again if we try the removal of EcMAMSEVER - Juniper
TREAT_MODEL_6 <- brm( #I use bf() to set up several formula
  # the cover model
  bf(PropCov ~ 0 + TreatmentGroup * GMYCPFUNC + (1 + TreatmentGroup | SiteSubsite:GMYCPFUNC) + (1|SiteSubsitePlot:GMYCPFUNC)+ (1 + TreatmentGroup|SPECIES_NAME:GMYCPFUNC),
     # the absence presence model
     zi~ 0 + TreatmentGroup * GMYCPFUNC + (1 + TreatmentGroup | SiteSubsite:GMYCPFUNC) + (1|SiteSubsitePlot:GMYCPFUNC)+ (1 + TreatmentGroup |SPECIES_NAME:GMYCPFUNC) ),
  data = MYC_ARCTIC_TREAT,
  family = zero_inflated_beta(),
  prior = my_prior,  
  chains = 3,
  iter = 1250, #in the exploratory phase Less iterations and 3 chains is enough
  warmup = 250,# again, in exploratory phase less burnin
  cores = 3, # A core is used per chain, compute three time faster !
  threads = threading(4), # hom many core per core (confusing I know) to speed up even more, use less core that parallel::detectCores()
  # control = list(adapt_delta = 0.99), # this increase the computation time like crazy, and often better model specification can fix the issues
  #  backend = "cmdstanr", # you can comment this or install cmdstanr, more stable API <- wouldn't install on my R
  file = "TREAT_MODEL_6"   # <-- autosaves every few minutes
)
