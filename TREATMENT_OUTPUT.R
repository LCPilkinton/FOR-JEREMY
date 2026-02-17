#Changes in PLant cover in the ITEX network- the role of my: corhizal type
#Lisa Pilkinton
#11/02/2026- Outputs of a treatment effect model


#clear the brain
rm (list=ls())
#setwd("M:/R_Lisa")



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
library(lmerTest)
library(tidybayes)
library(broom.mixed)
library(forcats)
library(viridis)
library(bayesplot)

#create a theme for plots
theme_clean <- function(){
  theme_bw()+
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "plain"),             
          axis.title.y = element_text(size = 14, face = "plain"),             
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 10, face = "italic"),          
          legend.position = "right")
}
#add in the data used


load("~/New_R_Lisa/MYC_ARCTIC_SP_AN.RData") #75470


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

#keep only the final year #3978 - we have only done this on TREAT_1 and TREAT_5 and TREAT_6
#MYC_ARCTIC_TREAT<- MYC_ARCTIC_TREAT %>%
 # group_by(SiteSubsite) %>%
  #filter(YEAR == max(YEAR, na.rm = TRUE)) %>%
  #ungroup()


#add adjusted year relative for interaction models TREAT_2 and TREAT_3
MYC_ARCTIC_TREAT <- MYC_ARCTIC_TREAT%>%group_by(SiteSubsitePlot) %>%
  mutate(
    YEAR_relative_2 = YEAR - mean(range(YEAR)),        # or use your existing YEAR_relative
    # scale YEAR_relative globally for numerical stability:
    YEAR_relative_sc = as.numeric(scale(YEAR_relative_2)) #Lisa editetd to YEAR_relative_2, rather than YEAR_relative
  )

#remove juniper for TREAT_3 and TREAT_5 and TREAT_6
#remove the group EcMAM:SEVER #25262
MYC_ARCTIC_TREAT<-MYC_ARCTIC_TREAT %>% filter(GMYCPFUNC !="EcMAM:SEVER") #13688
unique(MYC_ARCTIC_Test$SiteSubsite)

#also for TREAT_6
MYC_ARCTIC_AM_DECI<-MYC_ARCTIC_TREAT %>% filter(GMYCPFUNC =="AM:SDECI") #13688
unique(MYC_ARCTIC_Test$SiteSubsite)


#create propertion


MYC_ARCTIC_TREAT<-mutate(MYC_ARCTIC_TREAT,PropCov=Cover/100)#make a proportion
MYC_ARCTIC_TREAT$GMYCPFUNC <- as.factor(MYC_ARCTIC_TREAT$GMYCPFUNC)
MYC_ARCTIC_TREAT$MOISTURE <- as.factor(MYC_ARCTIC_TREAT$MOISTURE)


#### Adjust the proportional values

sum(MYC_ARCTIC_TREAT$PropCov>= 1)
(new_val <- max(MYC_ARCTIC_TREAT$PropCov[MYC_ARCTIC_TREAT$PropCov < 1]))
MYC_ARCTIC_TREAT$PropCov <- ifelse(MYC_ARCTIC_TREAT$PropCov>= 1,new_val,MYC_ARCTIC_TREAT$PropCov)
sum(MYC_ARCTIC_TREAT$PropCov>= 1) 

############################################################################################################################30/1 SPECIES as random effect

fitmod_TREAT <- readRDS("New_R_Lisa/TREAT_MODEL_1.rds") # we add the {R_Lisa/} for remote desktop <- Model originally in Treatment_models

summary(fitmod_TREAT) #nope not great- its the EcMSM:SEVER that is a problem- juniper


fitmod_TREAT_2 <- readRDS("New_R_Lisa/TREAT_MODEL_2.rds") # we add the {New_R_Lisa/} for remote desktop <- Model originally in Treatment_models

summary(fitmod_TREAT_2) #no- not good at all, see without juniper

fitmod_TREAT_3 <- readRDS("New_R_Lisa/TREAT_MODEL_3.rds") # we add the {New_R_Lisa/} for remote desktop <- Model originally in Treatment_model_minusjuniper

summary(fitmod_TREAT_3) #no- not good at all, see without juniper and without interaction

fitmod_TREAT_5 <- readRDS("New_R_Lisa/TREAT_MODEL_5.rds") # we add the {New_R_Lisa/} for remote desktop <- Model originally in Treatment_model_minusjuniper

summary(fitmod_TREAT_5) #no- not good at all, see without juniper - looks like AM deci is a problem- have a look at it


fitmod_TREAT_6<- readRDS("New_R_Lisa/TREAT_MODEL_6.rds") # we add the {New_R_Lisa/} for remote desktop <- Model originally in Treatment_modelsminusjuniper
summary(fitmod_TREAT_6) #This is much better- will look at the draws next. possibly be better with more iterations


