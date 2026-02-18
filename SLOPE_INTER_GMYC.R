#Changes in PLant cover in the ITEX network- the role of mycorhizal type
#Lisa Pilkinton
#5/2/26 Species included with a random slope and random intercept looking at mycorrhizal grouping only

#clear the brain
rm (list=ls())
#setwd("/ITEX/ITEX_data")

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

load("~/New_R_Lisa/MYC_ARCTIC_SP_AN.RData") #75470


#lets look just at controls

MYC_ARCTIC_CTL<-MYC_ARCTIC_ix %>%
  filter(TreatmentGroup !="Warming")#----------------------------#61782obs


str(MYC_ARCTIC_CTL)

# duration of monitoring
duration_per_plot <- MYC_ARCTIC_CTL %>%
  group_by(SiteSubsitePlot) %>%
  summarize(
    start_year = min(YEAR),
    end_year   = max(YEAR),
    duration   = end_year - start_year
  ) %>%
  ungroup()

# Add to the main dataset
MYC_ARCTIC_CTL <- MYC_ARCTIC_CTL %>%
  left_join(duration_per_plot, by = "SiteSubsitePlot")

#Remove plots with duration < 5 years#...........60357
MYC_ARCTIC_CTL <- MYC_ARCTIC_CTL %>%
  filter(duration >= 5)

str(MYC_ARCTIC_CTL)

#create propertion

MYC_ARCTIC_CTL<-mutate(MYC_ARCTIC_CTL,PropCov=Cover/100)#make a proportion
MYC_ARCTIC_CTL$GMYCPFUNC <- as.factor(MYC_ARCTIC_CTL$GMYCPFUNC)
MYC_ARCTIC_CTL$MOISTURE <- as.factor(MYC_ARCTIC_CTL$MOISTURE)


# create a YEAR_relative variable and scale it
MYC_ARCTIC_CTL <- MYC_ARCTIC_CTL %>%
  mutate(
    YEAR_relative = YEAR - start_year,        # or use your existing YEAR_relative
    # scale YEAR_relative globally for numerical stability:
    YEAR_relative_sc = as.numeric(scale(YEAR_relative))
  )


#### Start of Jeremy's suggestions

# use the middle year as the relative instead of start, so that 
# each interecept of each subsite is centered around the middle of the survey
# Same as in Bjorkman 2018
MYC_ARCTIC_CTL <- MYC_ARCTIC_CTL%>%group_by(SiteSubsitePlot) %>%
  mutate(
    YEAR_relative_2 = YEAR - mean(range(YEAR)),        # or use your existing YEAR_relative
    # scale YEAR_relative globally for numerical stability:
    YEAR_relative_sc = as.numeric(scale(YEAR_relative_2)) #Lisa editetd to YEAR_relative_2, rather than YEAR_relative
  )


sum(MYC_ARCTIC_CTL$PropCov>= 1)
# There is only 38 100% in the whole dataset, but they need 
# there own "one model" which increase the complexity a bunch
# I suggest turning ones into the biggest non-one cover value
(new_val <- max(MYC_ARCTIC_CTL$PropCov[MYC_ARCTIC_CTL$PropCov < 1]))
MYC_ARCTIC_CTL$PropCov <- ifelse(MYC_ARCTIC_CTL$PropCov>= 1,new_val,MYC_ARCTIC_CTL$PropCov)
sum(MYC_ARCTIC_CTL$PropCov>= 1) # worked like a charm

## I suggest setting more explciit prior, this will also help the model to fit
my_prior <- set_prior("normal(-2,0.5)",class = "Intercept",dpar = "mu")# Cover Intercept
my_prior <- c(my_prior,set_prior("normal(-2,1)",class = "Intercept",dpar = "zi")) # absence intercept

my_prior <- c(my_prior,prior(gamma(0.01, 0.02), class = phi))
my_prior <- c(prior(gamma(0.01, 0.02), class = phi))

# now we make the model

GMYC_SP_RAN_SL_INT <- brm( #I use bf() to set up several formula
  # the cover model
  bf(PropCov ~ 0 + YEAR_relative_sc * GMYC + (1 + YEAR_relative_sc | SiteSubsite:GMYC) + (1|SiteSubsitePlot:GMYC)+ (1 + YEAR_relative_sc |SPECIES_NAME:GMYC),
     # the absence presence model
     zi~ 0 + YEAR_relative_sc * GMYC + (1 + YEAR_relative_sc | SiteSubsite:GMYC) + (1|SiteSubsitePlot:GMYC)+ (1 + YEAR_relative_sc |SPECIES_NAME:GMYC) ),
  data = MYC_ARCTIC_CTL,
  family = zero_inflated_beta(),
  prior = my_prior,  
  chains = 3,
  iter = 1250, #in the exploratory phase Less iterations and 3 chains is enough
  warmup = 250,# again, in exploratory phase less burnin
  cores = 3, # A core is used per chain, compute three time faster !
  threads = threading(4), # hom many core per core (confusing I know) to speed up even more, use less core that parallel::detectCores()
  # control = list(adapt_delta = 0.99), # this increase the computation time like crazy, and often better model specification can fix the issues
  #  backend = "cmdstanr", # you can comment this or install cmdstanr, more stable API <- wouldn't install on my R
  file = "GMYC_SP_RAN_SL_INT"   # <-- autosaves every few minutes
)

