#Changes in PLant cover in the ITEX network- the role of mycorhizal type
#Lisa Pilkinton
#January- Script for Jeremy

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

library(lmerTest)


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

#########################################################################      SET UP THE DATA        ##########################################################
load("~/R_Lisa/ITEX/ITEX_data/MYC_ARCTIC_SP_AN.RData") #75470


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

#create proportion

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

############################################    SET UP MODEL, THIS IS FOR THE FULL GMYCPFUC GROUP, MODEL NAME : model3_brms               ####################################

#Make the model I want just GMYCPFUNC*time, this will give each GMYCPFUNC at each sitesubsite.

model3_brms <- brm(
  PropCov ~ 0 + YEAR_relative_sc * GMYCPFUNC + (1 + YEAR_relative_sc | SiteSubsite:GMYCPFUNC),
  data = MYC_ARCTIC_CTL,
  family = zero_one_inflated_beta(),
  prior = prior(gamma(0.01, 0.01), class = phi),  
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99),
  file = "model3_brms_autosave"   # <-- autosaves every few minutes
)

# Save final object 
write_rds(model3_brms, "model3_brms.rds")

##########################################    MODEL OUTPUT FOR: FULL GMYCPFUC GROUP, MODEL NAME : model3_brms   ############################################

#----Model: PropCov ~ 0 + YEAR_relative_sc * GMYCPFUNC + (1 + YEAR_relative_sc | SiteSubsite:GMYCPFUNC)

#we have model saved so can load it here
fitmod3 <- readRDS("R_Lisa/model3_brms.rds") # we add the {R_Lisa/} for remote desktop

summary(fitmod3) #this looks good so far

plot(fitmod3) # all good

pp_check(fitmod3, ndraws = 100) # this looks okay


#-----------runing things to get the slopes for each GMYCPFUNC at each subsite
library(tidybayes)
library(dplyr)
library(stringr)


#I need to get the fixed slope for the interaction of year and gmycpfunc, then tHE random efect of sitesubsite:Gmycpfunc and together they will give my gmycpfunc over time for each sitesubsite

fitmod3%>% get_variables() %>% grep("^b_", ., value = TRUE)
# "b_YEAR_relative_sc"                           "b_GMYCPFUNCAM:FORB"                           "b_GMYCPFUNCAM:GRAMMINOID"                    
#"b_GMYCPFUNCAM:SDECI"                          "b_GMYCPFUNCAM:SEVER"                          "b_GMYCPFUNCEcM:FORB"                         
#"b_GMYCPFUNCEcM:SDECI"                         "b_GMYCPFUNCEcM:SEVER"                         "b_GMYCPFUNCEcMAM:FORB"                       
# "b_GMYCPFUNCEcMAM:GRAMMINOID"                  "b_GMYCPFUNCEcMAM:SDECI"                       "b_GMYCPFUNCEcMAM:SEVER"                      
# "b_GMYCPFUNCEcMErM:SDECI"                      "b_GMYCPFUNCErM:SDECI"                         "b_GMYCPFUNCErM:SEVER"                        
# "b_GMYCPFUNCNM:FORB"                           "b_GMYCPFUNCNMAM:FORB"                         "b_GMYCPFUNCNMAM:GRAMMINOID"                  
# "b_YEAR_relative_sc:GMYCPFUNCAM:GRAMMINOID"    "b_YEAR_relative_sc:GMYCPFUNCAM:SDECI"         "b_YEAR_relative_sc:GMYCPFUNCAM:SEVER"        
# "b_YEAR_relative_sc:GMYCPFUNCEcM:FORB"         "b_YEAR_relative_sc:GMYCPFUNCEcM:SDECI"        "b_YEAR_relative_sc:GMYCPFUNCEcM:SEVER"       
# "b_YEAR_relative_sc:GMYCPFUNCEcMAM:FORB"       "b_YEAR_relative_sc:GMYCPFUNCEcMAM:GRAMMINOID" "b_YEAR_relative_sc:GMYCPFUNCEcMAM:SDECI"     
# "b_YEAR_relative_sc:GMYCPFUNCEcMAM:SEVER"      "b_YEAR_relative_sc:GMYCPFUNCEcMErM:SDECI"     "b_YEAR_relative_sc:GMYCPFUNCErM:SDECI"       
#"b_YEAR_relative_sc:GMYCPFUNCErM:SEVER"        "b_YEAR_relative_sc:GMYCPFUNCNM:FORB"          "b_YEAR_relative_sc:GMYCPFUNCNMAM:FORB"       
#"b_YEAR_relative_sc:GMYCPFUNCNMAM:GRAMMINOID" 


library(posterior)
#I need this to see the random effect interaction,
Row.names_3 <- as.data.frame(fitmod3)
#r_SiteSubsite:GMYCPFUNC[.........] #in the brackets site:subsite_gmyc:pfunc, year_rel_sc # this means we divide where there is a comma, then divide by underscore.




# 1) get posterior draws for every interaction slope
# here we extract all parameters starting with b_YEAR_relative_sc:GMYCPFUNC, and we also need b_YEAR_relative_sc, because it is the baseline we need to use
options(scipen = 999)
interaction_terms_3 <- names(as_draws_df(fitmod3)) %>% #17 terms
  grep(
    "^b_YEAR_relative_sc$|^b_YEAR_relative_sc:GMYCPFUNC",
    .,
    value = TRUE
  )

all_fixed_draws_3 <- fitmod3 %>%
  spread_draws(!!!syms(interaction_terms_3))

#The baseline "b_YEAR_relative_sc" is the b_YEAR_relative_sc:AM:FORB interaction group- we want all but that group.
INT_cols <- interaction_terms_3[
  interaction_terms_3!= "b_YEAR_relative_sc"
]

# 2) add baseline (b_YEAR_relative_sc) to each GMYCPFUNC-specific slope
added_fixed_draws_3 <- all_fixed_draws_3 %>%
  mutate(
    across(
      all_of(INT_cols),
      ~ .x + b_YEAR_relative_sc
    )
  )

#We can make these tidy by removing the stuuctrre from teh model output ans calling year_rel AM:FORB
tidy_draws <- added_fixed_draws_3 %>%
  rename_with(
    ~ str_remove(.x, "^b_YEAR_relative_sc:GMYCPFUNC"),
    -c(.chain, .iteration, .draw)
  ) %>%
  rename(
    `AM:FORB` = b_YEAR_relative_sc
  )

#W can summarise fixed effects:
#fixed_slope_summary_3 <- tidy_draws %>%
# pivot_longer(
#  cols = -c(.chain, .iteration, .draw),
#  names_to = "GMYCPFUNC",
#  values_to = "slope"
#) %>%
# group_by(GMYCPFUNC) %>%
#  summarise(
#    mean_slope = mean(slope),
#    sd_slope   = sd(slope),
#    ci_low     = quantile(slope, 0.025),
#   ci_high    = quantile(slope, 0.975),
#   .groups = "drop"
# )
#

#2 get the posterior draws for the random effects of each subsite, 
#this is what we need to then add to the interaction fixed effect for each GMYCPFUC:YEAR,
#We will end up with a slope value for each GMYCPFUC at each SiteSubsite.


subsite_terms_3 <- names(as_draws_df(fitmod3)) %>% #739 terms
  grep("^r_SiteSubsite:GMYCPFUNC\\[.*,YEAR_relative_sc\\]$", 
       ., value = TRUE)

all_rand_draws_3 <- fitmod3 %>%
  spread_draws(!!!syms(subsite_terms_3))

#we now add together the fixed effect and the random effect for each site

#we pivot longer prior to adding to the random draws
fixed_long_3 <- tidy_draws %>%
  pivot_longer(
    cols = -c(.chain, .iteration, .draw),
    names_to = "GMYCPFUNC",
    values_to = "fixed_slope"
  )

#pivot longer for these too, we also remove the syntax from the model output but put the information into their own columns. 
#This allows us to join by GMYCPFUC and add teh fixed effects to random for each of the GMYCPFUNC groups in each sitesubsite
rand_long_3 <- all_rand_draws_3 %>%
  pivot_longer(
    cols = -c(.chain, .iteration, .draw),
    names_to = "rand_term",
    values_to = "rand_slope"
  ) %>%
  mutate(
    # remove Stan syntax
    rand_term = str_remove(rand_term, "^r_SiteSubsite:GMYCPFUNC\\["),
    rand_term = str_remove(rand_term, ",YEAR_relative_sc\\]$"),
    
    # extract SITE (everything before the first ':')
    SITE = str_extract(rand_term, "^[^:]+"),
    
    # extract PFUNC (everything after the last '_')
    GMYCPFUNC = str_extract(rand_term, "[^_]+$"),
    
    # extract SiteSubsite (everything between first ':' and last '_')
    SiteSubsite = str_remove(rand_term, paste0("^", SITE, ":")) %>%
      str_remove(paste0("_", GMYCPFUNC, "$"))
  )


total_slopes_3<- rand_long_3 %>%
  left_join(fixed_long_3, by = c(".draw", "GMYCPFUNC")) %>%
  mutate(total_slope = fixed_slope + rand_slope)

#now we have the total slope, we can summarise the draws
total_slope_summary_3 <- total_slopes_3 %>%
  group_by(SiteSubsite, GMYCPFUNC) %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )

save(total_slopes_3, file = "R_Lisa/ITEX/ITEX_data/total_slopes_3.RData")

#summarise not by subsite totalty by GNYCPFUNC

total_slope_summary_GMYCPFUNC <- total_slopes_3 %>%
  group_by(GMYCPFUNC) %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )

# We use the total slope summery to show the data, fo rthe plot to look nice, we reorder teh GMYC group, so its shrubs, grammiods, forbs

total_slope_summary_3<- total_slope_summary_3 %>%
  mutate(GMYCPFUNC = fct_relevel(GMYCPFUNC,  
                                 "EcM:SDECI",
                                 "ErM:SDECI",
                                 "EcMErM:SDECI",
                                 "EcMAM:SDECI", 
                                 "AM:SDECI" ,
                                 "EcM:SEVER", 
                                 "ErM:SEVER",
                                 "EcMAM:SEVER",  
                                 "AM:SEVER",
                                 "EcMAM:GRAMMINOID",
                                 "AM:GRAMMINOID",
                                 "NMAM:GRAMMINOID",
                                 "EcM:FORB",
                                 "EcMAM:FORB",
                                 "AM:FORB",
                                 "NMAM:FORB",
                                 "NM:FORB"
  ))

#We make a custon Pallete, so that the mycorrhizal groups share the same colour
custom_palette <- c(
  "EcM:SDECI"="dodgerblue4",
  "ErM:SDECI"="chartreuse4",
  "EcMErM:SDECI"="purple3",
  "EcMAM:SDECI"="deeppink4", 
  "AM:SDECI"="slateblue1" ,
  "EcM:SEVER"="dodgerblue4", 
  "ErM:SEVER"="chartreuse4",
  "EcMAM:SEVER"="deeppink4",  
  "AM:SEVER"="slateblue1",
  "EcMAM:GRAMMINOID"="deeppink4",
  "AM:GRAMMINOID"="slateblue1",
  "NMAM:GRAMMINOID"="lightpink1",
  "EcM:FORB"="dodgerblue4",
  "EcMAM:FORB"="deeppink4",
  "AM:FORB"="slateblue1",
  "NMAM:FORB"="lightpink1",
  "NM:FORB"="deepskyblue2")


####---------------make the plot
(GMYC_AM_JAN <- ggplot(
  total_slope_summary_3,
  aes(
    x = GMYCPFUNC,
    y = mean_slope,
    fill = GMYCPFUNC
  )
) +
    geom_boxplot(
      alpha = 0.5,
      linewidth = 0.3,
      width = 0.7
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    scale_fill_manual(
      values = custom_palette,
      guide = guide_legend(
        override.aes = list(
          shape = NA,
          linewidth = 0,
          alpha = 0.5
        )
      )
    ) +
    theme_clean() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.key.height = unit(0.8, "cm"),
      legend.key.width  = unit(1.2, "cm")
    ) +
    labs(
      x = "Functional group (GMYCPFUNC)",
      y = "Mean site-level slope (PropCov per standardized year)",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("GMYC_AM_JAN.png", plot = GMYC_AM_JAN, width = 16, height = 6, dpi = 300)

############################################    SET UP MODEL, THIS IS FOR JUST PFUNC GROUP, MODEL NAME : model_Pfunc_2               ####################################

#Make the model I want just PFUNC*time, this will give each PFUNC at each sitesubsite.

model_Pfunc_2 <- brm(
  PropCov ~ 0 + PFUNC + YEAR_relative_sc * PFUNC + (1 + YEAR_relative_sc | SiteSubsite:PFUNC),
  data = MYC_ARCTIC_CTL,
  family = zero_one_inflated_beta(),
  prior = prior(gamma(0.01, 0.01), class = phi),  
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99),
  file = "model_Pfunc_2_autosave"   # <-- autosaves every few minutes
)

# Save final object cleanly
write_rds(model_Pfunc_2, "model_Pfunc_2.rds")

##########################################    MODEL OUTPUT FOR: FULL PFUNC GROUP, MODEL NAME : model_Pfunc_2   ############################################


fitmodPfunc_2 <- readRDS("R_Lisa/model_Pfunc_2.rds") # we add the {R_Lisa/} for remote desktop

summary(fitmodPfunc_2 ) #this looks good so far

plot(fitmodPfunc_2) # all good

pp_check(fitmodPfunc_2, ndraws = 100) # this looks okay


#-----------runiing things to get the slopes for each GMYCPFUNC at each subsite
library(tidybayes)
library(dplyr)
library(stringr)


#I need to get the fixed slope for the interaction of year and pfunc, then the random effect of sitesubsite:pfunc and together they will give my gmycpfunc over time for each sitesubsite

fitmodPfunc_2 %>% get_variables() %>% grep("^b_", ., value = TRUE)
# "b_PFUNCFORB"                        "b_PFUNCGRAMMINOID"                  "b_PFUNCSDECI"                       "b_PFUNCSEVER"                      
# "b_YEAR_relative_sc"                 "b_PFUNCGRAMMINOID:YEAR_relative_sc" "b_PFUNCSDECI:YEAR_relative_sc"      "b_PFUNCSEVER:YEAR_relative_sc"     


library(posterior)
#I need this to see therandom effect interaction,
Row.names_pf <- as.data.frame(fitmodPfunc_2)
#r_SiteSubsite:PFUNC[.........] #in the brackets site:subsite_gmyc:pfunc, year_rel_sc # this means we divide where there is a comma, then divide by underscore.




# 1) get posterior draws for every interaction slope
#  here we extract all parameters starting with b_PFUNC.*:YEAR_relative_sc, and we also need b_YEAR_relative_sc, because it is the baseline we need to use

interaction_terms_pf <- names(as_draws_df(fitmodPfunc_2)) %>%
  grep(
    "^b_YEAR_relative_sc$|^b_PFUNC.*:YEAR_relative_sc$",
    .,
    value = TRUE
  )
#spread all the draws for all the groups, including YEAR_rel, as that is FORB
all_fixed_draws_pf <- fitmodPfunc_2 %>%
  spread_draws(!!!syms(interaction_terms_pf))

#get terms but not the year_rel (FORB) we using this term to add and adjust values
pfunc_cols <- interaction_terms_pf[
  interaction_terms_pf != "b_YEAR_relative_sc"
]

# 2) add baseline year slope (FORB) to each PFUNC-specific slope
added_fixed_draws_pf <- all_fixed_draws_pf %>%
  mutate(
    across(
      all_of(pfunc_cols),
      ~ .x + b_YEAR_relative_sc
    )
  )

#tidy names
tidy_draws <- added_fixed_draws_pf %>%
  rename_with(
    ~ str_remove(.x, "^b_PFUNC") %>%
      str_remove(":YEAR_relative_sc$"),
    -c(.chain, .iteration, .draw)
  ) %>%
  rename(`FORB` = b_YEAR_relative_sc)



#summarise fixed effects:
fixed_slope_summary_pf <- tidy_draws %>%
  pivot_longer(
    cols = -c(.chain, .iteration, .draw),
    names_to = "PFUNC",
    values_to = "slope"
  ) %>%
  group_by(PFUNC) %>%
  summarise(
    mean_slope = mean(slope),
    sd_slope   = sd(slope),
    ci_low     = quantile(slope, 0.025),
    ci_high    = quantile(slope, 0.975),
    .groups = "drop"
  )


#2 get the posterior draws for the random effects of each subsite


subsite_terms_pf <- names(as_draws_df(fitmodPfunc_2)) %>%
  grep("^r_SiteSubsite:PFUNC\\[.*,YEAR_relative_sc\\]$", 
       ., value = TRUE)

all_rand_draws_pf <- fitmodPfunc_2 %>%
  spread_draws(!!!syms(subsite_terms_pf))

all_draws_pf <- added_fixed_draws_pf %>%
  left_join(all_rand_draws_pf, by = ".draw")


#we now add together the fixed effect and the random effect for each site
#5/1/26, trying a different adding approach

fixed_long_pf <- tidy_draws %>%
  pivot_longer(
    cols = -c(.chain, .iteration, .draw),
    names_to = "PFUNC",
    values_to = "fixed_slope"
  )


rand_long_pf <- all_rand_draws_pf %>%
  pivot_longer(
    cols = -c(.chain, .iteration, .draw),
    names_to = "rand_term",
    values_to = "rand_slope"
  ) %>%
  mutate(
    # remove Stan syntax
    rand_term = str_remove(rand_term, "^r_SiteSubsite:PFUNC\\["),
    rand_term = str_remove(rand_term, ",YEAR_relative_sc\\]$"),
    
    # extract SITE (everything before the first ':')
    SITE = str_extract(rand_term, "^[^:]+"),
    
    # extract PFUNC (everything after the last '_')
    PFUNC = str_extract(rand_term, "[^_]+$"),
    
    # extract SiteSubsite (everything between first ':' and last '_')
    SiteSubsite = str_remove(rand_term, paste0("^", SITE, ":")) %>%
      str_remove(paste0("_", PFUNC, "$"))
  )


total_slopes_pf <- rand_long_pf %>%
  left_join(fixed_long_pf, by = c(".draw", "PFUNC")) %>%
  mutate(total_slope = fixed_slope + rand_slope)


total_slope_summary_pf <- total_slopes_pf %>%
  group_by(SiteSubsite, PFUNC) %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )



#summarise not by subsite totalty by PFUNC

total_slope_summary_PFUNC <- total_slopes_pf %>%
  group_by(PFUNC) %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )


total_slope_summary_pf <- total_slope_summary_pf%>%
  mutate(PFUNC = fct_relevel(PFUNC,  
                             "SDECI",
                             "SEVER", 
                             "GRAMMINOID",
                             "FORB"
  ))


custom_palette <- c(
  "SDECI"="dodgerblue4",
  "SEVER"="chartreuse4",
  "GRAMMINOID"="lightpink1",
  "FORB"="slateblue1")


####---------------------This one to combine with the other for figure
(PFUNC_JAN <- ggplot(
  total_slope_summary_pf,
  aes(
    x = PFUNC,
    y = mean_slope,
    fill = PFUNC
  )
) +
    geom_boxplot(
      alpha = 0.5,
      linewidth = 0.3,
      width = 0.7
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    scale_fill_manual(
      values = custom_palette,
      guide = guide_legend(
        override.aes = list(
          shape = NA,
          linewidth = 0,
          alpha = 0.5
        )
      )
    ) +
    theme_clean() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.key.height = unit(0.8, "cm"),
      legend.key.width  = unit(1.2, "cm")
    ) +
    labs(
      x = "Functional group (PFUNC)",
      y = "Mean site-level slope (PropCov per standardized year)",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("PFUNC_JAN.png", plot = PFUNC_JAN, width = 16, height = 6, dpi = 300)

#Another way to visualise these is to count the slopes..
#################################################            SLOPE COUNTS FOR EACH SITESUBSITE AND GMYCPFUNC GROUPS
#add the model output

load("R_Lisa/ITEX/ITEX_data/total_slopes_3.RData")


total_slope_summary_3 <- total_slopes_3 %>%
  group_by(SITE,SiteSubsite, GMYCPFUNC) %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )


#add increase or decerease and add significant columns
total_slope_summary_3 <- total_slope_summary_3 %>%
  mutate(
    direction = case_when(
      mean_slope > 0 ~ "increase",
      mean_slope < 0 ~ "decrease",
      TRUE           ~ "no change"
    ),
    SIG = case_when(
      ci_low > 0 | ci_high < 0 ~ "significant",
      TRUE                    ~ "not significant"
    )
  )

total_slope_summary_3 <- total_slope_summary_3 %>%
  mutate(
    slope = case_when(
      direction == "increase" & SIG == "significant"     ~ "significant increase",
      direction == "increase" & SIG == "not significant" ~ "non-significant increase",
      direction == "decrease" & SIG == "significant"     ~ "significant decrease",
      direction == "decrease" & SIG == "not significant" ~ "non-significant decrease",
      TRUE                                               ~ "no change"
    )
  )



slope_counts <- total_slope_summary_3%>%
  group_by(GMYCPFUNC, slope, SIG) %>%
  summarise(count = n(), .groups = 'drop')


slope_counts <- slope_counts %>%
  mutate(
    total_count = if_else(
      str_detect(slope, "decrease"),
      -count,
      count
    )
  )

#we reorder intothe order we want
slope_counts <- slope_counts %>%
  mutate(GMYCPFUNC = fct_relevel(GMYCPFUNC,  
                                 "EcM:SDECI",
                                 "ErM:SDECI",
                                 "EcMErM:SDECI",
                                 "EcMAM:SDECI", 
                                 "AM:SDECI" ,
                                 "EcM:SEVER", 
                                 "ErM:SEVER",
                                 "EcMAM:SEVER",  
                                 "AM:SEVER",
                                 "EcMAM:GRAMMINOID",
                                 "AM:GRAMMINOID",
                                 "NMAM:GRAMMINOID",
                                 "EcM:FORB",
                                 "EcMAM:FORB",
                                 "AM:FORB",
                                 "NMAM:FORB",
                                 "NM:FORB"
  ))

#set the custonm pallete for out groups
custom_palette <- c(
  "EcM:SDECI"="dodgerblue4",
  "ErM:SDECI"="chartreuse4",
  "EcMErM:SDECI"="purple3",
  "EcMAM:SDECI"="deeppink4", 
  "AM:SDECI"="slateblue1" ,
  "EcM:SEVER"="dodgerblue4", 
  "ErM:SEVER"="chartreuse4",
  "EcMAM:SEVER"="deeppink4",  
  "AM:SEVER"="slateblue1",
  "EcMAM:GRAMMINOID"="deeppink4",
  "AM:GRAMMINOID"="slateblue1",
  "NMAM:GRAMMINOID"="lightpink1",
  "EcM:FORB"="dodgerblue4",
  "EcMAM:FORB"="deeppink4",
  "AM:FORB"="slateblue1",
  "NMAM:FORB"="lightpink1",
  "NM:FORB"="deepskyblue2")



#this plot....

(COUNTS_JAN <- ggplot(
  slope_counts,
  aes(
    x = GMYCPFUNC,
    y = total_count,
    fill = GMYCPFUNC,
    alpha = SIG #this makes the shading to lighter possible
  )
) +
    geom_col(
      width = 0.6, #how wide the bars are
      colour = "black"
    ) +
    scale_fill_manual(values = custom_palette) +
    scale_alpha_manual(
      values = c(
        "significant" = 1, #this is the level of alpoha- the shading
        "not significant" = 0.2
      )
    ) +
    scale_y_continuous(
      breaks = seq(
        -70,
        max(abs(slope_counts$total_count), na.rm = TRUE),
        by = 10
      ),
      labels = function(x) abs(x)
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    theme_clean() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      legend.position = "right",
      panel.border = element_blank(),
      plot.title = element_blank()
    ) +
    labs(
      y = "Number of sites",
      fill = "Mycorrhizal type with Functional Group"
    ) +
    guides(alpha = "none"))


ggsave("COUNTS_JAN.png", plot = COUNTS_JAN, width = 16, height = 6, dpi = 300)


############################################               WE DO SLOPE COUNTS EACH SITESUBSITE AND PFUNC            ##################################################
#add the model output

load("R_Lisa/ITEX/ITEX_data/total_slopes_pf.RData")


total_slope_summary_pf <- total_slopes_pf %>%
  group_by(SITE,SiteSubsite, GMYCPFUNC) %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )


#add increase or decerease and add significant columns
total_slope_summary_pf <- total_slope_summary_pf %>%
  mutate(
    direction = case_when(
      mean_slope > 0 ~ "increase",
      mean_slope < 0 ~ "decrease",
      TRUE           ~ "no change"
    ),
    SIG = case_when(
      ci_low > 0 | ci_high < 0 ~ "significant",
      TRUE                    ~ "not significant"
    )
  )

total_slope_summary_pf <- total_slope_summary_pf %>%
  mutate(
    slope = case_when(
      direction == "increase" & SIG == "significant"     ~ "significant increase",
      direction == "increase" & SIG == "not significant" ~ "non-significant increase",
      direction == "decrease" & SIG == "significant"     ~ "significant decrease",
      direction == "decrease" & SIG == "not significant" ~ "non-significant decrease",
      TRUE                                               ~ "no change"
    )
  )



slope_counts_pf <- total_slope_summary_pf%>%
  group_by(PFUNC, slope, SIG) %>%
  summarise(count = n(), .groups = 'drop')


slope_counts_pf <- slope_counts_pf %>%
  mutate(
    total_count = if_else(
      str_detect(slope, "decrease"),
      -count,
      count
    )
  )



slope_counts_pf <- slope_counts_pf%>%
  mutate(PFUNC = fct_relevel(PFUNC,  
                             "SDECI",
                             "SEVER", 
                             "GRAMMINOID",
                             "FORB"
  ))


custom_palette <- c(
  "SDECI"="dodgerblue4",
  "SEVER"="chartreuse4",
  "GRAMMINOID"="lightpink1",
  "FORB"="slateblue1")


(COUNTS_JAN_pf <- ggplot(
  slope_counts_pf,
  aes(
    x = PFUNC,
    y = total_count,
    fill = PFUNC,
    alpha = SIG #this makes the shading to lighter possible
  )
) +
    geom_col(
      width = 0.6, #how wide the bars are
      colour = "black"
    ) +
    scale_fill_manual(values = custom_palette) +
    scale_alpha_manual(
      values = c(
        "significant" = 1, #this is the level of alpoha- the shading
        "not significant" = 0.2
      )
    ) +
    scale_y_continuous(
      breaks = seq(
        -70,
        max(abs(slope_counts$total_count), na.rm = TRUE),
        by = 10
      ),
      labels = function(x) abs(x)
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    theme_clean() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      legend.position = "right",
      panel.border = element_blank(),
      plot.title = element_blank()
    ) +
    labs(
      y = "Number of sites",
      fill = "Mycorrhizal type with Functional Group"
    ) +
    guides(alpha = "none"))


ggsave("COUNTS_JAN_pf.png", plot = COUNTS_JAN_pf, width = 16, height = 6, dpi = 300)


#########################----                            -------LOOK AT THE GROUPS SEPERTELY--------------------#############################################################


# I made some seperate models where I divided the dataset into the GMYCPFUNC Groups, went through the same process to see if they had the same results... they didn't

###################################                      the group models                       ##########################

###########################################################################################-----------------------AM:FORB--------------------

AM_FORB_data <- MYC_ARCTIC_CTL %>%
  dplyr::filter(GMYCPFUNC == "AM:FORB")

#Maske the model I want just GMYCPFUNC*time, this will give each GMYCPFUNC at each sitesubsite.

Model_AM_FORB <- brm(
  PropCov ~ YEAR_relative_sc+
    (1 + YEAR_relative_sc | SiteSubsite),
  data = AM_FORB_data,
  family = zero_one_inflated_beta(),
  prior = prior(gamma(0.01, 0.01), class = phi),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99),
  file = "Model_AM_FORB_autosave"
)


# Save final object cleanly
write_rds(Model_AM_FORB, "Model_AM_FORB.rds")

###########################################################################################-----------------------ECMDECI--------------------

EcM_SDECI_data <- MYC_ARCTIC_CTL %>%
  dplyr::filter(GMYCPFUNC == "EcM:SDECI")

#Maske the model I want just GMYCPFUNC*time, this will give each GMYCPFUNC at each sitesubsite.

Model_EcM_SDECI<- brm(
  PropCov ~ YEAR_relative_sc+
    (1 + YEAR_relative_sc | SiteSubsite),
  data = EcM_SDECI_data,
  family = zero_one_inflated_beta(),
  prior = prior(gamma(0.01, 0.01), class = phi),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99),
  file = "EcM:SDECI_autosave"
)


# Save final object cleanly
write_rds(Model_EcM_SDECI, "Model_EcM_SDECI.rds")


###########################################################################################-----------------------ErM;SEVER--------------------

ErM_SEVER_data <- MYC_ARCTIC_CTL %>%
  dplyr::filter(GMYCPFUNC == "ErM:SEVER")

#Maske the model I want just GMYCPFUNC*time, this will give each GMYCPFUNC at each sitesubsite.

Model_ErM_SEVER<- brm(
  PropCov ~ YEAR_relative_sc+
    (1 + YEAR_relative_sc | SiteSubsite),
  data = ErM_SEVER_data,
  family = zero_one_inflated_beta(),
  prior = prior(gamma(0.01, 0.01), class = phi),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99),
  file = "ErM_SEVER_autosave"
)


# Save final object cleanly
write_rds(Model_ErM_SEVER, "Model_ErM_SEVER.rds")

###########################################################################################-----------------------AM:GRAMMINOID--------------------

AM_GRAMMINOID_data <- MYC_ARCTIC_CTL %>%
  dplyr::filter(GMYCPFUNC == "AM:GRAMMINOID")

#Maske the model I want just GMYCPFUNC*time, this will give each GMYCPFUNC at each sitesubsite.

Model_AM_GRAMMINOID<- brm(
  PropCov ~ YEAR_relative_sc+
    (1 + YEAR_relative_sc | SiteSubsite),
  data = AM_GRAMMINOID_data,
  family = zero_one_inflated_beta(),
  prior = prior(gamma(0.01, 0.01), class = phi),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99),
  file = "AM_GRAMMINOID_autosave"
)


# Save final object cleanly
write_rds(Model_AM_GRAMMINOID, "Model_AM_GRAMMINOID.rds")


##################################               group  model outputs               ########################################




fitmodAMFORB<- readRDS("R_Lisa/Model_AM_FORB.rds") # we add the {R_Lisa/} for remote desktop

summary(fitmodAMFORB) #this looks good so far- we have AM:FORB for seperate interction

#fixef(fitmod3) # to get more detailed values for estimates
#coef(fitmod3) # if you have group-level effects (hierarchical data)

plot(fitmodAMFORB) # all good

pp_check(fitmodAMFORB, ndraws = 100) # this looks okay

#I need to get the fixed slope for the interaction of year and gmycpfunc, then teh random efect of sitesubsite:Gmycpfunc and together they will give my gmycpfunc over time for each sitesubsite

fitmodAMFORB %>% get_variables() %>% grep("^b_", ., value = TRUE)
#"b_Intercept"        "b_YEAR_relative_sc"




#I need this to see therandom effect interation,
Row.names <- as.data.frame(fitmodAMFORB)
#r_SiteSubsite[.........] #in the brackets site:subsite, year_relative_sc # this means we divide where there is a comma, then divide by underscore.




# 1) get posterior draws for every interaction slope
# spread_draws with matches() pulls many columns; here we extract all parameters starting with b_YEAR_relative_sc:
AM_FORB_YEAR<- names(as_draws_df(fitmodAMFORB)) %>%
  grep("^b_YEAR_relative_sc", ., value = TRUE)

all_fixed_draws_AM_FORB <- fitmodAMFORB %>%
  spread_draws(!!!syms(AM_FORB_YEAR))

all_fixed_draws_AM_FORB$"YEAR" <- all_fixed_draws_AM_FORB$"b_YEAR_relative_sc"


#summarise fixed effects:
AM_fixed_slope_summary <- all_fixed_draws_AM_FORB %>%
  summarise(
    mean_slope = mean(YEAR),
    sd_slope   = sd(YEAR),
    ci_low     = quantile(YEAR, 0.025),
    ci_high    = quantile(YEAR, 0.975)
  )

#2 get the posterior draws for the random effects of each subsite


AM_subsite_terms <- names(as_draws_df(fitmodAMFORB)) %>%
  grep("^r_SiteSubsite\\[.*,YEAR_relative_sc\\]$", 
       ., value = TRUE)

AM_all_rand_draws <- fitmodAMFORB%>%
  spread_draws(!!!syms(AM_subsite_terms))

AM_all_draws <- all_fixed_draws_AM_FORB %>%
  left_join(AM_all_rand_draws, by = ".draw")


#we now add together the fixed effect and the random effect for each site
#5/1/26, trying a diffrernt adding approach


AM_rand_long <- AM_all_rand_draws %>%
  pivot_longer(
    cols = -c(.chain, .iteration, .draw),
    names_to = "rand_term",
    values_to = "rand_slope"
  ) %>%
  mutate(
    # remove Stan syntax
    rand_term = str_remove(rand_term, "^r_SiteSubsite\\["),
    rand_term = str_remove(rand_term, ",YEAR_relative_sc\\]$"),
    
    # split identifiers
    SITE = str_extract(rand_term, "^[^:]+"),
    SiteSubsite = str_remove(rand_term, paste0("^", SITE, ":")) 
  )

AM_total_slopes <- AM_rand_long %>%
  left_join(all_fixed_draws_AM_FORB, by = c(".draw")) %>%
  mutate(total_slope = YEAR + rand_slope)


AM_total_slope_summary <- AM_total_slopes %>%
  group_by(SiteSubsite) %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )



AM_total_summary <- AM_total_slopes %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )





fitmodEcM_SDECI<- readRDS("R_Lisa/Model_EcM_SDECI.rds") # we add the {R_Lisa/} for remote desktop

summary(fitmodEcM_SDECI) #this looks good so far- we have AM:FORB for seperate interction

#fixef(fitmod3) # to get more detailed values for estimates
#coef(fitmod3) # if you have group-level effects (hierarchical data)

plot(fitmodEcM_SDECI) # all good

pp_check(fitmodEcM_SDECI, ndraws = 100) # this looks okay

#I need to get the fixed slope, then the random effect of sitesubsite and together they will give my over time for each sitesubsite

fitmodEcM_SDECI %>% get_variables() %>% grep("^b_", ., value = TRUE)
#"b_Intercept"        "b_YEAR_relative_sc"




#I need this to see therandom effect interation,
Row.names <- as.data.frame(fitmodEcM_SDECI)
#r_SiteSubsite[.........] #in the brackets site:subsite, year_relative_sc # this means we divide where there is a comma, then divide by underscore.




# 1) get posterior draws for every interaction slope
# spread_draws with matches() pulls many columns; here we extract all parameters starting with b_YEAR_relative_sc:
EcM_SDECI_YEAR<- names(as_draws_df(fitmodEcM_SDECI)) %>%
  grep("^b_YEAR_relative_sc", ., value = TRUE)

all_fixed_draws_EcM_SDECI <- fitmodEcM_SDECI %>%
  spread_draws(!!!syms(EcM_SDECI_YEAR))

all_fixed_draws_EcM_SDECI$"YEAR" <- all_fixed_draws_EcM_SDECI$"b_YEAR_relative_sc"


#summarise fixed effects:
EcM_fixed_slope_summary <- all_fixed_draws_EcM_SDECI %>%
  summarise(
    mean_slope = mean(YEAR),
    sd_slope   = sd(YEAR),
    ci_low     = quantile(YEAR, 0.025),
    ci_high    = quantile(YEAR, 0.975)
  )

#2 get the posterior draws for the random effects of each subsite


EcM_subsite_terms <- names(as_draws_df(fitmodEcM_SDECI)) %>%
  grep("^r_SiteSubsite\\[.*,YEAR_relative_sc\\]$", 
       ., value = TRUE)

EcM_all_rand_draws <- fitmodEcM_SDECI%>%
  spread_draws(!!!syms(EcM_subsite_terms))

EcM_all_draws <- all_fixed_draws_AM_FORB %>%
  left_join(AM_all_rand_draws, by = ".draw")


#we now add together the fixed effect and the random effect for each site
#5/1/26, trying a diffrernt adding approach


EcM_rand_long <- EcM_all_rand_draws %>%
  pivot_longer(
    cols = -c(.chain, .iteration, .draw),
    names_to = "rand_term",
    values_to = "rand_slope"
  ) %>%
  mutate(
    # remove Stan syntax
    rand_term = str_remove(rand_term, "^r_SiteSubsite\\["),
    rand_term = str_remove(rand_term, ",YEAR_relative_sc\\]$"),
    
    # split identifiers
    SITE = str_extract(rand_term, "^[^:]+"),
    SiteSubsite = str_remove(rand_term, paste0("^", SITE, ":")) 
  )

EcM_total_slopes <- EcM_rand_long %>%
  left_join(all_fixed_draws_EcM_SDECI, by = c(".draw")) %>%
  mutate(total_slope = YEAR + rand_slope)


EcM_total_slope_summary <- EcM_total_slopes %>%
  group_by(SiteSubsite) %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )



EcM_total_summary <- EcM_total_slopes %>%
  summarise(
    mean_slope = mean(total_slope),
    sd_slope   = sd(total_slope),
    ci_low     = quantile(total_slope, 0.025),
    ci_high    = quantile(total_slope, 0.975),
    .groups = "drop"
  )

#