#Changes in PLant cover in the ITEX network- the role of mycorhizal type
#Lisa Pilkinton
#29/01/2026- making species a random effect
#Update and overview 17/02:

#Output from the full model that includes the random slope and intercept for species, 
#model: SPECIES_RANDOM_SLOPE <- in Script Models_Slope_and_Intercept
#Figure for GMYCPFUNC: RAN_GMYC_Lisa_JER_SLOPE
#Count figure: RAN_COUNTS_JAN_SLOPE

#Output from the models that were PFUNC and GMYC separatley
#
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


############################################################################################################################30/1 SPECIES as random effect

fitmod_RAN <- readRDS("New_R_Lisa/SPECIES_RANDOM_SLOPE.rds") # we add the {R_Lisa/} for remote desktop <- Model originally in SPECIES_RANDOM_ITEX

summary(fitmod_RAN) #this looks good so far

#catterpillar checks
plot(fitmod_RAN) # all good

pp_check(fitmod_RAN, ndraws = 100) # this looks okay

#following Jerme'y code for getting prosteriors
library(data.table)
MYC_ARCTIC_CTL_dt <- data.table(MYC_ARCTIC_CTL)

#making a dummy dataframe for teh output we want
RAN_dummy_dt <- MYC_ARCTIC_CTL_dt[,.(YEAR_relative_sc = range(YEAR_relative_sc)),by = .(SiteSubsite,GMYCPFUNC,GMYC,PFUNC)] #<-add sp.nm
RAN_dummy_dt[,diff_year := diff(YEAR_relative_sc),by = .(SiteSubsite,GMYCPFUNC)] #<- add sp. nm

RAN_preds <- fitted(fitmod_RAN,
                    newdata = RAN_dummy_dt,
                    re_formula  =  ~  (1 + YEAR_relative_sc | SiteSubsite:GMYCPFUNC) , #+ (1+ YEAR_relative_sc|SPECIES_NAME:GMYCPFUNC) , <- this taken out so its accounted for but when i want species slopes to present i put ot back in
                    seed = 0,allow_new_levels = F,summary = F)

RAN_preds <- RAN_preds*100 # just to get a cover 100
RAN_preds <- t(RAN_preds)
RAN_preds <- cbind(RAN_dummy_dt,RAN_preds)


RAN_preds <- data.table(melt(RAN_preds,id.vars = c("SiteSubsite","GMYCPFUNC","YEAR_relative_sc","GMYC","PFUNC","diff_year")))  # prob_clim_dist <- add species name for species slopes


RAN_preds_change <- RAN_preds[,.(change = diff(value)/ unique(diff_year)  ),by = .(SiteSubsite,GMYCPFUNC,GMYC,PFUNC,variable)] #<- add species_name for sp,. slopes


RAN_preds_change_summary <- RAN_preds_change[,.(mean = mean(change),
                                                Q05 = quantile(change,probs = c(0.025)),
                                                Q95 = quantile(change,probs = c(0.975))),by = .(SiteSubsite,GMYCPFUNC,GMYC,PFUNC)][order(-mean),] #<- add sp.name for sp. slopes
RAN_preds_change_summary[,signif := ifelse(sign(Q05) == sign(Q95) & sign(mean) == -1,"Decline",
                                           ifelse(sign(Q05) == sign(Q95) & sign(mean) == 1,"Increase","Non signif"))]
table(RAN_preds_change_summary$signif,RAN_preds_change_summary$GMYCPFUNC)
table(RAN_preds_change_summary$mean>0)

ggplot(preds_change_summary,aes( x = reorder(GMYC,PFUNC), y = mean,fill = PFUNC , group = GMYCPFUNC))+
  geom_boxplot()+theme_bw()

reduced_MYC_ARCTIC_CTL <- unique(
  MYC_ARCTIC_CTL[, c("SiteSubsite", "AZONE", "AlpArc", "Region", "LAT")]
)

# Merge with all teh other data about the sites
RAN_preds_change_summary <- left_join(RAN_preds_change_summary, reduced_MYC_ARCTIC_CTL, by = "SiteSubsite")
#now we use my code to make the figures I want. 

#First the global fig..
RAN_preds_change_summary<- RAN_preds_change_summary %>%
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
(RAN_GMYC_Lisa_JER_SLOPE <- ggplot(
  RAN_preds_change_summary,
  aes(
    x = GMYCPFUNC,
    y = mean,
    fill = GMYCPFUNC
  )
) +
    geom_boxplot(
      alpha = 1,
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
          alpha = 1
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
      y = "Mean site-level slope (Change in Percent of cover per year)",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("RAN_GMYC_Lisa_JER_SLOPE.png", plot = RAN_GMYC_Lisa_JER_SLOPE, width = 16, height = 6, dpi = 300)

#############################################################################################################################COUNT FIGURE INC RANDOM EFFECT

#make count figure

RAN_slope_summary <- RAN_preds_change_summary %>%
  mutate(
    direction = case_when(
      mean > 0 ~ "increase",
      mean < 0 ~ "decrease",
      TRUE           ~ "no change"
    ),
    SIG = case_when(
      Q05 > 0 | Q95 < 0 ~ "significant",
      TRUE                    ~ "not significant"
    )
  )

RAN_slope_summary <- RAN_slope_summary %>%
  mutate(
    slope = case_when(
      direction == "increase" & SIG == "significant"     ~ "significant increase",
      direction == "increase" & SIG == "not significant" ~ "non-significant increase",
      direction == "decrease" & SIG == "significant"     ~ "significant decrease",
      direction == "decrease" & SIG == "not significant" ~ "non-significant decrease",
      TRUE                                               ~ "no change"
    )
  )



RAN_slope_counts <- RAN_slope_summary%>%
  group_by(GMYCPFUNC, slope, SIG) %>%
  summarise(count = n(), .groups = 'drop')


RAN_slope_counts <- RAN_slope_counts %>%
  mutate(
    total_count = if_else(
      str_detect(slope, "decrease"),
      -count,
      count
    )
  )

#we reorder intothe order we want
RAN_slope_counts <- RAN_slope_counts %>%
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




#Use this plot....

(RAN_COUNTS_JAN_SLOPE <- ggplot(
  RAN_slope_counts,
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
        max(abs(RAN_slope_counts$total_count), na.rm = TRUE),
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


ggsave("RAN_COUNTS_JAN_SLOPE.png", plot = RAN_COUNTS_JAN_SLOPE, width = 16, height = 6, dpi = 300)


#####################################################################################################################FINISH FIG 1 &2


####---------------make the plot without ALPINE sites

RAN_preds_change_summary_HL <- RAN_preds_change_summary[
  RAN_preds_change_summary$AZONE %in% c("HIGH", "LOW"),
]

(RAN_GMYC_Lisa_JER_ALP <- ggplot(
  RAN_preds_change_summary_HL,
  aes(
    x = GMYCPFUNC,
    y = mean,
    fill = GMYCPFUNC
  )
) +
    geom_boxplot(
      alpha = 1,
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
          alpha = 1
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
      y = "Mean site-level slope (Change in Percent of cover per year)",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("RAN_GMYC_Lisa_JER_ALP.png", plot = RAN_GMYC_Lisa_JER_ALP, width = 16, height = 6, dpi = 300)

####---------------make the plot
(RAN_FACET_AZONE <- ggplot(
  RAN_preds_change_summary,
  aes(
    x = GMYCPFUNC,
    y = mean,
    fill = GMYCPFUNC
  )
) +
    geom_boxplot(
      alpha = 1,
      linewidth = 0.3,
      width = 0.7
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    facet_wrap(~AZONE, ncol = 1) +
    scale_fill_manual(
      values = custom_palette,
      guide = guide_legend(
        override.aes = list(
          shape = NA,
          linewidth = 0,
          alpha = 1
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
      y = "Mean site-level slope (Change in Percent of cover per year)",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("RAN_FACET_AZONE.png", plot = RAN_FACET_AZONE, width = 16, height = 18, dpi = 300)

#############FACET Moisture

####---------------make the plot
(RAN_FACET_MOIST <- ggplot(
  RAN_preds_change_summary,
  aes(
    x = GMYCPFUNC,
    y = mean,
    fill = GMYCPFUNC
  )
) +
    geom_boxplot(
      alpha = 1,
      linewidth = 0.3,
      width = 0.7
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    facet_wrap(~MOISTURE, ncol = 1) +
    scale_fill_manual(
      values = custom_palette,
      guide = guide_legend(
        override.aes = list(
          shape = NA,
          linewidth = 0,
          alpha = 1
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
      y = "Mean site-level slope (Change in Percent of cover per year)",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("RAN_FACET_MOIST.png", plot = RAN_FACET_MOIST, width = 16, height = 30, dpi = 300)




############################################################################################################################30/1 SPECIES as random effect PFUNC

PFUNC_fitmod_RAN <- readRDS("New_R_Lisa/PFUNC_SP_RAN_SL_INT.rds") # we add the {R_Lisa/} for remote desktop <- Model originally in SPECIES_RANDOM_ITEX

summary(PFUNC_fitmod_RAN) #this looks good so far

#catterpillar checks
plot(PFUNC_fitmod_RAN) # all good

pp_check(PFUNC_fitmod_RAN, ndraws = 100) # this looks okay

#following Jerme'y code for getting prosteriors
library(data.table)
MYC_ARCTIC_CTL_dt <- data.table(MYC_ARCTIC_CTL)

#making a dummy dataframe for teh output we want
PFUNC_RAN_dummy_dt <- MYC_ARCTIC_CTL_dt[,.(YEAR_relative_sc = range(YEAR_relative_sc)),by = .(SiteSubsite,PFUNC)] #<-add sp.nm
PFUNC_RAN_dummy_dt[,diff_year := diff(YEAR_relative_sc),by = .(SiteSubsite,PFUNC)] #<- add sp. nm

PFUNC_RAN_preds <- fitted(PFUNC_fitmod_RAN,
                    newdata = PFUNC_RAN_dummy_dt,
                    re_formula  =  ~  (1 + YEAR_relative_sc | SiteSubsite:PFUNC) , #+ (1+ YEAR_relative_sc|SPECIES_NAME:GMYCPFUNC) , <- this taken out so its accounted for but when i want species slopes to present i put ot back in
                    seed = 0,allow_new_levels = F,summary = F)

PFUNC_RAN_preds <- PFUNC_RAN_preds*100 # just to get a cover 100
PFUNC_RAN_preds <- t(PFUNC_RAN_preds)
PFUNC_RAN_preds <- cbind(PFUNC_RAN_dummy_dt,PFUNC_RAN_preds)


PFUNC_RAN_preds <- data.table(melt(PFUNC_RAN_preds,id.vars = c("SiteSubsite","PFUNC","YEAR_relative_sc","diff_year")))  # prob_clim_dist <- add species name for species slopes


PFUNC_RAN_preds_change <- PFUNC_RAN_preds[,.(change = diff(value)/ unique(diff_year)  ),by = .(SiteSubsite,PFUNC,variable)] #<- add species_name for sp,. slopes


PFUNC_RAN_preds_change_summary <- PFUNC_RAN_preds_change[,.(mean = mean(change),
                                                Q05 = quantile(change,probs = c(0.025)),
                                                Q95 = quantile(change,probs = c(0.975))),by = .(SiteSubsite,PFUNC)][order(-mean),] #<- add sp.name for sp. slopes
PFUNC_RAN_preds_change_summary[,signif := ifelse(sign(Q05) == sign(Q95) & sign(mean) == -1,"Decline",
                                           ifelse(sign(Q05) == sign(Q95) & sign(mean) == 1,"Increase","Non signif"))]
table(PFUNC_RAN_preds_change_summary$signif,PFUNC_RAN_preds_change_summary$PFUNC)
table(RAN_preds_change_summary$mean>0)

reduced_MYC_ARCTIC_CTL <- unique(
  MYC_ARCTIC_CTL[, c("SiteSubsite", "AZONE", "AlpArc", "Region")]
)

# Merge with all teh other data about the sites
PFUNC_RAN_preds_change_summary <- left_join(PFUNC_RAN_preds_change_summary, reduced_MYC_ARCTIC_CTL, by = "SiteSubsite")
#now we use my code to make the figures I want. 



#PFUNC fig..
PFUNC_RAN_preds_change_summary<- PFUNC_RAN_preds_change_summary %>%
  mutate(PFUNC = fct_relevel(PFUNC,  
                             "SDECI",
                             "SEVER", 
                             "GRAMMINOID",
                             "FORB",
  ))


custom_palette <- c(
  "SDECI"="blue4",
  "SEVER"="green",
  "GRAMMINOID"="red",
  "FORB"="purple")




####---------------make the plot
(PFUNC_6_2 <- ggplot(
  PFUNC_RAN_preds_change_summary,
  aes(
    x = PFUNC,
    y = mean,
    fill = PFUNC
  )
) +
    geom_boxplot(
      alpha = 1,
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
          alpha = 1
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
      y = "Mean change in percent of cover per year",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("PFUNC_6_2.png", plot = PFUNC_6_2, width = 16, height = 6, dpi = 300)


#merge teh two data summaries so we can make a plt with the sisde by side comaparison- call PFUC GMYC in new frame.. merge colour plateea nd prders..

MERGE_PFUNC_RAN_preds_change_summary <- PFUNC_RAN_preds_change_summary

MERGE_PFUNC_RAN_preds_change_summary <- 
PFUNC_RAN_preds_change_summary %>%
  rename(GMYCPFUNC = PFUNC)


# Merge with all teh other data about the sites
MIX_RAN_preds_change_summary <- bind_rows(
  RAN_preds_change_summary,
  MERGE_PFUNC_RAN_preds_change_summary
)

#we reorder intothe order we want
MIX_RAN_preds_change_summary <- MIX_RAN_preds_change_summary  %>%
  mutate(GMYCPFUNC = fct_relevel(GMYCPFUNC,  
                                 "EcM:SDECI",
                                 "ErM:SDECI",
                                 "EcMErM:SDECI",
                                 "EcMAM:SDECI", 
                                 "AM:SDECI" ,
                                 "SDECI",
                                 "EcM:SEVER", 
                                 "ErM:SEVER",
                                 "EcMAM:SEVER",  
                                 "AM:SEVER",
                                 "SEVER", 
                                 "EcMAM:GRAMMINOID",
                                 "AM:GRAMMINOID",
                                 "NMAM:GRAMMINOID",
                                 "GRAMMINOID",
                                 "EcM:FORB",
                                 "EcMAM:FORB",
                                 "AM:FORB",
                                 "NMAM:FORB",
                                 "NM:FORB",
                                 "FORB",
  ))

#set the custonm pallete for out groups
MIX_custom_palette <- c(
  "EcM:SDECI"="dodgerblue4",
  "ErM:SDECI"="chartreuse4",
  "EcMErM:SDECI"="purple3",
  "EcMAM:SDECI"="deeppink4", 
  "SDECI"="blue4",
  "AM:SDECI"="slateblue1" ,
  "EcM:SEVER"="dodgerblue4", 
  "ErM:SEVER"="chartreuse4",
  "EcMAM:SEVER"="deeppink4",  
  "AM:SEVER"="slateblue1",
  "SEVER"="green",
  "EcMAM:GRAMMINOID"="deeppink4",
  "AM:GRAMMINOID"="slateblue1",
  "NMAM:GRAMMINOID"="lightpink1",
  "GRAMMINOID"="red",
  "EcM:FORB"="dodgerblue4",
  "EcMAM:FORB"="deeppink4",
  "AM:FORB"="slateblue1",
  "NMAM:FORB"="lightpink1",
  "NM:FORB"="deepskyblue2",
  "FORB"="purple")



#Mixed PLot

(MIX_PLOT <- ggplot(
  MIX_RAN_preds_change_summary,
  aes(
    x = GMYCPFUNC,
    y = mean,
    fill = GMYCPFUNC
  )
) +
    geom_boxplot(
      alpha = 1,
      linewidth = 0.3,
      width = 0.7
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    scale_fill_manual(
      values = MIX_custom_palette,
      guide = guide_legend(
        override.aes = list(
          shape = NA,
          linewidth = 0,
          alpha = 1
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
      y = "Mean change in percent of cover per year",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("MIX_PLOT.png", plot = MIX_PLOT, width = 20, height = 6, dpi = 300)


#make count fig for PFUNC



PFUNC_COUNT_summary <- PFUNC_RAN_preds_change_summary %>%
  mutate(
    direction = case_when(
      mean > 0 ~ "increase",
      mean < 0 ~ "decrease",
      TRUE           ~ "no change"
    ),
    SIG = case_when(
      Q05 > 0 | Q95 < 0 ~ "significant",
      TRUE                    ~ "not significant"
    )
  )

PFUNC_COUNT_summary <- PFUNC_COUNT_summary%>%
  mutate(
    slope = case_when(
      direction == "increase" & SIG == "significant"     ~ "significant increase",
      direction == "increase" & SIG == "not significant" ~ "non-significant increase",
      direction == "decrease" & SIG == "significant"     ~ "significant decrease",
      direction == "decrease" & SIG == "not significant" ~ "non-significant decrease",
      TRUE                                               ~ "no change"
    )
  )



PFUNC_slope_counts <- PFUNC_COUNT_summary%>%
  group_by(PFUNC, slope, SIG) %>%
  summarise(count = n(), .groups = 'drop')


PFUNC_slope_counts <- PFUNC_slope_counts %>%
  mutate(
    total_count = if_else(
      str_detect(slope, "decrease"),
      -count,
      count
    )
  )
#Use this plot....

(PFUNC_FEB_SLOPE <- ggplot(
  PFUNC_slope_counts,
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
        max(abs(PFUNC_slope_counts$total_count), na.rm = TRUE),
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


ggsave("PFUNC_FEB_SLOPE.png", plot = PFUNC_FEB_SLOPE, width = 16, height = 6, dpi = 300)


############################################################################################################################30/1 SPECIES as random effect GMYC

GMYC_fitmod_RAN <- readRDS("New_R_Lisa/GMYC_SP_RAN_SL_INT.rds") # we add the {R_Lisa/} for remote desktop <- Model originally in SPECIES_RANDOM_ITEX

summary(GMYC_fitmod_RAN) #this looks good so far

#catterpillar checks
plot(GMYC_fitmod_RAN) # all good

pp_check(GMYC_fitmod_RAN, ndraws = 100) # this looks okay

#following Jerme'y code for getting prosteriors
library(data.table)
MYC_ARCTIC_CTL_dt <- data.table(MYC_ARCTIC_CTL)

#making a dummy dataframe for teh output we want
GMYC_RAN_dummy_dt <- MYC_ARCTIC_CTL_dt[,.(YEAR_relative_sc = range(YEAR_relative_sc)),by = .(SiteSubsite,GMYC)] #<-add sp.nm
GMYC_RAN_dummy_dt[,diff_year := diff(YEAR_relative_sc),by = .(SiteSubsite,GMYC)] #<- add sp. nm

GMYC_RAN_preds <- fitted(GMYC_fitmod_RAN,
                          newdata = GMYC_RAN_dummy_dt,
                          re_formula  =  ~  (1 + YEAR_relative_sc | SiteSubsite:GMYC) , #+ (1+ YEAR_relative_sc|SPECIES_NAME:GMYCPFUNC) , <- this taken out so its accounted for but when i want species slopes to present i put ot back in
                          seed = 0,allow_new_levels = F,summary = F)

GMYC_RAN_preds <- GMYC_RAN_preds*100 # just to get a cover 100
GMYC_RAN_preds <- t(GMYC_RAN_preds)
GMYC_RAN_preds <- cbind(GMYC_RAN_dummy_dt,GMYC_RAN_preds)


GMYC_RAN_preds <- data.table(melt(GMYC_RAN_preds,id.vars = c("SiteSubsite","GMYC","YEAR_relative_sc","diff_year")))  # prob_clim_dist <- add species name for species slopes


GMYC_RAN_preds_change <- GMYC_RAN_preds[,.(change = diff(value)/ unique(diff_year)  ),by = .(SiteSubsite,GMYC,variable)] #<- add species_name for sp,. slopes


GMYC_RAN_preds_change_summary <- GMYC_RAN_preds_change[,.(mean = mean(change),
                                                            Q05 = quantile(change,probs = c(0.025)),
                                                            Q95 = quantile(change,probs = c(0.975))),by = .(SiteSubsite,GMYC)][order(-mean),] #<- add sp.name for sp. slopes
GMYC_RAN_preds_change_summary[,signif := ifelse(sign(Q05) == sign(Q95) & sign(mean) == -1,"Decline",
                                                 ifelse(sign(Q05) == sign(Q95) & sign(mean) == 1,"Increase","Non signif"))]
table(GMYC_RAN_preds_change_summary$signif,GMYC_RAN_preds_change_summary$GMYC)
table(RAN_preds_change_summary$mean>0)

reduced_MYC_ARCTIC_CTL <- unique(
  MYC_ARCTIC_CTL[, c("SiteSubsite", "AZONE", "AlpArc", "Region")]
)

# Merge with all teh other data about the sites
GMYC_RAN_preds_change_summary <- left_join(GMYC_RAN_preds_change_summary, reduced_MYC_ARCTIC_CTL, by = "SiteSubsite")
#now we use my code to make the figures I want. 

GMYC_RAN_preds_change_summary <- GMYC_RAN_preds_change_summary %>%
  mutate(GMYC= fct_relevel(GMYC,  
                           "EcM",
                           "ErM",
                           "EcMErM",
                           "EcMAM", 
                           "AM" ,
                           "NMAM",
                           "NM"
  ))

#set the custonm pallete for out groups
custom_palette_3 <- c(
  "EcM"="dodgerblue4",
  "ErM"="chartreuse4",
  "EcMErM"="purple3",
  "EcMAM"="deeppink4", 
  "AM"="slateblue1" ,
  "NMAM"="lightpink1",
  "NM"="deepskyblue2")



####---------------make the plot
(GMYC_6_2 <- ggplot(
  GMYC_RAN_preds_change_summary,
  aes(
    x = GMYC,
    y = mean,
    fill = GMYC
  )
) +
    geom_boxplot(
      alpha = 1,
      linewidth = 0.3,
      width = 0.7
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    scale_fill_manual(
      values = custom_palette_3,
      guide = guide_legend(
        override.aes = list(
          shape = NA,
          linewidth = 0,
          alpha = 1
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
      x = "Functional group (GMYC)",
      y = "Mean change in percent of cover per year",
      title = "Variation in temporal trends among sites",
      subtitle = "Each box summarizes site-level slopes"
    ))

ggsave("GMYC_6_2.png", plot = GMYC_6_2, width = 16, height = 6, dpi = 300)


#make count fig for GMYC



GMYC_COUNT_summary <- GMYC_RAN_preds_change_summary %>%
  mutate(
    direction = case_when(
      mean > 0 ~ "increase",
      mean < 0 ~ "decrease",
      TRUE           ~ "no change"
    ),
    SIG = case_when(
      Q05 > 0 | Q95 < 0 ~ "significant",
      TRUE                    ~ "not significant"
    )
  )

GMYC_COUNT_summary <- GMYC_COUNT_summary%>%
  mutate(
    slope = case_when(
      direction == "increase" & SIG == "significant"     ~ "significant increase",
      direction == "increase" & SIG == "not significant" ~ "non-significant increase",
      direction == "decrease" & SIG == "significant"     ~ "significant decrease",
      direction == "decrease" & SIG == "not significant" ~ "non-significant decrease",
      TRUE                                               ~ "no change"
    )
  )



GMYC_slope_counts <- GMYC_COUNT_summary%>%
  group_by(GMYC, slope, SIG) %>%
  summarise(count = n(), .groups = 'drop')


GMYC_slope_counts <- GMYC_slope_counts %>%
  mutate(
    total_count = if_else(
      str_detect(slope, "decrease"),
      -count,
      count
    )
  )
#Use this plot....

(GMYC_FEB_SLOPE <- ggplot(
  GMYC_slope_counts,
  aes(
    x = GMYC,
    y = total_count,
    fill = GMYC,
    alpha = SIG #this makes the shading to lighter possible
  )
) +
    geom_col(
      width = 0.6, #how wide the bars are
      colour = "black"
    ) +
    scale_fill_manual(values = custom_palette_3) +
    scale_alpha_manual(
      values = c(
        "significant" = 1, #this is the level of alpoha- the shading
        "not significant" = 0.2
      )
    ) +
    scale_y_continuous(
      breaks = seq(
        -70,
        max(abs(GMYC_slope_counts$total_count), na.rm = TRUE),
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
      fill = "Mycorrhizal type"
    ) +
    guides(alpha = "none"))


ggsave("GMYC_FEB_SLOPE.png", plot = GMYC_FEB_SLOPE, width = 16, height = 6, dpi = 300)

