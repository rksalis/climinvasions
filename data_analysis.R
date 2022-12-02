# data visualisation and analysis of species March 16 2022
# R. K. Salis 
# Purpose: This document is for the plotting and analysis of the 18S data from weeks 0, 10, 12, 14, 18 and 22 of the mesocosm experiment. Bioinformatics and data processing is done in separate scripts.

# Includes:
# 1. Identification of potential and successful invaders
# calculate invasion success and get numbers for venn diagram Fig. 1
# 2. Calculation of invasion success
# 3. Fig. 2
# 4. Richness
# 4.1 plot Fig. 3
# 4.2 Evenness plot Fig. S3
# 4.3 RM-ANOVAs
# 4.4 2-way ANOVASs

#packages
install.packages("data.table")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("rlist")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("tidyr")
install.packages("tidytext")
library(data.table)
library(ggplot2)
library(ggpubr)
library(rlist)
library(tidyverse)
library(dplyr)
library(tidyr)
library(tidytext)

#load data
tax_table <- fread("tax_table_eDNA1_18S_pr2_mesocosms_vs_sp.csv",header=TRUE)
sample_data <- fread("sample_data_eDNA1_18S_pr2_mesocosms_vs_sp.csv")
IW_tax <- fread("tax_table_eDNA1_18Sf_invasion_vs_sp.csv",header=TRUE)
asv_table1 <- fread("asv_table_eDNA1_18S_pr2_mesocosms_vs_sp.csv")
asv_table <- dcast(melt(asv_table1, id.vars = "V1"), variable ~ V1)

# 1. Identification of potential and successful invaders
# invasion dataset
IW <-  IW_tax[, .(Species)]
IW_list <-  IW[, Species]
IW_list #238 species in invasion water dataset compared to 634 species in mesocosm dataset 

dt_ASVsamp <- sample_data[asv_table, on = .(V1 = variable)]
#convert to long format and retain only the variables want to keep
dt_long <- melt(dt_ASVsamp,
     id.vars = c("Sample.date","sampling.week","treatment","Mesocosm","Sample.number","rep","Heated","Invasion"),
     measure.vars = patterns("^ASV"),
     variable.name = "ASV",
     value.name = c("normalised_reads"))

#add taxonomy information
dt_long_sp <- dt_long[tax_table, on = .(ASV = V1)]

#convert back to wide
dt_wide_sp <- dcast(dt_long_sp, Kingdom+Supergroup+Division+Class+Order+Family+Genus+Species+treatment+Mesocosm+rep~Sample.date,
                              value.var =  c("normalised_reads"))
#sum sequence abundance before
dt_wide_sp[ , before :=rowSums(.SD, na.rm = TRUE), .SDcols = c("19/05/20", "28/07/20")]
#sum sequence abundance after
dt_wide_sp[ , after :=rowSums(.SD, na.rm = TRUE), .SDcols = c("11/08/20", "25/08/20", "22/09/20", "20/10/20")]
#means for species and treatment
dt_wide_sp[, mean_before:= mean(before),by = c("treatment","Species")] 
dt_wide_sp[, mean_after:= mean(after),by = c("treatment","Species")]
fwrite(dt_wide_sp, "dt_wide_sp_vs.csv") 

#sum sequence abundance w22
dt_wide_sp[ , w22 :=rowSums(.SD, na.rm = TRUE), .SDcols = c("20/10/20")]
dt_wide_sp[, mean_w22:= mean(w22),by = c("treatment","Species")]
fwrite(dt_wide_sp, "dt_wide_sp_vs_w22.csv") 

Species_mesocosms_before <- dt_wide_sp[mean_before > 0,Species]
Species_mesocosms_before_list <- unique(Species_mesocosms_before, by = Species)
Species_mesocosms_before_list #491 species in mesocosms before invasion
#POTENTIAL INVADERS - What is different in Invasion Water (those in Invasion Water but not in mesocosms before invasion)
potential_invaders <- setdiff(IW_list,Species_mesocosms_before_list)
potential_invaders #96 species
#then of these 96 potential invaders - how many found in the mesocosms afterwards??
setkey(dt_wide_sp,treatment)

dt_wide_sp[c("HI","I")]#invasion treatments (HI= Warming + Invasion, I = Invasion)
Species_mesocosms_after_HIxI <-  dt_wide_sp[c("HI","I")]
Species_mesocosms_after <- Species_mesocosms_after_HIxI[mean_after > 0,Species]
Species_mesocosms_after_list <- unique(Species_mesocosms_after, by = Species)
Species_mesocosms_after_list #461 species in mesocosms after invasion
Species<- intersect(potential_invaders,Species_mesocosms_after_list)
Species #13 of potential invaders present in mesocosms after invasion
dt_successful_invaders <- as.data.table(Species)
dt_successful_invaders[,succ_inv:= "TRUE"]
dt_long_sp_inv <- dt_successful_invaders [dt_long_sp, on = .(Species = Species )]
dt_long_sp_inv[is.na(succ_inv),succ_inv:= FALSE ]

#get numbers for venn diagram
IWafter<- intersect(IW_list,Species_mesocosms_after_list)
IWafter #139
IWbefore<- intersect(IW_list,Species_mesocosms_before_list)
IWbefore #142
IWab<- intersect(IWafter,IWbefore)
IWab #126
IWa_b<- setdiff(IWafter,IWbefore)
IWa_b #13
IWb_a<- setdiff(IWbefore,IWafter)
IWb_a #16
mesos_ba <- intersect(Species_mesocosms_before_list,Species_mesocosms_after_list)
mesos_ba #356
mesos_a_b <- setdiff(Species_mesocosms_after_list,Species_mesocosms_before_list)
mesos_a_b #105
mesos_b_a <- setdiff(Species_mesocosms_before_list,Species_mesocosms_after_list)
mesos_b_a #135

#create supplementary table (Dataset S1) of all species added in invasion
dt_potential_invaders <- as.data.table(potential_invaders)
dt_potential_invaders[,pot_inv:= "TRUE"]
dt_IW <- dt_potential_invaders [IW_tax, on = .(potential_invaders = Species )]
dt_IW_all <- dt_successful_invaders [dt_IW, on = .(Species = potential_invaders )]
#write csv
fwrite(dt_IW_all, "Invasion_allspeciesadded_vs.csv") 
dt_IW_all_inv <- dt_IW_all [dt_long_sp, on = .(Species = Species )]
fwrite(dt_IW_all_inv, "dt_IW_all_inv.csv") 

#subset dt_long_sp_inv
dt_succ_inv <- dt_long_sp_inv[succ_inv == TRUE, ]
dt_succ_inv2 <- dt_succ_inv[sampling.week != "0", ][sampling.week != "10", ] 
dt_succ_inv2_wide <- dcast(dt_succ_inv2, Sample.date+sampling.week+treatment+Heated+Invasion+Mesocosm+rep~Species,
                           value.var =  c("normalised_reads"))
fwrite(dt_succ_inv2_wide, "successful_invaders_normalised_vs.csv") 
dt_succ_inv2_wide2 <- dcast(dt_succ_inv2, treatment+Heated+Invasion+Mesocosm+rep~Species+sampling.week,
                           value.var =  c("normalised_reads"))
fwrite(dt_succ_inv2_wide2, "successful_invaders_normalisedwide_vs.csv")

#subset to retain rows that have normalised_reads > 0
dt_succ_inv3 <- dt_succ_inv2[normalised_reads > 0, ] 

#calculate overall invasion success - number of species out of the potential invaders that are detected in I and HI
dt_succ_invHI <- dt_succ_inv3[treatment=="HI", ] 
dt_succ_invHI2 <- dt_succ_invHI[normalised_reads > 0, ] 
dt_succ_invHI_wide <- dcast(dt_succ_invHI2, Species~treatment,
                           value.var =  c("normalised_reads")) # 12 species
dt_succ_invI <- dt_succ_inv3[treatment=="I", ]
dt_succ_invI2 <- dt_succ_invI[normalised_reads > 0, ] 
dt_succ_invI_wide <- dcast(dt_succ_invI2, Species~treatment,
                                  value.var =  c("normalised_reads")) # 6

# 1.3 calculate invasion success(%)
#number of mesocosms and sampling days a species is present in following the invasion 
# / 4sampling days*6 mesocosms = 24 (for each treatment)
dt_invasion_success <- dt_succ_inv3[, (count = .N), by = c("treatment","Species")]
dt_invasion_success[, invasion_succ:=V1/24*100 ]

dt_invasion_success_wide <- dcast(dt_invasion_success, Species~treatment,
                    value.var =  c("invasion_succ"))
#add taxonomy data
dt_invasion_success_wide <- IW_tax[dt_invasion_success_wide, on = .(Species = Species)]
#write csv
fwrite(dt_invasion_success_wide, "Invasion_success_normalised_vs.csv") 


dt_invasion_success2 <- melt(dt_invasion_success_wide,
                id.vars = c("V1","Kingdom", "Supergroup", "Division","Class","Order","Family","Genus","Species"),
                variable.name = "treatment",
                value.name = c("invasion_succ"))
dt_invasion_success2[is.na(invasion_succ),invasion_succ:= 0 ]


#Fig 2
dt_invasion_success_wide2 <- dt_invasion_success_wide[is.na(H),H:= 0 ]
dt_invasion_success_wide2[is.na(HI),HI:= 0 ]
dt_invasion_success_wide2[is.na(I),I:= 0 ]
dt_invasion_success_wide2[is.na(C),C:= 0 ]
dt_invasion_success_wide2[, HI_I:=HI-I ]
df_invasion_success_wide2 <- as.data.frame(dt_invasion_success_wide2)
df_invasion_success_wide2$Species <- factor(df_invasion_success_wide2$Species)
fct_reorder(df_invasion_success_wide2$Species,df_invasion_success_wide2$HI_I,.desc = TRUE)

dt_invasion_success_wide2 <- dt_invasion_success_wide2[ , Species := as.factor(Species)]  
setorder(dt_invasion_success_wide2,-HI_I)

df_invasion_success_wide3 <- 
  df_invasion_success_wide2 %>%
ungroup() %>% arrange(Supergroup,Division,desc(HI_I)) %>%
  mutate(order = row_number())

ggplot(df_invasion_success_wide3) +
  geom_col(aes(y=HI_I, x=order,group=Division, fill = factor(sign(HI_I))),position = position_dodge(1))+
  theme(legend.position = "none") +  theme_bw() +
  scale_fill_manual(values = c( "#0072B2","black","orange"),breaks=c("-1", "0", "1"),labels = c("I>HI", "HI=I", "HI>I"))+
  theme(axis.text.x = element_text(angle = 90, size = 5, hjust=1, vjust=-.05))+
  ylab("Effect of Climate Warming Scenario on Invasion Success (% Change HI-I)")  +
  xlab("Species")+
  facet_grid(rows=vars(Supergroup,Division),scales = "free_y", space="free",switch="y")+                                                                # Change font size
  theme(strip.text.x = element_text(size = 5))+
  theme(legend.position="none",panel.border = element_blank())+
  scale_x_continuous(breaks = df_invasion_success_wide3$order, labels = df_invasion_success_wide3$Species,expand = c(0,0) ) +
  coord_flip() +
  theme(strip.text.y.left = element_text(angle = 0)) +
  theme(panel.spacing.y=unit(0.1, "lines"))

# 4. Richness
library(ggplot2)
library(plyr)

richness_sp=read.csv("Richness_18S_vs_sp.csv")                      
richness_sp_AFTER=read.csv("Richness_18S_vs_sp_after.csv")
richness_sp_AFTER_w22=read.csv("Richness_18S_vs_sp_w22.csv")

### Create summary function
#Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).  data: a data frame.  measurevar: the name of a column that contains the variable to be summarized.  groupvars: a vector containing names of columns that contain grouping variables.  na.rm: a boolean that indicates whether to ignore NA's.  conf.interval: the percent range of the confidence interval (default is 95%)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
} 

# 4.1 plot Fig. 3
richness_sp_tot <- summarySE(richness_sp, measurevar="richness", groupvars=c("sampling.week","treatment"))

ggplot(richness_sp_tot, aes(x=sampling.week, y=richness, group=treatment)) + 
  geom_errorbar(aes(ymin=richness-se, ymax=richness+se), colour="dark grey", width=.1, position=position_dodge(0.1)) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(aes(shape=treatment), size = 3, position = position_dodge(0.1), fill="white") +
  scale_shape_manual(values = c(21, 22, 15, 19)) +
  xlab("Sampling date") +
  ylab("Species richness") +
  theme_bw()  + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

# 4.2 Evenness plot Fig. S3
evenness_sp_tot <- summarySE(richness_sp, measurevar="evenness", groupvars=c("sampling.week","treatment"))

ggplot(evenness_sp_tot, aes(x=sampling.week, y=evenness, group=treatment)) + 
  geom_errorbar(aes(ymin=evenness-se, ymax=evenness+se), colour="dark grey", width=.1, position=position_dodge(0.1)) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(aes(shape=treatment), size = 3, position = position_dodge(0.1), fill="white") +
  scale_shape_manual(values = c(21, 22, 15, 19)) +
  xlab("Sampling date") +
  ylab("Species evenness") +
  theme_bw()  + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 


# Statistics
library(rstatix)
library(heplots)

# 4.3 RM-ANOVAs
richness_sp_AFTER$sampling.week <- as.numeric(as.character(richness_sp_AFTER$sampling.week))

## Richness
M_richRM <- anova_test(
  data = richness_sp_AFTER, dv = richness, wid = Mesocosm,
  between = c(Heated,Invasion),
  within = c(sampling.week)
)
M_richRM
get_anova_table(M_richRM)

## Evenness
M_evenRM <- anova_test(
  data = richness_sp_AFTER, dv = evenness, wid = Mesocosm,
  between = c(Heated,Invasion),
  within = c(sampling.week)
)
M_evenRM
get_anova_table(M_evenRM)

# 4.4 2-way ANOVASs 
## Richness - final sampling
M_rich22 <- aov(richness ~ Heated*Invasion, data = richness_sp_AFTER_w22)
summary(M_rich22)
etasq(M_rich22)
## Evenness - final sampling
M_even22 <- aov(evenness ~ Heated*Invasion, data = richness_sp_AFTER_w22)
summary(M_even22)
etasq(M_even22)
