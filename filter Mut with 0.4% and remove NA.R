library(readxl)
library(readr)
library(dplyr)
library(tidyverse)
##readi dataframe##
df1 <- read_xlsx("/Users/chendd/LowFreMut/screening/baseline error rate/raw data copy 2.xlsx", sheet = "abs fre plasmid (lib1-3)")

###assign the library region based on amino acid position. Firstly parse amino acid position by scripts below:
parse_number(df1$Frame1)
df1$Lib_region <- ifelse(parse_number(df1$Frame1)< 81,  "lib1",
                         ifelse(parse_number(df1$Frame1) >80 & parse_number(df1$Frame1) < 161, "lib2", "lib3" ))
### group data by library region
df2 <- group_by(df1, Lib_region) 

## replace zero with NA
df2[df2 == 0] <- NA

##filter out rows if mutant frequency in WT = NA. 
WT_statistic <- df2[!is.na(df2$`lib1-9-wt abs fre`), ]
##count the mutants in WT plasmid
WT_total <- as.data.frame(table(WT_statistic$Lib_region)) %>% rename("WT_total" = "Freq")
##filter out rows if mutant frequency in mut library = NA.
Mut_statistics <- df2[!is.na(df2$`Lib1-9-mut-backup abs fre`), ]
##count the mutants in mut library.
Mut_total <- as.data.frame(table(Mut_statistics$Lib_region)) %>% rename("Mut_total" = "Freq")


##provide statistics after remove NA values in WT plasmid or Mut library.
Mut_Mean <- summarise(Mut_statistics, Mut_mean = mean(`Lib1-9-mut-backup abs fre`), Mut_std = sd(`Lib1-9-mut-backup abs fre`))
WT_Mean <- summarise(WT_statistic, WT_mean = mean(`lib1-9-wt abs fre`),WT_std = sd(`lib1-9-wt abs fre`))

###################  filter all mutant with 0.04% threshold   #################
##########################################################################
Mut_statistics_post_filter <- filter(df2, `Lib1-9-mut-backup abs fre` > 0.0004) ##filter Mut via 0.0004 threshold, thus no need to filter rows with NA values.
##count the mutants in mut library.
Mut_C_postfilter <- as.data.frame(table(Mut_statistics_post_filter$Lib_region)) %>% rename("Mut_C_postfilter" = "Freq")
##provide statistics after remove NA values in WT plasmid or Mut library.
Mut_Mean_postfilter <- summarise(Mut_statistics_post_filter, Mut_mean_postfilter = mean(`Lib1-9-mut-backup abs fre`),Mut_std_postfilter = sd(`Lib1-9-mut-backup abs fre`))


WT_statistic_post_filter <- filter(df2, `lib1-9-wt abs fre` > 0.0004) ##filter WT via 0.0004 threshold.
##count the mutants in WT plasmid.
WT_C_postfilter <- as.data.frame(table(WT_statistic_post_filter$Lib_region)) %>% rename ("WT_C_postfilter" = "Freq")
WT_Mean_postfilter <- summarise(WT_statistic_post_filter, WT_mean_postfilter = mean(`lib1-9-wt abs fre`), WT_std_postfilter = sd(`lib1-9-wt abs fre`))


list(WT_total, Mut_total, WT_C_postfilter, Mut_C_postfilter) %>% reduce (full_join, by = "Var1") -> summary_counts

list(WT_Mean, Mut_Mean, WT_Mean_postfilter, Mut_Mean_postfilter) %>% reduce (full_join, by = "Lib_region")  -> summary_mean





