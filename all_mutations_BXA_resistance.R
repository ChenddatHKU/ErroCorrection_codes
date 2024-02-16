library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(GGally)



####reanalyzed all mutations without fitering.
data1 <- read_xlsx("/Users/chendd/Desktop/P1_reanalysized_Feb152024.xlsx", sheet = "Sheet1")
colnames(data1)

points <- c("E23G", "A37T", "I38V", "I38T", "I38M", "L106R", "N136D", "E181D", "E198D", "N228I", "N239Y")
data1 %>% ggplot(aes(x = log10(`P1_withoutBXA_meanRF`), y = log10(`P1_withBXA_meanRF`), color = ifelse(`AA` %in% points, "blue", "black"))) + 
  geom_hline(yintercept=0, linetype ="dashed", color= "black")+ 
  geom_vline(xintercept=0, linetype ="dashed", color= "black") +
  geom_point(size =0.8, alpha=0.6 )+ 
  geom_text_repel(aes(label=ifelse( `AA` %in% points,`AA`,'')), color ="blue", size = 3, fontface='bold', hjust=-0.1, vjust=0.4, max.overlaps = 500)+
  scale_color_identity()+
  theme_classic()+
  ggtitle("all mutations") + 
  scale_x_continuous(limits= c(-2, 2)) + 
  scale_y_continuous(limits = c(-2,2)) +xlab(expression(bold("log"[10]*"RF(0 nM BXA)"))) + ylab(expression(bold("log"[10]*"RF(0.2 nM BXA)")))->p1
p1

####HCMuts###
data2 <- read_xlsx("/Users/chendd/Library/CloudStorage/OneDrive-Personal/Documents/SCREENING DATA analysis/fitness calculation_for plotting.xlsx", sheet = "HcMut")
colnames(data2)

data2 %>% ggplot(aes(y = log10((`P1W-1-fit` + `P1W-2-fit`)/2), x = log10((`P1NO-1-fit` + `P1NO-2-fit`)/2), color = ifelse(`AA` %in% points, "blue", "black"))) + 
  geom_point(size =0.8, alpha = 0.6)+ 
  geom_text_repel(aes(label=ifelse( `AA` %in% points,`AA`,'')), color ="blue", size = 3, fontface='bold',  hjust=-0.1, vjust=0.4, max.overlaps = 500)+
  geom_hline(yintercept=0, linetype ="dashed", color= "black")+ 
  geom_vline(xintercept=0, linetype ="dashed", color= "black")+
  scale_color_identity()+
  theme_classic()+
  ggtitle("fillted mutations")+ 
  scale_x_continuous(limits= c(-2, 2)) + 
  scale_y_continuous(limits = c(-2,2))+xlab(expression(bold("log"[10]*"RF(0 nM BXA)"))) + ylab(expression(bold("log"[10]*"RF(0.2 nM BXA)")))->p2
p2


ggsave("Fig_S1_all.png", plot = p1, width = 5, height = 5, dpi = 300)
ggsave("Fig_S1_filltered.png", plot = p2, width = 5, height = 5, dpi = 300)




##HCMut VS minigenome, Correlation###3
data2 %>% filter(`AA` %in% points) ->data3
data3$average_RF_P1_0nMBXA_ampliconSeq <- (data3$`P1NO-2-fit` + data3$`P2NO-1-fit`)/2
data3$average_RF_P1_0.2nMBXA_ampliconSeq <- (data3$`P1W-1-fit` + data3$`P1W-2-fit`)/2
data3_mean <- select(data3, c("AA", "average_RF_P1_0nMBXA_ampliconSeq", "average_RF_P1_0.2nMBXA_ampliconSeq"))


minigenome <- read_xlsx("/Users/chendd/Desktop/minigenome_data.xlsx", sheet = "Sheet1")
minigenome_mean <- select(minigenome, c("AA", "mean(0nM_BXA)_minigenome", "mean(2nM_BXA)_minigenome"))

CorrDF <- merge(data3_mean, minigenome_mean, by = "AA")

ggpairs(CorrDF, columns = c(2,4))


####all mutant VS minigenome, Correlation#######
data1 %>% filter(`AA` %in% points) %>% select("AA", "P1_withBXA_meanRF", "P1_withoutBXA_meanRF") %>% filter(! `AA` %in% c("I38M", "E198D", "E181D"))->data4
ggpairs(merge(data4, minigenome_mean, by = "AA"), columns = c(2,4))