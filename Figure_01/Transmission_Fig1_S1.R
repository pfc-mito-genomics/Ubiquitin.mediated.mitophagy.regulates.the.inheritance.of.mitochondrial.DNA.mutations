library(tidyverse)
library(ggpubr)

#function to make a qq plot
makeqq <- function(d1){
  md1 <- mean(d1)
  sd1 <- sd(d1)
  standardized <- (d1-md1)/sd1
  qqnorm( standardized )
  abline(0,1,col="red")
}

##############################################################################
################## Mother-to-offspring heteroplasmy analysis##################
##############################################################################
#run this to clear the environment
rm(list = ls())

#set the  Figure_01 folder. The path depends on which folder you currently are
setwd("")

#import data
All_Ear_data <- read.csv("All_Ear_data.csv")
All_Ear_data$Usp30_allele <- factor(All_Ear_data$Usp30_allele, levels = c("WT","Het", "KO"))


###check that Heteroplasmy  column is copied as value
data.class(All_Ear_data$Het)

#transform the heteroplasmy data into Heteroplasmic Fraction data
All_Ear_data <- All_Ear_data %>% 
  mutate(Mother_HF=Mother_Het/100, HF=Het/100)

#calculate heteroplasmy shift
All_Ear_data$TransfShift <- rep(0,nrow(All_Ear_data))
for(i in 1:nrow(All_Ear_data)){
  All_Ear_data$TransfShift[i] <- log((All_Ear_data$HF[i]*(1-All_Ear_data$Mother_HF[i]))/(All_Ear_data$Mother_HF[i]*(1-All_Ear_data$HF[i])))
}

All_Ear_data$CrossType <- factor(All_Ear_data$CrossType, levels = c("0","1","2","3"))

All_Ear_data <- All_Ear_data %>% 
  mutate(Usp30_genotype = case_when(Usp30_allele == "WT" ~ "+/+",
                                    Usp30_allele == "Het" ~ "+/-",
                                    Usp30_allele == "KO" ~ "-/-"))

All_Ear_data$Usp30_genotype <- factor(All_Ear_data$Usp30_genotype, levels = c("+/+","+/-","-/-"))

#check the Mendelian balance of F2
TotMice <- nrow(All_Ear_data[which(All_Ear_data$CrossType == "2" & All_Ear_data$Usp30_allele != "NA"),]) 

Mendel <- All_Ear_data %>% 
  filter(CrossType == "2") %>% 
  group_by(CrossType, Usp30_allele) %>% 
  summarise(Nmice = n(), Frequency = Nmice/TotMice)
Mendel

chisq.test(as.double(unlist(Mendel[,4])),c(0.25,0.5,0.25))

########Mother to offspring heteroplasmy - Fig 1B and S1
All_Ear_data %>% 
  filter(CrossType == "2") %>% 
  group_by(Mother_Het, Usp30_allele) %>%
  ggplot(aes(x=Mother_Het, y= Het, colour=Usp30_allele)) +
  geom_point(alpha=0.9) +
  theme_classic2(base_size = 20) + theme(text = element_text(face = "bold"),
                                         strip.text = element_text(size = 20, face = "bold")) +
  facet_wrap(~Usp30_allele, labeller = labeller(Usp30_allele = 
                                                  c("WT" = "+/+","Het" = "+/-","KO" = "-/-"))) +
  scale_colour_manual(values = c("darkgoldenrod","#aadce0","#72bcd5")) +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(), panel.grid.major.x = element_line(),
        panel.grid.minor.y = element_line(), panel.grid.minor.x = element_line(),
        axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) +
  scale_y_continuous(breaks = seq(0,90, 10), limits = c(10,85), minor_breaks = seq(0,90,5)) +
  scale_x_continuous(breaks = seq(0,90, 10), limits = c(20,80), minor_breaks = seq(0,90,5)) +
  geom_smooth(colour = "black", method='lm',se=TRUE) +
  geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
  labs(x= "Mother heteroplasmy %", y = "Offspring heteroplasmy %",
       title ="Cross C: Mothers - Usp30+/-") 

ggsave("Fig 1B.png", width = 9, height = 7)

All_Ear_data %>% 
  filter(CrossType != "2") %>% 
  group_by(Mother_Het, Usp30_genotype) %>%
  ggplot(aes(x=Mother_Het, y= Het, colour=Usp30_genotype)) +
  geom_point(alpha=0.9) +
  theme_classic2(base_size = 15) + theme(text = element_text(face = "bold")) +
  scale_colour_manual(values = c("darkgoldenrod","#aadce0", "#72bcd5")) +
  theme(legend.position = c(0.25,0.1), legend.background = element_blank(),
        panel.grid.major.y = element_line(),  panel.grid.major.x = element_line(),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) +
  scale_y_continuous(breaks = seq(0,95, 5), limits = c(10,85)) +
  scale_x_continuous(breaks = seq(0,90, 10), limits = c(50,80)) +
  #geom_smooth(colour = "black", method='lm',se=TRUE) +
  facet_wrap(~CrossType, labeller = labeller(CrossType = c("0" = "Cross A: Mothers - Usp30+/+",
                                                           "1" = "Cross B: Mothers - Usp30+/+",
                                                           "2" = "Cross C: Mothers - Usp30+/-",
                                                           "3" = "Cross D: Mothers - Usp30-/-"))) +
  geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
  labs(x= "Mother heteroplasmy %", y = "Offspring heteroplasmy %", color = "Offspring Usp30 genotype") 

ggsave("Fig S1E.png", width = 9, height = 7)


#####Stepwise AIC analysis
StepwiseAICdata_full <- All_Ear_data %>% 
  filter(CrossType == "2") %>% 
  dplyr::select(Het, Mother_Het, Usp30_allele)

#----- Build model
lm_full <- lm(Het ~ Mother_Het * Usp30_allele, data = StepwiseAICdata_full) 
summary(lm_full)
# Reduce variables in the model using stepwise AIC selection
library("MASS")
lm_reduced <- stepAIC(lm_full)
summary(lm_reduced) # it's the same model, meaning the genotype is significant

#repeat on sane dataset, but with WT removed. To test if the two KO genotypes are significantly different
StepwiseAICdata_noWT <- StepwiseAICdata_full %>% 
  filter(Usp30_allele != "WT")

lm_full <- lm(Het ~ Mother_Het * Usp30_allele, data = StepwiseAICdata_noWT) 
summary(lm_full)
lm_reduced <- stepAIC(lm_full)
summary(lm_reduced) #Usp30 allele is removed. So they are not different

################Figure 1C
#check that the mother heteroplasmies are similar and the number of mice per group
All_Ear_data %>% 
  mutate(bins = cut(c(Mother_Het),breaks=c(20,60,80))) %>%
  group_by(Usp30_allele, CrossType, bins) %>% 
  summarise(Nmice = n(),mean(TransfShift), mean(Mother_Het))

###make heteroplasmy bins to x and y plot - Fig 1C
#first below 60%
All_Ear_data %>% 
  mutate(bins = cut(c(Mother_Het),breaks=c(20,60,80))) %>%
  filter(CrossType == "2" & bins == "(20,60]") %>% 
  ggplot(aes(x=Usp30_genotype, y= TransfShift)) +
  theme_classic2(base_size = 20) + theme(text = element_text(face = "bold"),
                                         strip.text = element_text(size = 20, face = "bold")) +
  geom_violin() +
  theme(text = element_text(face = "bold"), legend.position = "none",
        axis.text.x = element_text(size=25), axis.text.y = element_text(size=25)) +
  geom_boxplot(width = 0.15) +
  geom_dotplot(aes(fill = Usp30_genotype),binaxis = "y", stackdir = "center",
               dotsize = 0.15, position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(-2.2,2.2)) +
  stat_summary(fun.y = median, geom = "crossbar", width = 0.15, colour = "black") +
  ggforce::facet_row(CrossType~bins, scales = "free_x", space = "free", 
                     labeller = labeller(bins =c("(20,60]" = "Mothers <60%", "(60,80]" = " Mothers >60%"),
                                         CrossType = c("0" = "Cross A",
                                                       "1" = "Cross B",
                                                       "2" = "Cross C : Mothers - Usp30+/-",
                                                       "3" = "Cross D"))) +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dotted") +
  scale_colour_manual(values = c("darkgoldenrod","#aadce0","#72bcd5")) +
  scale_fill_manual(values = c("darkgoldenrod","#aadce0","#72bcd5")) +
  labs(x= "Offspring Usp30 genotype", y = "Transformed Heteroplasmy shift")

ggsave("Fig 1C_lowHet.png", width = 6.5, height = 9)

All_Ear_data %>% 
  mutate(bins = cut(c(Mother_Het),breaks=c(20,60,80))) %>%
  filter(bins == "(60,80]") %>% 
  ggplot(aes(x=Usp30_genotype, fill= Usp30_genotype, y= TransfShift)) +
  theme_classic2(base_size = 20) +
  geom_violin() +
  theme(legend.position = "none", strip.text = element_text(size = 20, face = "bold"),
        axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=25), text = element_text(face = "bold")) +
  geom_boxplot(width = 0.15) +
  geom_dotplot(binaxis = "y", stackdir = "center",
               dotsize = 0.15, position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(-2.2,2.2)) +
  stat_summary(fun.y = median, geom = "crossbar", width = 0.15, colour = "black") +
  ggforce::facet_row(CrossType~bins, scales = "free_x", space = "free", 
                     labeller = labeller(bins =c("(20,60]" = "Mothers <60%", "(60,80]" = " Mothers >60%"),
                                         CrossType = c("0" = "Cross A",
                                                       "1" = "Cross B",
                                                       "2" = "Cross C : Mothers - Usp30+/-",
                                                       "3" = "Cross D"))) +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dotted") +
  scale_fill_manual(values = c("darkgoldenrod","#aadce0","#72bcd5")) +
  labs(x= "Offspring Usp30 genotype", y = "")

ggsave("Fig 1C_highHet.png", width = 14, height = 9)

 ###############stats#########
#define control and sample as the comparison you're checking each time
control <- All_Ear_data %>% 
  filter(CrossType == "0") %>% 
  mutate(bins = cut(c(Mother_Het),breaks=c(20,60,80))) %>%
  filter(bins=="(60,80]",  Usp30_allele == "WT") %>% 
  dplyr::select(TransfShift)

sample <- All_Ear_data %>% 
  filter(CrossType == "2") %>% 
  mutate(bins = cut(c(Mother_Het),breaks=c(20,60,80))) %>%
  filter(bins=="(60,80]", Usp30_allele == "WT") %>% 
  dplyr::select(TransfShift)

#check if normally distributed
count(control)
shapiro.test(unlist(control))
makeqq(unlist(control))
count(sample)
shapiro.test(unlist(sample))
makeqq(unlist(sample))
#do appropriate test 
wilcox.test(as.double(unlist(control)),as.double(unlist(sample)))
t.test(as.double(unlist(control)),as.double(unlist(sample)))


##############################################################################
######################## Litter size analysis - Fig S1########################
##############################################################################
setwd("../Figure_S01/")

All_Ear_data %>% 
  filter(CrossType == "2") %>%
  group_by(MotherNumber,LitterNumber,Mother_Het) %>% 
  summarise(LitterSize=n()) %>% 
  #print(n=100) #uncomment to check the tibble
  mutate(bins = cut(c(Mother_Het),breaks=c(0,60,100))) %>%
  ggplot(aes(x=bins,y=LitterSize, fill = bins)) +
  theme_classic2(base_size = 15) + 
  geom_violin() +
  geom_boxplot(width = 0.15) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
  scale_fill_manual(values = c("darkcyan","darkgoldenrod")) +
  stat_summary(fun.y = median, geom = "crossbar", width = 0.15, colour = "black") +
  theme(legend.position = "none", text = element_text(face = "bold"),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) +
  #stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", ref.group = "(0,60]",label.y=14) + #uncomment this for stats
  scale_y_continuous(breaks = c(1:15)) +
  scale_x_discrete(labels=c("Low","High")) +
  labs(x= "Mother heteroplasmy", y = "Number of offspring per litter at weaning")

ggsave("Fig S1A_pretty.png", width = 5, height = 5)

#############calculating the proportion of each genotype over total#############

#creating tibbles with the number of mice per litter of each genotype
LitterSizes <- All_Ear_data %>% 
  filter(CrossType == "2") %>%
  group_by(MotherNumber,LitterNumber,Mother_Het) %>% 
  summarise(LitterSize=n())

WT_LitterSizes <- All_Ear_data %>% 
  filter(CrossType == "2") %>%
  filter(Usp30_allele == "WT") %>% 
  group_by(MotherNumber,LitterNumber,Mother_Het) %>% 
  summarise(LitterSize=n()) 

Het_LitterSizes <- All_Ear_data %>% 
  filter(CrossType == "2") %>%
  filter(Usp30_allele == "Het") %>% 
  group_by(MotherNumber,LitterNumber,Mother_Het) %>% 
  summarise(LitterSize=n())

KO_LitterSizes <- All_Ear_data %>% 
  filter(CrossType == "2") %>%
  filter(Usp30_allele == "KO") %>% 
  group_by(MotherNumber,LitterNumber,Mother_Het) %>% 
  summarise(LitterSize=n()) 

#making a unique label for each litter
LitterSizes$Litter <- paste(LitterSizes$MotherNumber,LitterSizes$LitterNumber,sep = ".")
WT_LitterSizes$Litter <- paste(WT_LitterSizes$MotherNumber,WT_LitterSizes$LitterNumber,sep = ".")
Het_LitterSizes$Litter <- paste(Het_LitterSizes$MotherNumber,Het_LitterSizes$LitterNumber,sep = ".")
KO_LitterSizes$Litter <- paste(KO_LitterSizes$MotherNumber,KO_LitterSizes$LitterNumber,sep = ".")

#inserting genotype-specific litter sizes into final data.frame
LitterSizes$WTSize <- rep(0,nrow(LitterSizes))
for(i in 1:nrow(LitterSizes)){
  LitterSizes$WTSize[i] <- WT_LitterSizes$LitterSize[match(LitterSizes$Litter[i],WT_LitterSizes$Litter)]
}
LitterSizes$WTSize[which(is.na(LitterSizes$WTSize))] <- 0

LitterSizes$HetSize <- rep(0,nrow(LitterSizes))
for(i in 1:nrow(LitterSizes)){
  LitterSizes$HetSize[i] <- Het_LitterSizes$LitterSize[match(LitterSizes$Litter[i],Het_LitterSizes$Litter)]
}
LitterSizes$HetSize[which(is.na(LitterSizes$HetSize))] <- 0

LitterSizes$KOSize <- rep(0,nrow(LitterSizes))
for(i in 1:nrow(LitterSizes)){
  LitterSizes$KOSize[i] <- KO_LitterSizes$LitterSize[match(LitterSizes$Litter[i],KO_LitterSizes$Litter)]
}
LitterSizes$KOSize[which(is.na(LitterSizes$KOSize))] <- 0

rm(WT_LitterSizes,Het_LitterSizes,KO_LitterSizes)

#checking that the sum of the offspring for each litter is always like the total
LitterSizes %>% 
  mutate(WT_percent=(WTSize/LitterSize)*100,
         KO_percent=(KOSize/LitterSize)*100,
         Het_percent=(HetSize/LitterSize)*100,
         check = sum(WTSize,HetSize,KOSize)) %>% 
  filter(check == 0) 

LitterSizes %>% 
  mutate(WT_percent=(WTSize/LitterSize)*100,
         KO_percent=(KOSize/LitterSize)*100,
         Het_percent=(HetSize/LitterSize)*100,
         check = sum(WTSize,HetSize,KOSize)) %>%
  filter(check != 0) %>% 
  mutate(bins = cut(c(Mother_Het),breaks=c(0,60,100))) %>%
  ggplot(aes(x=bins,y=KO_percent, fill = bins)) + #change the KO to WT or Het for each plot
  geom_boxplot(width = 0.15) +
  geom_violin() +
  theme_classic2(base_size = 15) + 
  scale_fill_manual(values = c("darkcyan","darkgoldenrod")) +
  theme(legend.position = "none", text = element_text(face = "bold"),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
  stat_summary(fun.y = mean, geom = "crossbar", width = 0.5, colour = "black") +
  scale_y_continuous(breaks = seq(0,100, by=25), limits = c(0,125))+
  #stat_compare_means(aes(label = ..p.format..), method = "wilcox.test", ref.group = "(0,60]",label.y=120) + #uncomment for stats
  scale_x_discrete(labels=c("Low","High")) +
  geom_hline(aes(yintercept = 25), color = "red", linetype = "dotted") +
  labs(title="Usp30-/-", #change the label for each plot
       x= "Mother heteroplasmy", y = "allele prevalence in litter (%)")

ggsave("Fig S1J.png", width = 5, height = 5)

write.csv(x = LitterSizes, file = "LitterSizes.csv", row.names = FALSE)
#no special stats here as the distributions were either clearly non-parametric (WT and KO) or very similar
