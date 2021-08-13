##Read in packages
library(ggplot2); library(dplyr)

##Read in data
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Expression/Boxplots/")
FullSpeciesExpression<-read.csv("FullSpeciesExpression_boxplot.csv")

####DYM####
SpeciesOJPJ <- FullSpeciesExpression %>% 
  filter(TaxaID=="Labeotropheus" | TaxaID=="Tropheops") %>% 
  filter(GeneID=="DYM") %>% 
  select(TaxaID, TissueID, AllSpeciesData) 

# Calculates mean, sd, se and IC
my_sum_OJPJ <- SpeciesOJPJ %>%
  group_by(TaxaID, TissueID) %>%
  summarise( 
    n=n(),
    mean=mean(AllSpeciesData, na.rm = T),
    sd=sd(AllSpeciesData, na.rm = T)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# Standard error
OJPJChart_DYM <- ggplot(my_sum_OJPJ, aes(x=TissueID, y=mean, fill=factor(TaxaID))) +
  geom_bar(stat="identity", position='dodge', alpha=0.5) +
  geom_point(data=SpeciesOJPJ,aes(y=AllSpeciesData, x=TissueID, group = factor(TaxaID), color=factor(TaxaID)), position = position_dodge(1))+
  scale_color_manual(values=c("#0571b0","#ca0020"), name="Clade", labels=c("Labeotropheus", "Tropheops")) +
  scale_fill_manual(values=c("#0571b0","#ca0020"), name="Clade", labels=c("Labeotropheus", "Tropheops")) +
  geom_errorbar(aes(x=TissueID, ymin=mean-se, ymax=mean+se), position = position_dodge(.9), width=0.4, colour="orange", alpha=0.9, size=1.5) +
  ggtitle("Dym Expression Level (±se)") + labs(y = "Relative Expression", x = "Taxa", fill = "Taxa") +
  theme_bw() + theme(aspect.ratio=1, legend.position = "none") + coord_fixed()

OJPJChart_DYM


####SMAD7####
SpeciesOJPJ <- FullSpeciesExpression %>% 
  filter(TaxaID=="Labeotropheus" | TaxaID=="Tropheops") %>% 
  filter(GeneID=="SMAD7") %>% 
  select(TaxaID, TissueID, AllSpeciesData) 

# Calculates mean, sd, se and IC
my_sum_OJPJ <- SpeciesOJPJ %>%
  group_by(TaxaID, TissueID) %>%
  summarise( 
    n=n(),
    mean=mean(AllSpeciesData, na.rm = T),
    sd=sd(AllSpeciesData, na.rm = T)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# Standard error
OJPJChart_SMAD7 <- ggplot(my_sum_OJPJ, aes(x=TissueID, y=mean, fill=factor(TaxaID))) +
  geom_bar(stat="identity", position='dodge', alpha=0.5) +
  geom_point(data=SpeciesOJPJ,aes(y=AllSpeciesData, x=TissueID, group = factor(TaxaID), color=factor(TaxaID)), position = position_dodge(1))+
  scale_color_manual(values=c("#0571b0","#ca0020"), name="Clade", labels=c("Labeotropheus", "Tropheops")) +
  scale_fill_manual(values=c("#0571b0","#ca0020"), name="Clade", labels=c("Labeotropheus", "Tropheops")) +
  geom_errorbar(aes(x=TissueID, ymin=mean-se, ymax=mean+se), position = position_dodge(.9), width=0.4, colour="orange", alpha=0.9, size=1.5) +
  ggtitle("Smad7 Expression Level (±se)") + labs(y = "Relative Expression", x = "Taxa", fill = "Taxa") +
  theme_bw() + theme(aspect.ratio=1, legend.position = "none") + coord_fixed()

OJPJChart_SMAD7



####NOTCH1a####
SpeciesOJPJ <- FullSpeciesExpression %>% 
  filter(TaxaID=="Labeotropheus" | TaxaID=="Tropheops") %>% 
  filter(GeneID=="NOTCH1a") %>% 
  select(TaxaID, TissueID, AllSpeciesData) 

# Calculates mean, sd, se and IC
my_sum_OJPJ <- SpeciesOJPJ %>%
  group_by(TaxaID, TissueID) %>%
  summarise( 
    n=n(),
    mean=mean(AllSpeciesData, na.rm = T),
    sd=sd(AllSpeciesData, na.rm = T)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# Standard error
OJPJChart_NOTCH1a <- ggplot(my_sum_OJPJ, aes(x=TissueID, y=mean, fill=factor(TaxaID))) +
  geom_bar(stat="identity", position='dodge', alpha=0.5) +
  geom_point(data=SpeciesOJPJ,aes(y=AllSpeciesData, x=TissueID, group = factor(TaxaID), color=factor(TaxaID)), position = position_dodge(1))+
  scale_color_manual(values=c("#0571b0","#ca0020"), name="Clade", labels=c("Labeotropheus", "Tropheops")) +
  scale_fill_manual(values=c("#0571b0","#ca0020"), name="Clade", labels=c("Labeotropheus", "Tropheops")) +
  geom_errorbar(aes(x=TissueID, ymin=mean-se, ymax=mean+se), position = position_dodge(.9), width=0.4, colour="orange", alpha=0.9, size=1.5) +
  ggtitle("Notch1a Expression Level (±se)") + labs(y = "Relative Expression", x = "Taxa", fill = "Taxa") +
  theme_bw() + theme(aspect.ratio=1, legend.position = "none") + coord_fixed()

OJPJChart_NOTCH1a



####Combined####
library("cowplot")
plot_grid(OJPJChart_DYM, OJPJChart_SMAD7, OJPJChart_NOTCH1a,
          nrow = 1,
          ncol = 3)
