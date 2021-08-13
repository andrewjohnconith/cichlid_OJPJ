##Read in packages
library(ggplot2); library(dplyr); library(ggpmisc); library(car)

##Read in data
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Expression/Statistics")
ScatterSpeciesExpression<-read.csv("AllExpressionData.csv")

#Filter
SpeciesDYM <- ScatterSpeciesExpression %>%
  filter(GeneID=="DYM")
SpeciesSMAD <- ScatterSpeciesExpression %>%
  filter(GeneID=="SMAD7")
SpeciesNOTCH <- ScatterSpeciesExpression %>%
  filter(GeneID=="NOTCH1a")

####Statistics####
##T-Tests
#dym
SpeciesDYM %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops") %>%
{t.test(.$lOJ ~ .$TaxaID, var.equal = TRUE)}

SpeciesDYM %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops") %>%
{t.test(.$lPJ ~ .$TaxaID, var.equal = TRUE)}

#smad7
SpeciesSMAD %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops") %>%
{t.test(.$lOJ ~ .$TaxaID, var.equal = TRUE)}

SpeciesSMAD %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops") %>%
{t.test(.$lPJ ~ .$TaxaID, var.equal = TRUE)}

#notch1a
SpeciesNOTCH %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops") %>%
{t.test(.$lOJ ~ .$TaxaID, var.equal = TRUE)}

SpeciesNOTCH %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops") %>%
{t.test(.$lPJ ~ .$TaxaID, var.equal = TRUE)}



##Bartlett's K-squared for homogeneity of variances
ScatterSpeciesExpression<-read.csv("FullSpeciesExpression_boxplot.csv")

#Filter
SpeciesDYM <- ScatterSpeciesExpression %>%
  filter(GeneID=="DYM")
SpeciesSMAD <- ScatterSpeciesExpression %>%
  filter(GeneID=="SMAD7")
SpeciesNOTCH <- ScatterSpeciesExpression %>%
  filter(GeneID=="NOTCH1a")

##smad7
#Compare b/w tissues
#LF OJ and PJ exhibit differences in their variances, OJ is a little more constrained
SpeciesSMAD %>% filter(TaxaID=="Labeotropheus") %>%
{bartlett.test(.$AllSpeciesData ~ .$TissueID)}
#Bartlett's K-squared = 2.92, df = 1, p-value = 0.08749

SpeciesSMAD %>% filter(TaxaID=="Tropheops") %>%
{bartlett.test(.$AllSpeciesData ~ .$TissueID)}
#Bartlett's K-squared = 0.57529, df = 1, p-value = 0.4482

#Compare b/w taxa
#LF OJ is more constrained relative to TRC OJ
#PJ variances are similar between taxa
SpeciesSMAD %>%
  filter(TaxaID=="Labeotropheus" | TaxaID=="Tropheops") %>%
  filter(TissueID=="OJ") %>%
  {bartlett.test(.$AllSpeciesData ~ .$TaxaID)}
#Bartlett's K-squared = 4.7441, df = 1, p-value = 0.0294

SpeciesSMAD %>%
  filter(TaxaID=="Labeotropheus" | TaxaID=="Tropheops") %>%
  filter(TissueID=="PJ") %>%
  {bartlett.test(.$AllSpeciesData ~ .$TaxaID)}
#Bartlett's K-squared = 0.048121, df = 1, p-value = 0.8264

#Compare b/w taxa
#OJ and PJ variances are similar between taxa
SpeciesDYM %>%
  filter(TaxaID=="Labeotropheus" | TaxaID=="Tropheops") %>%
  filter(TissueID=="OJ") %>%
  {bartlett.test(.$AllSpeciesData ~ .$TaxaID)}
#Bartlett's K-squared = 0.2056, df = 1, p-value = 0.6502

SpeciesDYM %>%
  filter(TaxaID=="Labeotropheus" | TaxaID=="Tropheops") %>%
  filter(TissueID=="PJ") %>%
  {bartlett.test(.$AllSpeciesData ~ .$TaxaID)}
#Bartlett's K-squared = 0.004837, df = 1, p-value = 0.9446
