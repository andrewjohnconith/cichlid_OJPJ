##Read in packages
library(ggplot2); library(dplyr); library(ggpmisc); library(car)

##Read in data
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Expression/Regressions/")
ScatterSpeciesExpression<-read.csv("AllExpressionData.csv")

#Filter
SpeciesDYM <- ScatterSpeciesExpression %>%
  filter(Gene=="DYM")
SpeciesSMAD <- ScatterSpeciesExpression %>%
  filter(Gene=="SMAD7")
SpeciesNOTCH <- ScatterSpeciesExpression %>%
  filter(Gene=="NOTCH1a")


####Regressions####
#GGplot plotting
MyRatio <- SpeciesDYM %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops"|TaxaID=="Maylandia") %>% 
  with(diff(range(lOJ))/diff(range(lPJ)))

DYMplotLMall <- SpeciesDYM %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops"|TaxaID=="Maylandia") %>% 
  ggplot(aes(x=lOJ, y=lPJ)) +
  stat_smooth(method=lm, fullrange = TRUE, se = T, col="black") +
  geom_point(aes(colour = factor(TaxaID)), size = 2.5) +
  labs(colour = "Clade", x = "L. Oral Jaw Relative Expression", y = "L. Pharyngeal Jaw Relative Expression") +
  scale_color_manual(values=c("#0571b0","#af8dc3","#ca0020"), name="Clade", labels=c("Labeotropheus","Maylandia","Tropheops")) +
  theme_bw() + coord_fixed(ratio=MyRatio)

Apt <- DYMplotLMall + stat_fit_glance(method = "lm", label.x = "left", label.y = "top",method.args = list(formula = y ~ x),mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),parse = TRUE)

#GGplot plotting
MyRatio <- SpeciesSMAD %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops"|TaxaID=="Maylandia") %>% 
  with(diff(range(lOJ))/diff(range(lPJ)))

SMADplotLMall <- SpeciesSMAD %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops"|TaxaID=="Maylandia") %>% 
  ggplot(aes(x=lOJ, y=lPJ)) +
  stat_smooth(method=lm, fullrange = TRUE, se = T, col="black") +
  geom_point(aes(colour = factor(TaxaID)), size = 2.5) +
  labs(colour = "Clade", x = "L. Oral Jaw Relative Expression", y = "L. Pharyngeal Jaw Relative Expression") +
  scale_color_manual(values=c("#0571b0","#af8dc3","#ca0020"), name="Clade", labels=c("Labeotropheus","Maylandia","Tropheops")) +
  theme_bw() + coord_fixed(ratio=MyRatio)


Bpt <- SMADplotLMall + stat_fit_glance(method = "lm", label.x = "left", label.y = "top",method.args = list(formula = y ~ x),mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),parse = TRUE)

#Notch1a GGplot plotting
MyRatio <- SpeciesNOTCH %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops"|TaxaID=="Maylandia") %>% 
  with(diff(range(lOJ))/diff(range(lPJ)))

NOTCHplotLMall <- SpeciesNOTCH %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops"|TaxaID=="Maylandia") %>% 
  ggplot(aes(x=lOJ, y=lPJ)) +
  stat_smooth(method=lm, fullrange = TRUE, se = T, col="black") +
  geom_point(aes(colour = factor(TaxaID)), size = 2.5) +
  labs(colour = "Clade", x = "L. Oral Jaw Relative Expression", y = "L. Pharyngeal Jaw Relative Expression") +
  scale_color_manual(values=c("#0571b0","#af8dc3","#ca0020"), name="Clade", labels=c("Labeotropheus","Maylandia","Tropheops")) +
  theme_bw() + coord_fixed(ratio=MyRatio)

Cpt <- NOTCHplotLMall  + stat_fit_glance(method = "lm", label.x = "right", label.y = "top",method.args = list(formula = y ~ x),mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),parse = TRUE)

library(cowplot)
plot_grid(Apt, Bpt, Cpt,
          labels = c("Dym", "Smad7", "Notch1a"),
          nrow = 1)

##Smad7 Inset plotting
#GGplot plotting
MyRatio <- SpeciesSMAD %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops"|TaxaID=="Maylandia") %>% 
  with(diff(range(lOJ))/diff(range(lPJ)))

SMADplotLM <- SpeciesSMAD %>% filter(TaxaID=="Labeotropheus"|TaxaID=="Tropheops"|TaxaID=="Maylandia") %>% 
  ggplot(aes(x=lOJ, y=lPJ, group=as.factor(TaxaID))) +
  stat_smooth(method=lm, fullrange = TRUE, se = F, aes(colour = factor(TaxaID))) +
  geom_point(aes(colour = factor(TaxaID)), size = 2.5) +
  labs(colour = "Clade", x = "L. Oral Jaw Relative Expression", y = "L. Pharyngeal Jaw Relative Expression") +
  scale_color_manual(values=c("#0571b0","#af8dc3","#ca0020"), name="Clade", labels=c("Labeotropheus","Maylandia","Tropheops")) +
  theme_bw() + coord_fixed(ratio=MyRatio)

Ypt <- SMADplotLM + stat_fit_glance(method = "lm", label.x = "right", label.y = "top",method.args = list(formula = y ~ x),mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g',stat(r.squared), stat(p.value))),parse = TRUE)
Ypt
