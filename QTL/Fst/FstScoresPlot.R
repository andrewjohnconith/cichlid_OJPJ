####Fst Plot####
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/QTL/Fst/")
#Read in data
FullMap<-read.csv("FstScores.csv")

#Insert a user defined alpha value
#To adjust the opacity range play with the
#'scale_alpha_continuous' and the 0 value below
FullMap$FstAlpha<-FullMap$Fst
FullMap$FstAlpha[which(FullMap$FstAlpha<0.6)]<-0

Zrange <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
FullMap$FstAlpha[which(FullMap$FstAlpha>0.6)]<-Zrange(FullMap$FstAlpha[which(FullMap$FstAlpha>0.6)])

FineMap_Fst <- filter(FullMap, !is.na(Fst)) %>%
  ggplot(aes(x=RelativePOS, y=Fst))

FineMap_FstPoint <- FineMap_Fst + theme_bw() + annotate("rect", xmin = -1286590, xmax = 2195347, ymin = -Inf, ymax = Inf, fill = "#404040", alpha = 0.2) +
  geom_point(size=2, aes(alpha=FstAlpha)) +
  scale_alpha_continuous(range = c(0.2,1)) +
  labs(x = "Marker position on scaffold:0/21 (Mbs)", y = "Fst") +
  theme(legend.position = "none", aspect.ratio = 1, legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

FineMap_FstPointGene <- FineMap_FstPoint + 
  geom_hline(yintercept = 0.9, col="dark green", linetype="longdash") + #Z=2 position
  geom_hline(yintercept = 0.6, col="dark green", linetype="dashed") #Z=1 position

FineMap_FstPointGene
