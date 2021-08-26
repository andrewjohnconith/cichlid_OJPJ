library(geomorph); library(geiger); library(plotrix); library(caper); library(ggplot2)

####Data Engineering####
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Macroevolution")

LPJ_CichlidMeanShapes<-readland.tps(file = "MacroLPJ_Mean_20210826.tps", specID = "ID")
LOJ_CichlidMeanShapes<-readland.tps(file = "MacroLOJ_Mean_20210826.tps", specID = "ID")

TaxaIDData<-read.csv("TaxaData5.csv", row.names = 1)
Tree<-read.tree("CichlidTreeNature2020.tre")

Pruning<-treedata(Tree, TaxaIDData)
Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid<-dimnames(LPJ_CichlidMeanShapes)[[3]][dimnames(LPJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes[,,RM_Cichlid]

RM_Cichlid<-dimnames(LOJ_CichlidMeanShapes)[[3]][dimnames(LOJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes[,,RM_Cichlid]


LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LPJ_CichlidMeanShapes_Tr)[[3]]))]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LOJ_CichlidMeanShapes_Tr)[[3]]))]

##Group Testing
gp.Lake <- as.factor(FData$Lake)
gp.Clade <- as.factor(FData$Clade)
gp.Diet <- as.factor(FData$Diet)

#Sample size check
table(gp.Lake) 
#A = African Nile Basin; M = Malawi; T = Tanganyika; V = Victoria
  # A  M  T  V 
  #14 40 29  5

table(gp.Clade)
#A = African Nile Basin; M = Mbuna (Malawi); T = Tanganyika; U = Utaka (Malawi); V = Victoria
  # A  M  T  U  V 
  #14 20 29 20  5 

table(gp.Diet)
#A = Aufwuchs, F = Fish, O = Omnivore, S = Scale, Z = Zooplankton, Z1 = Zoobenthos1, Z2 = Zoobenthos2
  # A  F  O  S  Z Z1 Z2 
  #27 12 13  2  7  8 19



####Phylo-Covariation Analysis####
phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree) #r-PLS: 0.506; P-value: 0.004; Effect Size: 2.6578
plot(phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree), label = dimnames(LOJ_CichlidMeanShapes_Tr)[[3]])


#X-Y Coords
LOJLPJcoords<-phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree)

LOJLPJ<-cbind.data.frame(LOJLPJcoords$XScores[,1],LOJLPJcoords$YScores[,1], FData$Clade, FData$Lake, FData$Diet)
colnames(LOJLPJ)<-c("XBlock", "YBlock", "Clade", "Lake", "Diet")
#write.csv(LOJLPJ, "CorrData.csv")

####Integration Plotting####
library(ggplot2)
#LOJLPJ<-read.csv("CorrData.csv")

MyRatio <- with(LOJLPJ, diff(range(XBlock))/diff(range(YBlock)))
#Region-Basin
p <- ggplot(LOJLPJ, aes(x=XBlock, y=YBlock)) +
  geom_point(aes(colour = factor(Clade), shape = factor(Lake)), size = 2.5) +
  labs(colour = "Clade", shape = "Lake", x = "L. Oral Jaw Block", y = "L. Pharyngeal Jaw Block") +
  theme_bw() + coord_fixed(ratio=MyRatio)
p

#Diet
r <- ggplot(LOJLPJ, aes(x=XBlock, y=YBlock)) +
  geom_point(aes(colour = factor(Diet), shape = factor(Lake)), size = 2.5) +
  labs(colour = "Diet", shape = "Lake", x = "L. Oral Jaw Block", y = "L. Pharyngeal Jaw Block") +
  theme_bw() + coord_fixed(ratio=MyRatio)
r

####Extract pPLS Line####
##Extracts the coefficients directly from the pPLS analysis
#Run function below to extract the coefficients so you can plot it in the ggplot2 framework
#pPLS_Info is a global assignment
MYPLSFUNCTIONT(LOJLPJcoords)
pPLSCoefficients<-pPLS_Info

p + geom_abline(intercept = pPLSCoefficients$coefficients[1], slope = pPLSCoefficients$coefficients[2]) +
  scale_shape_manual(values = c(15, 19, 17, 18), name="Lake", labels=c("Africa","Malawi", "Tanganyika", "Victoria")) +
  scale_colour_manual(values=c("#bdbdbd","#f3766d","#7caf41","#1cbdc2", "#a781ba"), name="Clade", labels=c("Africa","Mbuna", "Tanganyika", "Utaka", "Victoria")) +
  #scale_shape_discrete(name="Lake", labels=c("Africa","Malawi", "Tanganyika", "Victoria")) +
  theme(legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1)) #+
  #geom_text(aes(label=rownames(LOJLPJ)),hjust=0, vjust=0)

r + geom_abline(intercept = pPLSCoefficients$coefficients[1], slope = pPLSCoefficients$coefficients[2]) +
  scale_shape_manual(values = c(15, 19, 17, 18), name="Lake", labels=c("Africa","Malawi", "Tanganyika", "Victoria")) +
  scale_colour_manual(values=c("#e41a1c", "#377eb8" ,"#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628"), name="Diet", labels=c("Aufwuchs", "Fish", "Omnivore", "Scale", "Zooplankton", "Zoobenthos1", "Zoobenthos2")) +
  #scale_shape_discrete(name="Lake", labels=c("Africa","Malawi", "Tanganyika", "Victoria")) +
  theme(legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1)) #+
  #geom_text(aes(label=rownames(LOJLPJ)),hjust=0, vjust=0)


#Extract out the pPLS line with this function
MYPLSFUNCTIONT<-function (x, label = NULL, warpgrids = TRUE, shapes = FALSE, 
                         ...) 
{
  A1 <- x$A1
  A2 <- x$A2
  XScores <- x$XScores
  YScores <- x$YScores
  if (is.matrix(XScores)) 
    XScores <- XScores[, 1]
  if (is.matrix(YScores)) 
    YScores <- YScores[, 1]
  Xmin <- min(XScores)
  Xmax <- max(XScores)
  Ymin <- min(YScores)
  Ymax <- max(YScores)
  plsRaw <- geomorph:::pls(x$A1.matrix, x$A2.matrix, verbose = TRUE)
  XScoresRaw <- plsRaw$XScores[, 1]
  YScoresRaw <- plsRaw$YScores[, 1]
  pc <- prcomp(cbind(XScores, YScores))$x[, 1]
  px <- predict(lm(XScores ~ pc))
  py <- predict(lm(YScores ~ pc))
  pxmax <- max(px)
  pxmin <- min(px)
  pymax <- max(py)
  pymin <- min(py)
  pcRaw <- prcomp(cbind(XScoresRaw, YScoresRaw))$x[, 1]
  pxRaw <- predict(lm(XScoresRaw ~ pcRaw))
  pyRaw <- predict(lm(YScoresRaw ~ pcRaw))
  if (length(dim(A1)) == 3) {
    A1.ref <- mshape(A1)
    preds <- shape.predictor(A1, x = XScores, method = "LS", 
                             Intercept = TRUE, pred1 = Xmin, pred2 = Xmax)
    pls1.min <- preds$pred1
    pls1.max <- preds$pred2
  }
  if (length(dim(A2)) == 3) {
    A2.ref <- mshape(A2)
    preds <- shape.predictor(A2, x = YScores, method = "LS", 
                             Intercept = TRUE, pred1 = Ymin, pred2 = Ymax)
    pls2.min <- preds$pred1
    pls2.max <- preds$pred2
  }
  if (length(dim(A1)) != 3 && length(dim(A2)) != 3) {
    plot(XScores, YScores, pch = 21, bg = "black", main = "PLS Plot", 
         xlab = "PLS1 Block 1", ylab = "PLS1 Block 2")
    abline(lm(py ~ px), col = "red")
    if (length(label != 0)) {
      text(XScores, YScores, label, adj = c(-0.7, -0.7))
    }
  }
  if (length(dim(A1)) == 3 || length(dim(A2)) == 3) {
    par(mar = c(1, 1, 1, 1) + 0.1)
    split.screen(matrix(c(0.22, 1, 0.22, 1, 0.19, 0.39, 0, 
                          0.19, 0.8, 1, 0, 0.19, 0, 0.19, 0.19, 0.39, 0, 0.19, 
                          0.8, 1), byrow = TRUE, ncol = 4))
    screen(1)
    plot(XScores, YScores, pch = 21, bg = "black", main = "PLS1 Plot: Block 1 (X) vs. Block 2 (Y) ", 
         xlab = "PLS1 Block 1", ylab = "PLS1 Block 2")
    abline(lm(py ~ px), col = "red")
    if (length(label != 0)) {
      text(XScores, YScores, label, adj = c(-0.7, -0.7))
    }
    if (warpgrids == TRUE) {
      if (length(dim(A1)) == 3 && dim(A1)[2] == 2) {
        screen(2)
        tps(A1.ref, pls1.min, 20, sz = 0.7)
        screen(3)
        tps(A1.ref, pls1.max, 20, sz = 0.7)
      }
      if (length(dim(A2)) == 3 && dim(A2)[2] == 2) {
        screen(4)
        tps(A2.ref, pls2.min, 20, sz = 0.7)
        screen(5)
        tps(A2.ref, pls2.max, 20, sz = 0.7)
      }
    }
    close.screen(all.screens = TRUE)
    par(mar = c(5.1, 4.1, 4.1, 2.1))
  }
  if (length(dim(A1)) == 3 && dim(A1)[2] == 3) {
    #plot(XScores, YScores, pch = 21, bg = "black", main = "PLS Plot", 
    #     xlab = "PLS1 Block 1", ylab = "PLS1 Block 2")
    if (length(label != 0)) {
    #  text(XScores, YScores, label, adj = c(-0.7, -0.7))
    }
    pPLS_Info<<-lm(py ~ px)
    #abline(lm(py ~ px), col = "red")
    #open3d()
   # mfrow3d(1, 2)
   # plot3d(pls1.min, type = "s", col = "gray", main = paste("PLS Block1 negative"), 
    #       size = 1.25, aspect = FALSE, xlab = "", ylab = "", 
     #      zlab = "", box = FALSE, axes = FALSE)
    #plot3d(pls1.max, type = "s", col = "gray", main = paste("PLS Block1 positive"), 
     #      size = 1.25, aspect = FALSE, xlab = "", ylab = "", 
      #     zlab = "", box = FALSE, axes = FALSE)
  }
  if (length(dim(A2)) == 3 && dim(A2)[2] == 3) {
   # open3d()
  #  mfrow3d(1, 2)
  #  plot3d(pls2.min, type = "s", col = "gray", main = paste("PLS Block2 negative"), 
   #        size = 1.25, aspect = FALSE, xlab = "", ylab = "", 
    #       zlab = "", box = FALSE, axes = FALSE)
    #plot3d(pls2.max, type = "s", col = "gray", main = paste("PLS Block2 positive"), 
    #       size = 1.25, aspect = FALSE, xlab = "", ylab = "", 
    #       zlab = "", box = FALSE, axes = FALSE)
  }
  layout(1)
  if (shapes == TRUE) {
    if (length(dim(A1)) == 3 || length(dim(A2)) == 3) {
      rtrn <- list()
      if (length(dim(A1)) == 3) {
        rtrn$pls1.min = pls1.min
        rtrn$pls1.max = pls1.max
      }
      if (length(dim(A2)) == 3) {
        rtrn$pls2.min = pls2.min
        rtrn$pls2.max = pls2.max
      }
    }
    if (length(dim(A1)) == 3 || length(dim(A2)) == 3) 
      return(rtrn)
  }
}



####Lake Basin / Region Morphospace####
#Basin - Region
gp.Clade <- as.factor(FData$Clade)
ColourVec<-c("#e66101", "#1a9641", "#2c7bb6", "#7b3294","#fdb863")
ShapeVec<-c(15,19,17,19,18)

#Plots
par(pty='s')
#LOJ
MandiblePCAs<-gm.prcomp(A = LOJ_CichlidMeanShapes_Tr, phy = FTree, align.to.phy = F)
pc.plot<-plot(MandiblePCAs, axis1 = 1, axis2 = 2, phylo = F, col=ColourVec[unclass(gp.Clade)], pch=ShapeVec[unclass(gp.Clade)])
shapeHulls(pc.plot, groups = as.character(gp.Clade), group.cols = c("#e66101", "#2c7bb6", "#fdb863", "#7b3294","#1a9641"))
pc.means1 <- aggregate(MandiblePCAs$x~gp.Clade, FUN=mean)
rownames(pc.means1) <- pc.means1[,1]
pc.means1 <- pc.means1[,-1]
points(pc.means1[,1:2], asp=1, pch=24, bg=ColourVec, cex=3)
text(MandiblePCAs$x, pos = 4, label = dimnames(LOJ_CichlidMeanShapes_Tr)[[3]])
#picknplot.shape(plot(MandiblePCAs, axis1 = 1, axis2 = 2, phylo = F, col=ColourVec[unclass(gp.Clade)], pch=ShapeVec[unclass(gp.Clade)]))

#LPJ
LPJPCAs<-gm.prcomp(A = LPJ_CichlidMeanShapes_Tr, phy = FTree, align.to.phy = F)
pc.plot<-plot(LPJPCAs, axis1 = 1, axis2 = 2, phylo = F, col=ColourVec[unclass(gp.Clade)], pch=ShapeVec[unclass(gp.Clade)])
shapeHulls(pc.plot, groups = as.character(gp.Clade), group.cols = c("#e66101", "#2c7bb6", "#fdb863", "#7b3294","#1a9641"))
pc.means1 <- aggregate(LPJPCAs$x~gp.Clade, FUN=mean)
rownames(pc.means1) <- pc.means1[,1]
pc.means1 <- pc.means1[,-1]
points(pc.means1[,1:2], asp=1, pch=24, bg=ColourVec, cex=3)
text(LPJPCAs$x, pos = 4, label = dimnames(LPJ_CichlidMeanShapes_Tr)[[3]])



####Diet Morphospace####
#Diet
gp.Diet <- as.factor(FData$Diet)
ColourVec<-c("#e41a1c", "#377eb8" ,"#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")
ShapeVec<-c(19,19,19,19,19,19,19)

#LOJ
MandiblePCAs<-gm.prcomp(A = LOJ_CichlidMeanShapes_Tr, phy = FTree, align.to.phy = F)
pc.plot<-plot(MandiblePCAs, axis1 = 1, axis2 = 2, phylo = F, col=ColourVec[unclass(gp.Diet)], pch=ShapeVec[unclass(gp.Diet)])
shapeHulls(pc.plot, groups = as.character(gp.Diet), group.cols = c("#4daf4a", "#377eb8", "#ff7f00", "#984ea3","#a65628", "#e41a1c", "#ffff33"))
pc.means1 <- aggregate(MandiblePCAs$x~gp.Diet, FUN=mean)
rownames(pc.means1) <- pc.means1[,1]
pc.means1 <- pc.means1[,-1]
points(pc.means1[,1:2], asp=1, pch=24, bg=ColourVec, cex=3)
text(MandiblePCAs$x, pos = 4, label = dimnames(LOJ_CichlidMeanShapes_Tr)[[3]])
legend("topleft", legend=levels(gp.Diet), col=ColourVec, pch = 19, cex=0.8)
#picknplot.shape(plot(MandiblePCAs, axis1 = 1, axis2 = 2, phylo = F, col=ColourVec[unclass(gp.Diet)], pch=ShapeVec[unclass(gp.Diet)]))

#LPJ
LPJPCAs<-gm.prcomp(A = LPJ_CichlidMeanShapes_Tr, phy = FTree, align.to.phy = F)
pc.plot<-plot(LPJPCAs, axis1 = 1, axis2 = 2, phylo = F, col=ColourVec[unclass(gp.Diet)], pch=ShapeVec[unclass(gp.Diet)], flip = T)
shapeHulls(pc.plot, groups = as.character(gp.Diet), group.cols = c("#4daf4a", "#377eb8", "#ff7f00", "#984ea3","#a65628", "#e41a1c", "#ffff33"))
pc.means1 <- aggregate(LPJPCAs$x~gp.Diet, FUN=mean)
rownames(pc.means1) <- pc.means1[,1]
pc.means1 <- pc.means1[,-1]
points(pc.means1[,1:2], asp=1, pch=24, bg=ColourVec, cex=3)
text(LPJPCAs$x, pos = 4, label = dimnames(LPJ_CichlidMeanShapes_Tr)[[3]])

#write.csv(LPJPCAs$x[,1:5], "LPJ_PCscores.csv")
#write.csv(MandiblePCAs$x[,1:5], "LOJ_PCscores.csv")


####Disparity - Malawi vs. Tanganyika####
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Macroevolution")

LPJ_CichlidMeanShapes<-readland.tps(file = "MacroLPJ_Mean_20210826.tps", specID = "ID")
LOJ_CichlidMeanShapes<-readland.tps(file = "MacroLOJ_Mean_20210826.tps", specID = "ID")

TaxaIDData<-read.csv("TaxaData5.csv", row.names = 1)
Tree<-read.tree("CichlidTreeNature2020.tre")

LT_LM<-TaxaIDData[TaxaIDData$Lake=="T" | TaxaIDData$Lake=="M",]
Pruning<-treedata(Tree, LT_LM)

Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid<-dimnames(LPJ_CichlidMeanShapes)[[3]][dimnames(LPJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes[,,RM_Cichlid]

RM_Cichlid<-dimnames(LOJ_CichlidMeanShapes)[[3]][dimnames(LOJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes[,,RM_Cichlid]

LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LPJ_CichlidMeanShapes_Tr)[[3]]))]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LOJ_CichlidMeanShapes_Tr)[[3]]))]


gp.Lake <- as.factor(FData$Lake)
#LOJ
LOJfit_Lake<-procD.pgls(f1 = LOJ_CichlidMeanShapes_Tr~gp.Lake, phy = FTree)
morphol.disparity(f1 = LOJfit_Lake)
#M vs. T: p=0.253; M proc var = 0.07404536; T proc var = 0.05672349

#LPJ
LPJfit_Lake<-procD.pgls(f1 = LPJ_CichlidMeanShapes_Tr~gp.Lake, phy = FTree)
morphol.disparity(f1 = LPJfit_Lake)
#M vs. T: p=0.012; M proc var = 0.01529984; T proc var = 0.02262777*


####Compare Disparity and Integration####
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Macroevolution")

#Set-up for Lake Malawi
LPJ_CichlidMeanShapes<-readland.tps(file = "MacroLPJ_Mean_20210826.tps", specID = "ID")
LOJ_CichlidMeanShapes<-readland.tps(file = "MacroLOJ_Mean_20210826.tps", specID = "ID")

TaxaIDData<-read.csv("TaxaData5.csv", row.names = 1)
Tree<-read.tree("CichlidTreeNature2020.tre")

LMal<-TaxaIDData[TaxaIDData$Lake=="M",]

Pruning<-treedata(Tree, LMal)

Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid<-dimnames(LPJ_CichlidMeanShapes)[[3]][dimnames(LPJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes[,,RM_Cichlid]

RM_Cichlid<-dimnames(LOJ_CichlidMeanShapes)[[3]][dimnames(LOJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes[,,RM_Cichlid]

LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LPJ_CichlidMeanShapes_Tr)[[3]]))]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LOJ_CichlidMeanShapes_Tr)[[3]]))]


LMalPI<-phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree) #r-PLS: 0.556; P-value: 0.011; Effect Size: 2.2174
plot(phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree), label = dimnames(LOJ_CichlidMeanShapes_Tr)[[3]])


#Set-up for Lake Tanganyika
LPJ_CichlidMeanShapes<-readland.tps(file = "MacroLPJ_Mean_20210826.tps", specID = "ID")
LOJ_CichlidMeanShapes<-readland.tps(file = "MacroLOJ_Mean_20210826.tps", specID = "ID")

TaxaIDData<-read.csv("TaxaData5.csv", row.names = 1)
Tree<-read.tree("CichlidTreeNature2020.tre")

LTan<-TaxaIDData[TaxaIDData$Lake=="T",]

Pruning<-treedata(Tree, LTan)

Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid<-dimnames(LPJ_CichlidMeanShapes)[[3]][dimnames(LPJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes[,,RM_Cichlid]

RM_Cichlid<-dimnames(LOJ_CichlidMeanShapes)[[3]][dimnames(LOJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes[,,RM_Cichlid]

LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LPJ_CichlidMeanShapes_Tr)[[3]]))]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LOJ_CichlidMeanShapes_Tr)[[3]]))]


LTanPI<-phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree) #r-PLS: 0.682; P-value: 0.003; Effect Size: 2.851
plot(phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree), label = dimnames(LOJ_CichlidMeanShapes_Tr)[[3]])


#Compare PLS
compare.pls(LMalPI,LTanPI)
