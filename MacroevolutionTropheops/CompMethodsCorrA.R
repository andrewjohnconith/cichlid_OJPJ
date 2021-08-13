library(geomorph); library(geiger); library(plotrix); library(caper); library(ggplot2)

####Data Engineering####
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/MacroevolutionTropheops/")

LPJ_CichlidMeanShapes<-readland.tps(file = "NaturalPopulation_LPJ_MeanData_20201214.tps", specID = "ID")
LOJ_CichlidMeanShapes<-readland.tps(file = "NaturalPopulation_LOJ_MeanData_20201214.tps", specID = "ID")

TaxaIDData<-read.csv("TaxaData.csv", row.names = 1)
Tree<-read.tree("Tree-CompMeth-Ultra_20181101.tre")

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

gp.depth <- as.factor(FData$Depth)


####Phylo-Covariation Analysis####
phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree) #r-PLS: 0.802; P-value: 0.004; Effect Size: 2.5716
plot(phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree), label = gp.depth)

#X-Y Coords
LOJLPJcoords<-phylo.integration(A = LOJ_CichlidMeanShapes_Tr, A2 = LPJ_CichlidMeanShapes_Tr, phy = FTree)

LOJLPJ<-cbind.data.frame(LOJLPJcoords$XScores[,1],LOJLPJcoords$YScores[,1], FData$Depth)
colnames(LOJLPJ)<-c("XBlock", "YBlock", "Depth")
#write.csv(LOJLPJ, "CorrDataTroph.csv")

#Plot
MyRatio <- with(LOJLPJ, diff(range(XBlock))/diff(range(YBlock)))
r <- ggplot(LOJLPJ, aes(x=XBlock, y=YBlock)) +
   geom_point(aes(colour = factor(Depth), shape = factor(Depth)), size = 2.5) +
   labs(colour = "Depth", x = "L. Oral Jaw Block", y = "L. Pharyngeal Jaw Block") +
   theme_bw() + coord_fixed(ratio=MyRatio)
r

####Extract pPLS Line####
##Extracts the coefficients directly from the pPLS analysis
#Run function below to extract the coefficients
#pPLS_Info is a global assignment
MYPLSFUNCTIONT(LOJLPJcoords)
pPLSCoefficients<-pPLS_Info

MyRatio <- with(LOJLPJ, diff(range(XBlock))/diff(range(YBlock)))
p <- ggplot(LOJLPJ, aes(x=XBlock, y=YBlock)) +
  geom_point(aes(colour = factor(Depth), shape = factor(Depth)), size = 2.5) +
  geom_abline(intercept = pPLSCoefficients$coefficients[1], slope = pPLSCoefficients$coefficients[2]) +
  scale_shape_manual(values = c(19, 17), name="Depth", labels=c("Deep","Shallow")) +
  scale_colour_manual(values=c("#F8766D","#00BFC4"), name="Depth", labels=c("Deep","Shallow")) +
  theme_bw() + coord_fixed(ratio=MyRatio)

p

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


