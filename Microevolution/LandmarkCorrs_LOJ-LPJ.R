#install.packages("/Users/home/Dropbox/Amherst/Courses/SOURCE/Undergraduate/PostDoc/Week3/geomorph_3.0.4.tar.gz", repos = NULL, type="source")
library(geomorph); library(abind)

####LOJ####
#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Microevolution/")

#Read in data and TPS file
LOJp<-readland.tps("LOJ_p_20201105.tps", specID = "imageID", readcurves = F)
LOJp[,1,]<-LOJp[,1,]*-1

#Procrustes
Y.gpa<-gpagen(LOJp[c(1:11,20:12),,]) #


###Allometry###
gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.lm(coords~Csize, data = gdf, iter = 999, RRPP=TRUE) 
summary(HybAnova) 

shape.resid <- arrayspecs(HybAnova$residuals, p=dim(Y.gpa$coords)[1], k=dim(Y.gpa$coords)[2]) # size-adjusted residuals

#allometry-free shapes
MandibleShape <- shape.resid + array(Y.gpa$consensus, dim(shape.resid))

#Set up symmetric structure
#ManSym<-cbind(3:12,22:13)
ManSym<-cbind(3:11,20:12)

#Analyze asymmetry
gdf <- geomorph.data.frame(shape = MandibleShape, ind=dimnames(MandibleShape)[[3]])
Mandible.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
                               land.pairs=ManSym, data = gdf, RRPP = TRUE, iter = 999)
summary(Mandible.sym)

PGroup <- c(rep("blue",17),rep("red",11))

####LPJ####
#Read in data and TPS file
LPJp<-readland.tps("LPJ_p_20201105.tps", specID = "imageID", readcurves = F)

#Procrustes
Y.gpa<-gpagen(LPJp)

###Allometry###
#Y.gpa$Csize<-AllData$centroid.size[c(5:144)] #HBs
gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.lm(coords~Csize, data = gdf, iter = 999, RRPP=TRUE) 
summary(HybAnova) 

shape.resid <- arrayspecs(HybAnova$residuals, p=dim(Y.gpa$coords)[1], k=dim(Y.gpa$coords)[2]) # size-adjusted residuals

#allometry-free shapes
LPJShape <- shape.resid + array(Y.gpa$consensus, dim(shape.resid))

#Set up symmetric structure
LPJSym<-cbind(6:10,11:15)

#Analyze asymmetry
gdf <- geomorph.data.frame(shape = LPJShape, ind=dimnames(LPJShape)[[3]])
LPJ.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
                          land.pairs=LPJSym, data = gdf, RRPP = TRUE, iter = 999)
summary(LPJ.sym)

PGroup <- c(rep("blue",17),rep("red",11))


####Plotting####
#LOJ
par(pty='s', mfrow=c(1,2))
LOJPCAs<-gm.prcomp(A = Mandible.sym$symm.shape)
plot(LOJPCAs, axis1 = 1, axis2 = 2, pch=19, col=PGroup)
plot(LOJPCAs, axis1 = 1, axis2 = 3, pch=19, col=PGroup)
#text(LOJPCAs$x, pos = 4, label = dimnames(Mandible.sym$symm.shape)[[3]])

#write.csv(LOJPCAs$x, "LOJPCAsParentals.csv")

#LPJ
LPJPCAs<-gm.prcomp(A = LPJ.sym$symm.shape)
plot(LPJPCAs, axis1 = 1, axis2 = 2, pch=19, col=PGroup)
plot(LPJPCAs, axis1 = 1, axis2 = 3, pch=19, col=PGroup)
#text(LPJPCAs$x, pos = 4, label = dimnames(Mandible.sym$symm.shape)[[3]])

#write.csv(LPJPCAs$x, "LPJPCAsParentals.csv")


####Correlation####
par(pty='s', mfrow=c(1,1))
##All
All_Corr<-two.b.pls(A1 = Mandible.sym$symm.shape, A2 = LPJ.sym$symm.shape)
plot(All_Corr)

#LF
LF_Corr<-two.b.pls(A1 = Mandible.sym$symm.shape[,,1:17], A2 = LPJ.sym$symm.shape[,,1:17], seed = 12345)
plot(LF_Corr)

#TRC
TRC_Corr<-two.b.pls(A1 = Mandible.sym$symm.shape[,,18:28], A2 = LPJ.sym$symm.shape[,,18:28],  seed = 12345)
plot(TRC_Corr)

#Compare integration between LOJ and LPJ
compare.pls(LF = LF_Corr, TRC = TRC_Corr)
?compare.pls
#Effect sizes
#  LF_Corr  TRC_Corr 
#0.7936856 3.0380033 

#Effect sizes for pairwise differences in PLS effect size
#          LF_Corr TRC_Corr
#LF_Corr  0.000000 1.678421
#TRC_Corr 1.678421 0.000000

#P-values
#            LF_Corr   TRC_Corr
#LF_Corr  1.00000000 0.04663243
#TRC_Corr 0.04663243 1.00000000

#TRC is more integrated than LF


####Plotting####
library(ggplot2)
##LF##
#X-Y Coords
LF_LOJLPJ<-cbind.data.frame(LF_Corr$XScores[,1],LF_Corr$YScores[,1])
colnames(LF_LOJLPJ)<-c("XBlock", "YBlock")

##Extracts the coefficients directly from the PLS analysis
#Run function below to extract the coefficients
#PLEASEWORK is a global assignment
MYPLSFUNCTIONT(LF_Corr)
PLSCoefficients<-pPLS_Info

MyRatio <- with(LF_LOJLPJ, diff(range(XBlock))/diff(range(YBlock)))
LFplot <- ggplot(LF_LOJLPJ, aes(x=XBlock, y=YBlock)) +
  geom_point(size = 2.5, color="#0571b0") +
  labs(x = "L. Oral Jaw Block", y = "L. Pharyngeal Jaw Block") +
  theme_bw() + coord_fixed(ratio=MyRatio)

LFplotFULL <- LFplot + geom_abline(intercept = PLSCoefficients$coefficients[1], slope = PLSCoefficients$coefficients[2]) +
  theme(legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))


##TRC##
#X-Y Coords
TRC_LOJLPJ<-cbind.data.frame(TRC_Corr$XScores[,1],TRC_Corr$YScores[,1])
colnames(TRC_LOJLPJ)<-c("XBlock", "YBlock")

##Extracts the coefficients directly from the PLS analysis
#Run function below to extract the coefficients
#PLEASEWORK is a global assignment
MYPLSFUNCTIONT(TRC_Corr)
PLSCoefficients<-pPLS_Info

MyRatio <- with(TRC_LOJLPJ, diff(range(XBlock))/diff(range(YBlock)))
TRCplot <- ggplot(TRC_LOJLPJ, aes(x=XBlock, y=YBlock)) +
  geom_point(size = 2.5, color="#ca0020") +
  labs(x = "L. Oral Jaw Block", y = "L. Pharyngeal Jaw Block") +
  theme_bw() + coord_fixed(ratio=MyRatio)

TRCplotFULL <- TRCplot + geom_abline(intercept = PLSCoefficients$coefficients[1], slope = PLSCoefficients$coefficients[2]) +
  theme(legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

library("cowplot")
plot_grid(TRCplotFULL, LFplotFULL)


####PLS Extract####
#Extract out the PLS line with this function
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
