library(geomorph)

#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Hybrids/")

####LPJ Hybrid####
#Read in allo/sym corrected TPS file
LPJ.adj.shape<-readland.tps("LPJ_20200719.tps", specID = "ID")

####LOJ Hybrid####
#Read in allo/sym corrected TPS file
LOJ.adj.shape<-readland.tps("LOJ_20200719.tps", specID = "ID")

####Covariation Analysis####
#Common Individuals
CommonInds<-c(dimnames(LPJ.adj.shape)[[3]], dimnames(LOJ.adj.shape)[[3]])[duplicated(c(dimnames(LPJ.adj.shape)[[3]], dimnames(LOJ.adj.shape)[[3]]))]

#LPJ setup
LPJnew<-LPJ.adj.shape[,,CommonInds]

#LOJ setup
LOJnew<-LOJ.adj.shape[,,CommonInds]

LOJLPJCorrs<-two.b.pls(A1 = LOJnew, A2 = LPJnew) #r-PLS: 0.491; P-value: 0.001; Effect Size: 6.1886
plot(LOJLPJCorrs, label = CommonInds)
plot(LOJLPJCorrs, pch=19)
summary(LOJLPJCorrs)


####Correlation####
#X-Y Coords
LOJLPJCorrsH<-cbind.data.frame(LOJLPJCorrs$XScores[,1],LOJLPJCorrs$YScores[,1])
colnames(LOJLPJCorrsH)<-c("XBlock", "YBlock")

##Extracts the coefficients directly from the PLS analysis
#Run function below to extract the coefficients
#PLEASEWORK is a global assignment
MYPLSFUNCTIONT(LOJLPJCorrs)
PLSCoefficients<-pPLS_Info

MyRatio <- with(LOJLPJCorrsH, diff(range(XBlock))/diff(range(YBlock)))
Hybplot <- ggplot(LOJLPJCorrsH, aes(x=XBlock, y=YBlock)) +
  geom_point(size = 2.5, color="#636363", alpha=0.8) +
  labs(x = "L. Oral Jaw Block", y = "L. Pharyngeal Jaw Block") +
  theme_bw() + coord_fixed(ratio=MyRatio)

Hybplot + geom_abline(intercept = PLSCoefficients$coefficients[1], slope = PLSCoefficients$coefficients[2]) +
  theme(legend.key=element_blank(), axis.ticks = element_blank(), panel.border = element_rect(colour = "dark gray", fill=NA, size=1))

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
