#install.packages("/Users/home/Dropbox/Amherst/Courses/SOURCE/Undergraduate/PostDoc/Week3/geomorph_3.0.4.tar.gz", repos = NULL, type="source")
library(geomorph); library(abind)

#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/HybridsSetUp")

#Read in data and TPS file
HeartInds<-readland.tps("LOJa_20190929.tps", specID = "imageID", readcurves = F)

#Reflect specific individuals
RefInd<-c("LFxTRC5620", "LFxTRC5614", "LFxTRC5613", "LFxTRC5606", "LFxTRC5603", "LFxTRC5601", "LFxTRC5595", "LFxTRC5591", "LFxTRC5589", "LFxTRC5583", "LFxTRC5582", "LFxTRC5578", "LFxTRC5569", "LFxTRC5563", "LFxTRC5552", "LFxTRC5547", "LFxTRC5527", "LFxTRC5526", "LFxTRC5510", "LFxTRC5509", "LFxTRC5115", "LFxTRC5113", "LFxTRC5110", "LFxTRC5108", "LFxTRC5098", "LFxTRC5095", "LFxTRC5087", "LFxTRC5086", "LFxTRC5009", "TRC17-0031", "LF17-N3")
HeartInds[,1,RefInd]<-HeartInds[,1,RefInd]*-1

LOJInds<-readland.tps("LOJb_20200619.tps", specID = "imageID", readcurves = F)
LOJInds[,1,]<-LOJInds[,1,]*-1

#Need to reverse the right  side of new individuals
AllInds<-abind(HeartInds,LOJInds[c(1:11,20:12),,],along = 3)

#Procrustes
Y.gpa<-gpagen(AllInds[,,c(-1:-4,-147:-150)]) #remove parentals

####Allometry####
gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.lm(coords~Csize, data = gdf, iter = 999, RRPP=TRUE) 
summary(HybAnova) 

shape.resid <- arrayspecs(HybAnova$residuals, p=dim(Y.gpa$coords)[1], k=dim(Y.gpa$coords)[2]) # size-adjusted residuals

#allometry-free shapes
MandibleShape <- shape.resid + array(Y.gpa$consensus, dim(shape.resid))

#Set up symmetric structure
ManSym<-cbind(3:11,20:12)

#Analyze asymmetry
gdf <- geomorph.data.frame(shape = MandibleShape, ind=dimnames(MandibleShape)[[3]])
Mandible.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
                             land.pairs=ManSym, data = gdf, RRPP = TRUE, iter = 999)
summary(Mandible.sym)

#Sym
MandiblePCAs<-gm.prcomp(A = Mandible.sym$symm.shape)

plot(MandiblePCAs, axis1 = 1, axis2 = 2, pch=19)
text(MandiblePCAs$x, pos = 4, label = dimnames(MandibleShape)[[3]])
picknplot.shape(plot(MandiblePCAs))

#writeland.tps(A = Mandible.sym$symm.shape, file = "LOJ_20200719.tps")

#Save PC Scores
#write.csv(MandiblePCAs$x, "ManPCScores_S.csv")
