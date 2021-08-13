library(geomorph); library(abind)

#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/HybridsSetUp")

#Read in data and TPS file
Hybrids<-readland.tps("LOJ_20200723.tps", specID = "ID", readcurves = F)

#Procrustes
Y.gpa<-gpagen(Hybrids) #

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

##Plot
par(pty='s')
#Sym
MandiblePCAs<-gm.prcomp(A = Mandible.sym$symm.shape)

plot(MandiblePCAs, axis1 = 1, axis2 = 2, pch=19)
text(MandiblePCAs$x, pos = 4, label = dimnames(MandibleShape)[[3]])
picknplot.shape(plot(MandiblePCAs))

#Save PC Scores
#writeland.tps(A = Mandible.sym$symm.shape, file = "LOJ_20200719.tps")
#write.csv(MandiblePCAs$x[,1:2], "ManPCScores_S.csv")
