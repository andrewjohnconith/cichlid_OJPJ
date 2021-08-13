library(geomorph); library(abind)

#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/HybridsSetUp")

#Read in data and TPS file
LPJInds<-readland.tps("LPJ_20200723.tps", specID = "ID", readcurves = F)

#Procrustes
Y.gpa<-gpagen(LPJInds) #

####Allometry####
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

##Plot
par(pty='s')
#Sym
LPJPCAs<-gm.prcomp(A = LPJ.sym$symm.shape)

plot(LPJPCAs, axis1 = 1, axis2 = 2, pch=19)
text(LPJPCAs$x, pos = 4, label = dimnames(LPJShape)[[3]])
picknplot.shape(plot(LPJPCAs))

#Save PC Scores
#writeland.tps(A = LPJ.sym$symm.shape, file = "LPJ_20200719.tps")
#write.csv(LPJPCAs$x, "LPJPCScores_S.csv")
