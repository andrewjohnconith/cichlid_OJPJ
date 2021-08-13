library(geomorph); library(abind); library(plyr)

#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/MacroevolutionTropheopsSetUp")

#Read in data and TPS file
NaturalInds<-readland.tps("NaturalPopulation_LOJ_20201103.tps", specID = "imageID", readcurves = F)
NaturalInds[,1,]<-NaturalInds[,1,]*-1

#Procrustes, there's a few extra specimens in there, this removes them
Y.gpa<-gpagen(NaturalInds[,,NaturalData$TPS_ID[NaturalData$Genus=="Tropheops"]])

####Allometry####
gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.lm(coords~Csize, data = gdf, iter = 999, RRPP=TRUE) 
summary(HybAnova) 

shape.resid <- arrayspecs(HybAnova$residuals, p=dim(Y.gpa$coords)[1], k=dim(Y.gpa$coords)[2]) # size-adjusted residuals

#allometry-free shapes
MandibleShape <- shape.resid + array(Y.gpa$consensus, dim(shape.resid))

#Set up symmetric structure
ManSym<-cbind(3:11,12:20)

#Analyze asymmetry
gdf <- geomorph.data.frame(shape = MandibleShape, ind=dimnames(MandibleShape)[[3]])
Mandible.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
                             land.pairs=ManSym, data = gdf, RRPP = TRUE, iter = 999)
summary(Mandible.sym)

#Plot
par(pty='s')
#Sym
MandiblePCAs<-gm.prcomp(A = Mandible.sym$symm.shape)

plot(MandiblePCAs, axis1 = 1, axis2 = 2, pch=19)
text(MandiblePCAs$x, pos = 4, label = dimnames(MandibleShape)[[3]])
picknplot.shape(plot(MandiblePCAs))

#Save PC Scores and allometric/symmetric corrected tps if needed
#writeland.tps(A = Mandible.sym$symm.shape, file = "LOJ_20210719.tps")
#write.csv(MandiblePCAs$x, "ManPCScores_S.csv")


####Species Means####
NaturalData<-read.csv("TropheopsPopulation_LOJ_Data_20201214.csv")
#Coordinates#
NaturalData$SpeciesLocID<-as.factor(NaturalData$SpeciesLocID)
#Create mean GPA from all specimens
CichlidMeanShapes<-vector("list", length(levels(NaturalData$SpeciesLocID)))

for(i in 1:length(levels(NaturalData$SpeciesLocID))){
  
  CichlidNames <- levels(NaturalData$SpeciesLocID)[i]
  NamesRows <- which(NaturalData$SpeciesLocID == CichlidNames)
  
  CatTPS <- Mandible.sym$symm.shape[,,paste0(as.character(NaturalData[NamesRows,"TPS_ID"]))]
  
  #If you have a single entry need to convert to array (R thinks it is a matrix)
  if(class(CatTPS) != "array"){
    CatTPS<-array(CatTPS, dim=c(20,3,1))
  }
  
  Cichlid.gpa<-gpagen(CatTPS)
  CichlidMeanShapes[[i]] <- mshape(Cichlid.gpa$coords)
  
}

#Convert from list to array
CichlidMeanShapes <-array(unlist(CichlidMeanShapes), dim=c(20,3,length(levels(NaturalData$SpeciesLocID))))

#Add species names
CichlidDimName<-vector(mode="list", 3)
CichlidDimName[[3]]<-levels(NaturalData$SpeciesLocID)
dimnames(CichlidMeanShapes)<-CichlidDimName
CichlidMeanShapesY<-gpagen(CichlidMeanShapes)

#writeland.tps(A = CichlidMeanShapesY$coords, file = "NaturalPopulation_LOJ_MeanData_20201214.tps")
