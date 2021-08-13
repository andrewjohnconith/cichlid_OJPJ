library(geomorph); library(abind); library(plyr)

#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/MacroevolutionTropheopsSetUp")

#Read in data and TPS file
LPJInds<-readland.tps("NaturalPopulation_LPJ_20201103.tps", specID = "imageID", readcurves = F)
NaturalData<-read.csv("NaturalPopulation_LPJ_Data_20201103.csv")

#Procrustes, there's a few extra specimens in there, this removes them
Y.gpa<-gpagen(LPJInds[,,NaturalData$TPS_ID[NaturalData$Genus=="Tropheops"]])

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

#Save PC Scores and allometric/symmetric corrected tps if needed
#writeland.tps(A = LPJ.sym$symm.shape, file = "LPJ_20210719.tps")
#write.csv(LPJPCAs$x, "LPJPCScores_S.csv")


####Species Means####
NaturalData<-read.csv("TropheopsPopulation_LPJ_Data_20201214.csv")
#Coordinates#
NaturalData$SpeciesLocID<-as.factor(NaturalData$SpeciesLocID)
#Create mean GPA from all specimens
CichlidMeanShapes<-vector("list", length(levels(NaturalData$SpeciesLocID)))

for(i in 1:length(levels(NaturalData$SpeciesLocID))){
  
  CichlidNames <- levels(NaturalData$SpeciesLocID)[i]
  NamesRows <- which(NaturalData$SpeciesLocID == CichlidNames)
  
  CatTPS <- LPJ.sym$symm.shape[,,paste0(as.character(NaturalData[NamesRows,"TPS_ID"]))]
  
  #If you have a single entry need to convert to array (R thinks it is a matrix)
  if(class(CatTPS) != "array"){
    CatTPS<-array(CatTPS, dim=c(15,3,1))
  }
  
  Cichlid.gpa<-gpagen(CatTPS)
  CichlidMeanShapes[[i]] <- mshape(Cichlid.gpa$coords)
  
}

#Convert from list to array
CichlidMeanShapes <-array(unlist(CichlidMeanShapes), dim=c(15,3,length(levels(NaturalData$SpeciesLocID))))

#Add species names
CichlidDimName<-vector(mode="list", 3)
CichlidDimName[[3]]<-levels(NaturalData$SpeciesLocID)
dimnames(CichlidMeanShapes)<-CichlidDimName
CichlidMeanShapesY<-gpagen(CichlidMeanShapes)

#writeland.tps(A = CichlidMeanShapesY$coords, file = "NaturalPopulation_LPJ_MeanData_20201214.tps")
