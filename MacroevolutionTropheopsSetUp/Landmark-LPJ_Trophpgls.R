library(geomorph); library(abind); library(plyr); library(geiger)

#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/MacroevolutionTropheopsSetUp")

#Read in data and TPS file
LPJInds<-readland.tps("NaturalPopulation_LPJ_20201103.tps", specID = "imageID", readcurves = F)
NaturalData<-read.csv("NaturalPopulation_LPJ_Data_20201103.csv")

#Procrustes, there's a few extra specimens in there, this removes them
Y.gpa<-gpagen(LPJInds[,,NaturalData$TPS_ID[NaturalData$Genus=="Tropheops"]])


####Species Means####
NaturalData<-read.csv("TropheopsPopulation_LPJ_Data_20201214.csv")
#Coordinates#
NaturalData$SpeciesLocID<-as.factor(NaturalData$SpeciesLocID)
#Create mean GPA from all specimens
CichlidMeanShapes<-vector("list", length(levels(NaturalData$SpeciesLocID)))

for(i in 1:length(levels(NaturalData$SpeciesLocID))){
  
  CichlidNames <- levels(NaturalData$SpeciesLocID)[i]
  NamesRows <- which(NaturalData$SpeciesLocID == CichlidNames)
  
  CatTPS <- Y.gpa$coords[,,paste0(as.character(NaturalData[NamesRows,"TPS_ID"]))]
  
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

morphmean <- aggregate(Csize~SpeciesLocID,data=NaturalData,FUN="mean",na.rm=TRUE,na.action=NULL)
colnames(morphmean)<-c("PhyloID","Csize")
CichlidMeanShapesY$Csize<-morphmean$Csize

####Allometry####
LPJ_CichlidMeanShapes<-CichlidMeanShapesY$coords

TaxaIDData<-read.csv("TaxaData.csv")
Tree<-read.tree("Tree-CompMeth-Ultra_20181101.tre")

TaxaIDData<-merge(TaxaIDData, morphmean, by = "PhyloID")
rownames(TaxaIDData)<-TaxaIDData$PhyloID

Pruning<-treedata(Tree, TaxaIDData)
Ftree<-Pruning$phy
FTree<-ladderize(Ftree)
FData<-Pruning$data

FData <- FData[match(FTree$tip.label, rownames(FData)),]
FData<-as.data.frame(FData)

RM_Cichlid<-dimnames(LPJ_CichlidMeanShapes)[[3]][dimnames(LPJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes[,,RM_Cichlid]

LPJ_CichlidMeanShapes_Tr<-LPJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LPJ_CichlidMeanShapes_Tr)[[3]]))]

gdf <- geomorph.data.frame(coords = LPJ_CichlidMeanShapes_Tr, species = dimnames(LPJ_CichlidMeanShapes_Tr)[[3]], Csize = as.numeric(FData$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.pgls(coords~Csize, phy = FTree, data = gdf) 
summary(HybAnova) 

shape.resid.LPJ <- arrayspecs(HybAnova$pgls.residuals, p=dim(LPJ_CichlidMeanShapes_Tr)[1], k=dim(LPJ_CichlidMeanShapes_Tr)[2]) # size-adjusted residuals
LPJShape <- shape.resid.LPJ + array(CichlidMeanShapesY$consensus, dim(shape.resid.LPJ))

####Symmetric####
#Set up symmetric structure
LPJSym<-cbind(6:10,11:15)

#Analyze asymmetry
gdf <- geomorph.data.frame(shape = LPJShape, ind=dimnames(LPJShape)[[3]])
LPJ.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
                             land.pairs=LPJSym, data = gdf, RRPP = TRUE, iter = 999)
summary(LPJ.sym)


#writeland.tps(A = LPJ.sym$symm.shape, file = "NaturalPopulation_LPJ_MeanData_20210826.tps")
