library(geomorph); library(abind); library(plyr); library(geiger)

#Set WD
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/MacroevolutionTropheopsSetUp")

#Read in data and TPS file
NaturalInds<-readland.tps("NaturalPopulation_LOJ_20201103.tps", specID = "imageID", readcurves = F)
NaturalInds[,1,]<-NaturalInds[,1,]*-1

NaturalData<-read.csv("NaturalPopulation_LOJ_Data_20201103.csv")

#Procrustes, there's a few extra specimens in there, this removes them
Y.gpa<-gpagen(NaturalInds[,,NaturalData$TPS_ID[NaturalData$Genus=="Tropheops"]])


####Species Means####
NaturalData<-read.csv("TropheopsPopulation_LOJ_Data_20201214.csv")
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

morphmean <- aggregate(Csize~SpeciesLocID,data=NaturalData,FUN="mean",na.rm=TRUE,na.action=NULL)
colnames(morphmean)<-c("PhyloID","Csize")
CichlidMeanShapesY$Csize<-morphmean$Csize

####Allometry####
LOJ_CichlidMeanShapes<-CichlidMeanShapesY$coords

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

RM_Cichlid<-dimnames(LOJ_CichlidMeanShapes)[[3]][dimnames(LOJ_CichlidMeanShapes)[[3]]%in%rownames(FData)]
LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes[,,RM_Cichlid]

LOJ_CichlidMeanShapes_Tr<-LOJ_CichlidMeanShapes_Tr[,,(match(FTree$tip.label, dimnames(LOJ_CichlidMeanShapes_Tr)[[3]]))]

gdf <- geomorph.data.frame(coords = LOJ_CichlidMeanShapes_Tr, species = dimnames(LOJ_CichlidMeanShapes_Tr)[[3]], Csize = as.numeric(FData$Csize))

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.pgls(coords~Csize, phy = FTree, data = gdf) 
summary(HybAnova) 

shape.resid.LOJ <- arrayspecs(HybAnova$pgls.residuals, p=dim(LOJ_CichlidMeanShapes_Tr)[1], k=dim(LOJ_CichlidMeanShapes_Tr)[2]) # size-adjusted residuals
LOJShape <- shape.resid.LOJ + array(CichlidMeanShapesY$consensus, dim(shape.resid.LOJ))


####Symmetric####
#Set up symmetric structure
ManSym<-cbind(3:11,12:20)

#Analyze asymmetry
gdf <- geomorph.data.frame(shape = LOJShape, ind=dimnames(LOJShape)[[3]])
Mandible.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
                             land.pairs=ManSym, data = gdf, RRPP = TRUE, iter = 999)
summary(Mandible.sym)


#writeland.tps(A = Mandible.sym$symm.shape, file = "NaturalPopulation_LOJ_MeanData_20210826.tps")
