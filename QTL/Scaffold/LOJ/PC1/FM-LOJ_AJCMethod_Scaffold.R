library(qtl); library(plyr); library(plotrix)
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/QTL/Scaffold/LOJ/PC1/")
Scaffold0 <- read.cross(format="csv",file="Scaffold0FULLPhenoMapEdnoNA.csv",
                         na.strings="NA",genotypes=c("AA","AB","BB"),
                         alleles=c("A","B"),convertXdata=TRUE)

#Check sample sizes
apply(Scaffold0$geno$`7`$data, 2, table)
  #

#ScafNames<-names(Scaffold0$geno$`7`$map)
ScafNames<-rownames(geno.table(Scaffold0))[as.logical(apply(X = geno.table(Scaffold0,7)[,3:5]>=15,MARGIN = 1,FUN = mean)==1)]
EffResult<-matrix(data = NA, nrow = length(ScafNames), ncol = 6)

for (i in 1:length(ScafNames)){
  ScafID<-Scaffold0$geno$`7`$data[,ScafNames[i]]
  TestingD<-as.data.frame(cbind(Scaffold0$pheno$HybridID,ScafID,Scaffold0$pheno$PC1))
  colnames(TestingD)<-c("ID","factor","value")
  r2<-ddply(TestingD, .(factor), summarize, mean=mean(value))
  r3<-ddply(TestingD, .(factor), summarize, mean=std.error(value))
  
  EffResult[i,]<-c(r2$mean[r2$factor=="1"][1],r2$mean[r2$factor=="2"][1],r2$mean[r2$factor=="3"][1],r3$mean[r3$factor=="1"][1],r3$mean[r3$factor=="2"][1],r3$mean[r3$factor=="3"][1])
}

CombinedDataEffectE<-cbind.data.frame(ScafNames,as.numeric(sub("^[^_]*_[^_]*_", "", ScafNames)), EffResult,EffResult[,3]-EffResult[,1])
colnames(CombinedDataEffectE)<-c("Scaff", "Pos", "AA_mean", "AB_mean", "BB_mean", "AA_se", "AB_se", "BB_se", "AA.BB")

CombinedDataEffectEo<-CombinedDataEffectE[order(CombinedDataEffectE$Pos),]

#write.csv(CombinedDataEffectEo,"Scaffold0FineMapDataT_LOJpc1_Stringent_noaug.csv", quote = F, row.names = F)


##Plot
#setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/QTL/New/LOJ/FineMap/PC1/Scaffold0/")
#CombinedDataEffectEo<-read.csv("Scaffold0FineMapDataT_LOJpc1_Stringent_noaug.csv")

MapPP<-CombinedDataEffectEo
MapP<-MapPP[complete.cases(MapPP$AA.BB) & complete.cases(MapPP$AA_se) & complete.cases(MapPP$BB_se),]

MapP<-MapP[39:56,]

MapPr<-MapP[rev(order(MapP$Pos)), ]

MapPr$Pos<-(MapPr$Pos-3785116)*-1




setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/QTL/Scaffold/LOJ/PC1/")
Scaffold21 <- read.cross(format="csv",file="Scaffold21FULLPhenoMapEdnoNA.csv",
                         na.strings="NA",genotypes=c("AA","AB","BB"),
                         alleles=c("A","B"),convertXdata=TRUE)

#Check sample sizes
apply(Scaffold21$geno$`7`$data, 2, table)
#scaffold_21_1081954; 7-BB: scaffold_21_1556666; 7-AA
ScafNames<-rownames(geno.table(Scaffold21))[as.logical(apply(X = geno.table(Scaffold21,7)[,3:5]>=15,MARGIN = 1,FUN = mean)==1)]
#ScafNames<-names(Scaffold21$geno$`7`$map)
EffResult<-matrix(data = NA, nrow = length(ScafNames), ncol = 6)

for (i in 1:length(ScafNames)){ 
  ScafID<-Scaffold21$geno$`7`$data[,ScafNames[i]]
  TestingD<-as.data.frame(cbind(Scaffold21$pheno$HybridID,ScafID,Scaffold21$pheno$PC1))
  colnames(TestingD)<-c("ID","factor","value")
  r2<-ddply(TestingD, .(factor), summarize, mean=mean(value))
  r3<-ddply(TestingD, .(factor), summarize, mean=std.error(value))
  
  EffResult[i,]<-c(r2$mean[r2$factor=="1"][1],r2$mean[r2$factor=="2"][1],r2$mean[r2$factor=="3"][1],r3$mean[r3$factor=="1"][1],r3$mean[r3$factor=="2"][1],r3$mean[r3$factor=="3"][1])
}

CombinedDataEffectE<-cbind.data.frame(ScafNames,as.numeric(sub("^[^_]*_[^_]*_", "", ScafNames)), EffResult,EffResult[,3]-EffResult[,1])
colnames(CombinedDataEffectE)<-c("Scaff", "Pos", "AA_mean", "AB_mean", "BB_mean", "AA_se", "AB_se", "BB_se", "AA.BB")

CombinedDataEffectEo<-CombinedDataEffectE[order(CombinedDataEffectE$Pos),]

#write.csv(CombinedDataEffectEo,"Scaffold21FineMapDataT_LOJpc1_Stringent_noaug.csv", quote = F, row.names = F)


##Plot
#setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/QTL/New/LOJ/FineMap/PC1/Scaffold/")
#CombinedDataEffectEo<-read.csv("Scaffold21FineMapDataT_LOJpc1_Stringent_noaug.csv")

MapPP<-CombinedDataEffectEo
MapP<-MapPP[complete.cases(MapPP$AA.BB) & complete.cases(MapPP$AA_se) & complete.cases(MapPP$BB_se),]
MapP<-MapP[1:25,]

MapPPP<-rbind(MapPr,MapP)

par(pty='s', mfrow=c(1,1))
plot(MapPPP$Pos, abs(MapPPP$AA.BB), ylim = range(abs(MapPPP$AA.BB)), type='l', lwd=2, axes=T, xlab=NA, ylab=NA, col="red")
abline(h=0)
