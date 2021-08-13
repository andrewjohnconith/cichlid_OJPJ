library(qtl); library(plyr); library(plotrix)
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/QTL/FineMap/LPJ/PC1/")
LG7 <- read.cross(format="csv",file="LG7FULLPhenoMapEd_noNA.csv",
                         na.strings="NA",genotypes=c("AA","AB","BB"),
                         alleles=c("A","B"),convertXdata=TRUE)

#Check sample sizes
apply(LG7$geno$`7`$data, 2, table)
  #

EffResult<-matrix(data = NA, nrow = 132, ncol = 6)
#ScafNames<-names(LG7$geno$`7`$map)
ScafNames<-rownames(geno.table(LG7))[as.logical(apply(X = geno.table(LG7,7)[,3:5]>=20,MARGIN = 1,FUN = mean)==1)]
         
for (i in 1:length(ScafNames)){
  effect<-effectplot(LG7, pheno.col=2, mname1=ScafNames[i], var.flag=c("group"), draw = F) #effect plot at qtl
  EffResult[i,]<-c(effect$Means, effect$SEs)
}


CombinedDataEffectE<-cbind.data.frame(ScafNames,as.numeric(sub("^[^_]*_", "", ScafNames)), EffResult,EffResult[,3]-EffResult[,1])
colnames(CombinedDataEffectE)<-c("Scaff", "Pos", "AA_mean", "AB_mean", "BB_mean", "AA_se", "AB_se", "BB_se", "AA.BB")

CombinedDataEffectEo<-CombinedDataEffectE[order(CombinedDataEffectE$Pos),]

#write.csv(CombinedDataEffectEo,"LG7FineMapDataT_LPJpc1_Stringent_noaug.csv", quote = F, row.names = F)


##Plot
#setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/QTL/New/LPJ/FineMap/PC1/NewLG/")
#CombinedDataEffectEo<-read.csv("LG7FineMapDataT_LPJpc1_Stringent_noaug.csv")

MapPP<-CombinedDataEffectEo
MapP<-MapPP[complete.cases(MapPP$AA.BB) & complete.cases(MapPP$AA_se) & complete.cases(MapPP$BB_se),]

par(pty='s')
plot(MapP$Pos, MapP$AA.BB+MapP$AA_se, ylim = range(c(MapP$AA.BB,MapP$AA.BB+MapP$AA_se,MapP$AA.BB+MapP$BB_se,MapP$AA.BB-MapP$AA_se,MapP$AA.BB-MapP$BB_se), na.rm = T), type='n', lwd=1, lty=2, xlab = 'LG7', ylab = 'Avg. Pheno. Effect')
polygon(c(MapP$Pos,rev(MapP$Pos)),c(MapP$AA.BB+MapP$AA_se,rev(MapP$AA.BB-MapP$AA_se)), col = "purple", border = "purple", lwd=2)
polygon(c(MapP$Pos,rev(MapP$Pos)),c(MapP$AA.BB+MapP$BB_se,rev(MapP$AA.BB-MapP$BB_se)), col = "dark blue", border = "dark blue", lwd=2)
lines(MapP$Pos, MapP$AA.BB, col="light blue", lwd=1.5)

abline(h=0)
