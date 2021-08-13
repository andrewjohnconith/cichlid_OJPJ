#####Calculate ddCT####
install.packages('pcr')
devtools::install_github('MahShaaban/pcr')
library(pcr); library(ggplot2)

setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/Expression/SetUp")
FullData<-read.csv("RAWcq.csv")

####DYM####
###Tropheops###
SubsetData<-FullData[FullData$Species=="Tropheops" & FullData$Gene %in% c("DYM", "BACTIN"),]

SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")
SubsetDataAGsd<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "sd")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="DYM",][,4])
colnames(GeneCq)<-c("BACTIN","DYM")

#Split by individual
Ind1<-GeneCq[c(1,7,13),]
Ind2<-GeneCq[c(2,8,14),]
Ind3<-GeneCq[c(3,9,15),]
Ind4<-GeneCq[c(4,10,16),]
Ind5<-GeneCq[c(5,11,17),]
Ind6<-GeneCq[c(6,12,18),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

Tr_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllTropeopsNorm<-rbind.data.frame(
Tr_res1$relative_expression[2:3],
Tr_res2$relative_expression[2:3],
Tr_res3$relative_expression[2:3],
Tr_res4$relative_expression[2:3],
Tr_res5$relative_expression[2:3],
Tr_res6$relative_expression[2:3]
)

colnames(AllTropeopsNorm)<-c("OJ","PJ")

par(pty='s')
plot(AllTropeopsNorm, pch=19, xlim=c(0,2.5), ylim=c(0,2.5), col="red")


###Labeotropheus###
SubsetData<-FullData[FullData$Species=="Lfuelleborni" & FullData$Gene %in% c("DYM", "BACTIN"),]

##Save this next line as a CSV file
SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")
SubsetDataAGsd<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "sd")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="DYM",][,4])
colnames(GeneCq)<-c("BACTIN","DYM")

#Split by individual
Ind1<-GeneCq[c(1,9,17),]
Ind2<-GeneCq[c(2,10,18),]
Ind3<-GeneCq[c(3,11,19),]
Ind4<-GeneCq[c(4,12,20),]
Ind5<-GeneCq[c(5,13,21),]
Ind6<-GeneCq[c(6,14,22),]
Ind7<-GeneCq[c(7,15,23),]
Ind8<-GeneCq[c(8,16,24),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

LF_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res7 <- pcr_analyze(Ind7, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res8 <- pcr_analyze(Ind8, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllLabeotropheusNorm<-rbind.data.frame(
  LF_res1$relative_expression[2:3],
  LF_res2$relative_expression[2:3],
  LF_res3$relative_expression[2:3],
  LF_res4$relative_expression[2:3],
  LF_res5$relative_expression[2:3],
  LF_res6$relative_expression[2:3],
  LF_res7$relative_expression[2:3],
  LF_res8$relative_expression[2:3]
)

colnames(AllLabeotropheusNorm)<-c("OJ","PJ")

points(AllLabeotropheusNorm, pch=19, col="blue")


###Maylandia###
SubsetData<-FullData[FullData$Species=="Mcallainos" & FullData$Gene %in% c("DYM", "BACTIN"),]

##Save this next line as a CSV file
SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="DYM",][,4])
colnames(GeneCq)<-c("BACTIN","DYM")

#Split by individual
Ind1<-GeneCq[c(1,9,17),]
Ind2<-GeneCq[c(2,10,18),]
Ind3<-GeneCq[c(3,11,19),]
Ind4<-GeneCq[c(4,12,20),]
Ind5<-GeneCq[c(5,13,21),]
Ind6<-GeneCq[c(6,14,22),]
Ind7<-GeneCq[c(7,15,23),]
Ind8<-GeneCq[c(8,16,24),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

MC_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res7 <- pcr_analyze(Ind7, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res8 <- pcr_analyze(Ind8, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllMaylandiaNorm<-rbind.data.frame(
  MC_res1$relative_expression[2:3],
  MC_res2$relative_expression[2:3],
  MC_res3$relative_expression[2:3],
  MC_res4$relative_expression[2:3],
  MC_res5$relative_expression[2:3],
  MC_res6$relative_expression[2:3],
  MC_res7$relative_expression[2:3],
  MC_res8$relative_expression[2:3]
)

colnames(AllMaylandiaNorm)<-c("OJ","PJ")

points(AllMaylandiaNorm, pch=19, xlim=c(0,4), ylim=c(0,4), col="purple")


###Collate###
AllSpecies<-rbind.data.frame(
  AllMaylandiaNorm, AllTropeopsNorm, AllLabeotropheusNorm
  
)

TaxaID<-c(rep('Maylandia', 8),
          rep('Tropheops', 6),
          rep('Labeotropheus', 8))

TissueID<-c(rep("OJ",22),
            rep("PJ",22))

AllSpeciesData<-as.numeric(unlist(AllSpecies))

DYMFullExpression<-cbind.data.frame(TaxaID,TissueID,AllSpeciesData)



####NOTCH1A####
###Tropheops###
SubsetData<-FullData[FullData$Species=="Tropheops" & FullData$Gene %in% c("NOTCH1A", "BACTIN"),]

##Save this next line as a CSV file
SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")
SubsetDataAGsd<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "sd")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="NOTCH1A",][,4])
colnames(GeneCq)<-c("BACTIN","NOTCH1A")

#Split by individual
Ind1<-GeneCq[c(1,7,13),]
Ind2<-GeneCq[c(2,8,14),]
Ind3<-GeneCq[c(3,9,15),]
Ind4<-GeneCq[c(4,10,16),]
Ind5<-GeneCq[c(5,11,17),]
Ind6<-GeneCq[c(6,12,18),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

Tr_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllTropeopsNorm<-rbind.data.frame(
  Tr_res1$relative_expression[2:3],
  Tr_res2$relative_expression[2:3],
  Tr_res3$relative_expression[2:3],
  Tr_res4$relative_expression[2:3],
  Tr_res5$relative_expression[2:3],
  Tr_res6$relative_expression[2:3]
)

colnames(AllTropeopsNorm)<-c("OJ","PJ")

plot(AllTropeopsNorm, pch=19, xlim=c(0,3), ylim=c(0,3), col="red")


###Labeotropheus###
SubsetData<-FullData[FullData$Species=="Lfuelleborni" & FullData$Gene %in% c("NOTCH1A", "BACTIN"),]

##Save this next line as a CSV file
SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")
SubsetDataAGsd<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "sd")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="NOTCH1A",][,4])
colnames(GeneCq)<-c("BACTIN","NOTCH1A")

#Split by individual
Ind1<-GeneCq[c(1,9,17),]
Ind2<-GeneCq[c(2,10,18),]
Ind3<-GeneCq[c(3,11,19),]
Ind4<-GeneCq[c(4,12,20),]
Ind5<-GeneCq[c(5,13,21),]
Ind6<-GeneCq[c(6,14,22),]
Ind7<-GeneCq[c(7,15,23),]
Ind8<-GeneCq[c(8,16,24),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

LF_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res7 <- pcr_analyze(Ind7, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res8 <- pcr_analyze(Ind8, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllLabeotropheusNorm<-rbind.data.frame(
  LF_res1$relative_expression[2:3],
  LF_res2$relative_expression[2:3],
  LF_res3$relative_expression[2:3],
  LF_res4$relative_expression[2:3],
  LF_res5$relative_expression[2:3],
  LF_res6$relative_expression[2:3],
  LF_res7$relative_expression[2:3],
  LF_res8$relative_expression[2:3]
)

colnames(AllLabeotropheusNorm)<-c("OJ","PJ")

#Remove clearly errounous individuals
AllLabeotropheusNorm<-AllLabeotropheusNorm[c(-1,-4),]

points(AllLabeotropheusNorm, pch=19, col="blue")


###Maylandia###
SubsetData<-FullData[FullData$Species=="Mcallainos" & FullData$Gene %in% c("NOTCH1A", "BACTIN"),]

##Save this next line as a CSV file
SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="NOTCH1A",][,4])
colnames(GeneCq)<-c("BACTIN","NOTCH1A")

#Split by individual
Ind1<-GeneCq[c(1,9,17),]
Ind2<-GeneCq[c(2,10,18),]
Ind3<-GeneCq[c(3,11,19),]
Ind4<-GeneCq[c(4,12,20),]
Ind5<-GeneCq[c(5,13,21),]
Ind6<-GeneCq[c(6,14,22),]
Ind7<-GeneCq[c(7,15,23),]
Ind8<-GeneCq[c(8,16,24),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

MC_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res7 <- pcr_analyze(Ind7, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res8 <- pcr_analyze(Ind8, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllMaylandiaNorm<-rbind.data.frame(
  MC_res1$relative_expression[2:3],
  MC_res2$relative_expression[2:3],
  MC_res3$relative_expression[2:3],
  MC_res4$relative_expression[2:3],
  MC_res5$relative_expression[2:3],
  MC_res6$relative_expression[2:3],
  MC_res7$relative_expression[2:3],
  MC_res8$relative_expression[2:3]
)

colnames(AllMaylandiaNorm)<-c("OJ","PJ")

#Remove clearly errounous individuals
AllMaylandiaNorm<-AllMaylandiaNorm[c(-3,-5,-8),]

points(AllMaylandiaNorm, pch=19, col="purple")


###Collate###
AllSpecies<-rbind.data.frame(
  AllMaylandiaNorm, AllTropeopsNorm, AllLabeotropheusNorm
  
)

TaxaID<-c(rep('Maylandia', 5),
          rep('Tropheops', 6),
          rep('Labeotropheus', 6))

TissueID<-c(rep("OJ",17),
            rep("PJ",17))

AllSpeciesData<-as.numeric(unlist(AllSpecies))

NOTCHFullExpression<-cbind.data.frame(TaxaID,TissueID,AllSpeciesData)



####SMAD7####
###Tropheops###
SubsetData<-FullData[FullData$Species=="Tropheops" & FullData$Gene %in% c("SMAD7", "BACTIN"),]

##Save this next line as a CSV file
SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")
SubsetDataAGsd<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "sd")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="SMAD7",][,4])
colnames(GeneCq)<-c("BACTIN","SMAD7")

#Split by individual
Ind1<-GeneCq[c(1,7,13),]
Ind2<-GeneCq[c(2,8,14),]
Ind3<-GeneCq[c(3,9,15),]
Ind4<-GeneCq[c(4,10,16),]
Ind5<-GeneCq[c(5,11,17),]
Ind6<-GeneCq[c(6,12,18),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

Tr_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
Tr_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllTropeopsNorm<-rbind.data.frame(
  Tr_res1$relative_expression[2:3],
  Tr_res2$relative_expression[2:3],
  Tr_res3$relative_expression[2:3],
  Tr_res4$relative_expression[2:3],
  Tr_res5$relative_expression[2:3],
  Tr_res6$relative_expression[2:3]
)

colnames(AllTropeopsNorm)<-c("OJ","PJ")

plot(AllTropeopsNorm, pch=19, xlim=c(0,3), ylim=c(0,3), col="red")


###Labeotropheus###
SubsetData<-FullData[FullData$Species=="Lfuelleborni" & FullData$Gene %in% c("SMAD7", "BACTIN"),]

##Save this next line as a CSV file
SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="SMAD7",][,4])
colnames(GeneCq)<-c("BACTIN","SMAD7")

#Split by individual
Ind1<-GeneCq[c(1,9,17),]
Ind2<-GeneCq[c(2,10,18),]
Ind3<-GeneCq[c(3,11,19),]
Ind4<-GeneCq[c(4,12,20),]
Ind5<-GeneCq[c(5,13,21),]
Ind6<-GeneCq[c(6,14,22),]
Ind7<-GeneCq[c(7,15,23),]
Ind8<-GeneCq[c(8,16,24),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

LF_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res7 <- pcr_analyze(Ind7, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
LF_res8 <- pcr_analyze(Ind8, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllLabeotropheusNorm<-rbind.data.frame(
  LF_res1$relative_expression[2:3],
  LF_res2$relative_expression[2:3],
  LF_res3$relative_expression[2:3],
  LF_res4$relative_expression[2:3],
  LF_res5$relative_expression[2:3],
  LF_res6$relative_expression[2:3],
  LF_res7$relative_expression[2:3],
  LF_res8$relative_expression[2:3]
)

colnames(AllLabeotropheusNorm)<-c("OJ","PJ")

points(AllLabeotropheusNorm, col="blue", pch=19)


###Maylandia###
SubsetData<-FullData[FullData$Species=="Mcallainos" & FullData$Gene %in% c("SMAD7", "BACTIN"),]

##Save this next line as a CSV file
SubsetDataAG<-aggregate(Cq ~ Individual + Tissue + Gene, data = SubsetData, FUN= "mean")

GeneCq<-cbind.data.frame(SubsetDataAG[SubsetDataAG$Gene=="BACTIN",][,4],SubsetDataAG[SubsetDataAG$Gene=="SMAD7",][,4])
colnames(GeneCq)<-c("BACTIN","SMAD7")

#Split by individual
Ind1<-GeneCq[c(1,9,17),]
Ind2<-GeneCq[c(2,10,18),]
Ind3<-GeneCq[c(3,11,19),]
Ind4<-GeneCq[c(4,12,20),]
Ind5<-GeneCq[c(5,13,21),]
Ind6<-GeneCq[c(6,14,22),]
Ind7<-GeneCq[c(7,15,23),]
Ind8<-GeneCq[c(8,16,24),]

TissueVecB<-rep(c('CF', 'OJ', 'PJ'), each = 1)

MC_res1 <- pcr_analyze(Ind1, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res2 <- pcr_analyze(Ind2, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res3 <- pcr_analyze(Ind3, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res4 <- pcr_analyze(Ind4, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res5 <- pcr_analyze(Ind5, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res6 <- pcr_analyze(Ind6, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res7 <- pcr_analyze(Ind7, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')
MC_res8 <- pcr_analyze(Ind8, group_var = TissueVecB, reference_gene = 'BACTIN', reference_group = 'CF')

AllMaylandiaNorm<-rbind.data.frame(
  MC_res1$relative_expression[2:3],
  MC_res2$relative_expression[2:3],
  MC_res3$relative_expression[2:3],
  MC_res4$relative_expression[2:3],
  MC_res5$relative_expression[2:3],
  MC_res6$relative_expression[2:3],
  MC_res7$relative_expression[2:3],
  MC_res8$relative_expression[2:3]
)

colnames(AllMaylandiaNorm)<-c("OJ","PJ")

points(AllMaylandiaNorm, pch=19, col="purple")


###Collate###
AllSpecies<-rbind.data.frame(
  AllMaylandiaNorm, AllTropeopsNorm, AllLabeotropheusNorm
  
)

TaxaID<-c(rep('Maylandia', 8),
          rep('Tropheops', 6),
          rep('Labeotropheus', 8))

TissueID<-c(rep("OJ",22),
            rep("PJ",22))

AllSpeciesData<-as.numeric(unlist(AllSpecies))

SMADFullExpression<-cbind.data.frame(TaxaID,TissueID,AllSpeciesData)




####Full Collate####
GeneID<-c(rep('DYM', 44),
          rep('SMAD7', 44),
          rep('NOTCH1a', 34))

AllSpeciesExpression<-rbind.data.frame(
  DYMFullExpression, SMADFullExpression, NOTCHFullExpression
  
)

FullSpeciesExpression<-cbind.data.frame(
  GeneID,
  AllSpeciesExpression
)

#Save
FullSpeciesExpression
write.csv(FullSpeciesExpression, "FullSpeciesExpression_boxplot.csv", row.names = F)
