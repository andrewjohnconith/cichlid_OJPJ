rm(list = ls(all = TRUE)) 
#load qtl library
library(qtl)

#load data
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/QTL/MQM/LOJ/")

F5All <- read.cross(format="csv",file="F5Map_wData_noNAIndsLGcNEWs.csv",
                    na.strings="NA",genotypes=c("AA","AB","BB"),
                    alleles=c("A","B"),convertXdata=TRUE)


gt <- geno.table(F5All)
todrop3 <- rownames(gt[gt$P.value < 1e-150,])
F5Omit <- drop.markers(F5All, todrop3)


augdata1<-fill.geno(F5Omit)


####PC1####
  #If you get the singular matrix, re-run the fill.geno step
scan1<-mqmscan(augdata1,pheno.col=2,plot=T,model="dominance")  #Step 1: scan without cofactors
summary(scan1)

cofactorsindex<-NULL
cofactorsindex<-c(224,119,346)
find.markerindex(augdata1,find.marker(augdata1,chr=7,pos=20))
find.markerindex(augdata1,find.marker(augdata1,chr=4,pos=50))
find.markerindex(augdata1,find.marker(augdata1,chr=11,pos=10))

homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex)) #Step 3: designate cofactors for multivariate scan.

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=2,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

bayesint(homescan1,7,qtl.index=224,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.
bayesint(homescan1,4,qtl.index=119,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.

find.marker(augdata1,chr=7,pos=20)
effectplot(F5All, pheno.col=2, mname1="scaffold_21_2195347",var.flag=c("group"))

find.marker(augdata1,chr=4,pos=50)
effectplot(F5All, pheno.col=2, mname1="scaffold_85_1747922",var.flag=c("group"))

results <- mqmpermutation(augdata1, pheno.col=2, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=1000, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(results)
summary(sig.results)
mqmplot.permutations(results)


####PC2####
scan1<-mqmscan(augdata1,pheno.col=3,plot=T,model="dominance")  #Step 1: scan without cofactors
summary(scan1)

cofactorsindex<-NULL
cofactorsindex<-c(224,84,57,514,392)
find.markerindex(augdata1,find.marker(augdata1,chr=7,pos=20))
find.markerindex(augdata1,find.marker(augdata1,chr=3,pos=55))
find.markerindex(augdata1,find.marker(augdata1,chr=2,pos=35))
find.markerindex(augdata1,find.marker(augdata1,chr=16,pos=5))
find.markerindex(augdata1,find.marker(augdata1,chr=12,pos=15))

homemadecofactors<-mqmsetcofactors(augdata1,cofactors=c(cofactorsindex)) #Step 3: designate cofactors for multivariate scan.

homescan1<-mqmscan(cross=augdata1,cofactors=homemadecofactors,pheno.col=3,cofactor.significance=0.02,verbose=T,plot=T,model="dominance") #using a dominance model tests for both additive effects and dominance, i.e., AA+AB vs BB and BB+AB vs AA, default is additive only
summary(homescan1)

bayesint(homescan1,3,qtl.index=84,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.
bayesint(homescan1,7,qtl.index=224,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.
bayesint(homescan1,12,qtl.index=392,prob=0.95,lodcolumn=1, expandtomarkers=T) #Use this one.

find.marker(augdata1,chr=7,pos=20)
effectplot(F5All, pheno.col=3, mname1="scaffold_21_2195347", var.flag=c("group"))

find.marker(augdata1,chr=3,pos=55)
effectplot(F5All, pheno.col=3, mname1="scaffold_339_46552", var.flag=c("group"))

find.marker(augdata1,chr=12,pos=15)
effectplot(F5All, pheno.col=3, mname1="scaffold_9_4999064", var.flag=c("group"))

results <- mqmpermutation(augdata1, pheno.col=3, scanfunction=mqmscan, cofactors=homemadecofactors, n.perm=1000, multicore=T, plot=F)
sig.results <- mqmprocesspermutation(results)
summary(sig.results)
mqmplot.permutations(results)
