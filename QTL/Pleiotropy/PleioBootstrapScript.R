#install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl"))

#library(devtools)
#install_github("rqtl/qtl2")

#devtools::install_github("fboehm/qtl2pleio")

#Import packages
library(qtl); library(qtl2); library(qtl2pleio); library(ggplot2); library(stringr); library(dplyr); library(readr)

#Read in map data
setwd("/Users/home/Dropbox/Amherst/Post-Doc/Muscles/OJ_PJ/Writing/MS/Writing/ReSubmission/FinalSubmission/GitHub/QTL/Pleiotropy/")
F5 <- read.cross(format="csv",file="F5Map_wData_noNAIndsLGcPleio.csv",
                    na.strings="NA",genotypes=c("AA","AB","BB"),
                    alleles=c("A","B"),convertXdata=TRUE)

#Convert map to r/qtl2 format
F5All <- convert2cross2(F5)

#Data engineering
pmap <- insert_pseudomarkers(map=F5All$gmap, step=1)
pr <- calc_genoprob(cross=F5All, map=pmap, error_prob=0.002)
pp <- genoprob_to_alleleprob(pr)
kinship__loco <- calc_kinship(pp, "loco")


#Pleiotropy scan# #LOJ_PC1 - LPJ_PC1
LOJPC1_LPJPC1 <- scan_pvl(probs = pp$`7`,
                pheno = F5All$pheno[,c(2,4)],
                kinship = kinship__loco$`7`,
                n_snp = 145)
save(LOJPC1_LPJPC1, file = "Ple_LOJPC1-LPJPC1.rda")

#Pleiotropy scan# #LOJ_PC2 - LPJ_PC1
LOJPC2_LPJPC1 <- scan_pvl(probs = pp$`7`,
                pheno = F5All$pheno[,c(3,4)],
                kinship = kinship__loco$`7`,
                n_snp = 145)
save(LOJPC2_LPJPC1, file = "Ple_LOJPC2-LPJPC1.rda")


##Plot two-dimensional scan
	#Calculate profile LOD for each trait
		#LOJPC1 - LPJPC1
out_lodsOJ1PJ1 <- LOJPC1_LPJPC1 %>%
  calc_profile_lods() %>%
  add_pmap(pmap = F5All$gmap$`7`)

LOJPC1_LPJPC1_Plot<-out_lodsOJ1PJ1 %>% 
  ggplot() + 
  geom_line(aes(x = marker_position, y = profile_lod, colour = trait), size=1.3, alpha = 0.90) +
  labs(x = "Marker position on LG7 (cM)", y = "LOD Score") +
  scale_colour_manual(values=c(pleiotropy="#404040",tr1="#ca0020",tr2="#0571b0"), name="Trace", labels=c("Pleiotropy", "L. Oral Jaw PC1", "L. Pharyngeal Jaw PC1")) +
  theme(legend.key=element_blank(), legend.position="bottom", axis.ticks = element_blank(),
        panel.border = element_rect(colour = "dark gray", fill=NA, size=1),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "dark gray"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "dark gray"),
        panel.background = element_rect(fill = "white", colour = "white", size = 1, linetype = "solid"))

MyRatio <- with(out_lodsOJ1PJ1, diff(range(marker_position))/diff(range(profile_lod)))
LOJPC1_LPJPC1_PlotL <- LOJPC1_LPJPC1_Plot + coord_fixed(ratio=MyRatio)

(mylrt_OJ1PJ1 <- calc_lrt_tib(LOJPC1_LPJPC1))
find_pleio_peak_tib(LOJPC1_LPJPC1, start_snp = 1)

##Plot two-dimensional scan
	#Calculate profile LOD for each trait
		#LOJPC2 - LPJPC1
out_lodsOJ2PJ1 <- LOJPC2_LPJPC1 %>%
  calc_profile_lods() %>%
  add_pmap(pmap = F5All$gmap$`7`)

LOJPC2_LPJPC1_Plot<-out_lodsOJ2PJ1 %>% 
  ggplot() + 
  geom_line(aes(x = marker_position, y = profile_lod, colour = trait), size=1.3, alpha = 0.90) +
  labs(x = "Marker position on LG7 (cM)", y = "LOD Score") +
  scale_colour_manual(values=c(pleiotropy="#404040",tr1="#f4a582",tr2="#0571b0"), name="Trace", labels=c("Pleiotropy", "L. Oral Jaw PC2", "L. Pharyngeal Jaw PC1")) +
  theme(legend.key=element_blank(), legend.position="bottom", axis.ticks = element_blank(),
        panel.border = element_rect(colour = "dark gray", fill=NA, size=1),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "dark gray"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "dark gray"),
        panel.background = element_rect(fill = "white", colour = "white", size = 1, linetype = "solid"))

MyRatio <- with(out_lodsOJ2PJ1, diff(range(marker_position))/diff(range(profile_lod)))
LOJPC2_LPJPC1_PlotR <- LOJPC2_LPJPC1_Plot + coord_fixed(ratio=MyRatio)

(mylrt_OJ2PJ1 <- calc_lrt_tib(LOJPC2_LPJPC1))
find_pleio_peak_tib(LOJPC2_LPJPC1, start_snp = 1)


library("cowplot")
  plot_grid(LOJPC1_LPJPC1_PlotL, LOJPC1_LPJPC1_LG4PlotR)

#LOJ_PC1 - LPJ_PC1
#(mylrt_OJ1PJ1 <- calc_lrt_tib(LOJPC1_LPJPC1))
#0.02218583
#find_pleio_peak_tib(LOJPC1_LPJPC1, start_snp = 1)
#log10lik44
#44

#LOJ_PC2 - LPJ_PC1
#(mylrt_OJ2PJ1 <- calc_lrt_tib(LOJPC2_LPJPC1))
#0.1992048
#find_pleio_peak_tib(LOJPC2_LPJPC1, start_snp = 1)
#log10lik46
#46



####LOJ_PC1 - LPJ_PC1 LG7 Bootstrap####
set.seed(20201029)
##
pleio_peak_index <- 44 #LOJ_PC1 - LPJ_PC1
##

###############

# insert pseudomarkers
pmap <- insert_pseudomarkers(map=F5All$gmap, step=1)
pr <- calc_genoprob(cross=F5All, map=pmap, error_prob=0.002)
pp <- genoprob_to_alleleprob(pr)
kinship__loco <- calc_kinship(pp, "loco")

## ------------------------------------------------------------------------
phe <- F5All$pheno[,c(2,4)] #LOJ-PC1, LPJ-PC1
k <- kinship__loco$`7`
## ------------------------------------------------------------------------

###############
# simulate a phenotype
X1 <- pp$`7`[ , , pleio_peak_index]
## remove subjects with missing values of phenotype
is.na(phe[ , 1]) | is.na(phe[ , 2]) -> missing_indic
phe_nona <- phe[!missing_indic, ]
Xpre_nona <- X1[!missing_indic, ]
k_nona <- k[!missing_indic, !missing_indic]

##
gemma2::stagger_mats(Xpre_nona,Xpre_nona) -> X

calc_covs(pheno = phe_nona, kinship = k_nona) -> cc_out
(cc_out$Vg -> Vg)
(cc_out$Ve -> Ve)
# calculate Sigma
calc_Sigma(Vg = Vg, Ve = Ve,  kinship =  k_nona) -> Sigma
solve(Sigma) -> Sigma_inv
# calc Bhat 
B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe_nona)
# Start loop
lrtOJpc1PJpc1 <- numeric()
for (i in 1:100){
  sim1(X = X, B = B, Sigma = Sigma) -> foo
  matrix(foo, ncol = 2, byrow = FALSE) -> Ysim
  rownames(Ysim) <- rownames(phe_nona)
  colnames(Ysim) <- c("L. Oral Jaw PC1", "L. Pharyngeal Jaw PC1")
  scan_pvl(probs = pp$`7`[!missing_indic, , ], pheno = Ysim, kinship = k_nona, start_snp = 1, n_snp = 145) -> loglik
  # in above call, s1 & nsnp come from command line args
  calc_lrt_tib(loglik) -> lrtOJpc1PJpc1[i]
}

#(pvalue <- mean(lrt >= mylrt))

run_num<-"OJpc1PJpc1"
proc_num<-20201029

fn_out <- paste0("recla-boot_", run_num, "_", proc_num, ".csv")
write.csv(lrtOJpc1PJpc1, fn_out)

save(lrtOJpc1PJpc1, save="OJpc1PJpc1_20201029.rda")



####LOJ_PC2 - LPJ_PC1 LG7 Bootstrap####
set.seed(20201029)
##
pleio_peak_index <- 46 #LOJ_PC2 - LPJ_PC1
##

###############

# insert pseudomarkers
pmap <- insert_pseudomarkers(map=F5All$gmap, step=1)
pr <- calc_genoprob(cross=F5All, map=pmap, error_prob=0.002)
pp <- genoprob_to_alleleprob(pr)
kinship__loco <- calc_kinship(pp, "loco")

## ------------------------------------------------------------------------
phe <- F5All$pheno[,c(3,4)] #LOJ-PC2, LPJ-PC1
k <- kinship__loco$`7`
## ------------------------------------------------------------------------

###############
# simulate a phenotype
X1 <- pp$`7`[ , , pleio_peak_index]
## remove subjects with missing values of phenotype
is.na(phe[ , 1]) | is.na(phe[ , 2]) -> missing_indic
phe_nona <- phe[!missing_indic, ]
Xpre_nona <- X1[!missing_indic, ]
k_nona <- k[!missing_indic, !missing_indic]

##
gemma2::stagger_mats(Xpre_nona,Xpre_nona) -> X

calc_covs(pheno = phe_nona, kinship = k_nona) -> cc_out
(cc_out$Vg -> Vg)
(cc_out$Ve -> Ve)
# calculate Sigma
calc_Sigma(Vg = Vg, Ve = Ve,  kinship =  k_nona) -> Sigma
solve(Sigma) -> Sigma_inv
# calc Bhat 
B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe_nona)
# Start loop
lrtOJpc2PJpc1 <- numeric()
for (i in 1:100){
  sim1(X = X, B = B, Sigma = Sigma) -> foo
  matrix(foo, ncol = 2, byrow = FALSE) -> Ysim
  rownames(Ysim) <- rownames(phe_nona)
  colnames(Ysim) <- c("L. Oral Jaw PC2", "L. Pharyngeal Jaw PC1")
  scan_pvl(probs = pp$`7`[!missing_indic, , ], pheno = Ysim, kinship = k_nona, start_snp = 1, n_snp = 145) -> loglik
  # in above call, s1 & nsnp come from command line args
  calc_lrt_tib(loglik) -> lrtOJpc2PJpc1[i]
}

#(pvalue <- mean(lrt >= mylrt))

run_num<-"OJpc2PJpc1"
proc_num<-20201029

fn_out <- paste0("recla-boot_", run_num, "_", proc_num, ".csv")
write.csv(lrtOJpc2PJpc1, fn_out)

save(lrtOJpc2PJpc1, save="OJpc2PJpc1_20201029.rda")





####LG4####
#Pleiotropy scan# #LOJ_PC1 - LPJ_PC1
LOJPC1_LPJPC1_LG4 <- scan_pvl(probs = pp$`4`,
                          pheno = F5All$pheno[,c(2,4)],
                          kinship = kinship__loco$`4`,
                          n_snp = 100)
save(LOJPC1_LPJPC1_LG4, file = "Ple_LOJPC1-LPJPC1_LG4.rda")

out_lodsOJ1PJ1_LG4 <- LOJPC1_LPJPC1_LG4 %>%
  calc_profile_lods() %>%
  add_pmap(pmap = F5All$gmap$`4`)

LOJPC1_LPJPC1_LG4Plot<-out_lodsOJ1PJ1_LG4 %>% 
  ggplot() + 
  geom_line(aes(x = marker_position, y = profile_lod, colour = trait), size=1.3, alpha = 0.90) +
  labs(x = "Marker position on LG4 (cM)", y = "LOD Score") +
  scale_colour_manual(values=c(pleiotropy="#404040",tr1="#ca0020",tr2="#0571b0"), name="Trace", labels=c("Pleiotropy", "L. Oral Jaw PC1", "L. Pharyngeal Jaw PC1")) +
  theme(legend.key=element_blank(), legend.position="bottom", axis.ticks = element_blank(),
        panel.border = element_rect(colour = "dark gray", fill=NA, size=1),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "dark gray"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "dark gray"),
        panel.background = element_rect(fill = "white", colour = "white", size = 1, linetype = "solid"))

MyRatio <- with(out_lodsOJ1PJ1_LG4, diff(range(marker_position))/diff(range(profile_lod)))
LOJPC1_LPJPC1_LG4PlotR<-LOJPC1_LPJPC1_LG4Plot + coord_fixed(ratio=MyRatio)

(mylrt_OJ1PJ1_LG4 <- calc_lrt_tib(LOJPC1_LPJPC1_LG4))
find_pleio_peak_tib(LOJPC1_LPJPC1_LG4, start_snp = 1)

#(mylrt_OJ1PJ1_LG4 <- calc_lrt_tib(LOJPC1_LPJPC1_LG4))
#1.85236
#find_pleio_peak_tib(LOJPC1_LPJPC1_LG4, start_snp = 1)
#log10lik59 
#59 



####LOJ_PC1 - LPJ_PC1 LG4 Bootstrap####
set.seed(20201029)
##
pleio_peak_index <- 59 #LOJ_PC1 - LPJ_PC1
##

###############

# insert pseudomarkers
pmap <- insert_pseudomarkers(map=F5All$gmap, step=1)
pr <- calc_genoprob(cross=F5All, map=pmap, error_prob=0.002)
pp <- genoprob_to_alleleprob(pr)
kinship__loco <- calc_kinship(pp, "loco")

## ------------------------------------------------------------------------
phe <- F5All$pheno[,c(2,4)] #LOJ-PC1, LPJ-PC1
k <- kinship__loco$`4`
## ------------------------------------------------------------------------

###############
# simulate a phenotype
X1 <- pp$`4`[ , , pleio_peak_index]
## remove subjects with missing values of phenotype
is.na(phe[ , 1]) | is.na(phe[ , 2]) -> missing_indic
phe_nona <- phe[!missing_indic, ]
Xpre_nona <- X1[!missing_indic, ]
k_nona <- k[!missing_indic, !missing_indic]

##
gemma2::stagger_mats(Xpre_nona,Xpre_nona) -> X

calc_covs(pheno = phe_nona, kinship = k_nona) -> cc_out
(cc_out$Vg -> Vg)
(cc_out$Ve -> Ve)
# calculate Sigma
calc_Sigma(Vg = Vg, Ve = Ve,  kinship =  k_nona) -> Sigma
solve(Sigma) -> Sigma_inv
# calc Bhat 
B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe_nona)
# Start loop
lrtOJpc1PJpc1_LG4 <- numeric()
for (i in 1:100){
  sim1(X = X, B = B, Sigma = Sigma) -> foo
  matrix(foo, ncol = 2, byrow = FALSE) -> Ysim
  rownames(Ysim) <- rownames(phe_nona)
  colnames(Ysim) <- c("L. Oral Jaw PC1", "L. Pharyngeal Jaw PC1")
  scan_pvl(probs = pp$`4`[!missing_indic, , ], pheno = Ysim, kinship = k_nona, start_snp = 1, n_snp = 100) -> loglik
  # in above call, s1 & nsnp come from command line args
  calc_lrt_tib(loglik) -> lrtOJpc1PJpc1_LG4[i]
}

#(pvalue <- mean(lrt >= mylrt))

run_num<-"OJpc1PJpc1_LG4"
proc_num<-20201029

fn_out <- paste0("recla-boot_", run_num, "_", proc_num, ".csv")
write.csv(lrtOJpc1PJpc1_LG4, fn_out)

save(lrtOJpc1PJpc1_LG4, save="lrtOJpc1PJpc1LG4_20201029.rda")
