# R01 GSH DO Mapping Code 
# Updated March 2021
# Jess Strosahl

#Heart-GSH-QTL 

#Load in Heart-GSH-QTL.Rdata
#Run RankZ Data Prep R Script before doing this**


#set wd

library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (pander)
library (writexl)
library (RSQLite)

####################################################
## Plot Genome Scans with Permutation Tests
####################################################

qtlscan_HeartGSHGSSGRatio <- scan1(genoprobs = probs, pheno = pheno["zHeart_GSHGSSGRatio"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_HeartGSHGSSGRatio <- scan1perm(genoprobs = probs, pheno = pheno["zHeart_GSHGSSGRatio"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Heart GSH/GSSG Ratio QTL Results.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_HeartGSHGSSGRatio = summary(perm_HeartGSHGSSGRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_HeartGSHGSSGRatio, map = control$gmap,  main = "Genome Scan for Heart GSH/GSSG Ratio", ylim = c(0,11))
  abline(h = threshold_HeartGSHGSSGRatio, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  

  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_HeartGSHGSSGRatio, map = control$gmap, threshold = summary(perm_HeartGSHGSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksHeartGSHGSSGRatio <- find_peaks(scan1_output = qtlscan_HeartGSHGSSGRatio, map = control$gmap, threshold = summary(perm_HeartGSHGSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksHeartGSHGSSGRatio)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_HeartGSHGSSGRatio, map = control$pmap, threshold = summary(perm_HeartGSHGSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  pmap_peaksHeartGSHGSSGRatio <- find_peaks(scan1_output = qtlscan_HeartGSHGSSGRatio, map = control$pmap, threshold = summary(perm_HeartGSHGSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksHeartGSHGSSGRatio)
  
  #save excel sheet with the suggestive (p < 0.2) and significant (p < 0.05) peaks
  write_xlsx(list("GSH/GSSG Ratio gmap (cM)" = gmap_peaksHeartGSHGSSGRatio,
                  "GSH/GSSG Ratio pmap (Mbp)" = pmap_peaksHeartGSHGSSGRatio),
             "GSH/GSSG Ratio Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

dev.off()

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "GSH/GSSG Ratio GWAS - RankZ sexgen.pdf")
out_gwas_HeartGSHGSSGRatio <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeart_GSHGSSGRatio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_HeartGSHGSSGRatio$lod, out_gwas_HeartGSHGSSGRatio$snpinfo, altcol="green4", gap=0, main = "Heart GSH/GSSG Ratio GWAS", ylim = c(0,6))
dev.off()


##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSH/GSSG Ratio Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_HeartGSHGSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_GSHGSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Gpx1 Position -- Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_HeartGSHGSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_GSHGSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Heart GSHGSSGRatio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Gclc Position -- Heart GSHGSSGRatio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_HeartGSHGSSGRatio_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_GSHGSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr3, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Heart GSHGSSGRatio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr3, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Gclm Position -- Heart GSH GSSG Ratio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_HeartGSHGSSGRatio_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_GSHGSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr2, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Heart GSHGSSGRatio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr2, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Gss Position -- Heart GSHGSSGRatio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_HeartGSHGSSGRatio_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_GSHGSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr8, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_HeartGSHGSSGRatio_chr8, map = control$gmap, scan1_output = qtlscan_HeartGSHGSSGRatio, main = "Gsr Position -- Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()



