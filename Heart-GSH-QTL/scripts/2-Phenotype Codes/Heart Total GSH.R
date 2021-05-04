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

qtlscan_HeartTotalGSH <- scan1(genoprobs = probs, pheno = pheno["zHeart_TotalGSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_HeartTotalGSH <- scan1perm(genoprobs = probs, pheno = pheno["zHeart_TotalGSH"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Heart GSSG QTL Results.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_HeartTotalGSH = summary(perm_HeartTotalGSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_HeartTotalGSH, map = control$gmap,  main = "Genome Scan for Heart Total GSH", ylim = c(0,11))
  abline(h = threshold_HeartTotalGSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  

  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_HeartTotalGSH, map = control$gmap, threshold = summary(perm_HeartTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksHeartTotalGSH <- find_peaks(scan1_output = qtlscan_HeartTotalGSH, map = control$gmap, threshold = summary(perm_HeartTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksHeartTotalGSH)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_HeartTotalGSH, map = control$pmap, threshold = summary(perm_HeartTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  pmap_peaksHeartTotalGSH <- find_peaks(scan1_output = qtlscan_HeartTotalGSH, map = control$pmap, threshold = summary(perm_HeartTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksHeartTotalGSH)
  
  #save excel sheet with the suggestive (p < 0.2) and significant (p < 0.05) peaks
  write_xlsx(list("Total GSH gmap (cM)" = gmap_peaksHeartTotalGSH,
                  "Total GSH pmap (Mbp)" = pmap_peaksHeartTotalGSH),
             "Total GSH Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

dev.off()

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "Total GSH GWAS - RankZ sexgen.pdf")
out_gwas_HeartTotalGSH <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeart_TotalGSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_HeartTotalGSH$lod, out_gwas_HeartTotalGSH$snpinfo, altcol="green4", gap=0, main = "Heart Total GSH GWAS", ylim = c(0,6))
dev.off()


##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "Total GSH Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_HeartTotalGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_TotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr9, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Heart Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr9, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Gpx1 Position -- Heart Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_HeartTotalGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_TotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_HeartTotalGSH_chr9, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Heart TotalGSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr9, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Gclc Position -- Heart TotalGSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_HeartTotalGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_TotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr3, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Heart TotalGSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr3, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Gclm Position -- Heart TotalGSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_HeartTotalGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_TotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr2, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Heart TotalGSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr2, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Gss Position -- Heart TotalGSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_HeartTotalGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_TotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr8, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Heart TotalGSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_HeartTotalGSH_chr8, map = control$gmap, scan1_output = qtlscan_HeartTotalGSH, main = "Gsr Position -- Heart TotalGSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()



