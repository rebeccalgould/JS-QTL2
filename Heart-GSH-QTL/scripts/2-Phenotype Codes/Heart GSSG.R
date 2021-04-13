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

qtlscan_HeartGSSG <- scan1(genoprobs = probs, pheno = pheno["zHeartGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=10)
perm_HeartGSSG <- scan1perm(genoprobs = probs, pheno = pheno["zHeartGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)

#set working directory
pdf(file = "Heart GSSG QTL Results.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_HeartGSSG = summary(perm_HeartGSSG, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_HeartGSSG, map = control$gmap,  main = "Genome Scan for Heart GSSG", ylim = c(0,11))
  abline(h = threshold_HeartGSSG, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  

  #using gmap (cM)
  find_peaks(scan1_output = qtlscan_HeartGSSG, map = control$gmap, threshold = summary(perm_HeartGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  gmap_peaksHeartGSSG <- find_peaks(scan1_output = qtlscan_HeartGSSG, map = control$gmap, threshold = summary(perm_HeartGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(gmap_peaksHeartGSSG)
  
  #using pmap (Mbp)
  find_peaks(scan1_output = qtlscan_HeartGSSG, map = control$pmap, threshold = summary(perm_HeartGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  pmap_peaksHeartGSSG <- find_peaks(scan1_output = qtlscan_HeartGSSG, map = control$pmap, threshold = summary(perm_HeartGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  print(peaksHeartGSSG)
  
  #save excel sheet with the suggestive (p < 0.2) and significant (p < 0.05) peaks
  write_xlsx(list("GSSG gmap (cM)" = gmap_peaksHeartGSSG,
                  "GSSG pmap (Mbp)" = pmap_peaksHeartGSSG),
             "GSSG Peaks - RankZ sexgen.xlsx")


####################################################
## Estimate QTL Effects (Coefficients) + Connect to SNP and Gene Databases
####################################################

#For Heart GSSG --- Chromosome 16
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #estimate QTL effects by founder strain
  #using gmap (cM)
  chr = 16
  coef_blup_HeartGSSG_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_HeartGSSG_chr16, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,58)
  plot_coefCC(x = coef_blup_HeartGSSG_chr16, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

  #using pmap (Mbp)
  chr = 16
  #could use ci_lo or ci_hi, but in this case, I want a specific chromosome 16 peak
  #start = pmap_peaksHeartGSSG[pmap_peaksHeartGSSG$chr ==  chr,"ci_lo"]
  #end = pmap_peaksHeartGSSG[pmap_peaksHeartGSSG$chr == chr, "ci_hi"] 
  
  pander(pmap_peaksHeartGSSG)
  #based on pmap_peaksHeartGSSG, peak of interest is ~96.99176 Mbp
  #Becca typically does +/- of the QTL interval
  variants_HeartGSSG_chr16 <- query_variants(chr, 95, 99)
  out_snps_HeartGSSG_chr16 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                         chr = chr, start = 95, end = 99, keep_all_snps = TRUE)
  plot_snpasso(out_snps_HeartGSSG_chr16$lod, out_snps_HeartGSSG_chr16$snpinfo, main = "Heart GSSG SNPs")
    
  HeartGSSG_Genes_MGI_chr16 <- query_genes_mgi(chr = chr, start = 95, end = 99)
  plot(out_snps_HeartGSSG_chr16$lod, out_snps_HeartGSSG_chr16$snpinfo, drop_hilit=1.5, genes = HeartGSSG_Genes_MGI_chr16, main = "Heart GSSG Genes MGI")
  
  
#For Heart GSSG --- Chromosome 19
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #estimate QTL effects by founder strain
  #using gmap (cM)
  chr = 19
  coef_blup_HeartGSSG_chr19 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_HeartGSSG_chr19, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,56)
  plot_coefCC(x = coef_blup_HeartGSSG_chr19, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 19
  #could use ci_lo or ci_hi, but in this case, I want a specific chromosome 16 peak
  #start = pmap_peaksHeartGSSG[pmap_peaksHeartGSSG$chr ==  chr,"ci_lo"]
  #end = pmap_peaksHeartGSSG[pmap_peaksHeartGSSG$chr == chr, "ci_hi"] 
  
  pander(pmap_peaksHeartGSSG)
  #based on pmap_peaksHeartGSSG, peak of interest is ~57.15 Mbp
  #Becca typically does +/- of the QTL interval
  variants_HeartGSSG_chr19 <- query_variants(chr, 56, 58.5)
  out_snps_HeartGSSG_chr19 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = 56, end = 58.5, keep_all_snps = TRUE)
  plot_snpasso(out_snps_HeartGSSG_chr19$lod, out_snps_HeartGSSG_chr19$snpinfo, main = "Heart GSSG SNPs")
  
  HeartGSSG_Genes_MGI_chr19 <- query_genes_mgi(chr = chr, start = 56, end = 58.5)
  plot(out_snps_HeartGSSG_chr19$lod, out_snps_HeartGSSG_chr19$snpinfo, drop_hilit=1.5, genes = HeartGSSG_Genes_MGI_chr19, main = "Heart GSSG Genes MGI")

  
dev.off()

####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

pdf(file = "GSSG GWAS - RankZ sexgen.pdf")
out_gwas_HeartGSSG <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeartGSSG"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_HeartGSSG$lod, out_gwas_HeartGSSG$snpinfo, altcol="green4", gap=0, main = "Heart GSSG GWAS", ylim = c(0,6))
dev.off()


##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSSG Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_HeartGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSSG_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_HeartGSSG_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Gpx1 Position -- Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_HeartGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_HeartGSSG_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_HeartGSSG_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Gclc Position -- Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_HeartGSSG_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSSG_chr3, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_HeartGSSG_chr3, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Gclm Position -- Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_HeartGSSG_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSSG_chr2, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_HeartGSSG_chr2, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Gss Position -- Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_HeartGSSG_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSSG_chr8, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(19,22.5)
plot_coefCC(x = coef_blup_HeartGSSG_chr8, map = control$gmap, scan1_output = qtlscan_HeartGSSG, main = "Gsr Position -- Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()



