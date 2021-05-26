#Helpful codes

#comparing two QTL plots - in this case, comparing the QTL scan for Heart GSH versus Kidney Redox Potential
plot_scan1(x = qtlscan_HeartGSH, map = control$gmap,  ylim = c(0,11))
plot_scan1(x = qtlscan_KidneyEh, map = control$gmap, col = "#009999", add = TRUE)
legend("topleft", lwd=2, col=c("darkslateblue", "#009999"), c("Heart GSH", "Kidney Redox Potential"), bg="gray90")



####################################################
## Heritability calculation - the ratio of genetic variance to total variance using a linear mixed model
####################################################

#example of heritability calculation for HeartGSH
herit_HeartGSH_sexgen <- est_herit(pheno["zHeartGSH"], kinship_lmm, sexgen, cores = 10)



####################################################
## How to get X chromosome-specific permutation tests and plot them
####################################################

Xcovar = get_x_covar(control)
perm_strata = mat2strata(Xcovar)
perm_X_KidneyGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSH"], addcovar = sexgen, n_perm = 1000, perm_Xsp = TRUE, perm_strata = perm_strata, chr_lengths = chr_lengths(control$gmap), cores=10)

summary(perm_X_KidneyGSH, alpha=c(0.2, 0.1, 0.05))
#Autosome LOD thresholds (1000 permutations)
#zKidneyGSH
#0.2        6.97
#0.1        7.33
#0.05       7.80

#X chromosome LOD thresholds (18090 permutations)
#zKidneyGSH
#0.2        6.50
#0.1        6.89
#0.05       7.40

#set working directory
pdf(file = "GSH QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#the normal process, defining the thresholds 
par(mar=c(4.1, 4.1, 2.6, 2.6))
threshold_KidneyGSH = summary(perm_KidneyGSH, alpha = c(0.2, 0.1, 0.05))
threshold_X_KidneyGSH = summary(perm_X_KidneyGSH, alpha = c(0.2, 0.1, 0.05))

#normal plot, showing autosome specific thresholds ONLY
plot_scan1(x = qtlscan_KidneyGSH, map = control$gmap,  main = "Genome Scan for Kidney GSH", ylim = c(0,11))
abline(h = threshold_KidneyGSH, col = c("purple", "red", "blue"), lwd = 2)
plot_scan1(x = qtlscan_KidneyGSH, map = control$gmap,  main = "Genome Scan for Kidney GSH", ylim = c(0,11))
abline(h = threshold_KidneyGSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")

#showing X chromosome significance thresholds ONLY
plot_scan1(x = qtlscan_KidneyGSH, map = control$gmap,  main = "Genome Scan for Kidney GSH (X Chrom)", ylim = c(0,11))
abline(h = c(6.50, 6.89, 7.40), col = c("purple", "red", "blue"), lwd = 2)
plot_scan1(x = qtlscan_KidneyGSH, map = control$gmap,  main = "Genome Scan for Kidney GSH (X Chrom)", ylim = c(0,11))
abline(h = c(6.50, 6.89, 7.40), col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
#could also do:
  #abline(h = threshold_X_KidneyGSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")

#plotting separate autosome versus X axis significance thresholds - WHAT YOU WANT, DASHED VERSION
plot_scan1(x = qtlscan_KidneyGSH, map = control$gmap,  main = "Genome Scan for Kidney GSH (Autosome vs X)", ylim = c(0,11))
segments(x0 = 0, y0 = threshold_X_KidneyGSH$A, x1 = 1695, y1 =   threshold_X_KidneyGSH$A, col = c("purple", "red", "blue"), lwd = 2)
segments(x0 = 1695, y0 = threshold_X_KidneyGSH$X, x1 = 2000, y1 = threshold_X_KidneyGSH$X, col = c("purple", "red", "blue"), lwd = 2)

#plotting separate autosome versus X axis significance thresholds - WHAT YOU WANT
plot_scan1(x = qtlscan_KidneyGSH, map = control$gmap,  main = "Genome Scan for Kidney GSH (Autosome vs X)", ylim = c(0,11))
segments(x0 = 0, y0 = threshold_X_KidneyGSH$A, x1 = 1695, y1 =   threshold_X_KidneyGSH$A, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
segments(x0 = 1695, y0 = threshold_X_KidneyGSH$X, x1 = 2000, y1 = threshold_X_KidneyGSH$X, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
