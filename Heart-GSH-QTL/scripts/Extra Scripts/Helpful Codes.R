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