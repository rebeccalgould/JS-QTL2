# R01 GSH DO Mapping Code 
# Updated March 2021
# Jess Strosahl


#Heart-GSH-QTL 

#Load in Heart-GSH-QTL.Rdata


#Creates Heat Maps that shows the individual CC founder allele contribution at a specific SNP from the GigaMUGA

#Created by Becca Gould with Sue McClatchy on 9/24/2020
#Based on Greg Keele's heat map code

#load the command line tools 
library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(RSQLite)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (RSQLite)
library (writexl)
library (pander)


####################################################
## Create the probability plotting function (made by Greg Keele)
####################################################

prob_plot <- function(pheno_vec,
                      pheno_name = NULL,
                      genoprobs,
                      qtl_chr,
                      qtl_marker,
                      cols = gray(10000:1/10000),
                      label_col = as.character(qtl2::CCcolors),
                      founders = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"),
                      main = "") {
  sorted_pheno <- sort(pheno_vec)
  image(genoprobs[[qtl_chr]][names(sorted_pheno), rev(LETTERS[1:8]), qtl_marker] * 2,
        yaxt = "n", xaxt = "n", col = cols)
  axis(2, at = seq(0, 8, 1 + 1/8)/8, labels = FALSE,
       lty = 0, srt = 90, las = 2)
  mtext(text = main, side = 3, padj = -1, cex = 1.25)
  mtext(text = rev(founders), side = 2, col = rev(label_col), at = seq(0, 1, length.out = 8),
        las = 1, cex = 1.25, adj = 1.25, font = 2)
  mtext(text = paste("lower", "<--", ifelse(is.null(pheno_name), "phenotype", pheno_name), "-->", "higher"), side = 1, padj = 1.25, cex = 1.25)
}

    
####################################################
## Alter the pheno file accordingly 
## pheno file needs to be a data frame for QTL analysis, but for these allele probability plots, it must be a matrix
####################################################

  #need to make the pheno file a matrix so that it runs in the code (currently a data frame)
  #first need to identify what specifically to make part of the phenotype matrix (only need transformed data!)
    names(pheno)
    #from this, I've identified I only need columns 24-34
    #pheno_mat is the matrix of outcomes (phenotypes)
    pheno_mat <- as.matrix(pheno[c(24:34)])
    
  #check rownames to make sure they are already set as the write row names (they are)
    rownames(pheno[c(24:34)])


####################################################
## ALLELE PLOTS CODE - LOOP
####################################################
  
#set working directory to store the plots
  pdf(file = "allele-plots_cM - RankZ sexgen.pdf") # create a file called allele-plots.pdf
  # loop through all qtl_gmap above lod threshold of 6 and create an individual plot
  for (i in 1:dim(qtl_gmap)[1]) {
    prob_plot(pheno_vec = pheno_mat[,qtl_gmap$lodcolumn[i]],
              genoprobs = probs,
              qtl_chr = qtl_gmap$chr[i],
              qtl_marker = qtl_gmap$marker.id[i],
              main = paste("lodindex", qtl_gmap$lodindex[i], "Chr", qtl_gmap$chr[i], qtl_gmap$marker.id[i], qtl_gmap$pos[i], qtl_gmap$lodcolumn[i]))
  }
  # be sure to turn the graphics output off at the end!
  dev.off()  
  
#set working directory to store the plots
  pdf(file = "allele-plots_Mbp - RankZ sexgen.pdf") # create a file called allele-plots.pdf
  # loop through all qtl_gmap above lod threshold of 6 and create an individual plot
  for (i in 1:dim(qtl_pmap)[1]) {
    prob_plot(pheno_vec = pheno_mat[,qtl_pmap$lodcolumn[i]],
              genoprobs = probs,
              qtl_chr = qtl_pmap$chr[i],
              qtl_marker = qtl_pmap$marker.id[i],
              main = paste(qtl_pmap$lodindex[i], "Chr", qtl_pmap$chr[i], qtl_pmap$marker.id[i], qtl_pmap$pos[i], qtl_pmap$lodcolumn[i]))
  }
  # be sure to turn the graphics output off at the end!
  dev.off()  
  

  
####################################################
## Printing out the results
####################################################
  
#INDIVIDUAL PLOTS CODE:
  prob_plot(pheno_vec = pheno_mat[,qtl$lodcolumn[1]],
            genoprobs = probs,
            qtl_chr = qtl$chr[1],
            qtl_marker = qtl$marker.id[1])

  #shows you the individual column for each lodindex
  qtl$lodcolumn[1]
  #zLiverGSH
  #[1] is a vector of values for lod index -- in this case, [1] refers to the first peak of zLiverGSH


  


  






