# R01 GSH DO Mapping Code 
# Updated March 2021
# Jess Strosahl


#Heart-GSH-QTL 

#Load in Heart-GSH-QTL.Rdata


#load the command line tools 
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

#set wd

#tells you all of the qtlscans that you have
ls(pattern = "qtl") 


####################################################
## Grab the QTL from all of the qtl scans 
####################################################
  
## I just set threshold to 6 (tells you all of the important qtl peaks with a LOD score > 6)
## map is the qtl2 map you want to use (gmap or pmap)


#tells you all of the qtlscans that you have
  ls(pattern = "qtlscan")

#use cbind to combine all of the qtlscans + take those 2D tables and combining them with another table over and over again
## scans is an R object containing your genome scans from scan1() that are loaded in to the R environment
  scans <- cbind(qtlscan_HeartGSH, qtlscan_HeartGSSG, qtlscan_HeartTotalGSH, qtlscan_HeartGSH_GSSGRatio, qtlscan_HeartRedoxPotentialGSSG2GSH)
  head(scans)
  
  
####################################################
## Review and print the QTL peaks from all of the QTL scans 
####################################################
  
  ## I just set threshold to 6 (tells you all of the important qtl peaks with a LOD score > 6)
  ## map is the qtl2 map you want to use (gmap or pmap)
  qtl_gmap <- find_peaks(scans, map = control$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  qtl_gmap
  
  qtl_pmap <- find_peaks(scans, map = control$pmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  qtl_pmap
  
  #Add marker information
  qtl_gmap$marker.id <- find_marker(map = control$gmap, chr = qtl_gmap$chr, pos = qtl_gmap$pos)
  qtl_gmap$marker.id
  qtl_gmap
  
  qtl_pmap$marker.id <- find_marker(map = control$pmap, chr = qtl_pmap$chr, pos = qtl_pmap$pos)
  qtl_pmap$marker.id
  qtl_pmap
  
  
#set wd
  write_xlsx(list("QTL List RankZ SexGen - cM" = qtl_gmap,
                  "QTL List RankZ SexGen - Mbp" = qtl_pmap),
             "LOD Scores 6.xlsx")
  #gives print out of all LOD peaks > 6
  
  
####################################################
## exporting genes in each interval
####################################################
#set working directory
write_xlsx(list(  "GSH chr16" = HeartGSH_Genes_MGI_chr16, 
                  "GSH chr19" = HeartGSH_Genes_MGI_chr19), 
              "GlutathioneGenesMGI.xlsx")
  
  
  