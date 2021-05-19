# R01 GSH DO Mapping Code - JS QTL2
# Updated March 2021
# Jess Strosahl

#Heart GSH QTL
#RankZ TRANSFORMATION AND DATA PREP

#load the command line tools - see https://github.com/Rdatatable/data.table/wiki/Installation for more information - must do every time you open up the Rproject!
library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (writexl)
library (RSQLite)

#setwd

####################################################
## Helpful Sites
####################################################

#Karl Broman info on R/qtl2: https://kbroman.org/qtl2/assets/vignettes/user_guide.html#QTL_analysis_in_Diversity_Outbred_mice
#DO mapping course created by Sue McClatchy: https://smcclatchy.github.io/mapping/ 
#link for R colors: http://www.sthda.com/english/wiki/colors-in-r


####################################################
## Read in the control file (gm.json)
####################################################

#Load in the control file to tell R/qtl2 to read the data and title it according to your project. For mine, it's R01_GSH_DO_QTLdata. 
#For this to work, all of the files in the control file need to be in the folder. Ex: if calling for the "genoprobs" file, it needs to actually be in the folder.
  control <- read_cross2(file = "~/JS-QTL2/data/R01_GSH_DO_control.json")

####################################################
## Genotype probabilities and allele probabilities - provided by Belinda and Vivek
####################################################

#read in the genoprobs file that is sorted by chromosomes in numerical order - the 8state.rds is the allele probabilities, the 36state.rds is the genotype probabilities
#^this is actually the ALlELE probabilities, but for simplicity, we will call it "probs"
  probs <- readRDS("~/JS-QTL2/data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")

  nrow(data.frame(R01_GSH_DO_QTLdata$gmap[1]))
  #should be 10415
  dim(probs[[1]])
  #should be 347 individuals, 8 alleles, 10415 markers


####################################################
## Variant files
####################################################

#Will need these for the final lesson episodes on SNP association mapping and QTL analysis in Diversity Outbred mice. Make sure they are the most updated versions!
  query_variants <- create_variant_query_func("~/JS-QTL2/data/cc_variants.sqlite")
  query_genes_mgi <- create_gene_query_func("~/JS-QTL2/data/mouse_genes_mgi.sqlite")
  query_genes <- create_gene_query_func("~/JS-QTL2/data/mouse_genes.sqlite")

####################################################
## Calculating kinship
####################################################

#calculate the kinship loco
#you can increase the cores amount if you have more cores in your computer. For mine, I have 18 cores available so to speed it up, I'll use 10 of them.
  kinship_loco <- calc_kinship(probs = probs, "loco", use_allele_probs = TRUE, cores = 10)

#Create the r plot of the kinship matrix
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  image(1:nrow(kinship_loco[[1]]), 1:ncol(kinship_loco[[1]]), kinship_loco[[1]][,ncol(kinship_loco[[1]]):1], xlab = "Samples", 
        ylab = "Samples", yaxt = "n", main = "Kinship between samples", 
        breaks = 0:100/100, col = heat.colors(length(0:100) - 1))


####################################################
## Editing the phenotype file to make it R/qtl2-friendly
####################################################
#to have the phenotype file for reference - can be used when plotting the data to see if it needs to be transformed
  pheno <- read.csv(file = "~/JS-QTL2/data/R01_GSH_DO_pheno_covar.csv", header = TRUE)

#make row names the ID of each sample
  rownames(pheno) <- pheno$id
#checking pheno file
  pheno[1:10,]

#change sexes to numeric variables
  pheno$sex[pheno$sex == "M"] <- 1
  pheno$sex[pheno$sex == "F"] <- 0

#both added covariates must be numeric, not characters 
  pheno$sex <- as.numeric(pheno$sex)
  pheno$generation <- as.numeric(pheno$generation)  

#check pheno file
  pheno[1:10,]
  str(pheno)


####################################################
## checking if data should be transformed
####################################################

#gives you the names of each phenotype
  control[["pheno"]]

#The rankZ function is a nonparametric method that replaces the values of a variable with their rank in 
#ascending order (e.g. values 0.2, 1.3, and 4.3 are replaced by 1, 2 and 3). This effectively forces the data into a normal distribution and eliminates outliers.
  rankZ <- function(x) {x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))}
  
#Rank Z transformations of each phenotype
  pheno$zHeartGSH = rankZ(pheno$Heart_GSH)
  pheno$zHeartGSSG = rankZ(pheno$Heart_GSSG)
  pheno$zHeartTotalGSH = rankZ(pheno$Heart_TotalGSH)
  pheno$zHeartGSHGSSGRatio = rankZ(pheno$Heart_GSHGSSGRatio)
  pheno$zHeartEh = rankZ(pheno$Heart_Eh)
  

#####Plot the transformations  

#set working directory 
pdf(file = "Box Plots and QQ Norm Plots - RankZ and Raw.pdf")
  ##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##
  
#For Heart GSH
  boxplot(pheno$Heart_GSH, main = "Heart GSH Box Plot")
  boxplot(pheno$Heart_GSH~pheno$generation, main = "Heart GSH Box Plot - by generation")
  boxplot(pheno$zHeartGSH, main = "Rank Z Heart GSH Box Plot")
  boxplot(pheno$zHeartGSH~pheno$generation, main = "Rank Z Heart GSH Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zHeartGSH, main = "Normal QQ Plot - Rank Z Heart GSH")

#For Heart GSSG
  boxplot(pheno$Heart_GSSG, main = "Heart GSSG Box Plot")
  boxplot(pheno$Heart_GSSG~pheno$generation, main = "Heart GSSG Box Plot - by generation")
  boxplot(pheno$zHeartGSSG, main = "Rank Z Heart GSSG Box Plot")
  boxplot(pheno$zHeartGSSG~pheno$generation, main = "Rank Z Heart GSSG Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zHeartGSSG, main = "Normal QQ Plot - Rank Z Heart GSSG")
  
#For Heart Total GSH
  boxplot(pheno$Heart_TotalGSH, main = "Heart TotalGSH Box Plot")
  boxplot(pheno$Heart_TotalGSH~pheno$generation, main = "Heart TotalGSH Box Plot - by generation")
  boxplot(pheno$zHeartTotalGSH, main = "Rank Z Heart TotalGSH Box Plot")
  boxplot(pheno$zHeartTotalGSH~pheno$generation, main = "Rank Z Heart TotalGSH Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zHeartTotalGSH, main = "Normal QQ Plot - Rank Z Heart TotalGSH")
  
#For Heart GSH/GSSG Ratio
  boxplot(pheno$Heart_GSHGSSGRatio, main = "Heart GSHGSSGRatio Box Plot")
  boxplot(pheno$Heart_GSHGSSGRatio~pheno$generation, main = "Heart GSHGSSGRatio Box Plot - by generation")
  boxplot(pheno$zHeartGSHGSSGRatio, main = "Rank Z Heart GSHGSSGRatio Box Plot")
  boxplot(pheno$zHeartGSHGSSGRatio~pheno$generation, main = "Rank Z Heart GSHGSSGRatio Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zHeartGSHGSSGRatio, main = "Normal QQ Plot - Rank Z Heart GSHGSSGRatio")
  
#For Heart Eh
  boxplot(pheno$Heart_Eh, main = "Heart Eh Box Plot")
  boxplot(pheno$Heart_Eh~pheno$generation, main = "Heart Eh Box Plot - by generation")
  boxplot(pheno$zHeartEh, main = "Rank Z Heart Eh Box Plot")
  boxplot(pheno$zHeartEh~pheno$generation, main = "Rank Z Heart Eh Box Plot - by generation")
  #check if it is normally distributed
  qqnorm(pheno$zHeartEh, main = "Normal QQ Plot - Rank Z Heart Eh")
  

dev.off()

####################################################
## add covariates
####################################################

#because of how large my cohorts were, it is important to account for generation + sex

#adding sex and generation as covariates
sexgen = model.matrix(~ sex + generation, data = pheno)[,-1]


#For heritability calculation, need a linear mixed model
#####make kinship function using linear mixed model, not loco
#####default type of kinship is "overall" aka "linear mixed model" -- did not need to specify a type
  
kinship_lmm <- calc_kinship(probs = probs, use_allele_probs = TRUE, cores = 10)

#adding sex as covariate to compare to sexgen
sex = model.matrix(~ sex, data = pheno)[,-1] 



