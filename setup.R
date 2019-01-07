# set up of the ssh
install.packages("ssh")
library("ssh", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

session <- ssh_connect("k1632479@login.rosalind.kcl.ac.uk", passwd="Hkdogs222*")


# loading out of date packages 
install.packages("remotes")
library(remotes)

# loading out of date but the experimental data file
install_bioc("airway", version = "3.8")
library("airway")

#install ggbeeswarm
install.packages("ggbeeswarm")

#install GTF file downloader 
install.packages("GenomicFeatures")

#loading STAR
# install.packages("STAR")

# install
source("https://bioconductor.org/biocLite.R")
biocLite("vsn")

install.packages("dplyr")
install.packages("ggplot2")
install.packages("hexbin")
install.packages("pheatmap")
