# set up of the ssh
install.packages("ssh")
library("ssh", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

ssh_connect("k1632479@login.rosalind.kcl.ac.uk", passwd="Hkdogs222*")

# loading out of date packages 
install.packages("remotes")
library(remotes)

# 
install_bioc("airway", version = "3.8")
library("airway")
library("shh")
