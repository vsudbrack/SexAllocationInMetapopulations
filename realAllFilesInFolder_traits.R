library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(DescTools)

# Goes through all mutInfo files in the same folder
folder = "/work/FAC/FBM/DEE/cmullon/default/sexAllocation/"
names = list.files(path=folder, 
                   pattern="Trait_*", 
                   full.names=TRUE, recursive=FALSE)

# Summary statistics is saved here
results = data.frame()

for (filename in names){
  mutpoints = read.table(filename, 
                         col.names=c("GEN", "MEAN.Z", "VARB.Z", "VARW.Z", "PI", "FIS", "FST", "FIT", "JOSHD", "GST"),
                         row.names = NULL)
  
  filename = substr(filename, 1, nchar(filename) - 4)
  row =  mutpoints
  
  pts = str_split(filename, "_")
  row$K  = as.numeric(gsub("[A-Z, a-z, =]", '', pts[[1]][str_detect(pts[[1]], 'K=' )]))
  row$E  = as.numeric(gsub("[A-Z, a-z, =]", '', pts[[1]][str_detect(pts[[1]], 'lifetimeDeme=' )]))
  row$MP  = as.numeric(gsub("[A-Z, a-z, =]", '', pts[[1]][str_detect(pts[[1]], 'dPolen=' )]))
  row$MS  = as.numeric(gsub("[A-Z, a-z, =]", '', pts[[1]][str_detect(pts[[1]], 'dSeeds=' )]))
  row$N  = as.numeric(gsub("[A-Z, a-z, =]", '', pts[[1]][str_detect(pts[[1]], 'N=' )]))
  row$SELF  = as.numeric(gsub("[A-Z, a-z, =]", '', pts[[1]][str_detect(pts[[1]], 'selfing=' )]))
  
  results = rbind(results, row)
}

write.csv(results, "results_SexAllocation.csv")
