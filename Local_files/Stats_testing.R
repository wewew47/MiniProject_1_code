options(scipen = 999)
setwd("~/Documents/PhD/Year1/MiniProject_1/Results/Tables")

Resistant <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/gyrA_FoldX_Averages_Resistant.csv')[,-c(1)]
Random <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/MutateXAveragesData.csv')[,-c(1)]

names <- list.files('~/Documents/PhD/Year1/MiniProject_1/Data/ProFitOutput/gyrA')
for(i in 1:nrow(Resistant)){
  for(j in 1:length(names)){
    print(paste('Resistant', Resistant$column_label[i]))
    print(paste('j', j))
    #print(names[j])
    if(Resistant$column_label[i] == j){
      print("TRUE")
      Resistant$column_label[i] = names[j]
    }
  }
}

Resistant <- Resistant[,c(1,4)]
Random <- Random[,c(1,2)]
colnames(Resistant) <- c('species', 'unfolding_energy')
colnames(Random) <- c('species', 'unfolding_energy')

Resistant$Dataset <- 'Resistant'
Random$Dataset <- 'Random'

x <- t.test(Resistant$unfolding_energy[Resistant$species==names[2]], Random$unfolding_energy[Random$species==names[2]])

pvalues <- c()
for(i in 1:length(names)){
  print(i)
  print(names[i])
  pvalue <- t.test(Resistant$unfolding_energy[Resistant$species==names[i]], Random$unfolding_energy[Random$species==names[i]])
  pvalues <- c(pvalues, pvalue$p.value)
}

significance1 <- as.data.frame(cbind(names, pvalues))

write.csv(significance1,'~/Documents/PhD/Year1/MiniProject_1/Results/Tables/gyrA_FoldX_Resistant_Vs_Random_Stats.csv')



###Profit analysis
Random <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/RMS_Random_Combined.csv')[,-c(1)]
Resistant <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/RMS_Resistant_Combined.csv')[,-c(1)]
colnames(Random) <- c('Sample', 'RMS','Species')
colnames(Resistant) <- c('Sample', 'RMS','Species')
Random$Species[Random$Species=="5_Clostridioides_difficule"] <- "5_Clostridioides_difficile"
names <- unique(sort(Random$Species))
x <- t.test(Resistant$RMS[Resistant$Species==names[2]], Random$RMS[Random$Species==names[2]])
names[2]
names

pvalues <- c()
for(i in 1:length(names)){
  print(i)
  print(names[i])
  pvalue <- t.test(Resistant$RMS[Resistant$Species==names[i]], Random$RMS[Random$Species==names[i]])
  pvalues <- c(pvalues, pvalue$p.value)
}
options(scipen=999)
significance2 <- as.data.frame(cbind(names, pvalues))

write.csv(significance2,'~/Documents/PhD/Year1/MiniProject_1/Results/Tables/gyrA_ProFit_Resistant_Vs_Random_Stats.csv')

