library(stringr)
library(ggplot2)

setwd('~/Documents/PhD/Year1/MiniProject_1/Results/Tables')
data1 <- read.csv('pbp_MutateXSubsetMeans.csv')[1:16, ]
data2 <- read.csv('pbp_MutateXSubsetAveragesData.csv')
colnames(data1)

species = unique(sort(data2$species))
species
truespecies = c('Staphylococcus aureus (pbp2)', 'Staphylococcus aureus (pbp4)', 'Enterococcus faecium (pbp5)', 'Streptococcus agalactiae (pbp2x)', 'Streptococcus pyogenes (pbp2x)')

truespecies
for(i in 1:570){
  data2$species[i] <- truespecies[1]
}

for(i in 571:1146){
  data2$species[i] <- truespecies[2]
}

for(i in 1147:1722){
  data2$species[i] <- truespecies[3]
}

for(i in 1723:2292){
  data2$species[i] <- truespecies[4]
}

for(i in 2293:2870){
  data2$species[i] <- truespecies[5]
}
# 
# for(i in 1:length(species)){
#   for(j in 1:nrow(data2)){
#     if(str_detect(data2$species[j], species[i]) == TRUE){
#       if(data2$species[j] == '10'){
#         print(species[i])
#         print(truespecies[i])
#         print(data2$species[j])
#       }
#       data2$species[j] <- truespecies[i]
#       if(data2$species[j] == 'Enterococcus faecium (pbp5)'){
#         print(data2$species[j])
#         print('')
#       }
#     }
#   }
# }

colnames(data2) <- c('X' ,'avgs_vec', 'std_vec','min_vec', 'max_vec','Species')
unique(sort(data2$Species))

ggplot(data2, aes(x=Species, y=`avgs_vec`)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle =45,vjust=1,hjust=1)) +
  ggtitle(paste('\U0394','\U0394', 'G by species',sep='')) +
  ylab(paste('\U0394','\U0394','G',sep=''))



setwd("~/Documents/PhD/Year1/MiniProject_1/Data/FoldXData/gyrA/AlphaFoldStabilityOutput/Susceptible")
myfiles <- list.files()
AllSus <- lapply(myfiles, read_tsv, c(header=FALSE))
for(i in 1:length(AllSus)){
  AllSus[[i]]$X1 <- myfiles[i]
}

setwd("~/Documents/PhD/Year1/MiniProject_1/Data/FoldXData/gyrA/AlphaFoldStabilityOutput/Resistant")
myfiles <- list.files()
AllRes <- lapply(myfiles, read_tsv, c(header=FALSE))
for(i in 1:length(AllRes)){
  AllRes[[i]]$X1 <- myfiles[i]
}

