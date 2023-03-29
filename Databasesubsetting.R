setwd("~/Documents/PhD/Year1/MiniProject_1/AMRFinderDatabase")
RefGeneTablelong <- read.csv('refgenecatalog.csv', header = TRUE)
#Dont need every column for what we want to do
TableNames <- names(RefGeneTablelong)
RefGeneTable <- subset.data.frame(RefGeneTablelong, select = c(TableNames[1:11]))

targetindices <- c()
for(i in 1:length(RefGeneTable$allele)){
  if (str_detect(RefGeneTable$subtype[i], 'POINT') == TRUE) {
    targetindices <- c(targetindices, i)
  }
}
#Make the shortened table containing only rows of the query protein variants
RefGeneTablePoints<- RefGeneTable[targetindices,]
length(unique(RefGeneTablePoints$gene_family))
RefGeneTablePoints <- RefGeneTablePoints[-(1:75),]
#parC and parE could be a useful combo to look at. treat them as the same ie. analyse at same time
#also the pbp proteins (penicillin binding proteins)
#gyrB is present in a few species
#rpoB
target_finder <- function(query){
  targetindices <- c()
  for(i in 1:length(RefGeneTablePoints$allele)){
    if (str_detect(RefGeneTablePoints$gene_family[i], query) == TRUE) {
      targetindices <- c(targetindices, i)
    }
  }
  assign(paste('targets_',query,sep=''), RefGeneTablePoints[targetindices,])
  targets_par <- RefGeneTablePoints[targetindices,]
  return(targets_par)
}

targets_par <- target_finder('par')
targets_pbp <- target_finder('pbp')
targets_rpoB <- target_finder('rpoB')
targets_gyrB <- target_finder('gyrB')

write.csv(targets_par, '~/Documents/PhD/Year1/MiniProject_1/Data/Databases/par.csv')
write.csv(targets_pbp, '~/Documents/PhD/Year1/MiniProject_1/Data/Databases/pbp.csv')
write.csv(targets_rpoB, '~/Documents/PhD/Year1/MiniProject_1/Data/Databases/rpoB.csv')
write.csv(targets_gyrB, '~/Documents/PhD/Year1/MiniProject_1/Data/Databases/gyrB.csv')

