library(stringr)
library(ggplot2)

setwd("~/Documents/PhD/Year1/MiniProject_1/Data/FoldXData/gyrA/AlphaFoldStabilityOutput/Resistant")
myfiles <- list.files()
allAvgs <- lapply(myfiles, read.delim, sep ='\t', header = FALSE, stringsAsFactors = FALSE)

names <- list.files()
colheads <- c('pdb', 'SD', 'total energy', 'Backbone Hbond', 'Sidechain Hbond', 
              'Van der Waals', 'Electrostatics', 'Solvation Polar', 'Solvation Hydrophobic', 
              'Van der Waals clashes', 'entropy sidechain', 'entropy mainchain', 'sloop_entropy', 'mloop_entropy', 
              'cis_bond', 'torsional clash', 'backbone clash', 'helix dipole', 'water bridge', 'disulfide',
              'electrostatic kon', 'partial covalent bonds', 'energy ionisation', 'Entropy Complex')
for(i in 1:length(allAvgs)){
  colnames(allAvgs[[i]]) <- colheads
  allAvgs[[i]]$pdb <- list.files()[i]
}

#Make combined dataframe
allAvgs <- bind_rows(allAvgs)
allAvgs <- cbind(allAvgs, species=NA)
#Add species column to combined dataframe
species = c(
  'Acinetobacter_baumannii', 'Burkholderia_cepacia', 'Neisseria_meningitidis',
  'Burkholderia_pseudomallei', 'Neisseria_gonorrhoeae', 'Pseudomonas_aeruginosa',
  'Staphylococcus_pseudintermedius', 'Campylobacter2', 'Enterococcus_faecium',
  'Clostridioides_difficile', 'Enterococcus_faecalis', 'Klebsiella_pneumoniae',
  'Campylobacter1', 'Staphylococcus_aureus', 'Escherichia1',
  'Salmonella1')
truespecies = c(
  'Acinetobacter baumannii', 'Burkholderia cepacia', 'Neisseria meningitidis',
  'Burkholderia pseudomallei', 'Neisseria gonorrhoeae', 'Pseudomonas aeruginosa',
  'Staphylococcus pseudintermedius', 'Campylobacter spp', 'Enterococcus faecium',
  'Clostridioides difficile', 'Enterococcus faecalis', 'Klebsiella pneumoniae',
  'Campylobacter jejuni', 'Staphylococcus aureus', 'Escherichia spp',
  'Salmonella spp')
for(i in 1:length(species)){
  for(j in 1:nrow(allAvgs)){
    if(str_detect(allAvgs$pdb[j], species[i]) == TRUE){
      allAvgs$species[j] <- truespecies[i]
    }
  }
}



#the object only contains averages for resistant runs
ResAvgs <- allAvgs

#Obtain susceptible averages
setwd("~/Documents/PhD/Year1/MiniProject_1/Data/FoldXData/gyrA/AlphaFoldStabilityOutput/Susceptible")
myfiles <- list.files()
allAvgs <- lapply(myfiles, read.delim, sep ='\t', header = FALSE, stringsAsFactors=FALSE)

names <- list.files()
colheads <- c('pdb', 'SD', 'total energy', 'Backbone Hbond', 'Sidechain Hbond', 
              'Van der Waals', 'Electrostatics', 'Solvation Polar', 'Solvation Hydrophobic', 
              'Van der Waals clashes', 'entropy sidechain', 'entropy mainchain', 'sloop_entropy', 'mloop_entropy', 
              'cis_bond', 'torsional clash', 'backbone clash', 'helix dipole', 'water bridge', 'disulfide',
              'electrostatic kon', 'partial covalent bonds', 'energy ionisation', 'Entropy Complex')
for(i in 1:length(allAvgs)){
  colnames(allAvgs[[i]]) <- colheads
  allAvgs[[i]]$pdb <- list.files()[i]
}

#Make combined dataframe
allAvgs <- bind_rows(allAvgs)
allAvgs <- cbind(allAvgs, species=NA)
#Add species column to combined dataframe
species = c(
  'Acinetobacter_baumannii', 'Burkholderia_cepacia', 'Neisseria_meningitidis',
  'Burkholderia_pseudomallei', 'Neisseria_gonorrhoeae', 'Pseudomonas_aeruginosa',
  'Staphylococcus_pseudintermedius', '2_Campylobacter', 'Enterococcus_faecium',
  'Clostridioides_difficile', 'Enterococcus_faecalis', 'Klebsiella_pneumoniae',
  '1_Campylobacter', 'Staphylococcus_aureus', '1_Escherichia',
  '1_Salmonella')
for(i in 1:length(species)){
  for(j in 1:nrow(allAvgs)){
    if(str_detect(allAvgs$pdb[j], species[i]) == TRUE){
      allAvgs$species[j] <- truespecies[i]
    }
  }
}
SusAvgs <- allAvgs



ResDiff <- ResAvgs
for(i in 1:nrow(ResDiff)){
  for(j in 1:nrow(SusAvgs)){
    if(str_detect(ResDiff$species[i], SusAvgs$species[j]) == TRUE){
      print('match')
      for(k in 2:(length(ResDiff)-1)){
        print(ResDiff[i,k] - SusAvgs[j, k])
        ResDiff[i,k] <- ResDiff[i,k] - SusAvgs[j, k]
      }
    }
  }
}

#Plot as a barplot
ggplot(ResDiff, aes(x=species, y=`total energy`)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle =45,vjust=1,hjust=1)) +
  ggtitle(paste('Mean ','\U0394','\U0394', 'G by species using AlphaFold Structures',sep='')) +
  ylab(paste('\U0394','\U0394','G',sep=''))


#Find mean difference per species and plot as bar plot
species <- SusAvgs$species
meanvec <- c()
class(meanvec)

for(i in 1:length(species)){
  tempvec <- c()
  for(j in 1:nrow(ResDiff)){
    if(str_detect(ResDiff$species[j], species[i]) == TRUE){
      tempvec <- c(tempvec, as.numeric(ResDiff$`total energy`[j]))
    }
  }
  meanvec <- c(meanvec, as.numeric(round(mean(tempvec), digits = 2)))
}

mean.df <- as.data.frame(cbind(as.numeric(meanvec), species), stringsAsFactors = FALSE)
mean.df$V1 <- as.numeric(mean.df$V1)


ggplot(mean.df,aes(x=species, y=meanvec))+
  geom_bar(fill='black',stat='identity') +theme(plot.margin = unit(c(1,1,1,1.7),'cm'),axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
  xlab('Species')+ylab('Difference of means (kcal/mol)')+ggtitle('Differences of means of stability between Susceptible and Resistant using AlphaFold structures')












