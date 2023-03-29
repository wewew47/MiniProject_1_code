library(stringr)
setwd("~/Documents/PhD/Year1/MiniProject_1/Data/FoldXData/pbp/Average")
myfiles <- list.files()
allAvgs <- lapply(myfiles, read.delim, sep ='\t', header = TRUE, skip =8)

##Next need to separate the allAvgs into a separate variable for each species. 
#Need averages of all resistant within species, of all resistant across species, and of all wildtype across species
#compare averages of resistants against averages of species. 

#Extract individual datasets for gyrA
# Acinetobacter_baumannii <- as.data.frame(allAvgs[1])
# Burkholderia_cepacia <- as.data.frame(allAvgs[2])
# Neisseria_meningitidis <- as.data.frame(allAvgs[3])
# Burkholderia_pseudomallei <- as.data.frame(allAvgs[4])
# Neisseria_gonorrhoeae <- as.data.frame(allAvgs[5])
# Pseudomonas_aeruginosa <- as.data.frame(allAvgs[6])
# Staphylococcus_pseudintermedius <- as.data.frame(allAvgs[7])
# Campylobacter2 <- as.data.frame(allAvgs[8])
# Enterococcus_faecium <- as.data.frame(allAvgs[9])
# Clostridioides_difficile <- as.data.frame(allAvgs[10])
# Enterococcus_faecalis <- as.data.frame(allAvgs[11])
# Klebsiella_pneumoniae <- as.data.frame(allAvgs[12])
# Campylobacter1 <- as.data.frame(allAvgs[13])
# Staphylococcus_aureus <- as.data.frame(allAvgs[14])
# Escherichia1 <- as.data.frame(allAvgs[15])
# Salmonella1 <- as.data.frame(allAvgs[16])

#Extract datasets for pbp
Staphylococcus_aureus1 <- as.data.frame(allAvgs[1])
Staphylococcus_aureus1 <- Staphylococcus_aureus1[-c(2),]
Staphylococcus_aureus2 <- as.data.frame(allAvgs[2])
Staphylococcus_aureus2 <- Staphylococcus_aureus2[-c(6:10),]
Enterococcus_faecium <- as.data.frame(allAvgs[3])
Enterococcus_faecium <- Enterococcus_faecium[-c(3,4),]
Streptococcus_agalactiae <- as.data.frame(allAvgs[4])
Streptococcus_agalactiae <- Streptococcus_agalactiae[-c(2),]
Streptococcus_pyogenes <- as.data.frame(allAvgs[5])
Streptococcus_pyogenes <- Streptococcus_pyogenes[-c(2),]
#2020 Nature paper by Kalman et al. used as basis for categorising changes as (de)stabilising or not, and to what degree
#highly stabilising =  < -1.84 kcal/mol
#stabilising =  >= -1.84 && < -0.92 kcal/mol
#slightly stabilising =  >= -0.92 && < -0.46 kcal/mol
#neutral = >= -0.46 && <= 0.46 kcal/mol
#slightly destabilising =  <= +0.92 && > +0.46 kcal/mol
#destabilising =  <= +1.84 && > +0.92 kcal/mol
#highly destabilising =  > +1.84 kcal/mol



#Could make plots or heatmaps of the total.energy of each species 
setwd("~/Documents/PhD/Year1/MiniProject_1/Data/FoldXData/pbp/Raw/")
myRawfiles <- list.files()
allRaw <- lapply(myRawfiles, read.delim, sep ='\t', header = TRUE, skip =8)
Raw_Susceptible <- allRaw
Raw_Resistant <- allRaw
for(i in 1:length(allRaw)){
  print(i)
  sublist <- allRaw[[i]]
  Resistantrowcounter <- c()
  Susceptiblerowcounter <- c()
  for(j in 1:nrow(sublist)){
    print(j)
    print(sublist$Pdb[j])
    if(str_detect(sublist$Pdb[j], 'WT') == TRUE){
      print('TRUE')
      Resistantrowcounter <- c(Resistantrowcounter, j)
    }
    else if(str_detect(sublist$Pdb[j], 'WT') == FALSE){
      print('FALSE')
      Susceptiblerowcounter <- c(Susceptiblerowcounter, j)
    }
  }
  Raw_Susceptible[[i]] <- Raw_Susceptible[[i]][-Susceptiblerowcounter,]
  Raw_Resistant[[i]] <- Raw_Resistant[[i]][-Resistantrowcounter,]
}

Raw_Resistant <- t(cbind(t(Raw_Resistant[[1]]),t(Raw_Resistant[[2]]),t(Raw_Resistant[[3]]),
              t(Raw_Resistant[[4]]), t(Raw_Resistant[[5]])))

Raw_Susceptible <- t(cbind(t(Raw_Susceptible[[1]]),t(Raw_Susceptible[[2]]),t(Raw_Susceptible[[3]]),
                           t(Raw_Susceptible[[4]]), t(Raw_Susceptible[[5]])))

resistantmean <- mean(as.numeric(Raw_Resistant[,2]))
susceptiblemean <- mean(as.numeric(Raw_Susceptible[,2]))
resistantSD <- sd(as.numeric(Raw_Resistant[,2]))
susceptibleSD <- mean(as.numeric(Raw_Susceptible[,2]))
resistantmean
susceptiblemean
resistantSD
susceptibleSD
boxplot(as.numeric(Raw_Resistant[,2]), as.numeric(Raw_Susceptible[,2]),names = c('Resistant','Susceptible'))

testfunc <- function(number){
  boxplot(as.numeric(Raw_Resistant[,number]), as.numeric(Raw_Susceptible[,number]),names = c('Resistant','Susceptible'),
          xlab='dataset',ylab='Kcal/mol',main='Van der Waals Clashes')
}
testfunc(9) #Potentially more vdw clashes in some resistant structures than susceptible ones, but means are still extremely close


#Averages of everything dont look good across species,
#Next step is to try and look at within species averages.
Smean <- c()
SSD <- c()
RSD <- c()
Rmean <- c()

for(i in 1:length(allRaw)){
  sublist <- allRaw
}
test <- bind_rows(allAvgs, .id = "column_label")
#write.csv(test, 'gyrA_FoldX_Averages_Resistant.csv')
for(i in 1:length(allRaw)){
  sublist <- allRaw[[i]] ##Should save sublist as a .csv
  Resistantvalues <- c()
  Susceptiblevalues <- c()
  for(j in 1:nrow(sublist)){
    if(str_detect(sublist$Pdb[j], 'WT') == TRUE){
      Resistantvalues <- c(Resistantvalues, sublist[j,2])
    }
    else if(str_detect(sublist$Pdb[j], 'WT') == FALSE){
      Susceptiblevalues <- c(Susceptiblevalues, sublist[j,2])
    }
  }
  Smean <- c(Smean, mean(Susceptiblevalues))
  Rmean <- c(Rmean, mean(Resistantvalues))
  SSD <- c(SSD, sd(Susceptiblevalues))
  RSD <- c(RSD, sd(Resistantvalues))
  print(i)
  print(sublist$Pdb[i])
}

withinspecies <- as.data.frame(cbind(species =list.files('~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/pbp'), Smean, SSD, Rmean, RSD))

for(i in 1:nrow(withinspecies)){
  print(i)
  diff <- as.numeric(withinspecies$Rmean[i]) - as.numeric(withinspecies$Smean[i])
  withinspecies$diff[i] <- diff
}
withinspecies$species
withinspecies$species <- c('Staphylococcus aureus1', 'Staphylococcus aureus2', 'Enterococcus faecium', 
                           'Streptococcus agalactiae', 'Streptococcus pyogenes')




boxplot(as.numeric(withinspecies$diff))
par(las=2)
par(mar=c(10,5,2,2))
barplot(withinspecies$diff, names.arg=withinspecies$species, cex.names=0.5,
        main='Differences of means between Susceptible and Resistant')

#write.csv(withinspecies, 'gyrA_unfolding_energy_within_species_means.csv')
library(ggplot2)
ggplot(withinspecies,aes(x=species, y=diff))+
  ylim(-2, 0.5) +
  geom_bar(fill='black',stat='identity') +theme(plot.margin = unit(c(1,1,1,1.7),'cm'),axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
  xlab('Species')+ylab('Difference of means (kcal/mol)')+ggtitle('Differences of means of stability between Susceptible and Resistant')

comparator <- function(colnumber){
  Smean <- c()
  SSD <- c()
  RSD <- c()
  Rmean <- c()
  for(i in 1:length(allRaw)){
    sublist <- allRaw[[i]]
    Resistantvalues <- c()
    Susceptiblevalues <- c()
    headers <- colnames(sublist)
    heading <- headers[colnumber]
    for(j in 1:nrow(sublist)){
      if(str_detect(sublist$Pdb[j], 'WT') == TRUE){
        Resistantvalues <- c(Resistantvalues, sublist[j,colnumber])
      }
      else if(str_detect(sublist$Pdb[j], 'WT') == FALSE){
        Susceptiblevalues <- c(Susceptiblevalues, sublist[j,colnumber])
      }
    }
    Smean <- c(Smean, mean(Susceptiblevalues))
    Rmean <- c(Rmean, mean(Resistantvalues))
    SSD <- c(SSD, sd(Susceptiblevalues))
    RSD <- c(RSD, sd(Resistantvalues))
  }
  
  withinspecies <- as.data.frame(cbind(species =list.files('~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/pbp'), Smean, SSD, Rmean, RSD))
  
  for(i in 1:nrow(withinspecies)){
    diff <- as.numeric(withinspecies$Rmean[i]) - as.numeric(withinspecies$Smean[i])
    withinspecies$diff[i] <- diff
  }
  par(las=2)
  par(mar=c(10,5,2,2))
  p <- barplot(withinspecies$diff, names.arg=withinspecies$species, cex.names=0.5,
          main=paste('Differences of means of', heading, 'between Susceptible and Resistant'))
  #ggplot(mean.df,aes(x=species, y=diff))+
  #  geom_bar(fill='black',stat='identity') +theme(plot.margin = unit(c(1,1,1,1.7),'cm'),axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
  #  xlab('Species')+ylab('Difference of means (kcal/mol)')+ggtitle(paste('Differences of means', heading  ,'between Susceptible and Resistant using AlphaFold structures'))
  
}

for(i in 2:23){
  comparator(i)
}

comparator(16)
