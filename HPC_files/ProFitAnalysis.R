library(stringr)
library(ggplot2)
library(dplyr)
library(tibble)
library(seqinr)
library(rentrez)


setwd("~/Documents/PhD/Year1/MiniProject_1/Data/ProFitOutputFull/pbp/ProFit_pbp/RandomAll")
names <- c()
dflist <- vector(mode = 'list', length = 5)
for(i in 1:length(list.files())){
  names <- c(names, list.files()[i])
  names <- c(  'Staphylococcus aureus (pbp2)', 'Staphylococcus aureus (pbp4)', 'Enterococcus faecium (pbp5)', 'Streptococcus agalactiae (pbp2x)', 'Streptococcus pyogenes (pbp2x)')
  print(paste("~/Documents/PhD/Year1/MiniProject_1/Data/ProFitOutputFull/pbp/ProFit_pbp/RandomAll/", list.files()[i], sep=''))
  setwd(paste("~/Documents/PhD/Year1/MiniProject_1/Data/ProFitOutputFull/pbp/ProFit_pbp/RandomAll/", list.files()[i], sep=''))
  RMSvector <- c()
  variantvector <- c()
  for(j in 1:length(list.files())){
    print(j)
    initialdata = readLines(list.files()[j])
    for(k in 1:length(initialdata)){
      print(initialdata[k])
      if(str_detect(initialdata[k], 'RMS') == TRUE){
        RMSline = initialdata[k]
      }
    }
    RMS = as.numeric(strsplit(RMSline, ' ')[[1]][5])
    RMSvector = c(RMSvector, RMS)
    variantvector <- c(variantvector, j)
    nam <- paste('RMS', names[i], sep ='')
    assign(nam, as.data.frame(cbind(variantvector, RMSvector)))
    dflist[[i]] <- as.data.frame(cbind(variantvector, RMSvector))
  }
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/ProFitOutputFull/pbp/ProFit_pbp/RandomAll")
}

meanvec <- c()
for(i in 1:length(dflist)){
  df <- dflist[[i]]
  meanvec <- c(meanvec, mean(df$RMSvector))
}

RMSmeans <- as.data.frame(cbind(names, meanvec))
colnames(RMSmeans) <- c('Species','Mean RMS')

#Also need to make boxplots of each dataframe

for(i in 1:length(dflist)){
  df <- dflist[[i]]
  ggplot(df, aes(y=RMSvector))+
    geom_boxplot()
}

combinedRMS <- bind_rows(dflist, .id='species')
for(i in 1:16){
  print(i)
  for(j in 1:nrow(combinedRMS)){
    print(j)
    if(combinedRMS$species[j] == i){
      combinedRMS$species[j] <- names[i]
    }
  }
}

colnames(combinedRMS) <- c('Species', 'Variant vector', 'RMS')

database <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Data/Databases/pbp.csv')[,-c(1)]
for(i in 1:length(database$allele)){
  database$mutation[i] <- substr(database$allele[i], start = 6, stop = nchar(database$allele[i]))
}

combinedRMS$mutation <- database$mutation





write.csv(combinedRMS, "~/Documents/PhD/Year1/MiniProject_1/Results/Tables/pbp_RandomRMS.csv")
write.csv(RMSmeans, "~/Documents/PhD/Year1/MiniProject_1/Results/Tables/pbp_mean_RandomRMS.csv")



is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

Resistant <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/pbp_ResistantRMS.csv')
Resistant$Species
combinedRMS %>%
  group_by(Species) %>%
  mutate(outlier=ifelse(is_outlier(RMS),as.character(mutation), as.character(NA))) %>%
  ggplot(., aes(x=Species, y=RMS)) +
    geom_boxplot() + theme(axis.text.x = element_text(angle =45,vjust=1,hjust=1)) +
    ggtitle('RMS deviation of resistant gyrA proteins by species') + geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.5, size = 3) +
  ylab('RMS')


#+ geom_text(aes(label=ifelse(RMSvector>0.15,variantvector,'')),hjust = 1.4, vjust = -0.1) +
#  geom_text(aes(label=ifelse(RMSvector<0.15,variantvector,'')),hjust = 1.4, vjust = -0.1) 



Random <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/gyrA_RMS_Random_Combined_Subset.csv')[,-c(1)]
Resistant <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/RMS_Resistant_Combined.csv')[,-c(1)]
colnames(Random) <- c('Sample', 'RMS','Species')
colnames(Resistant) <- c('Sample', 'RMS','Species')
Random$Species[Random$Species=="5_Clostridioides_difficule"] <- "5_Clostridioides_difficile"
names <- unique(sort(Random$Species))
names <- c('Acinetobacter baumannii', 'Neisseria gonorrhoeae', 'Pseudomonas aeruginosa', 'Salmonella spp.', 'Staphylococcus aureus',
           'Staphylococcus pseudintermedius', 'Burkholderia cepacia', 'Burkholderia pseudomallei', 'Campylobacter jejuni', 'Campylobacter spp.', 'Clostridioides difficile',
           'Enterococcus faecalis', 'Enterococcus faecium', 'Escherichia spp.', 'Klebsiella pneumoniae')
p <- ggplot(Random, aes(x=Species, y=RMS)) +
  geom_boxplot() + 
  theme(plot.margin = unit(c(1,1,1,1.7),'cm'),axis.text.x=element_text(angle=45,vjust=1,hjust=1))

combinedRMS <- bind_rows(Random, Resistant)

for(i in 1:length(combinedRMS$Species)){
  print(i)
  for(j in 1:length(names)){
    print(j)
    name <- names[j]
    sub <- substr(name, start=nchar(name)-6, stop=nchar(name))
    if(str_detect( combinedRMS$Species[i], sub) == TRUE){
      combinedRMS$Species[i] <- names[j]
    }
  }
}

for(i in 1:length(combinedRMS$Species)){
  sample = combinedRMS$Species[i]
  if (sample == '11_Neisseria_meningitidis'){
    combinedRMS$Species[i] <- 'Neisseria meningitidis'
  }
  if( sample == '13_1_Salmonella'){
    combinedRMS$Species[i] <- 'Salmonella spp'
  }
  if ( sample == '14_Staphylococcus_aureus'){
    combinedRMS$Species[i] <- 'Staphylococcus auerus'
  }
  if ( sample == '4_1_Campylobacter'){
    combinedRMS$Species[i] <- 'Campylobacter jejuni'
  }
  if (sample == '4_2_Campylobacter'){
    combinedRMS$Species[i] <- 'Campylobacter spp'
  }
  if ( sample == '8_1_Escherichia'){
    combinedRMS$Species[i] <- 'Escherichia spp'
  }
}

unique(sort(combinedRMS$Species))

combinedRMSnew <- combinedRMS[-c(8832:8841),]
combinedRMSnew %>%
  group_by(Species) %>%
  ggplot(., aes(x=Species, y=RMS)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle =45,vjust=1,hjust=1)) +
  ggtitle('RMS deviation of gyrA proteins by species') +
  ylab('RMS')











Random <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/pbp_RMS_Random_Combined_Subset.csv')[,-c(1)]
colnames(Random) <- c('Sample', 'RMS','Species')
unique(Random$Species)

names <- unique(sort(Random$Species))
names
truenames <- c('Staphylococcus aureus (pbp2)', 'Staphylococcus aureus (pbp4)', 'Enterococcus faecium (pbp5)', 'Streptococcus agalactiae (pbp2x)', 'Streptococcus pyogenes (pbp2x)')
p <- ggplot(Random, aes(x=Species, y=RMS)) +
  geom_boxplot() + 
  theme(plot.margin = unit(c(1,1,1,1.7),'cm'),axis.text.x=element_text(angle=45,vjust=1,hjust=1))
p
for(i in 1:length(Random$Species)){
  print(i)
  for(j in 1:length(names)){
    print(j)
    name <- names[j]
    if(str_detect( Random$Species[i], name) == TRUE){
      Random$Species[i] <- truenames[j]
    }
  }
}


unique(sort(Random$Species))

Random %>%
  group_by(Species) %>%
  ggplot(., aes(x=Species, y=RMS)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle =45,vjust=1,hjust=1)) +
  ggtitle('RMS deviation of gyrA proteins by species') +
  ylab('RMS')





Resistant <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Results/Tables/pbp_RMS_Resistant_Combined.csv')[,-c(1)]
colnames(Resistant) <- c('Sample', 'RMS','Species')

names <- unique(sort(Resistant$Species))
names <- c('Staphylococcus aureus (pbp2)', 'Staphylococcus aureus (pbp4)', 'Enterococcus faecium (pbp5)', 'Streptococcus agalactiae (pbp2x)', 'Streptococcus pyogenes (pbp2x)')
p <- ggplot(Resistant, aes(x=Species, y=RMS)) +
  geom_boxplot() + 
  theme(plot.margin = unit(c(1,1,1,1.7),'cm'),axis.text.x=element_text(angle=45,vjust=1,hjust=1))
p


for(i in 1:length(Resistant$Species)){
  print(i)
  for(j in 1:length(names)){
    print(j)
    name <- names[j]
    sub <- substr(name, start=nchar(name)-6, stop=nchar(name))
    if(str_detect( Resistant$Species[i], sub) == TRUE){
      Resistant$Species[i] <- names[j]
    }
  }
}

unique(sort(Resistant$Species))

Resistant %>%
  group_by(Species) %>%
  mutate(outlier=ifelse(is_outlier(RMS),as.character(Sample), as.character(NA))) %>%
  ggplot(., aes(x=Species, y=RMS)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle =45,vjust=1,hjust=1)) +
  ggtitle('RMS deviation of gyrA proteins by species') +
  ylab('RMS')
