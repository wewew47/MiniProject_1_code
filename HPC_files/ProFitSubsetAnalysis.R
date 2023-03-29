library(stringr)
library(ggplot2)
library(dplyr)
library(tibble)
library(seqinr)
library(rentrez)

#Code to get coordiantes for selecting mutatex files
targets <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Data/Databases/pbp.csv')
#Now need the range of mutations for each species and take a the mean position
means <- c()
for(i in 1:length(unique(targets$whitelisted_taxa))){
  species <- sort(unique(targets$whitelisted_taxa))
  mutationrange <- c()
  for(j in 1:nrow(targets)){
    if(targets$whitelisted_taxa[j] == species[i]){
      mutationcoord <- substr(targets$allele[j], start=7, stop=(nchar(targets$allele[j])-1))
      mutationrange <- c(mutationrange, as.numeric(mutationcoord))
    }
  }
  means <- c(means, round(mean(mutationrange)))
}
print(paste('mean mutations are:', means))


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
  resistantcoord <- means[i]
  for(j in 1:length(list.files())){
    currentfile <- read.delim(list.files()[j], header=FALSE, skip =1, sep= ' ')
    filename <- list.files()[j]
    #Get coords of that file
    coord_location <- str_locate(list.files()[j], '[:alnum:][:upper:][:digit:]{1,3}')
    coord_location <- as.numeric(coord_location)
    filecoords <- substr(filename, start=coord_location[1], stop =coord_location[2])
    print(class(filecoords))
    #filecoords <- as.numeric(filecoords)
    print(filecoords)
    #for each possible substitution at that coordinate, do the following
    if(filecoords > (resistantcoord-15) && filecoords < (resistantcoord+15)){
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
