library(dplyr)
library(stringr)

#Code to get coordiantes for selecting mutatex files
targets <- read.csv('/gpfs/home/wms/lfrkmq/MiniProject_1/AMRFinderDatabase/gyrAdatabase.csv')
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



setwd('/gpfs/home/wms/lfrkmq/MiniProject_1/ProFitOutput/ProFitOutput/gyrA/Resistant')

names <- c()
dflist <- vector(mode = 'list', length = 5)
for(i in 1:length(list.files())){
  names <- c(names, list.files()[i])
  print(list.files()[i])
  setwd(paste("/gpfs/home/wms/lfrkmq/MiniProject_1/ProFitOutput/ProFitOutput/gyrA/Resistant/", list.files()[i], sep=''))
  RMSvector <- c()
  variantvector <- c()
  resistantcoord <- means[i]
  for(j in 1:length(list.files())){
      initialdata = readLines(list.files()[j])
      identifier = FALSE
      for(k in 1:length(initialdata)){
        if(str_detect(initialdata[k], 'RMS') == TRUE){
          identifier = TRUE
          RMSline <- initialdata[k]
        }
      }
    if(identifier == TRUE){
      RMS = as.numeric(strsplit(RMSline, ' ')[[1]][5])
      RMSvector = c(RMSvector, RMS)
      variantvector <- c(variantvector, j)
      nam <- paste('RMS', names[i], sep ='')
      assign(nam, as.data.frame(cbind(variantvector, RMSvector, names[i])))
      dflist[[i]] <- as.data.frame(cbind(variantvector, RMSvector, names[i]))
      identifier = FALSE
    }
  }
  setwd("/gpfs/home/wms/lfrkmq/MiniProject_1/ProFitOutput/ProFitOutput/gyrA/Resistant")
}

RMS_Resistant_Combined <- bind_rows(dflist)
write.csv(RMS_Resistant_Combined, '/gpfs/home/wms/lfrkmq/MiniProject_1/pbp_RMS_Resistant_Combined.csv')

setwd('/gpfs/home/wms/lfrkmq/MiniProject_1/ProFitOutput/ProFitOutput/gyrA/RandomAll')

names <- c()
dflist <- vector(mode = 'list', length = 5)
for(i in 1:length(list.files())){
  names <- c(names, list.files()[i])
  print(list.files()[i])
  setwd(paste("/gpfs/home/wms/lfrkmq/MiniProject_1/ProFitOutput/ProFitOutput/gyrA/RandomAll/", list.files()[i], sep=''))
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
    #filecoords <- as.numeric(filecoords)
    print(filecoords)
    #for each possible substitution at that coordinate, do the following
    if(filecoords > (resistantcoord-15) && filecoords < (resistantcoord+15)){
      print(j)
      initialdata = readLines(list.files()[j])
      identifier = FALSE
      for(k in 1:length(initialdata)){
        if(str_detect(initialdata[k], 'RMS') == TRUE){
          identifier = TRUE
          RMSline <- initialdata[k]
        }
      }
    }
    if(identifier == TRUE){
      RMS = as.numeric(strsplit(RMSline, ' ')[[1]][5])
      RMSvector = c(RMSvector, RMS)
      variantvector <- c(variantvector, j)
      nam <- paste('RMS', names[i], sep ='')
      assign(nam, as.data.frame(cbind(variantvector, RMSvector, names[i])))
      dflist[[i]] <- as.data.frame(cbind(variantvector, RMSvector, names[i]))
      identifier = FALSE
    }
  }
  setwd("/gpfs/home/wms/lfrkmq/MiniProject_1/ProFitOutput/ProFitOutput/gyrA/RandomAll")
}

RMS_Random_Combined <- bind_rows(dflist)
write.csv(RMS_Random_Combined, '/gpfs/home/wms/lfrkmq/MiniProject_1/gyrA_RMS_Random_Combined_Subset.csv')

