MutateXAnalysis <- function(){
  #Code to get coordiantes for selecting mutatex files
  targets <- read.csv('~/Documents/PhD/Year1/MiniProject_1/Data/Databases/gyrAdatabase.csv')
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
  #setwd to be where all the species folders are located
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/MutateXInputFinal/results/mutation_ddgs")
  #number of species files is the directory size minus 1 as there is a final avg folder
  num <- (length(list.files()))
  #Initialise vector to store mean total.energy for each species
  speciesmeans <- c()
  speciesvec <- c()
  #For each species file do the following
  out = vector(mode = 'list', length = num+1)
  for(i in 1:num){
    #save species name so i can record the averages per species
    species <- list.files()[i]
    speciesvec <- c(speciesvec, species)
    resistantcoord <- means[i]
    print(means[i])
    print(means)
    #change to that species' folder
    setwd(list.files()[i])
    #Initialise vectors to store the data for all files of a species
    avgs_vec <- c()
    std_vec <- c()
    min_vec <- c()
    max_vec <- c()
    #For each file (i.e. each position in the protein), do the following
    for(j in 1:length(list.files())){
      #read file line by line
      currentfile <- read.delim(list.files()[j], header=FALSE, skip =1, sep= ' ')
      filename <- list.files()[j]
      #Get coords of that file
      filecoords <- substr(filename, start=3, stop =nchar(filename))
      #for each possible substitution at that coordinate, do the following
      if(filecoords > (resistantcoord-15) && filecoords < (resistantcoord+15)){
        for(k in 1:20){
          #one line in each file will be the original residue at that position, so ignore that to avoid skewing the mean
          if(currentfile$V1[k] != 0.00000 && currentfile$V1[k] != -0.00000){
            avgs_vec <- c(avgs_vec, currentfile$V1[k])
            std_vec <- c(std_vec, currentfile$V2[k])
            min_vec <- c(min_vec, currentfile$V3[k])
            max_vec <- c(max_vec, currentfile$V4[k])
          }
        }
      }
    } #End of run for that species
    df <- as.data.frame(cbind(avgs_vec, std_vec, min_vec, max_vec))#Bind each vector into a dataframe
    out[[i]]=df #Make list of vectors, one index per species 
    setwd("~/Documents/PhD/Year1/MiniProject_1/Data/MutateXInputFinal/results/mutation_ddgs")
  } #End of main for loop
  out[[num+1]] = speciesvec
  return(out)
} #End

x <- MutateXAnalysis()

namevec <- x[[3]]
names <- list.files('~/Documents/PhD/Year1/MiniProject_1/Data/ProFitOutput')
df_list <- x[1:(length(x)-1)]
#To plot boxplots ideally all data would be together in on dataframe
combineddata <- bind_rows(x[1:(length(x)-1)], .id='species')
for(i in 1:(length(namevec))){
  for(j in 1:length(combineddata$species)){
    if(combineddata$species[j] == i){
      combineddata$species[j] <- names[i]
    }
  }
}

ggplot(combineddata, aes(x=species, y=avgs_vec)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle =45,vjust=1,hjust=1)) +
  ggtitle(paste('Mean ','\U0394','\U0394', 'G by species',sep='')) +
  ylab(paste('\U0394','\U0394','G',sep=''))

#Obtain a dataframe of mean per species
meansvec <- c()
for(i in 1:length(df_list)){
  meansvec <- c(meansvec, mean(df_list[[i]]$avgs_vec))
}

df_mean <- as.data.frame(cbind(names[1:2], meansvec))
colnames(df_mean) <- c('Species', paste('Mean ', '\U0394', '\U0394', 'G', sep=''))

ggplot(df_mean,aes(x=Species, y=`Mean ΔΔG`))+
  geom_bar(fill='black',stat='identity') +theme(plot.margin = unit(c(1,1,1,1.7),'cm'),axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
  xlab('Species')+ylab('Difference of means (kcal/mol)')+ggtitle('Differences of means of stability between Susceptible and Resistant')

###Each line in each mutateX output result file is the total.energy of one mutation for one position for one wildtype protein
#Each reuslt ifle has 20 lines, for each amino acid that could be at that position in that structure
#Each species folder has as many files as there are coordinates in the protein





