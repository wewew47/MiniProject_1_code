library(stringr)
library(tibble)
library(seqinr)
library(rentrez)
library(dplyr)
setwd("~/Documents/PhD/Year1/MiniProject_1/Protein_fastas/Susceptible/gyrA")
RefGeneTablelong <- read.csv('refgenecatalog.csv', header = TRUE)


#Dont need every column for what we want to do
TableNames <- names(RefGeneTablelong)
RefGeneTable <- subset.data.frame(RefGeneTablelong, select = c(TableNames[1:11]))
#Make regex for the protein being queried.
refcheck <- function(rx,string) str_detect(string, rx) #function that returns true or false if the query is present
#Find indices of the rows containing the query string. Check if each $allele entry contains the string
targetindices <- c()
for(i in 1:length(RefGeneTable$allele)){
  if (refcheck("gyrA_[:upper:]\\d*[:upper:]", RefGeneTable$allele[i]) == TRUE) {
    targetindices <- c(targetindices, i)
  }
}
#Make the shortened table containing only rows of the query protein variants
targets <- RefGeneTable[targetindices,]
for(i in 1:length(targets$allele)){
  #if genus is found, not species, do the following
  if (str_detect(targets$whitelisted_taxa[i], '[:alnum:]*_[:alnum:]') == FALSE) {
    accessionlist <- c()
    genus = targets$whitelisted_taxa[i]
    genusindices <- str_which(targets$whitelisted_taxa, genus)
    temp <- targets[genusindices,]
    numbers = length(unique(temp$refseq_protein_accession))
    for(j in 1:numbers){
      accession <- unique(temp$refseq_protein_accession)[j]
      accessionlist <- c(accessionlist, accession)
    }
    for(k in 1:length(accessionlist)){
      if(targets$refseq_protein_accession[i] == accessionlist[k]){
        targets$whitelisted_taxa[i] = paste(targets$whitelisted_taxa[i],k,sep='')
      }
    }
  }
}
#add column for mutation identifier only (the allele col entries take the form 'genename_mutationidentifier')
targets <- add_column(targets, mutation = '', .after = 'allele')
print(length(targets$allele))
for(i in 1:length(targets$allele)){
  name <- targets$allele[i]
  newname <- ''
  newname = paste(substr(name, start=6,stop=6),'A', substr(name, start = 7, stop = nchar(name)),sep ='')
  targets$mutation[i] = newname 
}



mutantfilemaker <- function(species,j){
  setwd('~/Documents/PhD/Year1/MiniProject_1/Protein_fastas/Resistant/gyrA')
  files <- str_subset(list.files(), species)
  print(paste('species is',species))
  wildtype <- read.fasta(paste('~/Documents/PhD/Year1/MiniProject_1/Protein_fastas/Susceptible/gyrA/',j,'_',species,'.fasta',sep=''),as.string = TRUE, seqonly=TRUE)
  sequences <- c(wildtype)
  for(i in 1:length(files)){
    print(i)
    sequence <- read.fasta(files[i],as.string = TRUE,seqonly = TRUE)
    print(i)
    temp <- str_extract(files[i],'\\d+[:upper:]{1}[:digit:]{1,3}')
    removable <- str_extract(temp, '[:digit:]{1,3}[:upper:]{1}')
    coordinate <- as.numeric(gsub(removable, '',temp))
    insertable <- substr(sequence, coordinate-30, coordinate-1)
    sequences <- c(sequences,sequence)
  }
  write.table(sequences,file=paste("~/Documents/PhD/Year1/MiniProject_1/AMRFinderDatabase/mutant_file_",j,'_',species,'.pdb.txt',sep=''),sep="\n",dec = ";",row.names=FALSE,col.names =FALSE, quote = FALSE)
}

for(j in 1:length(unique(targets$whitelisted_taxa))){
  mutantfilemaker(unique(targets$whitelisted_taxa)[j],j)
}

mutantfilemaker('Staphylococcus_pseudintermedius',15)


mutantfilemakeralt <- function(species,species2,j,k){
  setwd('~/Documents/PhD/Year1/MiniProject_1/Protein_fastas/Resistant/gyrA')
  files <- str_subset(list.files(), species)
  print(paste('species is',species))
  wildtype <- read.fasta(paste('~/Documents/PhD/Year1/MiniProject_1/Protein_fastas/Susceptible/gyrA/',j,'_',k,'_',species2,'.fasta',sep=''),as.string = TRUE, seqonly=TRUE)
  sequences <- c(wildtype)
  for(i in 1:length(files)){
    print(i)
    sequence <- read.fasta(files[i],as.string = TRUE,seqonly = TRUE)
    print(i)
    temp <- str_extract(files[i],'\\d+[:upper:]{1}[:digit:]{1,3}')
    removable <- str_extract(temp, '[:digit:]{1,3}[:upper:]{1}')
    coordinate <- as.numeric(gsub(removable, '',temp))
    insertable <- substr(sequence, coordinate-30, coordinate-1)
    sequences <- c(sequences,sequence)
  }
  write.table(sequences,file=paste("~/Documents/PhD/Year1/MiniProject_1/AMRFinderDatabase/mutant_file_",j,'_',k,'_',species2,'.pdb.txt',sep=''),sep="\n",dec = ";",row.names=FALSE,col.names =FALSE, quote = FALSE)
}

mutantfilemakeralt('Salmonella1','Salmonella',13,1)
