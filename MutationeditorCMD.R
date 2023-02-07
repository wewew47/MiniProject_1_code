library(stringr)
library(tibble)
library(seqinr)
library(rentrez)
library(dplyr)
####Solution to the mutationeditor not reading genus file name problem:
####Change the column names in 'targets' to equal the filenames
#if species name lacks '_', then make the table consist of only as many rows as there are variants for that genus, so for campy its two rows
#then name the first instance genus1 and the seocnd genus2 etc
mutationeditor <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  proteinqueryfolder <- args[1]
  regularexpression <- args[2]
  print(paste('Regular expression is', regularexpression))
  #Get the database loaded as an object
  setwd("~/home/wms/lfrkmq/MiniProject_1/AMRFinderDatabase")
  RefGeneTablelong <- read.csv('refgenecatalog.csv', header = TRUE)
  #Dont need every column for what we want to do
  TableNames <- names(RefGeneTablelong)
  RefGeneTable <- subset.data.frame(RefGeneTablelong, select = c(TableNames[1:11]))
  #Make regex for the protein being queried.
  refcheck <- function(rx,string) str_detect(string, rx) #function that returns true or false if the query is present
  #Find indices of the rows containing the query string. Check if each $allele entry contains the string
  targetindices <- c()
  for(i in 1:length(RefGeneTable$allele)){
    if (refcheck(regularexpression, RefGeneTable$allele[i]) == TRUE) {
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
    newname = substr(name, start = 6, stop = nchar(name))
    targets$mutation[i] = newname 
  }
  #Using the entry in the mutation column, copy that species' wildtype sequence and introduce the point mutation. Saving as new file
  for(i in 1:length(targets$mutation)){
    species <- targets$whitelisted_taxa[i] #Identifier for future file organisation if i choose to organise files by species
    sequence <- read.fasta(paste("~/home/wms/lfrkmq/MiniProject_1/Protein_fastas/Susceptible/",proteinqueryfolder,"/", targets$whitelisted_taxa[i], '.fasta',sep = ''),
                           seqtype = 'AA', as.string = TRUE, seqonly = TRUE)
    
    #identify the original residue so i can check this base is present at the coordinate in the sequence
    original <- substr(targets$mutation[i],start=1,stop=1)
    #find coordinate in the sequence where the PM occurs
    coordinate <- as.numeric(substr(targets$mutation[i], start = 2, stop = (nchar(targets$mutation[i])-1)))
    
    #Check for original base not matching the base at the coordinate. Stop if this happens
    if(substr(sequence,start=coordinate,stop=coordinate) != original){   ###IF there are two identical adjacent bases this will mutate the wrong one!!!
      print(paste('Error position is', i))
      print('Error, sequence original position doesnt match the mutation original position')
      print('Increasing coordinate by 1')
      coordinate = coordinate + 1
      if(substr(sequence,start=coordinate,stop=coordinate) != original){
        print(paste('Error position is', i))
        print('Error, sequence original position doesnt match the new mutation position')
      }
    }
    else{
      #Obtain the residue that will replace original and make the mutant sequence
      pointmutant <- substr(targets$mutation[i],start = (nchar(targets$mutation[i])), stop = (nchar(targets$mutation[i])))
      mutantsequence <- paste(substr(sequence,start=1,stop=(coordinate-1)), pointmutant,substr(sequence,start=(coordinate+1),stop = nchar(sequence)), sep = '')
      write.fasta(sequences = mutantsequence, names = proteinqueryfolder, file.out = paste("~/home/wms/lfrkmq/MiniProject_1/Protein_fastas/Resistant/" ,proteinqueryfolder, '/', targets$mutation[i],'_', species, '.fasta', sep = ''))
    }
  }
  print('Finished')
}

mutationeditor()
