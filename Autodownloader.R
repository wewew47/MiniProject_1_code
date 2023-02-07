library(stringr)
library(tibble)
library(seqinr)
library(rentrez)
library(dplyr)

autodownloader <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  proteinqueryfolder <- args[1]
  regularexpression <- args[2]
  #Get the database loaded as an object
  setwd("~/Documents/PhD/Year1/MiniProject_1/AMRFinderDatabase")
  RefGeneTablelong <- read.csv('refgenecatalog.csv', header = TRUE)
  
  #Dont need every column for what we want to do
  TableNames <- names(RefGeneTablelong)
  TableNames
  RefGeneTable <-
    subset.data.frame(RefGeneTablelong, select = c(TableNames[1:11]))
  
  refcheck <-
    function(rx, string)
      str_detect(string, rx) #function that returns true or false if the query is present
  
  #Find indices of the rows containing the query string. Check if each $allele entry contains the string
  targetindices <- c()
  for (i in 1:length(RefGeneTable$allele)) {
    if (refcheck(regularexpression, RefGeneTable$allele[i]) == TRUE) {
      targetindices <- c(targetindices, i)
    }
  }
  
  #Make the shortened table containing only rows of the query protein variants
  targets <- RefGeneTable[targetindices, ]
  
  ###Check if i have the wildtype file downloaded, if not the code should download it into the susceptible folder
  icount <- 0
  for (i in 1:length(unique(targets$whitelisted_taxa))) {
    #should be one file per species
    icount <- icount + 1
    species <- unique(targets$whitelisted_taxa)[i]
    print(paste('species is', species))
    print(i)
    #check for cases where $whitelisted_taxa lists genuses instead of species.
    if (str_detect(unique(targets$whitelisted_taxa)[i],
                   '[:alnum:]*_[:alnum:]') == FALSE) {
      print(unique(targets$whitelisted_taxa)[i])
      print('look here') #temp stataments to troubleshoot
      #obtain each refseq associated with this genus and store as a dataframe
      temp <- filter(targets, whitelisted_taxa == species)
      refseqs <- unique(temp$refseq_protein_accession)
      
      #Special cases handled in this loop
      for (j in 1:length(refseqs)) {
        print(j)
        print(substr(refseqs[j], start = 1, stop = 12))
        #retrieve seqs from entrez and save
        newfasta <-
          entrez_fetch(
            db = 'protein',
            id = substr(refseqs[j], start = 1, stop = 12),
            rettype = 'fasta'
          )
        write(
          newfasta,
          file = paste(
            "~/Documents/PhD/Year1/MiniProject_1/Protein_fastas/Susceptible/",proteinqueryfolder,"/",
            species,
            j,
            '.fasta',
            sep = ''
          )
        )
        if (j != 1) {
          icount <- icount + 1
        }
      }
    } #if the entry doesnt contain an underscore it is a genus, not a single species, so needs to be run multiple times as there may be multiple refseqs
    
    #Handling default cases if no file already exists
    else if (!file.exists(paste(species, '.fasta', sep = ''))) {
      print(substr(
        unique(targets$refseq_protein_accession)[icount],
        start = 1,
        stop = 12
      ))
      newfasta <- entrez_fetch(
        db = 'protein',
        id = substr(
          unique(targets$refseq_protein_accession)[icount],
          start = 1,
          stop = 12
        ),
        rettype = 'fasta'
      )
      write(
        newfasta,
        file = paste(
          "~/Documents/PhD/Year1/MiniProject_1/Protein_fastas/Susceptible/",proteinqueryfolder,"/",
          species,
          '.fasta',
          sep = ''
        )
      )
    }
    
  }
}

autodownloader()