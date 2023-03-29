library(stringr)
library(tibble)
library(seqinr)
library(rentrez)
library(dplyr)


Dloader <- function(protein){
  setwd('~/Documents/PhD/Year1/MiniProject_1/Data/Databases/')
  data <- read.csv(paste(protein, '.csv', sep=''))
  ###Check if i have the wildtype file downloaded, if not the code should download it into the susceptible folder
  icount <- 0
  for (i in 1:length(unique(data$whitelisted_taxa))) {
    #should be one file per species
    icount <- icount + 1
    species <- unique(data$whitelisted_taxa)[i]
    print(paste('species is', species))
    print(i)
    #check for cases where $whitelisted_taxa lists genuses instead of species.
    if (str_detect(unique(data$whitelisted_taxa)[i],
                   '[:alnum:]*_[:alnum:]') == FALSE) {
      print(unique(data$whitelisted_taxa)[i])
      print('look here') #temp stataments to troubleshoot
      #obtain each refseq associated with this genus and store as a dataframe
      temp <- filter(data, whitelisted_taxa == species)
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
            "~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/",protein,"/",
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
        unique(data$refseq_protein_accession)[icount],
        start = 1,
        stop = 12
      ))
      newfasta <- entrez_fetch(
        db = 'protein',
        id = substr(
          unique(data$refseq_protein_accession)[icount],
          start = 1,
          stop = 12
        ),
        rettype = 'fasta'
      )
      write(
        newfasta,
        file = paste(
          "~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/",protein,"/",
          species,
          '.fasta',
          sep = ''
        )
      )
    }
    
  }
}

Dloader('par')

data <- read.csv('gyrB.csv')
