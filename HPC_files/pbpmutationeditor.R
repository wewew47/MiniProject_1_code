setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Databases")

targets <- read.csv('pbp.csv', header = TRUE)[,-1]
#add column for mutation identifier only (the allele col entries take the form 'genename_mutationidentifier')
targets <- add_column(targets, mutation = '', .after = 'allele')
print(length(targets$allele))
for(i in 1:length(targets$allele)){
  name <- targets$allele[i]
  newname <- ''
  newname = substr(name, start = 6, stop = nchar(name))
  targets$mutation[i] = newname 
}

targets$mutation[2] = 'Q557E'
targets$mutation[3] = 'P601L'
targets$mutation[4] = 'T266A'

#Using the entry in the mutation column, copy that species' wildtype sequence and introduce the point mutation. Saving as new file
for(i in 1:length(targets$mutation)){
  species <- targets$whitelisted_taxa[i] #Identifier for future file organisation if i choose to organise files by species
  sequence <- read.fasta(paste("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/pbp/", targets$whitelisted_taxa[i], '.fasta',sep = ''),
                         seqtype = 'AA', as.string = TRUE, seqonly = TRUE)
  
  #identify the original residue so i can check this base is present at the coordinate in the sequence
  original <- substr(targets$mutation[i],start=1,stop=1)
  #find coordinate in the sequence where the PM occurs
  coordinate <- as.numeric(substr(targets$mutation[i], start = 2, stop = (nchar(targets$mutation[i])-1)))
  print(paste('coordinate is', coordinate))
  print(paste('sequence at this coord is',substr(sequence,start=coordinate,stop=coordinate)))
  #Check for original base not matching the base at the coordinate. Stop if this happens
  if(substr(sequence,start=coordinate,stop=coordinate) != original){   ###IF there are two identical adjacent bases this will mutate the wrong one!!!
    print(paste('Error position is', i))
    print('Error, sequence original position doesnt match the mutation original position')
    print('Increasing coordinate by 1')
    coordinate = coordinate - 1
    if(substr(sequence,start=coordinate,stop=coordinate) != original){
      print(paste('Error position is', i))
      print('Error, sequence original position doesnt match the new mutation position')
    }
    else{
      #Obtain the residue that will replace original and make the mutant sequence
      pointmutant <- substr(targets$mutation[i],start = (nchar(targets$mutation[i])), stop = (nchar(targets$mutation[i])))
      mutantsequence <- paste(substr(sequence,start=1,stop=(coordinate-1)), pointmutant,substr(sequence,start=(coordinate+1),stop = nchar(sequence)), sep = '')
      write.fasta(sequences = mutantsequence, names = 'pbp', file.out = paste("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Resistant/pbp/", targets$mutation[i],'_', species, '.fasta', sep = ''))
    }
  }
}

