# all_seqs <- list()
# for(i in 1:16){
#   setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/gyrA")
#   Sus_seq <- unlist(read.fasta(list.files()[i], seqtype = 'AA', as.string = TRUE, seqonly = TRUE))
#   species <- list.files()[i]
#   species <- substr(species, start = 5, stop = nchar(species))
#   setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Resistant/gyrA")
#   resseqs <- c()
#   if(i == 5){
#     species = 'Salmonella1'
#   }
#   if(i == 10){
#     species = 'Campylobacter1'
#   }
#   if(i == 11){
#     species = 'Campylobacter2'
#   }
#   if(i == 15){
#     species = 'Escherichia1'
#   }
#   print(i)
#   print(species)
#   for(j in 1:length(list.files())){
#     if(str_detect(list.files()[j], species) == TRUE){
#       resseqs <- c(resseqs, substr(unlist(read.fasta(list.files()[j], seqtype = 'AA', as.string=TRUE, seqonly = TRUE)), start = 60, stop = 140))
#     }
#     print(class(resseqs))
#     combine_seq <- c(substr(Sus_seq, start = 60, stop = 140), resseqs)
#   }
#   all_seqs[[i]] <- unlist(combine_seq)
# }
# 
# `names<-`(all_seqs, c('Acinetobacter_baumannii', 'Neisseria_gonorrhoeae', 'Neisseria_meningitidis','Pseudomonas_aeruginosa',
#                       'Salmonella', 'Staphylococcus_aureus', 'Staphylococcus_pseudintermedius', 'Burkholderia_cepacia', 
#                       'Burkholderia_pseudomallei', 'Campylobacter1', 'Campylobacter2', 'Clostridioides_difficile', 'Enterococcus_faecalis',
#                       'Enterococcus_faecium', 'Escherichia', 'Klebsiella_pneumoniae'))

library(stringr)
setwd('~/Documents/PhD/Year1/MiniProject_1/Data/Databases')
results <- read.csv('gyrAdatabase.csv')[,-1]

results$mutation <- ''
for(i in 1:length(results$mutation)){
  string <- substr(results$allele[i], start=6, stop=nchar(results$allele[i])) 
  results$mutation[i] <- string
}

groups = data.frame(
  letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
  group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
  stringsAsFactors = F
)

cons_count <- 0
rad_count <- 0
adj_count <- 0
for(i in 1:length(results$mutation)){
  original <- substr(results$mutation[i], 1, 1)
  replacement <- substr(results$mutation[i], nchar(results$mutation[i]), nchar(results$mutation[i]))
  for(j in 1:length(groups$group)){
    if(groups$letter[j]==original){
      org_group <- groups$group[j]
    }
    if(groups$letter[j]==replacement){
      rpl_group <- groups$group[j]
    }
  }
  if(org_group == rpl_group){
    cons_count = cons_count+1
  }
  if(org_group != rpl_group){
    rad_count = rad_count+1
    setwd('~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Resistant/gyrA')
    for(k in 1:length(list.files())){
      if(str_detect(substr(list.files()[k],1,3), as.character(i)) == TRUE){
        sequence <- read.fasta(list.files()[k], seqtype = 'AA', seqonly = TRUE)
        coord <- substr(results$mutation[i], 2,4)
        if(str_detect(substr(coord, 3,3), '[:upper:]')==TRUE){
          coord <- substr(coord, 1,2)
        }
        coord <- as.numeric(coord)
        before <- substr(sequence,(coord-1), (coord-1))
        after <- substr(sequence, (coord+1), (coord+1))
      }
    }
    for(j in 1:length(groups$group)){
      if(groups$letter[j] == before){
        bef_group <- groups$group[j]
      }
      if(groups$letter[j] == after){
        aft_group <- groups$group[j]
      }
    }
    
    if(rpl_group == bef_group | rpl_group == aft_group){
      adj_count <- adj_count+1
      rad_count <- rad_count-1
    }
  }
}

cons_count #20
rad_count #39
adj_count #82

20/141*100 #14%
39/141*100 #28%
82/141*100 #58%



library(stringr)
setwd('~/Documents/PhD/Year1/MiniProject_1/Data/Databases')
results <- read.csv('pbp.csv')[,-1]

results$mutation <- ''
for(i in 1:length(results$mutation)){
  string <- substr(results$allele[i], start=6, stop=nchar(results$allele[i])) 
  results$mutation[i] <- string
}

results$mutation[2] <- 'Q557E'
results$mutation[3] <- 'P601L'
results$mutation[4] <- 'T266A'

groups = data.frame(
  letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
  group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
  stringsAsFactors = F
)

#Initialise counters
cons_count <- 0
rad_count <- 0
adj_count <- 0

for(i in 1:length(results$mutation)){
  print(i)
  
  #Find original and replacement residues
  original <- substr(results$mutation[i], 1, 1)
  replacement <- substr(results$mutation[i], nchar(results$mutation[i]), nchar(results$mutation[i]))
  
  #For each residue check which property group its in
  for(j in 1:length(groups$group)){
    #Check property group of original
    if(groups$letter[j]==original){
      org_group <- groups$group[j]
    }
    #Check property group of replacement
    if(groups$letter[j]==replacement){
      rpl_group <- groups$group[j]
    }
  }
  
  #If the two residues match add to the conservative count
  if(org_group == rpl_group){
    cons_count = cons_count+1
  }
  #If they do not match, check if radical, or if conservative to an adjacent residue
  if(org_group != rpl_group){
    rad_count = rad_count+1
    setwd('~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Resistant/pbp')
    
    for(k in 1:length(list.files())){
      if(str_detect(substr(list.files()[k],1,3), as.character(i)) == TRUE){
        sequence <- read.fasta(list.files()[k], seqtype = 'AA', seqonly = TRUE)
        print(sequence)
        coord <- substr(results$mutation[i], 2,4)
        if(str_detect(substr(coord, 3,3), '[:upper:]')==TRUE){
          coord <- substr(coord, 1,2)
        }
        coord <- as.numeric(coord)
        before <- substr(sequence,(coord-1), (coord-1))
        after <- substr(sequence, (coord+1), (coord+1))
      }
    }
    for(j in 1:length(groups$group)){
      if(groups$letter[j] == before){
        bef_group <- groups$group[j]
      }
      if(groups$letter[j] == after){
        aft_group <- groups$group[j]
      }
    }
    
    if(rpl_group == bef_group | rpl_group == aft_group){
      adj_count <- adj_count+1
      rad_count <- rad_count-1
    }
  }
}

cons_count #2
rad_count #5
adj_count #3

2/10*100 #
5/10*100 #
3/10*100 #
