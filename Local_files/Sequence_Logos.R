library(seqinr)
library(memes)
library(stringr)

all_seqs <- list()
for(i in 1:16){
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/gyrA")
  Sus_seq <- unlist(read.fasta(list.files()[i], seqtype = 'AA', as.string = TRUE, seqonly = TRUE))
  species <- list.files()[i]
  species <- substr(species, start = 5, stop = nchar(species))
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Resistant/gyrA")
  resseqs <- c()
  if(i == 5){
    species = 'Salmonella1'
  }
  if(i == 10){
    species = 'Campylobacter1'
  }
  if(i == 11){
    species = 'Campylobacter2'
  }
  if(i == 15){
    species = 'Escherichia1'
  }
  print(i)
  print(species)
  for(j in 1:length(list.files())){
    if(str_detect(list.files()[j], species) == TRUE){
      resseqs <- c(resseqs, substr(unlist(read.fasta(list.files()[j], seqtype = 'AA', as.string=TRUE, seqonly = TRUE)), start = 60, stop = 110))
    }
    print(class(resseqs))
    combine_seq <- c(substr(Sus_seq, start = 60, stop = 110), resseqs)
  }
  all_seqs[[i]] <- unlist(combine_seq)
}
par(mfrow=c(1,2))
plot_sequence_heatmap(all_seqs[[15]], alph = 'AA', legend = 'right', heights = c(0.3,0.7))
plot_sequence_heatmap(all_seqs[[2]], alph = 'AA', legend = 'right', heights = c(0.3,0.7))

all_seqs <- c()
for(i in 1:16){
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/gyrA")
  Sus_seq <- unlist(read.fasta(list.files()[i], seqtype = 'AA', as.string = TRUE, seqonly = TRUE))
  species <- list.files()[i]
  species <- substr(species, start = 5, stop = nchar(species))
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Resistant/gyrA")
  resseqs <- c()
  if(i == 5){
    species = 'Salmonella1'
  }
  if(i == 10){
    species = 'Campylobacter1'
  }
  if(i == 11){
    species = 'Campylobacter2'
  }
  if(i == 15){
    species = 'Escherichia1'
  }
  print(i)
  print(species)
  for(j in 1:length(list.files())){
    if(str_detect(list.files()[j], species) == TRUE){
      resseqs <- c(resseqs, substr(unlist(read.fasta(list.files()[j], seqtype = 'AA', as.string=TRUE, seqonly = TRUE)), start = 61, stop = 110))
    }
    print(class(resseqs))
    combine_seq <- c(substr(Sus_seq, start = 61, stop = 110), resseqs)
  }
  all_seqs <- c(all_seqs, combine_seq)
}

plot_sequence_heatmap(all_seqs, alph = 'AA', legend = 'right', heights = c(0.3,0.7)) 
x

ggplot(x = c(1:10), y = c(1:10))
###PBP Analysis

all_seqs <- list()
for(i in 1:5){
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/pbp")
  Sus_seq <- unlist(read.fasta(list.files()[i], seqtype = 'AA', as.string = TRUE, seqonly = TRUE))
  species <- list.files()[i]
  species <- substr(species, start = 5, stop = nchar(species))
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Resistant/pbp")
  resseqs <- c()
  if(i == 1){
    species = 'Staphylococcus_aureus1'
  }
  if(i == 2){
    species = 'Staphylococcus_aureus2'
  }
  print(i)
  print(species)
  for(j in 1:length(list.files())){
    if(str_detect(list.files()[j], species) == TRUE){
      print(j)
      resseqs <- c(resseqs, substr(unlist(read.fasta(list.files()[j], seqtype = 'AA', as.string=TRUE, seqonly = TRUE)), start = 180, stop = 270))
    }
    combine_seq <- c(substr(Sus_seq, start = 180, stop = 270), resseqs)
  }
  all_seqs[[i]] <- unlist(combine_seq)
}
par(mfrow=c(1,2))

plot_sequence_heatmap(all_seqs[[2]], alph = 'AA', legend = 'right', heights = c(0.3,0.7))


all_seqs <- c()
for(i in 1:5){
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Susceptible/pbp")
  Sus_seq <- unlist(read.fasta(list.files()[i], seqtype = 'AA', as.string = TRUE, seqonly = TRUE))
  species <- list.files()[i]
  species <- substr(species, start = 5, stop = nchar(species))
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_fastas/Resistant/pbp")
  resseqs <- c()
  if(i == 1){
    species = 'Staphylococcus_aureus1'
  }
  if(i == 2){
    species = 'Staphylococcus_aureus2'
  }
  print(i)
  print(species)
  for(j in 1:length(list.files())){
    if(str_detect(list.files()[j], species) == TRUE){
      resseqs <- c(resseqs, substr(unlist(read.fasta(list.files()[j], seqtype = 'AA', as.string=TRUE, seqonly = TRUE)), start = 180, stop = 610))
    }
    print(class(resseqs))
    combine_seq <- c(substr(Sus_seq, start = 180, stop = 610), resseqs)
  }
  all_seqs <- c(all_seqs, combine_seq)
}

plot_sequence_heatmap(all_seqs, alph = 'AA', legend = 'right', heights = c(0.3,0.7))

for(i in 1:length(all_seqs)){
  print(nchar(all_seqs[i]))
}





ggplot(mtcars, aes(x=mpg, y=disp, size=hp, col=as.factor(cyl), shape=as.factor(gear))) +
  geom_point() +
  labs(x='test1', y='test2', size='hp',
       col= 'test3', shape='test4')
