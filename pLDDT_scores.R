setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_structures/gyrA")
files <- list.files()
plddt_vec <- c()
for(i in 1:length(files)){
  setwd(files[i])
  plddt <- mean(read.pdb('ranked_0.pdb')[[1]]$b)
  plddt_vec <- c(plddt_vec, plddt)
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_structures/gyrA")
}

mean(plddt_vec)

setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_structures/pbp")
files <- list.files()
plddt_vec <- c()
for(i in 1:length(files)){
  setwd(files[i])
  plddt <- mean(read.pdb('ranked_0.pdb')[[1]]$b)
  plddt_vec <- c(plddt_vec, plddt)
  setwd("~/Documents/PhD/Year1/MiniProject_1/Data/Protein_structures/pbp")
}

mean(plddt_vec)


setwd('A_baumannii')
test <- read.pdb('ranked_0.pdb')[[1]]$b
  
  
#61 to 66