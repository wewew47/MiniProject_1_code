#!/bin/sh
#SBATCH --job-name=ProFit 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3700

#Clear modules, load needed modules
module purge

#Parent for loop to select each wildtype and change to relevant mutateX directory
#Sub for loop to change to each cooridnates directory
#Sub sub for loop to run profit on one of the 5 structures of each mutation at that coordinate
#Run ProFit as with array job formatting
cd /gpfs/home/wms/lfrkmq/MiniProject_1/MutateXInputFinal/mutations
filecount1=$(find *_* -maxdepth 0 -type d|wc -l)
PDB=.pdb
for i in $(eval echo {1..$filecount1})
do
   cd /gpfs/home/wms/lfrkmq/MiniProject_1/Protein/structures/Susceptible/gyrA/FoldXStructures
   wildtype=$(ls | sed -n ${i}p)
   cd /gpfs/home/wms/lfrkmq/MiniProject_1/MutateXInputFinal/mutations
   folder=$(ls | sed -n ${i}p)
   echo folder
   echo $folder
   cd $folder
   filecount2=$(find * -maxdepth 0 -type d|wc -l)
   for j in $(eval echo {1..$filecount2})
   do
      subfolder=$(ls | sed -n ${j}p)
      echo subfolder
      echo $subfolder
      cd $subfolder
      for k in {1..20}
      do
      resistant=$(ls *_0.pdb | sed -n ${k}p)
      cd /gpfs/home/wms/lfrkmq/MiniProject_1/ProFitOutput/gyrA/RandomAll
      profit -f /gpfs/home/wms/lfrkmq/software_tools/ProFit_V3.3/ProFitCommand.txt /gpfs/home/wms/lfrkmq/MiniProject_1/Protein_structures/Susceptible/gyrA/FoldXStructures/${wildtype} /gpfs/home/wms/lfrkmq/MiniProject_1/MutateXInputFinal/mutations/${folder}/${subfolder}/${resistant} > out_$folder$wildtype$resistant.txt
      cd ../${subfolder}
      done
   done
done
