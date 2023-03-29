#!/bin/sh
#SBATCH --job-name=FoldX_RepairPDB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3700

#Clear modules, load needed modules
module purge
cd /gpfs/home/wms/lfrkmq/MiniProject_1/Protein_structures/Susceptible/gyrA/FoldX_Repaired/
folder=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)
foldx --command=BuildModel --pdb=/gpfs/home/wms/lfrkmq/MiniProject_1/Protein_structures/Susceptible/gyrA/FoldX_Repaired/${folder} --mutant-file=/gpfs/home/wms/lfrkmq/MiniProject_1/AMRFinderDatabase/mutant_files/mutant_file_${folder}.txt --output-file=${folder} --output-dir=BuildModelOutput

