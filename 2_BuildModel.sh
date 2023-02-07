#!/bin/sh
#SBATCH --job-name=FoldX_RepairPDB
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3700

#Clear modules, load needed modules
module purge
cd /gpfs/home/wms/lfrkmq/MiniProject_1/Protein_structures/Susceptible/gyrA/
folder=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)
cd /gpfs/home/wms/lfrkmq/MiniProject_1/Protein_structures/Susceptible/gyrA/FoldX_Repaired/
FoldX --command=BuildModel --pdb=${folder}.pdb --mutant-file=/gpfs/home/wms/lfrkmq/MiniProject_1/Mutation_lists/gyrA/mutant_file_${folder}.txt --output-file=${folder} --output-dir=/gpfs/home/wms/lfrkmq/MiniProject_1/FoldX_Output/

