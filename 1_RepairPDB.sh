#!/bin/sh
#SBATCH --job-name=FoldX_RepairPDB
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=3700

#Clear modules, load needed modules
module purge
cd /gpfs/home/wms/lfrkmq/MiniProject_1/Protein_structures/Susceptible/gyrA/
folder=$(ls * | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo ${folder}
cd /gpfs/home/wms/lfrkmq/MiniProject_1/Protein_structures/Susceptible/gyrA/${folder}/
FoldX --command=RepairPDB --pdb=0_ranked.pdb --output-file=${folder}
