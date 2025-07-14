#!/bin/sh
#SBATCH -p blades7
#SBATCH --ntasks=1
#SBATCH --mem=128gb
#SBATCH --time=100:00:00
#SBATCH -J harmony_integrate
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL,BEGIN,END

#SBATCH -o /suffolk/WorkGenomics/cdr42/MF_selection/logs/slurm-%j.out
#SBATCH -e /suffolk/WorkGenomics/cdr42/MF_selection/logs/slurm-%j.err

module load conda
conda init
source activate /home/cdr42/.conda/envs/seurat
conda list

Rscript /suffolk/WorkGenomics/cdr42/MF_selection/scripts/seurat_manipulation.R
