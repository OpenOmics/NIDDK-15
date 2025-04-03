#! /bin/bash
#SBATCH --partition=norm
#SBATCH --job-name=libra
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=12g
#SBATCH --time=120:00:00
#SBATCH --mail-type=END
#SBATCH --output=libra_script_RNA_celltype.b.out

module load R/4.3.0

Rscript libra_script_RNA_celltype.R > libra_script_RNA_celltype.log

exit


