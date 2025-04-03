#! /bin/bash
#SBATCH --partition=largemem
#SBATCH --job-name=libra
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=350g
#SBATCH --time=24:00:00
#SBATCH --mail-type=END
#SBATCH --output=libra_script_RNA.b.out

module load R/4.3.0

Rscript libra_script_RNA.R > libra_script_RNA.log

exit


