#! /bin/bash
#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=200g
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=seurat_merge_bcells_subset.b.out

module load R/4.2.0

Rscript /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/scripts/seurat_merge_bcells_subset.R > seurat_merge_bcells_subset.log

exit
