#! /bin/bash
#SBATCH --job-name=NIDDK15_cellchat
#SBATCH --mail-user=chenv3@nih.gov
#SBATCH --mail-type=END,FAIL
#SBATCH --time=72:00:00
#SBATCH --partition=norm
#SBATCH --mem=240g
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=8
#SBATCH -D /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/cellchat_BD/logs

module purge
export R_LIBS_USER="/data/NHLBI_IDSS/rawdata/NIDDK-15_temp/R4.3.0_cellchat"
module load R/4.3.0

cd /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/cellchat_BD

echo " Running R script	"

Rscript 01_cellChat_splitBcells_trimester2.R

chmod -R 775 *

echo "R	script runs successfully"

