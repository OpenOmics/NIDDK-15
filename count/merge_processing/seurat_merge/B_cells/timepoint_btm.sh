#! /bin/bash
#SBATCH --partition=norm
#SBATCH --job-name=merge
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=150g
#SBATCH --gres=lscratch:30
#SBATCH --time=12:00:00
#SBATCH --mail-type=END
#SBATCH --output=run_enrichment.b.out

module load R/4.3.0


#for FILE in */FindAllMarkers_MAST_res*;
#do Rscript /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/scripts/20231019_enrich_seurat_merge_findallmarkers_thresholds_pval_methods-BTM_hallmark_kegg.R --pval 0.05 --logfold 0.5 --pct 0.2 --degfile $FILE --outputdir $(dirname $FILE)/ClusterProfiler
#done

Rscript /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/scripts/20231019_enrich_seurat_merge_findallmarkers_thresholds_pval_methods-BTM_hallmark_kegg.R --pval 0.05 --logfold 0.5 --pct 0.2 --degfile 0_Pre/FindAllMarkers_MAST_res0.2.csv --outputdir 0_Pre/ClusterProfiler
Rscript /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/scripts/20231019_enrich_seurat_merge_findallmarkers_thresholds_pval_methods-BTM_hallmark_kegg.R --pval 0.05 --logfold 0.5 --pct 0.2 --degfile 1_1Tri/FindAllMarkers_MAST_res0.2.csv --outputdir 1_1Tri/ClusterProfiler
Rscript /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/scripts/20231019_enrich_seurat_merge_findallmarkers_thresholds_pval_methods-BTM_hallmark_kegg.R --pval 0.05 --logfold 0.5 --pct 0.2 --degfile 2_2Tri/FindAllMarkers_MAST_res0.2.csv --outputdir 2_2Tri/ClusterProfiler
Rscript /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/scripts/20231019_enrich_seurat_merge_findallmarkers_thresholds_pval_methods-BTM_hallmark_kegg.R --pval 0.05 --logfold 0.5 --pct 0.2 --degfile 3_3Tri/FindAllMarkers_MAST_res0.4.csv --outputdir 3_3Tri/ClusterProfiler
Rscript /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/scripts/20231019_enrich_seurat_merge_findallmarkers_thresholds_pval_methods-BTM_hallmark_kegg.R --pval 0.05 --logfold 0.5 --pct 0.2 --degfile 4_Post2M/FindAllMarkers_MAST_res0.4.csv --outputdir 4_Post2M/ClusterProfiler
Rscript /data/NHLBI_IDSS/rawdata/NIDDK-15_temp/count/scripts/20231019_enrich_seurat_merge_findallmarkers_thresholds_pval_methods-BTM_hallmark_kegg.R --pval 0.05 --logfold 0.5 --pct 0.2 --degfile 5_Post6M/FindAllMarkers_MAST_res0.2.csv --outputdir 5_Post6M/ClusterProfiler

exit
