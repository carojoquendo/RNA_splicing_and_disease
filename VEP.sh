#!/bin/bash
#SBATCH --mem=15G
#SBATCH --time=48:00:00
#SBATCH --nodes=1

module load ensembl-vep/103.0

reference=/path/to/references
Splicefiles=/path/to/SpliceAIfiles

vep \
-i Final.vcf.gz \
-o VEPresults.txt \
--vcf \
--everything \
--dir_cache /mainfs/lyceum/wir1n21/.vep \
--plugin SpliceAI,snv=$Splicefiles/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$Splicefiles/spliceai_scores.raw.indel.hg38.vcf.gz \
--offline

