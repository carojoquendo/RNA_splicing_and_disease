#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=60GB

module load conda

source activate SALMON

# load variables folder
# this folder has variables for sampleID and lane pairs for fastq files
. variables

transcript="/RNA_refs/salmon_index/transcript_index"
fastq_dir="/splicing_disease/FASTQ"

salmon quant -p 16 -i $transcript --gcBias -l A -1 $fastq_dir/$lane1_pair1 $fastq_dir/$lane2_pair1 $fastq_dir/$lane3_pair1 $fastq_dir/$lane4_pair1 -2 $fastq_dir/$lane1_pair2 $fastq_dir/$lane2_pair2 $fastq_dir/$lane3_pair2 $fastq_dir/$lane4_pair2 -o "$sampleID"_quant


conda deactivate
