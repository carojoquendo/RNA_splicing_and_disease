#!/bin/bash
#SBATCH --ntasks-per-node=20
#SBATCH --time=15:00:00
#SBATCH --mem=40GB

#Load variable folder
. variables

## Create output directory structure
## CHECK PATHS ARE CORRECT !!
MAIN_out="/RNA_samples/$sampleID/MAIN_STAR"
fastq_dir="/RNA_samples/FASTQ/batch2"
sample_dir="/RNA_samples/$sampleID"

mkdir $MAIN_out

#Check fastq quality
module load biobuilds/2017.05
fastqc -o $sample_dir $fastq_dir/$lane1_pair1 $fastq_dir/$lane1_pair2 --threads 32

################################################################################
#Align with STAR
#quantMode GeneCounts > counts number of reads while mapping
#--outReadsUnmapped Fastx > outputs in separate file (fasta) unmapped reads
#--outFilterType BySJout > reduces the number of "spurious" junctions
#################################################################################

cd $MAIN_out
/STAR-2.7.9a/bin/Linux_x86_64_static/STAR \
        --genomeDir /path/to/RNA_refs \
        --readFilesCommand zcat \
        --readFilesIn $fastq_dir/$lane1_pair1 $fastq_dir/$lane1_pair2 \
        --runThreadN 16 \
        --twopassMode Basic \
        --twopass1readsN -1 \
        --outSAMmapqUnique 60 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSJfilterOverhangMin 12 12 12 12 \
        --outSJfilterCountUniqueMin 1 1 1 1 \
        --outSJfilterCountTotalMin 1 1 1 1 \
        --outSJfilterDistToOtherSJmin 0 0 0 0 \
        --quantMode GeneCounts \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM Unsorted

## Clean up STAR intermediates
rm -r _STARgenome
rm -r _STARpass1

####################################################
## Samtools (v1.9) sort and index the alignments
####################################################

module load samtools/1.9

## Sort main alignment
samtools sort -@ 16 Aligned.out.bam > "$sampleID".main.sorted.bam
samtools index "$sampleID".main.sorted.bam

## Clean up unsorted bams as they take up a lot of space
rm Aligned.out.bam

#####################################
#RSeQC QC on main alignment
#####################################

module load python/3.7.3
module load R/4.1.1

## CHECK PATH IS CORRECT !!
qc_outdir="/RNA_samples/$sampleID/RSeQC/"
sorted_main="$MAIN_out/"$sampleID".main.sorted.bam"

mkdir $qc_outdir
cd $qc_outdir

## CHECK PATH IS CORRECT !!
housekeeping="/RNA_refs/hg38.HouseKeepingGenes.bed"
bedfile="/RNA_refs/hg38_Gencode_V28.bed"

python /home/.local/bin/infer_experiment.py -r $bedfile -i $sorted_main > "$sampleID"_infer_exp

python /home/.local/bin/bam_stat.py -i $sorted_main > "$sampleID"_bamstat

python /home/.local/bin/geneBody_coverage.py -r $housekeeping -i $sorted_main -o "$sampleID"_gene_body

python /home/.local/bin/junction_annotation.py -r $bedfile -i $sorted_main -o "$sampleID"_junc_ann

python /home/.local/bin/junction_saturation.py -r $bedfile -i $sorted_main -o "$sampleID"_junc_sat
