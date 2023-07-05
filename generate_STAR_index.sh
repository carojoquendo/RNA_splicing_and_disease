#!/bin/bash
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH --time=01:00:00

/STAR-2.7.9a/bin/Linux_x86_64_static/STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v38.primary_assembly.annotation.gtf
