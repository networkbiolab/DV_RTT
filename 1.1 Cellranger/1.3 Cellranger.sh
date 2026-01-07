#!/bin/bash
#SBATCH --job-name=CellrangerLoop
#SBATCH -N 1
#SBATCH --mem=60G
#SBATCH --output=salida_%j.txt
#SBATCH --error=error_%j.txt

module load cellranger #Load cellranger
for f in [homedir]/raw/*/SRR* #Search for folder Samples and save them in object f
do
	for i in $f/*fastq #Save the Fastq in object i
		do
		j=`basename $i fastq.gz` #Save the name of the file in object
		echo cellranger count --id=$j --transcriptome [GRCh 38-3.0.0-2020-A] --fastqs $f --jobmode slurm #Process each object effectively looping and working in all the Samples.
		echo JJJ $j #Check object j
		echo III $i #Check Object i
		echo $f/$j\_2_trimmed.fq.gz #Check Object j      
		done
	done
