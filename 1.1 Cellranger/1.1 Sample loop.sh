#!/bin/bash
#SBATCH --job-name=
#SBATCH -
#SBATCH --mem=60G
#SBATCH --output=salida_%j.txt
#SBATCH --error=error_%j.txt

module load [Program]/[ver] #carga el programa que se va a utilizar 
for f in [Data folder]/([Sample_dir]/*) /[ID]* #Busca todas las carpetas nombradas con un ID (SRR*) guardando la dirección en el objeto "f"
do
	for i in $f/Sample. [extension] #Dentro de "f" busca guardando el directorio del archivo en el objeto "i"
		do
		j=`basename $i fastq.gz` #Save the name of the file in object
		[Program] [Parameters] --Working Dir $f --Read_files $i –Outputdir [Out Dir]/$j  --Outputname $j #los objetos "f", "i" y "j" reemplazan distintos requerimientos
		echo III $i #Check Object i
		echo JJJ $j #Check Object j      
		done
	done
