#!/bin/bash

#######################

# Este script está diseñado para obtener una tabla de equivalencias entre los distintos identificadores
# usados en las redes de regulación transcripcional y el ID único empleado por NCBI para
# E.coli, B.subtilis, P.aeruginosa, S.enterica, M.tuberculosis y S.aureus

#######################

for file in ../../genomes/*faa; do

	echo "Procesing ${file}"

	BaseName=$(basename $file)
	BaseName=$(echo "${BaseName%.*}.txt")
	File_name=$(printf "${BaseName}" )

	echo "    Saving into ${File_name}"
	grep ">" $file | cut -d"|" -f1,2 |  sed -r 's/>//g' | sed -r 's/\|/\t/g' > $File_name

done
