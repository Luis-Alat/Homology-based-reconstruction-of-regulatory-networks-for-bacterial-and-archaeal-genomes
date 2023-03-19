#!/bin/bash

#################################################################################

# Este script está deseñado para procesar y posteeriormente
# explorar la red de bacillus subtillus descargada de wikisubtillus a comparación
# de la usada por DBTBS

#################################################################################

set -e

GetReferenceNet() {

	# Removiendo elementos innecesarios
	printf "Removing unnecesary elements\n"
	pre_cleaned_net=$(grep -vi "sigma factor" bacsu_regulations-2023-01-13.csv | column -s',' -t |
				grep -v "box" | grep -vi "Riboswitch" |
				grep -v "stringent response")

	# Consiguiendo las equivalencias entre NCBI ID y nombre genérico
	printf "Getting NCBI ID\n"
	equival_pre_ready=$(grep ">" ../../IIMAS_V4/genomes/GCF_000009045.1_ASM904v1_B_subtilis_168_genomic.faa | cut -d"|" -f1,2 |
				sed -r 's/\|/\t/g' | sed -r 's/>//g' |
				cut -f1,3 | awk -F'\t' '{OFS=FS}{print $2,$1}')

	# Formateando equivalencias para hacer múltiples RegEx
	printf "Generating multiples RegEx rules\n"
	equival=$(printf "${equival_pre_ready}" | sed -r 's/^/s\/\\</g' |
				sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' |
				sed -r 's/\./\\./g')

	# Mapeando los ID a la red
	printf "Mapping\n"
	pre_cleaned_net_mapped=$(printf "${equival}" | sed -f - <(printf "${pre_cleaned_net}"))

	# Limpiando la red y obteniendo la versión de referencia
	printf "Retrieving network of reference\n"
	net_ref=$(printf "${pre_cleaned_net_mapped}" | grep -v "regulator locus" |
			sed -r 's/\"\s+\"/\t/g' | awk -F'\t' '{OFS=FS}{print $2,$5,$3,$1,$4}' |
			sed -r 's/\"//g' |
			perl -F"\t" -nae 'if(($F[0] =~ /NP_/) or ($F[0] =~ /YP_/)){
						if(($F[1] =~ /NP_/) or ($F[1] =~ /YP_/)){
							print $_
						}
					}' | sort | uniq | nl - | sed -r 's/^\s+//g')

	printf "${net_ref}" > bacillus_wiki_regulations-2023-01-13.tsv

}

ShowGeneralMetrics() {

	local network=$1

	read -r Edges FileName <<<$( wc -l "${network}" )
	TFs=$( cut -f2 "${network}" | sort | uniq | wc -l)
	TGs=$( cut -f3 "${network}" | sort | uniq | wc -l)

	echo ""
	printf "Network: ${FileName}\n"
	printf "Edges: ${Edges}\n"
	printf "TFs: ${TFs}\n"
	printf "TGs: ${TGs}\n"
	echo ""

}

CompareNets() {
	echo ""
	echo "Comparing networks"
	python ComparaBacillus.py
	echo ""
}

GetReferenceNet
ShowGeneralMetrics "bacillus_wiki_regulations-2023-01-13.tsv"
CompareNets
