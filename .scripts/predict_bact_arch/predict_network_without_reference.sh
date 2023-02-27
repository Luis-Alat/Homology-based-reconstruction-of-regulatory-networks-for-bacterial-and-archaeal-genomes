#!bin/bash

####################################################################
#### Script to predict interactions based on our organism models ###
####################################################################


############################ Functions ############################

############################################################################
########### THE FOLLOWING FUNCTION IS DESIGN TO RUN PROTEINORTHO ###########
############################################################################

runProteinortho() {

	#### Running proteinortho individually

	printf "\n\n${pgreen}Running proteinortho${NC}\n\n"

	### Creating a new directory if it doesn't exit
	if [ ! -d "${OutputPath}proteinortho" ]; then
        	mkdir "${OutputPath}proteinortho"
	fi

	for ((i=0;i<${#FolderNames[@]};i++)); do

        	if [ ! -d "${OutputPath}proteinortho/${FolderNames[$i]}" ];then
                	mkdir "${OutputPath}proteinortho/${FolderNames[$i]}"
        	fi

        	# Apparently proteinortho only saves results in the current folder
        	cd "${OutputPath}proteinortho/${FolderNames[$i]}"

	        for ((j=0;j<${#TargetGenomes[@]};j++)); do

        	        ProjectName=$(echo ${TargetGenomes[$j]}"_"${FolderNames[$i]})
			printf "$ProjectName"
                	proteinortho5 -verbose -cpus=5 --cov=70 -project=$ProjectName $GenomesPath${GenomeNames[$i]} $InputPath${TargetGenomes[$j]}

	        done
	done

	cd ${OutputPath}

}


########################################################################################
########### THE FOLLOWING FUNCTION IS DESIGN TO CREATE THE EXTENDED NETWORKS ###########
########################################################################################

makeExtention() {

	printf "\n\n${pgreen}Getting new networks${NC}\n\n"

	# Creating folder if it doesn't exist
	if [ ! -d "${OutputPath}net" ];then
		mkdir "${OutputPath}net"
	fi

	# In this "for" every ortolog is retrieved using the proteinortho output and a new file is created. That new file is used to create another one wich is useful to get the new interactions
        # All outputs are going to be saved inside of "net" directory in the respective organism folder

	# Run for in N batches
	Batches=700

	for ((i=0; i<${#FolderNames[@]}; i++)); do

        	# Creating folder if it doesn't exist
        	if [ ! -d "${OutputPath}net/${FolderNames[$i]}" ];then
                	mkdir "${OutputPath}net/${FolderNames[$i]}" "${OutputPath}net/${FolderNames[$i]}/ortoFiles/" "${OutputPath}net/${FolderNames[$i]}/PreNewInteractionsFiles/"
	        fi

		# Defining output path
		local PathOrtho=$(echo "${OutputPath}net/${FolderNames[$i]}/ortoFiles/")
		local PathNew=$(echo "${OutputPath}net/${FolderNames[$i]}/PreNewInteractionsFiles/")

		# Counting iterations
		count_batches=1

        	# Iterating in every "blast-graph" file found inside of ${FolderNames[@]} using the file names ("${TargetGenomes[@]}")
		for file in ${TargetGenomes[@]}; do

			# Getting important file of proteinortho results
			local ProtResult=$(ls ${OutputPath}proteinortho/${FolderNames[$i]}/${file}*blast-graph)

			# Setting new file names and outputPath
                	local OutFileOrtho=$(echo $file | sed -r 's/$/_ortho/g')
                	local OutFileNewInteractions=$(echo $file | sed -r 's/$/_new_interactions/g')

	                printf "(background) Making files from "$ProtResult"\n"

 	                # Given a specific file inside of the folder, search for a pattern (file name) and get the column number. Later it will be useful to reorder the file in the "correct" order (giving structure)
        	        ## First 5 lines contain the necessary information ## Search by using the name file our column of interest and print the right order separated by comma
                	local string=$(head -n 5 $ProtResult | awk -F$"\t" -v pattern=$file '$1 ~ pattern {print 2","1} $2 ~ pattern {print 1","2}')

			## The output before now will be an array
                	IFS="," read -r -a columnOrder <<< "$string"

		        # Reording the columns
        	        ## Pasting in this order the column described in ${columnOrder[0]} (model organism) and later ${columnOrder[1]} (model target) ## Replacing "|" by TAB ## Saving only the columns importants
	                paste <(grep -v "#" $ProtResult | cut -f${columnOrder[0]}) <(grep -v "#" $ProtResult | cut -f${columnOrder[1]})  | sed -r 's/\|/\t/g' | cut -f1,3 > $PathOrtho$OutFileOrtho

        	        # Given the new file, get a previous file useful to find the new interactions. For this purpose will be searched and replaced (in the modified networks) every match into a $OutFileOrtho format and add the organism where the interaction came from
                	sed -r 's/^/s\/\\</g' $PathOrtho$OutFileOrtho | sed -r 's/\t/\\>\//g' | sed -r 's/$/MATCH\/gI;/g' | sed -r 's/\./\\./g' | sed -f - <(cut -f1-3 $ModelPath${NetFiles[$i]}) | perl -F"\t" -nae 'if(($F[1] =~ /MATCH/) and ($F[2] =~ /MATCH/)){print $_}' | sed -r 's/MATCH//g' | awk -F $'\t' -v organism=${FolderNames[$i]} 'BEGIN{OFS = FS} {print organism,$0}' > $PathNew$OutFileNewInteractions &

			# Wait until end N batches of the previous command
			((++count_batches)); [ "${count_batches}" -eq "${Batches}" ] && count_batches=1 && wait
		done
	done

	wait

	# Creating folder if it doesn't exist
        if [ ! -d "${OutputPath}net/Merge/" ];then
                mkdir "${OutputPath}net/Merge/"
        fi

        # Merging results ("for" loop below) of every Folder (${FolderNames[@]}) using the name found inside of "${TargetGenomes[@]}"
        # Results will be saved in the folder "${OutputPath}net/Merge/"

	for file in ${TargetGenomes[@]}; do

        	# Defining a new file name and making a temporally empty file
        	MergeFileName=$(echo $file"_Merged")
		echo -n > "${OutputPath}net/Merge/tmp"

	        # Add a column specifying the organism where the new interaction came from
	        for folder in ${FolderNames[@]}; do

        	        printf "Creating a tmp file using "$file" found in "$folder"\n"

                	# Merging all results into a temporaly file
	                cat ${OutputPath}net/$folder/PreNewInteractionsFiles/$file* >> "${OutputPath}net/Merge/tmp"

	        done

	        printf "Making file "$MergeFileName"\n"

	        # This line calls an sql implementation to group the columns associated with regulatory interactions to describe (1) transcription factor, (2) target, (3) organism where that interaction was found and (4), line number associated with that interaction found in modified networks. In addition adds a number indicating in how many organism that interactions was found
        	bash "${SCRIPTPATH}/sqlite_join_net_without_reference.sh" "${OutputPath}net/Merge/tmp" | perl -F"\t" -nae 'chomp($_); @array=split(",",$F[2]); %seen=(); @unique = grep { ! $seen{$_} ++ } @array; $len = scalar(@unique); print ($_,"\t",$len,"\n")' > "${OutputPath}net/Merge/${MergeFileName}"

	        rm "${OutputPath}net/Merge/tmp"
	done

}

###########################################################################
########### THE FOLLOWING FUNCTION IS DESIGN FOR TRACKING ERRORS ###########
###########################################################################

set -eE -o functrace

failure() {
        local lineno=$1
        local msg=$2
        echo "Error at $lineno: $msg"
}

trap 'failure ${LINENO} "$BASH_COMMAND"' ERR


#####################################################################################
########### THE FOLLOWING FUNCTION IS DESIGN TO MAKE AN EXTENTION WITH TU ###########
#####################################################################################

makeExtentionTU() {

	printf "\n\n${pgreen}Using transcriptional units files${NC}\n\n"

        # Creating folder if it doesn't exist
        if [ ! -d "${OutputPath}net/Merge_plus_TU" ]; then
                mkdir "${OutputPath}net/Merge_plus_TU" "${OutputPath}net/Merge_plus_TU/tmp"
        fi

	# Count iterations
	count_batches=1
	# Run for in N batches
        Batches=500

	for file in ${TargetGenomes[@]}; do

		# Removing "faa" extention
		local file_cp=$(printf $file | sed -r 's/\.faa$//g')

        	# Verifying if the file exists
	        ( if [ -f $PathTu$file_cp* ]; then

			# Getting input files (full name)
			local InputExt=$(find "${OutputPath}net/Merge/${file_cp}"*)
			local InputTU=$(find "${PathTu}${file_cp}"*)
			local OutputF=$(printf "${OutputPath}net/Merge_plus_TU/tmp/${file}_net_predictions_TU")
			local OutputFinal=$(printf "${OutputPath}net/Merge_plus_TU/${file}_net_predictions_TU")

			# Maping "new" interactions
        	        python3 "${SCRIPTPATH}/ParserTUandExt.py" $InputExt $InputTU $OutputF 2

			# Giving a new format to the earlier output
			tmp_1=$(printf "${file}_tmp_1")
			awk -F$"\t" 'BEGIN{OFS=FS}{print $3,$4,"NOT_REFR_ORG","NOT_REFR_NUM","NULL",$2,$5}' $OutputF > "${OutputPath}net/Merge_plus_TU/tmp/${tmp_1}"

			# Giving a new format to the extended network
			tmp_2=$(printf "${file}_tmp_2")
			awk -F$"\t" 'BEGIN{OFS=FS}{print $0,"NOT_STRAND","NOT_TU_REFER"}' $InputExt > "${OutputPath}net/Merge_plus_TU/tmp/${tmp_2}"

			# Merging results and eliminating information not relevant (labels no longer used)
			tmp_3=$(printf "${file}_tmp_3")
			cat "${OutputPath}net/Merge_plus_TU/tmp/${tmp_2}" "${OutputPath}net/Merge_plus_TU/tmp/${tmp_1}" > "${OutputPath}net/Merge_plus_TU/tmp/${tmp_3}"
			bash "${SCRIPTPATH}/sqlite_group_tu_without_reference.sh" "${OutputPath}net/Merge_plus_TU/tmp/${tmp_3}" | sed -r 's/,NOT_REFR_ORG//g' | sed -r 's/,NOT_REFR_NUM//g' | sed -r 's/,NULL//g' | sed -r 's/NOT_STRAND,//g' | sed -r 's/NOT_TU_REFER,//g' | sed -r 's/,\w+_tu//g' > $OutputFinal

			rm "${OutputPath}net/Merge_plus_TU/tmp/${tmp_1}" "${OutputPath}net/Merge_plus_TU/tmp/${tmp_2}" "${OutputPath}net/Merge_plus_TU/tmp/${tmp_3}"

	        fi ) &

		((++count_batches)); [ "${count_batches}" -eq "${Batches}" ] && count_batches=1 && wait

	done

	wait

}


############################ End functions ############################

### Path of this script

SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname $SCRIPT)

### Defining default variables

### Names of model organism (to create folder, assign sufix, etc), genome and network of the organism that we know (In the same order)
FolderNames=('Bacsu' 'Ecoli' 'Myctu' 'Pseudo' 'Salmon' 'Staph')
#GenomeNames=("GCF_000009045.1_ASM904v1_B_subtilis_168_genomic.faa" "GCF_000005845.2_E_coli_K12_genomic.faa" "GCF_000195955.2_ASM19595v2_M_tuberculosis_H37Rv_genomic.faa" "GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_genomic.faa" "GCF_000006945.2_ASM694v2_S_enterica_LT2_genomic.faa" "GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa")
NetFiles=('GCF_000009045.1_ASM904v1_B_subtilis_168_modified_net.txt' 'GCF_000005845.2_E_coli_K12_modified_net.txt' 'GCF_000195955.2_ASM19595v2_M_tuberculosis_H37Rv_modified_net.txt' 'GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_modified_net.txt' 'GCF_000006945.2_ASM694v2_S_enterica_LT2_modified_net.txt' 'GCF_000009645.1_ASM964v1_S_aureus_N315_modified_net.txt')

### Names of model organism (old name fasta file)

GenomeNames=("B_subtilis_168_GCF_000009045.1_ASM904v1_genomic.faa" "E_coli_RegulonDB_U00096.3.faa" "Mycobacterium_tuberculosis_H37Rv_proteins_v3.faa" "P_aeruginosa_CF_000006765.1_ASM676v1_genomic.faa" "salmonella_enterica_serovar_LT2_RN_GCF_000006945.2_ASM694v2_genomic.faa" "S_aureusN315_GCA_000009645.1_ASM964v1_genomic.faa")


### Defining color messages on screen
pgreen=$'\e[1;32m' # Green color
NC='\033[0m' # No Color

### Help message if -h was provided or if there were no arguments provided

help="\n$(basename "$0") [-h] [-d PATH] [-o PATH] [-m PATH] [-g PATH] [-a STRING] [-b STRING] [-c STRING] [-t PATH] -- Script to predict transciptional networks based on orthologs

-d and -o commands must be provided

where:
        -h show this message
        -d directory where the FASTA files are (\".faa\" filename extension)
	-o output path where all files will be saved
	-m [OPTIONAL] directory where the modified transcriptional networks are (default at the same level as this script ./modifiedNetworks)
	-g [OPTIONAL] directory where the FASTA files of the model organisms are (default at the same level as this script ./genomes)
	-a [OPTIONAL] append new organism models. String with different fasta files names (must be found in -g) of the new organism models
           separated them by comma (e.g \"org_1.faa,org_2.faa\")
	-b Mandatory if -a was supplied. String with the different sufix of -a (useful, for example, to name folders).They must be in the same
           order as -a (e.g \"org_1,org_2\")
	-c Mandatory if -a was supplied. String with the different file names (must be found in -m) of the modified networks of -a. They must be
           in the same order as -a (e.g \"net_org_1.txt,net_org_2.txt\")
	-t [OPTIONAL] Path where the transcriptional units files (BZ2) are (same name as FASTA files without \".faa\" extention)

"

### Check if arguments were no provided and exit

if [ $# -eq 0 ]; then
        printf "$help"
        exit
fi

### Default parameters, getting absolute path and adding "/" at the end of the line

GenomesPath=$(realpath ./genomes | sed -r 's/($|\/$)/\//g')
ModelPath=$(realpath ./modifiedNetworks | sed -r 's/($|\/$)/\//g')

### Changing IFS if -a was supplied
IFS=","

### Retrieving parameters

while getopts ":hd:o:m:g:a:b:c:t:" opt; do
  case $opt in
    h) printf "$help"
       exit
    ;;
    o) OutputPath="$OPTARG"
    ;;
    :) printf "Missing argument for -$OPTARG\n" >&2
       exit 1
    ;;
    d) InputPath="$OPTARG"
    ;;
    :) printf "Missing argument for -$OPTARG\n" >&2
       exit 1
    ;;
    m) ModelPath="$OPTARG"
    ;;
    g) GenomesPath="$OPTARG"
    ;;
    a) AppendFasta="$OPTARG"
       AppendFasta=$(printf "${AppendFasta}" | sed -r 's/\s*,\s*/,/g')
       AppendFasta=($AppendFasta)
       GenomeNames=("${GenomeNames[@]}" "${AppendFasta[@]}")
    ;;
    b) AppendSufix="$OPTARG"
       AppendSufix=$(printf "${AppendSufix}" | sed -r 's/\s*,\s*/,/g')
       AppendSufix=($AppendSufix)
       FolderNames=("${FolderNames[@]}" "${AppendSufix[@]}")

    ;;
    c) AppendNet="$OPTARG"
       AppendNet=$(printf "${AppendNet}" | sed -r 's/\s*,\s*/,/g')
       AppendNet=($AppendNet)
       NetFiles=("${NetFiles[@]}" "${AppendNet[@]}")
    ;;
    t) PathTu="$OPTARG"
    ;;
    \?) printf "\nInvalid option -$OPTARG" >&2
        printf "$help" >&2
        exit 1
    ;;
  esac
done

### Returning default value IFS
IFS=" "

### Verifying sizes of arrays

if [ ${#FolderNames[@]} -ne ${#GenomeNames[@]} ] || [ ${#NetFiles[@]} -ne ${#GenomeNames[@]} ] || [ ${#FolderNames[@]} -ne ${#NetFiles[@]} ]; then

	printf "\nNumber of input files dont match\n\n" >&2
	exit 1
fi

### Verifying if input directory exist and adding "/" if it's necessary (Getting absolute path too)
if [ -d $InputPath ]; then
	InputPath=$(realpath $InputPath | sed -r 's/($|\/$)/\//g')
else
	printf "Input directory does not exist ${InputPath}\n" >&2
	exit 1
fi

### Verifying if output directory exist and adding "/" if it's necessary (Getting absolute path too)
if [ -d $OutputPath ]; then

	OutputPath=$(realpath $OutputPath | sed -r 's/($|\/$)/\//g')
	CheckFolder=$(ls $OutputPath)

	### Verifying if output directory is empty. Printing warning message
	if [ ! -z $CheckFolder ]; then
		printf "$OutputPath is not empty. Files could be overwritten\n"
		secs=10; while [ $secs -gt 0 ]; do echo -ne "Waiting $secs seconds to continue...\033[0K\r"; sleep 1; : $((secs--)); done
	fi

else
	printf "Output directory does not exist ${OutputPath}\n" >&2
	exit 1
fi

### Verifying if transcriptional unit directory exist (Getting absolute path) if this argument was supplied

if [ ! -z $PathTu ]; then

	if [ -d $PathTu ]; then

		PathTu=$(realpath $PathTu | sed -r 's/($|\/$)/\//g')

	else
		printf "Transcriptional units directory does not exist ${PathTu}\n" >&2
		exit 1

	fi

fi

### Parameters

printf "\n${pgreen}Parameters${NC}\n\n"

printf "Output path: ${OutputPath}\n"
printf "Input path: ${InputPath}\n"
printf "Path of the fasta files of the organism models: ${GenomesPath}\n"
printf "Path of the modified networks: ${ModelPath}\n"
printf "\nSufix, fasta and networks files:\n"

for ((i=0;i<${#FolderNames[@]};i++)); do
        printf "\t${FolderNames[$i]}\t${GenomeNames[$i]}\t${NetFiles[$i]}\n"
done

if [ ! -z $PathTu ]; then
	printf "\nPath of transcription units: ${PathTu}\n\n"
fi

### Making an array with the target genomes
shopt -s nullglob
TargetGenomes=($InputPath/*faa)
shopt -u nullglob

### Getting only base name

for ((i=0;i<${#TargetGenomes[@]};i++)); do
	TargetGenomes[$i]=$(basename "${TargetGenomes[$i]}")
done

### Running functions

#runProteinortho
makeExtention

if [ ! -z $PathTu ]; then
	makeExtentionTU
fi
