#!bin/bash

ShowHelp() {

	# This function shows on screen the help guide for the main script
	
	local help="

-- Script to predict transciptional networks based on orthologs ($(basename "$0")) 

Arguments Type: 

		[-h|--help] [-g|--genomes STR] [-n|--nets STR] [-l|--labels STR]
		[--proteinortho_output STR] [--batches INT]

NOTE: It is expected that the content of the files supplied matches in order. For example, the first line in the 
      genomes file describing an organism X must match with the first line in the network file describing also
      the organism X

Arguments:

        -h|--help -> Show this message.
        -g|--genomes -> File name containing both the path and name of the fasta files to be read. One organism per line.
        -n|--nets -> File name containing both the path and name of the networks in TSV format to be read. One organism per line.
	-l|--labels -> File name containing labels to identify the respective organism rather than using the file name in -g.
	--proteinortho_output -> Path Name to create and place proteinortho final outputs
	--batches -> Number of baches for running paralized process [Proteinortho]

"
	printf "${help}"
	exit 1

}


ShowArguments() {

	# This function shows on screen the arguments suppplied by the user in format "Name_of_argument: value_of_argument".
	# This functions is expecting a hash table

	local -n ARG=$1
	local KEY

	printf "\n  ${GREEN_COLOR}Arguments${RESET_COLOR}\n\n"

	for KEY in "${!ARG[@]}"; do
		echo "${KEY}: ${ARG[$KEY]}"
	done

	printf "\n"

}

# Color codes variables (global)
GREEN_COLOR=$'\e[1;32m'
RESET_COLOR=$'\033[0m'
