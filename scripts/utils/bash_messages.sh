#!bin/bash

ShowHelp() {

	# This function shows on screen the help guide for the main script
	
	local help="
$(basename "$0") [-h|--help] [-g|--genomes FILE] [-n|--nets FILE] [-l|--labels FILE] -- Script to predict transciptional networks based on orthologs

NOTE: It is expected that the content of the files supplied matches in order. For example, the first line in the 
      genomes file describing an organism X must match with the first line in the network file describing also
      the organism X

Arguments:

        -h|--help Show this message.
        -g|--genomes File containing both the path and name of the fasta files to be read. One organism per line.
        -n|--nets File containing both the path and name of the networks in TSV format to be read. One organism per line.
		-l|--labels File containing labels to identify the respective organism rather than using the file name in -g

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
