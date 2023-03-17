#!/bin/bash

set -e

RunProteinortho() {

    # This function runs proteinortho doing paired comparations between a reference 
    # organism and a target organism
    # This functions is expecting 4 arguments. In order:
    # 1. Array containing the path and name of reference organism of the fasta files (passed by reference)
    # 2. Array containing the labels identifying each fasta file of the reference organisms (passed by reference)
    # 3. Array containing the path and name of target organisms of the fasta files (passed by reference)
    # 4. String to put output files (passed by value)
    # 5. Integer to paralize proteinortho (launch it multiple times)

    local -n GENOMES_REFERENCE=$1
    local -n LABELS=$2
    local -n GENOMES_TARGET=$3
    local OUTPUT=$(echo $4 | sed 's/\/*$/\//g')

    # Defining arguments to paralyze
    [ -z $5 ] && local BATCHES=1 || local BATCHES=$5
    local COUNT_BATCHES=0

	printf "${GREEN_COLOR}  Running Function RunProteinortho${RESET_COLOR}\n\n"

	### Creating a new directory if it doesn't exit
	[ ! -d "${OUTPUT}" ] && mkdir -p "${OUTPUT}"

    local i
    local j
    local COUNTER=1

	for ((i=0; i < ${#GENOMES_REFERENCE[@]}; i++)); do

        [ ! -d "${OUTPUT}${LABELS}" ] && mkdir "${OUTPUT}${LABELS[$i]}"

        for ((j=0; j < ${#GENOMES_TARGET[@]}; j++)); do

            ( mkdir "tmp_${COUNTER}" && cd "tmp_${COUNTER}"

            TARGET_GENOME_BASENAME=$(basename "${GENOMES_TARGET[$j]}")
            PROJECT_NAME=$(echo "${TARGET_GENOME_BASENAME}_${LABELS[$i]}" )

            proteinortho6.pl -verbose=1 -cpus=6 --cov=70 -project=$PROJECT_NAME ../$GENOMES_REFERENCE[$i] ../$GENOMES_TARGET[$j]

            mv ${PROJECT_NAME}* "../${OUTPUT}${LABELS[$i]}" && cd ../ ) &

            ((++COUNTER))
            ((++COUNT_BATCHES)); [ "${COUNT_BATCHES}" -eq "${BATCHES}" ] && COUNT_BATCHES=0 && wait

        done

    done

    wait

    rm -r "tmp_"*

}

# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../../utils/tracking.sh
        source ../../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        RunProteinortho

fi