#!/bin/bash

set -e

RunProteinortho() {

        # This function runs proteinortho doing paired comparations
        # This functions is expecting 4 arguments. In order:
        # 1. String to put output files (passed by value)
        # 2. Array containing the labels identifying each fasta file or organism (passed by reference)
        # 3. Array containing the path and name of the fasta files (passed by reference)
        # 4. Integer specifying if proteinortho will run in batches

        printf "${GREEN_COLOR}  Running Function RunProteinortho${RESET_COLOR}\n\n"

        local i
        local j
        local COUNTER=1

        # Parsing arguments
        local OUTPUT=$(echo $1 | sed 's/\/*$/\//g')
        local -n GENOMES=$2
        local -n FOLDER_LABELS=$3

        # Defining arguments to paralyze
        [ -z $4 ] && local BATCHES=1 || local BATCHES=$4
        local COUNT_BATCHES=0

        echo "    Checking if ${OUTPUT} exits"
        if [ ! -d $OUTPUT ];then
                echo "      Creating ${OUTPUT}"
                mkdir $OUTPUT
        fi

        # Proteinortho was run in pairs. Every comparation has its own folder to save the results of the orthologous
        for ((i=0; i < ${#FOLDER_LABELS[@]}; i++)) ; do

                ## Creating a folder for each organism compared
                local FULL_OUTPUT_PATH=$(printf "${OUTPUT}${FOLDER_LABELS[$i]}")
                
                if [ ! -d $FULL_OUTPUT_PATH ]; then
                        mkdir $FULL_OUTPUT_PATH
                fi

                for ((j=0; j < ${#GENOMES[@]}; j++));do

                        # Do not search orthologues against itself (dah)
                        if [[ $i -ne $j ]]; then

                                ( mkdir "tmp_${COUNTER}" && cd "tmp_${COUNTER}" 
                                local BASE_NAME=$(basename ${GENOMES[$j]})
                                local PROJECT_NAME=$(printf "${BASE_NAME}_${FOLDER_LABELS[$i]}")

                                proteinortho6.pl --clean --verbose=2 -cpus=5 --cov=70 -project=$PROJECT_NAME ../${GENOMES[$i]} ../${GENOMES[$j]}
                                
                                # Apparently the proteinortho outputs only can be saved in the directory where it's ran
                                mv ${PROJECT_NAME}* ../$FULL_OUTPUT_PATH && cd ../ && rm -r "tmp_${COUNTER}" ) &

                                ((++COUNTER))
                                ((++COUNT_BATCHES)); [ "${COUNT_BATCHES}" -eq "${BATCHES}" ] && COUNT_BATCHES=0 && wait

                        fi

                done

        done

        wait

}

RunOperonMapper() {

        exit 1

}

# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        RunProteinortho
        RunOperonMapper

fi
