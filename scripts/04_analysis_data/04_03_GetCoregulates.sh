#!/bin/bash

set -e

RunCoreg() {

        # All the following output files are going to be saved in the current path and within of a tmp folder
        printf "${GREEN_COLOR}  Getting co-regulators${RESET_COLOR}\n\n"

        local -n GENOMES=$1
        local INPUT_PATH=$(echo $2 | sed 's/\/*$/\//g')
        local OUTPUT=$(echo $3 | sed 's/\/*$/\//g')
        local SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

        declare -a local BASE_NAME
        for ORG in ${GENOMES[@]}; do
            BASE_NAME+=($(basename $ORG))
        done

        # Creating folder if it doesn't exist
        [ ! -d "${OUTPUT}" ] && mkdir "${OUTPUT}" "${OUTPUT}tmp/"

        # Using CoReg to find co-regulators in the networks
        for FILE_NAME in ${BASE_NAME[@]}; do
                Rscript --vanilla $SCRIPT_DIR/04_03_01_CoReg.R -i "${INPUT_PATH}${FILE_NAME}"* -o $OUTPUT
        done

}

# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        RunCoreg

fi