#!/bin/bash

RunLoTo() {

        # This function runs the LoTo method to create and compare similarity between 
        # networks using triplets of graphlets
        # This function is expecting two arguments. In order
        # 1. (passed by value) String indicating path to the extended networks.
        #     Using this network, the reference network and the infered network will be
        #     obtained based on metadata expected to be found inside each file of
        #     the networks
        # 2. (passed by value) String indicating path to place the LoTo results

        printf "${GREEN_COLOR}  Running LoTo implementation${RESET_COLOR}\n\n"

        local INPUT_PATH=$(echo $1 | sed 's/\/*$/\//g')
        local OUTPUT_PATH=$(echo $2 | sed 's/\/*$/\//g')
        local SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
        local BASENAME INFERRED_NETWORK OUTPUT_NAME

        # Creating output folder if it doesn't exist
        [ ! -d $OUTPUT_PATH ] && mkdir "${OUTPUT_PATH}"
        [ ! -d "${OUTPUT_PATH}input_files" ] && mkdir "${OUTPUT_PATH}input_files"

        # Creating expected format to be read by LoTo and running it

        for INFERRED_NETWORK in "${INPUT_PATH}"*; do

                if [ -f $INFERRED_NETWORK ]; then

                        BASENAME=$(basename $INFERRED_NETWORK)
                        OUTPUT_NAME=$(printf "${BASENAME}_LoTo")

                        grep "Known" $INFERRED_NETWORK | cut -f2,3 |
                         perl -ne 'chomp($_); print($_, "\t", 1, "\n")' > "${OUTPUT_PATH}input_files/${BASENAME}_ref"

                        perl -ne 'chomp($_); print($_, "\t", 1, "\n")' <(grep "New" $INFERRED_NETWORK | grep -v "NOT_REFR_ORG" | cut -f2,3) | 
                         cat - "${OUTPUT_PATH}input_files/${BASENAME}_ref" > "${OUTPUT_PATH}input_files/${BASENAME}_inf"

                        perl $SCRIPT_DIR/04_07_01_DoLoto.pl -g "${OUTPUT_PATH}input_files/${BASENAME}_ref" -i "${OUTPUT_PATH}input_files/${BASENAME}_inf" -f "${OUTPUT_PATH}" -o "${OUTPUT_NAME}"
                
                fi

        done

}

# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        RunLoTo $1 $2

fi
