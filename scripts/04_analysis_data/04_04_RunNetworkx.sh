#!/bin/bash

set -e

GetHubs() {

    printf "${GREEN_COLOR}  Getting Hubs by networkx${RESET_COLOR}\n\n"

    local NETWORK_PATH=$(echo $1 | sed 's/\/*$/\//g')
    local OUTPUT_PATH=$(echo $2 | sed 's/\/*$/\//g')
    local SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

    [ ! -d $OUTPUT_PATH ] && mkdir $OUTPUT_PATH 

    local NETWORK_FILES
    readarray -t NETWORK_FILES < <(ls -dp "${NETWORK_PATH}"* | grep -vP "/$")

    local NET
    local BASE_NAME
    local OUTPUT_FILE_NAME

    for NET in ${NETWORK_FILES[@]}; do
        
        BASE_NAME=$(basename $NET)
        OUTPUT_FILE_NAME=$(printf "${OUTPUT_PATH}${BASE_NAME}_hits")

        python $SCRIPT_DIR/04_04_01_GetHubs.py --input_network $NET --arguments_read_net '{"header":None, "sep":"\t", "usecols":[1,2]}' --output_path $OUTPUT_FILE_NAME
    
    done
}


# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        GetHubs

fi