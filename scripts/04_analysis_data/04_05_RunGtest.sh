#!/bin/bash

set -e

RunGtest() {

    printf "${GREEN_COLOR}  Running G-test${RESET_COLOR}\n"

    local INPUT_PATH=$(echo $1 | sed 's/\/*$/\//g')
    local OUTPUT_PATH=$(echo $2 | sed 's/\/*$/\//g')
    local SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

    [ ! -d $OUTPUT_PATH ] && mkdir -p $OUTPUT_PATH

    local FILE_PATH_NAME=$(printf "${OUTPUT_PATH}g_test_results.txt")

    python $SCRIPT_DIR/04_05_01_RunGtest.py --input_path $INPUT_PATH --filter_files "^GCF" --output_file $FILE_PATH_NAME --arguments_load_files '{"sep":"\t", "header":None, "usecols":[1,2,3,4,8], "index_col":None, "names":["TF","TG","STATUS","ORG_REF", "TU_UNIT"]}'

}

# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        RunGtest

fi