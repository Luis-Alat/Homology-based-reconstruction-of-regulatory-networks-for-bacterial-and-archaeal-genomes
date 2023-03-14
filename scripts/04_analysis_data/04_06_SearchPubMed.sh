set -e 

SearchOnPubMed() {

    local INPUT=$1
    local OUTPUT=$(echo $2 | sed 's/\/*$/\//g')
    local SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
    
    local INPUT_BASENAME=$(basename $1)

    printf "  ${GREEN_COLOR}Searching on PubMed${RESET_COLOR}\n\n"

    [ ! -d $OUTPUT ] && mkdir -p $OUTPUT

    Rscript --vanilla $SCRIPT_DIR/04_06_01_SearchPubMed.R -i $INPUT -o $OUTPUT | tee "${OUTPUT}${INPUT_BASENAME}_stdoutput"

}


# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        SearchOnPubMed

fi