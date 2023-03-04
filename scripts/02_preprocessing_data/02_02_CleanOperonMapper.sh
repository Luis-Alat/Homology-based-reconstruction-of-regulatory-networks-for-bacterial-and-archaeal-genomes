#!/bin/bash

set -e

TransformOperonMapperFiles() {

    # This function parser the operons to be decribed in a single-line
    # This function transform the common name id into the NCBI id
    # This function is expecting 3 arguments. In order
    # 1. Path where to find the operons obtained from operon mapper
    # 2. Path to place operon results parsed
    # 3. Ṕath to find the equivalences table between the common name and the NCBI id

    # NOTE THAT THE FILES OF THE OPERONS AND THE EQUIVALENCES ARE EXPECTED TO HAVE THE SAME NAME

    local OPERON_RESULTS=$(echo $1 | sed 's/\/*$/\//g')
    local OUTPUT=$(echo $2 | sed 's/\/*$/\//g')
    local EQUIV_PATH_TABLE=$(echo $3 | sed 's/\/*$/\//g')
    local TXT
    local BATCHES=10
    local COUNTER_BATCHES=0
    local COUNTER_TMP=1

    [ -d $OUTPUT ] && rm -r $OUTPUT; mkdir $OUTPUT || mkdir $OUTPUT

    for TXT in ${OPERON_RESULTS}*; do
        
        (OUTPUT_FILE_NAME=$(basename $TXT)
        mkdir "${OUTPUT}tmp_${COUNTER_TMP}/"

        printf "Processing ${OUTPUT_FILE_NAME}\n"

        # Replace the common name by the ncbi id
        awk -F"\t" 'BEGIN{OFS=FS}{print $3,$1}' "${EQUIV_PATH_TABLE}${OUTPUT_FILE_NAME}" |
            sed -r 's/^/s\/\\</g; s/\t/\\>\//g; s/$/\/gI;/g; s/\./\\./g' | 
            sed -f - <(tail -n +2 $TXT | cut -f1,2,7) > "${OUTPUT}tmp_${COUNTER_TMP}/tmp_operon_mapper"

        # Removing names not matching with the NCBI id
        grep -Fvif <(tail -n +2 $TXT | grep -P "^\s+" | cut -f2) "${OUTPUT}tmp_${COUNTER_TMP}/tmp_operon_mapper" > "${OUTPUT}tmp_${COUNTER_TMP}/tmp"

        # Transform multiples lines describing elements in the same one operon into one describing the elements in the operon
        cut -f1,2 "${OUTPUT}tmp_${COUNTER_TMP}/tmp" | 
            sed ':r;$!{N;br};s/\n\t/,/g' | 
            sed -r 's/,/\t/1' | 
            grep -vP "^\d+$" > "${OUTPUT}tmp_${COUNTER_TMP}/operon_group_elements"
        
        cut -f1,3 "${OUTPUT}tmp_${COUNTER_TMP}/tmp" | 
            sed ':r;$!{N;br};s/\n\t/,/g' | 
            sed -r 's/,/\t/1'  | 
            grep -vP "^\d+$" > "${OUTPUT}tmp_${COUNTER_TMP}/operon_strand_type"

        # Join elements of the operon with the strand
        bash 02_02_01_sqlite_join.sh "${OUTPUT}tmp_${COUNTER_TMP}" | 
            sed -r 's/,(\+|\-)+//g' > "${OUTPUT}${OUTPUT_FILE_NAME}") &

        ((++COUNTER_TMP))
        ((++COUNTER_BATCHES)); [ $COUNTER_BATCHES -eq $BATCHES ] && wait || COUNTER_BATCHES=0 

    done

    wait
    
    rm -r "${OUTPUT}tmp_"*

}

# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        TransformOperonMapperFiles "../../tus_models_operon/" "../../tus_models_operon_processed/" "../../others/ids_equivalences/"

fi