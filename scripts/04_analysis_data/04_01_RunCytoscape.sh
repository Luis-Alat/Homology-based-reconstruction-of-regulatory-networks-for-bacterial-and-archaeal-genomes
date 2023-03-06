#!/bin/bash

set -e

AnalyzeByCytoscape() {

        local -n GENOMES=$1
        local OUTPUT=$(echo $2 | sed 's/\/*$/\//g')
        local PATH_NET_ORTHO=$(echo $3 | sed 's/\/*$/\//g')
        local PATH_NET_PLUS_TU=$(echo $4 | sed 's/\/*$/\//g')
        local SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

        printf "${GREEN_COLOR}  CytoScape block${RESET_COLOR}\n\n"

        # It useful to assign an output name and also to find files the following array
        declare -a local BASE_NAME
        readarray -t GENOMES_WITHOUT_EXTENTION < <(printf "%s\n" ${GENOMES[@]%.*})

        local FASTA
        for FASTA in ${GENOMES_WITHOUT_EXTENTION[@]}; do
                BASE_NAME+=( $(basename $FASTA) )
        done

        # Creating folder if it doesn't exist
        if [ ! -d $OUTPUT ]; then
            mkdir $OUTPUT "${OUTPUT}tmp"
        else
            rm -r $OUTPUT
            mkdir $OUTPUT "${OUTPUT}tmp"
        fi

        cd $OUTPUT

        for ((i=0; i < ${#GENOMES[@]}; i++)); do

                local INPUT_ORTHO_NETS=$(printf "${PATH_NET_ORTHO}${BASE_NAME[$i]}"*)
                local INPUT_PLUS_TU_NETS=$(printf "${PATH_NET_PLUS_TU}${BASE_NAME[$i]}"*)

                # Calling the R script to execute commands on CytoScape (API), generating and analyzing networks
                Rscript --vanilla $SCRIPT_DIR/04_01_01_CytoScapeBuildNetworks.R -i "${INPUT_ORTHO_NETS}" -t new -p $OUTPUT
                Rscript --vanilla $SCRIPT_DIR/04_01_01_CytoScapeBuildNetworks.R -i "${INPUT_PLUS_TU_NETS}" -t new -p $OUTPUT

                # Adding additional (inplace) information in file of statistic...
                # (cut -f2 gets the number of TFs in that network, cut -f3 gets the target)...
                # from modified network and the extended network

                local BASENAME_OUTPUT_CYTO_ORTHO=$(basename ${INPUT_ORTHO_NETS})
                local BASENAME_OUTPUT_CYTO_TUS=$(basename ${INPUT_PLUS_TU_NETS})

                local TFCounts=$(cut -f2 "${INPUT_ORTHO_NETS}" | sort | uniq | wc -l | sed -r 's/^/#TF_counts\t/g')
                sed -i "2i\\$TFCounts" "${OUTPUT}${BASENAME_OUTPUT_CYTO_ORTHO}"*"metrics"
                local TargetCounts=$(cut -f3 "${INPUT_ORTHO_NETS}" | sort | uniq | wc -l | sed -r 's/^/#Target_counts\t/g')
                sed -i "2i\\$TargetCounts" "${OUTPUT}${BASENAME_OUTPUT_CYTO_ORTHO}"*"metrics"

                TFCounts=$(cut -f2 "${INPUT_PLUS_TU_NETS}" | sort | uniq | wc -l | sed -r 's/^/#TF_counts\t/g')
                sed -i "2i\\$TFCounts" "${OUTPUT}${BASENAME_OUTPUT_CYTO_TUS}"*"metrics"
                TargetCounts=$(cut -f3 "${INPUT_PLUS_TU_NETS}" | sort | uniq | wc -l | sed -r 's/^/#Target_counts\t/g')
                sed -i "2i\\$TargetCounts" "${OUTPUT}${BASENAME_OUTPUT_CYTO_TUS}"*"metrics"

        done

}

# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        AnalyzeByCytoscape
        CreateCytoTables

fi
