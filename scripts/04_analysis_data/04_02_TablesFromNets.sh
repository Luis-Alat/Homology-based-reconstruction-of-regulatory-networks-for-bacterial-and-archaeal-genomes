#!/bin/bash

set -e

CreateNetworkTableMetrics() {

        # Based on the extended network with info from Tus. Create two general tables about
        # how many interactions the new network was extended and where that interactions came from and
        # how many interactions rebuild the modified network and where that interactions came from
        # This functions is expecting 5 arguments. In order
        # 1. Array of the labels of the organism (passed by reference)
        # 2. Array of the paths to find the genomes organis fasta files (passed by reference)
        # 3. String indicating the path to find the networks (passed by value)
        # 4. String indicating tha path to save the tables (passed by value)

        printf "${GREEN_COLOR}  Creating tables${RESET_COLOR}\n\n"

        local -n ORG_LABELS=$1
        local -n GENOMES=$2
        local INPUT=$(echo $3 | sed 's/\/*$/\//g')
        local OUTPUT=$(echo $4 | sed 's/\/*$/\//g')

        [ ! -d $OUTPUT ] && mkdir $OUTPUT

        for ((i=0; i < ${#GENOMES[@]}; i++)); do
                GENOMES[$i]=$(basename ${GENOMES[$i]})
        done

        # Creating headers for the first table
        printf "  Creating headers and table files\n"
        local header=","
        local ORG
        for ORG in ${ORG_LABELS[@]}; do
                header=$(printf "${header}${ORG},")
        done
        printf "${header}TU,Total of new interactions\n" > "${OUTPUT}Table_of_percentage_of_news.csv"
        printf "${header}TU,Total of reconstructions,Total Known interactions\n" > "${OUTPUT}Table_of_percentage_of_reconstruction.csv"

        for ((i=0; i < ${#ORG_LABELS[@]}; i++)); do

                local VALUES_ROW_NEW_TABLE=()
                local VALUES_ROW_REC_TABLE=()

                local COUNT_NEW_INTERACTIONS=$(grep "New" "${INPUT}${GENOMES[$i]}"* | wc -l)
                local TU=$(grep "New" "${INPUT}${GENOMES[$i]}"* | grep -v "NOT_TU_REFER" | wc -l)
                local COUNT_NEW_TU=$(echo "scale=2; (${TU}*100)/${COUNT_NEW_INTERACTIONS}" | bc)
                local NEW_TU_FORMAT=$(echo -n "${TU} (${COUNT_NEW_TU}%)")

                local COUNT_KNOWN_RECONSTRUCTION=$(grep "Known" "${INPUT}${GENOMES[$i]}"* | grep -vP "(?=.*NOT_REFR_ORG)(?=.*NOT_TU_REFER)" | wc -l)
                local COUNT_KNOWN_INTERACTIONS=$(grep "Known" "${INPUT}${GENOMES[$i]}"* | wc -l)
                local TU=$(grep "Known" "${INPUT}${GENOMES[$i]}"* | grep -v "NOT_TU_REFER" | wc -l)
                local COUNT_KNOWN_TU=$(echo "scale=2; (${TU}*100)/${COUNT_KNOWN_RECONSTRUCTION}" | bc)
                local KNOWN_TU_FORMAT=$(echo -n "${TU} (${COUNT_KNOWN_TU}%)")

                for ((j=0; j < ${#ORG_LABELS[@]}; j++)) do

                        # Don't search againts itself. Add '---' inplace
                        if [ $i -eq $j ]; then

                                local CURRENT_COUNT=$(echo -n "---,")
                                local VALUES_ROW_NEW_TABLE+=($CURRENT_COUNT)
                                local VALUES_ROW_REC_TABLE+=($CURRENT_COUNT)

                        else

                                # Getting new interactions per organism and number of it
                                local ORG_VALUE=$(grep "New" "${INPUT}${GENOMES[$i]}"* | grep ${ORG_LABELS[$j]} | wc -l)
                                
                                # Getting percentage
                                local PCG_VALUE=$(echo "scale=2; (${ORG_VALUE}*100)/${COUNT_NEW_INTERACTIONS}" | bc)
                                
                                # Giving format to save in array, useful later to save results into file
                                CURRENT_COUNT=$(echo -n "${ORG_VALUE} (${PCG_VALUE}%),")
                                VALUES_ROW_NEW_TABLE+=($CURRENT_COUNT)

                                # Similar as above, but now with the Known interactions
                                ORG_VALUE=$(grep "Known" "${INPUT}${GENOMES[$i]}"* | grep ${ORG_LABELS[$j]} | wc -l)
                                PCG_VALUE=$(echo "scale=2; (${ORG_VALUE}*100)/${COUNT_KNOWN_RECONSTRUCTION}" | bc)
                                
                                CURRENT_COUNT=$(echo -n "${ORG_VALUE} (${PCG_VALUE}%),")
                                VALUES_ROW_REC_TABLE+=($CURRENT_COUNT)

                        fi

                done

                # Append in files created above
                echo "${ORG_LABELS[$i]},${VALUES_ROW_NEW_TABLE[*]}${NEW_TU_FORMAT},${COUNT_NEW_INTERACTIONS}" >> "${OUTPUT}Table_of_percentage_of_news.csv"
                echo "${ORG_LABELS[$i]},${VALUES_ROW_REC_TABLE[*]}${KNOWN_TU_FORMAT},${COUNT_KNOWN_RECONSTRUCTION},${COUNT_KNOWN_INTERACTIONS}" >> "${OUTPUT}Table_of_percentage_of_reconstruction.csv"
        done

}


# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        CreateNetworkTableMetrics

fi
