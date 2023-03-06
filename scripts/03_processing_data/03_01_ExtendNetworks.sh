#!/bin/bash

set -e

ExtendNetworksByOtho() {

        # This function creates the extended networks by mapping the orthologs retrieved in proteinortho
        # This function is expecting 5 arguments. In order
        # 1. String to put output files (passed by value)
        # 2. Array containing the labels identifying each fasta file, network or organism (passed by reference)
        # 3. String indicating where to find the proteinortho output (passed by value)
        # 4. Array containing the path and name of the fasta files (passed by reference)
        # 5. Array containing the path and name of the network files (passed by reference)

        printf "${GREEN_COLOR}  Extended network block by orthologs${RESET_COLOR}\n\n"

        local i
        local FILE_PATH_PROTEINORTHO_BLAST_GRAPH
        local OUTPUT=$(echo $1 | sed 's/\/*$/\//g') # Before It was net/
        local -n FOLDER_LABELS=$2
        local INPUT_PROTEINORTHO=$(echo $3 | sed 's/\/*$/\//g')
        local -n GENOMES=$4
        local -n NETWORKS=$5

        echo "    Checking if ${OUTPUT} exits"
        if [ ! -d $OUTPUT ]; then
                echo "      Creating ${OUTPUT}"
                mkdir $OUTPUT
        else
                rm -r $OUTPUT
                mkdir $OUTPUT
        fi

        for ((i=0; i < ${#FOLDER_LABELS[@]}; i++)); do

                # Creating folder if it doesn't exist
                # In this case, we need two, one to save the relations among orthologs,
                # and the another one to save a draft of the interactions predicted

                local REFERENCE_OUTPUT_PATH_CURRENT_ORG=$(printf "${OUTPUT}${FOLDER_LABELS[$i]}")
                local REFERENCE_OUTPUT_PATH_ORTHO_FILES=$(printf "${REFERENCE_OUTPUT_PATH_CURRENT_ORG}/ortoFiles")
                local REFERENCE_OUTPUT_PATH_PRE_NEW_INT=$(printf "${REFERENCE_OUTPUT_PATH_CURRENT_ORG}/PreNewInteractionsFiles")

                if [ ! -d "${OUTPUT}${FOLDER_LABELS[$i]}" ]; then
                        mkdir $REFERENCE_OUTPUT_PATH_CURRENT_ORG $REFERENCE_OUTPUT_PATH_ORTHO_FILES $REFERENCE_OUTPUT_PATH_PRE_NEW_INT
                else
                        rm -r "${OUTPUT}${FOLDER_LABELS[$i]}"
                        mkdir $REFERENCE_OUTPUT_PATH_CURRENT_ORG $REFERENCE_OUTPUT_PATH_ORTHO_FILES $REFERENCE_OUTPUT_PATH_PRE_NEW_INT         
                fi

                local CURRENT_PROTEINORTHO_FOLDER_FILES=$(echo "${INPUT_PROTEINORTHO}${FOLDER_LABELS[$i]}/*blast-graph")

                # Making a tmp file to remove interactions at the end of the following for loop
                # Some interactions are described by multiples TFs, the reason of this file is to remove
                # those interactions if one of the TFs was not found in the target organism
                # at the moment of doing an extrapolation of the interactions among organisms
                
                grep ">" "${GENOMES[$i]}" |
                        cut -d"|" -f1 |
                        sed -r 's/>//g' > "${REFERENCE_OUTPUT_PATH_PRE_NEW_INT}/tmp"

                # Iterating in every "blast-graph" file found inside of ${FolderNames[@]}
                for FILE_PATH_PROTEINORTHO_BLAST_GRAPH in $CURRENT_PROTEINORTHO_FOLDER_FILES; do

                        local OUTPUT_FILE_BASENAME=$(basename $FILE_PATH_PROTEINORTHO_BLAST_GRAPH)
                        local OUTPUT_FILE_NAME_ORTHO=$(printf "${OUTPUT_FILE_BASENAME}_ortho")
                        local OUTPUT_FILE_NAME_PRE_NEW_INTER=$(echo "${OUTPUT_FILE_BASENAME}_new_interactions")

                        printf "    Making files from ${FILE_PATH_PROTEINORTHO_BLAST_GRAPH}\n"

                        # Getting an array for a correct order of the columns in the proteinortho files
                        # First column organism model (E.coli, M.tuberculosis with Known network) and the
                        # second column it the target organism

                        # Error es en la siguiente linea, no busca bien el patrón o nombre del archivo en el header
                        # primero ecoli luego el otro
                        local CORRECT_PROTEINORTHO_COLUMN_ORDER=$(head -n 5 $FILE_PATH_PROTEINORTHO_BLAST_GRAPH |
                                awk -F$"\t" -v pattern=$(basename ${GENOMES[$i]}) '$1 ~ pattern {print 1","2} $2 ~ pattern {print 2","1}')

                        # The output above now will be an array
                        local CORRECT_PROTEINORTHO_COLUMN_ORDER_ARR
                        IFS="," read -r -a CORRECT_PROTEINORTHO_COLUMN_ORDER_ARR <<< "${CORRECT_PROTEINORTHO_COLUMN_ORDER}"

                        # Reording the columns and removing comments and extra columns
                        paste <(grep -v "#" $FILE_PATH_PROTEINORTHO_BLAST_GRAPH | cut -f${CORRECT_PROTEINORTHO_COLUMN_ORDER_ARR[0]}) <(grep -v "#" $FILE_PATH_PROTEINORTHO_BLAST_GRAPH | cut -f${CORRECT_PROTEINORTHO_COLUMN_ORDER_ARR[1]}) |
                                sed -r 's/\|/\t/g' | 
                                cut -f1,3 > "${REFERENCE_OUTPUT_PATH_ORTHO_FILES}/${OUTPUT_FILE_NAME_ORTHO}"

                        # Given the new file above, we are going to use it to identify interactions
                        # Creating multiples regex commands to replace in the network of the organism
                        # models (mapping orthologues) by the NCBI_id of the target organism.
                        # Every regex command is going to have a format NCBI_ID_model,NCBI_ID_targetMATCH
                        # The keyword MATCH is going to be useful to retrieve what interactions were mapped
                        sed -r 's/^/s\/\\</g' "${REFERENCE_OUTPUT_PATH_ORTHO_FILES}/${OUTPUT_FILE_NAME_ORTHO}" | 
                                sed -r 's/\t/\\>\//g' | 
                                sed -r 's/$/MATCH\/gI;/g' | 
                                sed -r 's/\./\\./g' | 
                                sed -f - <(cut -f1-3 ${NETWORKS[$i]}) | 
                                perl -F"\t" -nae 'if(($F[1] =~ /MATCH/) and ($F[2] =~ /MATCH/)){print $_}' | 
                                sed -r 's/MATCH//g' | 
                                awk -F $'\t' -v organism=${FOLDER_LABELS[$i]} 'BEGIN{OFS = FS} {print organism,$0}' | 
                                grep -Fvwif "${REFERENCE_OUTPUT_PATH_PRE_NEW_INT}/tmp" - > "${REFERENCE_OUTPUT_PATH_PRE_NEW_INT}/$OUTPUT_FILE_NAME_PRE_NEW_INTER"

                        done
        done

        # Creating array of file names, It's useful to iterate
        declare -a local BASENAMES_GENOMES 
        local FASTA_FILE

        for FASTA_FILE in "${GENOMES[@]}"; do
                BASENAMES_GENOMES+=($(basename $FASTA_FILE))
        done

        # Creating folder if it doesn't exist
        [ -d "${OUTPUT}Merge" ] && rm -r "${OUTPUT}Merge"; mkdir "${OUTPUT}Merge" || mkdir "${OUTPUT}Merge"

        # Merging results ("for" loop below) of every Folder (${FolderNames[@]}) using the name base founds inside of "proteinortho/${FolderNames[@]}/"
        # Every final file generated in the for loop above contains ...
        # organism_name_where_interaction_came_from, organism_interaction_number_where_interaction_came_from, TF, TG by model and target organism spread ...
        # among different folders. This for loop merge the respective files
        # Results will be saved in the folder "net/Merge"
        local FASTA_FILE_BASENAME
        for FASTA_FILE_BASENAME in ${BASENAMES_GENOMES[@]}; do

                # Defining new file name and making empty file
                MERGE_FILE_NAME=$(echo "${FASTA_FILE_BASENAME}_Merged")
                echo -n > "${OUTPUT}/Merge/tmp"

                # Merging files into the "Merge" folder
                local ORGANISM_FOLDER
                for ORGANISM_FOLDER in ${FOLDER_LABELS[@]}; do

                        printf "    ${FASTA_FILE_BASENAME}\t${ORGANISM_FOLDER}\n"
                        INPUT_FILE=$(echo "${OUTPUT}${ORGANISM_FOLDER}/PreNewInteractionsFiles/${FASTA_FILE_BASENAME}*")

                        # Check if file exist in that folder
                        if [ -f $INPUT_FILE ]; then
                                # Merging all results into a temporaly file
                                cat $INPUT_FILE >> "${OUTPUT}Merge/tmp"
                        fi

                done

                printf "    Making file ${MERGE_FILE_NAME}\n"

                # This line calls an sql implementation to group the columns associated with regulatory interactions
                # It requires in order: Col names, GROUP BY col names, SELECT instruction and file with tha table (tsv format)

                bash utils/sqlite_grouper.sh "tf, tg, org_ref, line_num" "org_ref, line_num" "org_ref, line_num, GROUP_CONCAT(tf) , GROUP_CONCAT(tg)" "${OUTPUT}/Merge/tmp" | 
                        perl -F"\t" -nae '
                                chomp($_);
                                @array=split(",",$F[2]);
                                %seen=();
                                @unique = grep { ! $seen{$_} ++ } @array;
                                $len = scalar(@unique);
                                print ($_,"\t",$len,"\n")' > "${OUTPUT}Merge/${MERGE_FILE_NAME}"

                rm "${OUTPUT}/Merge/tmp"

        done

        # In the following "for" loop an "extended network" is going to be created for every model organism...
        # It has the interactions described in the modified network as well as the new interactions found
        # Results will be saved in the folder "net/results"

        # Creating folder if it doesn't exist        
        if [ ! -d "${OUTPUT}results" ]; then
                mkdir "${OUTPUT}results" "${OUTPUT}results/tmp"
        else
                rm -r "${OUTPUT}results"; mkdir "${OUTPUT}results" "${OUTPUT}results/tmp"
        fi

        for ((i=0; i < ${#GENOMES[@]}; i++)); do

                OUTPUT_FILE_TMP=$(echo "${BASENAMES_GENOMES[$i]}_TMP")

                printf "    Creating extended network of ${GENOMES[$i]} using ${NETWORKS[$i]}\n"

                # Giving format both $MergeFile and network and later they are concatenated...
                # Any missing info is added too
                ## Retrieving TF and TG columns, adding info about the state of the interaction (does the interactions is known?)
                cut -f1-3 ${NETWORKS[$i]} |
                        awk -F$"\t" 'BEGIN{OFS=FS}{print $0,"Known","NOT_REFR_ORG","NOT_REFR_NUM","NULL"}' > "${OUTPUT}results/tmp/${OUTPUT_FILE_TMP}"

                # Adding info about if the interaction is new and a new enumerate id
                awk -F$'\t' 'BEGIN{OFS=FS} $2=$2 FS "New?"' "${OUTPUT}Merge/${BASENAMES_GENOMES[$i]}"* |
                        perl -nse 'print($count,"_new?","\t",$_); $count+=1' -- -count=1 >> "${OUTPUT}results/tmp/${OUTPUT_FILE_TMP}"

                # Grouping (group_by sql) results and removing "innecesary" information (labels no longer used)
                OUTPUT_FILE_NAME_EXT_NET=$(echo "${BASENAMES_GENOMES[$i]}_extended_network")
                
                bash utils/sqlite_grouper.sh "line_number, tg, tf, status, org_ref, numeric_ref, total_ref" "tg, tf"  "GROUP_CONCAT(line_number), tg, tf, GROUP_CONCAT(status), GROUP_CONCAT(org_ref), GROUP_CONCAT(numeric_ref), GROUP_CONCAT(total_ref)" "${OUTPUT}results/tmp/${OUTPUT_FILE_TMP}" |
                        ## Removing innecesary labels, for example, an interaccion "Known" could be also described by the mapping approach ...
                        ## therefore, that interaction will have the label "Known" as "New", but we already know thats not new. "New" is deleted
                        sed -r 's/(,\w+_new\?)+//1' | sed -r 's/_new\?/_new/1' |
                        sed -r 's/,New\?//g' | sed -r 's/\tNew\?\t/\tNew\t/g' |
                        sed 's/NOT_REFR_ORG,//g' | sed 's/NOT_REFR_NUM,//g' |
                        sed -r 's/(NULL,)+//g' | sed -r 's/(Known,)+//g' | sort -k2,3 > "${OUTPUT}results/${OUTPUT_FILE_NAME_EXT_NET}"

        done

}

ExtendNetworksByTranscriptionUnit() {

        # This function integrates the TU information into an existing network
        # This function is expecting 3 arguments. In order
        # 1. String indicating where to save new networks; This function also will search the networks in this path (passed by value)
        # 2. Array indicating the path ad file names of the genomes or organism (passed by reference)
        # 3. String indicating where the transcription units are; Those files should have the same name as the genomes (passed by value)

        printf "${GREEN_COLOR}  Extended network block by Transcription Units (TU)${RESET_COLOR}\n\n"

        local OUTPUT=$(echo $1 | sed 's/\/*$/\//g')
        local -n GENOMES=$2
        local TU_PATH_FILES=$(echo $3 | sed 's/\/*$/\//g')
        local SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

        # Creating a bash array and get ID names
        declare -a local BASENAMES_GENOMES_NO_EXT
        declare -a local BASENAMES_GENOMES
        local BASE

        for FASTA_FILE in ${GENOMES[@]}; do
                BASE=$(basename $FASTA_FILE)
                BASENAMES_GENOMES+=($BASE)
                BASENAMES_GENOMES_NO_EXT+=(${BASE%.*})
        done

        # If folder doesn't exit, create it
        if [ ! -d "${OUTPUT}results_plus_TU" ]; then
                mkdir "${OUTPUT}results_plus_TU" "${OUTPUT}results_plus_TU/tmp"
        else
                rm -r "${OUTPUT}results_plus_TU" "${OUTPUT}results_plus_TU/tmp"
                mkdir "${OUTPUT}results_plus_TU" "${OUTPUT}results_plus_TU/tmp"
        fi

        local i
        for ((i=0; i < ${#BASENAMES_GENOMES[@]}; i++)); do

                # Verifying if the file exists (Does exist the TU file for this particular model organism?)
                printf "${TU_PATH_FILES}${BASENAMES_GENOMES_NO_EXT[$i]}"*
                if [ -f "${TU_PATH_FILES}${BASENAMES_GENOMES_NO_EXT[$i]}"* ]; then

                        # Defining an output final file name for the python script below
                        local OUTPUT_FILE_NAME_TUS=$(echo "${BASENAMES_GENOMES[$i]}_extended_network_plus_tu")

                        # Searching TUs where the first gene is equal to the target in the extended network
                        python ${SCRIPT_DIR}/03_01_01_ParserTUandExt.py --input_net "${OUTPUT}results/${BASENAMES_GENOMES[$i]}"* --input_tus "${TU_PATH_FILES}${BASENAMES_GENOMES_NO_EXT[$i]}"* --output_path_name "${OUTPUT}results_plus_TU/tmp/${BASENAMES_GENOMES[$i]}_tmp" --arguments_read_net '{"header":None, "sep":"\t", "usecols":[1,2]}' --arguments_read_tus '{"header":None, "sep":"\t"}'

                        printf "  Making ${OUTPUT}results_plus_TU/${OUTPUT_FILE_NAME_TUS}\n"

                        # Giving a new format to the earlier output to match with the extended network format. Saving in a tmp_1 file
                        awk -F$"\t" 'BEGIN{OFS=FS}{print $4,$2,$3,"New?","NOT_REFR_ORG","NOT_REFR_NUM","NULL",$1,$4}' "${OUTPUT}results_plus_TU/tmp/${BASENAMES_GENOMES[$i]}_tmp" > ${OUTPUT}results_plus_TU/tmp_1

                        # Giving a new format to the extended network to match with the output python file above. Saving in a tmp_2 file 
                        awk -F$"\t" 'BEGIN{OFS=FS}{print $0,"NOT_STRAND","NOT_TU_REFER"}' "${OUTPUT}results/${BASENAMES_GENOMES[$i]}"* > ${OUTPUT}results_plus_TU/tmp_2

                        # Merging results of both files
                        cat ${OUTPUT}results_plus_TU/tmp_2 ${OUTPUT}results_plus_TU/tmp_1 > ${OUTPUT}results_plus_TU/tmp_3

                        # Removing information not relevant (labels no longer used) useful during the group_by, but not later
                        bash utils/sqlite_grouper.sh "id,tf,tg,status,org_ref_name,org_ref_line,count,strand,id_tus" "tf,tg" "GROUP_CONCAT(id),tf,tg,GROUP_CONCAT(status),GROUP_CONCAT(org_ref_name),GROUP_CONCAT(org_ref_line),GROUP_CONCAT(count),GROUP_CONCAT(strand),GROUP_CONCAT(id_tus)" "${OUTPUT}results_plus_TU/tmp_3" |
                                sed -r 's/,New\?//g' | sed -r 's/,NOT_REFR_ORG//g' |
                                sed -r 's/,NOT_REFR_NUM//g' | sed -r 's/,NULL//g' |
                                sed -r 's/NOT_STRAND,//g' | sed -r 's/NOT_TU_REFER,//g' |
                                sed -r 's/,\w+_tu//g' | sed -r 's/New\?/New/g' > "${OUTPUT}results_plus_TU/${OUTPUT_FILE_NAME_TUS}"

                        printf "  All done\n\n"

                fi

        done

        # Removing tmp files generated above
        rm ${OUTPUT}results_plus_TU/tmp_1 ${OUTPUT}results_plus_TU/tmp_2 ${OUTPUT}results_plus_TU/tmp_3

}



# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        ExtendNetworksByOtho
        ExtendNetworksByTranscriptionUnit

fi
