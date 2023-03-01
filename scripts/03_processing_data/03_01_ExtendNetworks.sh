#!/bin/bash

set -e

ExtendNetworks() {

        # This function creates the extended networks by mapping the orthologs retrieved in proteinortho
        # This function is expecting 5 arguments. In order
        # 1. String to put output files (passed by value)
        # 2. Array containing the labels identifying each fasta file, network or organism (passed by reference)
        # 3. String indicating where to find the proteinortho output (passed by value)
        # 4. Array containing the path and name of the fasta files (passed by reference)
        # 5. Array containing the path and name of the network files (passed by reference)

        printf "${GREEN_COLOR}  Extended network block${RESET_COLOR}\n\n"

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

                ## This line calls an sql implementation to group the columns associated with regulatory interactions to describe
                ##(1) transcription factor,
                ##(2) target,
                ##(3) organism where that interaction was found and
                ##(4), line number associated with that interaction found in modified networks...
                ##In addition adds a number indicating in how many organism that interactions was found

                exit 1

                bash scripts/sqlite_join.sh "${OUTPUT}/Merge/tmp"  |
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
        if [ ! -d net/results ];then
                mkdir net/results net/results/tmp
        fi

        for ((i=0;i<${#GenomeNames[@]};i++)); do

                # Getting base name
                BaseName=$(echo ${GenomeNames[$i]} | sed -r 's/\.faa.*//g')
                OutPutFileTmp=$(echo $BaseName"TMP")

                printf "Creating result files of "${GenomeNames[$i]}" using the modified network "${NetFiles[$i]}"\n"

                # Giving format both $MergeFile and the modified network and later they are concatenated...
                # Any missing info is added too
                ## Retrieving TF and TG columns, adding info about the state of the interaction (does the interactions is known?)
                cut -f1-3 modifiedNetworks/${NetFiles[$i]} |
                        awk -F$"\t" 'BEGIN{OFS=FS}{print $0,"Known","NOT_REFR_ORG","NOT_REFR_NUM","NULL"}' > net/results/tmp/$OutPutFileTmp

                # Adding info about if the interaction is new and a new enumerate id
                awk -F$'\t' 'BEGIN{OFS=FS} $2=$2 FS "New?"' net/Merge/$BaseName* |
                        perl -nse 'print($count,"_new?","\t",$_); $count+=1' -- -count=1 >> net/results/tmp/$OutPutFileTmp

                # Grouping (group_by sql) results and removing "innecesary" information (labels no longer used)
                OutPutFile=$(echo $BaseName"_extended_network.txt")

                bash scripts/sqlite_groupe_results.sh net/results/tmp/$OutPutFileTmp |
                        ## Removing innecesary labels, for example, an interaccion "Known" could be also described by the mapping approach ...
                        ## therefore, that interaction will have the label "Known" as "New", but we already know thats not new. "New" is deleted
                        sed -r 's/(,\w+_new\?)+//1' | sed -r 's/_new\?/_new/1' |
                        sed -r 's/,New\?//g' | sed -r 's/\tNew\?\t/\tNew\t/g' |
                        sed 's/NOT_REFR_ORG,//g' | sed 's/NOT_REFR_NUM,//g' |
                        sed -r 's/(NULL,)+//g' | sed -r 's/(Known,)+//g' | sort -k2,3 > net/results/$OutPutFile

        done

}

# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following calling for now is expecting the values directly and not by command line
        ExtendNetworks

fi
