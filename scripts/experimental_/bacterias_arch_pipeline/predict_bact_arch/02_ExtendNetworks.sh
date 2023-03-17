#!/bin/bash

set -e 

ExtendNetworksByOtho() {

    # This function create the regulatory networks of organism without a
    # transcriptional regulatiory network by using as reference organisms
    # which have one
    # This function is expecting 5 arguments or 6 optionally. In order
    # 1. String indicating path to save results of this 
    #    function (passed by value)
    # 2. Array containing the labels to identify the organism
    #    used as reference (passed by reference)
    # 3. String indicating the path to find the proteinortho 
    #    results (passed by value)
    # 4. Array containing the path and name of the fasta file 
    #    of the organism without network (passed by reference)
    # 5. Array containing the path and name of the networks of
    #    the organism used as reference (passed by reference)
    # 6. Integer indicating how many processes in parallel to
    #    certain sections of this function will be done (passed by value)

    printf "${GREEN_COLOR}  Extended network block by orthologs${RESET_COLOR}\n\n"

    local OUTPUT=$(echo $1 | sed 's/\/*$/\//g')
    local -n LABELS=$2
    local INPUT_PROTEINORTHO=$(echo $3 | sed 's/\/*$/\//g')
    local -n TARGET_GENOMES=$4
    local -n NETWORKS=$5

    # Run for in N batches
    [ -z $6 ] && local BATCHES=1 || local BATCHES=$6
    local COUNT_BATCHES=1

    # Current script path
    SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

	# Creating folder if it doesn't exist
	[ ! -d "${OUTPUT}" ] && mkdir -p "${OUTPUT}"

	for ((i=0; i < ${#LABELS[@]}; i++)); do

        # Creating folder if it doesn't exist
        if [ ! -d "${OUTPUT}${LABELS[$i]}" ]; then
            mkdir "${OUTPUT}${LABELS[$i]}"
            mkdir "${OUTPUT}${LABELS[$i]}/ortoFiles/" "${OUTPUT}${LABELS[$i]}/PreNewInteractionsFiles/"
        fi

		# Defining paths to place files
        local REFERENCE_OUTPUT_PATH_ORTHO_FILES=$(echo "${OUTPUT}${LABELS[$i]}/ortoFiles/")
        local REFERENCE_OUTPUT_PATH_PRE_NEW_INT=$(echo "${OUTPUT}${LABELS[$i]}/PreNewInteractionsFiles/")

        # Iterating in every "blast-graph" file created by proteinortho using the names of the target organism
		local CURRENT_NUMBER=1
		for ORG in ${TARGET_GENOMES[@]}; do

            local TARGET_BASENAME=$(basename $ORG)

            if [ ! -f "${INPUT_PROTEINORTHO}${LABELS[$i]}/${TARGET_BASENAME}"*"blast-graph" ]; then

                printf "  ${YELLOW_COLOR} File ${INPUT_PROTEINORTHO}${LABELS[$i]}/${TARGET_BASENAME}*blast-graph was not found${RESET_COLOR}\n"
                continue

            fi

			local BLAST_FILE_PROTEINORTHO=$(ls ${INPUT_PROTEINORTHO}${LABELS[$i]}/${TARGET_BASENAME}*blast-graph)

			# New files to be created
            local OUTPUT_FILE_NAME_ORTHO=$(printf "${TARGET_BASENAME}_ortho")
            local OUTPUT_FILE_NAME_INTER=$(printf "${TARGET_BASENAME}_new_interactions")

            printf "(background) Making ${CURRENT_NUMBER} of ${#TARGET_GENOMES[@]} files: ${BLAST_FILE_PROTEINORTHO}\n"

            # Given a specific file inside of the folder, search for a pattern (file name) and get the column number. Later it will be useful to reorder the file in the "correct" order (giving structure)
            # First 5 lines contain the necessary information ## Search by using the name file our column of interest and print the right order separated by comma
            local COLUMN_RIGTH_ORDER=$(head -n 5 $BLAST_FILE_PROTEINORTHO | awk -F$"\t" -v pattern=$TARGET_BASENAME '$1 ~ pattern {print 2","1} $2 ~ pattern {print 1","2}')

            IFS="," read -r -a COLUMN_RIGTH_ORDER_ARRAY <<< $COLUMN_RIGTH_ORDER

            # Reording the columns
            # Pasting in this order the column described in ${COLUMN_RIGTH_ORDER_ARRAY[0]} (model organism) and later ${COLUMN_RIGTH_ORDER_ARRAY[1]} (model target)
            # Replacing "|" by TAB and saving the id of the interaction and TF-TG
            paste <( grep -v "#" $BLAST_FILE_PROTEINORTHO | cut -f${COLUMN_RIGTH_ORDER_ARRAY[0]} ) <( grep -v "#" $BLAST_FILE_PROTEINORTHO | cut -f${COLUMN_RIGTH_ORDER_ARRAY[1]} ) | 
                sed -r 's/\|/\t/g' |
                cut -f1,3 > "${REFERENCE_OUTPUT_PATH_ORTHO_FILES}${OUTPUT_FILE_NAME_ORTHO}"

            # Given the new file, get a previous file useful to find the new interactions. For this purpose will be searched and replaced (in the modified networks) every match into a $OutFileOrtho format and add the organism where the interaction came from
            sed -r 's/^/s\/\\</g' "${REFERENCE_OUTPUT_PATH_ORTHO_FILES}${OUTPUT_FILE_NAME_ORTHO}" | 
                sed -r 's/\t/\\>\//g; s/$/MATCH\/gI;/g; s/\./\\./g' |
                sed -f - <( cut -f1-3 "${NETWORKS[$i]}" ) | 
                perl -F"\t" -nae 'if(($F[1] =~ /MATCH/) and ($F[2] =~ /MATCH/)){print $_}' | 
                sed -r 's/MATCH//g' | 
                awk -F $'\t' -v organism=${LABELS[$i]} 'BEGIN{OFS = FS} {print organism,$0}' > "${REFERENCE_OUTPUT_PATH_PRE_NEW_INT}${OUTPUT_FILE_NAME_INTER}" &

            ((++CURRENT_NUMBER))

			# Wait if there are N processes in background
			((++COUNT_BATCHES)); [ "${COUNT_BATCHES}" -eq "${BATCHES}" ] && COUNT_BATCHES=1 && wait
	
        done
	
    done

	wait

	# Creating folder if it doesn't exist
    [ ! -d "${OUTPUT}Merge/" ] && mkdir "${OUTPUT}Merge/"

    # Merging results ("for" loop below) of every Folder (${FolderNames[@]}) using the name found inside of "${TargetGenomes[@]}"
    # Results will be saved in the folder "${OutputPath}net/Merge/"

	for ORG in ${TARGET_GENOMES[@]}; do

        # Defining a new file name and making a temporally empty file
        TARGET_BASENAME=$(basename $ORG)
        PREDICTED_NET_FILE_NAME=$(echo "${TARGET_BASENAME}_Merged")
		
        echo -n > "${OUTPUT}Merge/tmp"

        # Add a column specifying the organism where the new interaction came from
        for ((i=1; i < ${#LABELS[@]}; i++)); do

            printf "  Creating a tmp file using ${TARGET_BASENAME} found in ${LABELS[$i]}\n"

            # Merging all results into a temporaly file
            cat "${OUTPUT}${LABELS[$i]}/PreNewInteractionsFiles/${TARGET_BASENAME}"* >> "${OUTPUT}Merge/tmp"

        done

        printf "  Making file ${PREDICTED_NET_FILE_NAME}\n"

        # This line calls an sql implementation to group the columns associated with regulatory interactions to describe (1) transcription factor, (2) target, (3) organism where that interaction was found and (4), line number associated with that interaction found in modified networks. In addition adds a number indicating in how many organism that interactions was found
        bash "${SCRIPT_DIR}/../../../utils/sqlite_grouper.sh" "LINE_NUMBER_REF, ORG_REF, TF, TG" "TF, TG" " TF, TG, group_concat(LINE_NUMBER_REF) , group_concat(ORG_REF)" "${OUTPUT}Merge/tmp" |
            perl -F"\t" -nae 'chomp($_);
                            @array=split(",", $F[2]);
                            %seen=();
                            @unique = grep { ! $seen{$_} ++ } @array;
                            $len = scalar(@unique);
                            print ($_,"\t",$len,"\n")' > "${OUTPUT}Merge/${PREDICTED_NET_FILE_NAME}"

        rm "${OUTPUT}Merge/tmp"

	done

}

ExtendNetworksByTranscriptionUnit() {

    printf "${GREEN_COLOR}  Using transcriptional units files${RESET_COLOR}\n\n"

    # Creating folder if it doesn't exist
    if [ ! -d "${OutputPath}net/Merge_plus_TU" ]; then
        mkdir "${OutputPath}net/Merge_plus_TU" "${OutputPath}net/Merge_plus_TU/tmp"
    fi

	# Count iterations
	count_batches=1
	# Run for in N batches
        Batches=500

	for file in ${TargetGenomes[@]}; do

		# Removing "faa" extention
		local file_cp=$(printf $file | sed -r 's/\.faa$//g')

        	# Verifying if the file exists
	        ( if [ -f $PathTu$file_cp* ]; then

			# Getting input files (full name)
			local InputExt=$(find "${OutputPath}net/Merge/${file_cp}"*)
			local InputTU=$(find "${PathTu}${file_cp}"*)
			local OutputF=$(printf "${OutputPath}net/Merge_plus_TU/tmp/${file}_net_predictions_TU")
			local OutputFinal=$(printf "${OutputPath}net/Merge_plus_TU/${file}_net_predictions_TU")

			# Maping "new" interactions
        	        python3 "${SCRIPTPATH}/ParserTUandExt.py" $InputExt $InputTU $OutputF 2

			# Giving a new format to the earlier output
			tmp_1=$(printf "${file}_tmp_1")
			awk -F$"\t" 'BEGIN{OFS=FS}{print $3,$4,"NOT_REFR_ORG","NOT_REFR_NUM","NULL",$2,$5}' $OutputF > "${OutputPath}net/Merge_plus_TU/tmp/${tmp_1}"

			# Giving a new format to the extended network
			tmp_2=$(printf "${file}_tmp_2")
			awk -F$"\t" 'BEGIN{OFS=FS}{print $0,"NOT_STRAND","NOT_TU_REFER"}' $InputExt > "${OutputPath}net/Merge_plus_TU/tmp/${tmp_2}"

			# Merging results and eliminating information not relevant (labels no longer used)
			tmp_3=$(printf "${file}_tmp_3")
			cat "${OutputPath}net/Merge_plus_TU/tmp/${tmp_2}" "${OutputPath}net/Merge_plus_TU/tmp/${tmp_1}" > "${OutputPath}net/Merge_plus_TU/tmp/${tmp_3}"
			bash "${SCRIPTPATH}/sqlite_group_tu_without_reference.sh" "${OutputPath}net/Merge_plus_TU/tmp/${tmp_3}" | sed -r 's/,NOT_REFR_ORG//g' | sed -r 's/,NOT_REFR_NUM//g' | sed -r 's/,NULL//g' | sed -r 's/NOT_STRAND,//g' | sed -r 's/NOT_TU_REFER,//g' | sed -r 's/,\w+_tu//g' > $OutputFinal

			rm "${OutputPath}net/Merge_plus_TU/tmp/${tmp_1}" "${OutputPath}net/Merge_plus_TU/tmp/${tmp_2}" "${OutputPath}net/Merge_plus_TU/tmp/${tmp_3}"

	        fi ) &

		((++count_batches)); [ "${count_batches}" -eq "${Batches}" ] && count_batches=1 && wait

	done

	wait



}

# Run current script if it "${OutputPath}net/Merge/tmp"was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../../utils/tracking.sh
        source ../../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        ExtendNetworksByOtho

fi
