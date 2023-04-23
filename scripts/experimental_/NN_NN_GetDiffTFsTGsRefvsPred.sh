#!/bin/bash

set -e

GetValues() {

	local FILE_NAME=$1
	local COL_VALUE=$2
	local COL_NAME=$3

	local REF=$(grep "Known" "$FILE_NAME" | cut -f${COL_VALUE} | sort | uniq | wc -l)
	local EXT=$(cut -f${COL_VALUE} $FILE_NAME | sort | uniq | wc -l)
	local DIFF=$(echo "$(echo "${EXT} - ${REF}" | bc)")


	local OUTPUT_STRING=$(printf "${COL_NAME} REF:${REF} EXT:${EXT} DIFF:${DIFF}")

	printf "\t${OUTPUT_STRING}\n"

}

source ../utils/tracking.sh

trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

PATH_NETS="../../network/predicted_nets/models/results_plus_TU/"

for FILE in "${PATH_NETS}"*; do

	if [ -f $FILE ]; then

		echo $FILE

		GetValues $FILE "2" "TF"
		GetValues $FILE "3" "TG"

	fi

done
