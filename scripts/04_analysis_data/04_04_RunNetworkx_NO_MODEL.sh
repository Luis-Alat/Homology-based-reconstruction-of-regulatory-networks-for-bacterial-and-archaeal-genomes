#!/bin/bash

set -e

AnalyzeByNetworkx() {

	# This function will run networkx to get general topology metrics of predicted networks for organism
	# without a network reference
	# This function is expecting 3 arguments. In order:
	# 1. String describing path to find all networks (passed by value)
	# 2. String describing path to place all results (passed by value)
	# 3. Integer of how many processes in background will be executed (passed by value)

	printf "  ${GREEN_COLOR}Running block AnalyzeByNetworkx${RESET_COLOR}\n\n"

	# Current script path
	local SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)
	local NETWORK_PATH=$(echo $1 | sed 's/\/*$/\//g')
	local OUPUT_PATH=$(echo $2 | sed 's/\/*$/\//g')
	
	[ -z $3 ] && local BATCHES=1 || local BATCHES=$3
	local COUNT_BATCHES=0

	[ ! -d $OUPUT_PATH ] && mkdir -p $OUPUT_PATH

	for NET in "${NETWORK_PATH[@]}"*; do
		
		if [ -f $NET ]; then

			NETWORK_BASENAME=$(basename $NET)
			OUTPUT_FILE_NAME=$(printf "${OUPUT_PATH}${NETWORK_BASENAME}_topology")

			python $SCRIPT_DIR/04_04_02_RunNetworkx_NO_MODEL.py --inputFile $NET --outputFile $OUTPUT_FILE_NAME &

			((++COUNT_BATCHES)); [ "${COUNT_BATCHES}" -eq "${BATCHES}" ] && COUNT_BATCHES=0 && wait

		fi
	done

	wait

}


# Run current script if it was called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

        source ../utils/tracking.sh
        source ../utils/bash_messages.sh

        # Stop execution and show on screen line number and bash command if there is any error
        trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

        # The following functions for now are expecting the values directly and not by command line
        AnalyzeByNetworkx

fi
