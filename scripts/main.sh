#!/bin/bash

################################################
#### Homology based-reconstruction pipeline ####
################################################

set -eE +o functrace

# Importing bash_messages.sh [ShowHelp, ShowArguments], tracking.sh [TrackFailure] 
# Importing variables bash_messages GREEN_COLOR; RESET_COLOR
source ./utils/bash_messages.sh
source ./utils/tracking.sh

source ./01_acquisition_data/01_RunProteinortho_RunOperonMapper.sh

# Stop execution and show on screen line number and bash command if there is any error
trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

# Default arguments values

declare -A ARGUMENTS

GENOMES_FILE_PATH=""
NETWORK_FILE_PATH=""
LABELS=""

# Parsing argument values

if [[ $# -eq 0 ]]; then
    ShowHelp
fi

while [[ $# -gt 0 ]]; do

    ARGUMENTS[$1]=$2

    case "$1" in
        -g|--genomes)
            readarray -t GENOMES_FILE_PATH_VALUES < $2
            shift 2
            ;;
        -n|--nets)
            readarray -t NETWORK_FILE_PATH_VALUES < $2
            shift 2
            ;;
        -l|--labels)
            readarray -t LABELS_VALUES < $2
            shift 2
            ;;
        -h|--help)
            ShowHelp
            ;;
        *)
        echo "Invalid argument: $1"
        exit 1
        ;;
    esac

done

ShowArguments ARGUMENTS

# The preprocessing part is done manually for every organism used in this work
# and for that reason is commented and also more suitable to be ran separately
#preprocessing

RunProteinortho "../proteinortho/" GENOMES_FILE_PATH_VALUES LABELS_VALUES 2


#extendNet
#orthosumm
#extendNetPlusTU
#cytoScape
#coregulators
#hubs
#assignTypeTransc
