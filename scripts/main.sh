#!/bin/bash

################################################
#### Homology based-reconstruction pipeline ####
################################################

#bash main.sh -g Fasta_files_path.txt -n Nets_files_path.txt -l Labels_organism.txt --extended_nets_output ../net --proteinortho_output ../proteinortho

set -eE +o functrace

# Importing bash_messages.sh [ShowHelp, ShowArguments], tracking.sh [TrackFailure] 
# Importing variables bash_messages GREEN_COLOR; RESET_COLOR
source ./utils/bash_messages.sh
source ./utils/tracking.sh

source ./01_acquisition_data/01_01_RunProteinortho_RunOperonMapper.sh
source ./03_processing_data/03_01_ExtendNetworks.sh

# Stop execution and show on screen line number and bash command if there is any error
trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

# Default arguments values

declare -A ARGUMENTS

GENOMES_FILE_PATH_VALUES=""
NETWORK_FILE_PATH_VALUES=""
LABELS_VALUES=""
PROTEINORTHO_OUTPUT="../proteinortho/"
BATCHES_NUMBER=1
EXTENDED_NETWORKS_OUTPUT="../net/"

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
        -t|--transcription_units_path)
            TUS_PATH=$2
            shift 2
        ;;
        --proteinortho_output)
            PROTEINORTHO_OUTPUT=$2
            shift 2
            ;;
        --batches)
            BATCHES_NUMBER=$2
            shift 2
            ;;
        --extended_nets_output)
            EXTENDED_NETWORKS_OUTPUT=$2
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
# for a deep view about that, check 02_preprocessing folder

RunProteinortho $PROTEINORTHO_OUTPUT GENOMES_FILE_PATH_VALUES LABELS_VALUES $BATCHES_NUMBER

# Still working on it
#RunOperonMapper

ExtendNetworksByOtho $EXTENDED_NETWORKS_OUTPUT LABELS_VALUES $PROTEINORTHO_OUTPUT GENOMES_FILE_PATH_VALUES NETWORK_FILE_PATH_VALUES
ExtendNetworksByTranscriptionUnit $EXTENDED_NETWORKS_OUTPUT GENOMES_FILE_PATH_VALUES $TUS_PATH
#cytoScape
#coregulators
#hubs
