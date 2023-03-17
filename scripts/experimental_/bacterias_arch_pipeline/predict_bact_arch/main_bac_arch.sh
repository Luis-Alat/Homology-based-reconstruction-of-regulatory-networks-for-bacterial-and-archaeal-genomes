#!/bin/bash

#########################################################################
#### Homology based-reconstruction pipeline for bacteria and archael ####
#########################################################################

set -e

#bash main_bac_arch.sh -g GENOMES.txt -t TARGETS.txt -n NETS.txt -l LABELS.txt --proteinortho_output ../../../../proteinortho/bacteria --extended_nets_output ../../../../network/predicted_nets/bacteria--batches 20

# Current script path
SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

source $SCRIPT_DIR/../../../utils/bash_messages.sh
source $SCRIPT_DIR/../../../utils/tracking.sh
source $SCRIPT_DIR/01_RunProteinortho.sh
source $SCRIPT_DIR/02_ExtendNetworks.sh

trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

declare -A ARGUMENTS

BATCHES_NUMBER=1

# Parsing argument values

#if [[ $# -eq 0 ]]; then
#    ShowHelp
#fi

while [[ $# -gt 0 ]]; do

    ARGUMENTS[$1]=$2

    case "$1" in
        -g|--genomes_ref)
            readarray -t GENOMES_REF_FILE_PATH_VALUES < $2
            shift 2
            ;;
        -t|--genomes_tar)
            readarray -t GENOMES_TAR_FILE_PATH_VALUES < $2
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
        --networkx_output)
            NETX_OUTPUT=$2
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

# The following function has not been tested here. It was ran somewhere else
#RunProteinortho GENOMES_REF_FILE_PATH_VALUES LABELS_VALUES GENOMES_TAR_FILE_PATH_VALUES $PROTEINORTHO_OUTPUT

ExtendNetworksByOtho $EXTENDED_NETWORKS_OUTPUT LABELS_VALUES $PROTEINORTHO_OUTPUT GENOMES_TAR_FILE_PATH_VALUES NETWORK_FILE_PATH_VALUES $BATCHES_NUMBER     
#MakeExtentionTU
#AnalyzeNetx
