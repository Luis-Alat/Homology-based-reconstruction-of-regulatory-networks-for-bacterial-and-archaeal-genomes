#!bin/bash

################################################
#### Homology based-reconstruction pipeline ####
################################################

set -eE +o functrace

source ./utils/bash_messages.sh


# Stop execution and show on screen line number and bash command if there is any error

failure() {
        local LINE_NUMBER=$1
        local ERROR_MESSAGE=$2
        echo "Error at ${lINE_NUMBER}: ${ERROR_MESSAGE}"
}

trap ' failure ${LINENO} "$BASH_COMMAND" ' ERR

# Color code along the code

GREEN_COLOR=$'\e[1;32m'
RESET_COLOR=$'\033[0m'

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
            readarray -t GENOMES_FILE_PATH < $2
            shift 2
            ;;
        -n|--nets)
            readarray -t NETWORK_FILE_PATH < $2
            shift 2
            ;;
        -l|--labels)
            readarray -t LABELS < $2
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
# and for that reason is comment because is more like a separate step
#preprocessing

#proteinortho
#extendNet
#orthosumm
#extendNetPlusTU
#cytoScape
#coregulators
#hubs
#assignTypeTransc