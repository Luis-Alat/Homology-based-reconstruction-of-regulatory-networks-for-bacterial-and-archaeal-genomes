#!/bin/bash

################################################
#### Homology based-reconstruction pipeline ####
################################################

#bash main.sh -g Fasta_files_path.txt -n Nets_files_path.txt -l Labels_organism.txt -t ../tus_models_operon_processed/ --extended_nets_output ../net/ --proteinortho_output ../proteinortho/ --tables_output ../analysis/tables/ --cytoscape_output ../analysis/cytoscape/ --coreg_output ../analysis/coreg/ --networkx_output ../analysis/hits/ --g_test_output ../analysis/gtest

set -eE +o functrace

# Importing bash_messages.sh [ShowHelp, ShowArguments], tracking.sh [TrackFailure] 
# Importing variables bash_messages GREEN_COLOR; RESET_COLOR
source ./utils/bash_messages.sh
source ./utils/tracking.sh

source ./01_acquisition_data/01_01_RunProteinortho_RunOperonMapper.sh
source ./03_processing_data/03_01_ExtendNetworks.sh
source ./04_analysis_data/04_01_RunCytoscape.sh
source ./04_analysis_data/04_02_TablesFromNets.sh
source ./04_analysis_data/04_03_GetCoregulates.sh
source ./04_analysis_data/04_04_RunNetworkx.sh
source ./04_analysis_data/04_05_RunGtest.sh
#source ./04_analysis_data/04_06_SearchPubMed.sh

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
CYTO_OUTPUT=""
TABLE_OUTPUT=""
COREG_OUTPUT=""
NETX_OUTPUT=""
G_TEST_OUTPUT=""

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
        --cytoscape_output)
            CYTO_OUTPUT=$2
            shift 2
            ;;
        --tables_output)
            TABLE_OUTPUT=$2
            shift 2
            ;;
        --coreg_output)
            COREG_OUTPUT=$2
            shift 2
            ;;
        --networkx_output)
            NETX_OUTPUT=$2
            shift 2
            ;;
        --g_test_output)
            G_TEST_OUTPUT=$2
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

# Deprecated in the future
#RunOperonMapper

ExtendNetworksByOtho $EXTENDED_NETWORKS_OUTPUT LABELS_VALUES $PROTEINORTHO_OUTPUT GENOMES_FILE_PATH_VALUES NETWORK_FILE_PATH_VALUES
ExtendNetworksByTranscriptionUnit $EXTENDED_NETWORKS_OUTPUT GENOMES_FILE_PATH_VALUES $TUS_PATH
AnalyzeByCytoscape GENOMES_FILE_PATH_VALUES $CYTO_OUTPUT "${EXTENDED_NETWORKS_OUTPUT}results" "${EXTENDED_NETWORKS_OUTPUT}results_plus_TU"
CreateNetworkTableMetrics LABELS_VALUES GENOMES_FILE_PATH_VALUES "${EXTENDED_NETWORKS_OUTPUT}results_plus_TU" $TABLE_OUTPUT
RunCoreg GENOMES_FILE_PATH_VALUES "${EXTENDED_NETWORKS_OUTPUT}results_plus_TU" $COREG_OUTPUT
RunGtest "${EXTENDED_NETWORKS_OUTPUT}results_plus_TU" $G_TEST_OUTPUT

# Literature
