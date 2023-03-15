#!bin/bash

# Script to infer networks from aprox 12000 bacterias and 650 archaes
# This script was run in the parent folder

set -eE -o functrace

# Tracking errors

failure() {
        local lineno=$1
        local msg=$2
        echo "Error at $lineno: $msg"
}

#######################################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO RUN PROTEINORTHO AND PREDICT NETWORKS  ###############
#######################################################################################################

tool() {

    cd scripts/predict_bact_arch/
    # Predicting bacterias
    bash predict_network_without_reference.sh -d ../../bacterias_archaeans/bacteria_genomes/ -o ../../bacterias_archaeans/bacteria_nets_predicted/ -t ../../TUnits/bacteria/ -g genomes_old_names/
    # Predicting archae
    bash predict_network_without_reference.sh -d ../../bacterias_archaeans/archae_genomes/ -o ../../bacterias_archaeans/archae_nets_predicted/ -t ../../TUnits/archaea/ -g genomes_old_names/
    cd ../../

}


#####################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO GET NETWORK METRICS  ###############
#####################################################################################

NetMetrics() {

    cd  bacterias_archaeans/

    # Creating directories if they doesn't exist
    printf "Creating directories to save files if the don't exist..\n"

    if [ ! -d "archae_nets_predicted/metrics/" ]; then
        mkdir archae_nets_predicted/metrics/
    fi

    if [ ! -d "bacteria_nets_predicted/metrics/" ]; then
        mkdir bacteria_nets_predicted/metrics/
    fi

    # Main script of this function
    bash ../scripts/get_metrics_bact_archaes_pipeline.sh
    cd ../

}


## Running functions/main-program by order

tool
NetMetrics

