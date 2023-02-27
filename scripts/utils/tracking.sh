#!bin/bash

# Tracking errors
TrackFailure() {
        
        local LINE_NUMBER=$1
        local ERROR_MESSAGE=$2
        
        echo "Error at ${lINE_NUMBER}: ${ERROR_MESSAGE}"

        exit 1
}