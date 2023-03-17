#!/bin/bash

# Tracking errors
TrackFailure() {
        
        source ./bash_messages.sh

        local LINE_NUMBER=$1
        local ERROR_MESSAGE=$2
        
        printf "${RED_COLOR}Error at ${lINE_NUMBER}: ${ERROR_MESSAGE}${RESET_COLOR}\n"

        exit 1
}