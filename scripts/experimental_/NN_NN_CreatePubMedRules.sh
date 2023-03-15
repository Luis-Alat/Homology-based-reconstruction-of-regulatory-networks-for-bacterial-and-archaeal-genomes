#!/bin/bash

set -e

#bash NN_NN_CreatePubMedRules.sh ../../others/ids_equivalences/GCF_000005845.2_E_coli_K12_genomic.txt ../../network/predicted_nets/models/results_plus_TU/GCF_000005845.2_E_coli_K12_genomic.faa_extended_network_plus_tu "AND ((E. coli) OR (Escherichia coli))" > ../Ecoli_Rules_PubMed.txt

source ../utils/tracking.sh

trap ' TrackFailure ${LINENO} "$BASH_COMMAND" ' ERR

EQUIVALENCES_FILE=$1
NETWORK_FILE=$2
ORGANISM_STATEMENT=$3

GET_SED_RULES=$(cut -f1,3 $EQUIVALENCES_FILE |
                  sed -r 's/^/s\/\\</g; s/\t/\\>\//g; s/$/\/gI;/g; s/\./\\./g')

# All autoregulations will be removed
REPLACE_BY_RULES=$(sed "${GET_SED_RULES}" <(grep "New" $NETWORK_FILE | cut -f2,3 | perl -F"\t" -nae 'chomp($F[1]); if( $F[0] !~ /$F[1]/ ){ print $_ }' ) )

printf "%s\n" "${REPLACE_BY_RULES}" | 
  sed -r 's/^/\(\(/g' |
  sed -r 's/\t/\) AND \(/g' | 
  sed -r "s/$/\)\) /g" | 
  perl -slne 'chomp($_); print($_, $ORG)' -- -ORG="${ORGANISM_STATEMENT}";
