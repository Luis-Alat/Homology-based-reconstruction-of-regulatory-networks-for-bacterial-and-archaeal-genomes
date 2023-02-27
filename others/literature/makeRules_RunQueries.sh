## Make Ecoli queries/rules to search on pubmed

printf "Making E.coli queries file...\n"

# Get important columns from equiv tabs and make several RegEx commands
cut -f1,3 ../BAC9/others/equiv_tab/E_coli_RegulonDB_U00096.3 |
  sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
  # Replace NP by common name on the New interactions (predicted) with removed self-interactions
  sed -f - <(grep "New" ../../net/results_plus_TU/GCF_000005845.2_E_coli_K12_genomic_extended_network_plus_TU.txt |
    cut -f2,3 | perl -F"\t" -nae 'chomp($F[1]); if( $F[0] !~ /$F[1]/ ){print $_}') |
  grep -vE '(NP|YP)' |
  # Make queries
  sed -r 's/^/\(\(/g' | sed -r 's/\t/\) AND \(/g' | sed -r 's/$/\)\) AND ((E. coli) OR (Escherichia coli))/g' > E_coli_queries

printf "Making PubMed queries...\n"

Rscript --vanilla search_pubmed.R -i E_coli_queries -o E_coli_abstracts_new_search.txt
