# This script must be run inside of cytoscape directory
# Name to refer to the files

Prefix=("B_subtilis" "E_coli_K12" "P_aeruginosa_PA01" "S_enterica_LT2" "S_aureus_N315" "M_tuberculosis_H37Rv")
FolderNames=("Bacsu"  "Ecoli" "Pseudo" "Salmon" "Staph" "Myctu")

# Header of this new table
# Table describing how many interactions the new network was extended and where that interactions came from
printf ",B_subtilis,E_coli_K12,P_aeruginosa_PA01,S_enterica_LT2,S_aureus_N315,M_tuberculosis_H37Rv,Total de nuevas interacciones\n" > Table_percentage.txt

# Table describing how many interactions rebuild the modified network and where that interactions came from
printf ",B_subtilis,E_coli_K12,P_aeruginosa_PA01,S_enterica_LT2,S_aureus_N315,M_tuberculosis_H37Rv,Reconstrucción\n" > Table_percentage_reconstruction.txt

# Iterating over the number of files (6 files)
for ((i=0;i<${#Prefix[@]};i++)); do

        values_array_new=()
        values_array_rec=()

        new_interactions_number=$(grep "New" ../net/results_plus_TU/*"${Prefix[$i]}"*extended_network_plus_TU.txt | wc -l)
        reconstruction_number=$(grep "Known" ../net/results_plus_TU/*"${Prefix[$i]}"*extended_network_plus_TU.txt |
                grep -v "NOT_REFR_OR" | wc -l)

        # Iterating over every model organism (For example, Ecoli vs Bacsu, Ecoli vs Salmo, etc) ...
        # to get individually numbers
        for ((j=0;j<${#FolderNames[@]};j++)) do

                # Don't search againts itself. Add '---' inplace
                if [ $i -eq $j ]; then
                        value=$(echo -n "---,")
                        values_array_new+=$value
                        values_array_rec+=$value
                else

                        # Getting new interactions per organism and number of it
                        org_value=$(grep "New" ../net/results_plus_TU/*"${Prefix[$i]}"*extended_network_plus_TU.txt |
                                grep ${FolderNames[$j]} | wc -l)
                        # Getting percentage
                        percentage_value=$(echo "scale=2; ("$org_value"*"100")/"$new_interactions_number | bc)
                        # Giving format to save in array, useful later to save results into file
                        value=$(echo -n "${org_value} (${percentage_value}%),")
                        values_array_new+=$value

                        # Similar as above, but now with the Known interactions
                        org_value=$(grep "Known" ../net/results_plus_TU/*"${Prefix[$i]}"*extended_network_plus_TU.txt |
                                grep ${FolderNames[$j]} | wc -l)

                        percentage_value=$(echo "scale=2; ("$org_value"*"100")/"$reconstruction_number | bc)
                        
                        value=$(echo -n "${org_value} (${percentage_value}%),")
                        
                        values_array_rec+=$value

                fi

        done

        # Append in files created above
        echo "${Prefix[$i]},${values_array_new[*]}${new_interactions_number}" >> Table_percentage.txt
        echo "${Prefix[$i]},${values_array_rec[*]}${reconstruction_number}" >> Table_percentage_reconstruction.txt

done
