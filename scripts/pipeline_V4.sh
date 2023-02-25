#!bin/bash

# This script was run inside of the parent folder
# This is the fourth version to analyse the 6 models organisms against each other
# Comments along the script are available describing the processing steps

# Tracking errors

set -eE -o functrace

failure() {
        local lineno=$1
        local msg=$2
        echo "Error at $lineno: $msg"
}

trap 'failure ${LINENO} "$BASH_COMMAND"' ERR

# Printing in green each block

greenprintf=$'\e[1;32m' # Green color
NC='\033[0m' # No Color

printf "${greenprintf}Setting variables...${NC}\n"

# Creating an array containing the genome names to run proteinortho in "mode" 1vs1 and save the results in the corresponding folder

GenomeNames=("GCF_000009045.1_ASM904v1_B_subtilis_168_genomic.faa" "GCF_000005845.2_E_coli_K12_genomic.faa" "GCF_000195955.2_ASM19595v2_M_tuberculosis_H37Rv_genomic.faa" "GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_genomic.faa" "GCF_000006945.2_ASM694v2_S_enterica_LT2_genomic.faa" "GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa")

# Folder basename  where outputs are going to be saved (but different parent directory)

FolderNames=("Bacsu"  "Ecoli"  "Myctu"  "Pseudo"  "Salmon"  "Staph")

# Base name describing the modified network

NetFiles=('GCF_000009045.1_ASM904v1_B_subtilis_168_modified_net.txt' 'GCF_000005845.2_E_coli_K12_modified_net.txt' 'GCF_000195955.2_ASM19595v2_M_tuberculosis_H37Rv_modified_net.txt' 'GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_modified_net.txt' 'GCF_000006945.2_ASM694v2_S_enterica_LT2_modified_net.txt' 'GCF_000009645.1_ASM964v1_S_aureus_N315_modified_net.txt')


###################################### Start functions ####################################


###########################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO PREPOCESSED THE INPUT FILES ##############
###########################################################################################

preprocessing() {

        if [ ! -d modifiedNetworks/ ]; then
                mkdir modifiedNetworks
        fi

        printf "${greenprintf}Getting genomes...${NC}\n"

        # Parsing gbff files for getting faa genome files. Defining headers of sequences as NCBI_ID,id,name,product,organism separated by '|'

        python scripts/ParserGBK.py --file genomes/gbk/NC_000913.3.gbk --org "Escherichia coli str. K-12 substr. MG1655" --get "protein_id,locus_tag,gene,product,translation" --sep "|" > genomes/GCF_000005845.2_E_coli_K12_genomic.faa
        python scripts/ParserGBK.py --file genomes/gbk/GCF_000009645.1_ASM964v1_genomic.gbff --org "Staphylococcus aureus subsp. aureus N315" --get "protein_id,old_locus_tag,gene,product,translation" --sep "|" > genomes/GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa
        python scripts/ParserGBK.py --file genomes/gbk/GCF_000195955.2_ASM19595v2_genomic.gbff --org "Mycobacterium tuberculosis H37Rv" --get "protein_id,locus_tag,gene,product,translation" --sep "|" > genomes/GCF_000195955.2_ASM19595v2_M_tuberculosis_H37Rv_genomic.faa
        python scripts/ParserGBK.py --file genomes/gbk/B_subtilis_168_GCF_000009045.1_ASM904v1_genomic.gbff --org "Bacillus subtilis subsp. subtilis str. 168" --get "protein_id,old_locus_tag,gene,product,translation" --sep "|" > genomes/GCF_000009045.1_ASM904v1_B_subtilis_168_genomic.faa
        python scripts/ParserGBK.py --file genomes/gbk/P_aeruginosa_GCF_000006765.1_ASM676v1_genomic.gbff --org "Pseudomonas aeruginosa PAO1" --get "protein_id,locus_tag,gene,product,translation" --sep "|" > genomes/GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_genomic.faa
        python scripts/ParserGBK.py --file genomes/gbk/salmonella_enterica_serovar_LT2_RN_GCF_000006945.2_ASM694v2_genomic.gbff --org "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2" --get "protein_id,locus_tag,gene,product,translation" --sep "|" > genomes/GCF_000006945.2_ASM694v2_S_enterica_LT2_genomic.faa

        # Removing repeated sequences in S_aureus_n315 (Manually identfied and removed). Overwritting  the original file

        cat <(grep -A 1 -m 1 "WP_010959190.1" genomes/GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa) <(sed -r '/WP_010959190.1/ { N; d; }' genomes/GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa) > genomes/tmp && mv genomes/tmp genomes/GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa

        printf "${greenprintf}Getting modified networks...${NC}\n"

        # Getting the modified networks (cleaning and parsing)

        # Getting the modified network of E.coli
        ## Getting from the faa file NCBI_id and name to match and replace the data in the original network 
        grep ">" genomes/GCF_000005845.2_E_coli_K12_genomic.faa | cut -d"|" -f1,3 |
                ## Replacing '|' by tab and removing >. Resorting columns (id,NCBI_ID)
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                ## Creating multiples regex command to replace id in the original network by the chosen one (NCBI_ID)
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                ## Replacing ids in the original network
                sed -f - network/ecoli_tf_network.txt |
                ## Removing comments in the network. Retrieving only rows where a replaced were done
                grep -v "#" | grep -E "(^NP_|^YP_)" | grep -E "(\s+\S*NP_|\s+\S*YP_)" |
                ## Numerating row (interactions). Save file as modified
                nl | sed -r 's/^\s+//g' > modifiedNetworks/GCF_000005845.2_E_coli_K12_modified_net.txt

        # Getting the modified network of M.tuberculosis, but firts...
        ## Retrieving sigma factors as regex commands. Those elements will be removed from the network
        cut -f1 others/Sigma_factors_excluir/mtuberculosis_excluir | sed -r 's/#//g' | sed -r 's/$/\\s\+Rv/g' |
                grep -Evif - network/mtuberculosis_tf_network.txt > modifiedNetworks/tmp

        ## Pretty much as the steps above with Ecoli, but note some changes as getting NCBI_ID and ID instead of name due to the original network uses different elements...
        ## to describe the interactions (Ecoli uses common name, Mtuberculosis uses Rv number)
        grep ">" genomes/GCF_000195955.2_ASM19595v2_M_tuberculosis_H37Rv_genomic.faa | cut -d"|" -f1,2 |
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                ## Replacing ids on the network without sigma factors, not the original one. Removing empty lines, skipping first column (not important)..
                ## and replacing <spaces> by <tab>
                sed -f - <(grep -v "#" modifiedNetworks/tmp | grep -Ev "^$" | cut -f2- | sed -r 's/\s+/\t/g') |
                grep -E "(^NP_|^YP_)" | grep -E "(\s+\S*NP_|\s+\S*YP_)" |
                nl | sed -r 's/^\s+//g' > modifiedNetworks/GCF_000195955.2_ASM19595v2_M_tuberculosis_H37Rv_modified_net.txt

        # Getting the modified network of B.subtillis, but first (again)...
        ## Retrieving sigma factors to remove from the original network similar with Mtuberculosis
        cut -f1 others/Sigma_factors_excluir/bsubt_excluir | sed -r 's/#//g' | sed -r 's/$/\t\(yes\|no\)/g' |
                grep -Evwif - network/bacsu_tf_network.txt > modifiedNetworks/tmp

        # Bsubtillis network uses both a common name as a id (what a nightmare!). Here, the replacement is based on the id and changed by the NCBI_ID
        grep -v "#" modifiedNetworks/tmp | grep -v "Gene_Name" | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                sed -f - modifiedNetworks/tmp | grep -v "#" | grep -v "Gene_Name" | cut -f2- |
                awk -F"\t" 'BEGIN{OFS=FS}{print$2,$1,$3,$4,$5,$6,$7,$8}' > modifiedNetworks/tmp_1

        # Here we use the common name to replace it with the id (More interactions were mapped using this method)
        grep ">" genomes/GCF_000009045.1_ASM904v1_B_subtilis_168_genomic.faa | cut -d"|" -f2,3 |
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                sed -f - modifiedNetworks/tmp_1 > modifiedNetworks/tmp_2

        # Almost done. Here, we replace all the ids by the NCBI_id
                grep ">" genomes/GCF_000009045.1_ASM904v1_B_subtilis_168_genomic.faa | cut -d"|" -f1,2 |
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                sed -f - modifiedNetworks/tmp_2 |
                perl -F"\t" -nae 'if(($F[0] =~ /NP_/) or ($F[0] =~ /YP_/)){if(($F[1] =~ /NP_/) or ($F[1] =~ /YP_/)){print $_}}' > modifiedNetworks/tmp_3

        # Finally, elements described in the network like ribozyme or RNA are removed to give format and saving the modified network
        sed -r 's/( ribozyme|\*)//g' modifiedNetworks/tmp_3 | sed -r 's/(RNA - | riboswitch|as-)//g' |
                nl | sed -r 's/^\s+//g' > modifiedNetworks/GCF_000009045.1_ASM904v1_B_subtilis_168_modified_net.txt

        # Getting the modified network of P.aureus (process very similar like Bsubtillis)
        # Here, a new network is created to remove empty columns and, in general, anything that could be a problem later 
        grep -v "#" network/pseudomonas_tf_network.txt | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                sed -f - <(grep -v "#" network/pseudomonas_tf_network.txt | cut -f2-) |
                awk -F"\t" '{OFS=FS}{print $1,$3,$4,$2,$5,$6}' > modifiedNetworks/tmp

        # Replacing common name by id
        grep ">" genomes/GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_genomic.faa | cut -d"|" -f2,3 |
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                sed -f - modifiedNetworks/tmp > modifiedNetworks/tmp_2

        # Replacing id by NCBI_ID
        grep ">" genomes/GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_genomic.faa | cut -d"|" -f1,2 |
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                sed -f - modifiedNetworks/tmp_2 |
                perl -F"\t" -nae 'if(($F[0] =~ /NP_/) or ($F[0] =~ /YP_/)){if(($F[1] =~ /NP_/) or ($F[1] =~ /YP_/)){print $_}}' > modifiedNetworks/tmp_3

        # Tuning details and creating the modified network
        sed -r 's/\.11//g' modifiedNetworks/tmp_3 |
                sed -r 's/\.1\.1/\.1/g' |
                nl | sed -r 's/^\s+//g' > modifiedNetworks/GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_modified_net.txt

        # Getting the modified network of S.enterica (easier than the networks before)
        grep ">" genomes/GCF_000006945.2_ASM694v2_S_enterica_LT2_genomic.faa | cut -d"|" -f1,2 |
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}'|
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g'|
                sed -f - <(sed -r 's/,/\t/g' network/salmonella_enterica_serovar_Typhimurium_LT2_RN.csv |
                awk -F$"\t" 'BEGIN{OFS=FS}{print $1,$4,$3,$6}' | grep -v "Node_1_locus_tag") | sed -r 's/(\.S|\.1N|\.c)//g' |
                nl | sed -r 's/^\s+//g'  > modifiedNetworks/GCF_000006945.2_ASM694v2_S_enterica_LT2_modified_net.txt

        # Getting the modified network of S.aureus N315 (similar processing, removing sigma, removing inncesary elements, mapping etc..)
        cut -f1 others/Sigma_factors_excluir/staphy_excluir | sed -r 's/#//g' | grep -Fvwif - network/staphy_tf_network.txt > modifiedNetworks/tmp

        grep -v "#" modifiedNetworks/tmp | grep -vE "(Genomic Island|Phage region|T-box|_leader|nitrate reductase)" | cut -f2- > modifiedNetworks/tmp_1

        # Creating one unique file described by id as name to do a mapping of both
        cat <(cut -f1,4,6,7 modifiedNetworks/tmp_1) <(cut -f1,5,6,7 modifiedNetworks/tmp_1) <(cut -f 2,4,6,7 modifiedNetworks/tmp_1) <(cut -f2,5,6,7 modifiedNetworks/tmp_1) |
                sort | grep -Ev "^\s+" > modifiedNetworks/tmp_2

        grep ">" genomes/GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa | cut -d"|" -f1,3 |
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                sed -f - modifiedNetworks/tmp_2 > modifiedNetworks/tmp_3

        grep ">" genomes/GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa | cut -d"|" -f1,2 |
                sed -r "s/\|/\t/g" | sed -r 's/>//g' | awk -F$"\t" 'BEGIN{OFS=FS}{print $2,$1}' |
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                sed -f - modifiedNetworks/tmp_3 > modifiedNetworks/tmp_4

        grep "WP_" modifiedNetworks/tmp_4 | perl -F"\t" -nae 'if(($F[0] =~ /WP_/) and ($F[1] =~ /WP_/)){print $_}' |
                sort | uniq | nl |
                sed -r 's/^\s+//g' | sed -r 's/\.1,/\.1/g' > modifiedNetworks/GCF_000009645.1_ASM964v1_S_aureus_N315_modified_net.txt

        # All the tmp files are removed. We dont need them anymore! c:
        rm modifiedNetworks/tmp_4 modifiedNetworks/tmp_3 modifiedNetworks/tmp_2 modifiedNetworks/tmp_1 modifiedNetworks/tmp

}

################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO RUN PROTEINORTHO ##############
################################################################################

# In order to run proteinortho, we need to modify the headers of the fasta file in order to take ...
# from the headers the id given by NCBI and one additional id. This ensure the a standart format in the proteinortho output ...
# files

proteinortho() {

        printf "${greenprintf}Porteinortho block${NC}\n"

        for faa in ${GenomeNames[@]}; do
                # Changing the second '|' by <tab> inplace
                printf "Modifying faa file "$faa"\n"
                sed -r -i 's/\|/\t/2' genomes/$faa

        done

        # Running proteinortho to get orthologs
        ## Verifyng if the output directory doesn't exist. Create it if doesn't
        if [ ! -d proteinortho/ ];then
                mkdir proteinortho
        fi

        # Proteinortho was run like paired comparations. Every comparation has its own folder to save the results of the orthologous
        for ((i=0;i<${#FolderNames[@]};i++)) ; do

                ## Creating folder by each organism comparation
                if [ ! -d proteinortho/${FolderNames[$i]} ];then
                        mkdir proteinortho/${FolderNames[$i]}
                fi

                ## Apparently the proteinortho outputs only can be saved in the directory where it's ran (even if it has a parameter to output files)
                cd proteinortho/${FolderNames[$i]}

                for ((j=0;j<${#GenomeNames[@]};j++));do

                        # Do not search orthologues against itself (dah)
                        if [[ $i -eq $j ]]; then
                                continue
                        fi

                        ## Defining a name for the proyect
                        ProjectName=$(echo ${GenomeNames[$j]}"_"${FolderNames[$i]})
                        ## Running proteinortho
                        proteinortho6.pl --verbose=2 -cpus=5 --cov=70 -project=$ProjectName ../../genomes/${GenomeNames[$i]} ../../genomes/${GenomeNames[$j]}

                done

                ## Returning to the parent folder
                cd ../../

        done

}

############################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO CREATE THE EXTENDED NETWORKS ##############
############################################################################################

extendNet() {

        printf "${greenprintf}Extended network block${NC}\n"

        # In this "for" every ortolog is retrieved and a new file is created. That new file is used to create another one wich is useful..
        # to get the new interactions. The final file saved in PreNewInteractionsFiles will describe what interactions were mapped ...
        # among what organisms. For example, the folder and file net/Ecoli/*P_aeruginosa_PA01* contain based on the orthologue relation ...
        # what interactions were mapped from the network of Ecoli to P_aeruginosa_PA01

        # All the outputs are going to be saved inside of "net" directory in the respective organism folder

        # Creating folder if it doesn't exist
        if [ ! -d net/ ];then
                mkdir net
        fi

        for ((i=0; i<${#FolderNames[@]}; i++)); do

                # Creating folder if it doesn't exist (Remember a paired comparations were done, every organism is going to have its own folder! :o) ...
                # in this case, we need two, one to save the relations among orthologues, and the other one to save a draft (?) of the interactions predicted
                if [ ! -d net/${FolderNames[$i]} ];then
                        mkdir net/${FolderNames[$i]} net/${FolderNames[$i]}/ortoFiles/ net/${FolderNames[$i]}/PreNewInteractionsFiles/
                fi

                FolderFiles=$(echo "proteinortho/"${FolderNames[$i]}"/*blast-graph")

                # Making a tmp file to remove interactions at the end of the following for loop...
                # (some interactions are described by multiples TFs, the reason of this file is to remove those interactions ...
                # if one of the TFs was not found in the target organism at the moment of doing an extrapolation of the interactions among organisms)
                grep ">" genomes/${GenomeNames[$i]} | cut -d"|" -f1 | sed -e 's/>//g' > net/${FolderNames[$i]}/PreNewInteractionsFiles/tmp

                # Iterating in every "blast-graph" file found inside of ${FolderNames[@]}
                for filePath in $FolderFiles; do

                        # Defining new file names based in a sufix of the fasta files 
                        FileName=$(basename $filePath | sed -r 's/\.faa.*//g')
                        OutFileOrtho=$(echo $FileName'_ortho.txt')
                        OutFileNewInteractions=$(echo $FileName'_new_interactions.txt')

                        printf "Making files from "$filePath"\n"

                        # Given a specific file inside of the folder, search for a pattern and get the column number...
                        # Later it will be useful to reorder the file in the "correct" order (giving structure) of organism model...
                        # (same name as parent folder) and then target organism

                        ## First 5 lines contain the necessary information, then search by using the name file
                        ## our interes column and print the right order separated by comma
                        string=$(head -n 5 $filePath | awk -F$"\t" -v pattern=$FileName '$1 ~ pattern {print 2","1} $2 ~ pattern {print 1","2}')
        
                        ## The output above now will be an array
                        IFS="," read -r -a columnOrder <<< "$string"

                        # Reording the columns
                        ## Pasting in this order the column described in ${columnOrder[0]} (model organism) and later ${columnOrder[1]} (model target) ...
                        ## also removing the headers/rows written by proteinortho
                        paste <(grep -v "#" $filePath | cut -f${columnOrder[0]}) <(grep -v "#" $filePath | cut -f${columnOrder[1]}) |
                                ## Replacing "|" by TAB. Saving only the columns importants (NCBI_ID, not the id)
                                sed -r 's/\|/\t/g' | cut -f1,3 > net/${FolderNames[$i]}/ortoFiles/$OutFileOrtho

                        # Given the new file above, we are going to use it to identify interactions
                        ## Creating multiples regex commands to replace in the network of the model organism (mapping orthologues) ...
                        ## by the NCBI_id of the target organism. Every regex command is going to have a format NCBI_ID_model,NCBI_ID_targetMATCH ...
                        ## The keyword MATCH is going to be useful to retrieve what interactions were mapped
                        sed -r 's/^/s\/\\</g' net/${FolderNames[$i]}/ortoFiles/$OutFileOrtho | sed -r 's/\t/\\>\//g' | sed -r 's/$/MATCH\/gI;/g' | sed -r 's/\./\\./g' |
                                ## Replacing NCBI_ID_model by the NCBI_ID_targetMATCH
                                sed -f - <(cut -f1-3 modifiedNetworks/${NetFiles[$i]}) |
                                ## Identifying what interactions were replaced, retrieving those and removing keyword no longer used
                                perl -F"\t" -nae 'if(($F[1] =~ /MATCH/) and ($F[2] =~ /MATCH/)){print $_}' | sed -r 's/MATCH//g' |
                                ## Retrieving the organism name where that new interaction came from
                                awk -F $'\t' -v organism=${FolderNames[$i]} 'BEGIN{OFS = FS} {print organism,$0}' |
                                ## Filtering those interactions that contain valid NCBI_IDs
                                grep -Fvwif net/${FolderNames[$i]}/PreNewInteractionsFiles/tmp - > net/${FolderNames[$i]}/PreNewInteractionsFiles/$OutFileNewInteractions

                        done
        done

        # Creating array of file names, It's useful to iterate
        mapfile -t FileNewInterNames < <(find proteinortho/*/*blast-graph -exec basename {} \; | sed -r 's/\.faa.*//g' | sort | uniq)

        # Creating folder if it doesn't exist
        if [ ! -d net/Merge/ ];then
                mkdir net/Merge
        fi

        # Merging results ("for" loop below) of every Folder (${FolderNames[@]}) using the name base founds inside of "proteinortho/${FolderNames[@]}/"
        # Every final file generated in the for loop above contains ...
        # organism_name_where_interaction_came_from, organism_interaction_number_where_interaction_came_from, TF, TG by model and target organism spread ...
        # among different folders. This for loop merge the respective files

        # Results will be saved in the folder "net/Merge"
        for file in ${FileNewInterNames[@]}; do

                # Defining new file name and making empty file
                MergeFileName=$(echo $file"_Merged.txt")
                echo -n > net/Merge/tmp

                # Merging files into the "Merge" folder
                for folder in ${FolderNames[@]}; do

                        printf $file"\t"$folder"\n"
                        InputFile=$(echo "net/"$folder"/PreNewInteractionsFiles"/$file"*")

                        # Check if file exist in that folder
                        if [ -f $InputFile ]; then
                                # Merging all results into a temporaly file
                                cat $InputFile >> net/Merge/tmp
                        fi

                done

                printf "Making file "$MergeFileName"\n"

                ## This line calls an sql implementation to group the columns associated with regulatory interactions to describe
                ##(1) transcription factor,
                ##(2) target,
                ##(3) organism where that interaction was found and
                ##(4), line number associated with that interaction found in modified networks...
                ##In addition adds a number indicating in how many organism that interactions was found
                bash scripts/sqlite_join.sh net/Merge/tmp  |
                        perl -F"\t" -nae '
                                chomp($_);
                                @array=split(",",$F[2]);
                                %seen=();
                                @unique = grep { ! $seen{$_} ++ } @array;
                                $len = scalar(@unique);
                                print ($_,"\t",$len,"\n")' > net/Merge/$MergeFileName

                rm net/Merge/tmp
        done

        # In the following "for" loop an "extended network" is going to be created for every model organism...
        # It has the interactions described in the modified network as well as the new interactions found
        # Results will be saved in the folder "net/results"

        # Creating folder if it doesn't exist
        if [ ! -d net/results ];then
                mkdir net/results net/results/tmp
        fi

        for ((i=0;i<${#GenomeNames[@]};i++)); do

                # Getting base name
                BaseName=$(echo ${GenomeNames[$i]} | sed -r 's/\.faa.*//g')
                OutPutFileTmp=$(echo $BaseName"TMP")

                printf "Creating result files of "${GenomeNames[$i]}" using the modified network "${NetFiles[$i]}"\n"

                # Giving format both $MergeFile and the modified network and later they are concatenated...
                # Any missing info is added too
                ## Retrieving TF and TG columns, adding info about the state of the interaction (does the interactions is known?)
                cut -f1-3 modifiedNetworks/${NetFiles[$i]} |
                        awk -F$"\t" 'BEGIN{OFS=FS}{print $0,"Known","NOT_REFR_ORG","NOT_REFR_NUM","NULL"}' > net/results/tmp/$OutPutFileTmp

                # Adding info about if the interaction is new and a new enumerate id
                awk -F$'\t' 'BEGIN{OFS=FS} $2=$2 FS "New?"' net/Merge/$BaseName* |
                        perl -nse 'print($count,"_new?","\t",$_); $count+=1' -- -count=1 >> net/results/tmp/$OutPutFileTmp

                # Grouping (group_by sql) results and removing "innecesary" information (labels no longer used)
                OutPutFile=$(echo $BaseName"_extended_network.txt")

                bash scripts/sqlite_groupe_results.sh net/results/tmp/$OutPutFileTmp |
                        ## Removing innecesary labels, for example, an interaccion "Known" could be also described by the mapping approach ...
                        ## therefore, that interaction will have the label "Known" as "New", but we already know thats not new. "New" is deleted
                        sed -r 's/(,\w+_new\?)+//1' | sed -r 's/_new\?/_new/1' |
                        sed -r 's/,New\?//g' | sed -r 's/\tNew\?\t/\tNew\t/g' |
                        sed 's/NOT_REFR_ORG,//g' | sed 's/NOT_REFR_NUM,//g' |
                        sed -r 's/(NULL,)+//g' | sed -r 's/(Known,)+//g' | sort -k2,3 > net/results/$OutPutFile

        done

}

############################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO CREATE A SUMMARY OF THE OTHOLOGS ##########
############################################################################################

orthosumm() {

        printf "${greenprintf}Summary of orthologs block${NC}\n"

        # It useful to assign an output name and also to find files the following array
        mapfile -t BaseName < <(printf "%s\n" ${GenomeNames[@]%.faa*})

        # The following for will create a file by each organism describing the orthologues...
        # and the equivalences found by proteinortho among one organism vs all the other organism and ...
        # a number that will describe in how many organism that gene has orthologues  

        for ((i=0;i<${#FolderNames[@]};i++)); do

                # Creating and empty temporal file
                echo -n > net/results/tmp/tmp

                for ((j=0;j<${#GenomeNames[@]};j++)); do

                        # Do not search orthologues against itself
                        if [[ $i -eq $j ]]; then
                                continue
                        fi

                        # Given a specific folder, add a column specifying where the orthologue came from
                        printf "Retriving ortologs of "${FolderNames[$i]}" using "${BaseName[$j]}"\n"
                        awk -F $'\t' -v organism=${FolderNames[$j]} 'BEGIN{OFS = FS} {print  $0,organism}' net/${FolderNames[$i]}/ortoFiles/${BaseName[$j]}* >> net/results/tmp/tmp

                done

                # Adding information from modifed network, grouping results, removing labels no longer used and adding a new line about where the ortholog came from
                # Concept very similar seen above at the moment of retrieving the new interactions
                OutPutName=$(echo ${BaseName[$i]}"_Orto_Summary.txt")
                cut -f 2,3 modifiedNetworks/${NetFiles[$i]} |
                        sed -r 's/\t/\n/g' | sort | uniq |
                        awk -F$"\t" 'BEGIN{OFS=FS} {print $0, "NULL", "NULL"}' >> net/results/tmp/tmp

                bash scripts/sqlite_group_ortho.sh net/results/tmp/tmp |
                        sed -r 's/(,NULL)+//g' |
                        perl -F"\t" -nae '
                                chomp($_);
                                chomp($F[2]);
                                if($F[2] !~ /NULL/){
                                        @array=split(",",$F[2]);
                                        %seen=();
                                        @unique = grep { ! $seen{$_} ++ } @array;
                                        $len = scalar(@unique);
                                        print ($_,"\t",$len,"\n")
                                } else {
                                        print ($_,"\t0\n")
                                }' > net/results/$OutPutName

                rm net/results/tmp/tmp

        done
        
}

######################################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO PROCESS TRANSCRIPTION UNIT (TU) FILES ###############
######################################################################################################

extendNetPlusTU() {

        # Creating a bash array and get ID names
        mapfile -t FileNameID < <(find net/results/*extended_network.txt -exec basename {} \; | sed -r 's/\..*//g')

        printf "${greenprintf}Transcription units block${NC}\n"

        # If folder doesn't exit, create it

        if [ ! -d net/results_plus_TU/ ]; then
                mkdir net/results_plus_TU/ net/results_plus_TU/tmp/
        fi

        # Using FileNameID array

        for file in ${FileNameID[@]}; do

                # Verifying if the file exists (Does exist the TU file for this particular model organism?)
                if [ -f TUnits/bacteria/$file* ]; then

                        # Defining an output final file name for the python script below
                        OutPutName=$(find net/results/$file*extended_network.txt -exec basename {} \; |
                                sed -r 's/\.txt/_plus_TU\.txt/g')

                        # Searching TUs where the first gene is equal to the target in the extended network
                        python scripts/ParserTUandExt.py --inputNet net/results/$file*extended_network.txt --inputTus TUnits/bacteria/$file* --outputPath net/results_plus_TU/tmp/

                        printf "Making net/results_plus_TU/$OutPutName\n"

                        # Giving a new format to the earlier output to match with the extended network format. Saving in a tmp_1 file
                        awk -F$"\t" 'BEGIN{OFS=FS}{print $4,$2,$3,"New?","NOT_REFR_ORG","NOT_REFR_NUM","NULL",$1,$4}' net/results_plus_TU/tmp/$file* > net/results_plus_TU/tmp_1

                        # Giving a new format to the extended network to match with the output python file above. Saving in a tmp_2 file 
                        awk -F$"\t" 'BEGIN{OFS=FS}{print $0,"NOT_STRAND","NOT_TU_REFER"}' net/results/$file*extended_network.txt > net/results_plus_TU/tmp_2

                        # Merging results of both files
                        cat net/results_plus_TU/tmp_2 net/results_plus_TU/tmp_1 > net/results_plus_TU/tmp_3

                        # Removing information not relevant (labels no longer used) useful during the group_by, but not later
                        bash scripts/sqlite_group_TU.sh net/results_plus_TU/tmp_3 |
                                sed -r 's/,New\?//g' | sed -r 's/,NOT_REFR_ORG//g' |
                                sed -r 's/,NOT_REFR_NUM//g' | sed -r 's/,NULL//g' |
                                sed -r 's/NOT_STRAND,//g' | sed -r 's/NOT_TU_REFER,//g' |
                                sed -r 's/,\w+_tu//g' | sed -r 's/New\?/New/g' > net/results_plus_TU/$OutPutName
                fi

        done

        # Removing tmp files generated above
        rm net/results_plus_TU/tmp_1 net/results_plus_TU/tmp_2 net/results_plus_TU/tmp_3

}

#########################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO RUN CYTOSCAPE ##########
#########################################################################

cytoScape() {

        printf "${greenprintf}CytoScape block${NC}\n"

        # It useful to assign an output name and also to find files the following array
        mapfile -t BaseName < <(printf "%s\n" ${GenomeNames[@]%.faa*})

        # Creating folder if it doesn't exist
        if [ ! -d cytoscape/ ];then
                mkdir cytoscape/ cytoscape/tmp/
        fi

        # All the following output files are going to be saved in the current path
        # User must lauch CytoScape manually
        # If exist previous files, CytoScape will ask if you want to overwrite them (even if you modify the R API to overwrite)...
        # for simplicity the existent file could be erased to avoid this
        # Six files are going to be generated by every model organism: One a .cys file creating a network with...
        # the relation among orthologues where nodes are colored based on how many orthologues of each gene was found in the others...
        # organisms, a second and third .cys file with the extended network and extended network plus TUs info...
        # (new interactions colored in red) and a fourth, fith and sixth one describing general topological metrics of each network

        cd cytoscape/
        path=$(realpath . | sed -r 's/(\/$|$)/\//g')

        for ((i=0;i<${#GenomeNames[@]};i++)); do

                # Using the baseName (type array variable) previously defined
                newInputFile=$(echo ${BaseName[$i]}"_extended_network.txt")
                orthoInputFile=$(echo ${BaseName[$i]}"_Orto_Summary.txt")
                newInputFileTU=$(echo ${BaseName[$i]}"_extended_network_plus_TU.txt")

                TMPnewInputFile=$newInputFile
                TMPnewInputFileTU=$newInputFileTU

                # Calling the R script to execute commands on CytoScape (API), generating and analyzing networks
                Rscript --vanilla ../scripts/cytoScapeBuildNetworks.R -i ../net/results_plus_TU/$newInputFileTU -t new -p $path
                Rscript --vanilla ../scripts/cytoScapeBuildNetworks.R -i ../net/results/$newInputFile -t new -p $path
                Rscript --vanilla ../scripts/cytoScapeBuildNetworks.R -i ../net/results/$orthoInputFile -m ../modifiedNetworks/${NetFiles[$i]} -t ortho -c "#89D0F5,#9D6AC7,#97FA4A,#EFFA4A,#FAAD4A,#FA554A" -p $path


                # Retrieving output file names
                newInputFile=$(echo ${BaseName[$i]}"_extended_networkStatistics.txt")
                orthoInputFile=$(echo ${BaseName[$i]}"_Orto_SummaryStatistics.txt")
                newInputFileTU=$(echo ${BaseName[$i]}"_extended_network_plus_TUStatistics.txt")

                # Adding additional (inplace) information in file of statistic...
                # (cut -f2 gets the number of TFs in that network, cut -f3 gets the target)...
                # from modified network and the extended network

                TFCounts=$(cut -f2 ../modifiedNetworks/${NetFiles[$i]} | sort | uniq | wc -l | sed -r 's/^/#TF_counts\t/g')
                sed -i "2i\\$TFCounts" $orthoInputFile
                TargetCounts=$(cut -f3 ../modifiedNetworks/${NetFiles[$i]} | sort | uniq | wc -l | sed -r 's/^/#Target_counts\t/g')
                sed -i "2i\\$TargetCounts" $orthoInputFile

                TFCounts=$(cut -f2 ../net/results/$TMPnewInputFile | sort | uniq | wc -l | sed -r 's/^/#TF_counts\t/g')
                sed -i "2i\\$TFCounts" $newInputFile
                TargetCounts=$(cut -f3 ../net/results/$TMPnewInputFile | sort | uniq | wc -l | sed -r 's/^/#Target_counts\t/g')
                sed -i "2i\\$TargetCounts" $newInputFile

                TFCounts=$(cut -f2 ../net/results_plus_TU/$TMPnewInputFileTU | sort | uniq | wc -l | sed -r 's/^/#TF_counts\t/g')
                sed -i "2i\\$TFCounts" $newInputFileTU
                TargetCounts=$(cut -f3 ../net/results_plus_TU/$TMPnewInputFileTU | sort | uniq | wc -l | sed -r 's/^/#Target_counts\t/g')
                sed -i "2i\\$TargetCounts" $newInputFileTU

        done

        # The following script must be run inside the cytoscape folder. Copy that file there and run it as below
        if [ ! -f makeGeneralTables.sh ]; then
                cp ../scripts/makeGeneralTables_extendedTU_nets.sh .
        fi

        # Based on the extended network with info from Tus. Create two general tables about ...
        # how many interactions the new network was extended and where that interactions came from and...
        # how many interactions rebuild the modified network and where that interactions came from
        printf "\n\nMaking tables of results...\n"
        bash makeGeneralTables_extendedTU_nets.sh
        printf "Done\n\n"

        cd ../

}

##################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO GET CO-REGULATORS ###############
##################################################################################

coregulators() {

        # All the following output files are going to be saved in the current path and within of a tmp folder
        printf "${greenprintf}Co-regulatos block${NC}\n"

        mapfile -t FileNameID < <(find net/results/*extended_network.txt -exec basename {} \; | sed -r 's/\..*//g')

        # Creating folder if it doesn't exist
        if [ ! -d co_regulators/ ]; then
                mkdir co_regulators/ co_regulators/tmp/ co_regulators/modules_comparation/
        fi

        cd co_regulators/

        # Using FileNameID again
        # Using CoReg to find co-regulators in the modified network, extended network and extended with Tus info for
        # each organism model
        for file in ${FileNameID[@]}; do
                Rscript --vanilla ../scripts/CoReg.R -i ../net/results/$file*_extended_network*
                Rscript --vanilla ../scripts/CoReg.R -i ../net/results_plus_TU/$file*
                Rscript --vanilla ../scripts/CoReg.R -i ../modifiedNetworks/$file*
        done

        cd modules_comparation/
        # The following script must be run in t
        if [ ! -f compare_modules.sh ]; then
                cp ../../scripts/compare_modules.sh .
        fi

        # Comparing modules among organism to get an overview of those. It could take a while...
        # Every comparation is made for the modified network, extended network and extended plus Tu info; therefore...
        # the following scrip will make 2 files for every network for every organism, in other words, 6 files for each...
        # organism saved in co_regulators/modules_comparation path
        bash compare_modules.sh

        cd ../../
}

##########################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO GET HUBS  ###############
##########################################################################

hubs() {

        printf "${greenprintf}Hub scores block${NC}\n"

        # Using again FileNameID array
        mapfile -t FileNameID < <(find net/results/*extended_network.txt -exec basename {} \; | sed -r 's/\..*//g')

        for file in ${FileNameID[@]}; do

                # Getting hubs from extended networks
                OutPut=$(find net/results/$file*_extended_network* -exec basename {} \; | sed -r 's/_extended.*/_extended_network_hubs_scores.txt/g' )
                printf "Making $OutPut\n"
                python scripts/getHubs.py net/results/$file*_extended_network* "2,3" | sort -n -k2,2 > cytoscape/$OutPut

                # Getting hubs from modified networks
                OutPut=$(find modifiedNetworks/$file* -exec basename {} \; | sed -r 's/\.txt/_hubs_scores.txt/g' )
                printf "Making $OutPut\n"
                python scripts/getHubs.py modifiedNetworks/$file* "2,3" | sort -n -k2,2 > cytoscape/$OutPut

                # Getting hubs from extended networks with TU information
                OutPut=$(find net/results_plus_TU/$file* -exec basename {} \; | sed -r 's/_extended.*/_extended_network_plus_TU_hubs_scores.txt/g' )
                printf "Making $OutPut\n"
                python scripts/getHubs.py net/results_plus_TU/$file* "2,3" | sort -n -k2,2 > cytoscape/$OutPut

        done

}

###################################################################################################################
############### THE FOLLOWING BLOCK IS DESIGN TO ASSIGN TRANSCRIPTION TYPE (convergent/diovergent)  ###############
###################################################################################################################

assignTypeTransc() {

        printf "${greenprintf}Divergent/convergent transcription block${NC}\n"
        # Five parameters must be supplied
        bash scripts/div_conver_transcript_pipeline.sh genomes/ genomes/gbk/ modifiedNetworks/ TUnits/bacteria/ others/transcription_type/

}


###################################### End functions ####################################



## Running functions/main-program by order

preprocessing
proteinortho
extendNet
orthosumm
extendNetPlusTU
cytoScape
coregulators
hubs
assignTypeTransc
