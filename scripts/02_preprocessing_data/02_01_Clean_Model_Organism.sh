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

