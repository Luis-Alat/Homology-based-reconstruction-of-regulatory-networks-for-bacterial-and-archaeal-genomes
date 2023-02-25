#set -eE
#bash div_conver_transcript_pipeline.sh ../genomes/ ../genomes/gbk/ ../modifiedNetworks/ ../TUnits/bacteria/ tests/

# Defining names of files manually
fastaFiles=("GCF_000005845.2_E_coli_K12_genomic.faa" "GCF_000006765.1_ASM676v1_P_aeruginosa_PA01_genomic.faa" "GCF_000006945.2_ASM694v2_S_enterica_LT2_genomic.faa" "GCF_000009045.1_ASM904v1_B_subtilis_168_genomic.faa" "GCF_000009645.1_ASM964v1_S_aureus_N315_genomic.faa" "GCF_000195955.2_ASM19595v2_M_tuberculosis_H37Rv_genomic.faa")
gbkFiles=("NC_000913.3.gbk" "P_aeruginosa_GCF_000006765.1_ASM676v1_genomic.gbff" "salmonella_enterica_serovar_LT2_RN_GCF_000006945.2_ASM694v2_genomic.gbff" "B_subtilis_168_GCF_000009045.1_ASM904v1_genomic.gbff" "GCF_000009645.1_ASM964v1_genomic.gbff" "GCF_000195955.2_ASM19595v2_genomic.gbff")

# 5 paths must be passed
pathFasta=$1
pathGBK=$2
pathNetwork=$3
pathTus=$4
pathOutput=$5

# Iterating over all the organism model
for ((i=0;i<${#fastaFiles[@]};i++)); do

    # Defining input/output files
    inputGBK=$(printf "${pathGBK}${gbkFiles[$i]}")
    
    inputTUS=$(printf "${fastaFiles[$i]}" | cut -d"_" -f1,2 | sed -r 's/\.\w+//g')
    inputTUS=$(grep $inputTUS <(ls $pathTus))
    inputTUS=$(printf "${pathTus}${inputTUS}")

    inputFasta=$(printf "${pathFasta}${fastaFiles[$i]}")

    inputNetwork=$(printf "${fastaFiles[$i]}" | cut -d"_" -f1,2 | sed -r 's/\.\w+//g')
    inputNetwork=$(grep $inputNetwork <(ls $pathNetwork))
    inputNetwork=$(printf "${pathNetwork}${inputNetwork}")
    
    # This file will be overwriten later
    outputTMP=$(printf "${pathOutput}${fastaFiles[$i]}")

    # Getting convergent/divergent genes. Also get a table with the equivalences ids and type (is it a gene or TF?) based on the network
    python3 scripts/div_conver_transcript.py --inputPathGBK $inputGBK --inputPathTUS $inputTUS --inputPathFasta $inputFasta --inputPathNetwork $inputNetwork --outputPath $outputTMP

    # Defining a regex to know if the current iteration (below) is the divergent or convergent file
    regex="(convergent|divergent)"
    tmp_file=$(printf "${pathOutput}tmp")

    for file in $outputTMP*; do

        # Working over the convergent/divergent file. Skip (iteration) id table file generated above with the python script
        if [[ $file =~ $regex ]]; then

            table_ids=$(printf "${outputTMP}eq_ids")

            # Getting id ncbi and type (gene-TF)
            grep "TF" $table_ids | cut -f1,4 |
                # Creating multiple regex patterns to replace NCBI with type
                sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' |
                # Replacing over the convergent/divergent file (Only using important columns). Remove header #
                sed -f - <(cut -f1,3 $file) | grep -v "#" |
                # Replacing everything is not TF by gene
                sed -r 's/(\w+\.\w+)/gene/g' |
                # Merging the new columns with the original file without header
                paste <(grep -v "#" $file) - |
                # Adding the new header
                cat <(grep "#" $file | perl -nae 'chomp($_); print($_,"\ttype_element_1\ttype_element_2\n")') - |
                # Re-sorting columns. Save in a tmp file, later overwritten using an existent file
                awk -F"\t" 'BEGIN{OFS=FS}{print $1,$6,$2,$3,$7,$4,$5}' > $tmp_file

            mv $tmp_file $file
        
        fi

    done

done