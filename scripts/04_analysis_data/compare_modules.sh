set -e

###### Main line of this script #######################
# Removing headers

#grep -v "module" ../GCF_000005845.2_E_coli_K12_genomic_extended_network_Modules.txt |

# Sorting by column number

#sort -k2,2 -n |

# Getting first module (and then in a for loop for every module)

#grep -E "\s+0" |

# Cutting the IDs of the module

#cut -f1 |

# Searching and getting orthologs of that module; then search them in the modules of others organism (Bacsu in this case)

#grep -Fwif - ../../net/Bacsu/ortoFiles/GCF_000005845.2_E_coli_K12_genomic_ortho.txt | cut -f1 |

# Now, we are searching for ortologs in bacsu and ordering by module, in this way we are able to get ortologs of ecoli of the first module in all the others modules in Bacsu

#grep -Fwif - ../GCF_000009045.1_ASM904v1_B_subtilis_168_genomic_extended_network_Modules.txt | sort -k2,2 -n

########################################################

GenomeNames=("GCF_000009045.1_ASM904v1_B_subtilis_168" "GCF_000005845.2_E_coli_K12" "GCF_000195955.2_ASM19595v2_M_tuberculosis_H37R" "GCF_000006765.1_ASM676v1_P_aeruginosa_PA01" "GCF_000006945.2_ASM694v2_S_enterica_LT2" "GCF_000009645.1_ASM964v1_S_aureus_N315")
FolderNames=("Bacsu"  "Ecoli"  "Myctu"  "Pseudo"  "Salmon"  "Staph")

for prefix in "extended_network_Modules.txt" "extended_network_plus_TU_Modules.txt" "net_Modules.txt"; do

	for ((i=0;i<${#GenomeNames[@]};i++)); do

		# Counting how many times will be necesary iterate (number of modules)
		CountMod=$(grep -v "module" "../${GenomeNames[$i]}"*"${prefix}" | cut -f2 | sort -n | uniq | wc -l)

		# Finding the size of each module in the organism ${GenomeNames[$i]} and later using it as a pattern to find each module and add that value in the final file
		grep -v "module" "../${GenomeNames[$i]}"*"${prefix}" | cut -f2 | sort -n | uniq -c | sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' | awk -F$"\t" '{OFS=FS}{print $2,$2,$1}' | sed -r 's/^/s\/^/g' | sed -r 's/\t/$\//1' | sed -r 's/$/\/g;/g' > patterns_source_organism.txt

		# Creating an empty final file
		printf "" > "${FolderNames[$i]}_${prefix}_summary_compared_modules.txt"

		# Creating and second empty final file
		printf "" > "${FolderNames[$i]}_${prefix}_ids_compared_modules.txt"

		for ((j=0;j<${#FolderNames[@]};j++)); do

			# Do not iterate over itself
			if [[ $i -eq $j  ]]; then
				continue
			fi

			# Generating an empty temporaly file. It useful to manage every output
			printf "" > tmp_2

			for ((m=0;m<$CountMod;m++)); do

				# Information about the iteration status
				printf "MODULE $m of ${FolderNames[$i]} against ${FolderNames[$j]} modules\n"

				# Here, the equivalent (proteinortho processed ouput is used) of one module of the organism $i is going to be search in all modules in $j organism
				grep -v "module" "../${GenomeNames[$i]}"*"${prefix}" | sort -k2,2 -n | grep -E "\s+${m}$" | cut -f1 | grep -Fwif - "../../net/${FolderNames[$j]}/ortoFiles/${GenomeNames[$i]}"* | cut -f1 | grep -Fwif - "../${GenomeNames[$j]}"*"${prefix}" | sort -k2,2 -n > tmp_1

				# A new file is created based on the previous results adding information about how many matches were found in every module as well as organism information
				cut -f2 tmp_1 | uniq -c | sed -r 's/^\s+//g' | awk '{print $2,$1}' | sed -r 's/\s+/\t/g' | awk -F$"\t" -v sourceMod=${FolderNames[$i]} -v targetMod=${FolderNames[$j]} -v module=$m '{OFS=FS}{print sourceMod,targetMod,module,$0}' >> tmp_2

				# A new file is created based on the tmp_1 file. In this case this new file will describe the id where a match was found between both organismin every module

				cut -f1 tmp_1 | sort | uniq | grep -Fwif - ../../net/${FolderNames[$j]}/ortoFiles/${GenomeNames[$i]}* | awk -F$"\t" '{OFS=FS}{print $1,$1,$2}' | sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//1' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' | sed -f - tmp_1 | perl -F"\t" -nae 'if(scalar(@F) > 3){$string=""; for($i=1;$i<(scalar(@F) - 1);$i++){ chomp($F[$i]); $string = $string . $F[$i] . ","  } print($F[0],"\t",$string,"\t",$F[$#F]) } else{print $_}' | sed -r 's/,\t/\t/g' | awk -F$"\t" -v mod=$m -v source=${FolderNames[$i]} -v target=${FolderNames[$j]} '{OFS=FS}{print source,mod,$2,target,$3,$1}' > tmp_3

				# Removing wrong retrieved id of the modules

				grep -v "module" "../${GenomeNames[$i]}"*"${prefix}" | grep -E "\s+${m}$" | cut -f1 | sed -r 's/\//\\\//g' | perl -nae 'chomp($_); print($_,"\t",$_,"\n")' | sed -r 's/$/_MATCHES/g' | sed -r 's/^/s\/\\</g' | sed -r 's/\t/\\>\//g' | sed -r 's/$/\/gI;/g' | sed -r 's/\./\\./g' | sed -f - tmp_3 | perl -F"\t" -nae '@Items=split(",",$F[2]); $string=""; for($i=0;$i<scalar(@Items);$i++){ if($Items[$i] =~ /_MATCHES/){ $string=$string . $Items[$i] . "," }  } print($F[0],"\t",$F[1],"\t",$string,"\t",$F[3],"\t",$F[4],"\t",$F[5])' | sed -r 's/_MATCHES//g' | sed -r 's/,\t/\t/g' | sed -r 's/,,/,/g' >> "${FolderNames[$i]}_${prefix}_ids_compared_modules.txt"

			done

			# Modules are searched and sizes of each module are added using the previous file. Output is going to be saved in a new file of one column
			sed -f patterns_source_organism.txt <(cut -f3 tmp_2) | cut -f2 > size_mod_source.txt
			grep -v "module" ../${GenomeNames[$j]}*"${prefix}" | cut -f2 | sort -n | uniq -c | sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' | awk -F$"\t" '{OFS=FS}{print $2,$2,$1}' | sed -r 's/^/s\/^/g' | sed -r 's/\t/$\//1' | sed -r 's/$/\/g;/g' | sed -f - <(cut -f4 tmp_2) | cut -f2 > size_mod_target.txt

			# Merging into a file organismo "source", organism "target", module name of "source", module name of "target", number of matches, size of module "source", size of module "target"
			paste tmp_2 size_mod_source.txt size_mod_target.txt >> "${FolderNames[$i]}_${prefix}_summary_compared_modules.txt"

		done
	done
done

rm size_mod_source.txt size_mod_target.txt tmp_1 tmp_2 tmp_3 patterns_source_organism.txt

for file in *summary_compared*; do

	printf "Getting co-regulation value of ${file}\n"
	perl -F"\t" -nae 'chomp($_); chomp($F[6]); $metric = ($F[4]) / ($F[5] + $F[6] - $F[4]);  print($_,"\t",$metric,"\n")' $file > tmp_1 && mv tmp_1 $file

done
