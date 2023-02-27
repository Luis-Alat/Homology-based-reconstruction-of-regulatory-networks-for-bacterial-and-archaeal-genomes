# Script to get topological and general metrics of networks using networkx (python)...
# for every predicted network 

printf "Getting metrics of bacterias...\n"

# Defining count to run for in parallel (every 40 iterations)
count_batches=1
batches=40

for file in bacteria_nets_predicted/net/Merge_plus_TU/*TU; do

	( baseName=$(basename $file)
	outputfilename=$(printf "bacteria_nets_predicted/metrics/${baseName}_metrics")
	python3 ../scripts/get_topological_metrics_Bact_Arche.py --inputFile $file --outputFile $outputfilename ) &
	((++count_batches)); [ "${count_batches}" -eq "${batches}" ] && count_batches=1 && wait

done

printf "Getting metrics of archaes...\n"

count_batches=1
for file in archae_nets_predicted/net/Merge_plus_TU/*TU; do

	( baseName=$(basename $file)
	outputfilename=$(printf "archae_nets_predicted/metrics/${baseName}_metrics")
    python3 ../scripts/get_topological_metrics_Bact_Arche.py --inputFile $file --outputFile $outputfilename ) &
	((++count_batches)); [ "${count_batches}" -eq "${batches}" ] && count_batches=1 && wait

done
