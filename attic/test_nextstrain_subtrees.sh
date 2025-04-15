i=2
number_of_clusters=$(wc -l bmtree_cluster_extraction.tsv | awk '{ print $1 }')
if [ "true" = "true" ]; then printf "while %s < %s\n" "$i" "$number_of_clusters"; fi
# We are going to copy groups.tsv within this loop, but that shouldn't cause issues.
# Older versions of this task mv'd groups.txt and were fine, so this should be extra-safe.
while [ $i -le $number_of_clusters ]
do
	head -$i "bmtree_cluster_extraction.tsv" | tail -n 1 > groups.tsv
	if [ "true" = "true" ]
	then
		printf "line %s grouped clusters now in groups.tsv as:\n" "$i"
		cat groups.tsv
		printf "\n"
	fi
	while IFS="	" read -r cluster samples
	do
		echo $samples > this_cluster_samples.txt
		sed -i 's/,/\n/g' this_cluster_samples.txt
		number_of_samples_in_cluster=$(wc -l this_cluster_samples.txt | awk '{ print $1 }')
		minimum_tree_size=$((number_of_samples_in_cluster+0))
		if [ "true" = "true" ]
		then
			printf "%s in cluster + 0 context: expecting %s samples in output" "$number_of_samples_in_cluster" "$minimum_tree_size"
			printf "passing this_cluster_samples.txt:\n"
			cat this_cluster_samples.txt
			printf "matUtils extract -i bmtree_raw.pb -j %s -s this_cluster_samples.txt -N %s -M bmtree_cluster_annotation.tsv\n" "$cluster" "$minimum_tree_size"
		fi
		matUtils extract -i "bmtree_raw.pb" -j "prefix$cluster" -s this_cluster_samples.txt -N $minimum_tree_size -M "bmtree_cluster_annotation.tsv"
		mv subtree-assignments.tsv "$cluster-subtree-assignments.tsv"
		cp groups.tsv "$cluster-groups.tsv"
		i=$((i+1))
	done < groups.tsv
	rm groups.tsv
done