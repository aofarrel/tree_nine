version 1.0
# Tasks included:
# * annotate: matUtils annotate
# * convert_to_nextstrain_single: matUtils extract -j
# * convert_to_nextstrain_subtrees: matUtils extract -j -s
# * convert_to_newick: matUtils extract -t
# * convert_to_taxonium: usher_to_taxonium
# * extract: matUtils extract -s
# * mask: matUtils mask
# * matrix: python3 distancematrix_nwk.py
# * reroot: matUtils extract -y 
# * summary: matUtils summary
# * usher_sampled_diff: usher-sampled --diff

task extract {
	input {
		File input_mat
		File samples
		Int? nearest_k
		
		Int addldisk = 10
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
	}
	Int disk_size = ceil(size(input_mat, "GB")) + addldisk
	String output_mat = basename(input_mat, ".pb") + ".subtree" + ".pb"
	
	command <<<
	if [[ "~{nearest_k}" = "" ]]
	then
		matUtils extract -i ~{input_mat} -s ~{samples} -o ~{output_mat}
	else
		matUtils extract -i ~{input_mat} -K ~{samples}:~{nearest_k} -o ~{output_mat}
	fi
	>>>
	
	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: preempt
	}
	
	output {
		File subtree = output_mat
	}
}

task mask {
	input {
		File input_mat
		File mask_tsv
		
		Int addldisk = 10
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
	}
	Int disk_size = ceil(size(input_mat, "GB")) + addldisk
	String output_mat = basename(input_mat, ".pb") + ".masked" + ".pb"
	
	command <<<
	matUtils mask -i ~{input_mat} --mask-mutations ~{mask_tsv} -o ~{output_mat}
	>>>
	
	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: preempt
	}
	
	output {
		File masked_tree = output_mat
	}
}

task reroot {
	input {
		File input_mat
		String reroot_to_this_node

		Int addldisk = 10
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
	}
	Int disk_size = ceil(size(input_mat, "GB")) + addldisk
	String output_mat = basename(input_mat, ".pb") + ".reroot_to_~{reroot_to_this_node}" + ".pb"

	command <<<
	if [[ "~{reroot_to_this_node}" = "" ]]
	then
		echo "You need to specify the node to reroot upon"
		exit 1
	fi
	matUtils extract -i "~{input_mat}" -y "~{reroot_to_this_node}" -o "~{output_mat}"
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: preempt
	}
	
	output {
		File rerooted_tree = output_mat
	}
}

task summarize {
	# Generates most of the possible outputs of matUtils summarize:
	#
	# --samples (-s): Write a two-column tsv listing all samples in the tree and their parsimony score (terminal branch length). Auspice-compatible.
	# --clades (-c): Write a tsv listing all clades and the count of associated samples in the tree.
	# --mutations (-m): Write a tsv listing all mutations in the tree and their occurrence count.
	# --aberrant (-a): Write a tsv listing potentially problematic nodes, including duplicates and internal nodes with no mutations and/or branch length 0.
	# --haplotype (-H): Write a tsv listing haplotypes represented by comma-delimited lists of mutations and their count across the tree.
	# --sample-clades (-C): Write a tsv listing all samples and their clades.
	# --calculate-roho (-R): Write a tsv listing, for each mutation occurrence that is valid, the number of offspring and other numbers for RoHo calculation.
	#
	# Two outputs are not generated:
	# * expanded_roho: this slows things down too much
	# * translate: this would require taking in a gtf and ref genome

	input {
		File? input_mat
		String? prefix_outs

		Int addldisk = 10
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
	}
	Int disk_size = if defined(input_mat) then ceil(size(input_mat, "GB")) + addldisk else addldisk
	String prefix = select_first([prefix_outs, ""])

	command <<< 
	if [[ "~{input_mat}" = "" ]]
	then
		i="/HOME/usher/example_tree/for_debugging_only__tb_7K_noQC_diffs_mask2ref.L.fixed.pb"
	else
		i="~{input_mat}"
	fi
	
	matUtils summary -i "$i" > "~{prefix}.summary.txt"
	matUtils summary -i "$i" -A # samples, clades, mutations, aberrant
	matUtils summary -i "$i" -H haplotypes.tsv
	matUtils summary -i "$i" -C sample_clades.tsv
	matUtils summary -i "$i" -R roho.tsv
	for file in *.tsv
	do
		mv -- "$file" "~{prefix}.${file}"
	done
	
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File summary = prefix + ".summary.txt"
		File samples = prefix + ".samples.tsv"
		File clades = prefix + ".clades.tsv"
		File mutations = prefix + ".mutations.tsv"
		File aberrant = prefix + ".aberrant.tsv"
		File haplotype = prefix + ".haplotypes.tsv"
		File sample_clades = prefix + ".sample_clades.tsv"
		File calculate_roho = prefix + ".roho.tsv"
	}
}

task annotate {
	input {
		File? input_mat
		File metadata_tsv # only can annotate one column at a time
		String outfile_annotated

		Int addldisk = 10
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
	}
	Int disk_size = ceil(size(input_mat, "GB")) + ceil(size(metadata_tsv, "GB")) + addldisk

	command <<< 
	matUtils annotate -i "~{input_mat}" -P "~{metadata_tsv}" -o "~{outfile_annotated}"
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File annotated_tree = outfile_annotated
	}
}

task convert_to_newick_subtrees_by_cluster {
	input {
		File input_mat
		File metadata_tsv
		File grouped_clusters
		Int context_samples
		Int memory = 32
		Boolean debug = true
		String prefix = ""
	}

	command <<<
		i=2
		number_of_clusters=$(wc -l ~{grouped_clusters} | awk '{ print $1 }')
		if [ ~{debug} = "true" ]; then printf "while %s < %s\n" "$i" "$number_of_clusters"; fi
		# shellcheck disable=SC2086
		# We are going to copy groups.tsv within this loop, but that shouldn't cause issues.
		# Older versions of this task mv'd groups.txt and were fine, so this should be extra-safe.
		while [ $i -le $number_of_clusters ]
		do
			head -$i "~{grouped_clusters}" | tail -n 1 > groups.tsv
			if [ ~{debug} = "true" ]
			then
				printf "line %s grouped clusters now in groups.tsv as:\n" "$i"
				cat groups.tsv
				printf "\n"
			fi
			# shellcheck disable=SC2094
			while IFS="	" read -r cluster samples
			do
				# shellcheck disable=SC2086
				echo $samples > this_cluster_samples.txt
				sed -i 's/,/\n/g' this_cluster_samples.txt
				number_of_samples_in_cluster=$(wc -l this_cluster_samples.txt | awk '{ print $1 }')
				minimum_tree_size=$((number_of_samples_in_cluster+~{context_samples}))
				if [ ~{debug} = "true" ]
				then
					printf "%s in cluster + ~{context_samples} context: expecting %s samples in output" "$number_of_samples_in_cluster" "$minimum_tree_size"
					printf "passing this_cluster_samples.txt:\n"
					cat this_cluster_samples.txt
					printf "matutils extract -i ~{input_mat} -t %s -s this_cluster_samples.txt -N %s \n" "$cluster" "$minimum_tree_size"
				fi
				matUtils extract -i "~{input_mat}" -t "~{prefix}$cluster" -s this_cluster_samples.txt -N $minimum_tree_size -M "~{metadata_tsv}"
				if [ ~{debug} = "true" ]
				# for some reason, subtrees seem to end up with .nw as their extension
				for tree in *.nw; do
					mv -- "$tree" "${tree%.nw}.nwk"
				done
				mv subtree-assignments.tsv "$cluster-subtree-assignments.tsv"
				cp groups.tsv "$cluster-groups.tsv"
				then
					printf "Finished %s and moving files.\n" "$cluster"
					ls -lha
				fi
				i=$((i+1))
			done < groups.tsv
			rm groups.tsv
		done
	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: 1
	}

	output {
		Array[File] newick_subtrees = glob("*.nwk")
		Array[File] subtree_assignments = glob("*-subtree-assignments.tsv")
		Array[File] groups = glob("*groups.tsv")
		Array[File] metadata_tsvs = glob("*.tsv")  # for auspice.us, which supports nwk
	}
}

task convert_to_nextstrain_subtrees_by_cluster {
	input {
		File input_mat
		File metadata_tsv
		File grouped_clusters
		Int context_samples
		Int memory = 32
		Boolean debug = true
		String prefix = ""
	}

	command <<<
		i=2
		number_of_clusters=$(wc -l ~{grouped_clusters} | awk '{ print $1 }')
		if [ ~{debug} = "true" ]; then printf "while %s < %s\n" "$i" "$number_of_clusters"; fi
		# shellcheck disable=SC2086
		# We are going to copy groups.tsv within this loop, but that shouldn't cause issues.
		# Older versions of this task mv'd groups.txt and were fine, so this should be extra-safe.
		while [ $i -le $number_of_clusters ]
		do
			head -$i "~{grouped_clusters}" | tail -n 1 > groups.tsv
			if [ ~{debug} = "true" ]
			then
				printf "line %s grouped clusters now in groups.tsv as:\n" "$i"
				cat groups.tsv
				printf "\n"
			fi
			# shellcheck disable=SC2094
			while IFS="	" read -r cluster samples
			do
				# shellcheck disable=SC2086
				echo $samples > this_cluster_samples.txt
				sed -i 's/,/\n/g' this_cluster_samples.txt
				number_of_samples_in_cluster=$(wc -l this_cluster_samples.txt | awk '{ print $1 }')
				minimum_tree_size=$((number_of_samples_in_cluster+~{context_samples}))
				if [ ~{debug} = "true" ]
				then
					printf "%s in cluster + ~{context_samples} context: expecting %s samples in output" "$number_of_samples_in_cluster" "$minimum_tree_size"
					printf "passing this_cluster_samples.txt:\n"
					cat this_cluster_samples.txt
					printf "matutils extract -i ~{input_mat} -j %s -s this_cluster_samples.txt -N %s -M ~{metadata_tsv}\n" "$cluster" "$minimum_tree_size"
				fi
				matUtils extract -i "~{input_mat}" -j "~{prefix}$cluster" -s this_cluster_samples.txt -N $minimum_tree_size -M "~{metadata_tsv}"
				mv subtree-assignments.tsv "$cluster-subtree-assignments.tsv"
				cp groups.tsv "$cluster-groups.tsv"
				i=$((i+1))
			done < groups.tsv
			rm groups.tsv
		done
	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: 1
	}

	output {
		Array[File] nextstrain_subtrees = glob("*.json")
		Array[File] subtree_assignments = glob("*-subtree-assignments.tsv")
		Array[File] groups = glob("*groups.tsv")
	}
}


task convert_to_nextstrain_subtrees {
	# based loosely on Marc Perry's version
	input {
		File         input_mat # aka tree_pb

		Int          memory = 32
		Array[File?] metadata_files
		Int?         nearest_k
		String       outfile_nextstrain = "nextstrain"
		File?        selected_samples
		Int          treesize = 0
	}

	String metadata = if length(metadata_files) > 0 then "-M" else ""

	command <<<
		METAFILES_OR_EMPTY="~{sep=',' metadata_files}"
		if [[ "~{selected_samples}" == "" ]]
		then
			matUtils extract -i	~{input_mat} -S sample_paths.txt
			cut -f1 sample_paths.txt | tail -n +2 > sample.ids
			if [[ "~{nearest_k}" == "" ]]
			then
				matUtils extract -i ~{input_mat} -j ~{outfile_nextstrain} -s sample.ids -N ~{treesize} ~{metadata} $METAFILES_OR_EMPTY
			else
				matUtils extract -i ~{input_mat} -j ~{outfile_nextstrain} -K sample.ids:~{nearest_k} -N ~{treesize} ~{metadata} $METAFILES_OR_EMPTY
			fi
		else
			if [[ "~{nearest_k}" == "" ]]
			then
				matUtils extract -i ~{input_mat} -j ~{outfile_nextstrain} -s ~{selected_samples} -N ~{treesize} ~{metadata} $METAFILES_OR_EMPTY
			else
				matUtils extract -i ~{input_mat} -j ~{outfile_nextstrain} -K ~{selected_samples}:~{nearest_k} -N ~{treesize} ~{metadata} $METAFILES_OR_EMPTY
			fi
		fi
		ls -lha
	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: 1
	}

	output {
		Array[File] nextstrain_subtrees = glob("*.json")
	}
}

task convert_to_nextstrain_single {
	input {
		File input_mat # aka tree_pb
		Int memory = 32
		String outfile_nextstrain
		Array[File?] metadata_files
	}
	
	String metadata = if length(metadata_files) > 0 then "-M" else ""

	command <<<
		METAFILES_OR_EMPTY="~{sep=',' metadata_files}"
		matUtils extract -i ~{input_mat} ~{metadata} $METAFILES_OR_EMPTY -j ~{outfile_nextstrain}
	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: 1
	}

	output {
		File nextstrain_singular_tree = outfile_nextstrain
	}
}

task convert_to_nextstrain_single_terra_compatiable {
	input {
		File input_mat # aka tree_pb
		Int memory = 32
		String outfile_nextstrain
		File? one_metadata_file
	}

	command <<<
		if [ "~{one_metadata_file}" != "" ]
		then
			matUtils extract -i ~{input_mat} -M ~{one_metadata_file} -j ~{outfile_nextstrain}
		else
			matUtils extract -i ~{input_mat} -j ~{outfile_nextstrain}
		fi
	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: 1
	}

	output {
		File nextstrain_singular_tree = outfile_nextstrain
	}
}

task convert_to_newick {
	input {
		File input_mat
		String outfile_nwk
	}

	command <<<
		matUtils extract -i ~{input_mat} -t ~{outfile_nwk}
	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 8
		disks: "local-disk " + 100 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: 8 + " GB"
		preemptible: 1
	}

	output {
		File newick_tree = outfile_nwk
	}
}

task usher_sampled_diff {
	input {
		Int batch_size_per_process = 5
		Boolean detailed_clades
		File diff
		File? input_mat
		Int optimization_radius = 0
		Int max_parsimony_per_sample = 1000000
		Int max_uncertainty_per_sample = 1000000
		String output_mat
		File? ref_genome

		# WDL specific -- note that cpu does not directly set usher's
		# threads argument, but it does affect the number of cores
		# available for use (by default usher uses all available)
		Int addldisk = 10
		Int cpu = 40      # needed for CPDH but overkill for small numbers of samples -- 8 (yes, eight!) would do fine
		Int memory = 32   # needed for CDPH but overkill for small numbers of samples -- 16 would do fine
		Int preempt = 1
	}

	Int disk_size = ceil(size(diff, "GB")) + ceil(size(ref_genome, "GB")) +  ceil(size(input_mat, "GB")) + addldisk
	String D = if !(detailed_clades) then "" else "-D "
	#String ref = select_first([ref_genome, "/HOME/usher/ref/Ref.H37Rv/ref.fa"])
	#String i = select_first([input_mat, "/HOME/usher/example_tree/for_debugging_only__tb_7K_noQC_diffs_mask2ref.L.fixed.pb"])

	command <<<
		if [[ "~{input_mat}" = "" ]]
		then
			i="/HOME/usher/example_tree/for_debugging_only__tb_7K_noQC_diffs_mask2ref.L.fixed.pb"
		else
			i="~{input_mat}"
		fi

		if [[ "~{ref_genome}" = "" ]]
		then
			ref="/HOME/usher/ref/Ref.H37Rv/ref.fa"
		else
			ref="~{ref_genome}"
		fi
		
		echo "~{input_mat}"
		echo $i
		echo "~{ref_genome}"
		echo $ref
		echo "------------------"
		ls -lha
		echo "------------------"

		usher-sampled ~{D} --optimization_radius=~{optimization_radius} \
			-e ~{max_uncertainty_per_sample} \
			-E ~{max_parsimony_per_sample} \
			--batch_size_per_process ~{batch_size_per_process} \
			--diff "~{diff}" \
			-i "$i" \
			--ref "$ref" \
			-o "~{output_mat}"
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4_4"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File usher_tree = output_mat
		File? clades = "clades.txt" # only if detailed_clades = true
	}

}

task matOptimize {
	input {
		File input_mat
		Int max_hours = 1
		Float min_improvement = 0.00000001

		Int addldisk = 100
		Int cpu = 24
		Int memory = 32
		Int preempt = 1
	}

	Int disk_size = ceil(size(input_mat, "GB")) + addldisk
	String outfile = basename(input_mat) + "_optimized.pb"

	command <<<
	matOptimize -i "~{input_mat}" --max-hours ~{max_hours} --min-improvement ~{min_improvement} -o "~{outfile}"
	>>> 

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4ash_1"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File optimized_tree = outfile
	}
}

task convert_to_taxonium {
	input {
		File input_mat
		String outfile_taxonium

		Int addldisk = 100
		Int cpu = 12
		Int memory = 16
		Int preempt = 1
	}

	Int disk_size = ceil(size(input_mat, "GB")) + addldisk

	command <<<
		echo "booted into Docker successfully"
		echo "input file: ~{input_mat}"
		ls -lha ~{input_mat}
		echo "running usher_to_taxonium..."
		usher_to_taxonium -i "~{input_mat}" -o "~{outfile_taxonium}"
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/sranwrp:1.1.6"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File taxonium_tree = outfile_taxonium
	}
}

task matrix {
	input {
		File input_nwk
		Boolean only_matrix_special_samples
		File? special_samples
		String? outfile_matrix
	}
	
	command <<<
	wget https://raw.githubusercontent.com/aofarrel/parsevcf/1.3.1/distancematrix_nwk.py
	if [[ "~{only_matrix_special_samples}" = "true" ]]
	then
		samples=$(< "~{special_samples}" tr -s '\n' ',' | head -c -1)
		echo "Samples that will be in the distance matrix: $samples"
		python3 distancematrix_nwk.py "~{input_nwk}" --samples "$samples" -v
	else
		python3 distancematrix_nwk.py "~{input_nwk}" -v
	fi
	if [[ "~{outfile_matrix}" != "" ]]
	then
		icky_filename=$(find . -maxdepth 1 -name "*.tsv")
		mv "$icky_filename" "~{outfile_matrix}"
	fi
	>>>
	
	runtime {
		cpu: 8
		disks: "local-disk " + 100 + " SSD"
		docker: "ashedpotatoes/sranwrp:1.1.15"
		memory: 8 + " GB"
		preemptible: 1
	}

	output {
		File out_matrix = glob("*.tsv")[0]
	}
	
}

task matrix_and_find_clusters {
	input {
		File input_nwk
		File? persistent_cluster_tsv
		Boolean only_matrix_special_samples
		File? special_samples
		Int distance
		Int cpu = 8
		Int memory = 8
	}
	
	command <<<
	# temporarily overwriting DM script in the docker image until it's more stable
	rm /scripts/distancematrix_nwk.py
	wget https://raw.githubusercontent.com/aofarrel/parsevcf/refs/heads/slight-refactor/distancematrix_nwk.py
	mv distancematrix_nwk.py /scripts/distancematrix_nwk.py

	# TODO: persistent cluster IDs will kind of break on the "lonely" cluster. modify marc's script so
	# that if a sample in previous run was "lonely" and in most recent run is not "lonely", give it a
	# new name (that hasn't been used before) and mark it as a previously-run-newly-clustered sample.
	if [[ "~{only_matrix_special_samples}" = "true" ]]
	then
		samples=$(< "~{special_samples}" tr -s '\n' ',' | head -c -1)
		echo "Samples that will be in the distance matrix: $samples"
		python3 /scripts/distancematrix_nwk.py "~{input_nwk}" --samples "$samples" -d ~{distance}
	else
		python3 /scripts/distancematrix_nwk.py "~{input_nwk}" -d ~{distance}
	fi
	if [[ "~{persistent_cluster_tsv}" != "" ]]
	then
	perl /scripts/marcs_incredible_script.pl ""
		awk 'NR==FNR {keys[$1]; next} $1 in keys' "$(find . -name '*_cluster_annotation.tsv' | head -n 1)" ~{persistent_cluster_tsv} > filtered_latest_clusters.tsv
		awk 'NR==FNR {keys[$1]; next} $1 in keys' ~{persistent_cluster_tsv} "$(find . -name '*_cluster_annotation.tsv' | head -n 1)" > filtered_persistent_clusters.tsv
		perl /scripts/marcs_incredible_script.pl filtered_persistent_clusters.tsv filtered_latest_clustes.tsv
	fi

	>>>
	
	runtime {
		cpu: cpu
		disks: "local-disk " + 100 + " SSD"
		docker: "ashedpotatoes/sranwrp:1.1.26"
		memory: memory + " GB"
		preemptible: 1
	}

	output {
		Array[File] out_matrices = glob("*_dmtrx.tsv")
		File samp_cluster = glob("*_cluster_annotation.tsv")[0]
		File cluster_samps = glob("*_cluster_extraction.tsv")[0]
		File samp_UUID = glob("*_cluster_UUIDs.tsv")[0]
		File? persistent_cluster_translator = "mapped_persistent_cluster_ids_to_new_cluster_ids.tsv"
		Int n_clusters = read_int("n_clusters")
		Int n_samples_in_clusters = read_int("n_samples_in_clusters")
		Int total_samples_processed = read_int("total_samples_processed")
	}
	
}

task cluster_CDPH_method {
	# find_clusters.py: Generates 20-10-5 clusters and distance matrices (normal and backmasked)
	# process_clusters.py: Persistent cluster IDs, subtrees, and MR upload
	# Any clusters that have at least one sample without a diff file will NOT be backmasked
	input {
		File input_mat_with_new_samples
		String today
		Boolean upload_clusters_to_microreact = true
		File? persistent_denylist

		# Not actually optional, just marked as such due to WDL limitations when calling via Tree Nine
		File? persistent_ids
		File? persistent_cluster_meta
		File combined_diff_file           # used for local masking
		File? previous_run_cluster_json   # for comparisons

		# keep these files in the workspace bucket for now
		File? microreact_update_template_json # must be called REALER_template.json for now
		File? microreact_blank_template_json # must be called BLANK_template.json
		File? microreact_key

		Boolean only_matrix_special_samples
		Boolean inteight = false
		File? special_samples
		
		String? shareemail
		Int preempt = 0 # only set if you're doing a small test run
		Int memory = 50
		Boolean debug = true
		
		# temporary overrides
		File? find_clusters_script_override
		File? process_clusters_script_override
		File? summarize_changes_script_override

		#Int batch_size_per_process = 5
		#Boolean detailed_clades
		#Int optimization_radius = 0
		#Int max_parsimony_per_sample = 1000000
		#Int max_uncertainty_per_sample = 1000000
		#String output_mat
		
	}
	Array[Int] cluster_distances = [20, 10, 5] # CHANGING THIS WILL BREAK SECOND SCRIPT!
	String arg_denylist = if defined(persistent_denylist) then "--dl ~{persistent_denylist}" else ""
	String arg_shareemail = if defined(shareemail) then "-s ~{shareemail}" else ""
	String arg_microreact = if upload_clusters_to_microreact then "--yes_microreact" else ""
	String arg_token = if upload_clusters_to_microreact then "--token" else "" # cannot include microreact_key or else it will be gs://
	String arg_ieight = if inteight then "--int8" else ""

	command <<<
		matUtils extract -i ~{input_mat_with_new_samples} -t A_big.nwk
		cp ~{input_mat_with_new_samples} .

		if [[ "~{find_clusters_script_override}" == '' ]]
		then
			wget https://raw.githubusercontent.com/aofarrel/tree_nine/refs/heads/main/find_clusters.py
			mv find_clusters.py /scripts/find_clusters.py
		else
			mv "~{find_clusters_script_override}" /scripts/find_clusters.py
		fi

		if [[ "~{process_clusters_script_override}" == '' ]]
		then
			wget https://raw.githubusercontent.com/aofarrel/tree_nine/refs/heads/main/process_clusters.py
			mv process_clusters.py /scripts/process_clusters.py
		else
			mv "~{process_clusters_script_override}" /scripts/process_clusters.py
		fi

		if [[ "~{summarize_changes_script_override}" == '' ]]
		then
			wget https://raw.githubusercontent.com/aofarrel/tree_nine/refs/heads/main/summarize_changes.py
			mv summarize_changes.py /scripts/summarize_changes.py
		else
			mv "~{summarize_changes_script_override}" /scripts/summarize_changes.py
		fi

		wget https://gist.githubusercontent.com/aofarrel/6a458634abbca4eb16d120cc6694d5aa/raw/d6f5466e04394ca38f1a92b1580a9a5bd436bbc8/marcs_incredible_script_update.pl
		mv marcs_incredible_script_update.pl /scripts/marcs_incredible_script_update.pl

		wget https://raw.githubusercontent.com/aofarrel/tsvutils/refs/heads/main/extract_long_rows_and_truncate.sh
		mv extract_long_rows_and_truncate.sh /scripts/strip_tsv.sh
		wget https://raw.githubusercontent.com/aofarrel/tsvutils/refs/heads/main/equalize_tabs.sh
		mv equalize_tabs.sh /scripts/equalize_tabs.sh

		# never ever ever put this in the docker image (okay not really but like. for now.)
		mv ~{microreact_update_template_json} .
		mv ~{microreact_blank_template_json} .

		echo "Finished moving stuff, here is workdir"
		tree

		CLUSTER_DISTANCES="~{sep=',' cluster_distances}"
		FIRST_DISTANCE="${CLUSTER_DISTANCES%%,*}"
		OTHER_DISTANCES="${CLUSTER_DISTANCES#*,}"
		echo "cluster distances $CLUSTER_DISTANCES"
		echo "First distance $FIRST_DISTANCE"
		echo "Other distances $OTHER_DISTANCES"

		# TODO: on very large runs, the size of $/samples may eventually cause issues with ARG_MAX
		# should be fine for our purposes though

		# shellcheck disable=SC2086
		if [[ "~{only_matrix_special_samples}" = "true" ]]
		then
			samples=$(< "~{special_samples}" tr -s '\n' ',' | head -c -1)
			echo "Samples that will be in the distance matrix: $samples"

			python3 /scripts/find_clusters.py \
				~{input_mat_with_new_samples} \
				--samples $samples \
				--collection-name big \
				-t NB \
				-d "$FIRST_DISTANCE" \
				-rd "$OTHER_DISTANCES" \
				-v ~{arg_ieight}
		else
			python3 /scripts/find_clusters.py \
				~{input_mat_with_new_samples} \
				--collection-name big \
				-t NB \
				-d "$FIRST_DISTANCE" \
				-rd "$OTHER_DISTANCES" \
				-v ~{arg_ieight}
		fi

		set -eux pipefail  # now we can set pipefail since we are no longer returning non-0s
		echo "Current sample information:"
		cat latest_samples.tsv
		echo "Contents of workdir:"
		tree
		# A_big.nwk									big tree, nwk format
		# LONELY-subtree-n.nwk (n as variable)		subtrees (usually multiple) of unclustered samples
		# lonely-subtree-assignments.tsv			which subtree each unclustered sample ended up in
		# cluster_annotation_workdirIDs.tsv			can be used to annotate by nonpersistent cluster (but isn't, at least not yet)
		# latest_samples.tsv						used by persistent ID script
		# n_big_clusters (n as constant)			# of 20SNP clusters
		# n_samples_in_clusters (n as constant)		# of samples that clustered
		# n_samples_processed (n as constant)		# of samples processed by find_clusters.py
		# n_unclustered (n as constant)				# of samples that failed to cluster
		# ...and one distance matrix per cluster, and also one(?) subtree per cluster. Later, there will be two of each per cluster, once backmasking works!

		if [ "~{persistent_ids}" != "" ]
		then
			mkdir logs
			echo "Running second script"

			python3 /scripts/process_clusters.py --latestsamples latest_samples.tsv --persistentids ~{persistent_ids} -pcm ~{persistent_cluster_meta} ~{arg_token} ~{microreact_key} -mat ~{input_mat_with_new_samples} -cd ~{combined_diff_file} ~{arg_denylist} ~{arg_shareemail} ~{arg_microreact} --today ~{today} --allsamples $samples

			zip -r logs.zip ./logs
		fi

		if [ -f "rosetta_stone_20_merges.tsv" ]
		then
			cat rosetta_stone_20_merges.tsv
		fi
		if [ -f "rosetta_stone_10_merges.tsv" ]
		then
			cat rosetta_stone_10_merges.tsv
		fi
		if [ -f "rosetta_stone_5_merges.tsv" ]
		then
			cat rosetta_stone_5_merges.tsv
		fi

		#if [ "~{previous_run_cluster_json}" != "" ]
		#then
		#		echo "Running third script"
		#		python3 /scripts/summarize_changes.py ~{previous_run_cluster_json} all_cluster_information.json
		#fi
		if [ ~{debug} = "true" ]; then ls -lha; fi
		
		rm REALER_template.json # avoid globbing with the subtrees

	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.4ash_1"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File? logs = "logs.zip"
		File? change_report_full = "change_report_full"+today+".txt"
		File? change_report_cdph = "change_report_cdph"+today+".txt"

		# IMPORTANT FILES THAT SHOULD ALWAYS GO INTO SUBSEQUENT RUNS IF THEY EXIST
		# Try to avoid globbing where possible to make finding outs in Terra bucket easier
		# since globs create a folder with a randomized name, which is annoying
		File? clusterid_denylist = "clusterid_denylist.txt"
		File? new_persistent_ids = "persistentIDS" + today + ".tsv"
		File? new_persistent_meta = "persistentMETA" + today + ".tsv"
		File? final_cluster_information_json = "all_cluster_information" + today + ".json"
		File? change_report_json = "change_report" + today + ".json"

		# brand new samples list, not fully finished processing but here ya go
		File? new_samples = "new_samples" + today + ".tsv"

		# trees, all in nwk format for now
		# we have to comment some of these out because we're hitting the limit of what Terra can write to a data table
		# A = not internally masked
		# B = internally masked
		File? abig_tree = "A_big.nwk"
		File? bbig_tree = "b000000.nwk"
		#Array[File]? unclustered_subtrees = glob("LONELY*.nwk") # !UnnecessaryQuantifier
		#Array[File]? acluster_trees = glob("a*.nwk")            # !UnnecessaryQuantifier
		Array[File]? bcluster_trees = glob("b*.nwk")            # !UnnecessaryQuantifier
		Array[File]? subtree_assignments = glob("*subtree-assignments.tsv")  # !UnnecessaryQuantifier

		# distance matrices
		File? abig_matrix = "a000000.tsv"
		File? bbig_matrix = "b000000.tsv"
		Array[File]? acluster_matrices = glob("a*_dmtrx.tsv")  # !UnnecessaryQuantifier
		Array[File]? bcluster_matrices = glob("b*_dmtrx.tsv")  # !UnnecessaryQuantifier

		# cluster information
		File? unclustered_neighbors = "unclustered_neighbors.txt"
		#File rosetta_stone_20 = "rosetta_stone_20.tsv"
		#File rosetta_stone_10 = "rosetta_stone_10.tsv"
		#File rosetta_stone_5 = "rosetta_stone_5.tsv"
		#File nearest_and_furtherst_info = "all_neighbors.tsv"
		Int n_big_clusters = read_int("n_big_clusters")
		Int n_samples_in_clusters = read_int("n_samples_in_clusters")
		Int n_samples_processed = read_int("n_samples_processed")
		Int n_unclustered = read_int("n_unclustered")
		
		# old, maybe restore later?
		File latest_samples_temp = "latest_samples.tsv"
		#Array[File] abig_subtrees = glob("abig-subtree-*.nwk")
		File? samp_cluster = "samp_persiscluster" + today + ".tsv" # for nextstrain conversion
		#File? persistent_cluster_translator = "mapped_persistent_cluster_ids_to_new_cluster_ids.tsv"
		#Array[File] cluster_trees_json = glob("*.json")
		#Array[File] metadata_tsvs = glob("*.tsv")  # for auspice.us, which supports nwk
		#Array[String] unclustered_samples = read_lines("LONELY.txt")
	}
}

