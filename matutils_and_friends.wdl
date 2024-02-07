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
		
	    Int addldisk = 10
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
    }
    Int disk_size = ceil(size(input_mat, "GB")) + addldisk
	String output_mat = basename(input_mat, ".pb") + ".subtree" + ".pb"
	
	command <<<
	matUtils extract -i ~{input_mat} -s ~{samples} -o ~{output_mat}
	>>>
	
	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.0.2"
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
		docker: "ashedpotatoes/usher-plus:0.0.2"
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
		docker: "ashedpotatoes/usher-plus:0.0.2"
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
		i="/HOME/usher/example_tree/tb_alldiffs_mask2ref.L.fixed.pb"
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
		docker: "ashedpotatoes/usher-plus:0.0.2"
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
		docker: "ashedpotatoes/usher-plus:0.0.2"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File annotated_tree = outfile_annotated
	}

}

task convert_to_nextstrain_subtrees {
	# based loosely on Marc Perry's version
	input {
		File input_mat # aka tree_pb
		File? new_samples
		Int treesize = 0
		Int nearest_k = 250
		Int memory = 32
		Boolean new_samples_only
		String outfile_nextstrain = "nextstrain"
		Array[File?] metadata_files
	}

	String metadata = if defined(metadata_files) then "-M" else ""

	command <<<
		METAFILES_OR_EMPTY="~{sep=',' metadata_files}"
		if [[ "~{new_samples_only}" = "false" ]]
		then
			matUtils extract -i	~{input_mat} -S sample_paths.txt
			cut -f1 sample_paths.txt | tail -n +2 > sample.ids
			matUtils extract -i ~{input_mat} -j ~{outfile_nextstrain}.json -s sample.ids -N ~{treesize} ~{metadata} $METAFILES_OR_EMPTY
		else
			if [[ "~{new_samples}" == "" ]]
			then
				echo "Error -- new_samples_only is true, but no new_samples files was provided."
				exit 1
			else
				matUtils extract -i ~{input_mat} -j ~{outfile_nextstrain}.json -s ~{new_samples} -N ~{nearest_k} ~{metadata} $METAFILES_OR_EMPTY
			fi
		fi
		ls -lha
		
	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "yecheng/usher:latest"
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
		Array[File]? metadata_files
	}
	
	String metadata = if defined(metadata_files) then "-M" else ""

	command <<<
		METAFILES_OR_EMPTY="~{sep=',' metadata_files}"
		matUtils extract -i ~{input_mat} ~{metadata} $METAFILES_OR_EMPTY -j ~{outfile_nextstrain}
	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "yecheng/usher:latest"
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
		docker: "yecheng/usher:latest"
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
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
	}

	Int disk_size = ceil(size(diff, "GB")) + ceil(size(ref_genome, "GB")) +  ceil(size(input_mat, "GB")) + addldisk
	String D = if !(detailed_clades) then "" else "-D "
	#String ref = select_first([ref_genome, "/HOME/usher/ref/Ref.H37Rv/ref.fa"])
	#String i = select_first([input_mat, "/HOME/usher/example_tree/tb_alldiffs_mask2ref.L.fixed.pb"])

	command <<<
		if [[ "~{input_mat}" = "" ]]
		then
			i="/HOME/usher/example_tree/tb_alldiffs_mask2ref.L.fixed.pb"
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
		docker: "ashedpotatoes/usher-plus:0.0.2"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File usher_tree = output_mat
		File? clades = "clades.txt" # only if detailed_clades = true
	}

	meta {
		volatile: true
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
		Boolean only_matrix_special_samples
		File? special_samples
	}
	
	command <<<
	wget https://raw.githubusercontent.com/aofarrel/parsevcf/find-clusters-recursion/distancematrix_nwk.py
	if [[ "~{only_matrix_special_samples}" = "true" ]]
	then
		samples=$(< "~{special_samples}" tr -s '\n' ',' | head -c -1)
		echo "Samples that will be in the distance matrix: $samples"
		python3 distancematrix_nwk.py "~{input_nwk}" --samples "$samples" -v
	else
		python3 distancematrix_nwk.py "~{input_nwk}" -v
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
		Array[File] out_matrices = glob("*distance_matrix.tsv")
		File out_clusters = glob("*clusters.tsv")[0] # in auspice format
	}
	
}
