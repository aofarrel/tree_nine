version 1.0

# !! This calls the Google Cloud metadata server and WILL NOT WORK on other backends, including local !!

import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.29/tasks/processing_tasks.wdl" as processing
import "https://raw.githubusercontent.com/aofarrel/dropkick/1.0.0/dropkick.wdl" as dropkick

workflow Tree_Nine {
	input {
		Array[File] diffs
		File? input_tree
		File? existing_diffs
		File? existing_samples
		
		# matUtils/UShER options
		Boolean detailed_clades          = false
		Float?  max_low_coverage_sites

		# output file names and prefixes, extension not included
		String? comment
		Array[String]? rename_samples
		String out_prefix              = "bigtree"
		String out_diffs               = "_combined"
		
		# testing functions
		Boolean concat_files_then_exit = false
	}

	call processing.cat_files as cat_diff_files {
		input:
			new_files_to_concat = diffs,
			out_concat_file = out_prefix + out_diffs,
			keep_only_unique_lines = false,
			keep_only_unique_files = true,                             # STRICTLY NECESSARY UNLESS YOUR DATA *AND* SAMPLE IDS ARE DEDUPLICATED
			quality_report_removal_threshold = max_low_coverage_sites,
			out_sample_names = "samples_added",
			new_files_override_sample_names = rename_samples,
			king_file = existing_diffs,
			king_file_sample_names = existing_samples,
			and_then_exit_1 = concat_files_then_exit,
			datestamp_main_files = true,  # does not datestamp diffs
			out_concat_extension = ".diff"
	}
	String optional_datestamp           = cat_diff_files.today
	String outfile_usher_tree_raw       = out_prefix + optional_datestamp + "_raw.pb"

	call usher_sampled_diff as usher_sampled_diff_one {
		input:
			detailed_clades = detailed_clades,
			diff = cat_diff_files.outfile,
			input_mat = input_tree,
			output_mat = outfile_usher_tree_raw
	}

	call usher_sampled_diff as usher_sampled_diff_two {
		input:
			detailed_clades = detailed_clades,
			diff = cat_diff_files.outfile,
			input_mat = input_tree,
			output_mat = outfile_usher_tree_raw
	}
	
	output {
		String? out_comment = comment
		File  raw_tree_one = usher_sampled_diff_one.usher_tree
		File  raw_tree_two = usher_sampled_diff_two.usher_tree
	}
}


task usher_sampled_diff {
	input {
		# main files -- for TB, do not include ref_genome, it's already baked in!
		File diff
		File? input_mat
		File? ref_genome

		# usher options
		Int batch_size_per_process = 5
		Boolean detailed_clades
		Int optimization_radius = 0
		Int max_parsimony_per_sample = 1000000
		Int max_uncertainty_per_sample = 1000000
		String output_mat = basename(select_first([input_mat, "debugtree"]), ".pb") + "_new.pb"

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
		echo "----- CURRENT WORKDIR -----"
		ls -lha
		echo "---------- USHER ----------"

		usher-sampled ~{D} --optimization_radius=~{optimization_radius} \
			-e ~{max_uncertainty_per_sample} \
			-E ~{max_parsimony_per_sample} \
			--batch_size_per_process ~{batch_size_per_process} \
			--diff "~{diff}" \
			-i "$i" \
			--ref "$ref" \
			-o "~{output_mat}" >/dev/null 2>&1

		echo "--- CLOUD VM DIAGNOSTICS ---"
		apt-get update
		apt-get install -y curl

		echo -n "CPU Platform: "
		curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/cpu-platform || echo "Unknown"

		echo -n ""
		echo -n "Machine Type: "
		curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/machine-type | awk -F'/' '{print $NF}'

		lscpu | grep -E "Model name|flags" | head -n 2

		echo "------- DISK I/O TEST -------"
		# write 1 GB of 0s as a test of disk speed
		dd if=/dev/zero of=test_speed_file bs=1M count=1000 oflag=dsync 2>&1 | grep -E "copied|MB/s"
		rm test_speed_file
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