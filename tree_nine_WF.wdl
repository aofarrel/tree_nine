version 1.0

import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.1.18/tasks/processing_tasks.wdl" as processing
import "./matutils_and_friends.wdl"

# Anything marked !ForwardReference is using a bogus fallback value with select_first().

workflow Tree_Nine {
	input {
		Array[File] diffs
		
		# optional input - SNP distance matrix
		Boolean matrix_only_new_samples = false

		# optional inputs - filtering by coverage
		Array[File]? coverage_reports
		Float? max_low_coverage_sites
		
		# optional inputs - building trees
		Boolean detailed_clades = false
		File? input_tree                     # equivalent to UShER's i argument, if not defined, falls back to an SRA tree
		Boolean make_nextstrain_subtrees = true
		Boolean subtree_only_new_samples = true
		Boolean summarize_input_mat = true
		String? reroot_to_this_node          # equivalent to matUtils extract's y argument
		File? ref_genome                     # equivalent to USHER's ref argument
		File? metadata_tsv

		# output file names, extension not included
		Array[String]? rename_samples
		String out_prefix              = "tree"
		String out_prefix_summary      = out_prefix + "_"
		String in_prefix_summary       = basename(select_first([input_tree, "tb_alldiffs_mask2ref.L.fixed.pb"]))
		String out_diffs               = "_combined"
		String out_matrix              = "_matrix"
		String out_tree_annotated_pb   = "_annotated"
		String out_tree_nextstrain     = "_auspice"
		String out_tree_nwk            = "_nwk"
		String out_tree_taxonium       = "_taxonium"
		String out_tree_raw_pb         = "_raw"
		
	}

	parameter_meta {
		coverage_reports: "Single line text files generated by Lily's vcf to diff script, used to filter samples with low overall coverage"
		diffs: "Array of diff files"
		input_tree: "Input tree, equivalent to UShER's i argument"
		metadata_tsv: "TSV with one column of metadata"

		detailed_clades: "If true, run usher sampled diff with -D"
		make_nextstrain_subtrees: "If true, make nextstrain subtrees instead of one big nextstrain tree"
		matrix_only_new_samples: "If true, limit SNP distance matrix to only newly added samples"
		max_low_coverage_sites: "Maximum percentage of low coverage sites a sample can have before throwing it out"
		ref_genome: "Reference genome, equivalent to UShER's ref argument, default is H37Rv (M tuberculosis)"
		reroot_to_this_node: "Reroot the output tree relative to this node, leave blank to not reroot"
		out_prefix: "Prefix for all output files"
		subtree_only_new_samples: "If true and if make_nextstrain_subtrees true, nextstrain subtrees will only be focused on newly samples (ie samples added by your diffs)"
		summarize_input_mat: "If true and if an input tree is passed in, summarize that input tree"
	}

	call processing.cat_files as cat_diff_files {
		input:
			files = diffs,
			out_filename = out_prefix + out_diffs + ".diff",
			keep_only_unique_lines = false,
			removal_candidates = coverage_reports,
			removal_threshold = max_low_coverage_sites,
			first_lines_out_filename = "samples_added",
			overwrite_first_lines = rename_samples
	}

	File new_samples_added = select_first([cat_diff_files.first_lines, usher_sampled_diff.usher_tree]) #!ForwardReference

	if((summarize_input_mat)) {
		if (defined(input_tree)) {
			String basename_input_mat = basename(select_first([input_tree, ""]))

			# iff there is a metadata tsv, annotate input tree with it before summarizing
			if (defined(metadata_tsv)) {
				String annotated = "annotated_"
				call annotate as annotate_input_tree {
					input:
						input_mat = select_first([input_tree, usher_sampled_diff.usher_tree]), #!ForwardReference
						metadata_tsv = select_first([metadata_tsv, usher_sampled_diff.usher_tree]), #!ForwardReference
						outfile_annotated = "input_annotated_" + basename_input_mat + ".pb"
				}
			}

			File possibly_annotated_input_tree = select_first([annotate_input_tree.annotated_tree, input_tree])
			String annotated_or_blank = select_first([annotated, ""])

			call summarize as summarize_input_tree {
				input:
					input_mat = possibly_annotated_input_tree,
					prefix_outs = in_prefix_summary + annotated_or_blank
			}
		}
	}

	call usher_sampled_diff {
		input:
			detailed_clades = detailed_clades,
			diff = cat_diff_files.outfile,
			input_mat = input_tree,
			output_mat = out_prefix + out_tree_raw_pb + ".pb",
			ref_genome = ref_genome
	}

	if (defined(metadata_tsv)) {
		call annotate as annotate_output_tree {
			input:
				input_mat = usher_sampled_diff.usher_tree,
				metadata_tsv = select_first([metadata_tsv, usher_sampled_diff.usher_tree]), # bogus fallback
				outfile_annotated = out_prefix + out_tree_annotated_pb + ".pb"
		}
	}

	File possibly_annotated_output_tree = select_first([annotate_output_tree.annotated_tree, usher_sampled_diff.usher_tree])

	if(defined(reroot_to_this_node)) {

		call summarize as summarize_output_tree_before_reroot {
			input:
				input_mat = possibly_annotated_output_tree,
				prefix_outs = "before_reroot"
		}

		call reroot {
			input:
				input_mat = possibly_annotated_output_tree,
				reroot_to_this_node = select_first([reroot_to_this_node, ""])
		}
	}

	File final_output_tree = select_first([reroot.rerooted_tree, possibly_annotated_output_tree])

	# defined(metadata_tsv)   defined(reroot_to_this_node)          final_output_tree
	# ----------------------------------------------------------------------------------------------------------------
	#       true                       true                    annotated and rerooted
	#       true                      false                    annotated
	#       false                      true                    rerooted
	#       false                     false                    neither, just the output usher_sampled_diff.usher_tree

	call convert_to_newick {
		input:
			input_mat = final_output_tree,
			outfile_nwk = out_prefix + out_tree_nwk + ".nwk"
	}

	call convert_to_taxonium {
		input:
			input_mat = final_output_tree,
			outfile_taxonium = out_prefix + out_tree_taxonium + ".jsonl.gz"
	}

	if (make_nextstrain_subtrees) {
		call convert_to_nextstrain_subtrees {
			input:
				input_mat = final_output_tree,
				outfile_nextstrain = out_prefix + out_tree_nextstrain + ".json",
				new_samples = cat_diff_files.first_lines,
				new_samples_only = subtree_only_new_samples
		}
	}
	if (!make_nextstrain_subtrees) {
		call convert_to_nextstrain_single {
			input:
				input_mat = final_output_tree,
				outfile_nextstrain = out_prefix + out_tree_nextstrain + ".json"
		}
	}

	call matrix as make_distance_matrix {
		input:
			input_nwk = convert_to_newick.newick_tree,
			special_samples = new_samples_added,
			only_matrix_special_samples = matrix_only_new_samples,
			outfile_matrix = out_prefix + out_matrix + ".tsv"
	}

	call summarize as summarize_output_tree {
		input:
			input_mat = final_output_tree,
			prefix_outs = out_prefix_summary
	}

	output {
		# trees - protobuff
		#
		# note that tree_usher_rerooted is annotated if defined(metadata_tsv), but tree_usher_annotated is NOT rerooted
		# even if defined(reroot_to_this_node) -- this was done on purpose so people can get two annotated trees if they
		# want to easily compare the tree before and after rerooting
		File tree_usher_raw = usher_sampled_diff.usher_tree                                   # always
		File? tree_usher_rerooted = reroot.rerooted_tree                                      # iff defined(reroot_to_this_node)
		File? tree_usher_annotated = annotate_output_tree.annotated_tree                      # iff defined(metadata_tsv)

		# trees - other formats
		#
		# iff defined(reroot_to_this_node), these are based on usher_tree_rerooted
		# else, these are based on usher_tree_raw (and usher_tree_rerooted doesn't exist
		File tree_nwk = convert_to_newick.newick_tree                                         # always
		File tree_taxonium = convert_to_taxonium.taxonium_tree                                # always
		File? tree_nextstrain = convert_to_nextstrain_single.nextstrain_singular_tree         # mutually exclusive with nextstrain_subtrees
		Array[File]? subtrees_nextstrain = convert_to_nextstrain_subtrees.nextstrain_subtrees # mutually exclusive with nextstrain_tree

		# summaries
		File? summary_input = summarize_input_tree.summary                                    # iff summarize_input_mat
		File summary_output = summarize_output_tree.summary                                   # always
		File? summary_output_before_reroot = summarize_output_tree_before_reroot.summary      # iff defined(reroot_to_this_node)

		# sample information
		File? samples_input_tree = summarize_input_tree.samples                                # iff summarize_input_mat
		File samples_output_tree = summarize_output_tree.samples                               # always
		File? samples_output_tree_before_reroot = summarize_output_tree_before_reroot.samples  # iff defined(reroot_to_this_node)
		Array[String] samples_added = read_lines(new_samples_added)                            # always
		Array[String] samples_dropped = cat_diff_files.removed_files                           # always
		File distance_matrix = make_distance_matrix.out_matrix                                 # always
	}
}