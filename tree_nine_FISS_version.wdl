version 1.0

# This version of Tree Nine is almost identical to the main one. It is designed for use with the FISS API.
# Unless you intend on running Tree Nine in an automated fashion, you should use tree_nine.wdl, not this!

import "https://raw.githubusercontent.com/aofarrel/SRANWRP/v1.3.1/tasks/processing_tasks.wdl" as processing
import "https://raw.githubusercontent.com/aofarrel/dropkick/1.2.0/dropkick.wdl" as dropkick
import "https://raw.githubusercontent.com/aofarrel/microreact_WDLs/1.0.0/share_projects_with_team_via_file.wdl"
import "https://raw.githubusercontent.com/aofarrel/diffdiff/0.2.2/diffdiff.wdl" as diffdiff

# User notes:
# * If user doesn't define input_tree, a rudimentary 7K sample tree will serve as the base tree. This base tree
#   doesn't represent genetic diversity of MTBC well and should not be used for anything besides quick testing.
# * Should be run with --copy-input-files on miniwdl (required if clustering, may work w/o it if not clustering)

# Dev notes:
# * Anything marked !ForwardReference is using a bogus fallback value with select_first() to coerce to not-optional

workflow Tree_Nine {
	input {
		Array[File] diffs

		# Most important "modes" of running Tree Nine
		Boolean adhoc = false
		Boolean identify_clusters = false
		Boolean restart_clusters = false # WARNING: Will generate brand new cluster IDs and Microreact Projects,
										 # but Tree Nine intentionally CANNOT DELETE EXISTING MICROREACT PROJECTS.
										 # If you need to delete MR projects, use the API, or this WDL:
										 # https://github.com/aofarrel/microreact_WDLs/blob/main/delete_project.wdl
		Boolean upload_clusters_to_microreact  = false

		# Recommendation: Use the same base tree every time (do not pass in previous Tree Nine run's BIG_tree_usher)
		# Regardless, do not leave this undefined unless doing very quick tests; the fallback base tree includes low-quality samples
		File? input_tree

		# Recommendation: Pass in previous Tree Nine run's updated_diff_file and updated_diff_contents
		File? existing_diffs
		File? existing_samples

		String? listener_bucket
		
		# matUtils/UShER options
		Boolean detailed_clades          = false
		Float?  max_low_coverage_sites
		File?   matutils_clade_annotations
		Boolean optimize = true
		String? reroot_to_this_node
		Boolean summarize_tree_before_placing_samples   = false 
		Boolean summarize_tree_after_placing_samples    = false

		# Options related to clustering/distance matrix
		Boolean cluster_entire_tree            = false  # strongly recommended to leave as false or else crashing is likely
		File? cluster_these_samples_override
		File? cluster_ids_to_never_generate

		# metadata file; expected to be pulled via the FISS API but not strictly required
		File? sample_metadata_tsv
		Boolean strictly_check_metadata = true

		# if you are running with pre-existing clusters, all three of these must be filled in
		# if you are identifying clusters ad-hoc, all three of these must be undefined
		File? persistent_cluster_meta    # vital for persistent clusters
		File? persistent_cluster_ids     # vital for persistent clusters
		File? previous_run_cluster_json  # only used to generate a change report but required by input validation
		
		# related to putting clusters on Microreact
		File? microreact_blank_template_json
		File? microreact_decimated_template_json
		File? microreact_key
		File? microreact_update_template_json

		# non-exclusive ways of sharing your microreact projects
		String? microreact_share_email
		String? microreact_share_team
		
		# rarely used files (see parameter_meta)
		Array[File]? coverage_reports
		File? ref_genome               # do not define this if you're using H37Rv!
		
		# output file names and prefixes, extension not included
		String? comment
		Array[String]? rename_samples
		Boolean datestamp_outs         = true
		String out_prefix              = "bigtree"
		String out_diffs               = "_combined"
		
		# testing functions
		Boolean DEBUG_concat_files_then_exit = false
		File?   DEBUG_override_latest_samples
		File?   DEBUG_override_latest_clusters
		Boolean DEBUG_generate_debug_mr_jsons = false
	}

	parameter_meta {
		comment: "String that gets copied directly to output (useful for Terra data tables)"
		
		diffs: "Array of single-sample MAPLE-formatted diff files from myco"

		input_tree: "The base MAT tree (.pb) that samples will be placed upon; will fall back to a test tree on ~7K TB samples from SRA if not defined (test tree should ONLY be used for quick debugging; it is not representative of MTBC diversity nor proper sample QC)"

		existing_diffs: "A pre-concatenated multi-sample .diff file which will be concatenated with the single-sample diffs input (requires existing_samples)"
		
		existing_samples: "A file listing all samples within the multiple-sample existing_diffs file (requires existing_diffs)"
		
		matutils_clade_annotations: "Two column TSV for clade annotation via matUtils"
		
		cluster_entire_tree: "If true, matrix and cluster all samples on tree; if false, only matrix and cluster cluster_these_samples_override (if defined) or newly added samples."
		
		cluster_max_distance: "Soft-maximum SNP distance between two samples for them to be in the same cluster. NOTE if cluster_max_distance=10, A:B=10, B:C=5, and A:C=15, then all three will still be in a cluster even though A:C is above cluster_max_distance, since both are within 10 of another sample in that cluster."
		
		coverage_reports: "NOT USUALLY NEEDED - Single line text files generated by Lily's vcf to diff script, used to filter samples with low overall coverage. By default, myco's vcf_to_diff.py filters sites per site coverage, so this isn't typically needed here in Tree Nine."
		
		detailed_clades: "usher_sampled_diff -D"
		
		max_low_coverage_sites: "Maximum percentage of low coverage sites a sample can have before throwing it out (requires coverage_reports, does not apply to backmasked diffs)"
		
		cluster_these_samples_override: "Provide an override file containing names of the only samples to consider for matrix and clustering. If this isn't defined, matrixing and clustering is done on either entire tree (if cluster_entire_tree) or all samples with a diff file (if not cluster_entire_tree)."

		cluster_ids_to_never_generate: "Newline delimited text file of cluster IDs to never generate. For example, if you don't want a cluster to be assigned the ID 000013 because that feels unlucky, add 000013 to this file. Does not affect cluster IDs that already exist. If you are just trying to track cluster IDs persistently without them being reassigned, don't worry about this input, focus on the carryover files instead."
		
		ref_genome: "Reference genome, equivalent to UShER's ref argument, default is H37Rv (M tuberculosis)"
		
		rename_samples: "For file at index i in diffs[i], rename it to the corresponding string at rename_samples[i]."
		
		reroot_to_this_node: "matUtils extract -y (Reroot the output tree relative to this node, leave blank to not reroot)"
		
		out_prefix: "Prefix for all output files"
		
		upload_clusters_to_microreact: "If you know, you know"
	}

	# Metadata notes:
	#
	# 1) Metadata columns are currently hardcoded as they need to be in the Microreact template too. Newer versions of the clustering script attempt
	# to handle this on the fly, but just to be safe...
	# 2) Previously, the plan was to build Array[Pair[String, String]] "dictionaries" to replace TBProfiler lineage descriptions, and rename columns,
	# per CDPH request. However, due to https://github.com/broadinstitute/cromwell/issues/7883, I cannot actually pass Array[Pair[String, String]] into 
	# process_metadata without Cromwell crashing, even though it works perfectly on miniwdl, and even though womtool (Cromwell's checker, includes type
	# checking) does not have any issues with Array[Pair[String, String]]. As such, I am skipping the requested "La1.2 -> BCG" etc renames, and hardcoding
	# column names within process_CDPH_clusters.
	# 3) Due to how FISS works, if you download the sample level data table via FISS and then re-upload it to create your metadata table, FISS will probably
	# drop any columns that are 100% null. You'll want to make sure all columns are present before reuploading.
	Array[String]? microreact_metadata_columns = ["Epi_Duplication","Year_Collected","Patient_County","State","Country","20_Cluster_Date","10_Cluster_Date","5_Cluster_Date","tbd_strain_per_tbprof","tbd_resistance","Submitter_Facility","Submitter_Facility_Sample_ID","Sequencing_Facility","Latitude","Longitude"] #!UnnecessaryQuantifier
	

	call validate_treenine_inputs as validate_inputs {
		input:
			input_tree = input_tree,
			existing_diffs = existing_diffs,
			existing_samples = existing_samples,
			persistent_cluster_meta = persistent_cluster_meta,
			persistent_cluster_ids = persistent_cluster_ids,
			previous_run_cluster_json = previous_run_cluster_json,
			microreact_blank_template_json = microreact_blank_template_json,
			microreact_decimated_template_json = microreact_decimated_template_json,
			microreact_key = microreact_key,
			microreact_update_template_json = microreact_update_template_json,
			ref_genome = ref_genome,
			DEBUG_generate_debug_mr_jsons = DEBUG_generate_debug_mr_jsons,
			identify_clusters = identify_clusters,
			upload_clusters_to_microreact = upload_clusters_to_microreact,
			restart_clusters = restart_clusters,
			adhoc = adhoc
	}

	if (adhoc) {
		call diffdiff.diffdiff_usher_mask as diffdiff_usher {
			input:
				diffs = diffs
		}
	}

	if (defined(sample_metadata_tsv)) {

		# Surely if you reference an optional value within a defined() block, you can input it as a non-optional, right?
		# Not so! You will see here (and elsewhere) bogus select_first() fallbacks because of this quirk of WDL.

		call processing.process_metadata_table as process_metadata {
			input:
				table = select_first([sample_metadata_tsv, diffs[0]]),
				desired_columns = microreact_metadata_columns,
				strict = strictly_check_metadata
		}
	}

	call processing.cat_files as cat_diff_files {
		input:
			new_files_to_concat = diffs,
			out_concat_file = out_prefix + out_diffs,
			keep_only_unique_lines = false,
			keep_only_unique_files = true,                             # STRICTLY NECESSARY UNLESS YOUR DATA *AND* SAMPLE IDS ARE DEDUPLICATED
			new_files_quality_reports = coverage_reports,
			quality_report_removal_threshold = max_low_coverage_sites,
			out_sample_names = "samples_added",
			new_files_override_sample_names = rename_samples,
			king_file = existing_diffs,
			king_file_sample_names = existing_samples,
			and_then_exit_1 = DEBUG_concat_files_then_exit,
			datestamp_main_files = true,  # does not datestamp diffs
			out_concat_extension = ".diff"
	}

	File samples_considered_for_clustering = select_first([cluster_these_samples_override, cat_diff_files.first_lines, usher_sampled_diff.usher_tree]) #!ForwardReference

	# Tree Nine attempts to use a clear naming scheme to make its large number of output files unambigious, but you might have a better
	# system than I do, so I'm going to define all remaining major outfile-controlling variables here so you can edit it easily.

	String empty_string = ""
	if(!(datestamp_outs)) { String no_datestamp = "" }
	String optional_datestamp = select_first([no_datestamp, cat_diff_files.today])
	String presumed_input_mat_basename  = basename(select_first([input_tree, "default_debug-only_basetree"]))
	
	String outfile_annotated_input_tree = "input_" + presumed_input_mat_basename + optional_datestamp + ".pb"
	String outfile_usher_tree_raw       = out_prefix + optional_datestamp + "_raw.pb"
	String outfile_usher_tree_optimized = basename(outfile_usher_tree_raw, "_raw.pb") + "_optimized.pb"
	String outfile_usher_tree_annotated = basename(outfile_usher_tree_raw, "_raw.pb") + "_annotated.pb"
	String outfile_usher_tree_rerooted  = basename(outfile_usher_tree_raw, "_raw.pb") + "_reroot_to_" + select_first([reroot_to_this_node, empty_string]) + ".pb"
	String outfile_taxonium_tree        = basename(outfile_usher_tree_raw, "_raw.pb") + "_taxonium.jsonl.gz"
	String outfile_nextstrain_tree      = basename(outfile_usher_tree_raw, "_raw.pb") + ".json"
	#String outfile_nwk_matutils_tree    = basename(outfile_usher_tree_raw, "_raw.pb") + ".nwk"
	# There is also a nwk tree generated by the clustering script, that one is always called a000000.nwk and should be identical to the matutils one
	
	String outfile_input_tree_summaries               = "input_" + presumed_input_mat_basename + "_"
	String outfile_usher_tree_summaries_before_reroot = out_prefix + optional_datestamp + "_before_reroot_"
	String outfile_usher_tree_summaries_final         = out_prefix + optional_datestamp + "_final_"

	
	if(summarize_tree_before_placing_samples) {
		if (defined(input_tree)) {

			# iff there is a metadata tsv, annotate input tree with it before summarizing
			if (defined(matutils_clade_annotations)) {

				call annotate as annotate_input_tree {
					input:
						input_mat = select_first([input_tree, usher_sampled_diff.usher_tree]), #!ForwardReference
						metadata_tsv = select_first([matutils_clade_annotations, usher_sampled_diff.usher_tree]), #!ForwardReference
						outfile_mat = outfile_annotated_input_tree
				}
			}

			File possibly_annotated_input_tree = select_first([annotate_input_tree.annotated_tree, input_tree])

			call summarize as summarize_input_tree {
				input:
					input_mat = possibly_annotated_input_tree,
					prefix_outs = outfile_input_tree_summaries
			}
		}
	}

	call usher_sampled_diff as usher_sampled_diff {
		input:
			detailed_clades = detailed_clades,
			diff = cat_diff_files.outfile,
			input_mat = input_tree,
			output_mat = outfile_usher_tree_raw,
			ref_genome = ref_genome,
			noop_boolean = validate_inputs.didnt_crash
	}

	if (optimize) {
		call matOptimize as matOptimize_usher {
			input:
				input_mat = usher_sampled_diff.usher_tree,
				output_mat = outfile_usher_tree_optimized
		}
	}

	File optimized_or_raw_tree = select_first([matOptimize_usher.optimized_tree, usher_sampled_diff.usher_tree])
	

	if (defined(matutils_clade_annotations)) {
		call annotate as annotate_usher {
			input:
				input_mat = optimized_or_raw_tree,
				metadata_tsv = select_first([matutils_clade_annotations, usher_sampled_diff.usher_tree]), # bogus fallback
				outfile_mat = outfile_usher_tree_annotated
		}
	}

	File possibly_annotated_maximal_output_tree = select_first([annotate_usher.annotated_tree, optimized_or_raw_tree])

	if(defined(reroot_to_this_node)) {

		if(summarize_tree_after_placing_samples) {
			call summarize as summarize_before_reroot {
				input:
					input_mat = possibly_annotated_maximal_output_tree,
					prefix_outs = outfile_usher_tree_summaries_before_reroot
			}
		}

		call reroot as reroot_usher {
			input:
				input_mat = possibly_annotated_maximal_output_tree,
				reroot_to_this_node = select_first([reroot_to_this_node, ""]),
				output_mat = outfile_usher_tree_rerooted
		}
	}

	File final_maximal_output_tree = select_first([reroot_usher.rerooted_tree, possibly_annotated_maximal_output_tree])

	# defined(matutils_clade_annotations)   defined(reroot_to_this_node)          final_maximal_output_tree
	# ----------------------------------------------------------------------------------------------------------------
	#       true                       true                    annotated and rerooted
	#       true                      false                    annotated
	#       false                      true                    rerooted
	#       false                     false                    neither, just the output matOptimize_usher.optimized_tree

	call convert_to_taxonium as to_taxonium {
		input:
			input_mat = final_maximal_output_tree,
			outfile_taxonium = outfile_taxonium_tree
	}

	if (identify_clusters) {
		if (!defined(DEBUG_override_latest_samples)) {
			call find_CDPH_clusters as find_clusters {
				input:
					input_mat_with_new_samples = final_maximal_output_tree,
					special_samples = samples_considered_for_clustering,
					only_matrix_special_samples = !(cluster_entire_tree),
					datestamp = cat_diff_files.today
			}
		}
		
		call process_CDPH_clusters as process_clusters {
			input:
				shareemail = microreact_share_email,
				input_mat_with_new_samples = final_maximal_output_tree,
				special_samples = samples_considered_for_clustering,
				combined_diff_file = cat_diff_files.outfile,
				only_matrix_special_samples = !(cluster_entire_tree),
				persistent_ids = persistent_cluster_ids,
				persistent_cluster_meta = persistent_cluster_meta,
				previous_run_cluster_json = previous_run_cluster_json,
				microreact_key = microreact_key,
				microreact_update_template_json = microreact_update_template_json,
				microreact_blank_template_json = microreact_blank_template_json,
				microreact_decimated_template_json = microreact_decimated_template_json,
				persistent_denylist = cluster_ids_to_never_generate,
				upload_clusters_to_microreact = upload_clusters_to_microreact,
				datestamp = cat_diff_files.today,
				sample_metadata_tsv = process_metadata.processed_metadata_table,
				microreact_metadata_columns = microreact_metadata_columns,
				latest_samples_tsv = select_first([find_clusters.latest_samples_tsv, DEBUG_override_latest_samples]),
				latest_clusters_tsv = select_first([find_clusters.latest_clusters_tsv, DEBUG_override_latest_clusters]),
				cluster_matrices_randomIDs_tarball = find_clusters.cluster_matrices_randomIDs,
				cluster_subtrees_randomIDs_tarball = find_clusters.cluster_subtrees_randomIDs,
				DEBUG_generate_debug_mr_jsons = DEBUG_generate_debug_mr_jsons
		}

		# This is some trickery to prevent Cromwell from complaining about us putting an "optional" output
		# into a non-optional task input. We can do this because the "optional" output actually gets created
		# in non-error cases (at least, this is the case with how we call the task here in Tree Nine)
		# so it will never actually fall back on optimized_or_raw_tree -- and if it does error, the entire
		# pipeline crashes so what happens here is moot.
		File coerced_unclustered_txt = select_first([find_clusters.unclustered_samples, optimized_or_raw_tree])

		if (defined(listener_bucket)) {
			# More coercion workarounds here, this one is even sillier because we're in a defined() block, alas
			# this is required.
			String coerced_destination_bucket = select_first([listener_bucket, "nonsense fallback value"])
			call dropkick.Dropkick_Curl as upload_cluster_json {
				input:
					destination_bucket = coerced_destination_bucket,
					files_to_upload = [process_clusters.final_cluster_information_json]
			}

			call dropkick.Dropkick_Curl as upload_unclustered_txt {
				input:
					destination_bucket = coerced_destination_bucket,
					files_to_upload = [coerced_unclustered_txt]
			}
		}

		if (defined(microreact_share_team)) {
			if (defined(microreact_key)) {
				if (defined(process_clusters.updated_mr_URIs_file)) { # must explictly check as it could be undefined if nothing got updated 
					File coerced_microeract_key = select_first([microreact_key, optimized_or_raw_tree])
					String coerced_microreact_share_team = select_first([microreact_share_team, "nonsense fallback value"])
					File coerced_updated_mr_URIs_file = select_first([process_clusters.updated_mr_URIs_file, optimized_or_raw_tree])
					call share_projects_with_team_via_file.Microreact_Share_Projects_With_Team {
						input:
							token = coerced_microeract_key,
							team_uri = coerced_microreact_share_team,
							project_uris = coerced_updated_mr_URIs_file
					}
				}
			}
		}

		call convert_to_nextstrain_single_terra_compatiable as to_nextstrain_cluster {
			input:
				input_mat = final_maximal_output_tree,
				outfile_nextstrain = outfile_nextstrain_tree,
				one_metadata_file = process_clusters.samp_cluster_ten
		}
	}

	if (!(identify_clusters)) {
		call convert_to_nextstrain_single_terra_compatiable as to_nextstrain {
				input:
					input_mat = final_maximal_output_tree,
					outfile_nextstrain = outfile_nextstrain_tree
			}
	}
	

	if(summarize_tree_after_placing_samples) {
		call summarize as summarize_final {
			input:
				input_mat = final_maximal_output_tree,
				prefix_outs = outfile_usher_tree_summaries_final
		}
	}
	
	output {
		String? out_comment = comment
		File?   unclusted_samples = find_clusters.unclustered_samples

		# big trees - protobuff
		#
		# note that tree_usher_rerooted is annotated if defined(matutils_clade_annotations), but tree_usher_annotated is NOT rerooted
		# even if defined(reroot_to_this_node) -- this was done on purpose so people can get two annotated trees if they
		# want to easily compare the tree before and after rerooting
		#
		File  BIG_tree_usher = optimized_or_raw_tree
		File  BIG_tree_usher_raw_dont_use = usher_sampled_diff.usher_tree
		File? BIG_tree_reroot = reroot_usher.rerooted_tree
		File? BIG_tree_ushanno = annotate_usher.annotated_tree

		# big trees - other formats
		#
		# iff defined(reroot_to_this_node), these are based on usher_tree_rerooted
		# else, these are based on usher_tree_raw (and usher_tree_rerooted doesn't exist)
		#
		File?  BIG_tree_nwk_raw = find_clusters.bigtree_raw
		File?  BIG_tree_nwk_gen = find_clusters.bigtree_gen
		File   BIG_tree_taxonium = to_taxonium.taxonium_tree
		File?  BIG_tree_json_noanno = to_nextstrain.nextstrain_singular_tree
		File?  BIG_tree_json_clusteranno = to_nextstrain_cluster.nextstrain_singular_tree

		# cluster subtrees/matrices -- the ones from process_clusters have the expected persistent IDs, the ones
		# from find_clusters DO NOT, so we'll only include the process_clusters tarballs here. the reason why
		# these are tarballs now is because Terra has a hard limit on the number of workflow outputs.
		File? cluster_subtrees = process_clusters.cluster_trees_persisIDs                          # formerly Array[File] CLUSTER_trees_nwk
		File? cluster_matrices = process_clusters.cluster_matrices_persisIDs                       # formerly Array[File] CLUSTER_dmatrices
		File? cluster_matrices_backmasked= process_clusters.cluster_matrices_persisIDs             # formerly Array[File] BM_CLUSTER_dmatrices
		File? cluster_subtrees_backmasked = process_clusters.cluster_matrices_persisIDs_backmasked # formerly Array[File] BM_CLUSTER_trees_nwk

		# if you want to do persistent clustering in the future, you need all five of these files
		File updated_diff_file = cat_diff_files.outfile
		File updated_diff_contents = samples_considered_for_clustering
		File? updated_persistent_ids = process_clusters.new_persistent_ids
		File? updated_persistent_meta = process_clusters.new_persistent_meta
		File? updated_cluster_information_json = process_clusters.final_cluster_information_json

		# diffdiff outputs
		File? diffdiff_full_alignment = diffdiff_usher.full_alignment
        File? diffdiff_noteworthy_alignment = diffdiff_usher.noteworthy_alignment
        File? diffdiff_usher_mask = diffdiff_usher.usher_mask

		#### "stats for the nerds" section, most users don't need these but they're good context ####

		# cat_diff_files
		Int   n_new_samps_input = cat_diff_files.files_input
		Int   n_new_samps_skipped = cat_diff_files.files_removed
		Array[String] samples_dropped = cat_diff_files.removed_files

		# cluster-related
		File? BIG_matrix_nb = find_clusters.bigtree_matrix    # nb as in "not backmasked" although there is no backmasked version
		File? all_samples_nearest_relatives = find_clusters.all_nearest_relatives
		File? all_samples_that_clustered = process_clusters.all_samples_cluster_information
		File? new_samples_that_clustered = process_clusters.new_samples_cluster_information
		Int?  n_20SNP_clusters = find_clusters.n_big_clusters
		Int?  n_samps_unclustered = find_clusters.n_unclustered
		Int?  n_samps_clustered = find_clusters.n_samples_in_clusters
		Int?  n_samps_processed = find_clusters.n_samples_processed
		File? unclustered_subtrees_and_info = find_clusters.unclustered_subtrees_etc
		File? mr_uris_updated = process_clusters.updated_mr_URIs_file  # awkward name because not required for subsequent runs
		
		# tree summary tasks
		File? in_summary = summarize_input_tree.summary
		File? nb_summary_preroot = summarize_before_reroot.summary      # iff defined(reroot_to_this_node)
		File? nb_summary_final = summarize_final.summary
		File? in_list_samples = summarize_input_tree.samples
		File? nb_list_samples_preroot = summarize_before_reroot.samples # iff defined(reroot_to_this_node)
		File? nb_list_samples_final = summarize_final.samples
		

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
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
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
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
		memory: memory + " GB"
		preemptible: 1
	}

	output {
		File nextstrain_singular_tree = outfile_nextstrain
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
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File taxonium_tree = outfile_taxonium
	}
}

task reroot {
	input {
		File input_mat
		String reroot_to_this_node
		String output_mat = basename(input_mat, ".pb") + ".reroot_to_~{reroot_to_this_node}" + ".pb"

		Int addldisk = 10
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
	}
	Int disk_size = ceil(size(input_mat, "GB")) + addldisk

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
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
		memory: memory + " GB"
		preemptible: preempt
	}
	
	output {
		File rerooted_tree = output_mat
	}
}

task annotate {
	input {
		File? input_mat
		File metadata_tsv # only can annotate one column at a time
		String outfile_mat

		Int addldisk = 10
		Int cpu = 8
		Int memory = 16
		Int preempt = 1
	}
	Int disk_size = ceil(size(input_mat, "GB")) + ceil(size(metadata_tsv, "GB")) + addldisk

	command <<< 
	matUtils annotate -i "~{input_mat}" -P "~{metadata_tsv}" -o "~{outfile_mat}"
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File annotated_tree = outfile_mat
	}
}

task matOptimize {
	input {
		File input_mat
		Int max_hours = 1
		Int? max_iterations 
		Float min_improvement = 0.00000001
		String output_mat = basename(input_mat, ".pb") + "_optimized.pb"

		Int addldisk = 100
		Int cpu = 24
		Int memory = 32
		Int preempt = 1
	}

	Int disk_size = ceil(size(input_mat, "GB")) + addldisk

	command <<<
	if [[ ! "~{max_iterations}" = "" ]]
	then
		MAX_ITERATIONS="--max-iterations ~{max_iterations}"
	else
		MAX_ITERATIONS=""
	fi
	# shellcheck disable=SC2086
	matOptimize -i "~{input_mat}" --max-hours ~{max_hours} --min-improvement ~{min_improvement} $MAX_ITERATIONS -o "~{output_mat}"
	>>> 

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File optimized_tree = output_mat
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

		# No-op to force this task to run downstream of the validate_inputs task in Tree Nine
		Boolean noop_boolean = true #!UnusedDeclaration

		# Prevent "got unrecognized trialing" logs which grinds GCP to a halt
		Boolean silence_usher = true
	}

	Int disk_size = ceil(size(diff, "GB")) + ceil(size(ref_genome, "GB")) +  ceil(size(input_mat, "GB")) + addldisk
	String D = if !(detailed_clades) then "" else "-D "

	command <<<
		if [[ "~{input_mat}" = "" ]]
		then
			i="/HOME/ash/example_tree/for_debugging_only__tb_7K_noQC_diffs_mask2ref.L.fixed.pb"
		else
			i="~{input_mat}"
		fi

		if [[ "~{ref_genome}" = "" ]]
		then
			ref="/HOME/ash/ref/Ref.H37Rv/ref.fa"
		else
			ref="~{ref_genome}"
		fi
		
		echo "input mat: ~{input_mat}"
		echo "mat we will use: $i"
		echo "reference genome per user: ~{ref_genome}"
		echo "reference genome we will use: $ref"
		echo "------------------"
		tree
		echo "------------------"

		if [[ ~{silence_usher} = true ]]
		then
			echo "WARNING: All UShER prints will be silenced"
			usher-sampled ~{D} --optimization_radius=~{optimization_radius} \
				-e ~{max_uncertainty_per_sample} \
				-E ~{max_parsimony_per_sample} \
				--batch_size_per_process ~{batch_size_per_process} \
				--diff "~{diff}" \
				-i "$i" \
				--ref "$ref" \
				-o "~{output_mat}" >/dev/null 2>&1
		else
			usher-sampled ~{D} --optimization_radius=~{optimization_radius} \
				-e ~{max_uncertainty_per_sample} \
				-E ~{max_parsimony_per_sample} \
				--batch_size_per_process ~{batch_size_per_process} \
				--diff "~{diff}" \
				-i "$i" \
				--ref "$ref" \
				-o "~{output_mat}"
		fi

		tree
	>>>

	runtime {
		cpu: cpu
		disks: "local-disk " + disk_size + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File usher_tree = output_mat
		File? clades = "clades.txt"                   # only if detailed_clades = true
		File? mutation_paths = "mutation-paths.txt"
		File? placement_stats = "placement_stats.tsv"
	}

}

task find_CDPH_clusters {
	# Formally cluster_CDPH_method, now split into two tasks for easier debugging
	# find_clusters.py: Generates 20-10-5 clusters and distance matrices (normal and backmasked)
	# This might not work properly if any sample IDs contain a space
	input {
		File input_mat_with_new_samples
		String datestamp # has to be defined here for non-glob delocalization to work properly

		# If not provided, it is assumed you want to matrix the ENTIRE tree
		Boolean only_matrix_special_samples
		File? special_samples

		Int memory = 50

		# EXPERIMENTAL: Store the distance matrix in memory as eight-bit unsigned integers. Actual calculations are
		# done in 64 bit and anything that would overflow is set to 255. This can resolve out-of-memory issues on
		# limited hardware, but it's not recommended.
		Boolean inteight = false
		
		# these should only be set for test runs/debugging
		File?   override_find_clusters_script
		Int     preempt = 0
		
	}

	Array[Int] cluster_distances = [20, 10, 5] # CHANGING THIS MIGHT BREAK THINGS!
	String arg_ieight = if inteight then "--int8" else ""
	
	command <<<
		set -eux pipefail
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting task"
		
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting newick (A_big.nwk) from input PB"
		matUtils extract -i ~{input_mat_with_new_samples} -t A_big.nwk
		cp ~{input_mat_with_new_samples} .

		if [[ ! "~{override_find_clusters_script}" == '' ]]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Overwriting find_clusters.py script with user-defined input"
			rm /HOME/ash/scripts/find_clusters.py
			mv "~{override_find_clusters_script}" /HOME/ash/scripts/find_clusters.py
		fi

		CLUSTER_DISTANCES="~{sep=',' cluster_distances}"
		FIRST_DISTANCE="${CLUSTER_DISTANCES%%,*}"
		OTHER_DISTANCES="${CLUSTER_DISTANCES#*,}"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cluster distances: $CLUSTER_DISTANCES"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] First distance: $FIRST_DISTANCE"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Other distances: $OTHER_DISTANCES"

		# Turn off pipefail because find_clusters.py can return not-0 in non-error cases
		# TODO: this isn't great practice; there's better ways to handle the recursion;
		# should probably make the matrix generator its own script
		set +eo pipefail 

		# TODO: on very large runs, the size of $/samples may eventually cause issues with ARG_MAX
		# should be fine for our purposes though

		# shellcheck disable=SC2086
		if [[ "~{only_matrix_special_samples}" = "true" ]]
		then
			samples=$(< "~{special_samples}" tr -s '\n' ',' | head -c -1)
			echo "Samples that will be in the distance matrix: $samples"
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running find_clusters.py"
			python3 /HOME/ash/scripts/find_clusters.py \
				"~{input_mat_with_new_samples}" \
				--samples $samples \
				--collection-name big \
				-t NB \
				-d "$FIRST_DISTANCE" \
				-rd "$OTHER_DISTANCES" \
				-v ~{arg_ieight}
		else
			echo "No sample selection file passed in, will matrix the entire tree (WARNING: THIS MAY BE VERY SLOW)"
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running find_clusters.py"
			python3 /HOME/ash/scripts/find_clusters.py \
				"~{input_mat_with_new_samples}" \
				--collection-name big \
				-t NB \
				-d "$FIRST_DISTANCE" \
				-rd "$OTHER_DISTANCES" \
				-v ~{arg_ieight}
		fi
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished running find_clusters.py"

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Contents of workdir before processing outputs:"
		tree

		# A_big.nwk									big tree, nwk format (will be renamed later)
		# LONELY-subtree-n.nwk (n as variable)		subtrees (usually multiple) of unclustered samples
		# unclustered_samples.txt					what it says on the tin
		# lonely-subtree-assignments.tsv			which subtree each unclustered sample ended up in
		# cluster_annotation_workdirIDs.tsv			can be used to annotate by nonpersistent cluster (but isn't, at least not yet)
		# latest_samples.tsv						used by persistent ID script (will be renamed later)
		# n_big_clusters (n as constant)			# of 20SNP clusters
		# n_samples_in_clusters (n as constant)		# of samples that clustered
		# n_samples_processed (n as constant)		# of samples processed by find_clusters.py
		# n_unclustered (n as constant)				# of samples that failed to cluster
		# ...and one distance matrix per cluster, and also one(?) subtree per cluster. Later, there will be two of each per cluster thanks to backmasking

		mv A_big.nwk "BIGTREE~{datestamp}.nwk"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed A_big.nwk to BIGTREE~{datestamp}.nwk"
		mv all_closest_relatives.txt "all_nearest_relatives~{datestamp}.txt"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed all_closest_relatives.txt to all_nearest_relatives~{datestamp}.txt"
		mv latest_samples.tsv "latest_samples~{datestamp}.tsv"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed latest_samples.tsv to latest_samples~{datestamp}.tsv"
		mv latest_clusters.tsv "latest_clusters~{datestamp}.tsv"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed latest_clusters.tsv to latest_clusters~{datestamp}.tsv"
		mv unclustered_samples.txt "unclustered_samples~{datestamp}.txt"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Renamed unclustered_samples.txt to unclustered_samples~{datestamp}.txt"

		# copy stuff that will go into an archive but we also want as a task-level output
		cp aworkdir000000.nwk bigtree_gen.temp
		cp "unclustered_samples~{datestamp}.txt" nonclustered_samples.temp

		find . -maxdepth 1 \( -name "LONELY*.nwk" -o -name "lonely-subtree-assignments.tsv" -o -name "unclustered*.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > unclustered_subtrees_etc.tar.gz

		echo "The IDs of these clusters are random and DO NOT account for persistent cluster IDs. " > readme.txt
		echo "You'll need to run process_clusters.py to get your persistent cluster IDs!" >> readme.txt

		find . -maxdepth 1 \( -name "a*.nwk" -o -name "a*.pb" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > randomID_cluster_trees.tar.gz
		find . -maxdepth 1 \( -name "a*_dmtrx.tsv" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > randomID_cluster_matrices.tar.gz

		# for output matching (will include "workdir" for consistency)
		mv bigtree_gen.temp aworkdir000000.nwk
		mv nonclustered_samples.temp "unclustered_samples~{datestamp}.txt"

		# shellcheck disable=SC2317
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Contents of workdir after processing outputs:"
		tree
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished finding clusters"

	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		File all_nearest_relatives    = "all_nearest_relatives" + datestamp + ".txt" # includes unclustered + clustered samples
		File latest_samples_tsv       = "latest_samples"+datestamp+".tsv"            # formerly intermediate_samplewise
		File latest_clusters_tsv      = "latest_clusters"+datestamp+".tsv"           # formerly intermediate_clusterwise
		File unclustered_samples      = "unclustered_samples" + datestamp + ".txt"
		File unclustered_subtrees_etc = "unclustered_subtrees_etc.tar.gz"            # contains subtree assignment information

		# trees and matrices
		File      bigtree_gen                           = "aworkdir000000.nwk"               # generated by cluster script (should match bigtree_raw)
		File      bigtree_matrix                        = "aworkdir000000_dmtrx.tsv"
		File      bigtree_raw                           = "BIGTREE"+datestamp+".nwk"         # generated by matUtils (should match bigtree_gen)
		File      cluster_matrices_randomIDs            = "randomID_cluster_matrices.tar.gz" # formerly Array[File]? acluster_matrices
		File      cluster_subtrees_randomIDs            = "randomID_cluster_trees.tar.gz"    # formerly Array[File]? acluster_trees
		
		# stats for the nerds
		Int n_big_clusters        = read_int("n_big_clusters")
		Int n_samples_in_clusters = read_int("n_samples_in_clusters")
		Int n_samples_processed   = read_int("n_samples_processed")
		Int n_unclustered         = read_int("n_unclustered")
	}
}

task process_CDPH_clusters {
	# Formally cluster_CDPH_method, now split into two tasks for easier debugging and call cache savings.
	# This task runs DOWNSTREAM of find_CDPH_clusters.
	# find_clusters.py: Used to generate distance matrices
	# process_clusters.py: Persistent cluster IDs, subtrees, and MR upload
	# Any clusters that have at least one sample without a diff file will NOT be backmasked
	# This might not work properly if any sample IDs contain a space
	input {
		File input_mat_with_new_samples
		String datestamp # has to be defined here for non-glob delocalization to work properly

		# These come from find_CDPH_clusters WDL task/find_clusters.py
		File latest_samples_tsv
		File latest_clusters_tsv
		File? cluster_matrices_randomIDs_tarball
		File? cluster_subtrees_randomIDs_tarball

		Array[String]? microreact_metadata_columns

		Boolean upload_clusters_to_microreact  = true
		Boolean no_dropped_sample_failsafe     = false
		Boolean only_matrix_special_samples    # arg is assumed to be passed in from Tree Nine
		File? special_samples

		Boolean force_microreact_update        = false
		
		File? persistent_denylist
		File? persistent_ids
		File? persistent_cluster_meta
		File combined_diff_file           # used for local masking
		File? previous_run_cluster_json   # for comparisons -- currently we do this another way so this is unused

		# keep these files in the workspace bucket for now
		File? microreact_decimated_template_json
		File? microreact_update_template_json
		File? microreact_blank_template_json  # hardcoded to expect a file named BLANK_template.json
		File? microreact_key

		# actually optional
		File? sample_metadata_tsv
		String? shareemail
		
		Int preempt = 0 # only set if you're doing a small test run
		Int memory = 50
		Boolean verbose = true
		Boolean DEBUG_generate_debug_mr_jsons = false
		
		# temporary overrides
		File? override_find_clusters_script
		File? override_process_clusters_script
		File? override_summarize_changes_script
		File? override_mass_rename_script
		
	}
	# We cannot `String arg_token = if upload_clusters_to_microreact then "--token ~{microreact_key}" else "" ` or else the literal gs:// will
	# instead of the delocalized version, so some args will need to be handled in the command section itself

	Array[Int] cluster_distances = [20, 10, 5] # CHANGING THIS WILL BREAK SECOND SCRIPT!
	String arg_denylist = if defined(persistent_denylist) then "--dl ~{persistent_denylist}" else ""
	String arg_shareemail = if defined(shareemail) then "-s ~{shareemail}" else ""
	String arg_microreact = if upload_clusters_to_microreact then "--upload_to_microreact" else ""
	String arg_disable_dropped_sample_failsafe = if no_dropped_sample_failsafe then "--no_dropped_sample_failsafe" else ""
	String arg_force_mr_update = if force_microreact_update then "--force_mr_update" else ""
	String arg_verbose = if verbose then "--verbose" else ""
	String arg_debug_MR_jsons = if DEBUG_generate_debug_mr_jsons then "--debug_mr_json" else ""

	# naturally, this doesn't work on Cromwell
	#String? microreact_columns_csv = if defined(microreact_metadata_columns) then sep(",", microreact_metadata_columns) else ""
	
	command <<<
		set -eux pipefail
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting task"
		MICROREACT_COLUMNS_CSV=~{sep="," microreact_metadata_columns}
		
		# additional input validation now handled with a separate WDL task
		if [[ -f "~{microreact_key}" ]]
		then
			TOKEN_ARG="--token ~{microreact_key}"
		else
			TOKEN_ARG=""
		fi

		if [[ -f "~{microreact_update_template_json}" ]]
		then
			MR_UPDATE_JSON_ARG="--mr_update_template ~{microreact_update_template_json}"
		else
			MR_UPDATE_JSON_ARG=""
		fi

		if [[ -f "~{microreact_blank_template_json}" ]]
		then
			MR_BLANK_JSON_ARG="--mr_blank_template ~{microreact_blank_template_json}"
		else
			MR_BLANK_JSON_ARG=""
		fi

		if [[ -f "~{microreact_decimated_template_json}" ]]
		then
			MR_DECIMATED_JSON_ARG="--mr_decimated_template ~{microreact_decimated_template_json}"
		else
			MR_DECIMATED_JSON_ARG=""
		fi

		if [[ -f "~{cluster_subtrees_randomIDs_tarball}" && -f "~{cluster_matrices_randomIDs_tarball}" ]]
		then
			echo "Found cluster subtree and matrix tarballs; can upload to Microreact"
		elif [[ "~{upload_clusters_to_microreact}" = "true" ]]
		then
			echo -n "Upload to microreact is true, but either cluster_subtrees_randomIDs_tarball or cluster_matrices_randomIDs_tarball "
			echo -n "(or both) is missing. Although technically optional if just finding persistent IDs, these files are necessary to show "
			echo -n "subtrees and distance matrices on Microreact. Crashing!"
			exit 1
		fi

		# we do similar logic within process_clusters.py too, but if we can crash before find_clusters.py that'd be ideal
		if [[ -f "~{persistent_ids}" ]]
		then
			if [[ -f "~{persistent_cluster_meta}" ]]
			then
				PERSISTENTIDS_ARG="--persistentids ~{persistent_ids}"
				PERSISTENTMETA_ARG="--persistentclustermeta ~{persistent_cluster_meta}"
			else
				echo "Found persistent IDs file but no persistent cluster meta. You need neither or both. Crashing!"
				exit 1
			fi
		else
			if [[ -f "~{persistent_cluster_meta}" ]]
			then
				echo "Found persistent cluster meta file but no persistent IDs. You need neither or both. Crashing!"
				exit 1
			else
				echo "Found neither persistent IDs file nor persistent cluster meta, will be running without persistent IDs"
				PERSISTENTIDS_ARG=""
				PERSISTENTMETA_ARG=""
			fi
		fi

		if [[ -f "~{sample_metadata_tsv}" ]]
		then
			SAMPLEMETADATA_ARG="--samplemeta ~{sample_metadata_tsv}"
		else
			SAMPLEMETADATA_ARG=""
		fi

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting newick (A_big.nwk) from input PB"
		matUtils extract -i ~{input_mat_with_new_samples} -t A_big.nwk
		cp ~{input_mat_with_new_samples} .

		if [[ -f "~{cluster_matrices_randomIDs_tarball}" ]]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expanding cluster matrix tarball..."
			pigz -dc "~{cluster_matrices_randomIDs_tarball}" | tar xf -
		fi
		
		if [[ -f "~{cluster_subtrees_randomIDs_tarball}" ]]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Expanding cluster subtree tarball..."
			pigz -dc "~{cluster_subtrees_randomIDs_tarball}" | tar xf -
		fi

		if [[ ! "~{override_find_clusters_script}" == '' ]]
		then
			rm /HOME/ash/scripts/find_clusters.py
			mv "~{override_find_clusters_script}" /HOME/ash/scripts/find_clusters.py
		fi

		if [[ ! "~{override_process_clusters_script}" == '' ]]
		then
			rm /HOME/ash/scripts/process_clusters.py
			mv "~{override_process_clusters_script}" /HOME/ash/scripts/process_clusters.py
		fi

		if [[ ! "~{override_summarize_changes_script}" == '' ]]
		then
			rm /HOME/ash/scripts/summarize_changes_alt.py
			mv "~{override_summarize_changes_script}" /HOME/ash/scripts/summarize_changes_alt.py
		fi

		if [[ ! "~{override_mass_rename_script}" == '' ]]
		then
			touch /HOME/ash/scripts/mass_rename_to_persistent_id.py
			rm /HOME/ash/scripts/mass_rename_to_persistent_id.py
			mv "~{override_mass_rename_script}" /HOME/ash/scripts/mass_rename_to_persistent_id.py
		fi

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Files moved if necessary"
		if [[ "~{verbose}" == "true" ]]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Workdir after moving files (disable this print with !verbose)"
			tree
		fi

		CLUSTER_DISTANCES="~{sep=',' cluster_distances}"
		FIRST_DISTANCE="${CLUSTER_DISTANCES%%,*}"
		OTHER_DISTANCES="${CLUSTER_DISTANCES#*,}"
		echo "cluster distances $CLUSTER_DISTANCES"
		echo "First distance $FIRST_DISTANCE"
		echo "Other distances $OTHER_DISTANCES"

		# Turn off pipefail at this point for a few reasons
		# 1) find_clusters.py can return not-0 in non-error cases
		# 2) process_clusters.py writes a lot of logs to disk and we need them if it fails
		set +eo pipefail 

		# TODO: on very large runs, the size of $/samples may eventually cause issues with ARG_MAX
		# should be fine for our purposes though

		if [[ "~{only_matrix_special_samples}" = "true" ]]
		then
			samples=$(< "~{special_samples}" tr -s '\n' ',' | head -c -1)
			ALLSAMPLES_ARG_1="--allsamples"
			ALLSAMPLES_ARG_2="$samples"
		else
			ALLSAMPLES_ARG_1=""
			ALLSAMPLES_ARG_2=""
		fi
		mkdir logs
		mkdir logs/_microreact_jsons_

		echo "Contents of workdir:"
		tree

		# There seems to sometimes be inconsistent behavior r/e handling of whitespace, so
		# this is an overkill "print everything" dump (since stuff defined at runtime will
		# not appear in the WDL "command" file).
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generated these args for process_clusters.py:"
		echo "COMBINED_DIFF_FILE:--combineddiff ~{combined_diff_file}"
		echo "LATEST_SAMPLES_TSV:--latestsamples ~{latest_samples_tsv}"
		echo "LATEST_CLUSTER_TSV:--latestclustermeta ~{latest_clusters_tsv}"
		echo "INPUT_MAT_WITH_NEW_SAMPLES:--mat_tree ~{input_mat_with_new_samples}"
		echo "DATESTAMP:--today ~{datestamp}"
		echo "ARG_DENYLIST:~{arg_denylist}"
		echo "ARG_DISABLE_DROPPED_SAMPLE_FAILSAFE:~{arg_disable_dropped_sample_failsafe}"
		echo "ARG_VERBOSE:~{arg_verbose}"
		echo "PERSISTENTIDS_ARG:$PERSISTENTIDS_ARG"
		echo "PERSISTENTMETA_ARG:$PERSISTENTMETA_ARG"
		echo "SAMPLEMETADATA_ARG:$SAMPLEMETADATA_ARG"
		echo "MICROREACT_COLUMNS_CSV:--mr_metadata_columns $MICROREACT_COLUMNS_CSV"
		echo "ARG_FORCE_MR_UPDATE:~{arg_force_mr_update}"
		echo "ARG_MICROREACT:~{arg_microreact}"
		echo "ARG_SHAREEMAIL:~{arg_shareemail}"
		echo "ARG_DEBUG_MR_JSONS:~{arg_debug_MR_jsons}"
		echo "MR_UPDATE_JSON_ARG:$MR_UPDATE_JSON_ARG"
		echo "MR_BLANK_JSON_ARG:$MR_BLANK_JSON_ARG"
		echo "MR_DECIMATED_JSON_ARG:$MR_DECIMATED_JSON_ARG"
		echo "TOKEN_ARG:$TOKEN_ARG"
		echo "ALL_SAMPLES:$ALLSAMPLES_ARG_1 $ALLSAMPLES_ARG_2"
		
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running process_clusters.py"

		# shellcheck disable=SC2086
		python3 /HOME/ash/scripts/process_clusters.py \
			--combineddiff "~{combined_diff_file}" \
			--latestsamples "~{latest_samples_tsv}" \
			--latestclustermeta "~{latest_clusters_tsv}" \
			--mat_tree "~{input_mat_with_new_samples}" \
			--today ~{datestamp} \
			~{arg_denylist} \
			~{arg_disable_dropped_sample_failsafe} \
			~{arg_verbose} \
			$PERSISTENTIDS_ARG \
			$PERSISTENTMETA_ARG \
			$SAMPLEMETADATA_ARG \
			--mr_metadata_columns $MICROREACT_COLUMNS_CSV \
			~{arg_force_mr_update} \
			~{arg_microreact} \
			~{arg_shareemail} \
			~{arg_debug_MR_jsons} \
			$MR_UPDATE_JSON_ARG \
			$MR_BLANK_JSON_ARG \
			$MR_DECIMATED_JSON_ARG \
			$TOKEN_ARG \
			$ALLSAMPLES_ARG_1 $ALLSAMPLES_ARG_2

		PY_EXIT_CODE=$? # might intermittently fails on Terra (it happened once and now I'm scared)

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Zipping process_clusters.py's logs"
		zip -r logs.zip ./logs
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Logs zipped"

		if [[ "~{verbose}" == "true" ]]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Workdir after process_clusters.py (disable this print with !verbose)"
			tree
		fi

		# if process_clusters.py errored, NOW we should crash, since we have logs and such
		if [ "$PY_EXIT_CODE" -ne 0 ]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Crashing with rc $PY_EXIT_CODE because that's what process_clusters.py returned"
			exit $PY_EXIT_CODE
		fi

		if [ "~{previous_run_cluster_json}" != "" ]
		then
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running summarize_changes_alt.py"
			python3 /HOME/ash/scripts/summarize_changes_alt.py "all_cluster_information~{datestamp}.json"
			echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished summarize_changes_alt.py"
		fi

		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running mass_rename_to_persistent_id.py"
		python3 /HOME/ash/scripts/mass_rename_to_persistent_id.py "~{arg_verbose}" --json "all_cluster_information~{datestamp}.json"
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished mass_rename_to_persistent_id.py"

		echo "The IDs of these clusters were processed by process_clusters.py on ~{datestamp}(ish) and DO account for persistent cluster IDs. " > readme.txt
		echo "Note datestamp is set at execution of first task to ensure all outputs have same datestamp per workflow run, hence -ish. " >> readme.txt

		if [ ~{verbose} = "true" ]; then tree; fi

		find . -maxdepth 1 \( -name "a*.nwk" -o -name "a*.pb" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > "persisID_cluster_trees~{datestamp}.tar.gz"
		find . -maxdepth 1 \( -name "b*.nwk" -o -name "b*.pb" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > "persisID_cluster_trees_backmasked~{datestamp}.tar.gz"
		find . -maxdepth 1 \( -name "a*_dmtrx" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > "persisID_cluster_matrices~{datestamp}.tar.gz"
		find . -maxdepth 1 \( -name "b*_dmtrx" -o -name "readme.txt" \) -print0 | tar -cf - --null -T - | pigz -1 > "persisID_cluster_matrices_backmasked~{datestamp}.tar.gz"

		# shellcheck disable=SC2317
		echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished task"

	>>>

	runtime {
		bootDiskSizeGb: 15
		cpu: 12
		disks: "local-disk " + 150 + " SSD"
		docker: "ashedpotatoes/usher-plus:0.6.6_rev16"
		memory: memory + " GB"
		preemptible: preempt
	}

	output {
		# The amount of outputs we originally had was overloading Terra, so some of these are commented out now.
		# Also, we try to avoid globbing where possible to make finding outs in Terra bucket easier since globs
		# create a folder with a randomized name, which is annoying!

		###### IMPORTANT FILES THAT SHOULD ALWAYS GO INTO SUBSEQUENT RUNS IF THEY EXIST ######
		File new_persistent_ids = "persistentIDS" + datestamp + ".tsv"
		File new_persistent_meta = "persistentMETA" + datestamp + ".tsv"
		File final_cluster_information_json = "all_cluster_information" + datestamp + ".json"

		# cluster-specific subtrees and nwks
		# there is no internally masked big tree because masking is done per-cluster
		File cluster_trees_persisIDs = "persisID_cluster_trees" + datestamp + ".tar.gz"                        
		File cluster_trees_persisIDs_backmasked = "persisID_cluster_trees_backmasked" + datestamp + ".tar.gz"
		File cluster_matrices_persisIDs = "persisID_cluster_matrices" + datestamp + ".tar.gz"
		File cluster_matrices_persisIDs_backmasked = "persisID_cluster_matrices_backmasked" + datestamp + ".tar.gz"

		# stats for the nerds, currently unused downstream
		File new_samples_cluster_information = "new_samples" + datestamp + ".tsv"
		File all_samples_cluster_information = "all_samples" + datestamp + ".tsv"
		File? decimated_clusters = "decimated" + datestamp + ".tsv"
		File? updated_mr_URIs_file = "updated_mr_URIs" + datestamp + ".txt"

		# debug
		File? logs = "logs.zip"
		File? change_report_json       = "change_report" + datestamp + ".json"
		File? change_report_full       = "change_report_full"+datestamp+".txt"  # all clusters
		File? change_report_cdph       = "change_report_cdph"+datestamp+".txt"  # excludes 20-clusters
		#File? input_metadata_tsv       = sample_metadata_tsv                    # because metadata is mutable on Terra

		# can be used to annotate the Big Tree by cluster; uses Nextstrain (Auspice) metadata format
		# we currently don't annotate with Nextstrain as CDPH prefers Microreact, but the files are here if you want 'em
		File? samp_cluster_twn = "samp_persis20cluster" + datestamp + ".tsv"
		File? samp_cluster_ten = "samp_persis10cluster" + datestamp + ".tsv"
		File? samp_cluster_fiv = "samp_persis5cluster"  + datestamp + ".tsv"
	}
}

task validate_treenine_inputs {
	input {
		File? input_tree
		File? existing_diffs          # in theory we could just localize the string name on Terra but in practice it's iffy
		File? existing_samples
		File? persistent_cluster_meta
		File? persistent_cluster_ids
		File? previous_run_cluster_json

		Boolean adhoc
		Boolean identify_clusters
		Boolean restart_clusters
		Boolean upload_clusters_to_microreact

		File? microreact_blank_template_json
		File? microreact_decimated_template_json
		File? microreact_key
		File? microreact_update_template_json
		Boolean DEBUG_generate_debug_mr_jsons
		
		File? ref_genome

		# not checked
		#Array[File] diffs
		#String? listener_bucket
		#Boolean detailed_clades
		#Float?  max_low_coverage_sites
		#File?   matutils_clade_annotations
		#Boolean optimize = true
		#String? reroot_to_this_node
		#Boolean summarize_tree_before_placing_samples
		#Boolean summarize_tree_after_placing_samples
		#Boolean identify_clusters
		#Boolean cluster_entire_tree
		#File? special_samples
		#File? persistent_denylist
		#File? sample_metadata_tsv
		#Boolean strictly_check_metadata
		#String? microreact_share_email
		#String? microreact_share_team
		#Array[File]? coverage_reports
		#String? comment
		#Array[String]? rename_samples
		#Boolean datestamp_outs
		#String out_prefix
		#String out_diffs
		#Boolean DEBUG_concat_files_then_exit
		#File?   DEBUG_override_latest_samples
		#File?   DEBUG_override_latest_clusters
		
	}

	command <<<
		set -eux pipefail

		echo "CHECKING INPUT_TREE: ~{input_tree}"
		if [[ -f "~{input_tree}" ]]
		then
			echo "Found an input_tree to use as our base tree"
			if [[ "~{input_tree}" != *.pb ]]
			then
				echo "ERROR: input_tree does not end in .pb ergo is likely the wrong file format"
				exit 1
			fi
		else
			echo "WARNING: No input tree, will use a hardcoded fallback"
		fi

		echo "Inputs:"
		echo "adhoc = ~{adhoc}"
		echo "identify_clusters = ~{identify_clusters}"
		echo "restart_clusters = ~{restart_clusters}"
		echo "upload_clusters_to_microreact = ~{upload_clusters_to_microreact}"

		if [[ "~{adhoc}" = "true" && ~{upload_clusters_to_microreact} = "true" ]]
		then
			echo "ERROR: adhoc is true, but upload_clusters_to_microreact also true, not clear what user wants to do"
			echo "If you're attempting to restart cluster IDs, turn off adhoc and set restart_clusters to true"
			exit 1
		fi 

		if [[ "~{restart_clusters}" = "true" && ~{identify_clusters} = "false" ]]
		then
			echo "ERROR: restart_clusters is true, identify_clusters is false. can't restart clusters if no clusters!"
			exit 1
		fi 


		echo "CHECKING THE ZERO, TWO, OR FIVE PERSISTENT FILES"
		EXISTING_DIFFS="~{existing_diffs}"
		EXISTING_SAMPLES="~{existing_samples}"
		PERSISTENT_IDS="~{persistent_cluster_ids}"
		PERSISTENT_META="~{persistent_cluster_meta}"
		PREVOUS_CLUSTER_JSON="~{previous_run_cluster_json}"  # not technically required but needed for manual change reporting

		EXISTING_DIFFS_exists=$([[ -f "$EXISTING_DIFFS" ]] && echo 1 || echo 0)
		EXISTING_SAMPLES_exists=$([[ -f "$EXISTING_DIFFS" ]] && echo 1 || echo 0)
		sum_kingfiles=$(( EXISTING_DIFFS_exists + EXISTING_SAMPLES_exists))

		if [[ $sum_kingfiles = 1 ]]
		then
			echo "ERROR: EXISTING_DIFFS (.diff) and EXISTING_SAMPLES (no ext) must either both exist or both be missing"
			exit 1
		fi
		if [[ "~{adhoc}" = "true" && $sum_kingfiles -ne 0 ]]
		then
			echo "ERROR: adhoc is true, but EXISTING_DIFFS or EXISTING_SAMPLES (or both) was defined"
			exit 1
		fi

		PERSISTENT_IDS_exists=$([[ -f "$PERSISTENT_IDS" ]] && echo 1 || echo 0)
		PERSISTENT_META_exists=$([[ -f "$PERSISTENT_META" ]] && echo 1 || echo 0)
		PREVOUS_CLUSTER_JSON_exists=$([[ -f "$PREVOUS_CLUSTER_JSON" ]] && echo 1 || echo 0)
		sum_persistent_files=$(( PERSISTENT_IDS_exists + PERSISTENT_META_exists + PREVOUS_CLUSTER_JSON_exists ))

		if [[ "~{adhoc}" = "true" && $sum_persistent_files -ne 0 ]]
		then
			echo "ERROR: adhoc is true, but PERSISTENT_IDS or PERSISTENT_META or PREVIOUS_CLUSTER_JSON (or some combo thereof) was defined, see docs on adhoc runs"
			exit 1
		fi

		if [[ $sum_persistent_files -ne 0 && $sum_persistent_files -ne 3 ]]
		then
			echo "ERROR: PERSISTENT_IDS (.tsv), PERSISTENT_META (.tsv), and PREVIOUS_CLUSTER_JSON (.json/.ndjson) must either ALL exist or ALL be missing"
			exit 1
		fi

		# extensions
		if [[ -f "$EXISTING_DIFFS" && "$EXISTING_DIFFS" != *.diff ]]
		then
			echo "ERROR: EXISTING_DIFFS exists but does not have .diff after its datestamp: $EXISTING_DIFFS"
			exit 1
		fi
		if [[ -f "$PERSISTENT_IDS" && "$PERSISTENT_IDS" != *.tsv ]]
		then
			echo "ERROR: PERSISTENT_IDS exists but does not have .tsv after its datestamp: $PERSISTENT_IDS"
			exit 1
		fi
		if [[ -f "$PERSISTENT_META" && "$PERSISTENT_META" != *.tsv ]]
		then
			echo "ERROR: PERSISTENT_META exists but does not have .tsv after its datestamp: $PERSISTENT_META"
			exit 1
		fi
		if [[ -f "$PREVOUS_CLUSTER_JSON" && "$PREVOUS_CLUSTER_JSON" != *.json && "$PREVOUS_CLUSTER_JSON" != *.ndjson ]]
		then
			echo "ERROR: PREVOUS_CLUSTER_JSON exists but does not have .json or .ndjson after its datestamp: $PREVOUS_CLUSTER_JSON"
			exit 1
		fi

		# datestamp consistency
		target_date=""
		for f in "$EXISTING_DIFFS" "$EXISTING_SAMPLES" "$PERSISTENT_IDS" "$PERSISTENT_META" "$PREVOUS_CLUSTER_JSON"
		do
			if [[ -f "$f" ]]
			then
				filename="${f##*/}"
				no_ext="${filename%.*}"
				
				if [[ "$no_ext" =~ ([0-9-]+)$ ]]
				then
					current_date="${BASH_REMATCH[1]}"
				else
					current_date="${no_ext##*.}"
				fi

				if [[ -z "$target_date" ]]
				then
					target_date="$current_date"
					echo "Baseline datestamp set to: ${target_date}"
				elif [[ "$current_date" != "$target_date" ]]
				then
					echo "ERROR: Datestamp mismatch found!"
					echo "Expected: '${target_date}'"
					echo "Found '${current_date}' in file: ${f}"
					exit 1
				fi
			fi
		done
		
		echo "CHECKING MICROREACT STUFF"
		if [[ "~{upload_clusters_to_microreact}" = "true" ]]
		then
			echo "Upload to Microreact (upload_clusters_to_microreact) is true"
			if [[ "~{DEBUG_generate_debug_mr_jsons}" = "true" ]]
			then
				# breaking Microreact clusters is a destructive action so we have to be kinda strict on this
				echo "ERROR: Do not set DEBUG_generate_debug_mr_jsons to true if upload_clusters_to_microreact is true! Either both should be true or both should be false!"
				exit 1
			elif [[ ! -f "~{microreact_key}" ]]
			then
				echo "ERROR: Upload to microreact is true, but no token provided"
				exit 1
			fi
		fi
		if [[ "~{upload_clusters_to_microreact}" = "true" || "~{DEBUG_generate_debug_mr_jsons}" = "true" ]]
		then
			if [[ ! -f "~{microreact_update_template_json}" ]]
			then
				echo "ERROR: Upload to microreact or DEBUG_generate_debug_mr_jsons is true, but no microreact_update_template_json provided"
				exit 1
			fi
			if [[ ! -f "~{microreact_blank_template_json}" ]]
			then
				echo "ERROR: Upload to microreact or DEBUG_generate_debug_mr_jsons is true, but no microreact_blank_template_json provided"
				exit 1
			fi
			if [[ ! -f "~{microreact_decimated_template_json}" ]]
			then
				# not strictly required but should be done to avoid footguns
				echo "ERROR: Upload to microreact or DEBUG_generate_debug_mr_jsons is true, but no microreact_decimated_template_json provided"
				exit 1
			fi
		fi

		echo "CHECKING REF_GENOME: ~{ref_genome}"
		if [[ -f "~{ref_genome}" ]]
		then
			echo "WARNING: Found a ref_genome. If you're running this for H37Rv tuberculosis, don't do that, just use the default fallback!"
		else
			echo "No ref_genome provided, will use hardcoded H37Rv (this is recommended)"
		fi
	>>>

	runtime {
		cpu: 2
		disks: "local-disk " + 200 + " SSD"
		docker: "debian:trixie-20260713-slim"
		memory: 4 + " GB"
		preemptible: 2
	}

	output {
		Boolean didnt_crash = true
	}
}