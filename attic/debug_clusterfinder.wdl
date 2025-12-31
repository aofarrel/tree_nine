version 1.0

import "../matutils_and_friends.wdl" as matWDLlib

workflow DebugClusterScript {

	input {
		File input_mat
		File firstlines_from_cat_diff
		File cat_diff
		File find_clusters_script_override
		File process_clusters_script_override
		File summarize_changes_script_override
		File microreact_blank_template_json
		File microreact_update_template_json
		File microreact_key
		File persistent_ids
		File persistent_cluster_meta
		File previous_run_cluster_json
		String today = "1960-01-01"
		Boolean cluster_everything = true
	}

	
	call matWDLlib.cluster_CDPH_method as do_everything {
		input:
			input_mat_with_new_samples = input_mat,
			special_samples = firstlines_from_cat_diff,
			combined_diff_file = cat_diff,
			only_matrix_special_samples = !(cluster_everything),
			microreact_key = microreact_key,
			find_clusters_script_override = find_clusters_script_override,
			microreact_blank_template_json = microreact_blank_template_json,
			microreact_update_template_json = microreact_update_template_json,
			persistent_ids = persistent_ids,
			process_clusters_script_override = process_clusters_script_override,
			summarize_changes_script_override = summarize_changes_script_override,
			persistent_cluster_meta = persistent_cluster_meta,
			previous_run_cluster_json = previous_run_cluster_json,
			datestamp = today
	}
}