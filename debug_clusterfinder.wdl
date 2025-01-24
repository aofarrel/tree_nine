version 1.0

import "./matutils_and_friends.wdl" as matWDLlib

workflow DebugClusterScript {

	input {
		File input_mat
		File firstlines_from_cat_diff 
		File find_clusters_script_override
		File process_clusters_script_override
		File persistent_ids
		File bogus
		Boolean cluster_everything = false
	}

	File persistent_cluster_tsv = bogus
	File microreact_key = bogus
	File microreact_template_json = bogus

	
	call matWDLlib.nwk_json_cluster_matrix_microreact as normal_clusters {
		input:
			input_mat = input_mat,
			special_samples = firstlines_from_cat_diff,
			only_matrix_special_samples = !(cluster_everything),
			microreact_key = microreact_key,
			microreact_template_json = microreact_template_json,
			find_clusters_script_override = find_clusters_script_override,
			persistent_ids = persistent_ids,
			process_clusters_script_override = process_clusters_script_override
	}
}