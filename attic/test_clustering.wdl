version 1.0

import "../matutils_and_friends.wdl" as matWDLlib

workflow Test_Clustering {
	input {
		File input_tree
		Boolean cluster_entire_tree = true

		File? find_clusters_script_override
		File? process_clusters_script_override
		String? date = "1970-01-01"

		File microreact_update_template_json
		File microreact_blank_template_json
	}

	call matWDLlib.cluster_CDPH_method as cluster {
		input:
			input_mat_with_new_samples = input_tree,
			combined_diff_file = input_tree,
			only_matrix_special_samples = !(cluster_entire_tree),
			persistent_ids = input_tree,
			persistent_cluster_meta = input_tree,
			upload_clusters_to_microreact = false,
			find_clusters_script_override = find_clusters_script_override,
			process_clusters_script_override = process_clusters_script_override,
			microreact_update_template_json = microreact_update_template_json,
			microreact_blank_template_json = microreact_blank_template_json
	}
}