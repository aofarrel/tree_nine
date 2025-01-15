version 1.0

import "./matutils_and_friends.wdl" as matWDLlib

workflow DebugClusterScript {

	input {
		File bogus 
		File cluster_main_script_debug_override
		Boolean cluster_everything = false
	}

	File final_maximal_output_tree = bogus
	File persistent_cluster_tsv = bogus
	File microreact_key = bogus
	File microreact_template_json = bogus
	File special_samples_added = bogus

	
	call matWDLlib.nwk_json_cluster_matrix_microreact as normal_clusters {
		input:
			input_mat = final_maximal_output_tree,
			special_samples = special_samples_added,
			only_matrix_special_samples = !(cluster_everything),
			persistent_cluster_tsv = persistent_cluster_tsv,
			microreact_key = microreact_key,
			microreact_template_json = microreact_template_json,
			cluster_main_script_debug_override = cluster_main_script_debug_override
	}
}