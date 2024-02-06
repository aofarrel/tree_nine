version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Summarize_UShER_Tree {
    input {
        File pb_tree
    }

    String basename_input_mat = basename(pb_tree, ".pb")

    call treenine.summarize as get_info {
        input:
            input_mat = pb_tree,
            prefix_outs = basename_input_mat
    }

    output {
        File summary = get_info.summary
		File samples = get_info.samples
		File clades = get_info.clades
		File mutations = get_info.mutations
		File aberrant = get_info.aberrant
		File haplotype = get_info.haplotype
		File sample_clades = get_info.sample_clades
		File calculate_roho = get_info.calculate_roho
	}
}