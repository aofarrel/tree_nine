version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Convert_to_Nextstrain {
	input {
        File pb_tree
        Array[File?] metadata_files
        String? outfile = "nextstrain"

        Boolean subtrees = true
		Int subtree_treesize = 10
		Int subtree_nearest_k = 5
        File? subtree_only_these_samples
    }

    String out = select_first([outfile, basename(pb_tree, ".pb")])

    if(subtrees) {
        call treenine.convert_to_nextstrain_subtrees as many_trees {
            input:
                input_mat = pb_tree,
                selected_samples = subtree_only_these_samples,
                treesize = subtree_treesize,
                nearest_k = subtree_nearest_k,
                outfile_nextstrain = out,
                metadata_files = metadata_files
        }
    }
    if(!subtrees) {
        call treenine.convert_to_nextstrain_single as one_tree {
            input:
                input_mat = pb_tree,
                outfile_nextstrain = out,
                metadata_files = metadata_files
        }
    }
}