version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Nextstrainify {
	input {
        File pb_tree
        Array[File]? metadata_files
        Boolean subtrees = true

		Int treesize = 10
		Int nearest_k = 5
        File? only_these_samples
		String outfile = "nextstrain"
    }

    if(defined(only_these_samples)) { Boolean new_samples_only_1 = true }
    if(!defined(only_these_samples)) { Boolean new_samples_only_2 = false }

    if(subtrees) {
        call treenine.convert_to_nextstrain_subtrees as many_trees {
            input:
                input_mat = pb_tree,
                new_samples = only_these_samples,
                new_samples_only = select_first([new_samples_only_1, new_samples_only_2]),
                treesize = treesize,
                nearest_k = nearest_k,
                outfile_nextstrain = outfile,
                metadata_files = metadata_files
        }
    }
    if(!subtrees) {
        call treenine.convert_to_nextstrain_single as one_tree {
            input:
                input_mat = pb_tree,
                outfile_nextstrain = outfile,
                metadata_files = metadata_files
        }
    }
}