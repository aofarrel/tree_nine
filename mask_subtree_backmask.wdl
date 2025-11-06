version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Mask_Subtree_by_Shared_Masks {
    input {
        File pb_tree
        File combined_diff_of_desired_samples
        File desired_samples
        Boolean calc_distance_matrix = true
    }
    
    String basename = basename(pb_tree, ".pb")
    
    call treenine.extract as extract {
        input:
            input_mat = pb_tree,
            samples = desired_samples
    }
    
    if(calc_distance_matrix) {
        call treenine.convert_to_newick as convert_subtree_to_nwk {
            input:
                input_mat = extract.subtree,
                outfile_nwk = basename + ".nwk"
        }
        call treenine.matrix as matrix_subtree {
            input:
                input_nwk = convert_subtree_to_nwk.newick_tree,
                only_matrix_special_samples = false               # since it's already a subtree
        }
    }

    call treenine.backmask as backmask {
        input:
            input_mat = extract.subtree,
            combined_diff = combined_diff_of_desired_samples
    }
    
    if(calc_distance_matrix) {
        call treenine.convert_to_newick as convert_masked_to_nwk {
            input:
                input_mat = backmask.backmasked_tree,
                outfile_nwk = basename + ".nwk"
        }
        call treenine.matrix as matrix_masked {
            input:
                input_nwk = convert_masked_to_nwk.newick_tree,
                only_matrix_special_samples = false               # since it's already a subtree
        }
    }

    output {
        File unmasked_subtree = extract.subtree
        File masked_tree = backmask.backmasked_tree
        File? matrix_of_subtree = matrix_subtree.out_matrix
        File? matrix_of_masked_tree = matrix_masked.out_matrix
    }
}