version 1.0
import "./matutils_and_friends.wdl" as treenine

workflow Mask_Subtree_by_Position {
    input {
        File pb_tree
        File mask_tsv
        File highly_clonal_samples
        Boolean calc_distance_matrix = true
    }
    
    String basename = basename(pb_tree, ".pb")
    
    call treenine.extract as extract {
        input:
            input_mat = pb_tree,
            samples = highly_clonal_samples
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
                only_matrix_special_samples = false
        }
    }

    call treenine.mask as mask {
        input:
            input_mat = extract.subtree,
            mask_tsv = mask_tsv
    }
    
    if(calc_distance_matrix) {
        call treenine.convert_to_newick as convert_masked_to_nwk {
            input:
                input_mat = mask.masked_tree,
                outfile_nwk = basename + ".nwk"
        }
        call treenine.matrix as matrix_masked {
            input:
                input_nwk = convert_masked_to_nwk.newick_tree,
                only_matrix_special_samples = false
        }
    }

    output {
        File unmasked_subtree = extract.subtree
        File masked_tree = mask.masked_tree
        File? matrix_of_subtree = matrix_subtree.out_matrix
        File? matrix_of_masked_tree = matrix_masked.out_matrix
    }
}